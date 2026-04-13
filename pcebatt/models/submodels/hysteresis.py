import numpy as np
from numpy.typing import NDArray
from abc import abstractmethod
from collections.abc import Iterable
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
from ...utils import create_logger
from ..base import Model
from .ocv import qOCV_POLY
from types import SimpleNamespace


class HYSTERESIS(Model):
    """
    Template for an OCV object. Must have a solve method, which must return the simulated voltage array and the standard deviation of its uncertainty. If the model is deterministic, return zeros for the std. 
    
    """

    def __init__(self, config = None):
        super().__init__(config = config)
        self.model_type = 'hysteresis'
        self.voltage = SimpleNamespace()

    @abstractmethod
    def solve(self, *args, **kwargs) -> Iterable[NDArray, NDArray]:
        pass


class Plett_Hysteresis(HYSTERESIS):
    """
    Simulates hysteresis effects with the help of a one-state decay model. 
    Credit: G.L.Plett "Battery Management Systems Vol. 1, Battery Modeling" 
    """

    def __init__(self, g: float = 0.05, logger_level: str = 'DEBUG', config = None) -> None:
        super().__init__(config=config)

        self.ocv_model = qOCV_POLY(logger_level="CRITICAL", config=config) # model needs ocv charge and discharge
        self.g = g  # normalized hysteresis decay rates (g = gamma/C_Nom)
        self.max_solver_step = 20.
        self.logger = create_logger(__class__.__name__, logger_level)

        self.logger.info(f'[{__class__.__name__}] Initialized!')


    def solve(self, **kwargs) -> Iterable[NDArray, NDArray]:
        """
        Parameters
        ----------
            **kwargs : dict of NDArray
                Input arrays. Keys: 'soc', 'T', 'time', 'current'.
        """

        current = kwargs['current']
        time = kwargs['time'] 
        soc = kwargs['soc']
        T = kwargs['T']

        self.ocv_model.alpha = 1 # set ocv to charge mode
        ocv_cha = self.ocv_model.solve(**dict(soc=soc, T=T))[0]
        self.ocv_model.alpha = 0 # set ocv to discharge mode
        ocv_dch = self.ocv_model.solve(**dict(soc=soc, T=T))[0]
        M = np.sign(current)*(ocv_cha-ocv_dch)/2 # maximum ocv polarization

        # interpolators for input arrays
        current_ip = interp1d(time, current, kind='linear', bounds_error=False, fill_value=(current[0], current[-1]))
        M_ip = interp1d(time, M, kind='linear', bounds_error=False, fill_value=(M[0], M[-1]))

        # Right-hand-side for the differential equation
        def rhs(t:float, h:NDArray) -> NDArray:
            # Get interpolated soc, current and temperature
            current_t = current_ip(t)
            M_t = M_ip(t)

            a = np.abs(current_t*self.g)
            
            return a*(-h + M_t)

        voltage = solve_ivp(
                        fun=rhs, max_step=self.max_solver_step, y0=np.array([0.]),
                        t_span=(time[0], time[-1]), t_eval=time
                        ).y.T
        voltage = voltage.flatten() 
        
        self.voltage.mean, self.voltage.std = voltage, np.zeros_like(voltage)
        return self.voltage.mean, self.voltage.std
    