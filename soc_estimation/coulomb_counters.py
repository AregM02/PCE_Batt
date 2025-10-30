import sys
import chaospy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from collections.abc import Iterable
from utils import create_logger


class CoulombCounterGalerkin():
    """
    Calculates SoC and its uncertainty with Galerkin PCE.

    """

    def __init__(self, initial_soc: float, capacity: float, 
                 initial_soc_unc: float, capacity_unc: float, 
                 logger_level: str = 'DEBUG') -> Iterable[np.array, np.array]:

        self.logger = create_logger(__class__.__name__, level = logger_level)
        self.max_solver_step = 1. # set appropriate step size

        # ----PCE----
        # propagate uncertainty
        mu_c_inv = 1/capacity
        sigma_c_inv = capacity_unc / capacity**2  
        self.joint = chaospy.J(chaospy.Normal(mu_c_inv, sigma_c_inv),
                               chaospy.Normal(initial_soc, initial_soc_unc))
        expansion, norms = chaospy.generate_expansion(3, self.joint, retall=True)
        self.expansion = expansion
        c_inv, z_0 = chaospy.variable(2)

        self.b = chaospy.E(c_inv*expansion, self.joint) / norms
        self.initial_val = chaospy.E(z_0*expansion, self.joint) / norms
        
        self.logger.info(f"[{__class__.__name__}] Initialized!")


    def __call__(self, current: np.array, time: np.array) -> np.ndarray:

        current_ip = interp1d(time, current, kind='linear', bounds_error=False, fill_value=(current[0], current[-1]))

        def rhs(t:float, x:np.ndarray) -> np.ndarray:
            current_t = current_ip(t)
            return current_t * self.b / 3600
        
        self.logger.info(f"[{__class__.__name__}] Computing SoC sequence...")
        coefficients = solve_ivp(fun=rhs, max_step=self.max_solver_step, y0=self.initial_val,
                                 t_span=(time[0], time[-1]), t_eval=time).y.T
        
        z_approx = chaospy.sum(self.expansion * coefficients, -1)
        mean = chaospy.E(z_approx, self.joint)
        variance = np.abs(chaospy.Var(z_approx, self.joint))
        sigma = np.sqrt(variance)

        self.logger.info(f"[{__class__.__name__}] Computation of SoC completed!")

        return mean, sigma