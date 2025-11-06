import math
from pathlib import Path
import pandas as pd
import numpy as np
from abc import abstractmethod
from collections import defaultdict
from collections.abc import Iterable
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import chaospy
from utils import create_logger
from ..base import Model


PATH_QOCV = Path(__file__).parent.parent.parent / 'data' / 'measurements' / 'qocv20'


class OCV(Model):
    """
    Template for an OCV object. Must have a solve method, which must return the simulated voltage array and the standard deviation of its uncertainty. If the model is deterministic, return zeros for the std. 
    
    """

    @abstractmethod
    def solve(self, *args, **kwargs) -> Iterable[np.ndarray, np.ndarray]:
        pass


class qOCV_POLY(OCV):
    """
    Simulates a polynomial relationship between qOCV and SOC. Looks for qOCV data at <data/measurements/qocv20>
    
    """

    def __init__(self, order: int = 6, alpha: int = 0.5, logger_level: str = 'DEBUG', logger_name=None) -> None:
        super().__init__()
        
        if logger_name is None:
            self.logger = create_logger(__class__.__name__, logger_level)
        else: 
            self.logger = create_logger(logger_name, logger_level)

        self.functions = defaultdict(list)
        self.temperatures = np.array([15, 25, 35, 45])
        self.alpha = alpha # 0.5 = charge/discharge OCVs are averaged; 1 = only charge; 0 = only discharge

        for temp in self.temperatures:
            df_c = pd.read_parquet(PATH_QOCV/f'qocv_{temp}deg_charge.parquet')
            df_d = pd.read_parquet(PATH_QOCV/f'qocv_{temp}deg_discharge.parquet')
            
            poly_c = np.poly1d(np.polyfit(df_c.SOC, df_c.Spannung, order))
            poly_d = np.poly1d(np.polyfit(df_d.SOC, df_d.Spannung, order))
            
            self.functions[temp] = [poly_c, poly_d]

        self.logger.info(f"[{self.logger.name}] Initialized!")


    def solve(self, **kwargs) -> Iterable[np.ndarray, np.ndarray]:
        """
        Parameters
        ----------
            **kwargs : dict of np.ndarray
                Input arrays. Keys: 'soc', 'T'.
        """
        soc = kwargs['soc'][0]
        T = kwargs['T']

        # calculate OCV curves for all defined temperatures
        ocv_all = np.stack([
            (self.alpha*self.functions[temp][0](soc) + (1-self.alpha)*self.functions[temp][1](soc))
            for temp in self.temperatures
        ], axis=0)  # shape: (n_temps, n_soc)


        # clip T to bounds
        T_clip = np.clip(T, self.temperatures[0], self.temperatures[-1])

        # Find interval indices for T
        idx = np.searchsorted(self.temperatures, T_clip, side='right')
        idx = np.clip(idx, 1, len(self.temperatures)-1)  # ensure valid interval
        idx0 = idx - 1
        idx1 = idx

        # Gather OCV values at the interval endpoints
        ocv0 = ocv_all[idx0, np.arange(len(soc))]
        ocv1 = ocv_all[idx1, np.arange(len(soc))]
        t0 = self.temperatures[idx0]
        t1 = self.temperatures[idx1]

        # Linear interpolation    
        ocv_interp = ocv0 + (ocv1 - ocv0) * (T_clip - t0) / (t1 - t0)
  
        return ocv_interp, np.zeros_like(ocv_interp)


class Plett_Hysteresis(OCV):
    """
    Simulates hysteresis effects with the help of a one-state decay model. 
    Credit: G.L.Plett "Battery Management Systems Vol. 1, Battery Modeling" 
    """

    def __init__(self, g: float = 0.05, logger_level: str = 'DEBUG') -> None:
        super().__init__()

        self.ocv_model = qOCV_POLY(logger_level="CRITICAL") # model needs ocv charge and discharge
        self.g = g  # normalized hysteresis decay rates (g = gamma/C_Nom)
        self.max_solver_step = 20.
        self.logger = create_logger(__class__.__name__, logger_level)

        self.logger.info(f'[{__class__.__name__}] Initialized!')


    def solve(self, **kwargs) -> Iterable[np.ndarray, np.ndarray]:
        """
        Parameters
        ----------
            **kwargs : dict of np.ndarray
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
        def rhs(t:float, h:np.ndarray) -> np.ndarray:
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
        
        return voltage, np.zeros_like(voltage)
    

class qOCV_PCE(qOCV_POLY):
    """
    A polynomial chaos expanded variant of qOCV_POLY. Looks for qOCV data at <data/measurements/qocv20>
    
    """


    def __init__(self, order: int = 6, alpha: int = 0.5, logger_level: str = 'DEBUG') -> None:
        super().__init__(order=order, alpha=alpha, logger_level=logger_level, logger_name=__class__.__name__)
        
        # ----- Set up the A matrix -----
        H, norms = chaospy.expansion.hermite(order, retall=True)
        H = H/norms

        A = np.zeros((order+1, order+1))
        for n, Hn in enumerate(H):
            coeffs, exponents = Hn.coefficients, Hn.exponents.flatten()
            A[exponents, n] = coeffs

        # check matrix conditioning
        cond = np.linalg.cond(A)
        if cond > 1e12:
            # warn; use pseudo-inverse for stability
            self.A_inv = np.linalg.pinv(A)
            self.logger.warning("Warning: A is ill-conditioned (cond={cond:.3e}), using pinv")
        else:
            self.A_inv = np.linalg.inv(A)
        


    def get_beta(self, alpha:np.ndarray, mu_t: np.ndarray, sigma_t:np.ndarray) -> np.ndarray:
        """
        alpha: array shape (N+1,) with alpha[n] for x^n (ascending)
        mu_t, sigma_t: arrays shape (T,)
        returns beta: ndarray shape (N+1, T) with beta[k,t]
        """
        alpha = np.asarray(alpha, dtype=float)
        mu_t = np.asarray(mu_t, dtype=float)
        sigma_t = np.asarray(sigma_t, dtype=float)
        N = alpha.size - 1
        T = mu_t.size
        assert mu_t.shape == sigma_t.shape

        # precompute binomial coefficients C(n,k) for n,k in [0 ... N]
        binom = np.zeros((N+1, N+1), dtype=float)
        for n in range(N+1):
            for k in range(n+1):
                binom[n, k] = math.comb(n, k)

        beta = np.zeros((N+1, T), dtype=float)

        # loop over k (output power) and n (terms in alpha) - triangular
        for k in range(N+1):
            for n in range(k, N+1):
                # compute for all t: C(n,k) * alpha[n] * mu_t**(n-k) * sigma_t**k
                # use broadcasting to form a (T,) vector
                term = binom[n, k] * alpha[n] * (mu_t ** (n - k)) * (sigma_t ** k)
                beta[k, :] += term

        return beta


    def solve(self, **kwargs) -> Iterable[np.ndarray, np.ndarray]:
        """
        Parameters
        ----------
            **kwargs : dict of np.ndarray
                Input arrays. Keys: 'soc', 'T'
        """

        mu_soc, sigma_soc = kwargs['soc']
        temp = kwargs['T']

        _ix = np.argmin(np.abs(self.temperatures-temp.mean())) # pick closest from availabe temperatures
        T = self.temperatures[_ix]
        coeff_c = np.flip(self.functions[T][0].coeffs)
        coeffs_d = np.flip(self.functions[T][1].coeffs)
        alpha = self.alpha*coeff_c + (1-self.alpha)*coeffs_d

        beta = self.get_beta(alpha, mu_soc, sigma_soc)
        gamma = self.A_inv@beta

        return gamma[0, :], np.sqrt(np.sum(gamma[1:, :]**2, axis= 0))

