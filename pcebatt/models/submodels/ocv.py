import math
from pathlib import Path
import pandas as pd
import numpy as np
import chaospy
from numpy.typing import NDArray
from abc import abstractmethod
from collections import defaultdict
from collections.abc import Iterable
from scipy.interpolate import interp1d
from ...utils import create_logger, default_config
from ..base import Model
from types import SimpleNamespace

class OCV(Model):
    """
    Template for an OCV object. Must have a solve method, which must return the simulated voltage array and the standard deviation of its uncertainty. If the model is deterministic, return zeros for the std. 
    
    """

    def __init__(self, config = None):
        super().__init__(config = config)
        self.model_type = 'ocv'
        self.voltage = SimpleNamespace()

        self.PATH_QOCV = Path(self.config['simulation']['paths']['ocv']['ocv_folder'])
        self.NAME_PATTERN = {
        'cha': self.config['simulation']['paths']['ocv']['file_naming_scheme']['charge'], 
        'dch': self.config['simulation']['paths']['ocv']['file_naming_scheme']['discharge']}

    @abstractmethod
    def solve(self, *args, **kwargs) -> Iterable[NDArray, NDArray]:
        pass


class qOCV_POLY(OCV):
    """
    Simulates a polynomial relationship between qOCV and SOC. Looks for qOCV data at <data/measurements/qocv20>
    
    """

    def __init__(self, order: int = 6, alpha: int = 0.5, logger_level: str = 'DEBUG',
                 logger_name=None, config = None) -> None:
        super().__init__(config=config)
        
        if logger_name is None:
            self.logger = create_logger(__class__.__name__, logger_level)
        else: 
            self.logger = create_logger(logger_name, logger_level)

        self.functions = defaultdict(list)
        self.temperatures = np.array([15, 25, 35, 45])
        self.alpha = alpha # 0.5 = charge/discharge OCVs are averaged; 1 = only charge; 0 = only discharge

        for temp in self.temperatures:
            df_c = pd.read_parquet(self.PATH_QOCV/ self.NAME_PATTERN['cha'].format(temp=temp))
            df_d = pd.read_parquet(self.PATH_QOCV/ self.NAME_PATTERN['dch'].format(temp=temp))
            
            poly_c = np.poly1d(np.polyfit(df_c.SOC, df_c.Spannung, order))
            poly_d = np.poly1d(np.polyfit(df_d.SOC, df_d.Spannung, order))
            
            self.functions[temp] = [poly_c, poly_d]

        self.logger.info(f"[{self.logger.name}] Initialized!")


    def solve(self, **kwargs) -> Iterable[NDArray, NDArray]:
        """
        Parameters
        ----------
            **kwargs : dict of NDArray
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
  
        self.voltage.mean, self.voltage.std = ocv_interp, np.zeros_like(ocv_interp)
        return self.voltage.mean, self.voltage.std


class qOCV_PCE(qOCV_POLY):
    """
    A polynomial chaos expanded variant of qOCV_POLY. Looks for qOCV data at <data/measurements/qocv20>
    
    """

    def __init__(self, order: int = 6, alpha: int = 0.5, logger_level: str = 'DEBUG', config=None) -> None:
        super().__init__(order=order, alpha=alpha, logger_level=logger_level,
                         logger_name=__class__.__name__, config=config)
        
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
        

    def get_beta(self, alpha:NDArray, mu_t: NDArray, sigma_t:NDArray) -> NDArray:
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


    def solve(self, **kwargs) -> Iterable[NDArray, NDArray]:
        """
        Parameters
        ----------
            **kwargs : dict of NDArray
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

        self.voltage.mean, self.voltage.std = gamma[0, :], np.sqrt(np.sum(gamma[1:, :]**2, axis= 0))
        return self.voltage.mean, self.voltage.std


class qOCV_Differential(OCV):
    def __init__(self, alpha=0.5, logger_level='DEBUG', logger_name=None, config=None):
        super().__init__(config=config)
        self.alpha = alpha
        self.temperatures = np.array([15, 25, 35, 45])
        self.logger = create_logger(logger_name or __class__.__name__, logger_level)

        # PRE-LOAD DATA AND CREATE INTERPOLATORS
        self.interpolators = {}
        for T in self.temperatures:
            df_c = pd.read_parquet(self.PATH_QOCV / self.NAME_PATTERN['cha'].format(temp=T))
            df_d = pd.read_parquet(self.PATH_QOCV / self.NAME_PATTERN['dch'].format(temp=T))

            # Numerical derivatives
            dvdz_c = df_c.Spannung.diff(periods=100) / df_c.SOC.diff(periods=100)
            dvdz_d = df_d.Spannung.diff(periods=100) / df_d.SOC.diff(periods=100)

            valid_c = ~(dvdz_c.isna() | df_c.Spannung.isna())
            valid_d = ~(dvdz_d.isna() | df_d.Spannung.isna())

            self.interpolators[T] = {
                'v_c': interp1d(df_c.SOC.values[valid_c], df_c.Spannung.values[valid_c], kind='linear', fill_value="extrapolate"),
                'v_d': interp1d(df_d.SOC.values[valid_d], df_d.Spannung.values[valid_d], kind='linear', fill_value="extrapolate"),
                'dvdz_c': interp1d(df_c.SOC.values[valid_c], dvdz_c.values[valid_c], kind='linear', fill_value="extrapolate"),
                'dvdz_d': interp1d(df_d.SOC.values[valid_d], dvdz_d.values[valid_d], kind='linear', fill_value="extrapolate")
            }
        self.logger.info(f"[{self.logger.name}] Initialized with pre-loaded interpolators!")


    def solve(self, **kwargs):
        mu_soc, sigma_soc = kwargs['soc']
        T_val = kwargs['T']
        
        # Pick closest temperature
        T_nearest = self.temperatures[np.argmin(np.abs(self.temperatures - np.mean(T_val)))]
        ips = self.interpolators[T_nearest]

        # Quick lookups in RAM
        v_mean = self.alpha * ips['v_c'](mu_soc) + (1 - self.alpha) * ips['v_d'](mu_soc)
        dvdz_mean = self.alpha * ips['dvdz_c'](mu_soc) + (1 - self.alpha) * ips['dvdz_d'](mu_soc)
        
        self.voltage.mean, self.voltage.std = v_mean, np.abs(dvdz_mean) * sigma_soc
        return self.voltage.mean, self.voltage.std
