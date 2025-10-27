from pathlib import Path
import pandas as pd
import numpy as np
from collections import defaultdict
from abc import abstractmethod
from collections.abc import Iterable
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


class qOCV_POLY(Model):
    """
    Simulates a polynomial relationship between qOCV and SOC. Looks for qOCV data at <data/measurements/qocv20>
    
    """

    def __init__(self, order: int = 12, alpha: int = 0.5, logger_level: str = 'DEBUG'):
        super().__init__()
        
        self.logger = create_logger(__class__.__name__, logger_level)
        self.functions = defaultdict(list)
        self.temperatures = np.array([15, 25, 35, 45])
        self.alpha = alpha # 0.5 = charge/discharge OCVs are averaged; 1 = only charge; 0 = only discharge

        for temp in self.temperatures:
            df_c = pd.read_parquet(PATH_QOCV/f'qocv_{temp}deg_charge.parquet')
            df_d = pd.read_parquet(PATH_QOCV/f'qocv_{temp}deg_discharge.parquet')
            
            poly_c = np.poly1d(np.polyfit(df_c.SOC, df_c.Spannung, order))
            poly_d = np.poly1d(np.polyfit(df_d.SOC, df_d.Spannung, order))
            
            self.functions[temp] = [poly_c, poly_d]

        self.logger.info("[qOCV_POLY] Initialized!")


    def solve(self, **kwargs):
        """
        Parameters
        ----------
            **kwargs : dict of np.ndarray
                Input arrays. Keys: 'soc', 'T'.
        """
        soc = kwargs['soc']
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
    
