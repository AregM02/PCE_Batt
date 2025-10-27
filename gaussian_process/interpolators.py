import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from pathlib import Path
from typing import Dict, List
import logging
import textwrap

from utils import load_vars # custom dataloader

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)

class Interpolator:
    """Handles 2D interpolation across SOC and temperature for a single parameter."""

    def __init__(self, soc_points: np.ndarray, temps:np.ndarray, data:List[np.ndarray]):
        """
        Args:
            soc_points: Array of SOC points for each temperature
            temps: Array of temperature values
            data: Parameter arrays at each (SOC, temp) combination
        """
        # dictionary containing values at all SoC and Temperature combinations
        self.raw_data = {temps[i]: (soc_points[i], data[i]) for i in range(len(data))}

        # Create interpolators for each temperature
        self.temp_interpolators = {
            temp: interp1d(soc_points, values,
                          kind='linear',
                          bounds_error=False,
                          fill_value=(values[0], values[-1]))
            for temp, (soc_points, values) in self.raw_data.items()
        }

        # Extract sorted temperature points
        self.temp_points = np.array(sorted(self.raw_data.keys()))

    def __call__(self, T:float, SoC:float)->float:
        """Interpolate parameter value at given temperature and SOC."""

        # Clip temperature to available range
        T_clipped = np.clip(T, self.temp_points.min(), self.temp_points.max())

        # Get values at this SoC for all temperatures
        val_at_SoC = []
        for temp in self.temp_points:
            val = self.temp_interpolators[temp](SoC)
            val_at_SoC.append(val)

        # Handle single temperature case
        if len(self.temp_points) == 1:
            return float(val_at_SoC[0])

        # Interpolate between temperatures
        return float(np.interp(T_clipped, self.temp_points, val_at_SoC))


class BatteryParameterInterpolator:
    """Handles loading and interpolation of battery parameter distributions across SOC and temperature."""
    
    def __init__(self, temps=(15, 25, 35, 45)):
        self.temps = temps # default test temperatures; feel free to change, but make sure to rename files accordingly
        self.data = {} # {temperature: {parameter: value array}}
        self.soc_points = {} # soc arrays at each temp
        self.interpolators = None # interpolator object for each parameter

        logger.info('[BatteryParameterInterpolator] Setup in progress...')
        self._set_data()
        self._setup_interpolators()
        logger.info('[BatteryParameterInterpolator] Initialized!')

    
    def _load_distribution(self, temp: int) -> Dict[str, np.ndarray]:
        """
        Loads parameter distributions for a specific temperature at all SoCs. Expects data to be stored under <data/distributions/[temperature].parquet> 
        
        """

        parquet_path = Path(__file__).parent.parent / "data" / "distributions" / f"{temp}deg.parquet"
        try:
            dist = pd.read_parquet(parquet_path)
        except:
            logger.critical(f'[BatteryParameterInterpolator] Expected to find "{temp}deg.parquet" at {parquet_path} but failed. Check if data is missing or set up improperly.')
            quit()
        
        return load_vars(dist)
    
    def _set_data(self) -> None:
        """Sets 'data' and 'soc_points' attributes."""
        
        for temp in self.temps:
            params = self._load_distribution(temp)
            self.soc_points[temp] = params.pop('SOC')
            self.data[temp] = params
        logger.info('[BatteryParameterInterpolator] Parameter distributions set!')

    def _setup_interpolators(self) -> None:
        """Sets up vecrorized interpolator objects for all parameters."""

        # Get common parameter names (excluding SOC and T)
        param_names = [k for k in self.data[self.temps[0]].keys() 
                       if not k.startswith(('SOC', 'T'))]
        
        # Prepare data for each parameter across temperatures
        param_groups = []
        for name in param_names:
            group = [self.data[temp][name] for temp in self.temps] # [array@15C, array@25C, ...]
            param_groups.append(group)
        
        # Create interpolators
        self.interpolators = [
            np.vectorize(
                Interpolator(
                [self.soc_points[temp] for temp in self.temps],
                [self.data[temp]['T'] for temp in self.temps], # use the measured temperature
                param_group)
                        )
            for param_group in param_groups
        ]

        param_names = textwrap.fill(str(param_names), width=80, subsequent_indent=" ")
        logger.info(f'[BatteryParameterInterpolator] Interpolators created for the following variables: \n{param_names}.')
    
    def get_interpolated_params(self, T:float, SoC:float)->Dict[str, float]:
        """Gets all interpolated parameters for given temperature and SOC.
        
        Returns:
            Dictionary of interpolated parameters {mu_r0: value, sigma_r0: value, ...}
        """
        param_names = [k for k in self.data[self.temps[0]].keys() 
                       if not k.startswith(('SOC', 'T'))]
        
        return {
            name: interpolator(T, SoC)
            for name, interpolator in zip(param_names, self.interpolators)
        }
    

    def test_sample(self) -> Dict[str, float]:
        """Draws a test sample."""
        T = np.random.choice(self.temps)
        SoC = np.random.uniform(low=0.05, high=0.95, size=(1,))
        return self.get_interpolated_params(T, SoC)
    

    def __repr__(self):
        keys = list(self.test_sample().keys())
        keys_str = textwrap.fill(str(keys), width=80, subsequent_indent=" ")
        return f"BATTERY PARAMETER INTERPOLATOR at <{hex(id(self))}>\nReturns:\n{keys_str}"
    

