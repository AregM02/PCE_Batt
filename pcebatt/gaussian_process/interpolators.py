from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import logging
from typing import Dict
import matplotlib.pyplot as plt
import yaml
from ..utils import create_logger


# ---------- MultiIndex BatteryParameterInterpolator ----------
class BatteryParameterInterpolator:
    """
    Interpolates scalar parameters and correlation matrices over (Temperature, SOC) 
    using a unified MultiIndex grid.

    Args:
        temps: List of temperature values corresponding to available data files.
        data_dir: Directory containing the parquet files with parameter distributions.

    """

    def __init__(self, config=None):

        DIST_FOLDER = Path(config['simulation']['paths']['parameters']['dist_folder']).resolve()
        NAME_PATTERN = config['simulation']['paths']['parameters']['file_naming_scheme']
        TEMPS = list(config['simulation']['temperatures'])
        PARAM_NAMES = list(config['simulation']['paths']['parameters']['param_labels'])
        LOG_LEVEL = config['defaults']['logger_level']
        
        self.temps = TEMPS
        self.data_dir = DIST_FOLDER
        self.pattern = NAME_PATTERN

        # Define scalar parameter names (change as needed)
        self.scalar_names = PARAM_NAMES
        self.scalar_interpolator = None
        self.corr_interpolator = None
        self.D = None

        self.logger = create_logger(__class__.__name__, level = LOG_LEVEL)

        self.logger.info(f"[{self.__class__.__name__}] Setup in progress...")
        self._load_and_prepare_data()
        self._setup_interpolators()
        self.logger.info(f"[{self.__class__.__name__}] Initialized!")


    def _load_and_prepare_data(self):
        """Load parquet files, merge into a MultiIndex DataFrame, flatten CORR matrices."""
        dfs = []
        for temp in self.temps:
            file_path = self.data_dir / self.pattern.format(temp=temp)
            if not file_path.exists():
                raise FileNotFoundError(f"{file_path} not found.")
            df = pd.read_parquet(file_path)
            # Add temperature column
            df['T'] = temp
            dfs.append(df)

        # Merge into one DataFrame with MultiIndex (T, SOC)
        self.df_all = pd.concat(dfs)
        self.df_all.index = pd.MultiIndex.from_arrays([self.df_all['T'], self.df_all.index/100.], names=['T', 'SOC'])
        self.df_all.drop(columns='T', inplace=True)
        self.df_all.sort_index(inplace=True)

        # Flatten CORR matrices for interpolation
        if isinstance(self.df_all['CORR'].iloc[0], list):
            self.D = int(len(self.df_all['CORR'].iloc[0])**0.5)
        else: 
            self.D = int(self.df_all['CORR'].iloc[0].shape[0] ** 0.5)
        
        # Prepare scalar array and CORR array for grid
        self.scalar_array = self.df_all[self.scalar_names].to_numpy() # shape: (num_points, num_scalars)
        self.corr_array = np.stack(self.df_all['CORR'].to_numpy(), axis=0) # shape: (num_points, D*D)
        
        # Extract unique grid values
        self.T_unique = np.sort(self.df_all.index.get_level_values('T').unique())
        self.SOC_unique = np.sort(self.df_all.index.get_level_values('SOC').unique())

        # Reshape into 2D grid: (T, SOC, ...)
        self.scalar_grid = self.scalar_array.reshape(len(self.T_unique), len(self.SOC_unique), -1)
        self.corr_grid = self.corr_array.reshape(len(self.T_unique), len(self.SOC_unique), self.D*self.D)

    def _setup_interpolators(self):
        """Create RegularGridInterpolator for scalars and CORR."""
        self.scalar_interpolator = RegularGridInterpolator(
            (self.T_unique, self.SOC_unique), self.scalar_grid, bounds_error=False, fill_value=None
        )
        self.corr_interpolator = RegularGridInterpolator(
            (self.T_unique, self.SOC_unique), self.corr_grid, bounds_error=False, fill_value=None
        )


    def get_interpolated_params(self, T, SoC)-> dict:
        """
        Returns interpolated parameters at given T and SoC.
        Supports:
        - Scalars: T, SoC can be floats or 1D arrays of same length
        - Correlations: returned as D x D (single) or N x D x D (vectorized)

        Returns:
            dict: keys are parameter names (including 'CORR'), values are np.arrays
        """
        # Ensure inputs are arrays for uniform processing
        T = np.atleast_1d(T)
        SoC = np.atleast_1d(SoC)
        if T.shape != SoC.shape:
            raise ValueError("T and SoC must have the same shape")

        # Combine into query points
        query_points = np.column_stack([T, SoC])  # shape (N, 2)

        # Scalars: shape (N, num_scalars)
        scalar_vals = self.scalar_interpolator(query_points)  # (N, num_scalars)

        # Correlations: shape (N, D*D)
        corr_flat = self.corr_interpolator(query_points)
        corr_vals = corr_flat.reshape(-1, self.D, self.D)  # shape (N, D, D)

        # Build a dictionary of arrays
        result = {name: scalar_vals[:, j] for j, name in enumerate(self.scalar_names)}
        result['CORR'] = corr_vals

        # If scalar input (length 1), return 1D arrays or single matrix
        if len(T) == 1:
            result = {k: (v[0] if k != 'CORR' else v[0, :, :]) for k, v in result.items()}

        return result
    

    def get_deterministic_params(self, T=None, SoC=None) -> dict:
        """
        Returns a deterministic parameter dictionary. 
        If T and SoC are provided, the scalars are broadcasted to match their shape.
        
        - All 'mu_' parameters are averaged over the entire (T, SOC) range.
        - All 'sigma_' parameters are set to 0.0.
        - 'CORR' is returned as an identity matrix.
        
        Args:
            T (float or ndarray, optional): Temperature values to define output shape.
            SoC (float or ndarray, optional): SoC values to define output shape.

        Returns:
            dict: Deterministic parameter set matching the input shape.
        """

        # find output length
        if T is not None and SoC is not None:
            T = np.atleast_1d(T)
            SoC = np.atleast_1d(SoC)
            if T.shape != SoC.shape:
                raise ValueError("T and SoC must have the same shape")
            out_len = len(T)
        else:
            out_len = 1

        # calculate global averages for 'mu_' and other scalars
        mu_cols = [c for c in self.df_all.columns if c.startswith('mu_')]
        sigma_cols = [c for c in self.df_all.columns if c.startswith('sigma_')]
        global_means = self.df_all[mu_cols].mean()

        result = {}

        for col in mu_cols:
            val = global_means[col]
            result[col] = np.full(out_len, val) if out_len > 1 else val

        for col in sigma_cols:
            result[col] = np.zeros(out_len) if out_len > 1 else 0.0

        if out_len > 1:
            # Return shape (N, D, D)
            eye_stack = np.eye(self.D)[np.newaxis, :, :] # (1, D, D)
            result['CORR'] = np.repeat(eye_stack, out_len, axis=0)
        else:
            # Return shape (D, D)
            result['CORR'] = np.eye(self.D)

        return result


    def test_sample(self):
        return self.get_interpolated_params(0, 0)


if __name__ == "__main__":
    temps = [15, 25, 35, 45]
    interpolator = BatteryParameterInterpolator(temps=temps)

    # Test interpolation
    param_name = 'sigma_R0'
    T_test = 25.0
    SoC_test = np.linspace(0, 1, 20)
    param = []
    for soc in SoC_test:
        params = interpolator.get_interpolated_params(T_test, soc)
        param.append(params[param_name])
        # print(f"T={T_test}, SoC={soc} => Params: {params}")

    plt.plot(SoC_test, param, marker='o')
    plt.show()