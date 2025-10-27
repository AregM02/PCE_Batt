import pandas as pd
import numpy as np
from pathlib import Path


def load_vars(dist: pd.DataFrame)->dict['str', np.array]:
    """
    Parses data required for gaussian process interpolators. You are free to define your own variable names here, but keep in mind the specific naming scheme requirements for each model. See README.md for details.

    """
    
    return {
            'mu_r0': dist['mu_R0'].to_numpy(),
            'sigma_r0': dist['sigma_R0'].to_numpy(),
            'mu_tau1_inv': dist['mu_tau1_inv'].to_numpy(),
            'sigma_tau1_inv': dist['sigma_tau1_inv'].to_numpy(),
            'mu_c1_inv': dist['mu_c1_inv'].to_numpy(),
            'sigma_c1_inv': dist['sigma_c1_inv'].to_numpy(),
            'mu_tau2_inv': dist['mu_tau2_inv'].to_numpy(),
            'sigma_tau2_inv': dist['sigma_tau2_inv'].to_numpy(),
            'mu_c2_inv': dist['mu_c2_inv'].to_numpy(),
            'sigma_c2_inv': dist['sigma_c2_inv'].to_numpy(),
            'SOC': dist.index.to_numpy() / 100., # mandatory
            'T': int(round(dist['mu_T'].mean())), # mandatory
            'rho1': dist['rho1'].to_numpy(),
            'rho2': dist['rho2'].to_numpy(),
            }


# tesing
if __name__ == "__main__":
    TEMP = 15
    path = Path(__file__).parent.parent/'data'/'distributions'/f'{TEMP}deg.parquet'
    dist = pd.read_parquet(path.resolve())
    print(load_vars(dist).keys())