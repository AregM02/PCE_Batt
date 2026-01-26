"""
Some functions for caching results and combining with data. 

"""
import pickle
import numpy as np
from src.utils.paths import get_project_root

SAVE_LOC = get_project_root() / 'data' / 'cache'

def save_cache(data, name):
    with open(SAVE_LOC / f'{name}.pkl', 'wb') as file:
        pickle.dump(data, file)

def load_cached(name):
    with open(SAVE_LOC / f'{name}.pkl', 'rb') as file:
        x = pickle.load(file)
    return x

def unify_measurements(mean1, std1, mean2, std2):
    """
    Combines means and standard deviations of two independent measurements.
    total_variance = sum of all variances
     -> Assumes stochastic independence!
    """
    return mean1+mean2, np.sqrt(std1**2 + std2**2)