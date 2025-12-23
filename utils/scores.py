# Some numerical scores for assessing validation performace

import numpy as np


def RMSE(arr1, arr2):
    return np.sqrt(np.mean((arr1 - arr2)**2))