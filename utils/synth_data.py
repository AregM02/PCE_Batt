# SOME USEFUL FUNCTIONS FOR CREATING TEST PROFILES
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def make_data(short=False, plot=False):
    """Creates a pulse profile"""
    
    if short:
        dt = np.ones(1000)
        time = np.arange(dt.shape[0])
        ib = np.zeros(dt.shape)
        ib[200:300] = -1.
    else:
        dt = np.ones(30000)
        time = np.arange(dt.shape[0])
        ib = np.zeros(dt.shape)
        dist = int(len(dt) / 10)
        for n in range(10):
            ib[dist * n + 200:dist * n + 300] = 1.
            ib[dist * n + 600:dist * n + 700] = -1.
            ib[dist * n + 2000:dist * n + 2540] = -2.

    SOC = 1 + np.cumsum(dt * ib / (3600 * 3))
    ocv, vb, T = np.zeros_like(ib), np.zeros_like(ib), 25 * np.ones_like(ib)

    if plot:        
        plt.plot(ib)
        plt.plot(SOC)
        plt.show()

    return dt, time, ib, vb, SOC, ocv, T, 


make_data()