import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import random


FILE_PATH = Path(__file__).parent.parent / 'data' / 'measurements' / 'validation'

random.seed(1234)

def load_validation():
    fnames = [name for name in FILE_PATH.rglob('*.parquet')]
    
    dfs = []
    for name in fnames:
        df = pd.read_parquet(name)
        df = df.drop_duplicates(subset='time')
        df = df.reset_index(drop=True)

        # Cnom
        df_c = df.loc[df.status == 'CC DChg']
        jump_ix = df_c.index[df_c.index.diff()>1].values[0]
        df_c = df_c.loc[:jump_ix-5, :]
        cnom = -((df_c.current*df_c.time.diff()).cumsum()).iloc[-1]/3600.
        df = df.assign(cnom = cnom*np.ones_like(df.time))

        # SoC
        soc_100_ix = df.index[df.status=='CCCV Chg'].values[-1]
        soc = (df.current*df.time.diff()).cumsum()/(3600*cnom)
        soc.loc[0] = soc.loc[1] 
        df = df.assign(soc = 1.-soc.loc[soc_100_ix] + soc)

        # plt.plot(df.time - df.time.iloc[0], df.voltage)
        # break

        df = df.loc[df.status == 'Sim']
        df.time = df.time - df.time.iloc[0]
        df = df[['time', 'voltage', 'current', 'soc', 'cnom']]
        dfs.append(df)

        # plt.plot(df.time - df.time.iloc[0], df.voltage, c='b', alpha=0.04)
    
    # plt.show()
    # quit()


    which_df = np.argmin([len(frame) for frame in dfs])
    df = dfs[which_df]
    min_len = len(df)

    measurements = np.asarray([frame.voltage.to_numpy()[:min_len] for frame in dfs])
    n = len(measurements)
    measured_mean = measurements.mean(axis = 0)
    measured_std = measurements.std(axis = 0)*n/(n-1) # sample correction

    T = 25.*np.ones_like(df.time)
    cap = df.cnom.iloc[0]
    time = df.time.to_numpy()
    current = df.current.to_numpy()
    soc = df.soc.to_numpy()

    measurements = [measured_mean, measured_std]

    return (time, current, measurements, soc, T, cap)
