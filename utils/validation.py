import pandas as pd
import numpy as np
from pathlib import Path


FILE_PATH = Path(__file__).parent.parent / 'data' / 'measurements' / 'validation'


def load_validation():
    fnames = [name for name in FILE_PATH.rglob('*.parquet')]
    name = np.random.choice(fnames)
    df = pd.read_parquet(name)

    df_kap = df[(df.Prozedur == 'jri_KapTest_C2') & (df.Zustand == 'CHA')]
    cap = (df_kap.Strom*df_kap.Zeit.diff()).cumsum().iloc[-1]/3600 # determine the exact capacity

    df = df.drop_duplicates(subset='Programmdauer')
    df.Programmdauer = df.Programmdauer/1000
    dt = df.Programmdauer.diff()

    ref = df[df.Prozedur == 'jri_KapTest_C2'].index[0]
    soc = np.cumsum(dt * df.Strom / (3600 * cap))
    soc = soc + 1. - soc.loc[ref]
    df['soc'] = soc

    start = np.where(df.Prozedur == 'jri_Charge_C2')[0][-1] + 1
    start = df[df.Strom.abs() > 3.5].index[0] - 100
    end = df[df.Strom.abs() > 3.5].index[-1] + 1500

    df = df.loc[start:end]
    time, current, voltage, soc, T = df[['Programmdauer', 'Strom', 'Spannung', 'soc', 'Temp']].values.T

    return (time-time[0], current, voltage, soc, T, cap)
