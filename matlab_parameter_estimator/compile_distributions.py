from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation, norm
import matplotlib.pyplot as plt
from collections.abc import Iterable
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s:%(message)s')
# handler = logging.StreamHandler()
handler = logging.FileHandler((Path(__file__).parent/'outputs'/'info.log').resolve())
handler.setFormatter(formatter)
logger.addHandler(handler)


def filter_MAD(arr: np.array, ret_index=False, thresh=5)->Iterable[np.array, np.array]:
    """Returns: 
            arr: MAD filtered array
            indices: original indices of the removed elements(optional) 
    """

    med = np.median(arr)
    mad = median_abs_deviation(arr)
    new_arr =  arr[np.abs(arr - med) < thresh * mad]
    indices = np.where(np.abs(arr - med) > thresh * mad)[0]
    
    if ret_index:
        return new_arr, indices
    else:
        return new_arr


def sync_filter(df: pd.DataFrame, cols: Iterable[str])->None:
    """Applies synced MAD filtering along passed columns. 
       Useful if you need to caculate correlations, as separate column-by-column
       filtering would mess with the array order."""

    for row in df.index:
        idx_total = np.empty(0, dtype=int)
        for col in cols:
            _, ix = filter_MAD(df.at[row, col], ret_index=True)
            idx_total = np.append(idx_total, ix)
        for col in cols:
            df.at[row, col] = np.delete(df.at[row, col], idx_total)


def fit_dist(arr: np.array, plot=False):

    p = norm.fit(arr)
    if plot:
        print(p)
        x = np.linspace(0.9 * np.min(arr), np.max(arr) * 1.1, 1000)
        plt.plot(x, norm.pdf(x, loc=p[0], scale=p[1]), c='black')
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.hist(arr, bins=15, color='b', alpha=0.5, density=True, label='Data')
        plt.legend()
        # plt.savefig('dist_example.png', dpi=200)
        plt.show()

    return p


def snap_to_grid(df:pd.DataFrame, column_name:str, lower=0, upper=101, step=5, radius=5)->pd.DataFrame:
    """Snap a DataFrame column to nearest grid values"""

    grid = np.arange(lower, upper, step)
    values = df[column_name].values
    values = values.astype(int)
    
    # Find closest grid point for each value
    differences = np.abs(values[:, None] - grid)
    closest_indices = np.argmin(differences, axis=1)
    snapped_values = grid[closest_indices]
    
    # Apply radius constraint
    within_radius = np.abs(values - snapped_values) <= radius
    snapped_values = np.where(within_radius, snapped_values, values)
    snapped_values = snapped_values.astype(int)
    
    return df.assign(**{column_name: snapped_values})


def compile_data(temp:int)->pd.DataFrame:
    "Returns a combined dataframe, where each entry holds a list of parameter values from all batteries"

    if temp not in [15, 25, 35, 45]:
        raise ValueError("Temperature must be one of: 15, 25, 35, 45")
    
    param_path = (Path(__file__).parent.parent / "data"/ "parametrization" / f"{temp}deg").resolve() # full path
    try: 
        pd.read_parquet(next(param_path.rglob('*.parquet')))
    except:
        raise FileNotFoundError(f"Expected parametrization data at {param_path} but found nothing. Please check if data is set up properly.")
    
    dfs = [] # list of dataframes to be combined

    for fname in param_path.rglob('*.parquet'):
        df = pd.read_parquet(fname)

        # specific corrections to our data, change if yours doesn't have these issues
        if len(df) < 20: # separate treatment of initial measurements (10 SOCs, 300)
            if np.linalg.norm(np.mean(df['T'])-temp) < 1: # filter out wrong temperatures
                df = df.drop(0)
            else:
                continue
        if len(df) > 40: # the data in some of the files has been duplicated
            df = df.loc[0:20, :]

        df = snap_to_grid(df, 'SOC') # synchronize SOC values
        df = df.set_index('SOC') # set SOC as index
        df = df[~df.index.duplicated(keep='first')]
        # full_index = np.arange(100,-1, -5)
        # df = df.reindex(full_index) # index padding (unnecessary)
        df = df.sort_index()
        dfs.append(df)

    # group all dataframes into one
    df_grouped = pd.concat(dfs).groupby(level=0).agg(list)
    df_grouped = df_grouped.map(np.array) 

    return df_grouped


def main():
    for temp in [15, 25, 35, 45]:
        # load data
        data = compile_data(temp)

        data = data[['R0', 'R1', 'tau1', 'R2', 'tau2', 'T']]

        # change of variables
        data[['tau1', 'tau2']] = data[['tau1', 'tau2']].apply(lambda x: 1/x)
        data['R1'], data['R2'] = data['R1'] * data['tau1'], data['R2'] * data['tau2']
        data = data.rename(columns={'R1':'c1_inv', 'tau1':'tau1_inv',
                                    'R2':'c2_inv', 'tau2':'tau2_inv'})
        
        # clean outliers
        data['R0'] = data['R0'].apply(filter_MAD)
        sync_filter(data, ['c1_inv', 'tau1_inv'])
        sync_filter(data, ['c2_inv', 'tau2_inv'])
        
        # fit distributions
        data_fit = data.map(fit_dist)

        # split columns into mean and sigma
        data_mean = data_fit.map(lambda x: x[0])
        data_sigma = data_fit.map(lambda x: x[1])
        data_mean.columns = [f"mu_{col}" for col in data_fit.columns]
        data_sigma.columns = [f"sigma_{col}" for col in data_fit.columns]
        out = pd.concat([data_mean, data_sigma], axis=1)

        # calculate and add correlations to the output file
        rho1 = data.apply(lambda row: np.corrcoef(row['tau1_inv'], row['c1_inv'])[0, 1], axis=1)
        rho2 = data.apply(lambda row: np.corrcoef(row['tau2_inv'], row['c2_inv'])[0, 1], axis=1)
        out['rho1'], out['rho2'] = rho1, rho2

        # dump to file
        save_path = Path(__file__).parent/'outputs'
        file_name = f'{temp}deg.parquet'
        out.to_parquet((save_path/file_name).resolve())

        # log the creation date and content
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            logger.info(f' {file_name} created: contains \n {out} \n\n')

if __name__=="__main__":
    main()