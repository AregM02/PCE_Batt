import pytest
import numpy as np
from pcebatt.utils import load_validation, RMSE
from pcebatt.core import Cell


def test_cell_simulation_accuracy():
    time, current, measurement, soc, temperature, C_nom, C_nom_unc = load_validation(short=True)
    cell_params = dict(
        initial_soc=soc[0], 
        capacity=C_nom, 
        capacity_unc=C_nom * C_nom_unc
    )

    cell = Cell(**cell_params)
    cell.solve(current=current, time=time, temperature=temperature)

    # Calculate RMSE in mV
    mean_rmse = 1000 * RMSE(cell.voltage.mean, measurement[0])
    std_rmse = 1000 * RMSE(cell.voltage.std, measurement[1])

    assert mean_rmse < 50, f"Mean Voltage RMSE too high: {mean_rmse:.2f} mV"
    assert not np.isnan(cell.voltage.mean).any(), "Simulation resulted in NaN values"
    
    print(f"\nTest passed with Mean RMSE: {mean_rmse:.2f} mV")