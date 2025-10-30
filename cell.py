from models.composite import BatteryModel
from soc_estimation import CoulombCounterGalerkin
from types import SimpleNamespace
import numpy as np


class Cell:
    def __init__(self, capacity: float, initial_soc: float, 
                 initial_soc_unc: float = 1e-3, capacity_unc: float = 1e-3):
        
        self.model = BatteryModel()
        self.soc_counter = CoulombCounterGalerkin(initial_soc=initial_soc, capacity=capacity,
                                                  initial_soc_unc=initial_soc_unc, capacity_unc=capacity_unc)

        # save profiles in case the user wants to access them later (mean and std)
        self.voltage = SimpleNamespace()
        self.soc = SimpleNamespace()
        self.temperature = SimpleNamespace()


    def solve(self, current, time, soc=None, temperature=None):
        """
        Simulates the cell and sets the results as attributes. 
        Parameters <soc> and <temperature>, if passed, will override calculations.
        """

        if soc is None:
            soc_mean, soc_sigma = self.soc_counter(current, time)
        else:
            soc_mean, soc_sigma = soc, np.zeros_like(soc)

        if temperature is None:
            ## Thermal model here
            # temperature_mean, temperature_sigma = ThermalModel(...)
            pass
        else: 
            temperature_mean, temperature_sigma = temperature, np.zeros_like(temperature)

        voltage_mean, voltage_sigma = self.model.solve(current, time, soc_mean, temperature)

        self.voltage.mean, self.voltage.std = voltage_mean, voltage_sigma
        self.temperature.mean, self.temperature.std = temperature_mean, temperature_sigma
        self.soc.mean, self.soc.std = soc_mean, soc_sigma
