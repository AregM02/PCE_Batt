from ..models.composite import BatteryModel
from ..models.base import Model
from ..soc_estimation import CoulombCounterGalerkin
from types import SimpleNamespace
import numpy as np


class Cell:
    def __init__(self, capacity: float, initial_soc: float, 
                 initial_soc_unc: float = 1e-3, capacity_unc: float = 1e-3,
                 config = None, load_default_config = True):
        
        if config is None: # user has not specified their own config
            if not load_default_config: # user wants to configure the model in a script, keep it empty
                print(f"[{self.__class__.__name__}] Initializing empty cell model.")
            else: # load the default config
                print(f"[{self.__class__.__name__}] No configuration path specified. Loading defaults...")
        else:
            print(f"[{self.__class__.__name__}] Loading configuration from {config}")
        
        self.model = BatteryModel(config = config, load_default=load_default_config)
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

        voltage_mean, voltage_sigma = self.model.solve(current, time, [soc_mean, soc_sigma], temperature)

        self.voltage.mean, self.voltage.std = voltage_mean, voltage_sigma
        self.temperature.mean, self.temperature.std = temperature_mean, temperature_sigma
        self.soc.mean, self.soc.std = soc_mean, soc_sigma


    def get_model_dict(self):
        return self.model.get_submodels()
    

    def __str__(self):
        return f"Battery Cell with Submodels:\n{self.model.__str__()}"
    

    def set_model(self, model_obj: "Model"):
        """
        Add a desired model object to the full model.
        Params:
            model (obj): a reference to the instantiated model object
        """
        type = model_obj.model_type # the object itself should provide the type
        self.model.submodels = [m for m in self.model.submodels if m.model_type != type]
        self.model.submodels.append(model_obj)
    
    def get_current_config(self):
        return self.model.get_current_config()