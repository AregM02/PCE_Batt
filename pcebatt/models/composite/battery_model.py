import pprint
import numpy as np
import yaml
import importlib
from pathlib import Path
from ...utils import create_logger, default_config


class BatteryModel():
    """
    Composite battery model that combines multiple submodels.

    """

    def __init__(self, config = None, load_default = False):
        self.logger = create_logger(__class__.__name__)

        if config is None: # user has not specified their own config
            if not load_default: # user wants to configure the model in a script, keep it empty
                self.submodels = []
                return
            else: # load the default config
                self.config_dict, config = default_config()
        else:
            config = Path(config).resolve()
            try:
                with open(config, "r") as f:
                    self.config_dict = yaml.safe_load(f)

            except FileNotFoundError:
                self.logger.fatal(f"Error while fetching config at {config}. Check if the file is available.")
                raise SystemExit
        

        # import required submodels and instantiate
        instances = []
        for m in self.config_dict["models"]:
            mod = importlib.import_module(m["module"])
            cls = getattr(mod, m["class"])
            instances.append(cls(**m.get("kwargs", {}),
                                config=self.config_dict | {"abs_config_path": config})) # append the current config dict with its absolute path
        
        if not instances: # empty list, possible error with config
            self.logger.fatal(f"No submodels found! Check for proper configuration at {config}")
            raise SystemExit
        self.submodels = instances
        self.config_path = config


    def solve(self, current, time, soc, T):
        inputs = dict(current=current, time=time, soc=soc, T=T)
        combined_mean = None
        combined_var = None

        for model in self.submodels:
            mean, std = model.solve(**inputs)

            if combined_mean is None:
                combined_mean = mean.copy()
                combined_var = std**2
            else:
                combined_mean += mean
                combined_var += std**2

        combined_std = np.sqrt(combined_var)

        return combined_mean, combined_std
    
    def __str__(self):
        return pprint.pformat(self.get_submodels())
    
    def get_submodels(self):
        return {k.model_type : k for k in self.submodels}
    
    def get_current_config(self):
        return self.config_path, self.config_dict
    