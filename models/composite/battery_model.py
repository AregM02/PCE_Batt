import numpy as np
from collections.abc import Iterable
import yaml
import importlib
from pathlib import Path
from ..base import Model
from utils import create_logger

CONFIG_PATH = Path(__file__).parent.parent.parent / 'config' / 'models.yaml'

class BatteryModel(Model):
    """
    Composite battery model that combines multiple submodels.

    """

    def __init__(self, submodels: Iterable[Model] = None):
        super().__init__()

        self.logger = create_logger(__class__.__name__)

        # get the config (which models to compose)
        try:
            with open(CONFIG_PATH, "r") as f:
                config = yaml.safe_load(f)
        except FileNotFoundError:
            self.logger.fatal(f"Error while fetching config at {CONFIG_PATH}. Check if the file is available.")
            raise SystemExit
        
        # import required submodels and instantiate
        instances = []
        for m in config["submodels"]:
            mod = importlib.import_module(m["module"])
            cls = getattr(mod, m["class"])
            instances.append(cls(**m.get("kwargs", {})))
        
        if not instances: # empty list, possible error with config
            self.logger.fatal(f"No submodels found! Check for proper configuration at {CONFIG_PATH}")
            raise SystemExit
        self.submodels = instances


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

