from abc import ABC, abstractmethod
from ..utils import default_config

class Model(ABC):
    """Base template of a Model class. Must implement a solve method!"""

    def __init__(self, config = None):
        if config is None:
            dict_cfg, path_cfg = default_config()
            self.config = dict_cfg | {"abs_config_path": path_cfg}
        else:
            self.config = config
        pass

    @abstractmethod
    def solve(self, *args, **kwargs):
        pass
    