import yaml
from ..utils import get_project_root

DEFAULT_CONFIG_PATH = get_project_root() / 'config' / 'cfg.yaml'

def default_config():
    with open(DEFAULT_CONFIG_PATH, "r") as f:
        config_dict = yaml.safe_load(f)
        paths = config_dict['simulation']['paths']
        pre_path = DEFAULT_CONFIG_PATH.parent
        paths['parameters']['dist_folder'] = (pre_path / (paths['parameters']['dist_folder'])).resolve()
        paths['ocv']['ocv_folder'] = (pre_path / paths['ocv']['ocv_folder']).resolve()
        paths['validation']['validation_path'] = (pre_path / paths['validation']['validation_path']).resolve()
    return config_dict, DEFAULT_CONFIG_PATH