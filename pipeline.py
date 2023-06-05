import data_preprocessing
import model
import model_selection

import hydra
from hydra.core.config_store import ConfigStore
from hydra.core.hydra_config import HydraConfig
from config import MLConfig

cs = ConfigStore.instance()
cs.store(name='ml_config', node=MLConfig)


@hydra.main(config_path='conf', config_name='config_pipeline_start', version_base=None)
def main(cfg: MLConfig):
    data_preprocessing.main(cfg)
    model_selection.main(cfg)
    model.main(cfg)


if __name__ == "__main__":
    main()