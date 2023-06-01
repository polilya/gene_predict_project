import logging

import hydra
from hydra.core.config_store import ConfigStore
from hydra.core.hydra_config import HydraConfig
from config import MLConfig

from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression

from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
import scipy.stats as stats

import numpy as np
import pandas as pd
import tools

np.random.seed(42)

cs = ConfigStore.instance()
cs.store(name='ml_config', node=MLConfig)


@hydra.main(config_path='conf', config_name='config_start', version_base=None)
def main(cfg: MLConfig):

    logging.info('Script starts')
    hydra_cfg = HydraConfig.get()
    logging.info(f'Used config: {hydra_cfg.job.config_name}')

    data = pd.read_csv(f'{cfg.paths.data}/{cfg.files.features}')

    Y = data['target'].values
    X = data.drop(['target'], axis=1).values

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.33, random_state=42, stratify=Y)
    logging.info(f'X_train shape: {(X_train.shape[0], X_train.shape[1])}')

    X_train, X_test, dim_reduction, x_scaler = tools.apply_lda(X_train, X_test, Y_train)
    logging.info(f'X_train shape after {dim_reduction}: {(X_train.shape[0], X_train.shape[1])}')

    model = SVC()
    logging.info(f'Model: {model}')

    distributions = dict(C=stats.uniform(loc=1, scale=300),
                         kernel=['linear', 'rbf', 'sigmoid'],
                         gamma=['scale', 'auto'],
                         tol=stats.uniform(loc=0, scale=1e-1),
                         decision_function_shape=['ovo', 'ovr'],
                         class_weight=[{0: 1, 1: 1}, {0: 1, 1: 2},
                                       {0: 1, 1: 3}, {0: 2, 1: 1}])

    clf = RandomizedSearchCV(model, distributions, scoring='f1', random_state=0, verbose=2, n_iter=20, n_jobs=2)
    search = clf.fit(X_train, Y_train)
    model_params = search.best_params_
    model_params['probability'] = True

    logging.info(f'Model parameters: {model_params}')
    logging.info(f'Model best score: {search.best_score_:.3f}')
    logging.info('Script done!')


if __name__ == "__main__":
    main()
