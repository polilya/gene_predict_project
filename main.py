import logging

from Bio import SeqIO

import hydra
from hydra.core.config_store import ConfigStore, OmegaConf
from config import MLConfig

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import re

from tqdm import tqdm

import tools
from tools import calculate_kmer_features, signs

np.random.seed(42)

from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score

cs = ConfigStore.instance()
cs.store(name='ml_config', node=MLConfig)


@hydra.main(config_path='conf', config_name='config_start', version_base=None)
def main(cfg: MLConfig):

    logging.info('Script start')
    
    sequence_data = f'{cfg.paths.data}/{cfg.files.sequence}'
    with open(sequence_data) as file:
        fasta_sequences = SeqIO.parse(file, 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)

    data = pd.read_csv(f'{cfg.paths.data}/{cfg.files.features}')
    cds_location_data = pd.read_csv(f'{cfg.paths.data}/{cfg.files.cds_location}')
    cds_location_data = np.array(cds_location_data)

    Y = data['target'].values
    X = data.drop(['target'], axis=1).values

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.33, random_state=42, stratify=Y)

    model = SVC(**OmegaConf.to_container(cfg.params))

    X_train, X_test, LDA, x_scaler = tools.apply_lda(X_train, X_test, Y_train)

    logging.info('Start training...')
    model.fit(X_train, Y_train)
    logging.info('Training done...')
    prediction = model.predict(X_test)

    logging.info(f'f1 score: {f1_score(Y_test, prediction):.3f}')

    #loc_start_codon = [m.start() for m in re.finditer('ATG', sequence)]
    logging.info('Start evaluation...')

    step = cfg.add_params.step
    metrcis_dist = []
    plot_status = cfg.add_params.plot

    for graph in tqdm(range(0, cds_location_data.shape[0], 1)):

        i, j = cds_location_data[graph]
        length = j - i

        if plot_status:
            fig, ax = plt.subplots(dpi=120)
            ax.axvspan(i, j, alpha=0.6, color='red')

        scope = np.arange(i - 2*length, j + 2*length, step).astype('int')

        res = []
        
        for w in cfg.add_params.windows:
            for point in scope:
                scores = calculate_kmer_features(sequence[point - w // 2: point + w // 2], signs)['Entropy']
                scores = scores.values.reshape(1, -1)
                scores = x_scaler.transform(scores)
                scores = LDA.transform(scores)
                res.append(model.predict_proba(scores)[0][1])

        res = np.array(res).reshape(len(cfg.add_params.windows), -1).mean(axis=0)
        res_metric = tools.dice_metric(res, length, step)
        metrcis_dist.append(res_metric)

        # local_start_codon = [c for c in loc_start_codon if ((c >= i - lb) and (c <= j + rb))]
        # for codon in local_start_codon:
        #     ax.axvspan(codon, codon + 3, alpha=0.5, color='green')

        if plot_status:
            logging.info(f'prob_metric for sample {graph:3}: {res_metric:.3f}')
            plt.title(f'Example №{graph} [prob_metric: {res_metric:.3f}]')
            plt.plot(scope, res)
            plt.show()

    mu = np.array(metrcis_dist).mean()
    sigma = np.array(metrcis_dist).std()

    dist = np.random.normal(mu, sigma, 1000)
    density = stats.gaussian_kde(dist)

    n, x, _ = plt.hist(metrcis_dist, density=True, bins='sturges', linewidth=1, color='skyblue', edgecolor='black')
    x_d = np.linspace(x[0], x[-1], 500)
    plt.plot(x_d, density(x_d))
    plt.title(f'Dice metric distribution [μ={mu:.2f} std={sigma:.2f}]')
    plt.show()

    logging.info(f'Mean of Dice metric distribution: {mu:.2f}')
    logging.info(f'Std of Dice metric distribution: {sigma:.2f}')

    logging.info('Evaluation done')

    logging.info('Script done!')


if __name__ == "__main__":
    main()
