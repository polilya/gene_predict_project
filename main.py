import logging

from Bio import SeqIO

import hydra
from hydra.core.config_store import ConfigStore, OmegaConf
from config import MLConfig

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
    seq_before_cds_location = np.array(cds_location_data)

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

    lb = cfg.add_params.left_boundary
    rb = cfg.add_params.right_boundary
    step = cfg.add_params.step
    
    for graph in range(1, 251, 50):

        fig, ax = plt.subplots(dpi=120)
        i, j = seq_before_cds_location[graph]
        ax.axvspan(i, j, alpha=0.6, color='red')
        scope = np.arange(i - lb, j + rb, step).astype('int')

        res = []
        
        for w in cfg.add_params.windows:
            for point in scope:
                scores = calculate_kmer_features(sequence[point - w // 2: point + w // 2], signs)['Entropy']
                scores = scores.values.reshape(1, -1)
                scores = x_scaler.transform(scores)
                scores = LDA.transform(scores)
                res.append(model.predict_proba(scores)[0][1])

        res = np.array(res).reshape(len(cfg.add_params.windows), -1).mean(axis=0)

        prob_metric = res[lb // step:-rb // step].mean() / res.mean()

        logging.info(f'prob_metric for sample {graph}: {prob_metric:.3f}')

        # local_start_codon = [c for c in loc_start_codon if ((c >= i - lb) and (c <= j + rb))]
        # for codon in local_start_codon:
        #     ax.axvspan(codon, codon + 3, alpha=0.5, color='green')

        plt.title(f'Example №{graph} [prob_metric: {prob_metric:.3f}]')
        plt.plot(scope, res)
        plt.show()

    logging.info('Script done!')


if __name__ == "__main__":
    main()
