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

#np.random.seed(42)

from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score

cs = ConfigStore.instance()
cs.store(name='ml_config', node=MLConfig)


@hydra.main(config_path='conf', config_name='config_start', version_base=None)
def main(cfg: MLConfig):

    logging.info('Model.py starts')
    
    sequence_data = f'{cfg.paths.data}/{cfg.files.sequence}'
    with open(sequence_data) as file:
        fasta_sequences = SeqIO.parse(file, 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)

    data = pd.read_csv(f'{cfg.paths.data}/{cfg.files.features}')
    cds_location_data = pd.read_csv(f'{cfg.paths.data}/{cfg.files.cds_location}')
    cds_location_data = np.array(cds_location_data)
    codon_location = pd.read_csv(f'{cfg.paths.data}/{cfg.files.codon_loc}').values.squeeze()

    Y = data['target']
    X = data.drop(['target'], axis=1)
    type_of_data = cfg.add_params.data

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3, random_state=42, stratify=Y, shuffle=True)

    if type_of_data == 'test':
        indices = X_test.index[X_test.index < cds_location_data.shape[0]].values
    elif type_of_data == 'train':
        indices = X_train.index[X_train.index < cds_location_data.shape[0]].values

    X_train, X_test, dim_reduction, x_scaler = tools.apply_lda(X_train.values, X_test.values, Y_train.values)

    model = SVC(**OmegaConf.to_container(cfg.params))

    logging.info('Start training...')
    model.fit(X_train, Y_train)
    logging.info('Training done')
    prediction = model.predict(X_test)

    logging.info(f'f1 score: {f1_score(Y_test, prediction, average="binary"):.3f}')
    logging.info('Start evaluation...')

    step = cfg.add_params.step
    plot_status = cfg.add_params.plot

    num_of_plots = cfg.add_params.num_of_plots

    for graph in list(np.random.choice(indices, size=num_of_plots)):

        i, j = cds_location_data[graph]
        length = j - i

        scope = np.arange(i - 2*length, j + 2*length, step).astype('int')

        res = []
        
        for w in cfg.add_params.windows:
            for point in scope:
                scores = calculate_kmer_features(sequence[point - w // 2: point + w // 2], signs)['Entropy']
                scores = scores.values.reshape(1, -1)
                scores = x_scaler.transform(scores)
                scores = dim_reduction.transform(scores)
                #res.append(model.predict_proba(scores)[0][1])
                res.append(model.predict(scores)[0])

        res = np.array(res).reshape(len(cfg.add_params.windows), -1).mean(axis=0)

        border_codons_location = [c for c in codon_location if ((c >= i - 2*length) and (c <= j + 2*length))]

        if plot_status:
            fig, ax = plt.subplots(dpi=120)
            ax.axvspan(i, j, alpha=0.6, color='red', label='Positive area')
            plt.title(f'Example №{graph} [{type_of_data}]')
            plt.ylim(-0.1, 1.1)
            plt.ylabel('Predicted class')
            plt.xlabel('Sequence number of the nucleotide')
            plt.plot(scope, res, color='darkblue', linestyle='--')
            plt.scatter(scope, res, color='blue', s=20, label='Predicted class')
            for loc in border_codons_location:
                ax.axvspan(loc, loc + 3, alpha=0.8, color='darkgreen', label='Start codon location')
            plt.show()

    logging.info('Evaluation done!')


if __name__ == "__main__":
    main()
