import logging
from tqdm import tqdm

from Bio import SeqIO
import numpy as np
import pandas as pd
from bisect import bisect

import hydra
from hydra.core.config_store import ConfigStore
from config import MLConfig

import tools

np.random.seed(17)

cs = ConfigStore.instance()
cs.store(name='ml_config', node=MLConfig)


@hydra.main(config_path='conf', config_name='config_ppr', version_base=None)
def main(cfg: MLConfig):

    logging.info('Preprocessing starts')

    sequence_data = f'{cfg.paths.data}/{cfg.files.sequence}'
    with open(sequence_data) as file:
        fasta_sequences = SeqIO.parse(file, 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)

    cds = pd.read_csv(f'{cfg.paths.data}/{cfg.files.cds}')
    exons = pd.read_csv(f'{cfg.paths.data}/{cfg.files.exons}')

    logging.info(f'Total number of cds: {len(cds)}')
    logging.info(f'Total number of exons: {len(exons)}')

    cds_starts = tools.get_values(cds['Start'])
    cds_ends = tools.get_values(cds['End'])

    exons_starts = tools.get_values(exons['Start'])
    exons_ends = tools.get_values(exons['End'])

    cds_starts_new = []
    cds_ends_new = []

    for start, end in zip(cds_starts, cds_ends):

        if sequence[start - 1:start + 2] in ['ATG']:
            cds_starts_new.append(start)
        if sequence[end - 3:end] in ['TAA', 'TAG', 'TGA']:
            cds_ends_new.append(end)

    logging.info(f'Number of first cds: {len(cds_starts_new)}')
    logging.info(f'Number of last cds: {len(cds_ends_new)}')

    seq_before_cds_location = []
    for cds_start in cds_starts_new:
        p = bisect(exons_starts, cds_start)
        seq_before_cds_location.append((exons_starts[p - 1], cds_start))
    seq_before_cds_location = np.array(seq_before_cds_location)

    seq_after_cds_location = []
    for cds_end in cds_ends_new:
        p = bisect(exons_ends, cds_end)
        seq_after_cds_location.append((cds_end, exons_ends[p]))
    seq_after_cds_location = np.array(seq_after_cds_location)

    start_cds_location_data = tools.limit(seq_before_cds_location,
                                          min_length=cfg.preprocess_params.start_min_length,
                                          max_length=cfg.preprocess_params.start_max_length)
    start_cds_location_data.to_csv(f'{cfg.paths.data}/{cfg.files.first_cds_location}', index=False)

    end_cds_location_data = tools.limit(seq_after_cds_location,
                                        min_length=cfg.preprocess_params.end_min_length,
                                        max_length=cfg.preprocess_params.end_max_length)
    end_cds_location_data.to_csv(f'{cfg.paths.data}/{cfg.files.last_cds_location}', index=False)

    logging.info(f'Number of first cds after limiting: {len(start_cds_location_data)}')
    logging.info(f'Number of last cds after limiting: {len(end_cds_location_data)}')

    seq_before_cds = tools.get_sequence(start_cds_location_data, sequence)
    seq_after_cds = tools.get_sequence(end_cds_location_data, sequence)

    n = cfg.preprocess_params.num_of_samples
    seq_zero_location = []
    for _ in tqdm(range(n)):
        b = np.random.randint(1e4, 2e8)
        scale = np.random.randint(400, 800)
        L = int(np.random.exponential(scale=scale) + 50)
        bounds = (b, b + L)
        flag = tools.check_intersection(bounds, start_cds_location_data, end_cds_location_data)
        if flag:
            seq_zero_location.append(bounds)

    seq_zero_location_l = len(seq_zero_location)
    logging.info(f'Number of excluded samples [1-0 intersection]: {n - seq_zero_location_l}')

    seq_zero_location_new = tools.check_self_intersection(seq_zero_location)
    seq_zero_location_new_l = len(seq_zero_location_new)
    logging.info(f'Number of excluded samples [0-0 intersection]: {seq_zero_location_l - seq_zero_location_new_l}')

    logging.info(f'Total number of excluded samples: {n - seq_zero_location_new_l}')

    seq_zero = tools.get_sequence(seq_zero_location_new, sequence)



    start_data = tools.calculate_features(seq_before_cds)
    end_data = tools.calculate_features(seq_after_cds)
    zero_data = tools.calculate_features(seq_zero)

    logging.info(f'Number of samples: first cds: {start_data.shape[0]}, last cds: {end_data.shape[0]}, zero_class: {zero_data.shape[0]}')

    start_data.loc[:, 'target'] = 1
    end_data.loc[:, 'target'] = 0
    zero_data.loc[:, 'target'] = 0

    start_data_final = pd.concat([start_data, end_data, zero_data], axis=0)

    start_data.loc[:, 'target'] = 0
    end_data.loc[:, 'target'] = 1
    zero_data.loc[:, 'target'] = 0

    end_data_final = pd.concat([end_data, start_data, zero_data], axis=0)

    start_data_final.to_csv(f'{cfg.paths.data}/{cfg.files.start_features}', index=False)
    end_data_final.to_csv(f'{cfg.paths.data}/{cfg.files.end_features}', index=False)

    logging.info('Preprocessing done!')


if __name__ == "__main__":
    main()
