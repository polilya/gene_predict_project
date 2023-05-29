from Bio import SeqIO

import hydra
from hydra.core.config_store import ConfigStore
import numpy as np
import pandas as pd
from bisect import bisect

from preprocess_config import MLConfig
from tools import calculate_kmer_features, signs

np.random.seed(42)

cs = ConfigStore.instance()
cs.store(name='ml_config', node=MLConfig)

@hydra.main(config_path='conf', config_name='preprocess_config')
def main(cfg: MLConfig):
    sequence_data = f'{cfg.paths.data}/{cfg.files.sequence}'
    with open(sequence_data) as file:
        fasta_sequences = SeqIO.parse(file, 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)

    cds = pd.read_csv(f'{cfg.paths.data}/{cfg.files.cds}')
    exons = pd.read_csv(f'{cfg.paths.data}/{cfg.files.exons}')

    cds_starts = np.sort(list(set(cds['Start'].values)))
    cds_ends = np.sort(list(set(cds['End'].values)))

    exons_starts = np.sort(list(set(exons['Start'].values)))
    #exons_ends = np.sort(list(set(exons['Start'].values)))

    cds_starts_new = []
    cds_ends_new = []

    for start, end in zip(cds_starts, cds_ends):

        if sequence[start - 1:start + 2] in ['ATG']:
            cds_starts_new.append(start)
        if sequence[end - 3:end] in ['TAA', 'TAG', 'TGA']:
            cds_ends_new.append(end)

    seq_before_cds_location = []
    for cds_start in cds_starts_new:
        p = bisect(exons_starts, cds_start)
        seq_before_cds_location.append((exons_starts[p - 1], cds_start))
    seq_before_cds_location = np.array(seq_before_cds_location)

    mask = np.diff(seq_before_cds_location) > 50
    seq_before_cds_location = seq_before_cds_location[np.tile(mask, (1, 2))].reshape(-1, 2) - 1

    cds_location_data = pd.DataFrame(seq_before_cds_location, columns=['start', 'end'])
    cds_location_data.to_csv(f'{cfg.paths.data}/{cfg.files.cds_location}', index=False)

    seq_before_cds = []
    for i, j in seq_before_cds_location:
        seq_before_cds.append(sequence[i:j])

    seq_zero_location = []
    for sample in range(int(cfg.params.num_of_samples * seq_before_cds_location.shape[0])):
        b = np.random.randint(1e4, 2e8)
        bounds = (b, b + int(np.random.exponential(scale=200.0) + 50))

        flag = 0
        for bounds2 in seq_before_cds_location:
            r = max(bounds[0], bounds2[0]), min(bounds[1] + 1, bounds2[1] + 1)
            if set(range(*r)) != set():
                flag = 1
                break

        if flag == 0:
            seq_zero_location.append(bounds)

    seq_zero = []
    for i, j in seq_zero_location:
        seq_zero.append(sequence[i:j])

    start_data = pd.DataFrame(columns=signs)
    zero_data = pd.DataFrame(columns=signs)

    for cnt, seq in enumerate(seq_before_cds):
        start_data.loc[cnt] = calculate_kmer_features(seq, signs)['Entropy'].values

    for cnt, seq in enumerate(seq_zero):
        zero_data.loc[cnt] = calculate_kmer_features(seq, signs)['Entropy'].values

    start_data.loc[:, 'target'] = 1
    zero_data.loc[:, 'target'] = 0
    data = pd.concat([start_data, zero_data], axis=0)
    data.to_csv(f'{cfg.paths.data}/{cfg.files.final_data}', index=False)


if __name__ == "__main__":
    main()
