#!/usr/bin/env python
# coding: utf-8

import collections
import math
import pandas as pd
import numpy as np
from tqdm import tqdm

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis 


signs = ['A', 'T', 'G', 'C', 'AT', 'GC', 'GA', 'GT', 'GG', 'CG', 'CT', 'CA', 'CC', 
         'TG', 'TT', 'TA', 'TC', 'AC', 'AG', 'AA', 'ATG', 'TGA', 'TAA', 'TAG', 'TTT', 
         'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TTG', 'TCG', 
         'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA', 
         'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT', 'ATC', 
         'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ACG', 'AAG', 'AGG', 'GTT', 
         'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 
         'GTG', 'GCG', 'GAG', 'GGG']


def calculate_kmer_features(sequence, signs):
    # получаем количество всех k-мер в последовательности
    counter = {sign: sequence.count(sign) for sign in signs}
    
    # вычисляем общее количество k-мер
    total = sum(counter.values())
    
    # вычисляем частоту каждого k-мера
    frequency = {kmer: count / total if count > 0 else 0 for kmer, count in counter.items()}
    
    # вычисляем энтропию
    entropy = {kmer: -p_kmer * math.log(p_kmer, 4) if p_kmer > 0 else 0 for kmer, p_kmer in frequency.items()}

    return pd.DataFrame({
            'Kmer': list(frequency.keys()),
            'Frequency': list(frequency.values()),
            'Entropy': list(entropy.values()),
                        })


def get_values(data):
    return np.sort(list(set(data.values)))


def check_intersection(bounds, seq_before_location, seq_after_location):
    
    flag = 0
    
    for bounds2 in np.array(seq_before_location):
        r = max(bounds[0], bounds2[0]), min(bounds[1]+1, bounds2[1]+1)
        if set(range(*r)) != set():
            flag = 1
            break
         
    if flag == 1:
        return False
    
    for bounds2 in np.array(seq_after_location):
        r = max(bounds[0], bounds2[0]), min(bounds[1]+1, bounds2[1]+1)
        if set(range(*r)) != set():
            flag = 1
            break
    
    if flag == 1:
        return False
    
    elif flag == 0:
        return True


def check_self_intersection(seq_zero_location):
    seq_zero_location_new = []

    for bounds1 in tqdm(np.array(seq_zero_location)):
        seq_zero_location = seq_zero_location[1:]
        flag = 0
        for bounds2 in np.array(seq_zero_location):
            r = max(bounds1[0], bounds2[0]), min(bounds1[1] + 1, bounds2[1] + 1)
            if set(range(*r)) != set():
                flag = 1
                break
        if flag == 0: seq_zero_location_new.append(bounds1)

    return pd.DataFrame(seq_zero_location_new, columns=['start', 'end'])


def get_sequence(location_data, sequence):
    sequence_data = []
    for i, j in np.array(location_data):
        sequence_data.append(sequence[i:j])
    return sequence_data


def limit(location, min_length, max_length):
    
    mask_l = np.diff(location) > min_length
    mask_h = np.diff(location) < max_length
    mask = mask_l & mask_h
    
    limited_location = location[np.tile(mask, (1, 2))].reshape(-1, 2)-1
    return pd.DataFrame(limited_location, columns=['start', 'end'])


def calculate_features(in_data):
    data = pd.DataFrame(columns=signs)
    for cnt, seq in enumerate(in_data):
        data.loc[cnt] = calculate_kmer_features(seq, signs)['Entropy'].values
        
    return data


def apply_pca(X_train, X_test, Y_train=0):
    
    x_scaler = StandardScaler()
    X_train_scaled = x_scaler.fit_transform(X_train)
    X_test_scaled = x_scaler.transform(X_test)

    pca = PCA(n_components=0.9)
    pca.fit(X_train_scaled)

    PC_train, PC_test = pca.transform(X_train_scaled), pca.transform(X_test_scaled)

    return PC_train, PC_test, pca, x_scaler


def apply_lda(X_train, X_test, Y_train):
     
    
    x_scaler = StandardScaler()
    X_train = x_scaler.fit_transform(X_train)
    X_test = x_scaler.transform(X_test)
    
    LDA = LinearDiscriminantAnalysis() 
    LDA.fit(X_train, Y_train)
    
    X_train_new = LDA.transform(X_train)
    X_test_new = LDA.transform(X_test)
    
    return X_train_new, X_test_new, LDA, x_scaler
