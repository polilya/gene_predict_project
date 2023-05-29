# полезные функции в отдельном пакете

import collections
import math
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

# In[3]:


signs = ['A', 'T', 'G', 'C', 'AT', 'GC', 'GA', 'GT', 'GG', 'CG', 'CT', 'CA', 'CC',
         'TG', 'TT', 'TA', 'TC', 'AC', 'AG', 'AA', 'ATG', 'TGA', 'TAA', 'TAG', 'TTT',
         'TCT', 'TAT', 'TGT', 'TTC', 'TCC', 'TAC', 'TGC', 'TTA', 'TCA', 'TTG', 'TCG',
         'TGG', 'CTT', 'CCT', 'CAT', 'CGT', 'CTC', 'CCC', 'CAC', 'CGC', 'CTA', 'CCA',
         'CAA', 'CGA', 'CTG', 'CCG', 'CAG', 'CGG', 'ATT', 'ACT', 'AAT', 'AGT', 'ATC',
         'ACC', 'AAC', 'AGC', 'ATA', 'ACA', 'AAA', 'AGA', 'ACG', 'AAG', 'AGG', 'GTT',
         'GCT', 'GAT', 'GGT', 'GTC', 'GCC', 'GAC', 'GGC', 'GTA', 'GCA', 'GAA', 'GGA', 'GTG', 'GCG', 'GAG', 'GGG']


# In[4]:


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


def apply_pca(X_train, X_test):
    x_scaler = StandardScaler()
    X_train = x_scaler.fit_transform(X_train)
    X_test = x_scaler.transform(X_test)

    pca = PCA(n_components=0.9)
    pca.fit(X_train_scaled)

    PC_train, PC_test = pca.transform(X_train_scaled), pca.transform(X_test_scaled)

    return PC_train, PC_test, pca, x_scaler


def apply_lda(X_train, Y_train, X_test):
    x_scaler = StandardScaler()
    X_train = x_scaler.fit_transform(X_train)
    X_test = x_scaler.transform(X_test)

    LDA = LinearDiscriminantAnalysis()
    LDA.fit(X_train, Y_train)

    X_train_new = LDA.transform(X_train)
    X_test_new = LDA.transform(X_test)

    return X_train_new, X_test_new, LDA, x_scaler


def st_scale(X):
    return StandardScaler().fit_transform(X)

