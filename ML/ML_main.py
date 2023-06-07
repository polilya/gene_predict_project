# This is a sample Python script.

import os
import random

import gff3_parser
import keras
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
import torch

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from Bio import SeqIO
from hmmlearn import hmm

# from tensorflow.keras.preprocessing.sequence import pad_sequences
# from keras.preprocessing.sequence import pad_sequences
# from tensorflow.keras.preprocessing.sequence import pad_sequences
# from keras.preprocessing.sequence import pad_sequences
from keras.utils import pad_sequences
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from tensorflow.keras.layers import LSTM, Dense, Dropout, Embedding
from tensorflow.keras.models import Sequential
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.preprocessing.text import Tokenizer
from tensorflow.keras.utils import to_categorical
from torch.autograd import Function
from torch.nn import GRU, Conv1d, Dropout, Linear, Module, ModuleList, ReLU, Sigmoid
from torch.utils.data import DataLoader, Dataset

import classes
import dataparser
import model_utils
from classes import (
    ConvNet,
    Data,
    DNASequenceClassifier_N3,
    DNASequenceClassifier_N5,
    DNASequenceClassifier_N6,
    GuidedBackpropReLU,
    GuidedBackpropReLU_N6,
    GuidedGradCam,
    GuidedGradCam2,
    GuidedGradCam_N6,
)
from dataparser import (
    compute_segment_lengths,
    create_mask_array,
    extract_random_segment,
    extract_segment,
    parse_data,
    parse_sequence,
    preprocess_data,
    print_sequence_with_mask,
)
from model_utils import evaluate, evaluate_model, train, train_model, train_model_N6

# ML_main.py


# Пути к данным
gene_data = os.path.join("..", "data", "genomic.gff")
sequence_data = os.path.join("..", "data", "chromosome1.fasta")

# Используйте ваши функции из utils.py
sequence = parse_sequence(sequence_data)
data = parse_data(gene_data)

# Обработка и получение данных
cds = preprocess_data(data, sequence, 1000)
del data

# Вычисление и печать длины сегментов
max_length, min_length, average_length, lengths = compute_segment_lengths(cds)
print("Максимальная длина:", max_length)
print("Минимальная длина:", min_length)
print("Средняя длина:", average_length)

# Построение гистограммы
plt.hist(lengths, bins=30)
plt.xlabel("Длина")
plt.ylabel("Количество элементов")
plt.title("Распределение длин элементов")
# plt.show()

# Построение гистограммы с логарифмической осью X
plt.hist(lengths, bins=1000)
plt.xscale("log")
plt.xlabel("Длина")
plt.ylabel("Количество элементов")
plt.title("Распределение длин элементов (логарифмическая ось X)")
# plt.show()

# Ограничение длины цепочки
cds_short = cds[(cds["segments"].str.len() >= 9) & (cds["segments"].str.len() <= 300)]

# Получение значения 'End' из датафрейма cds_short
end_value = cds_short["End"]

# Вычисление длины отрезков
after_segment_length = 1000 - cds_short["segment_length"]
# print(after_segment_length)

after_segment_length_ch = after_segment_length.apply(
    lambda x: x + 1 if x % 2 != 0 else x
)
after_segment_length_ch = after_segment_length_ch * 0.5
print(after_segment_length_ch)

pre_segment_length = after_segment_length
pre_segment_length_ch = pre_segment_length.apply(lambda x: x - 1 if x % 2 != 0 else x)
pre_segment_length_ch = pre_segment_length_ch * 0.5
print(pre_segment_length_ch)

cds_short["pre_s_l"] = pre_segment_length_ch
cds_short["after_s_l"] = after_segment_length_ch

result = []
for _, row in cds_short.iterrows():
    start = int(row["Start"]) - int(row["pre_s_l"])  # Позиция начала отрезка
    end = int(row["Start"])  # Позиция конца отрезка
    substring = sequence[start:end]
    result.append(substring)

segments = []
masks = []
for _, row in cds_short.iterrows():
    start = row["Start"]
    stop = row["End"]
    segment, mask = extract_segment(sequence, start, stop, 401)
    # segment = extract_segment(sequence, start, stop, 1001)
    segments.append(segment)
    masks.append(mask)

seg = segments[0]
ma = masks[0]

print_sequence_with_mask(seg, ma)

cds_short["markered"] = segments
cds_short["masks"] = masks

rand_segments = []
for _ in range(len(cds_short)):
    segment = extract_random_segment(sequence, len(segments[5]))
    rand_segments.append(segment)

true_segments = cds_short["markered"].tolist()
true_masks = cds_short["masks"].tolist()

false_masks = create_mask_array(true_masks)
print(false_masks[924])

df = pd.DataFrame(
    {
        "Предложение": true_segments + rand_segments,
        "Метка": ["1"] * len(cds_short) + ["0"] * len(cds_short),
        "Маски": true_masks + false_masks,
    }
)
print(df)

# Проверка на равенство длин массивов
if len(cds_short["markered"]) != len(rand_segments):
    raise ValueError("Количество предложений в массивах не совпадает")

# Условие фильтрации для символа "N"
condition = df["Предложение"].str.contains("N")


# Фильтрация датафрейма
filtered_df = df[~condition]
filtered_df.shape

# Условие фильтрации для символа "N"
condition1 = filtered_df["Предложение"].str.contains("R")

# Фильтрация датафрейма
filtered_df1 = filtered_df[~condition1]
filtered_df1.shape

df_fin = filtered_df1

df_fin.to_csv("data_masks_rand.csv", index=False)

# Разделение данных на тренировочный и тестовый наборы
train_df, test_df = train_test_split(df, test_size=0.2, random_state=42)

# Вывод размерности тренировочного и тестового наборов
print("Размер тренировочного набора:", train_df.shape)
print("Размер тестового набора:", test_df.shape)

##Сеть 1
print("Начало сети 1")
# Создание модели HMM
model = hmm.MultinomialHMM(n_components=2)  # Количество скрытых состояний равно 2

# Обучающие данные
# sentences = ["This is a positive sentence",
#              "This is a negative sentence",
#              "Another positive sentence"]
sentences = train_df["Предложение"]

# Создание словаря уникальных токенов (отрезков по 3 символа)
tokens = set()
for sentence in sentences:
    for i in range(len(sentence) - 2):
        token = sentence[i : i + 3].lower()
        tokens.add(token)

# Создание индекса токенов в словаре
token2idx = {token: idx for idx, token in enumerate(tokens)}

# Создание словаря уникальных токенов и подсчет повторений
token_counts = {}
for sentence in sentences:
    for i in range(len(sentence) - 2):
        token = sentence[i : i + 3].lower()
        if token in token_counts:
            token_counts[token] += 1
        else:
            token_counts[token] = 1

# Вывод словаря уникальных токенов и количества повторений
for token, count in token_counts.items():
    print(f"Token: {token}, Count: {count}")

# Загрузка предобработанного датасета


# Разделение признаков и меток
X = df["Предложение"].values
y = df["Метка"].values

# Преобразование меток в числовой формат
label_encoder = LabelEncoder()
y = label_encoder.fit_transform(y)

# Преобразование геномных последовательностей в числовые последовательности
tokenizer = Tokenizer(
    num_words=500, char_level=True
)  # num_words - это максимальное количество слов/символов, которые должны быть учтены (в зависимости от выбранного уровня: char_level или word_level)
tokenizer.fit_on_texts(X)
X = tokenizer.texts_to_sequences(X)

# Разделение на обучающую и тестовую выборки
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Заполнение и выравнивание последовательностей до одинаковой длины
max_sequence_length = max(len(seq) for seq in X_train)
X_train = pad_sequences(X_train, maxlen=max_sequence_length)
X_test = pad_sequences(X_test, maxlen=max_sequence_length)


# Создание модели RNN
model = Sequential()
model.add(Embedding(input_dim=500, output_dim=32, input_length=max_sequence_length))
model.add(LSTM(64))
model.add(Dropout(0.5))
model.add(Dense(1, activation="sigmoid"))

# Компиляция модели
model.compile(loss="binary_crossentropy", optimizer="adam", metrics=["accuracy"])

# Обучение модели
model.fit(X_train, y_train, validation_data=(X_test, y_test), epochs=3, batch_size=32)

# Оценка модели
loss, accuracy = model.evaluate(X_test, y_test)
print("Потери:", loss)
print("Точность:", accuracy)
print("Конец сети 1")
##СЕТЬ 2
print("Начало сети 2")
# Гиперпараметры для настройки модели
EPOCHS = 15
BATCH_SIZE = 128
LR = 1e-3
# EMBED_DIM = 32
# HIDDEN_DIM = 32
TARGET_DIM = 2
# N_LAYER = 2
DROPOUT_PROB = 0.2

X = filtered_df1["Предложение"].values
y = filtered_df1["Метка"].values

# Сохранение данных в файлы
np.save("X.npy", X)
np.save("y.npy", y)

# Загрузка данных из файлов с разрешением pickle
X = np.load("X.npy", allow_pickle=True)
y = np.load("y.npy", allow_pickle=True)

batch_size = BATCH_SIZE
# Разделение данных на обучающий и тестовый наборы
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Создание объектов Dataset для обучающего и тестового наборов
train_dataset = Data(X_train, y_train)
test_dataset = Data(X_test, y_test)

# Создание DataLoader для обучающего и тестового наборов
train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

# Определяем модель
# Задаем параметры модели и обучения
input_size = 4  # Размерность входных данных
hidden_sizes = [128, 64, 32]  # Список размерностей скрытых слоев
output_size = 1  # Размерность выходных данных (бинарная классификация)
window_sizes = [3, 6, 9]  # Размеры окон
paddings = [0, 1, 2]  # Отступы от первого символа
dropout_rate = 0.2  # Коэффициент dropout
learning_rate = 0.01  # Скорость обучения
num_epochs = 10  # Количество эпох


# dataset = MyDataset(data, labels)

# Создаем DataLoader
batch_size = 64
shuffle = True
# dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)


# Создаем модель
model = ConvNet(
    input_size, hidden_sizes, output_size, window_sizes, paddings, dropout_rate
)

# Оптимизатор и функция потерь
optimizer = optim.SGD(model.parameters(), lr=learning_rate)
criterion = nn.BCELoss()

# Цикл обучения
for epoch in range(num_epochs):
    # Обучение
    train_loss, train_accuracy = train(model, train_dataloader, optimizer, criterion)

    # Проверка
    val_loss, val_accuracy = evaluate(model, test_dataloader, criterion)

    # Вывод результатов
    print(
        f"Epoch: {epoch+1}/{num_epochs}, Train Loss: {train_loss:.4f}, Train Accuracy: {train_accuracy:.4f}, Val Loss: {val_loss:.4f}, Val Accuracy: {val_accuracy:.4f}"
    )

# СЕТЬ 3
print("Начало сети 3")
# Загрузка данных из файлов с разрешением pickle
X = np.load("X.npy", allow_pickle=True)
y = np.load("y.npy", allow_pickle=True)

# Последовательности ДНК и метки классов
sequences = X
labels = y

# Создание словаря для кодирования
dictionary = {"A": 0, "C": 1, "G": 2, "T": 3}

# Преобразование последовательностей ДНК в векторы one-hot
max_sequence_length = max(len(seq) for seq in sequences)
num_classes = 2  # Количество классов (0 и 1)

encoded_sequences = np.zeros(
    (len(sequences), max_sequence_length, len(dictionary)), dtype=np.float32
)
for i, seq in enumerate(sequences):
    for j, nucleotide in enumerate(seq):
        encoded_sequences[i, j, dictionary[nucleotide]] = 1.0

# Преобразование меток классов в числовой формат
encoded_labels = np.array([int(label) for label in labels])

# Разделение данных на тренировочную и тестовую выборки
train_data, test_data, train_labels, test_labels = train_test_split(
    encoded_sequences, encoded_labels, test_size=0.2, random_state=42
)

# Преобразование в тензоры PyTorch
train_data = torch.tensor(train_data)
train_labels = torch.tensor(train_labels)
test_data = torch.tensor(test_data)
test_labels = torch.tensor(test_labels)

"""# Создание модели
input_size = train_data.size(1)
hidden_size = 128
num_classes = 2
model = DNASequenceClassifier_N3(input_size, hidden_size, num_classes)

# Обучение модели
num_epochs = 10
batch_size = 32
train_model(
    model,
    train_data,
    train_labels,
    test_data,
    test_labels,
    num_epochs,
    batch_size,
)
evaluate_model(model,
    train_data,
    train_labels,
    test_data,
    test_labels,
    num_epochs,
    batch_size,)"""

print("Конец сети 3")
# СЕТЬ 4
print("Начало сети 4")
print("Конец сети 4")
# СЕТЬ 5
print("Начало сети 5")
# Загрузка данных из файлов с разрешением pickle
"""X = np.load("X.npy", allow_pickle=True)
y = np.load("y.npy", allow_pickle=True)

# train_df, test_df = train_test_split(df, test_size=0.2, random_state=42)

X_small = X[:1000]
y_small = y[:1000]
X_small.shape

# Load DNA sequence data and labels (example using DeepBind dataset)
sequences = X_small
labels = y_small

# Create a dictionary for encoding nucleotides
dictionary = {"A": 0, "C": 1, "G": 2, "T": 3}

# Preprocess the DNA sequences
max_sequence_length = max(len(seq) for seq in sequences)
num_classes = 2
encoded_sequences = np.zeros(
    (len(sequences), max_sequence_length, len(dictionary)), dtype=np.float32
)
for i, seq in enumerate(sequences):
    for j, nucleotide in enumerate(seq):
        encoded_sequences[i, j, dictionary[nucleotide]] = 1.0

# Convert labels to numerical format
encoded_labels = np.array([int(label) for label in labels])

# Split the data into training and testing sets
train_data, test_data, train_labels, test_labels = train_test_split(
    encoded_sequences, encoded_labels, test_size=0.2, random_state=42
)

# Convert data to PyTorch tensors
train_data = torch.tensor(train_data)
train_labels = torch.tensor(train_labels)
test_data = torch.tensor(test_data)
test_labels = torch.tensor(test_labels)

# Create an instance of the DNASequenceClassifier model
input_size = train_data.size(2)
hidden_size = 32
num_classes = 2
model = DNASequenceClassifier_N5(input_size, hidden_size, num_classes)

# Define training parameters
num_epochs = 15
batch_size = 32

# Train the model
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

model.train()
for epoch in range(num_epochs):
    optimizer.zero_grad()
    outputs = model(train_data)
    loss = criterion(outputs.squeeze(), train_labels)
    loss.backward()
    optimizer.step()

    print(f"Epoch [{epoch+1}/{num_epochs}], Training Loss: {loss.item()}")

# Create an instance of GuidedGradCam
guided_grad_cam = GuidedGradCam(model, target_layer="fc")

# Select a random test sample
sample_idx = np.random.randint(len(test_data))
sample_data = test_data[sample_idx].unsqueeze(0)
sample_label = test_labels[sample_idx]


# Generate the Guided Grad-CAM visualization
cam = guided_grad_cam.generate(target_class=sample_label)


# Plot the DNA sequence and its corresponding Guided Grad-CAM
plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.imshow(sample_data.squeeze().T, cmap="binary", aspect="auto")
plt.title("DNA Sequence")
plt.xlabel("Position")
plt.ylabel("Nucleotide")
plt.subplot(1, 2, 2)
plt.imshow(cam, cmap="jet", aspect="auto")
plt.title("Guided Grad-CAM")
plt.xlabel("Position")
plt.ylabel("Importance")
plt.colorbar()
plt.tight_layout()
plt.show()


# Create an instance of GuidedGradCam
guided_grad_cam = GuidedGradCam2(model, target_layer="fc")

# ...

# Generate the Guided Grad-CAM visualization
cam = guided_grad_cam.generate(target_class=sample_label)"""

print("Конец сети 5")
# СЕТЬ 6
print("Начало сети 6")
# Load DNA sequence data and labels (example using DeepBind dataset)
X_small = X[:1000]
y_small = y[:1000]
X_small.shape

sequences = X_small
labels = y_small

# Create a dictionary for encoding nucleotides
dictionary = {"A": 0, "C": 1, "G": 2, "T": 3}

# Preprocess the DNA sequences
max_sequence_length = max(len(seq) for seq in sequences)
num_classes = 2
encoded_sequences = np.zeros(
    (len(sequences), max_sequence_length, len(dictionary)), dtype=np.float32
)
for i, seq in enumerate(sequences):
    for j, nucleotide in enumerate(seq):
        encoded_sequences[i, j, dictionary[nucleotide]] = 1.0

# Convert labels to numerical format
encoded_labels = np.array([int(label) for label in labels])

# Split the data into training and testing sets
train_data, test_data, train_labels, test_labels = train_test_split(
    encoded_sequences, encoded_labels, test_size=0.2, random_state=42
)

# Convert data to PyTorch tensors
train_data = torch.tensor(train_data)
train_labels = torch.tensor(train_labels)
test_data = torch.tensor(test_data)
test_labels = torch.tensor(test_labels)

# Create an instance of the DNASequenceClassifier model
input_size = train_data.size(2)
hidden_size = 32
num_classes = 1
model = DNASequenceClassifier_N6(input_size, hidden_size, num_classes)

# Define training parameters
num_epochs = 10
batch_size = 32


# Train the model
train_model_N6(
    model, train_data, train_labels, test_data, test_labels, num_epochs, batch_size
)

print("Конец сети 6")
# СЕТЬ 7
print("Начало сети 7")
print("Конец сети 7")
