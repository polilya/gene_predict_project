import random

import gff3_parser
import pandas as pd
from Bio import SeqIO


def parse_sequence(sequence_data):
    with open(sequence_data) as file:
        fasta_sequences = SeqIO.parse(file, "fasta")
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
    return sequence


def parse_data(gene_data):
    return gff3_parser.parse_gff3(gene_data)


def preprocess_data(data, sequence, sequence_length):
    cds = data[
        (data["Type"] == "CDS")
        & (data["Seqid"] == "NC_000001.11")
        & (data["Source"] == "BestRefSeq")
    ]

    cds.loc[:, "Start"] = cds.loc[:, "Start"].astype("int")
    cds.loc[:, "End"] = cds.loc[:, "End"].astype("int")
    cds = cds[:sequence_length]

    segments = []
    for _, row in cds.iterrows():
        start = row["Start"]
        end = row["End"]
        segment = sequence[
            start - 1 : end
        ]  # Извлечение отрезка, индексы начинаются с 0
        segments.append(segment)

    cds = cds.copy()
    cds.loc[:, "segments"] = segments
    cds["segment_length"] = cds["segments"].str.len()

    return cds


def compute_segment_lengths(cds):
    lengths = [len(segment) for segment in cds["segments"]]

    max_length = max(lengths)
    min_length = min(lengths)
    average_length = sum(lengths) / len(lengths)

    return max_length, min_length, average_length, lengths


def extract_segment(sequence, start, stop, length):
    # Вычисление максимально возможного отступа
    max_offset = length - (stop - start)

    # Генерация случайного отступа
    offset = random.randint(0, max_offset)

    # Вычисление границ внутреннего отрезка
    inner_start = start - offset
    inner_end = inner_start + length

    # Извлечение отрезка
    segment = sequence[inner_start:inner_end]

    # Создание маски
    mask = [0] * (length)
    mask[start - inner_start - 1 : stop - inner_start] = [1] * (stop - start + 1)

    return segment, mask


def print_sequence_with_mask(sequence, mask):
    # Проверка длины последовательности и маски
    if len(sequence) != len(mask):
        raise ValueError("Длина последовательности и маски должна быть одинаковой.")

    # Вывод последовательности с обведенной в квадратные скобки последовательностью символов
    is_prev_masked = False  # Флаг для отслеживания предыдущего символа
    for char, is_masked in zip(sequence, mask):
        if is_masked:
            if (
                not is_prev_masked
            ):  # Если предыдущий символ не был обведен, начинаем новую последовательность
                print("_____[", end="")
            print(char, end="")
            is_prev_masked = True
        else:
            if (
                is_prev_masked
            ):  # Если предыдущий символ был обведен, заканчиваем текущую последовательность
                print("]_____", end="")
            print(char, end="")
            is_prev_masked = False

    if (
        is_prev_masked
    ):  # Если последний символ в последовательности был обведен, заканчиваем последовательность
        print("]______", end="")


def extract_random_segment(sequence, length):
    # Генерация случайной позиции начала отрезка
    start = random.randint(0, len(sequence) - length)

    # Вычисление позиции конца отрезка
    end = start + length

    # Извлечение отрезка
    segment = sequence[start:end]

    return segment


def create_mask_array(true_masks):
    mask_length = len(true_masks[5])  # Длина каждой маски
    num_masks = len(true_masks)  # Количество масок
    mask_array = [
        [0] * mask_length for _ in range(num_masks)
    ]  # Создание массива масок из нулей

    return mask_array
