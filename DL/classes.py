# Импортируем необходимые библиотеки
import numpy as np
import torch
import torch.nn as nn
from torch.autograd import Function
from torch.nn import GRU, Conv1d, Dropout, Linear, Module, ModuleList, ReLU, Sigmoid
from torch.utils.data import Dataset


# Класс Data предназначен для подготовки данных для обучения.
# Он преобразует последовательности символов в one-hot вектора.
class Data(Dataset):
    def __init__(self, X, y):
        self.X = self._process_sequences(X)
        self.y = torch.from_numpy(y.astype(np.float32))
        self.len = self.X.shape[0]

    def __getitem__(self, index):
        return self.X[index], self.y[index]

    def __len__(self):
        return self.len

    def _process_sequences(self, sequences):
        alphabet = "ACGT"
        char_to_index = {char: index for index, char in enumerate(alphabet)}

        one_hot_vectors = []
        for sequence in sequences:
            sequence_one_hot = []
            for char in sequence:
                vector = [0] * 4
                vector[char_to_index[char]] = 1
                sequence_one_hot.append(vector)
            one_hot_vectors.append(sequence_one_hot)

        return torch.tensor(one_hot_vectors, dtype=torch.float32)


# Класс ConvNet представляет собой сверточную нейронную сеть.
class ConvNet(Module):
    def __init__(
        self,
        input_size,
        hidden_sizes,
        output_size,
        window_sizes,
        paddings,
        dropout_rate,
    ):
        super(ConvNet, self).__init__()
        self.convs = ModuleList(
            [
                Conv1d(
                    in_channels=input_size,
                    out_channels=64,
                    kernel_size=window_size,
                    padding=padding,
                )
                for window_size, padding in zip(window_sizes, paddings)
            ]
        )

        layers = []
        prev_size = len(window_sizes) * 64
        for hidden_size in hidden_sizes:
            layers.append(Linear(prev_size, hidden_size))
            layers.append(Sigmoid())
            layers.append(Dropout(dropout_rate))
            prev_size = hidden_size
        layers.append(Linear(prev_size, output_size))
        layers.append(Sigmoid())

        self.fc = nn.Sequential(*layers)

    def forward(self, x):
        x = x.permute(0, 2, 1)
        out = []
        for conv in self.convs:
            conv_out = torch.relu(conv(x))
            conv_out, _ = torch.max(conv_out, dim=2)
            out.append(conv_out)
        out = torch.cat(out, dim=1)
        out = out.view(out.size(0), -1)
        out = self.fc(out)
        return out


# Класс DNASequenceClassifier_N3 представляет собой модель нейронной сети для классификации последовательностей ДНК.
# Нейронная сеть включает в себя сверточный слой, нелинейность (ReLU) и полносвязный слой.
class DNASequenceClassifier_N3(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        # input_size - размерность входных данных
        # hidden_size - размер скрытого слоя
        # num_classes - количество классов для классификации
        super(DNASequenceClassifier_N3, self).__init__()
        self.conv1 = nn.Conv1d(input_size, hidden_size, kernel_size=3, padding=1)
        self.relu = nn.ReLU()
        self.fc = nn.Linear(hidden_size, num_classes)

    def forward(self, x):
        # Применяются слои сети: свертка, нелинейность и полносвязный слой
        x = self.conv1(x)
        x = self.relu(x)
        x = torch.mean(x, dim=2)  # Усреднение по временной оси
        x = self.fc(x)
        return x


# Определяем модель DNASequenceClassifier
class DNASequenceClassifier_N5(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super().__init__()
        self.hidden_size = hidden_size
        self.rnn = nn.GRU(
            input_size, hidden_size, batch_first=True
        )  # Определение слоя GRU
        self.fc = nn.Linear(hidden_size, num_classes)  # Определение полносвязного слоя

    def forward(self, x):
        _, h_n = self.rnn(x)  # Пропускание входных данных через слой GRU
        out = self.fc(
            h_n.squeeze(0)
        )  # Пропускание выходных данных GRU через полносвязный слой
        return out  # Возврат результата


class GuidedBackpropReLU(Function):
    @staticmethod
    def forward(ctx, input):
        positive_mask = (input > 0).type_as(
            input
        )  # Создание маски для положительных входных данных
        output = torch.addcmul(
            torch.zeros(input.size()).type_as(input), input, positive_mask
        )  # Применение функции ReLU к входным данным
        ctx.save_for_backward(
            output
        )  # Сохранение выходных данных для обратного прохода
        return output

    @staticmethod
    def backward(ctx, grad_output):
        (output,) = ctx.saved_tensors  # Извлечение сохраненных тензоров
        grad_input = None

        positive_mask_1 = (output > 0).type_as(
            grad_output
        )  # Создание маски для положительных выходных данных
        positive_mask_2 = (grad_output > 0).type_as(
            grad_output
        )  # Создание маски для положительных градиентов
        grad_input = torch.addcmul(
            torch.zeros(output.size()).type_as(grad_output),
            torch.addcmul(
                torch.zeros(output.size()).type_as(grad_output),
                grad_output,
                positive_mask_1,
            ),
            positive_mask_2,
        )  # Нулирование градиентов там, где либо вход, либо градиенты были отрицательными

        return grad_input  # Возврат результата обратного прохода


class GuidedGradCam:
    def __init__(self, model, target_layer):
        self.model = model  # Сохраняем модель
        self.target_layer = (
            target_layer  # Сохраняем целевой слой, который мы хотим визуализировать
        )
        self.gradients = None  # Инициализируем градиенты

        self.model.eval()  # Переводим модель в режим оценки
        self.register_hooks()  # Регистрируем хуки

    def register_hooks(self):
        target_layer = self.model._modules.get(
            self.target_layer
        )  # Получаем целевой слой

        # Функция, которая будет вызываться при прямом проходе
        def forward_hook(module, input, output):
            self.feature_maps = output  # Сохраняем карты признаков

        # Функция, которая будет вызываться при обратном проходе
        def backward_hook(module, grad_input, grad_output):
            self.gradients = grad_output[0]  # Сохраняем градиенты

        # Регистрируем функции хуков для прямого и обратного проходов
        target_layer.register_forward_hook(forward_hook)
        target_layer.register_backward_hook(backward_hook)

    def forward(self, x):
        return self.model(x)  # Производим прямой проход через модель

    def generate(self, target_class):
        self.model.zero_grad()  # Обнуляем градиенты модели
        output = self.forward(input)  # Вызываем прямой проход с входным тензором
        one_hot = torch.zeros((1, output.size()[-1]), dtype=torch.float32)
        one_hot[0][target_class] = 1.0  # Создаем вектор классов
        one_hot = torch.sum(one_hot * output)

        one_hot.backward(retain_graph=True)  # Вызываем обратный проход
        grads = self.gradients.cpu().data.numpy()  # Получаем градиенты
        feature_maps = self.feature_maps.cpu().data.numpy()[
            0
        ]  # Получаем карты признаков

        weights = np.mean(grads, axis=2)  # Вычисляем средние веса по каналам
        cam = np.zeros(feature_maps.shape[1:], dtype=np.float32)

        # Вычисляем взвешенную сумму карт признаков
        for i, w in enumerate(weights):
            cam += w * feature_maps[i, :]

        # Нормализуем полученную карту активации
        cam = np.maximum(cam, 0)
        cam = cam - np.min(cam)
        cam = cam / np.max(cam)

        return cam


class GuidedGradCam2:
    # ...

    def forward(self, x):
        return self.model(x)  # Use the model instance directly

    def generate(self, target_class):
        self.model.zero_grad()
        output = self.forward(x)  # Call forward with the input tensor
        one_hot = torch.zeros((1, output.size()[-1]), dtype=torch.float32)
        one_hot[0][target_class] = 1.0
        one_hot = torch.sum(one_hot * output)

        one_hot.backward(retain_graph=True)
        grads = self.gradients.cpu().data.numpy()
        feature_maps = self.feature_maps.cpu().data.numpy()[0]

        weights = np.mean(grads, axis=2)
        cam = np.zeros(feature_maps.shape[1:], dtype=np.float32)

        for i, w in enumerate(weights):
            cam += w * feature_maps[i, :]

        cam = np.maximum(cam, 0)
        cam = cam - np.min(cam)
        cam = cam / np.max(cam)

        return cam


# Define the DNASequenceClassifier model
class DNASequenceClassifier_N6(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(DNASequenceClassifier_N6, self).__init__()
        self.hidden_size = hidden_size
        self.rnn = nn.GRU(input_size, hidden_size, batch_first=True)
        self.fc = nn.Linear(hidden_size, num_classes)

    def forward(self, x):
        _, h_n = self.rnn(x)
        out = self.fc(h_n.squeeze(0))
        return out


class GuidedBackpropReLU_N6(Function):
    @staticmethod
    def forward(ctx, input):
        positive_mask = (input > 0).type_as(input)
        output = torch.addcmul(
            torch.zeros(input.size()).type_as(input), input, positive_mask
        )
        ctx.save_for_backward(output)
        return output

    @staticmethod
    def backward(ctx, grad_output):
        (output,) = ctx.saved_tensors
        grad_input = None

        positive_mask_1 = (output > 0).type_as(grad_output)
        positive_mask_2 = (grad_output > 0).type_as(grad_output)
        grad_input = torch.addcmul(
            torch.zeros(output.size()).type_as(grad_output),
            torch.addcmul(
                torch.zeros(output.size()).type_as(grad_output),
                grad_output,
                positive_mask_1,
            ),
            positive_mask_2,
        )

        return grad_input


class GuidedGradCam_N6:
    def __init__(self, model, target_layer):
        self.model = model
        self.target_layer = target_layer
        self.gradients = None

        self.model.eval()
        self.register_hooks()

    def register_hooks(self):
        target_layer = self.model._modules.get(self.target_layer)

        def forward_hook(module, input, output):
            self.feature_maps = output

        def backward_hook(module, grad_input, grad_output):
            self.gradients = grad_output[0]

        target_layer.register_forward_hook(forward_hook)
        target_layer.register_backward_hook(backward_hook)

    def forward(self, x):
        return self.model(x)  # Use the model instance directly

    def generate(self, target_class):
        self.model.zero_grad()
        output = self.forward(input)  # Call forward with the input tensor
        one_hot = torch.zeros((1, output.size()[-1]), dtype=torch.float32)
        one_hot[0][target_class] = 1.0
        one_hot = torch.sum(one_hot * output)

        one_hot.backward(retain_graph=True)
        grads = self.gradients.cpu().data.numpy()
        feature_maps = self.feature_maps.cpu().data.numpy()[0]

        weights = np.mean(grads, axis=2)
        cam = np.zeros(feature_maps.shape[1:], dtype=np.float32)

        for i, w in enumerate(weights):
            cam += w * feature_maps[i, :]

        cam = np.maximum(cam, 0)
        cam = cam - np.min(cam)
        cam = cam / np.max(cam)

        return cam
