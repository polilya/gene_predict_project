import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from torch.autograd import Function
from torch.nn import GRU, Conv1d, Dropout, Linear, Module, ModuleList, ReLU, Sigmoid
from torch.utils.data import DataLoader, Dataset


# Вспомогательная функция для обучения модели
def train(model, dataloader, optimizer, criterion):
    model.train()
    running_loss = 0.0
    correct = 0
    total = 0

    for inputs, labels in dataloader:
        optimizer.zero_grad()
        outputs = model(inputs)
        labels = labels.unsqueeze(1)  # Изменение размерности labels
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        running_loss += loss.item()
        predicted = (outputs >= 0.5).squeeze().long()
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

    accuracy = correct / total
    return running_loss, accuracy


# Вспомогательная функция для проверки модели
def evaluate(model, dataloader, criterion):
    model.eval()
    running_loss = 0.0
    correct = 0
    total = 0

    with torch.no_grad():
        for inputs, labels in dataloader:
            outputs = model(inputs)
            labels = labels.unsqueeze(1)
            loss = criterion(outputs, labels)

            running_loss += loss.item()
            predicted = (outputs >= 0.5).squeeze().long()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()

    accuracy = correct / total
    return running_loss, accuracy


def train_model(
    model, train_data, train_labels, test_data, test_labels, num_epochs, batch_size
    ):
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters())

    model.train()  # Переводим модель в режим обучения

    train_dataset = torch.utils.data.TensorDataset(train_data, train_labels)
    train_dataloader = torch.utils.data.DataLoader(
        train_dataset, batch_size=batch_size, shuffle=True
    )

    test_dataset = torch.utils.data.TensorDataset(test_data, test_labels)
    test_dataloader = torch.utils.data.DataLoader(
        test_dataset, batch_size=batch_size, shuffle=False
    )

    for epoch in range(num_epochs):
        for batch_data, batch_labels in train_dataloader:
            optimizer.zero_grad()
            outputs = model(batch_data)
            loss = criterion(outputs, batch_labels)
            loss.backward()
            optimizer.step()

        print(f"Epoch [{epoch+1}/{num_epochs}], Training Loss: {loss.item()}")

        # Оценка модели на тестовой выборке
        accuracy = evaluate_model(model, test_dataloader)
        print(f"Epoch [{epoch+1}/{num_epochs}], Test Accuracy: {accuracy:.2%}")

    print("Training finished.")


def evaluate_model(model, dataloader):
    model.eval()  # Переводим модель в режим оценки (необязательно, если только оцениваем точность)

    correct = 0
    total = 0

    with torch.no_grad():  # Отключаем вычисление градиентов для оценки
        for batch_data, batch_labels in dataloader:
            outputs = model(batch_data)
            _, predicted = torch.max(outputs.data, 1)
            total += batch_labels.size(0)
            correct += (predicted == batch_labels).sum().item()

    accuracy = correct / total
    return accuracy


def train_model_N6(
    model, train_data, train_labels, test_data, test_labels, num_epochs, batch_size
):
    # Define loss function and optimizer
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    for epoch in range(num_epochs):
        for i in range(0, len(train_data), batch_size):
            batch_data = train_data[i : i + batch_size]
            batch_labels = train_labels[i : i + batch_size]
            batch_labels = batch_labels.unsqueeze(1).float()  # Add this line

            # Zero the gradients
            optimizer.zero_grad()

            # Perform forward pass
            outputs = model(batch_data)

            # Compute loss
            loss = criterion(outputs, batch_labels)

            # Perform backward pass and optimization
            loss.backward()
            optimizer.step()

        # Compute training accuracy
        _, predicted = torch.max(model(train_data).data, 1)
        train_acc = (predicted == train_labels).float().mean()

        # Compute test accuracy
        _, predicted = torch.max(model(test_data).data, 1)
        test_acc = (predicted == test_labels).float().mean()

        print(
            f"Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}, Train Acc: {train_acc:.4f}, Test Acc: {test_acc:.4f}"
        )


def objective(trial):
    # Определяем гиперпараметры, которые хотим оптимизировать
    lr = trial.suggest_loguniform("lr", 1e-5, 1e-1)
    batch_size = trial.suggest_categorical("batch_size", [64, 128, 256])
    dropout_rate = trial.suggest_uniform("dropout_rate", 0, 1)
    hidden_sizes = [trial.suggest_int("hidden_size_{}".format(i), 64, 512) for i in range(3)]  # Пример для трех слоев

    # Создаем DataLoader с выбранным batch_size
    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    # Создаем модель с выбранными гиперпараметрами
    model = ConvNet(input_size, hidden_sizes, output_size, window_sizes, paddings, dropout_rate)

    # Оптимизатор и функция потерь
    optimizer = optim.SGD(model.parameters(), lr=lr)
    criterion = nn.BCELoss()

    for epoch in range(num_epochs):
        # Обучение
        train_loss, train_accuracy = train(model, train_dataloader, optimizer, criterion)

        # Проверка
        val_loss, val_accuracy = evaluate(model, test_dataloader, criterion)

        trial.report(val_loss, epoch)

        # Если обучение не улучшается, можно остановить его раньше
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    return val_loss