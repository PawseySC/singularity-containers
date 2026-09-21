#!/usr/bin/env python3

# A small single-GPU PyTorch training example: a fully-connected network
# classifying FashionMNIST images.  Adapted from Pawsey's PyTorch
# documentation for use with a Singularity/Apptainer container.

import os

import torch
from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets
from torchvision.transforms import ToTensor


# Program parameters

device = "cuda" if torch.cuda.is_available() else "cpu"
batch_size = 64
n_epochs = 10

# Where to cache the (small) FashionMNIST dataset.  Sourced from the
# environment so the Slurm batch script controls it; falls back to a
# local directory when run interactively.
data_path = os.environ.get("DATA_DIR", os.path.join(os.getcwd(), "mnist_data"))


class NeuralNetwork(nn.Module):
    def __init__(self):
        super().__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(28 * 28, 512),
            nn.ReLU(),
            nn.Linear(512, 512),
            nn.ReLU(),
            nn.Linear(512, 10)
        )

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits


def train(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    model.train()
    for batch, (X, y) in enumerate(dataloader):
        X, y = X.to(device), y.to(device)

        # Compute prediction error
        pred = model(X)
        loss = loss_fn(pred, y)

        # Backpropagation
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()

        if batch % 100 == 0:
            loss, current = loss.item(), (batch + 1) * len(X)
            print(f"loss: {loss:>7f} [{current:>5d}/{size:>5d}]")


if __name__ == "__main__":
    print("Device:", device)

    training_data = datasets.FashionMNIST(root=data_path, train=True, download=True, transform=ToTensor())
    test_data = datasets.FashionMNIST(root=data_path, train=False, download=True, transform=ToTensor())
    train_dataloader = DataLoader(training_data, batch_size=batch_size)
    test_dataloader = DataLoader(test_data, batch_size=batch_size)

    model = NeuralNetwork().to(device)
    print(model)

    loss_fn = nn.CrossEntropyLoss()
    optimizer = torch.optim.SGD(model.parameters(), lr=1e-3)

    for epoch in range(n_epochs):
        print(f"\n### Epoch {epoch}/{n_epochs} ###")
        train(train_dataloader, model, loss_fn, optimizer)
