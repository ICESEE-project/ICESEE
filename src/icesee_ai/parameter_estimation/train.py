import torch
from torch.utils.data import DataLoader
import torch.nn as nn
import torch.optim as optim

from model import BasalFrictionMLP
from dataset import BasalFrictionDataset


def train():

    dataset = BasalFrictionDataset(n_samples=5000)

    loader = DataLoader(
        dataset,
        batch_size=64,
        shuffle=True
    )

    model = BasalFrictionMLP()

    criterion = nn.MSELoss()

    optimizer = optim.Adam(
        model.parameters(),
        lr=1e-3
    )

    epochs = 50

    for epoch in range(epochs):

        running_loss = 0.0

        for X, y in loader:

            optimizer.zero_grad()

            predictions = model(X)

            loss = criterion(predictions, y)

            loss.backward()

            optimizer.step()

            running_loss += loss.item()

        print(
            f"Epoch {epoch+1}/{epochs} "
            f"Loss: {running_loss/len(loader):.6f}"
        )

    torch.save(
        model.state_dict(),
        "basal_friction_model.pt"
    )

    print("\nModel saved: basal_friction_model.pt")


if __name__ == "__main__":
    train()