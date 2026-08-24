import numpy as np
import torch
from torch.utils.data import Dataset


class BasalFrictionDataset(Dataset):
    """
    Synthetic dataset for basal friction estimation.

    Features:
        - surface velocity
        - ice thickness
        - surface elevation

    Target:
        - basal friction coefficient
    """

    def __init__(self, n_samples=1000):

        velocity = np.random.uniform(0, 1000, n_samples)
        thickness = np.random.uniform(100, 3000, n_samples)
        elevation = np.random.uniform(0, 4000, n_samples)

        friction = (
            0.002 * velocity
            + 0.0005 * thickness
            + 0.0002 * elevation
            + np.random.normal(0, 0.1, n_samples)
        )

        self.X = torch.tensor(
            np.vstack([velocity, thickness, elevation]).T,
            dtype=torch.float32
        )

        self.y = torch.tensor(
            friction.reshape(-1, 1),
            dtype=torch.float32
        )

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]