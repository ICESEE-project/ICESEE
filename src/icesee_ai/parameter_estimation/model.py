import torch
import torch.nn as nn


class BasalFrictionMLP(nn.Module):
    """
    Simple multilayer perceptron for estimating
    basal friction coefficients from surface observations.
    """

    def __init__(self, input_dim=3, hidden_dim=64, output_dim=1):
        super().__init__()

        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),

            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),

            nn.Linear(hidden_dim, output_dim)
        )

    def forward(self, x):
        return self.network(x)