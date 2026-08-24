import torch

from model import BasalFrictionMLP


def predict():

    model = BasalFrictionMLP()

    model.load_state_dict(
        torch.load("basal_friction_model.pt")
    )

    model.eval()

    sample = torch.tensor(
        [[500.0, 1500.0, 1200.0]],
        dtype=torch.float32
    )

    prediction = model(sample)

    print(
        f"Predicted basal friction: "
        f"{prediction.item():.4f}"
    )


if __name__ == "__main__":
    predict()