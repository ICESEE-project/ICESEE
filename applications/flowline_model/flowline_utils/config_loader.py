# =============================================================================
# src/project_name/config_loader.py
# @author: Brian Kyanjo
# @date: 2024-09-24
# @description: This file is used to load the parameters from the YAML file
# =============================================================================

# import libraries ========================
import yaml
import numpy as np

class ParamsLoader:
    def __init__(self, config_path):
        # Load the YAML file
        with open(config_path, 'r') as file:
            self.config = yaml.safe_load(file)['parameters']

        # Initialize the parameters dictionary
        self.icesee_kwargs = self.config.copy()

        # Ensure necessary values are cast to the correct types
        self._cast_parameters()

        # Perform derived calculations
        self._compute_derived_parameters()

        # Generate the grid dictionary
        self._generate_grid()

    def _cast_parameters(self):
        """Ensure specific parameters are of the correct type."""
        self.icesee_kwargs["A"] = float(self.icesee_kwargs["A"])
        self.icesee_kwargs["n"] = int(self.icesee_kwargs["n"])
        self.icesee_kwargs["C"] = float(self.icesee_kwargs["C"])
        self.icesee_kwargs["rho_i"] = float(self.icesee_kwargs["rho_i"])
        self.icesee_kwargs["rho_w"] = float(self.icesee_kwargs["rho_w"])
        self.icesee_kwargs["g"] = float(self.icesee_kwargs["g"])
        self.icesee_kwargs["transient"] = int(self.icesee_kwargs["transient"])
        self.icesee_kwargs["tcurrent"] = int(self.icesee_kwargs["tcurrent"])
        self.icesee_kwargs["accum"] = float(self.icesee_kwargs["accum"]) / self.icesee_kwargs["year"]  # Convert to per second
        self.icesee_kwargs["facemelt"] = float(self.icesee_kwargs["facemelt"]) / self.icesee_kwargs["year"]  # Convert to per second

    def _compute_derived_parameters(self):
        """Compute the derived scaling and other parameters."""
        self.icesee_kwargs["m"] = 1 / self.icesee_kwargs["n"]
        self.icesee_kwargs["B"] = self.icesee_kwargs["A"] ** (-1 / self.icesee_kwargs["n"])

        # Scaling parameters
        self.icesee_kwargs["hscale"] = 1000
        self.icesee_kwargs["ascale"] = 1.0 / self.icesee_kwargs["year"]
        self.icesee_kwargs["uscale"] = (self.icesee_kwargs["rho_i"] * self.icesee_kwargs["g"] * self.icesee_kwargs["hscale"] * self.icesee_kwargs["ascale"] / self.icesee_kwargs["C"]) ** (1 / (self.icesee_kwargs["m"] + 1))
        self.icesee_kwargs["xscale"] = self.icesee_kwargs["uscale"] * self.icesee_kwargs["hscale"] / self.icesee_kwargs["ascale"]
        self.icesee_kwargs["tscale"] = self.icesee_kwargs["xscale"] / self.icesee_kwargs["uscale"]
        self.icesee_kwargs["eps"] = self.icesee_kwargs["B"] * ((self.icesee_kwargs["uscale"] / self.icesee_kwargs["xscale"]) ** (1 / self.icesee_kwargs["n"])) / (2 * self.icesee_kwargs["rho_i"] * self.icesee_kwargs["g"] * self.icesee_kwargs["hscale"])
        self.icesee_kwargs["lambda"] = 1 - (self.icesee_kwargs["rho_i"] / self.icesee_kwargs["rho_w"])

        # Compute NX after ensuring N1 and N2 are cast to integers
        self.icesee_kwargs["N1"] = int(self.icesee_kwargs["N1"])
        self.icesee_kwargs["N2"] = int(self.icesee_kwargs["N2"])
        self.icesee_kwargs["NX"] = self.icesee_kwargs["N1"] + self.icesee_kwargs["N2"]  # Calculate NX

        # Grid time parameters
        self.icesee_kwargs["TF"] = self.icesee_kwargs["year"]  # 1 year in seconds
        self.icesee_kwargs["dt"] = self.icesee_kwargs["TF"] / self.icesee_kwargs["NT"]  # Time step

    def _generate_grid(self):
        """Generate sigma grid values."""
        self.icesee_kwargs["sigGZ"] = float(self.icesee_kwargs["sigGZ"])
        sigma1 = np.linspace(self.icesee_kwargs["sigGZ"] / (self.icesee_kwargs["N1"] + 0.5), self.icesee_kwargs["sigGZ"], int(self.icesee_kwargs["N1"]))
        sigma2 = np.linspace(self.icesee_kwargs["sigGZ"], 1, int(self.icesee_kwargs["N2"] + 1))
        sigma = np.concatenate((sigma1, sigma2[1:self.icesee_kwargs["N2"] + 1]))

        # Create the grid dictionary
        self.icesee_kwargs["grid"] = {
            "sigma": sigma,
            "sigma_elem": np.concatenate(([0], (sigma[:-1] + sigma[1:]) / 2)),
            "dsigma": np.diff(sigma)
        }

    def get_params(self):
        """Return the full parameters dictionary."""
        return self.icesee_kwargs
