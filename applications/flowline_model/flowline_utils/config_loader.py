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
        self.params = self.config.copy()
        
        # Ensure necessary values are cast to the correct types
        self._cast_parameters()

        # Perform derived calculations
        self._compute_derived_parameters()

        # Generate the grid dictionary
        self._generate_grid()

    def _cast_parameters(self):
        """Ensure specific parameters are of the correct type."""
        self.params["A"] = float(self.params["A"])
        self.params["n"] = int(self.params["n"])
        self.params["C"] = float(self.params["C"])
        self.params["rho_i"] = float(self.params["rho_i"])
        self.params["rho_w"] = float(self.params["rho_w"])
        self.params["g"] = float(self.params["g"])
        self.params["transient"] = int(self.params["transient"])
        self.params["tcurrent"] = int(self.params["tcurrent"])
        self.params["accum"] = float(self.params["accum"]) / self.params["year"]  # Convert to per second
        self.params["facemelt"] = float(self.params["facemelt"]) / self.params["year"]  # Convert to per second

    def _compute_derived_parameters(self):
        """Compute the derived scaling and other parameters."""
        self.params["m"] = 1 / self.params["n"]
        self.params["B"] = self.params["A"] ** (-1 / self.params["n"])

        # Scaling parameters
        self.params["hscale"] = 1000
        self.params["ascale"] = 1.0 / self.params["year"]
        self.params["uscale"] = (self.params["rho_i"] * self.params["g"] * self.params["hscale"] * self.params["ascale"] / self.params["C"]) ** (1 / (self.params["m"] + 1))
        self.params["xscale"] = self.params["uscale"] * self.params["hscale"] / self.params["ascale"]
        self.params["tscale"] = self.params["xscale"] / self.params["uscale"]
        self.params["eps"] = self.params["B"] * ((self.params["uscale"] / self.params["xscale"]) ** (1 / self.params["n"])) / (2 * self.params["rho_i"] * self.params["g"] * self.params["hscale"])
        self.params["lambda"] = 1 - (self.params["rho_i"] / self.params["rho_w"])

        # Compute NX after ensuring N1 and N2 are cast to integers
        self.params["N1"] = int(self.params["N1"])
        self.params["N2"] = int(self.params["N2"])
        self.params["NX"] = self.params["N1"] + self.params["N2"]  # Calculate NX

        # Grid time parameters
        self.params["TF"] = self.params["year"]  # 1 year in seconds
        self.params["dt"] = self.params["TF"] / self.params["NT"]  # Time step

    def _generate_grid(self):
        """Generate sigma grid values."""
        self.params["sigGZ"] = float(self.params["sigGZ"])
        sigma1 = np.linspace(self.params["sigGZ"] / (self.params["N1"] + 0.5), self.params["sigGZ"], int(self.params["N1"]))
        sigma2 = np.linspace(self.params["sigGZ"], 1, int(self.params["N2"] + 1))
        sigma = np.concatenate((sigma1, sigma2[1:self.params["N2"] + 1]))

        # Create the grid dictionary
        self.params["grid"] = {
            "sigma": sigma,
            "sigma_elem": np.concatenate(([0], (sigma[:-1] + sigma[1:]) / 2)),
            "dsigma": np.diff(sigma)
        }

    def get_params(self):
        """Return the full parameters dictionary."""
        return self.params
