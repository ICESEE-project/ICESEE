# =============================================================================
# @author: Brian Kyanjo
# @date: 2025-01-17
# @description: Add (register) or remove models to be supported by the application here.
#               Currently supported models include:
#               - icepack
#               - Lorenz96
#               - flowline (integration still underway)
#               - ISSM    (development still underway)
# =============================================================================

# --- Imports ---
import sys
import os
import re
import importlib

class SupportedModels:
    """
    Class to call the supported models in the application.
    Easily add/remove models by updating MODEL_CONFIG dictionary
    """

    # Dictionary mapping model names to their respective import paths and states
    MODEL_CONFIG = {
        "icepack": {
            "module": "icepack_model.icepack_utils._icepack_enkf",
            "description": "Icepack model",
            "status": "supported",
        },
        "issm": {
            "module": "issm_model.issm_utils._issm_enkf",
            "description": "ISSM model",
            "status": "supported",
        },
        "lorenz96": {
            "module": "lorenz_model.lorenz_utils._lorenz_enkf",
            "description": "Lorenz96 model",
            "status": "supported",
        },
        "flowline": {
            "module": "flowline_model.flowline_utils._flowline_enkf",
            "description": "Flowline model",
            "status": "development",
        },
    }

    def __init__(self, model=None, model_config=None,comm=None,verbose=False):
        self.model = model
        self.MODEL_CONFIG = model_config or self.MODEL_CONFIG
        self.comm = comm
        self.verbose = verbose

        self.rank = self.comm.Get_rank() if self.comm else 0

    def _get_project_root(self):
        """Automatically determines the root of the project."""
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # Traverse upwards until we reach the root of the project (assuming 'src' folder exists at root)
        while not os.path.exists(os.path.join(current_dir, 'src')):
            current_dir = os.path.dirname(current_dir)
        return current_dir
    
    # Globally add application directory to the path
    project_root = _get_project_root(None)
    application_dir = os.path.join(project_root, 'applications')
    sys.path.insert(0, application_dir)

    def list_models(self):
        """List all supported models with their descriptions and statuses."""
        for model, info in self.MODEL_CONFIG.items():
            if self.rank == 0:
                if self.verbose:
                    print(f"[ICESEE] {model.capitalize()}: {info['description']} (Status: {info['status']})")

    def call_model(self):
        """
        Dynamically import and return the modules for the specified model.
        Tries to load the _<model>_enkf module from the current working directory's example folder
        (e.g., <model>_model.examples.<dir_name>._<model>_enkf), falling back to the default module.
        """
        if not self.model:
            raise ValueError("No model specified. Please provide a model name.")

        # Normalize model name for case-insensitive matching
        normalized_model = self.model.lower()

        if normalized_model not in self.MODEL_CONFIG:
            raise ValueError(f"Model '{self.model}' is not supported or implemented.")

        model_info = self.MODEL_CONFIG[normalized_model]

        # Check model status
        if model_info["status"] != "supported":
            raise ValueError(f"Model '{self.model}' is still under development: {model_info['description']}.")

        # Dynamic discovery based on current directory
        current_dir = os.getcwd()
        # Extract the directory name (e.g., 'synthetic' or 'shallowIce')
        dir_name = os.path.basename(current_dir)
        # Construct the dynamic module path: <model>_model.examples.<dir_name>._<model>_enkf
        dynamic_module_path = f"{normalized_model}_model.examples.{dir_name}._{normalized_model}_enkf"

        # Validate directory name to ensure it's a valid Python module name
        if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', dir_name):
            print(f"Invalid directory name '{dir_name}' for module import. Falling back to default module.")
        else:
            try:
                model_module = importlib.import_module(dynamic_module_path)
                if self.rank == 0:
                    if self.verbose:
                        # Print a message indicating the successful loading of the dynamic module
                        print(f"[ICESEE] Successfully loaded example-specific {model_info['description']} from {dynamic_module_path}.")
                return model_module
            except ImportError as e:
                if self.rank == 0:
                    if self.verbose:
                        # Print a message indicating the failure to load the dynamic module
                        print(f"[ICESEE] Dynamic example module {dynamic_module_path} not found: {e}. Falling back to default module.")

        # Fall back to the default module
        try:
            model_module = importlib.import_module(model_info["module"])
            if self.rank == 0:
                if self.verbose:
                    print(f"[ICESEE] Successfully loaded {model_info['description']} from {model_info['module']}.")
            
            return model_module
        except ImportError as e:
            raise ImportError(f"Failed to import module for model '{self.model}': {e}")