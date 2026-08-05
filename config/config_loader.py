# =============================================================================
# src/project_name/config_loader.py
# @author: Brian Kyanjo
# @date: 2025-01-10
# @description: This file is used to load the parameters from the YAML file
# =============================================================================

# Import the required modules
import yaml
import os


def _deep_merge(base, override):
    """Recursively merge a YAML override onto a base configuration."""
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged

def load_yaml_to_dict(file_path, _seen=None):
    """
    Load a YAML file and store its contents in a dictionary.

    Parameters:
        file_path (str): Path to the YAML file.

    Returns:
        dict: A dictionary containing the parsed YAML entries.
    """
    # Check if the file exists before attempting to read it
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Error: The file '{file_path}' was not found.")

    try:
        with open(file_path, "r") as yaml_file:
            # Use safe_load for security and robust parsing
            data = yaml.safe_load(yaml_file) or {}

        base_name = data.pop("extends", None)
        if base_name is None:
            return data

        absolute_path = os.path.abspath(file_path)
        seen = set() if _seen is None else set(_seen)
        if absolute_path in seen:
            raise ValueError(f"Cyclic YAML inheritance detected at '{file_path}'")
        seen.add(absolute_path)

        base_path = base_name
        if not os.path.isabs(base_path):
            base_path = os.path.join(os.path.dirname(absolute_path), base_path)
        base = load_yaml_to_dict(base_path, _seen=seen)
        return _deep_merge(base, data)
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML file '{file_path}': {e}")


def get_section(data, section_name):
    """
    Safely retrieve a section from the loaded YAML dictionary.

    Parameters:
        data (dict): Dictionary containing YAML data.
        section_name (str): Section name to retrieve.

    Returns:
        dict: The requested section as a dictionary or an empty dictionary if not found.
    """
    return data.get(section_name, {})
