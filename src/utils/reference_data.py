"""Path safety for immutable reference data and fresh simulation output."""
from pathlib import Path


def validate_reference_paths(data_path, reference_path):
    if reference_path is None:
        return
    data, reference = Path(data_path).resolve(), Path(reference_path).resolve()
    if data == reference or data in reference.parents or reference in data.parents:
        raise ValueError("reference_data_path and writable data_path must be separate, non-nested directories")
