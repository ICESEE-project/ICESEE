# ==============================================================================
# @des: This file contains helper functions that are used in the main script.
# @date: 2024-10-4
# @author: Brian Kyanjo
# ==============================================================================

import os
import sys
import re
import time
import subprocess
import h5py
import numpy as np
import logging
from mpi4py import MPI

# Function to safely change directory
def safe_chdir(main_directory,target_directory):
    # Get the absolute path of the target directory
    target_path = os.path.abspath(target_directory)

    # Check if the target path starts with the main directory path
    if target_path.startswith(main_directory):
        os.chdir(target_directory)
    # else:
    #     print(f"[ICESEE] Error: Attempted to leave the main directory '{main_directory}'.")


def install_requirements(force_install=False, verbose=False):
    """
    Install dependencies listed in the requirements.txt file if not already installed,
    or if `force_install` is set to True.
    """
    # Check if the `.installed` file exists to determine if installation is needed
    if os.path.exists(".installed") and not force_install:
        print("[ICESEE] Dependencies are already installed. Skipping installation.")
        return
    
    try:
        # Run the command to install the requirements from requirements.txt
        print("[ICESEE] Installing dependencies from requirements.txt...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "../requirements.txt"])
        
        # Create a `.installed` marker file to indicate successful installation
        with open(".installed", "w") as f:
            f.write("Dependencies installed successfully.\n")

        print("[ICESEE] All dependencies are installed and verified.")
    except subprocess.CalledProcessError as e:
        # Print the error and raise a more meaningful exception
        print(f"[ICESEE] Error occurred while installing dependencies: {e}")
        raise RuntimeError("Failed to install dependencies from requirements.txt. Please check the file and try again.")

# ==== saves arrays to h5 file
def save_arrays_to_h5(filter_type=None, model=None, parallel_flag=None, commandlinerun=None, **datasets):
    """
    Save multiple arrays to an HDF5 file, optionally in a parallel environment (MPI).

    Parameters:
        filter_type (str): Type of filter used (e.g., 'ENEnKF', 'DEnKF').
        model (str): Name of the model (e.g., 'icepack').
        parallel_flag (str): Flag to indicate if MPI parallelism is enabled. Default is 'MPI'.
        commandlinerun (bool): Indicates if the function is triggered by a command-line run. Default is False.
        **datasets (dict): Keyword arguments where keys are dataset names and values are arrays to save.

    Returns:
        dict: The datasets if not running in parallel, else None.
    """
    output_dir = "results"
    output_file = f"{output_dir}/{filter_type}-{model}.h5"

    if parallel_flag == "MPI" or commandlinerun:
        # Create the results folder if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print("[ICESEE] Creating results folder")

        # Remove the existing file, if any
        if os.path.exists(output_file):
            os.remove(output_file)
            print(f"[ICESEE] Existing file {output_file} removed.")

        print(f"[ICESEE] Writing data to {output_file}")
        with h5py.File(output_file, "w") as f:
            for name, data in datasets.items():
                f.create_dataset(name, data=data, compression="gzip")
                print(f"[ICESEE] Dataset '{name}' written to file")
        print(f"[ICESEE] Data successfully written to {output_file}")
    else:
        print("[ICESEE] Non-MPI or non-commandline run. Returning datasets.")
        return datasets

# Routine extracts datasets from a .h5 file
def extract_datasets_from_h5(file_path):
    """
    Extracts all datasets from an HDF5 file and returns them as a dictionary.

    Parameters:
        file_path (str): Path to the HDF5 file.

    Returns:
        dict: A dictionary where keys are dataset names and values are numpy arrays.

    Raises:
        FileNotFoundError: If the specified HDF5 file does not exist.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    datasets = {}
    print(f"[ICESEE] Reading data from {file_path}...")

    with h5py.File(file_path, "r") as f:
        def extract_group(group, datasets):
            for key in group.keys():
                item = group[key]
                if isinstance(item, h5py.Dataset):
                    datasets[key] = np.array(item)
                    print(f"[ICESEE] Dataset '{key}' extracted with shape {item.shape}")
                elif isinstance(item, h5py.Group):
                    extract_group(item, datasets)

        extract_group(f, datasets)

    print("[ICESEE] Data extraction complete.")
    return datasets

# --- best for saving all data to h5 file in parallel environment
def save_all_data(enkf_params=None, nofilter=None, **kwargs):
    """
    General function to save datasets based on the provided parameters.
    """
    # Update filter_type only if nofilter is provided
    filter_type = "true-wrong" if nofilter else enkf_params["filter_type"]

    # --- Local MPI implementation ---
    if re.match(r"\AMPI\Z", enkf_params["parallel_flag"], re.IGNORECASE) or re.match(r"\AMPI_model\Z", enkf_params["parallel_flag"], re.IGNORECASE):
        from mpi4py import MPI
        comm = MPI.COMM_WORLD  # Initialize MPI
        rank = comm.Get_rank()  # Get rank of current MPI process
        size = comm.Get_size()  # Get total number of MPI processes

        comm.Barrier()
        if rank == 0:
            save_arrays_to_h5(
                filter_type=filter_type,  # Use updated or original filter_type
                model=enkf_params["model_name"],
                parallel_flag=enkf_params["parallel_flag"],
                commandlinerun=enkf_params["commandlinerun"],
                **kwargs
            )
        else:
            None
    else:
        save_arrays_to_h5(
            filter_type=filter_type,  # Use updated or original filter_type
            model=enkf_params["model_name"],
            parallel_flag=enkf_params["parallel_flag"],
            commandlinerun=enkf_params["commandlinerun"],
            **kwargs
        )

# ---- function to get the index of the variables in the vector dynamically
def icesee_get_index(vec, **kwargs):
    """
    Get the index of the variables in the vector dynamically.

    Parameters:
        - vec: The vector to be distributed
        - vec_inputs: List of input variable names, e.g., ['h', 'u', 'v', 'smb']
        - hdim: Size of each variable in vec_inputs
        - dim_list: List of sizes of each rank; ordered through the ranks using MPI_gather on root and broadcast
        - comm: MPI communicator containing the rank and size of the processors

    Returns:
        - var_indices: Dictionary where keys are variable names and values are their respective slices from `vec`
        - index_map: Dictionary where keys are variable names and values are the indices corresponding to their slices
    """
    # -- get parameters
    vec_inputs = kwargs.get("vec_inputs", None)
    params = kwargs.get("params", None)
    if params["default_run"]:
        comm = kwargs.get("subcomm", None)
    else:
        comm = kwargs.get("comm_world", None)
    
    # get size of input vector based on user inputs
    len_vec = params["total_state_param_vars"]

    # print(f"[ICESEE] dim_list: {kwargs['dim_list']}")
    dim_list_param = np.array(kwargs.get('dim_list', None)) // len_vec  # Get the size of each variable slice
    hdim = vec.shape[0] // len_vec  # Compute the size of each variable in vec_inputs

    if comm is None:
        # Non-MPI case
        rank = 0
        dim = dim_list_param[rank]
        offsets = [0]  # No offsets needed
    else:
        # MPI case
        size_world = kwargs.get("comm_world").Get_size()  # Get the total number of processors
        # if params["even_distribution"] or (params["default_run"] and size_world <= params["Nens"]):
        if params["even_distribution"]:
            rank = 0 # Set rank to 0 for even distribution
            dim = dim_list_param[rank]
            offsets = [0]
        else:
            rank = comm.Get_rank()  # Get the rank of the current processor
            dim = dim_list_param[rank]
            offsets = np.cumsum(np.insert(dim_list_param, 0, 0))  # Compute offsets per processor

    start_idx = offsets[rank]  # Get the start index of the current processor
   
    # Dynamically determine start indices for each variable
    var_indices = {}
    index_map = {}
    var_start = 0  # Initial start index

    for var in vec_inputs:
        start = var_start + start_idx
        end = start + dim
        var_indices[var] = vec[start:end]
        index_map[var] = np.arange(start, end)  # Store index range for easy fetching
        var_start += hdim  # Move to the next variable slice

    local_size_per_rank = kwargs.get('dim_list', None)
    return var_indices, index_map, local_size_per_rank[rank]
# ==============================================================================

# Refined ANSI color codes
COLORS = {
    "GRAY": "\033[90m",    # Subtle gray for borders
    "CYAN": "\033[36m",    # Calm cyan for title
    "GREEN": "\033[32m",   # Muted green for computational time
    "MAGENTA": "\033[35m", # Soft magenta for wall-clock time
    "RESET": "\033[0m"
}

def format_time_(seconds: float) -> str:
    """Convert seconds to a formatted HR:MIN:SEC string with milliseconds."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    millis = int((seconds % 1) * 1000)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}.{millis:03d}"

def format_time(seconds: float) -> str:
    """Convert seconds to a formatted DAY:HR:MIN:SEC string with milliseconds."""
    days = int(seconds // 86400)  # 86400 seconds in a day
    hours = int((seconds % 86400) // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    millis = int((seconds % 1) * 1000)
    return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02d}.{millis:03d}"

def setup_logger(log_file: str = "icesee_timing.log"):
    """Set up a logger for timing output."""
    logger = logging.getLogger("ICESEE_Timing")
    logger.setLevel(logging.INFO)
    
    # Avoid duplicate handlers
    if not logger.handlers:
        # File handler for logging to a file
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(file_handler)
        
        # Optional: Stream handler for console output (only for root process)
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            stream_handler = logging.StreamHandler(sys.stderr)  # Use stderr to avoid stdout issues
            stream_handler.setFormatter(logging.Formatter("%(message)s"))
            logger.addHandler(stream_handler)
    
    return logger

def display_timing(computational_time: float, wallclock_time: float) -> None:
    """Display computational and wall-clock times with perfectly aligned formatting using logging."""
    # Set up logger
    logger = setup_logger()
    
    # Only log from the root MPI process
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank != 0:
        return  # Non-root processes exit silently

    # Formatted time strings
    comp_time_str = format_time(computational_time)
    wall_time_str = format_time(wallclock_time)
    
    # Content lines (no trailing spaces after emojis)
    title = "[ICESEE] Performance Metrics"
    comp_line = f"Computational Time (Σ): {comp_time_str} (DAY:HR:MIN:SEC.ms) ⏱️"
    wall_line = f"Wall-Clock Time (max):  {wall_time_str} (DAY:HR:MIN:SEC.ms) 🕒"
    
    # Calculate max width based on plain text length (excluding ANSI codes)
    max_content_width = max(len(title), len(comp_line), len(wall_line))
    box_width = max_content_width + 12  # 2 for '║' on each side + 2 for padding
    
    # Box drawing
    header = f"{COLORS['GRAY']}╔{'═' * box_width}╗{COLORS['RESET']}"
    footer = f"{COLORS['GRAY']}╚{'═' * box_width}╝{COLORS['RESET']}"
    
    # Pad lines to exact width, ensuring no extra spaces
    def pad_line(text: str) -> str:
        padding = " " * (max_content_width - len(text) + 6 + 4)
        return f"{COLORS['GRAY']}║ {text}{padding} ║{COLORS['RESET']}"
    
    def pad_line_comp(text: str) -> str:
        padding = " " * (max_content_width - len(text) + 7 + 4)
        return f"{COLORS['GRAY']}║ {text}{padding} ║{COLORS['RESET']}"
    
    def pad_line_wall(text: str) -> str:
        padding = " " * (max_content_width - len(text) + 5 + 4)
        return f"{COLORS['GRAY']}║ {text}{padding} ║{COLORS['RESET']}"
    
    # Log with strict alignment
    logger.info(f"\n{header}")
    logger.info(f"{COLORS['CYAN']}{pad_line(title)}{COLORS['RESET']}")
    logger.info(f"{COLORS['GREEN']}{pad_line_comp(comp_line)}{COLORS['RESET']}")
    logger.info(f"{COLORS['MAGENTA']}{pad_line_wall(wall_line)}{COLORS['RESET']}")
    logger.info(footer)