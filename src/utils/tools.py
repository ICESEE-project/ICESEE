# ==============================================================================
# @des: This file contains helper functions that are used in the main script.
# @date: 2024-10-4
# @author: Brian Kyanjo
# ==============================================================================

import os
import sys
import re
import glob
import time
import subprocess
import h5py
import numpy as np
import logging
import traceback
from mpi4py import MPI
import json, glob, tempfile, hashlib

CKPT_DIRNAME = "_checkpoints"
CKPT_BASENAME = "icesee_ckpt.json"
FNAME_PATTERN = r'icesee_enkf_ens_(\d+)\.h5$'  # matches ..._0000.h5, ..._12.h5, etc.

def _extract_time(fname: str) -> int:
    m = re.search(FNAME_PATTERN, os.path.basename(fname))
    if not m:
        raise ValueError(f"Bad filename (no time index): {fname}")
    return int(m.group(1))

def _list_sorted_files(input_dir: str):
    files = glob.glob(os.path.join(input_dir, "icesee_enkf_ens_*.h5"))
    if not files:
        raise RuntimeError(f"No input files found in {input_dir}")
    files.sort(key=_extract_time)
    return files

def _infer_dataset_name(h5path: str, prefer="states") -> str:
    with h5py.File(h5path, "r") as f:
        if prefer in f:
            return prefer
        keys = [k for k in f.keys() if isinstance(f[k], h5py.Dataset)]
        if len(keys) == 1:
            return keys[0]
        raise RuntimeError(
            f"Could not infer dataset name in {h5path}. "
            f"Found: {keys}. Pass dset_name explicitly."
        )

def h5py_has_mpi():
    return bool(getattr(h5py.get_config(), "mpi", False))

# ---------- Option A: VDS (no data copy) ----------
def build_vds(input_dir: str,
              dset_name: str | None = None,
              out_file: str | None = None,
              fillvalue=np.nan) -> str:
    files = _list_sorted_files(input_dir)
    if dset_name is None:
        dset_name = _infer_dataset_name(files[0], prefer="states")
    if out_file is None:
        out_file = os.path.join(input_dir, "icesee_ensemble_data.h5")

    # Probe meta
    with h5py.File(files[0], "r") as f0:
        nd, nens = f0[dset_name].shape
        dtype = f0[dset_name].dtype

    nt = len(files)

    # Build VDS layout
    layout = h5py.VirtualLayout(shape=(nd, nens, nt), dtype=dtype)
    for t, f in enumerate(files):
        vsrc = h5py.VirtualSource(f, dset_name, shape=(nd, nens))
        layout[:, :, t] = vsrc

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    with h5py.File(out_file, "w", libver="latest") as fout:
        # dset_name='ensemble'
        fout.create_virtual_dataset(dset_name, layout, fillvalue=fillvalue)
        fout.attrs.update({
            "nd": nd, "nens": nens, "nt": nt,
            "stack_type": "VDS",
            "source_dir": os.path.abspath(input_dir),
            "dataset_name": dset_name
        })
        # Compute ensemble mean (iterate time slices lazily)
        mean_dset = fout.create_dataset(
            "ensemble_mean", shape=(nd, nt), dtype=np.float64,
            chunks=(nd, 1), fillvalue=np.nan
        )
        for t in range(nt):
            arr = fout[dset_name][:, :, t]
            mean_dset[:, t] = np.nanmean(arr, axis=1)

    return out_file

# ---------- Option B: materialized 3-D HDF5 ----------
def consolidate_h5(input_dir: str,
                   dset_name: str | None = None,
                   out_file: str | None = None,
                   compression: str = "gzip",
                   compression_opts: int = 4,
                   chunks: tuple[int,int,int] | None = None,
                   allow_missing: bool = False) -> str:
    files = _list_sorted_files(input_dir)
    if dset_name is None:
        dset_name = _infer_dataset_name(files[0], prefer="states")
    if out_file is None:
        out_file = os.path.join(input_dir, "icesee_ensemble_data.h5")

    # Probe meta
    with h5py.File(files[0], "r") as f0:
        nd, nens = f0[dset_name].shape
        dtype = f0[dset_name].dtype

    nt = len(files)
    if chunks is None:
        # Good default for time-wise appends and time-slice reads
        chunks = (nd, nens, 1)

    os.makedirs(os.path.dirname(out_file) or ".", exist_ok=True)
    with h5py.File(out_file, "w") as fout:
        # dset = fout.create_dataset(
        #     dset_name, shape=(nd, nens, nt), dtype=dtype,
        #     chunks=chunks, compression=compression,
        #     compression_opts=compression_opts,
        #     shuffle=True, fletcher32=True
        # )
        # dset_name='ensemble'
        dset = fout.create_dataset(
            'ensemble', shape=(nd, nens, nt), dtype=dtype,
            chunks=chunks, compression=compression,
            compression_opts=compression_opts,
            shuffle=True, fletcher32=True
        )
        mean_dset = fout.create_dataset(
            "ensemble_mean", shape=(nd, nt), dtype=np.float64,
            chunks=(nd, 1), compression=compression,
            compression_opts=compression_opts,
            shuffle=True, fletcher32=True
        )
        fout.attrs.update({
            "nd": nd, "nens": nens, "nt": nt,
            "stack_type": "materialized",
            "source_dir": os.path.abspath(input_dir),
            "dataset_name": dset_name
        })

        # Copy time slices, one file at a time (low memory)
        for t, fpath in enumerate(files):
            try:
                with h5py.File(fpath, "r") as fi:
                    arr = fi[dset_name][...]
            except Exception as e:
                if allow_missing:
                    arr = np.full((nd, nens), np.nan, dtype=dtype)
                else:
                    raise RuntimeError(f"Failed reading {fpath}: {e}") from e

            if arr.shape != (nd, nens):
                raise ValueError(f"Shape mismatch at {fpath}: {arr.shape} != {(nd, nens)}")

            dset[:, :, t] = arr
            mean_dset[:, t] = np.nanmean(arr, axis=1)  # (nd,) → store column

    return out_file

# ---------- Convenience entry point for your pipeline ----------
def finalize_stack(output_dir: str,
                   mode: str = "vds",
                   dset_name: str | None = "states",
                   **kwargs) -> str:
    """
    mode: 'vds' (no copy) or 'h5' (materialized).
    kwargs are passed to the underlying builder (e.g., allow_missing=True).
    """
    if mode.lower() == "vds":
        return build_vds(output_dir, dset_name=dset_name, **kwargs)
    elif mode.lower() in ("h5", "materialized"):
        return consolidate_h5(output_dir, dset_name=dset_name, **kwargs)
    else:
        raise ValueError("mode must be 'vds' or 'h5'")


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
def icesee_get_index(vec=None, **kwargs):
    """
    If var_nd is provided: variables in vec_inputs may have different global sizes.
    In this branch we DO NOT use dim_list, because dim_list is typically packed under
    equal-size assumptions elsewhere in the codebase.

    We compute a deterministic block distribution per variable across ranks:
      - each variable is split (almost) evenly across nranks
      - rank-local ownership is contiguous within each variable
    """
    try:
        var_nd = kwargs.get('var_nd', None)

        if var_nd is not None:
            vec_inputs = kwargs.get("vec_inputs", None)
            params = kwargs.get("params", None)
            if vec_inputs is None or params is None:
                raise ValueError("vec_inputs and params must be provided")

            # communicator selection
            if params["default_run"]:
                comm = kwargs.get("subcomm", None)
            else:
                comm = kwargs.get("comm_world", None)

            # rank/size
            if comm is None or params.get("even_distribution", False):
                rank = 0
                nranks = 1
            else:
                rank = comm.Get_rank()
                nranks = comm.Get_size()

            # var_nd: dict or list aligned with vec_inputs
            if isinstance(var_nd, dict):
                nd_vars = np.array([int(var_nd[v]) for v in vec_inputs], dtype=int)
            else:
                nd_vars = np.asarray(var_nd, dtype=int)
                if nd_vars.shape[0] != len(vec_inputs):
                    raise ValueError("var_nd must be dict keyed by vec_inputs or list aligned with vec_inputs")

            if np.any(nd_vars < 0):
                raise ValueError("var_nd contains negative sizes")

            # base offsets in concatenated global vector
            base_offsets = np.cumsum(np.insert(nd_vars, 0, 0))[:-1]

            def block_decomp(n, p, r):
                """
                Split n items into p contiguous blocks.
                Returns (start, count) for rank r.
                """
                q, rem = divmod(n, p)
                # first 'rem' ranks get q+1
                if r < rem:
                    count = q + 1
                    start = r * (q + 1)
                else:
                    count = q
                    start = rem * (q + 1) + (r - rem) * q
                return start, count

            index_map = {}
            local_size_total = 0

            for i, var in enumerate(vec_inputs):
                n = int(nd_vars[i])
                start_in_var, local_n = block_decomp(n, nranks, rank)

                g0 = int(base_offsets[i] + start_in_var)
                g1 = g0 + int(local_n)

                index_map[var] = np.arange(g0, g1, dtype=int) if local_n > 0 else np.array([], dtype=int)
                local_size_total += int(local_n)

            # Keep return signature compatible:
            # third return is "local size on this rank" (for this new branch, sum of locals)
            return None, index_map, local_size_total

        # ============================
        # Case 2: original equal-size logic (unchanged)
        # ============================
        else:
            vec_inputs = kwargs.get("vec_inputs", None)
            nd = kwargs.get("nd")
            # print(f"[ICESEE-debug] vec_inputs: {vec_inputs}, nd: {nd}, kwargs: {kwargs}\n")
            if kwargs["default_run"]:
                comm = kwargs.get("subcomm", None)
            else:
                comm = kwargs.get("comm_world", None)
            
            # len_vec = kwargs["total_state_param_vars"]
            len_vec = len(vec_inputs)
            dim_list_param = np.array(kwargs.get('dim_list', None)) // len(kwargs.get('vec_inputs_old', None))
            dim_list_param = dim_list_param[:len_vec]
            hdim = nd // len_vec
            # print(f"[ICESEE-debug] len_vec: {len_vec}, dim_list_param: {dim_list_param}, hdim: {hdim}\n")

            if comm is None:
                rank = 0
                dim = dim_list_param[rank]
                offsets = [0]
            else:
                if kwargs["even_distribution"]:
                    rank = 0
                    dim = dim_list_param[rank]
                    offsets = [0]
                else:
                    rank = comm.Get_rank()
                    dim = dim_list_param[rank]
                    offsets = np.cumsum(np.insert(dim_list_param, 0, 0))

            start_idx = offsets[rank]
            index_map = {}
            var_start = 0

            for var in vec_inputs:
                start = var_start + start_idx
                end = start + dim
                index_map[var] = np.arange(start, end)
                var_start += hdim

            local_size_per_rank = kwargs.get('dim_list', None)
            return None, index_map, local_size_per_rank[rank]
    except Exception as e:
        print(f"Error occurred in icesee_get_index: {e}")
        tb_str = "".join(traceback.format_exception(*sys.exc_info()))
        print(f"Traceback details:\n{tb_str}")
        # self.mpi_comm.Abort(1)
    
# ==============================================================================

# # Refined ANSI color codes
# COLORS = {
#     "GRAY": "\033[90m",    # Subtle gray for borders
#     "CYAN": "\033[36m",    # Calm cyan for title
#     "GREEN": "\033[32m",   # Muted green for computational time
#     "MAGENTA": "\033[35m", # Soft magenta for wall-clock time
#     "RESET": "\033[0m"
# }

# def format_time_(seconds: float) -> str:
#     """Convert seconds to a formatted HR:MIN:SEC string with milliseconds."""
#     hours = int(seconds // 3600)
#     minutes = int((seconds % 3600) // 60)
#     secs = int(seconds % 60)
#     millis = int((seconds % 1) * 1000)
#     return f"{hours:02d}:{minutes:02d}:{secs:02d}.{millis:03d}"

# def format_time(seconds: float) -> str:
#     """Convert seconds to a formatted DAY:HR:MIN:SEC string with milliseconds."""
#     days = int(seconds // 86400)  # 86400 seconds in a day
#     hours = int((seconds % 86400) // 3600)
#     minutes = int((seconds % 3600) // 60)
#     secs = int(seconds % 60)
#     millis = int((seconds % 1) * 1000)
#     return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02d}.{millis:03d}"

# def setup_logger(log_file: str = "icesee_timing.log"):
#     """Set up a logger for timing output."""
#     logger = logging.getLogger("ICESEE_Timing")
#     logger.setLevel(logging.INFO)
    
#     # Avoid duplicate handlers
#     if not logger.handlers:
#         # File handler for logging to a file
#         file_handler = logging.FileHandler(log_file)
#         file_handler.setFormatter(logging.Formatter("%(message)s"))
#         logger.addHandler(file_handler)
        
#         # Optional: Stream handler for console output (only for root process)
#         comm = MPI.COMM_WORLD
#         rank = comm.Get_rank()
#         if rank == 0:
#             stream_handler = logging.StreamHandler(sys.stderr)  # Use stderr to avoid stdout issues
#             stream_handler.setFormatter(logging.Formatter("%(message)s"))
#             logger.addHandler(stream_handler)
    
#     return logger

# def display_timing(computational_time: float, wallclock_time: float) -> None:
#     """Display computational and wall-clock times with perfectly aligned formatting using logging."""
#     # Set up logger
#     logger = setup_logger()
    
#     # Only log from the root MPI process
#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     if rank != 0:
#         return  # Non-root processes exit silently

#     # Formatted time strings
#     comp_time_str = format_time(computational_time)
#     wall_time_str = format_time(wallclock_time)
    
#     # Content lines (no trailing spaces after emojis)
#     title = "[ICESEE] Performance Metrics"
#     comp_line = f"Computational Time (Σ): {comp_time_str} (DAY:HR:MIN:SEC.ms) ⏱️"
#     wall_line = f"Wall-Clock Time (max):  {wall_time_str} (DAY:HR:MIN:SEC.ms) 🕒"
    
#     # Calculate max width based on plain text length (excluding ANSI codes)
#     max_content_width = max(len(title), len(comp_line), len(wall_line))
#     box_width = max_content_width + 12  # 2 for '║' on each side + 2 for padding
    
#     # Box drawing
#     header = f"{COLORS['GRAY']}╔{'═' * box_width}╗{COLORS['RESET']}"
#     footer = f"{COLORS['GRAY']}╚{'═' * box_width}╝{COLORS['RESET']}"
    
#     # Pad lines to exact width, ensuring no extra spaces
#     def pad_line(text: str) -> str:
#         padding = " " * (max_content_width - len(text) + 6 + 4)
#         return f"{COLORS['GRAY']}║ {text}{padding} ║{COLORS['RESET']}"
    
#     def pad_line_comp(text: str) -> str:
#         padding = " " * (max_content_width - len(text) + 7 + 4)
#         return f"{COLORS['GRAY']}║ {text}{padding} ║{COLORS['RESET']}"
    
#     def pad_line_wall(text: str) -> str:
#         padding = " " * (max_content_width - len(text) + 5 + 4)
#         return f"{COLORS['GRAY']}║ {text}{padding} ║{COLORS['RESET']}"
    
#     # Log with strict alignment
#     logger.info(f"\n{header}")
#     logger.info(f"{COLORS['CYAN']}{pad_line(title)}{COLORS['RESET']}")
#     logger.info(f"{COLORS['GREEN']}{pad_line_comp(comp_line)}{COLORS['RESET']}")
#     logger.info(f"{COLORS['MAGENTA']}{pad_line_wall(wall_line)}{COLORS['RESET']}")
#     logger.info(footer)


# Refined ANSI color codes
COLORS = {
    "GRAY": "\033[10m",    # Uniform gray for all text and borders
    "RESET": "\033[0m"
}

def format_time(seconds: float) -> str:
    """Convert seconds to a formatted DAY:HR:MIN:SEC string with milliseconds."""
    days = int(seconds // 86400)
    hours = int((seconds % 86400) // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    millis = int((seconds % 1) * 1000)
    return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02d}.{millis:03d}"

def setup_logger(log_file: str = "icesee_timing.log"):
    """Set up a logger for timing output."""
    import logging
    import sys
    from mpi4py import MPI
    
    logger = logging.getLogger("ICESEE_Timing")
    logger.setLevel(logging.INFO)
    
    if not logger.handlers:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(file_handler)
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            stream_handler = logging.StreamHandler(sys.stderr)
            stream_handler.setFormatter(logging.Formatter("%(message)s"))
            logger.addHandler(stream_handler)
    
    return logger

def display_timing_default(computational_time: float, wallclock_time: float) -> None:
    """Display computational and wall-clock times with perfectly aligned formatting using logging."""
    # Set up logger
    logger = setup_logger()
    
    # Only log from the root MPI process
    # comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank != 0:
        return

    # Formatted time strings
    comp_time_str = format_time(computational_time)
    wall_time_str = format_time(wallclock_time)
    
    # Content lines (no trailing spaces after emojis)
    # title = "[ICESEE] Performance Metrics"
    title = f"[ICESEE] Metrics on {MPI.COMM_WORLD.Get_size()} ranks"
    comp_line = f"Computational Time (Σ): {comp_time_str} (DAY:HR:MIN:SEC.ms) ⏱️"
    wall_line = f"Wall-Clock Time (max):  {wall_time_str} (DAY:HR:MIN:SEC.ms) 🕒"
    
    # Calculate max width based on the longest metric label and value
    max_label_width = max(len(entry[0]) for entry in time_entries)
    max_value_width = max(len(entry[1]) for entry in time_entries[1:])  # Skip header for value width
    total_width = max_label_width + max_value_width - 14  # 2 for '║' + 2 for padding
    
    # Box drawing
    header = f"{COLORS['GRAY']}╔{'═' * total_width}╗{COLORS['RESET']}"
    footer = f"{COLORS['GRAY']}╚{'═' * total_width}╝{COLORS['RESET']}"
    
    # Pad lines to exact width with strict alignment
    def pad_line(label: str, value: str = "") -> str:
        if not value:  # Header
            padding = " " * (total_width -10 - len(label))
            return f"{COLORS['GRAY']}║ \033[1m{label}{COLORS['RESET']}{padding}{COLORS['GRAY']}║{COLORS['RESET']}"
        else:  # Metric with value
            label_padding = " " * (max_label_width -17 - len(label))  # +1 for space
            value_padding = " " * (max_value_width -17 - len(value))  # +1 for space
            return f"{COLORS['GRAY']}║ {label}{label_padding}{value}{value_padding}{COLORS['RESET']}{COLORS['GRAY']}  ║{COLORS['RESET']}"
    
    # Log with strict alignment
    logger.info(f"{header}")
    for entry in time_entries:
        if len(entry) == 1:  # Header
            logger.info(pad_line(entry[0]))
        else:  # Metric with value
            logger.info(pad_line(entry[0], entry[1]))
    logger.info(footer)


# Refined ANSI color codes
COLORS = {
    "GRAY": "\033[10m",    # Uniform gray for all text and borders
    "RESET": "\033[0m"
}

def format_time(seconds: float) -> str:
    """Convert seconds to a formatted DAY:HR:MIN:SEC string with milliseconds."""
    days = int(seconds // 86400)
    hours = int((seconds % 86400) // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    millis = int((seconds % 1) * 1000)
    return f"{days:02d}:{hours:02d}:{minutes:02d}:{secs:02d}.{millis:03d}"

def setup_logger(log_file: str = "icesee_timing.log"):
    """Set up a logger for timing output."""
    import logging
    import sys
    from mpi4py import MPI
    
    logger = logging.getLogger("ICESEE_Timing")
    logger.setLevel(logging.INFO)
    
    if not logger.handlers:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(file_handler)
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        if rank == 0:
            stream_handler = logging.StreamHandler(sys.stderr)
            stream_handler.setFormatter(logging.Formatter("%(message)s"))
            logger.addHandler(stream_handler)
    
    return logger

def display_timing_verbose(
    computational_time: float,
    wallclock_time: float,
    true_wrong_time: float,
    assimilation_time: float,
    forecast_step_time: float,
    analysis_step_time: float,
    ensemble_init_time: float,
    init_file_time: float,
    forecast_file_time: float,
    analysis_file_time: float,
    total_file_time: float,
    forecast_noise_time: float,
    time_init_ensemble_mean_computation: float,
    time_forecast_ensemble_mean_computation: float,
    time_analysis_ensemble_mean_computation: float,
    comm: MPI.Comm = None,
    model_nprocs: int = 1
) -> None:
    """Display all timing metrics in a table with strict aligned formatting using logging, all in gray."""
    # from mpi4py import MPI
    
    # Set up logger
    logger = setup_logger()
    
    # Only log from the root MPI process
    # comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank != 0:
        return

    # Formatted time strings with metrics and values
    time_entries = [
        (f"[ICESEE] Performance Metrics ({MPI.COMM_WORLD.Get_size()*(model_nprocs+1)} ranks)     (DAY:HR:MIN:SEC.ms)",),  # Bold header
        ("Computational Time (Σ)", format_time(computational_time)),
        ("Wall-Clock Time (max)", format_time(wallclock_time)),
        ("True/Wrong State Time", format_time(true_wrong_time)),
        ("Ensemble Init Time", format_time(ensemble_init_time)),
        ("Forecast Step Time", format_time(forecast_step_time)),
        ("Analysis Step Time", format_time(analysis_step_time)),
        ("Assimilation Time", format_time(assimilation_time)),
        ("Init file I/O Time", format_time(init_file_time)),
        ("Forecast File I/O Time", format_time(forecast_file_time)),
        ("Analysis File I/O Time", format_time(analysis_file_time)),
        ("Total File I/O Time", format_time(total_file_time)),
        ("Forecast Noise Time", format_time(forecast_noise_time)),
        ("Init Ensemble Mean Computation", format_time(time_init_ensemble_mean_computation)),
        ("Forecast Ensemble Mean Computation", format_time(time_forecast_ensemble_mean_computation)),
        ("Analysis Ensemble Mean Computation", format_time(time_analysis_ensemble_mean_computation)),
    ]
    
    # Calculate max width based on the longest metric label and value
    max_label_width = max(len(entry[0]) for entry in time_entries)
    max_value_width = max(len(entry[1]) for entry in time_entries[1:])  # Skip header for value width
    total_width = max_label_width + max_value_width - 14  # 2 for '║' + 2 for padding
    
    # Box drawing
    header = f"{COLORS['GRAY']}╔{'═' * total_width}╗{COLORS['RESET']}"
    footer = f"{COLORS['GRAY']}╚{'═' * total_width}╝{COLORS['RESET']}"
    
    # Pad lines to exact width with strict alignment
    def pad_line(label: str, value: str = "") -> str:
        if not value:  # Header
            padding = " " * (total_width -10 - len(label))
            return f"{COLORS['GRAY']}║ \033[1m{label}{COLORS['RESET']}{padding}{COLORS['GRAY']}║{COLORS['RESET']}"
        else:  # Metric with value
            label_padding = " " * (max_label_width -17 - len(label))  # +1 for space
            value_padding = " " * (max_value_width -17 - len(value))  # +1 for space
            return f"{COLORS['GRAY']}║ {label}{label_padding}{value}{value_padding}{COLORS['RESET']}{COLORS['GRAY']}  ║{COLORS['RESET']}"
    
    # Log with strict alignment
    logger.info(f"{header}")
    for entry in time_entries:
        if len(entry) == 1:  # Header
            logger.info(pad_line(entry[0]))
        else:  # Metric with value
            logger.info(pad_line(entry[0], entry[1]))
    logger.info(footer)

def get_grid_dimensions(nx, ny, ndim):
    """
    Calculate grid dimensions mx and my based on physical dimensions and total points.
    
    Parameters:
    nx (int): Number of elements in x-direction
    ny (int): Number of elements in y-direction
    ndim (int): Total number of grid points (mx * my)
    
    Returns:
    tuple: (mx, my) - number of grid points in x and y directions
    """
    # Calculate aspect ratio from physical dimensions
    alpha = nx / ny
    
    # Initial estimate based on aspect ratio and ndim
    # mx/my = alpha and mx*my = ndim
    # mx = sqrt(ndim * alpha), my = sqrt(ndim / alpha)
    mx = np.sqrt(ndim * alpha)
    my = np.sqrt(ndim / alpha)
    
    # Initial rounding
    if mx - int(mx) > 0.5:
        mx = int(np.ceil(mx))
        my = int(np.floor(my))
    elif my - int(my) > 0.5:
        my = int(np.ceil(my))
        mx = int(np.floor(mx))
    else:
        mx, my = int(mx), int(my)
    
    # Quick adjustment to reach ndim
    current_product = mx * my
    if current_product != ndim:
        # Calculate scale factor
        scale = np.sqrt(ndim / current_product)
        mx = int(round(mx * scale))
        my = int(round(my * scale))
        
        # Fast fine-tuning with minimal iterations
        product = mx * my
        if product < ndim:
            while product < ndim:
                if mx/my < alpha:
                    mx += 1
                else:
                    my += 1
                product = mx * my
        elif product > ndim:
            while product > ndim:
                if mx/my > alpha:
                    mx -= 1
                else:
                    my -= 1
                product = mx * my
    
    return mx, my

def midpoint_rect(mx, my):
    return mx/2.0, my/2.0

def midprofiles_coords(mx, my, n=100):
    x_mid = np.full(n, mx/2.0)
    yv = np.linspace(0.0, my, n)         # vertical profile (x fixed)
    y_mid = np.full(n, my/2.0)
    xv = np.linspace(0.0, mx, n)         # horizontal profile (y fixed)
    return (x_mid, yv), (xv, y_mid)

def midindices(Nx, Ny):
    # choose the "left/bottom" center for even sizes; adjust if you prefer the right/top
    ix = (Nx-1)//2
    iy = (Ny-1)//2
    return ix, iy

def _extract_time_from_name(fname: str) -> int:
    m = re.search(FNAME_PATTERN, os.path.basename(fname))
    if not m:
        raise ValueError(f"Bad filename (no time index): {fname}")
    return int(m.group(1))

def _sorted_step_files(base_dir: str) -> list[str]:
    files = glob.glob(os.path.join(base_dir, "icesee_enkf_ens_*.h5"))
    files.sort(key=_extract_time_from_name)
    return files

def _last_completed_step(base_dir: str) -> int | None:
    """Return last completed time index (int) or None if none exist."""
    files = _sorted_step_files(base_dir)
    if not files:
        return None
    return _extract_time_from_name(files[-1])

def _ckpt_path(base_dir: str) -> str:
    return os.path.join(base_dir, CKPT_DIRNAME, CKPT_BASENAME)

def _atomic_write_json(path: str, payload: dict):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    # atomic via temp + rename
    fd, tmppath = tempfile.mkstemp(prefix=".tmp_ckpt_", dir=os.path.dirname(path))
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(payload, f, indent=2)
            f.flush(); os.fsync(f.fileno())
        os.replace(tmppath, path)
    except Exception:
        try: os.remove(tmppath)
        except Exception: pass
        raise

def save_checkpoint(base_dir: str, **state):
    """Rank 0 only: persist minimal restart info."""
    _atomic_write_json(_ckpt_path(base_dir), state)

def load_checkpoint(base_dir: str) -> dict | None:
    path = _ckpt_path(base_dir)
    if not os.path.exists(path):
        return None
    with open(path, "r") as f:
        return json.load(f)

def compute_km_from_tobserve(tobserve, k_start, m_obs=None):
    import numpy as np

    # coerce tobserve to a flat int64 1-D array
    tobserve = np.asarray(tobserve).astype(np.int64, copy=False).ravel()

    # coerce m_obs to an int (or default to len(tobserve))
    if m_obs is None:
        m_obs_i = tobserve.size
    else:
        try:
            m_obs_i = int(m_obs)  # handles python int or numpy scalar
        except Exception:
            m_obs_i = int(np.asarray(m_obs).reshape(-1)[0])
    # clamp to valid range
    m_obs_i = max(0, min(m_obs_i, tobserve.size))

    # count how many obs times have occurred at start (remember your check uses k+1)
    # k1 = int(k_start) + 1
    k1 = int(k_start)
    return int(np.count_nonzero(tobserve[:m_obs_i] <= k1))

def step_already_done(base_dir: str, k: int) -> bool:
    # accept both zero-padded and plain
    p1 = os.path.join(base_dir, f"icesee_enkf_ens_{k:04d}.h5")
    p2 = os.path.join(base_dir, f"icesee_enkf_ens_{k}.h5")
    return os.path.exists(p1) or os.path.exists(p2)

def reseed_for_step(base_seed: int, rank_world: int, k: int):
    """Deterministic per-step seeding (optional, keeps behavior reproducible on restart)."""
    # Any scheme is fine as long as it's stable. This mixes step + rank.
    seed = (base_seed * 1315423911 + 2654435761 * (k + 1) + rank_world) % (2**31 - 1)
    np.random.seed(seed)
    return seed

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Stack ICESEE DA slices into (nd, nens, nt)")
    ap.add_argument("output_dir", help="Directory containing icesee_enkf_ens_*.h5")
    ap.add_argument("--mode", choices=["vds", "h5"], default="vds",
                    help="vds = virtual dataset (no copy), h5 = materialized 3D file")
    ap.add_argument("--dset-name", default="states",
                    help="Dataset name inside each file (default: states)")
    ap.add_argument("--out-file", default=None, help="Path of output file")
    ap.add_argument("--allow-missing", action="store_true",
                    help="Fill missing/unreadable slices with NaN for materialized mode")
    ap.add_argument("--compression", default="gzip", help="HDF5 compression (h5 mode)")
    ap.add_argument("--compression-opts", type=int, default=4, help="Compression level (h5 mode)")
    ap.add_argument("--chunk-nd", type=int, default=None)
    ap.add_argument("--chunk-nens", type=int, default=None)
    ap.add_argument("--chunk-nt", type=int, default=1)
    args = ap.parse_args()

    chunks = None
    if args.chunk_nd or args.chunk_nens or args.chunk_nt:
        chunks = (args.chunk_nd, args.chunk_nens, args.chunk_nt)

    out = finalize_stack(
        args.output_dir, mode=args.mode, dset_name=args.dset_name,
        out_file=args.out_file, compression=args.compression,
        compression_opts=args.compression_opts, chunks=chunks,
        allow_missing=args.allow_missing
    )
    print(f"[Finalize] Stacked dataset written: {out}")

def icesee_fingerprint(params: dict, keys=("model_name","nd","nt","Nens","base_seed")) -> str:
    sub = {k: params.get(k) for k in keys}
    blob = json.dumps(sub, sort_keys=True, separators=(",",":"))
    return hashlib.sha1(blob.encode("utf-8")).hexdigest()

def h5_has_dataset_with_shape(path: str, dset: str, shape: tuple[int,...]) -> bool:
    if not os.path.exists(path):
        return False
    try:
        with h5py.File(path, "r") as f:
            if dset not in f: return False
            return tuple(f[dset].shape) == tuple(shape)
    except Exception:
        return False

def h5_attr_equals(path: str, attr: str, expected: str) -> bool:
    try:
        with h5py.File(path, "r") as f:
            return str(f.attrs.get(attr, "")) == str(expected)
    except Exception:
        return False

def mark_h5_with_fingerprint(path: str, attr="icesee_fingerprint", value: str | None = None, extra: dict | None = None):
    with h5py.File(path, "a") as f:
        if value is not None:
            f.attrs[attr] = value
        if extra:
            for k,v in extra.items():
                f.attrs[k] = v

def env_flag(name: str, default: bool = False) -> bool:
    """Interpret environment variable flags like 0/1, true/false, on/off."""
    val = os.environ.get(name)
    if val is None:
        return default
    val = str(val).strip().lower()
    if val in ("1", "true", "yes", "on"):
        return True
    if val in ("0", "false", "no", "off"):
        return False
    # fallback: any non-empty string means True
    return True


def load_bed_masks_from_h5(f):
    """
    Supports BOTH:
      - new format: /bed_masks/static/* and /bed_masks/cols/*
      - old format: /bed_mask_map (legacy)

    Returns:
      bed_mask_map_static: dict[str, np.ndarray(bool)]  shape (n_bed,)
      bed_mask_map_cols:   dict[str, np.ndarray(bool)]  shape (n_bed, m_obs)
      bed_snap_cols:       list[int]
      obs_model_to_col:    dict[int,int]
    """
    bed_mask_map_static = {}
    bed_mask_map_cols = {}
    bed_snap_cols = []
    obs_model_to_col = {}

    # ---------- NEW FORMAT ----------
    if "bed_masks" in f:
        # static masks
        if "static" in f["bed_masks"]:
            for k in f["bed_masks/static"].keys():
                bed_mask_map_static[k] = f["bed_masks/static"][k][:].astype(bool)

        # column/time-dependent masks
        if "cols" in f["bed_masks"]:
            for k in f["bed_masks/cols"].keys():
                bed_mask_map_cols[k] = f["bed_masks/cols"][k][:].astype(bool)

        if "bed_snap_cols" in f:
            bed_snap_cols = f["bed_snap_cols"][:].astype(int).tolist()

        # mapping model step -> obs column
        if "obs_model_to_col_keys" in f and "obs_model_to_col_vals" in f:
            keys = f["obs_model_to_col_keys"][:].astype(int)
            vals = f["obs_model_to_col_vals"][:].astype(int)
            obs_model_to_col = {int(k): int(v) for k, v in zip(keys, vals)}

        return bed_mask_map_static, bed_mask_map_cols, bed_snap_cols, obs_model_to_col

    # ---------- LEGACY FORMAT ----------
    if "bed_mask_map" in f:
        legacy = f["bed_mask_map"][:]
        # We can only interpret legacy as "static" (no per-km gating available)
        # If legacy stored (m_obs,) or something else, this is inherently ambiguous.
        bed_mask_map_static["bed"] = np.asarray(legacy, dtype=bool)
        # no cols gating
        bed_mask_map_cols = {}
        bed_snap_cols = []
        obs_model_to_col = {}

    return bed_mask_map_static, bed_mask_map_cols, bed_snap_cols, obs_model_to_col



def icesee_savefig(fig, name="results.png", dpi=300, show=True):
    """
    ICESEE-OnLINE helper:
    Always save plots into ./figures/ so the GUI can display them.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to save.
    name : str
        Output filename inside figures/.
    dpi : int
        Resolution.
    show : bool
        Whether to call plt.show() after saving.
    """

    from pathlib import Path
    import matplotlib.pyplot as plt
    # Ensure figures folder exists
    Path("figures").mkdir(exist_ok=True)

    # Full output path
    out_path = Path("figures") / name

    # Save figure
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")

    print(f"[ICESEE] Figure saved: {out_path}")

    # Optionally display inline
    if show:
        plt.show()