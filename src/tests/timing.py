# test_script.py
from mpi4py import MPI
import time

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
    time_forecast_ensemble_mean_generation: float,
    time_analysis_ensemble_mean_generation: float,
    comm: MPI.Comm = None
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
        ("[ICESEE] Performance Metrics       (DAY:HR:MIN:SEC.ms)",),  # Bold header
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
        ("Forecast Ensemble Mean Computation", format_time(time_forecast_ensemble_mean_generation)),
        ("Analysis Ensemble Mean Computation", format_time(time_analysis_ensemble_mean_generation)),
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

# Test script
if __name__ == "__main__":
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Simulate timing data
    if rank == 0:
        # Sample timing values in seconds
        computational_time = 60100.156  # ~16 hours, 10 minutes, 156 ms
        wallclock_time = 0.058          # ~58 ms
        true_wrong_time = 0.0           # Zero for testing
        ensemble_init_time = 0.0        # Zero for testing
        forecast_step_time = 0.009      # ~9 ms
        analysis_step_time = 0.021      # ~21 ms
        assimilation_time = 0.03        # Sum of above steps
        init_file_time = 0.008          # ~8 ms
        forecast_file_time = 0.008      # ~8 ms
        analysis_file_time = 0.0        # Zero for testing
        total_file_time = 0.017         # Sum of file times
        forecast_noise_time = 0.003     # ~3 ms

        # Call display_timing with simulated data
        display_timing_verbose(
            computational_time=computational_time,
            wallclock_time=wallclock_time,
            true_wrong_time=true_wrong_time,
            assimilation_time=assimilation_time,
            forecast_step_time=forecast_step_time,
            analysis_step_time=analysis_step_time,
            ensemble_init_time=ensemble_init_time,
            init_file_time=init_file_time,
            forecast_file_time=forecast_file_time,
            analysis_file_time=analysis_file_time,
            total_file_time=total_file_time,
            forecast_noise_time=forecast_noise_time,
            time_init_ensemble_mean_computation=0.005,
            time_forecast_ensemble_mean_generation=0.007,
            time_analysis_ensemble_mean_generation=0.006,
            comm=comm
        )

    # Synchronize processes
    comm.Barrier()