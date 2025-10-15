# Scripts

This directory contains utility scripts for data management, visualization, HPC job submission, and MATLAB process management.

## Subdirectories

### data_management/
Scripts for handling ICESEE data files:
- `get_data.py` - Data retrieval utilities
- `inspect_files.py` - File inspection tools
- `post_processing.py` - Post-processing workflows
- `stack_icesee_data.py` - Data stacking and aggregation

### matlab/
MATLAB-related utilities:
- `kill_matlab_processes.py` - Clean up MATLAB processes

### plotting/
Visualization and plotting utilities:
- `scaling_plots.py` - Performance scaling visualizations
- `scaling_plots_csv_details.py` - Detailed CSV-based scaling plots

### slurm/
SLURM job submission scripts for HPC environments:
- `run_da_issm.py` - Run ISSM data assimilation jobs
- `run_job*.sbatch` - Various SLURM batch scripts
- `run_model.m` - MATLAB model execution script

## Usage

These scripts are typically invoked from the command line or imported as utilities in other parts of the ICESEE framework. Refer to individual script documentation for specific usage patterns.
