# SLURM Scripts

This directory contains scripts and batch files for running ICESEE workflows on HPC systems using the SLURM workload manager.

## Batch Scripts

### run_job.sbatch
Basic SLURM job submission script for ICESEE data assimilation runs.

**Usage**:
```bash
sbatch run_job.sbatch
```

### run_job_in.sbatch
Interactive job submission with custom parameters.

### run_job_strong.sbatch
Configuration for strong scaling studies (fixed problem size, increasing processors).

### run_job_weak.sbatch
Configuration for weak scaling studies (problem size scales with processors).

## Python Scripts

### run_da_issm.py
Python script for orchestrating ISSM data assimilation runs on SLURM clusters.

**Features**:
- Job submission automation
- Parameter sweeps
- Dependency management
- Output organization

**Usage**:
```bash
python run_da_issm.py --config <config.yaml>
```

## MATLAB Scripts

### run_model.m
MATLAB script for executing ISSM model runs within SLURM jobs.

## Customization

To adapt these scripts for your HPC system:
1. Update partition/queue names
2. Adjust time limits
3. Modify memory requirements
4. Configure environment modules
5. Set correct paths

## Typical Workflow

1. Prepare configuration files
2. Customize batch script for your system
3. Submit job: `sbatch run_job.sbatch`
4. Monitor: `squeue -u $USER`
5. Check results in output directory

## Environment

These scripts typically require:
- SLURM workload manager
- Environment modules for ICESEE dependencies
- Proper MPI configuration
- Adequate scratch space

## Resources

For cluster-specific information, consult your HPC center's documentation on:
- SLURM configuration
- Queue policies
- File systems
- Module environment
