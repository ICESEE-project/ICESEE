# ISSM MATLAB Apptainer Image Usage Guide

This README provides instructions for using the Apptainer (Singularity) image built for running the **Ice Sheet System Model (ISSM)** with Python and MATLAB integration. The image is optimized for cluster environments and includes all necessary dependencies, such as MATLAB R2024b, ISSM, and Python libraries with MPI support.

---

## Prerequisites

- **Cluster Environment**: Access to a high-performance computing (HPC) cluster with Apptainer/Singularity installed.
- **Definition File**: The Apptainer definition file `issm_ap.def` used to build the image.
- **Input Data**: Ensure you have the necessary input files (e.g., model data) to populate the examples and execution directories.
- **Slurm Scheduler**: The command uses `srun`, indicating a Slurm-based cluster. Ensure Slurm is configured on your system.
- **MATLAB License**: A valid MATLAB license accessible via the license server (`1711@matlablic.edu`).
- **Build Permissions**: Ability to run `apptainer build`, which may require `fakeroot` or root privileges depending on your cluster configuration.

---

## Setup

### 1. Build the Apptainer Image

Build the image from the definition file using the following command:

```bash
apptainer build issm_matlab.sif issm_matlab_container.def
```

2. Create Directories

Create the directories required for input and output:

```bash
mkdir -p examples execution
```

3. Populate Directories
	•	Place ISSM example data and scripts (e.g., model input files) in the examples/ directory.
	•	The execution/ directory will be used for model outputs (initially empty, but must exist).

---

### Directory Structure

The following directories are bound from the host to the container:
	•	examples/ → /opt/ISSM/execution (container): Contains example data and model scripts.
	•	execution/ → /opt/execution (container): Stores output and execution results.

---

### Running the ISSM Script

Execute the run_da_issm.py script using:

```bash
srun -n 4 apptainer exec \
  -B examples:/opt/ISSM/execution,execution:/opt/execution \
  issm_matlab.sif python run_da_issm.py --Nens=2 --model_nprocs=2
```

### Command Breakdown
	•	srun -n 4: Runs with 4 MPI tasks via Slurm.
	•	apptainer exec: Executes a command inside the container.
	•	-B examples:/opt/ISSM/execution,execution:/opt/execution: Binds host directories into the container.
	•	issm_matlab.sif: The Apptainer image file.
	•	python run_da_issm.py: The Python script to be executed inside the container.
	•	--Nens=2: Number of ensemble members for data assimilation.
	•	--model_nprocs=2: Number of processors for the model simulation.

---

### Expected Behavior
	•	The script run_da_issm.py will execute within the Apptainer container.
	•	Input files are read from /opt/ISSM/examples (host: examples/).
	•	Execution files are written to /opt/execution (host: execution/).
	•	MATLAB is called from Python using the integrated environment within the image.
