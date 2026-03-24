#!/bin/bash -f
#SBATCH -J ISSM-ICESEESEESEE           # Job name
#SBATCH --account=gts-arobel3-atlas   # Charge account
#SBATCH -N 1                          # Number of nodes
#SBATCH --ntasks-per-node=24          # Tasks per node (3 nodes * 24 tasks = 72 total tasks)
#SBATCH -t 30:00:00                  # Job duration (120 hours)
#SBATCH -q inferno                    # QOS Name (inferno, embers)
#SBATCH --partition=cpu-large         # Partition name
#SBATCH --mem=10G                    # Total memory per node
#SBATCH -o issm-icsee-%j.out          # Output/error log file
#SBATCH --mail-type=BEGIN,END,FAIL    # Mail preferences
#SBATCH --mail-user=bkyanjo3@gatech.edu # Notification email

# Important variables
icesee_nprocs=1
issm_nprocs=2
Nens=4

# Clear module environment and load required modules
module purge
module load gcc/12
module load mvapich2

# Enable RDMA fork safety
export RDMAV_FORK_SAFE=1

# Define working directories
ICESEE_DIR="$HOME/p-arobel3-0/ICESEE"
CONTAINER_PATH="$ICESEE_DIR/applications/issm_model/issm_utils/containers/issm-icesee.sif"
EXAMPLES_DIR="$ICESEE_DIR/applications/issm_model/examples/ISMIP_Choi/examples"
EXECUTION_DIR="$ICESEE_DIR/applications/issm_model/examples/ISMIP_Choi/execution"

# Check if directories exist before creating them
[ -d "$EXAMPLES_DIR" ] || mkdir -p "$EXAMPLES_DIR"
[ -d "$EXECUTION_DIR" ] || mkdir -p "$EXECUTION_DIR"

# Run ICESEE and ISSM with Apptainer
mpirun -np "$icesee_nprocs" apptainer exec -e \
    --bind "$ICESEE_DIR:/opt/ICESEE" \
    --bind "$EXAMPLES_DIR:/opt/ISSM/examples" \
    --bind "$EXECUTION_DIR:/opt/ISSM/execution" \
    "$CONTAINER_PATH" python3 run_da_issm.py --Nens="$Nens" --model_nprocs="$issm_nprocs"