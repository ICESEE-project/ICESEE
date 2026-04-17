#!/bin/bash
#SBATCH -t40:00:00
#SBATCH -JISMIP
#SBATCH -N 2
#SBATCH -n 48
#SBATCH --ntasks-per-node=24
#SBATCH -p cpu-large                 # pick ONE partition; or switch to: -p inferno
#SBATCH --mem=256G                   # ensure this is allowed per node on that partition
#SBATCH -A gts-arobel3-atlas
#SBATCH -ostrong_scaling-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

# set -euo pipefail
cd $SLURM_SUBMIT_DIR

# Spack env (numpy, mpi4py, h5py, etc.)
module purge || true
source /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ICESEE-Spack/scripts/activate.sh
MPI_DIR="$(spack -e /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ICESEE-Spack/.spack-env/icesee location -i openmpi)"
export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:$LD_LIBRARY_PATH"

# --- Test Matrix ------------------------------------------------------
np=8
NENS=80
MODEL_NPROCS=1

run_case () {
  local np="$1"
  local tag="np${np}_Nens${NENS}_mp${MODEL_NPROCS}"

  /usr/bin/time -v \
    srun --mpi=pmix -n "$np" \
    python run_da_issm.py \
      --Nens=${NENS} \
      --model_nprocs=${MODEL_NPROCS} \
      --verbose

  echo "=== Finished: $tag ==="
}

# call the function with the specified number of processes
run_case "$np"

