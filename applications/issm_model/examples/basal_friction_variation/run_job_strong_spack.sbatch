#!/bin/bash
#SBATCH -t250:00:00
#SBATCH -JISMIP
#SBATCH -N 11
#SBATCH -n 264
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

# --- Thread pinning ---------------------------------------------------
export ICESEE_PERFORMANCE_TEST=1

# --- Test Matrix ------------------------------------------------------
NENS=64
MODEL_NPROCS=4
# NP_LIST="64 32 16 8 4 2 1"
NP_LIST="8 64 32 16 4 2 1"

OUTROOT="strong_scaling_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTROOT"

 # "$MPI_DIR/bin/mpirun"  -np "$np" \
 #  srun --mpi=pmix -n "$np" \

run_case () {
  local np="$1"
  local tag="np${np}_Nens${NENS}_mp${MODEL_NPROCS}"
  echo "=== Running case: $tag ==="
  mkdir -p "${OUTROOT}/${tag}"

  /usr/bin/time -v \
    srun --mpi=pmix -n "$np" \
    python run_da_issm.py \
      --Nens=${NENS} \
      --model_nprocs=${MODEL_NPROCS} \
      --verbose \
    > "${OUTROOT}/${tag}/stdout.log" \
    2> "${OUTROOT}/${tag}/stderr.log"

  echo "=== Finished: $tag ==="
}

for np in $NP_LIST; do
  (( np > NENS )) && { echo "Skip np=${np} (NENS=${NENS})"; continue; }
  run_case "$np"
done

echo "All runs complete. Results in: ${OUTROOT}"
