#!/bin/bash
#SBATCH -t250:00:00
#SBATCH -JISMIP
#SBATCH -N 14                         # 1080 / 24 = 45 nodes
#SBATCH -n 336                        # request enough tasks for the max np
#SBATCH --ntasks-per-node=24
#SBATCH -p cpu-large                  # pick ONE partition; or switch to: -p inferno
#SBATCH --mem=256G                    # per-node; ensure this is allowed on the partition
#SBATCH -A gts-arobel3-atlas
#SBATCH -oweak_scaling-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bankyanjo@gmail.com

# set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# Spack env (numpy, mpi4py, h5py, etc.)
module purge || true
source /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ICESEE-Spack/scripts/activate.sh
MPI_DIR="$(spack -e /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ICESEE-Spack/.spack-env/icesee location -i openmpi)"
export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:$LD_LIBRARY_PATH"

# --- Performance flags ------------------------------------------------
export ICESEE_PERFORMANCE_TEST=1

# --- Weak scaling: np = Nens; model_nprocs fixed ---------------------
MODEL_NPROCS=4
NP_LIST="4 8 16 32 64 128 256"

OUTROOT="weak_scaling_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTROOT"

run_case () {
  local np="$1"
  local nens="$np"   # enforce np = Nens
  local tag="np${np}_Nens${nens}_mp${MODEL_NPROCS}"
  echo "=== Running case: $tag ==="
  mkdir -p "${OUTROOT}/${tag}"

  # Compute total ranks used = np * (MODEL_NPROCS + 1)
  local rank_size=$(( np * (MODEL_NPROCS + 1) ))
  local data_path="_modelrun_datasets"

  /usr/bin/time -v \
    "$MPI_DIR/bin/mpirun"  -np "$np" \
      python run_da_issm.py \
      --Nens="${nens}" \
      --model_nprocs="${MODEL_NPROCS}" \
      --verbose \
    > "${OUTROOT}/${tag}/stdout.log" \
    2> "${OUTROOT}/${tag}/stderr.log"

  echo "=== Finished: $tag ==="
}

# Safety: skip if requested np exceeds SLURM allocation
for np in $NP_LIST; do
  if (( np > SLURM_NTASKS )); then
    echo "Skip np=${np} (exceeds allocation SLURM_NTASKS=${SLURM_NTASKS})"
    continue
  fi
  run_case "$np"
  # pause for 2 minutes between runs to allow for cleanup and resource availability
  sleep 120
done

echo "All runs complete. Results in: ${OUTROOT}"
