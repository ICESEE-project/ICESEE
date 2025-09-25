#!/bin/bash
#SBATCH -t250:00:00
#SBATCH -JISMIP
#SBATCH -N 45                         # 1080 / 24 = 45 nodes
#SBATCH -n 1080                       # request enough tasks for the max np
#SBATCH --ntasks-per-node=24
#SBATCH -p cpu-large                  # pick ONE partition; or switch to: -p inferno
#SBATCH --mem=256G                    # per-node; ensure this is allowed on the partition
#SBATCH -A gts-arobel3-atlas
#SBATCH -oweak_scaling-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

# set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# --- Modules & env ----------------------------------------------------
module purge
module load matlab/r2024b
module load gcc/13

# Python path (ICESEE)
export PYTHONPATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0:$PYTHONPATH

# ISSM root
export ISSM_DIR=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ISSM
export MATLAB_ADDPATH_CMD="addpath('${ISSM_DIR}/bin','${ISSM_DIR}/lib');"

# Your custom OpenMPI (must match mpi4py/h5py builds)
export PATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/openmpi/lib:$LD_LIBRARY_PATH

# Spack env (numpy, mpi4py, h5py, etc.)
source /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/spack_icesee_env.sh

# --- Performance flags ------------------------------------------------
export ICESEE_PERFORMANCE_TEST=1

# --- Weak scaling: np = Nens; model_nprocs fixed ---------------------
MODEL_NPROCS=4
NP_LIST="8 16 32 64 128 256"

OUTROOT="weak_scaling_npEqNens_mp${MODEL_NPROCS}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTROOT"

run_case () {
  local np="$1"
  local nens="$np"   # enforce np = Nens
  local tag="np${np}_Nens${nens}_mp${MODEL_NPROCS}"
  echo "=== Running case: $tag ==="
  mkdir -p "${OUTROOT}/${tag}"

  # Compute total ranks used = np * (MODEL_NPROCS + 1)
  local rank_size=$(( np * (MODEL_NPROCS + 1) ))
  local data_path="_modelrun_datasets_weak_${rank_size}"

  /usr/bin/time -v \
  mpirun \
    -mca coll ^hcoll \
    -np "$np" \
    python3 run_da_issm.py \
      --Nens="${nens}" \
      --model_nprocs="${MODEL_NPROCS}" \
      --verbose \
      --data_path="${data_path}" \
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
done

echo "All runs complete. Results in: ${OUTROOT}"