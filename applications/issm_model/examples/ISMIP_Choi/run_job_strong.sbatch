#!/bin/sh
#SBATCH -t150:00:00
#SBATCH -JISMIP
#SBATCH -N 11
#SBATCH -n 264
#SBATCH --ntasks-per-node=24
#SBATCH --partition=cpu-large         # Partition name
#SBATCH --mem=256G                     # Total memory per node
#SBATCH -pinferno
#SBATCH -A gts-arobel3-atlas
#SBATCH -ostrong_scaling-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# --- Modules & env ----------------------------------------------------
module purge
module load matlab/r2024b
module load gcc/13

# Python path (ICESEE)
export PYTHONPATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0:$PYTHONPATH

# ISSM root
export ISSM_DIR=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ISSM

alias matlab='matlab -r "addpath ${ISSM_DIR}/bin ${ISSM_DIR}/lib"'

# MPI from your custom OpenMPI toolchain (keep consistent with Python build)
export LD_LIBRARY_PATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/openmpi/lib:$LD_LIBRARY_PATH
export PATH=/storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/openmpi/bin:$PATH

# Your spack/env activations (numpy, mpi4py, h5py, etc.)
source /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/spack_icesee_env.sh

# --- Thread pinning to avoid oversubscription -------------------------
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export ICESEE_PERFORMANCE_TEST=1

# --- Test Matrix ------------------------------------------------------
NENS=64
MODEL_NPROCS=4                     # ISSM side; used inside MATLAB via mpiexec -np 4
NP_LIST="64 32 16 8 4 2 1"         # ICESEE MPI ranks (top-level parallelism)

# Output root
OUTROOT="runs_perf_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTROOT"

# --- Helper: run one case ---------------------------------------------
run_case () {
  local np="$1"
  local tag="np${np}_Nens${NENS}_mp${MODEL_NPROCS}"

  echo "=== Running case: $tag ==="
  mkdir -p "${OUTROOT}/${tag}"
 
  /usr/bin/time -v \
  mpirun -mca coll ^hcoll "$np" \
       bash -lc "set -euo pipefail; \
         echo [\$(hostname)] starting ICESEE with np=${np}, Nens=${NENS}, model_nprocs=${MODEL_NPROCS}; \
         python3 run_da_issm.py \
            --Nens=${NENS} \
            --model_nprocs=${MODEL_NPROCS} \
            --verbose" \
       >  "${OUTROOT}/${tag}/stdout.log" \
       2> "${OUTROOT}/${tag}/stderr.log"

  echo "=== Finished: $tag ==="
}

# --- Cold vs warm passes (optional) -----------------------------------
# Cold: fresh MATLAB servers (default behavior)
for np in $NP_LIST; do
  # Respect Nens ≥ np
  if (( np > NENS )); then
    echo "Skip np=${np} (NENS=${NENS})"
    continue
  fi
  run_case "$np"
done

# Warm: if your Python↔MATLAB layer supports persistent servers, uncomment:
# export ICESEE_REUSE_MATLAB=1
# for np in $NP_LIST; do
#   (( np > NENS )) && continue
#   run_case "$np"
# done

echo "All runs complete. Results in: ${OUTROOT}"