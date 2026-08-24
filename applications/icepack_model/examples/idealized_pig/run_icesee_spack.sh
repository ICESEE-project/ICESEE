#!/bin/bash

NENS="${1:-80}"
NPROCS="${2:-1}"
NP="${3:-${SLURM_NTASKS:-8}}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
while [[ "$SCRIPT_DIR" != "/" && ! -f "$SCRIPT_DIR/scripts/activate.sh" ]]; do
  SCRIPT_DIR="$(dirname "$SCRIPT_DIR")"
done
REPO_ROOT="$SCRIPT_DIR"

module purge || true
source "$REPO_ROOT/scripts/activate.sh"

MPI_DIR="$(spack location -i openmpi)"
export PATH="$(echo "$PATH" | tr ':' '\n' | grep -Ev 'openmpi|mvapich|mpich' | paste -sd: -)"
export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:$LD_LIBRARY_PATH"
hash -r

echo "NENS      = $NENS"
echo "NPROCS    = $NPROCS"
echo "NP        = $NP"
echo "python3   = $(which python3)"
echo "srun      = $(which srun)"
echo "mpirun    = $(which mpirun)"
echo "prte      = $(which prte 2>/dev/null || echo not-found)"

srun --mpi=pmix -n "$NP" python3 run_da_icepack.py \
  --Nens="$NENS" 
  #--model_nprocs="$NPROCS"