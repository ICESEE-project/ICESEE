#!/bin/bash

# --- Find repo root (looks for scripts/activate.sh) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
while [[ "$SCRIPT_DIR" != "/" && ! -f "$SCRIPT_DIR/scripts/activate.sh" ]]; do
  SCRIPT_DIR="$(dirname "$SCRIPT_DIR")"
done
REPO_ROOT="$SCRIPT_DIR"

# --- Activate environment ---
module purge || true
source "$REPO_ROOT/scripts/activate.sh"

# --- Fix MPI environment ---
MPI_DIR="$(spack location -i openmpi)"
export PATH="$(echo "$PATH" | tr ':' '\n' | grep -Ev 'openmpi|mvapich|mpich' | paste -sd: -)"
export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:$LD_LIBRARY_PATH"
hash -r

# --- Run ICESEE ---
mpirun -np 8 python3 run_da_issm.py --Nens=80 --model_nprocs=1
