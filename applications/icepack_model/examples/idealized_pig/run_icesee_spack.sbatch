#!/bin/sh
#SBATCH -t 50:00:00
#SBATCH -J ISMIP
#SBATCH -N 2
#SBATCH -n 48
#SBATCH --ntasks-per-node=24
#SBATCH --partition=cpu-large
#SBATCH --mem=128G
#SBATCH -p inferno
#SBATCH -A gts-arobel3-atlas
#SBATCH -o steadystate-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bankyanjo@gmail.com

cd $SLURM_SUBMIT_DIR

# --- User parameters ---
NENS=32
#MODEL_NPROCS=1
NP=16   # match -n above

echo "SLURM JOB INFO"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks: $SLURM_NTASKS"

# --- Call run script ---
#./run_icesee_spack.sh $NENS $MODEL_NPROCS $NP


#NENS="${1:-80}"
#NPROCS="${2:-1}"
#NP="${3:-${SLURM_NTASKS:-8}}"

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

srun -n "$NP" python3 run_da_icepack.py \
  --Nens="$NENS" 
  #--model_nprocs="$NPROCS"