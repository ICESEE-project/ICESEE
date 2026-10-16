#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -J PIG-DA40
#SBATCH -N 2
#SBATCH -n 40
#SBATCH --ntasks-per-node=20
#SBATCH -p cpu-medium
#SBATCH --mem=256G
#SBATCH -A gts-arobel3-atlas
#SBATCH -o pig_da40_Ne30_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

set -euo pipefail

# Slurm executes a spooled copy: resolve from submission directory or explicit root.
ROOT_HINT="${ICESEE_SPACK_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [[ "$ROOT_HINT" != / && ! -f "$ROOT_HINT/scripts/activate.sh" ]]; do ROOT_HINT="$(dirname "$ROOT_HINT")"; done
HELPER="$ROOT_HINT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
[[ -f "$HELPER" ]] || HELPER="$ROOT_HINT/ICESEE/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
[[ -f "$HELPER" ]] || { echo "Set ICESEE_SPACK_ROOT to the live Spack installation" >&2; exit 2; }
module purge || true
source "$HELPER"
set -euo pipefail

# ----------------------------------------------------------------------
# Run configuration
# ----------------------------------------------------------------------
NENS=20
MODEL_NPROCS=2
NP=$((NENS * MODEL_NPROCS))

NUM_YEARS=5
DT=0.05

OBS_START=1
OBS_MAX=5
FREQ_OBS=1

OUTDIR="pig_da40_Ne${NENS}_${SLURM_JOB_ID}"
mkdir -p "$OUTDIR"

# Preserve the user-facing dataset directory. Config import recreates data_path,
# so refuse an existing target and serialize launchers outside that target.
DATA_DIR="${ICESEE_DATA_PATH:-$PIG_DIR/_modelrun_datasets}"
if [[ -e "$DATA_DIR" || -L "$DATA_DIR" ]]; then
  echo "ERROR: $DATA_DIR already exists; preserve it before a fresh run, or set ICESEE_DATA_PATH." >&2
  exit 2
fi
mkdir -p "$(dirname "$DATA_DIR")"
DATA_LOCK="${DATA_DIR}.launch-lock"
mkdir "$DATA_LOCK" || { echo "ERROR: output target is locked: $DATA_LOCK" >&2; exit 2; }
trap 'rmdir "$DATA_LOCK" 2>/dev/null || true' EXIT


echo "============================================================"
echo "Idealized PIG Mode-3 DA"
echo "============================================================"
echo "Job ID             : $SLURM_JOB_ID"
echo "Nodes              : $SLURM_JOB_NUM_NODES"
echo "MPI ranks          : $NP"
echo "Ensembles          : $NENS"
echo "Ranks/model        : $MODEL_NPROCS"
echo "Years              : $NUM_YEARS"
echo "dt                 : $DT"
echo "Dataset directory  : $DATA_DIR"
echo "Log directory      : $OUTDIR"
echo "============================================================"

export ICESEE_PERFORMANCE_TEST=1
export ICESEE_DEBUG_INIT="${ICESEE_DEBUG_INIT:-1}"
export ICESEE_HANG_DUMP_S="${ICESEE_HANG_DUMP_S:-900}"
TPN=20
(( TPN % MODEL_NPROCS == 0 )) || { echo "Tasks/node must divide into whole model groups" >&2; exit 2; }
NODES=$(((NP + TPN - 1) / TPN))
(( NP % TPN == 0 && NODES <= SLURM_JOB_NUM_NODES && NP <= SLURM_NTASKS )) || {
  echo "Allocation cannot provide the requested node-aligned model groups" >&2; exit 2;
}

# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------

echo
echo "Starting ICESEE at $(date)"
echo

/usr/bin/time -v \
srun --mpi=pmix \
     -N "$NODES" -n "$NP" \
     --ntasks-per-node="$TPN" --distribution=block:block --cpu-bind=cores \
     "$PIG_DIR/pace/launch_rank.sh" \
     "$PIG_DIR/run_da_icepack.py" \
       --Nens="$NENS" \
       --execution_mode=3 \
       --model_nprocs="$MODEL_NPROCS" \
       --num_years="$NUM_YEARS" \
       --dt="$DT" \
       --obs_start_time="$OBS_START" \
       --obs_max_time="$OBS_MAX" \
       --freq_obs="$FREQ_OBS" \
       --data_path="$DATA_DIR" \
     > "${OUTDIR}/stdout.log" \
     2> "${OUTDIR}/stderr.log"

echo
echo "ICESEE completed at $(date)"
echo "Results: $OUTDIR"
