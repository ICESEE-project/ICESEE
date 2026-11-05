#!/bin/bash
#SBATCH --job-name=PIG-mode3
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=pig-mode3-%j.out

set -euo pipefail

# Select scientific source explicitly; Slurm executes a spooled script copy.
: "${ICESEE_SOURCE_ROOT:?Set ICESEE_SOURCE_ROOT to the updated ICESEE checkout}"
: "${ICESEE_SPACK_ROOT:?Set ICESEE_SPACK_ROOT to the installed dependency environment}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONUNBUFFERED=1
unset SLURM_OVERCOMMIT SLURM_OVERLAP
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail

# Defaults match the successful 32-member, five-year Mode-3 validation.
NENS="${ICESEE_NENS:-32}"
MODEL_NPROCS="${ICESEE_MODEL_NPROCS:-2}"
NODES="${SLURM_JOB_NUM_NODES:?Run this script with sbatch}"
NP="${SLURM_NTASKS:?Missing Slurm task count}"
TPN="${ICESEE_TASKS_PER_NODE:-${SLURM_NTASKS_PER_NODE:-8}}"
for value in "$NENS" "$MODEL_NPROCS" "$NODES" "$NP" "$TPN"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "Counts must be positive integers" >&2; exit 2; }
done
(( ${SLURM_CPUS_PER_TASK:-1} == 1 )) || { echo "Use one CPU per MPI rank" >&2; exit 2; }
(( NP == NODES * TPN && TPN % MODEL_NPROCS == 0 )) || {
  echo "Use a uniform tasks/node allocation divisible by ICESEE_MODEL_NPROCS" >&2; exit 2;
}
GROUPS_COUNT=$((NP / MODEL_NPROCS))
(( GROUPS_COUNT <= NENS )) || { echo "More model groups than ensemble members" >&2; exit 2; }
# Firedrake model groups are the MPI ranks; there is no extra coordinator per member.

OUTDIR="$PIG_DIR/pace_runs/mode3-${SLURM_JOB_ID:?}"
[[ ! -e "$OUTDIR" ]] || { echo "Log directory exists: $OUTDIR" >&2; exit 2; }
DATA_DIR="${ICESEE_DATA_PATH:-$PIG_DIR/_modelrun_datasets}"
[[ "$DATA_DIR" == /* ]] || DATA_DIR="$PIG_DIR/$DATA_DIR"
if [[ -e "$DATA_DIR" || -L "$DATA_DIR" ]]; then
  echo "Preserve existing results first, or set a new ICESEE_DATA_PATH: $DATA_DIR" >&2
  exit 2
fi
mkdir -p "$(dirname "$DATA_DIR")"
DATA_LOCK="${DATA_DIR}.launch-lock"
mkdir "$DATA_LOCK" || { echo "Output is locked: $DATA_LOCK" >&2; exit 2; }
trap 'rmdir "$DATA_LOCK" 2>/dev/null || true' EXIT
# Recheck after acquiring the lock in case another launcher just finished.
[[ ! -e "$DATA_DIR" && ! -L "$DATA_DIR" ]] || { echo "Output appeared: $DATA_DIR" >&2; exit 2; }
mkdir -p "$OUTDIR"
{ git -C "$ICESEE_ROOT" rev-parse HEAD; git -C "$ICESEE_ROOT" status --short; } > "$OUTDIR/source.txt"

export ICESEE_DEBUG_INIT="${ICESEE_DEBUG_INIT:-1}"
export ICESEE_HANG_DUMP_S="${ICESEE_HANG_DUMP_S:-900}"
ARGS=(--execution_mode=3 --mode3_experimental=True --Nens="$NENS"
      --model_nprocs="$MODEL_NPROCS" --num_years=5 --dt=0.05
      --obs_start_time=2 --obs_max_time=4 --freq_obs=1
      '--save_steps=[2,3,4]' --data_path="$DATA_DIR")
CMD=(srun --mpi=pmix --exact --exclusive --nodes="$NODES" --ntasks="$NP"
     --ntasks-per-node="$TPN" --cpus-per-task=1 --hint=nomultithread
     --distribution=block:block --cpu-bind=cores
     "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" "${ARGS[@]}")
{
  printf 'Nodes=%s MPI_ranks=%s tasks_per_node=%s members=%s spatial_ranks=%s model_groups=%s\n' \
    "$NODES" "$NP" "$TPN" "$NENS" "$MODEL_NPROCS" "$GROUPS_COUNT"
  printf 'Datasets=%s\nLogs=%s\nCommand:' "$DATA_DIR" "$OUTDIR"
  printf ' %q' "${CMD[@]}"
  printf '\n'
} | tee "$OUTDIR/launch.txt"

set +e
/usr/bin/time -v "${CMD[@]}" > "$OUTDIR/stdout.log" 2> "$OUTDIR/stderr.log"
RESULT=$?
set -e
printf '%s\n' "$RESULT" > "$OUTDIR/exit-status"
printf 'Mode 3 exit=%s; results: %s; logs: %s\n' "$RESULT" "$DATA_DIR" "$OUTDIR"
exit "$RESULT"
