#!/bin/bash
#SBATCH --job-name=PIG-reference
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=pig-reference-%j.out
set -euo pipefail
: "${ICESEE_SOURCE_ROOT:?Set ICESEE_SOURCE_ROOT to the validated source checkout}"
PIG="$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig"
PARAMS="$PIG/params.yaml"
if (( $# )); then
  [[ $# == 2 && ( "$1" == -F || "$1" == --force-params ) ]] || { echo 'Usage: sbatch pace/run_reference.sbatch [-F params.yaml]' >&2; exit 2; }
  PARAMS="$2"
fi
[[ "$PARAMS" == /* ]] || PARAMS="$PIG/$PARAMS"
[[ -f "$PARAMS" ]] || { echo "Missing parameters: $PARAMS" >&2; exit 2; }
source "$PIG/pace/resolve_environment.sh" >/dev/null
set -euo pipefail
(( SLURM_NTASKS == 1 && SLURM_JOB_NUM_NODES == 1 && ${SLURM_CPUS_PER_TASK:-1} == 1 )) || { echo 'Reference generation requires one MPI rank/core' >&2; exit 2; }
OUT="${ICESEE_LOG_PATH:-$PIG/pace_runs/reference-${SLURM_JOB_ID:?}}"
mkdir "$OUT"
cp "$PARAMS" "$OUT/params.yaml"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
set +e
srun --mpi=pmix --exact --exclusive --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores \
  "$PIG/pace/launch_rank.sh" "$PIG/tools/generate_reference.py" -F "$OUT/params.yaml" \
  > "$OUT/stdout.log" 2> "$OUT/stderr.log"
RESULT=$?
set -e
printf '%s\n' "$RESULT" > "$OUT/exit-status"
exit "$RESULT"
