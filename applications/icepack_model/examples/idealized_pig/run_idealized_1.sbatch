#!/bin/bash
#SBATCH --job-name=PIG-DA
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=pig-da-%j.out

set -euo pipefail

RUN_PARAMETERS=""
while (( $# )); do
  case "$1" in
    -F|--force-params)
      (( $# >= 2 )) && [[ -n "$2" && "$2" != -* ]] || { echo "Missing YAML path after $1" >&2; exit 2; }
      [[ -z "$RUN_PARAMETERS" ]] || { echo "Specify -F only once" >&2; exit 2; }
      RUN_PARAMETERS="$2"; shift 2 ;;
    --force-params=*)
      [[ -z "$RUN_PARAMETERS" && -n "${1#*=}" ]] || { echo "Specify one nonempty YAML path" >&2; exit 2; }
      RUN_PARAMETERS="${1#*=}"; shift ;;
    -h|--help)
      echo "Usage: sbatch run_idealized_1.sbatch [-F /path/to/params.yaml]"
      echo "Defaults to this example's params.yaml. Set run parameters in that YAML."
      exit 0 ;;
    *) echo "Unsupported argument: $1. Set run parameters in params.yaml; use -F to select another file." >&2; exit 2 ;;
  esac
done

# Slurm executes a spooled copy; discover source from the submission directory.
if [[ -z "${ICESEE_SOURCE_ROOT:-}" ]]; then
  ROOT_HINT="${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}"
  ROOT_HINT="$(cd "$ROOT_HINT" && pwd -P)"
  while [[ "$ROOT_HINT" != / && ! -f "$ROOT_HINT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh" ]]; do
    ROOT_HINT="$(dirname "$ROOT_HINT")"
  done
  [[ "$ROOT_HINT" != / ]] || { echo "Submit from the ICESEE checkout or set ICESEE_SOURCE_ROOT" >&2; exit 2; }
  export ICESEE_SOURCE_ROOT="$ROOT_HINT"
fi
ICESEE_SOURCE_ROOT="$(cd "$ICESEE_SOURCE_ROOT" && pwd -P)"
export ICESEE_SOURCE_ROOT
LAB_ROOT="$(dirname "$(dirname "$ICESEE_SOURCE_ROOT")")"
if [[ -z "${ICESEE_SPACK_ROOT:-}" ]]; then
  for DEP_ROOT in "$(dirname "$ICESEE_SOURCE_ROOT")" "$ICESEE_SOURCE_ROOT" "$LAB_ROOT/ICESEE-Spack"; do
    if [[ -f "$DEP_ROOT/scripts/activate.sh" && -d "$DEP_ROOT/.spack-env/icesee" ]]; then
      export ICESEE_SPACK_ROOT="$DEP_ROOT"
      break
    fi
  done
fi
: "${ICESEE_SPACK_ROOT:?Set ICESEE_SPACK_ROOT to the installed dependency environment}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONUNBUFFERED=1
unset SLURM_OVERCOMMIT SLURM_OVERLAP
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail

# The local example is authoritative unless -F selects another YAML.
RUN_PARAMETERS="${RUN_PARAMETERS:-$PIG_DIR/params.yaml}"
[[ "$RUN_PARAMETERS" == /* ]] || RUN_PARAMETERS="$PIG_DIR/$RUN_PARAMETERS"
[[ -f "$RUN_PARAMETERS" ]] || { echo "Missing parameter file: $RUN_PARAMETERS" >&2; exit 2; }
for old_key in ICESEE_PARAMS_FILE ICESEE_NENS ICESEE_MODEL_NPROCS ICESEE_DT ICESEE_NUM_YEARS; do
  if [[ -n "${!old_key:-}" ]]; then
    echo "Ignoring $old_key: edit the selected YAML; use -F to select another file." >&2
  fi
done
PREPARED_PARAMS="$(mktemp)"
DATA_LOCK=""
LOCK_HELD=0
cleanup() {
  rm -f -- "$PREPARED_PARAMS"
  if (( LOCK_HELD )); then rmdir "$DATA_LOCK" 2>/dev/null || true; fi
}
trap cleanup EXIT
# Resolve inheritance once without importing the model or its runtime configuration.
RESOLVED="$(python - "$RUN_PARAMETERS" "$ICESEE_ROOT" "$PREPARED_PARAMS" <<'PYTHON'
import math, sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[2]) / "config"))
from config_loader import load_yaml_to_dict
import yaml
cfg = load_yaml_to_dict(sys.argv[1])
model, enkf = cfg["modeling-parameters"], cfg["enkf-parameters"]
def count(value, name, minimum=1):
    number = float(value)
    if not math.isfinite(number) or not number.is_integer() or number < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(number)
mode = count(enkf.get("execution_mode", 1), "execution_mode")
if mode not in (1, 2, 3):
    raise ValueError("This batch script supports execution_mode 1, 2 or 3")
nens = count(enkf["Nens"], "Nens")
spatial = count(enkf.get("model_nprocs", 0), "model_nprocs", 0) or 1
years, dt = float(model["num_years"]), float(model["timesteps_per_year"])
if not all(math.isfinite(v) and v > 0 for v in (years, dt)):
    raise ValueError("num_years and timesteps_per_year (dt) must be finite and positive")
steps = round(years / dt)
if steps < 1:
    raise ValueError("Duration and dt must produce at least one forecast step")
Path(sys.argv[3]).write_text(yaml.safe_dump(cfg, sort_keys=False))
for value in (mode, nens, spatial, years, dt, steps, enkf.get("data_path", "_modelrun_datasets")):
    print(value)
PYTHON
)" || exit 2
VALUES=()
while IFS= read -r value; do VALUES+=("$value"); done <<< "$RESOLVED"
MODE="${VALUES[0]}" NENS="${VALUES[1]}" MODEL_NPROCS="${VALUES[2]}"
NUM_YEARS="${VALUES[3]}" DT="${VALUES[4]}" STEP_COUNT="${VALUES[5]}"
printf 'Source: %s\nDependencies: %s\nParameters: %s\nExecution mode: %s\n' "$ICESEE_ROOT" "$SPACK_REPO" "$RUN_PARAMETERS" "$MODE"
NODES="${SLURM_JOB_NUM_NODES:?Run this script with sbatch}"
NP="${SLURM_NTASKS:?Missing Slurm task count}"
TPN="${ICESEE_TASKS_PER_NODE:-${SLURM_NTASKS_PER_NODE:-8}}"
for value in "$NENS" "$MODEL_NPROCS" "$NODES" "$NP" "$TPN"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "Counts must be positive integers" >&2; exit 2; }
done
(( ${SLURM_CPUS_PER_TASK:-1} == 1 )) || { echo "Use one CPU per MPI rank" >&2; exit 2; }
(( NP == NODES * TPN && (MODE != 3 || TPN % MODEL_NPROCS == 0) )) || {
  echo "Use a uniform tasks/node allocation divisible by the YAML model_nprocs" >&2; exit 2;
}
GROUPS_COUNT=$((NP / MODEL_NPROCS))
(( MODE != 3 || (NP % MODEL_NPROCS == 0 && GROUPS_COUNT <= NENS) )) || { echo "More model groups than ensemble members" >&2; exit 2; }
# Native groups use the allocated MPI ranks without extra coordinators.
if (( MODE != 3 )); then GROUPS_COUNT="n/a"; fi

OUTDIR="$PIG_DIR/pace_runs/mode${MODE}-${SLURM_JOB_ID:?}"
[[ ! -e "$OUTDIR" ]] || { echo "Log directory exists: $OUTDIR" >&2; exit 2; }
DATA_DIR="${ICESEE_DATA_PATH:-${VALUES[6]}}"
[[ "$DATA_DIR" == /* ]] || DATA_DIR="$PIG_DIR/$DATA_DIR"
if [[ -e "$DATA_DIR" || -L "$DATA_DIR" ]]; then
  echo "Preserve existing results first, or set a new ICESEE_DATA_PATH: $DATA_DIR" >&2
  exit 2
fi
mkdir -p "$(dirname "$DATA_DIR")"
DATA_LOCK="${DATA_DIR}.launch-lock"
mkdir "$DATA_LOCK" || { echo "Output is locked: $DATA_LOCK" >&2; exit 2; }
LOCK_HELD=1
# Recheck after acquiring the lock in case another launcher just finished.
[[ ! -e "$DATA_DIR" && ! -L "$DATA_DIR" ]] || { echo "Output appeared: $DATA_DIR" >&2; exit 2; }
mkdir -p "$OUTDIR"
{ git -C "$ICESEE_ROOT" rev-parse HEAD; git -C "$ICESEE_ROOT" status --short; } > "$OUTDIR/source.txt"

cp -- "$RUN_PARAMETERS" "$OUTDIR/params.source.yaml"
python - "$PREPARED_PARAMS" "$OUTDIR/params.yaml" "$DATA_DIR" "$MODEL_NPROCS" "$STEP_COUNT" <<'PYTHON'
import sys, yaml
from pathlib import Path
cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
cfg["enkf-parameters"]["data_path"] = sys.argv[3]
cfg["enkf-parameters"]["model_nprocs"] = int(sys.argv[4])
cfg["modeling-parameters"]["nt"] = int(sys.argv[5])
Path(sys.argv[2]).write_text(yaml.safe_dump(cfg, sort_keys=False))
PYTHON

export ICESEE_DEBUG_INIT="${ICESEE_DEBUG_INIT:-1}"
export ICESEE_HANG_DUMP_S="${ICESEE_HANG_DUMP_S:-900}"
ARGS=(--execution_mode="$MODE" --force-params="$OUTDIR/params.yaml"
      --Nens="$NENS" --model_nprocs="$MODEL_NPROCS"
      --num_years="$NUM_YEARS" --dt="$DT" --data_path="$DATA_DIR")
CMD=(srun --mpi=pmix --exact --exclusive --nodes="$NODES" --ntasks="$NP"
     --ntasks-per-node="$TPN" --cpus-per-task=1 --hint=nomultithread
     --distribution=block:block --cpu-bind=cores
     "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" "${ARGS[@]}")
{
  printf 'Nodes=%s MPI_ranks=%s tasks_per_node=%s members=%s yaml_model_nprocs=%s native_model_groups=%s\n' \
    "$NODES" "$NP" "$TPN" "$NENS" "$MODEL_NPROCS" "$GROUPS_COUNT"
  printf 'Requested_years=%s dt=%s forecast_steps=%s history_columns=%s\n' \
    "$NUM_YEARS" "$DT" "$STEP_COUNT" "$((STEP_COUNT + 1))"
  printf 'Selected_parameters=%s\n' "$RUN_PARAMETERS"
  printf 'Parameters=%s\nDatasets=%s\nLogs=%s\nCommand:' "$OUTDIR/params.yaml" "$DATA_DIR" "$OUTDIR"
  printf ' %q' "${CMD[@]}"
  printf '\n'
} | tee "$OUTDIR/launch.txt"

set +e
/usr/bin/time -v "${CMD[@]}" > "$OUTDIR/stdout.log" 2> "$OUTDIR/stderr.log"
RESULT=$?
set -e
printf '%s\n' "$RESULT" > "$OUTDIR/exit-status"
if (( RESULT != 0 && MODE == 3 )); then
  # MPI has stopped: recover only snapshots with committed manifests.
  set +e
  python "$ICESEE_ROOT/src/parallelization/distributed_history_export.py" \
    "$DATA_DIR" --allow-partial --simulation-exit-code "$RESULT" \
    > "$OUTDIR/finalizer.log" 2>&1
  FINALIZER_RESULT=$?
  set -e
  printf '%s\n' "$FINALIZER_RESULT" > "$OUTDIR/finalizer-exit-status"
  printf 'History recovery exit=%s; log: %s/finalizer.log\n' "$FINALIZER_RESULT" "$OUTDIR"
fi
printf 'Mode %s exit=%s; results: %s; logs: %s\n' "$MODE" "$RESULT" "$DATA_DIR" "$OUTDIR"
exit "$RESULT"
