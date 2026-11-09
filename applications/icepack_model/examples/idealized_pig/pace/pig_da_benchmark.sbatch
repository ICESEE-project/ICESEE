#!/bin/bash
# ==============================================================================
# Idealized PIG Mode-3 DA benchmark on PACE Phoenix.
#
# Production timestep (dt=0.05) and physics; one truth trajectory; compact
# initialization (data/extended_beta1000yrs_compact_idx20000.h5); wrong/
# reference trajectory disabled (params.yaml defaults). Only the run length,
# observation window, ensemble size and rank layout are set here.
#
# Layout: NODES x TASKS_PER_NODE ranks, block-distributed, MODEL_NPROCS ranks
# per ensemble member. TASKS_PER_NODE must be a multiple of MODEL_NPROCS so
# no model group spans nodes. Measured on Phoenix (2026-10-03): one member per
# rank (MODEL_NPROCS=1) gives ~1.65x the ensemble throughput per node of
# MODEL_NPROCS=2 for this ~13.6k-node mesh. MODEL_NPROCS=2 is supported;
# icepack refuses model groups that span nodes.
#
# JIT caches (node-local, per job):
#   PyOP2 .so + TSFC kernels : $TMPDIR/icesee-<job>/{pyop2,tsfc}, shared by the
#                               ranks of one node (never NFS, never shared
#                               across nodes or jobs)
#   Loopy/pytools (sqlite)    : $TMPDIR/icesee-<job>/xdg/rank<N>, per rank
# Shared TMPDIR falls back to /tmp. Caches survive for the job lifetime.
# Model groups with more than one rank stay consistent through
# ICESEE/src/utils/pyop2_compile_cache.py (installed by run_da_icepack.py).
#
# Submit from anywhere:  sbatch pace/pig_da_benchmark.sbatch
# Overrides (environment): NENS NUM_YEARS OBS_START OBS_MAX FREQ_OBS
#                          MODEL_NPROCS RUN_ROOT TASKS_PER_NODE
#                          SRUN_EXTRA_ARGS (e.g. --overlap when rehearsing
#                          this script with bash inside an interactive
#                          allocation)
# ==============================================================================
#SBATCH -J PIG-DA-bench
#SBATCH -A gts-arobel3-atlas
#SBATCH -p cpu-medium
#SBATCH -N 2
#SBATCH --ntasks-per-node=15
#SBATCH --mem=96G
#SBATCH -t 08:00:00
#SBATCH -o pig_da_bench_%j.out
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

# --- Experiment (dt fixed at the production value) ---
NENS="${NENS:-30}"
MODEL_NPROCS="${MODEL_NPROCS:-1}"
NUM_YEARS="${NUM_YEARS:-40}"
DT=0.05
OBS_START="${OBS_START:-10}"
OBS_MAX="${OBS_MAX:-40}"
FREQ_OBS="${FREQ_OBS:-1}"
NODES="${SLURM_JOB_NUM_NODES}"
TPN="${TASKS_PER_NODE:-${SLURM_NTASKS_PER_NODE%%(*}}"
NP=$((NODES * TPN))

if (( TPN % MODEL_NPROCS != 0 )); then
  echo "ERROR: tasks per node ($TPN) must be a multiple of MODEL_NPROCS ($MODEL_NPROCS)" >&2
  exit 2
fi
if (( NP != NENS * MODEL_NPROCS )); then
  echo "NOTE: $NP ranks for NENS*MODEL_NPROCS=$((NENS * MODEL_NPROCS)); ICESEE will plan rounds/spares"
fi

# --- Unique run directory ---
RUN_ROOT="${RUN_ROOT:-$PIG_DIR/pace_runs}"
RUN_DIR="$RUN_ROOT/pig_da_bench_Ne${NENS}_y${NUM_YEARS}_${SLURM_JOB_ID}"
mkdir -p "$RUN_DIR"

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


export ICESEE_PERFORMANCE_TEST=1
export ICESEE_DEBUG_INIT="${ICESEE_DEBUG_INIT:-1}"
export ICESEE_HANG_DUMP_S="${ICESEE_HANG_DUMP_S:-900}"
# Rank placement and cache paths are recorded in stdout.log.

# --- Run record ---
{
  echo "dataset_directory : $DATA_DIR"
  echo "job_id            : $SLURM_JOB_ID"
  echo "submitted_from    : ${SLURM_SUBMIT_DIR:-}"
  echo "nodes             : $SLURM_JOB_NODELIST ($NODES)"
  echo "ranks             : $NP ($TPN per node)"
  echo "ensemble members  : $NENS"
  echo "model_nprocs      : $MODEL_NPROCS"
  echo "years / dt        : $NUM_YEARS / $DT"
  echo "observations      : years $OBS_START..$OBS_MAX every $FREQ_OBS"
  echo "initialization    : data/extended_beta1000yrs_compact_idx20000.h5"
  echo "git               : $(git -C "$ICESEE_ROOT" rev-parse HEAD 2>/dev/null) $(git -C "$ICESEE_ROOT" status --porcelain --untracked-files=no 2>/dev/null | wc -l) modified tracked files"
  echo "python            : $(which python)"
  echo "mpirun            : $(which mpirun)"
  echo "start             : $(date -Is)"
} | tee "$RUN_DIR/run_config.txt"
env | grep -E '^(SLURM_|PYTHONPATH|OMP_|ICESEE_|OMPI_MCA|PMIX_MCA|VIRTUAL_ENV)' | sort > "$RUN_DIR/env.txt"
python "$ICESEE_ROOT/scripts/benchmarks/pace_run_manifest.py" --out "$RUN_DIR/manifest.json" \
  > /dev/null 2>&1 || echo "WARNING: manifest generation failed (non-fatal)"

# --- Data assimilation ---
START=$(date +%s)
set +e
srun ${SRUN_EXTRA_ARGS:-} --mpi=pmix -N "$NODES" -n "$NP" --ntasks-per-node="$TPN" --distribution=block:block --cpu-bind=cores \
  "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" \
    --Nens="$NENS" \
    --execution_mode=3 \
    --model_nprocs="$MODEL_NPROCS" \
    --num_years="$NUM_YEARS" \
    --dt="$DT" \
    --obs_start_time="$OBS_START" \
    --obs_max_time="$OBS_MAX" \
    --freq_obs="$FREQ_OBS" \
    --data_path="$DATA_DIR" \
  > "$RUN_DIR/stdout.log" 2> "$RUN_DIR/stderr.log"
STATUS=$?
set -e
ELAPSED=$(( $(date +%s) - START ))


[[ -f "$DATA_DIR/performance.json" ]] && cp "$DATA_DIR/performance.json" "$RUN_DIR/performance.json"
{
  echo "end               : $(date -Is)"
  echo "exit status       : $STATUS"
  echo "srun elapsed [s]  : $ELAPSED"
} | tee -a "$RUN_DIR/run_config.txt"
[[ -f "$RUN_DIR/performance.json" ]] || echo "WARNING: performance.json was not produced"

sacct -j "$SLURM_JOB_ID" --format=JobID,Elapsed,TotalCPU,MaxRSS,NNodes,NTasks,State 2>/dev/null \
  > "$RUN_DIR/sacct.txt" || true
exit "$STATUS"
