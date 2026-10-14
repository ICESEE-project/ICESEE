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
#   PyOP2 .so + TSFC kernels : /tmp/icesee-<job>/{pyop2,tsfc}, shared by the
#                               ranks of one node (never NFS, never shared
#                               across nodes or jobs)
#   Loopy/pytools (sqlite)    : /tmp/icesee-<job>/xdg/rank<N>, per rank
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

set -uo pipefail

# --- Locations (independent of the submission directory) ---
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/run_da_icepack.py" ]]; then
  PIG_DIR="$SLURM_SUBMIT_DIR"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/../run_da_icepack.py" ]]; then
  PIG_DIR="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
else
  PIG_DIR="$HOME/r-arobel3-0/ICESEE-Spack/ICESEE/applications/icepack_model/examples/idealized_pig"
fi
ICESEE_ROOT="$(cd "$PIG_DIR/../../../.." && pwd)"          # .../ICESEE
SPACK_REPO="$(cd "$ICESEE_ROOT/.." && pwd)"                 # .../ICESEE-Spack
cd "$PIG_DIR"

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
DATA_DIR="$RUN_DIR/icesee_data"   # ICESEE recreates data_path; keep logs outside it
mkdir -p "$RUN_DIR"

# --- Environment ---
module purge >/dev/null 2>&1 || true
# shellcheck disable=SC1091
source "$SPACK_REPO/scripts/activate.sh"
MPI_DIR="$(spack -e "$SPACK_REPO/.spack-env/icesee" location -i openmpi)"
export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$SPACK_REPO:$ICESEE_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=1
export ICESEE_PERFORMANCE_TEST=1
export ICESEE_DEBUG_INIT="${ICESEE_DEBUG_INIT:-1}"        # rank-labelled init stage timings
export ICESEE_HANG_DUMP_S="${ICESEE_HANG_DUMP_S:-0}"      # set >0 to dump stacks periodically
hash -r

# --- Per-rank launcher: node-local JIT caches and placement record ---
JOB_CACHE="/tmp/icesee-${SLURM_JOB_ID}"
cat > "$RUN_DIR/launch_rank.sh" <<EOF
#!/bin/bash
export PYOP2_CACHE_DIR="$JOB_CACHE/pyop2"
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR="$JOB_CACHE/tsfc"
export XDG_CACHE_HOME="$JOB_CACHE/xdg/rank\${SLURM_PROCID}"
mkdir -p "\$PYOP2_CACHE_DIR" "\$FIREDRAKE_TSFC_KERNEL_CACHE_DIR" "\$XDG_CACHE_HOME"
echo "rank=\${SLURM_PROCID} node=\${SLURM_NODEID} local=\${SLURM_LOCALID} host=\$(hostname -s) cpus=\$(taskset -cp \$\$ 2>/dev/null | awk -F': ' '{print \$2}') group=\$((SLURM_PROCID / $MODEL_NPROCS))" >> "$RUN_DIR/placement.txt"
exec python "\$@"
EOF
chmod +x "$RUN_DIR/launch_rank.sh"
: > "$RUN_DIR/placement.txt"

# --- Run record ---
{
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
srun ${SRUN_EXTRA_ARGS:-} --mpi=pmix -N "$NODES" -n "$NP" --ntasks-per-node="$TPN" --distribution=block --cpu-bind=cores \
  "$RUN_DIR/launch_rank.sh" run_da_icepack.py \
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
ELAPSED=$(( $(date +%s) - START ))

sort -t= -k2 -n "$RUN_DIR/placement.txt" -o "$RUN_DIR/placement.txt"
[[ -f "$DATA_DIR/performance.json" ]] && cp "$DATA_DIR/performance.json" "$RUN_DIR/performance.json"
{
  echo "end               : $(date -Is)"
  echo "exit status       : $STATUS"
  echo "srun elapsed [s]  : $ELAPSED"
} | tee -a "$RUN_DIR/run_config.txt"
[[ -f "$RUN_DIR/performance.json" ]] || echo "WARNING: performance.json was not produced"

# Remove only this job's node-local JIT caches.
srun ${SRUN_EXTRA_ARGS:-} --ntasks-per-node=1 -N "$NODES" -n "$NODES" rm -rf "$JOB_CACHE" >/dev/null 2>&1 || true
sacct -j "$SLURM_JOB_ID" --format=JobID,Elapsed,TotalCPU,MaxRSS,NNodes,NTasks,State 2>/dev/null \
  > "$RUN_DIR/sacct.txt" || true
exit "$STATUS"
