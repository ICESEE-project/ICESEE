#!/bin/bash
# Short native PIG calibration using params_calibration.yaml.
#SBATCH -J icesee-run1-pig-calib
#SBATCH -A gts-arobel3-atlas          # ACCOUNT
#SBATCH -p inferno                    # PARTITION/QOS
#SBATCH --partition=cpu-small
#SBATCH -t 00:30:00                   # WALLTIME -- generous for a 40-step, Ne=4 calibration
#SBATCH -N 1
#SBATCH -n 4                          # NP == world_size below (ensemble_groups * model_nprocs)
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH -o run1-calib-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=CHANGE_ME@gatech.edu

set -euo pipefail

# Resolve the example directory from this script path.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IDEALIZED_PIG_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$IDEALIZED_PIG_DIR"

# --- Calibration target (item B) ---
NENS=4
MODEL_NPROCS=2       # P_model
# ensemble_groups is implied: world_size / model_nprocs = 4/2 = 2 groups,
# so Ne=4 over 2 groups -> rounds=2 (2 members per rank, sequential).

echo "=== RUN1 startup diagnostics ==="
echo "hostname: $(hostname)"
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unset}"
echo "working directory: $IDEALIZED_PIG_DIR"
echo "git branch: $(git -C "$IDEALIZED_PIG_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
echo "git revision: $(git -C "$IDEALIZED_PIG_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"

# --- Environment (same pattern as RUN0 / ../run_icesee_spack.sh) ---
# See run0_env_smoke.sbatch's comment: ICESEE-Spack's activate.sh is a
# separate repo, searched for independently of this ICESEE checkout's root.
REPO_SEARCH="$SCRIPT_DIR"
while [[ "$REPO_SEARCH" != "/" && ! -f "$REPO_SEARCH/scripts/activate.sh" ]]; do
  REPO_SEARCH="$(dirname "$REPO_SEARCH")"
done
ICESEE_REPO_ROOT="$(cd "$IDEALIZED_PIG_DIR/../../../.." && pwd)"
# ICESEE is imported as a package named after its checkout directory, and
# run_da_icepack.py relies on PYTHONPATH rather than adjusting sys.path.
if [[ "$(basename "$ICESEE_REPO_ROOT")" != "ICESEE" ]]; then
  echo "ERROR: the repository must be checked out in a directory named ICESEE (found $ICESEE_REPO_ROOT)" >&2
  exit 1
fi
export PYTHONPATH="$(dirname "$ICESEE_REPO_ROOT"):$ICESEE_REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"
if [[ -f "$REPO_SEARCH/scripts/activate.sh" ]]; then
  module purge || true
  # shellcheck disable=SC1091
  source "$REPO_SEARCH/scripts/activate.sh"
  MPI_DIR="$(spack location -i openmpi)"
  export PATH="$(echo "$PATH" | tr ':' '\n' | grep -Ev 'openmpi|mvapich|mpich' | paste -sd: -)"
  export PATH="$MPI_DIR/bin:$PATH"
  export LD_LIBRARY_PATH="$MPI_DIR/lib:$LD_LIBRARY_PATH"
  hash -r
else
  echo "WARNING: scripts/activate.sh not found -- assuming environment already active."
fi
echo "python3: $(which python3)"

# --- Unique, never-shared output/store location (item D) ---
export ICESEE_RUN1_OUTPUT="$IDEALIZED_PIG_DIR/pace_runs/run1_${SLURM_JOB_ID}"
mkdir -p "$ICESEE_RUN1_OUTPUT"
echo "output directory: $ICESEE_RUN1_OUTPUT"
# ICESEE removes and recreates data_path when its configuration is imported,
# so it gets its own fresh subdirectory; the manifest, log, and summary stay
# in ICESEE_RUN1_OUTPUT.
ICESEE_RUN1_DATA="$ICESEE_RUN1_OUTPUT/icesee_data"
echo "ICESEE data_path (member store, checkpoints): $ICESEE_RUN1_DATA"

# --- Reproducibility manifest (item C) ---
python3 "$ICESEE_REPO_ROOT/scripts/benchmarks/pace_run_manifest.py" \
  --out "$ICESEE_RUN1_OUTPUT/manifest.json" || echo "WARNING: manifest generation failed (non-fatal)"

# --- Fail-fast pre-checks (mirrors RUN0's checks, cheap to repeat here so
#     RUN1 never burns allocation on an environment that RUN0 didn't
#     actually validate for THIS job's node/rank layout) ---
srun --mpi=pmix -n "${SLURM_NTASKS}" python3 "$SCRIPT_DIR/run0_env_smoke.py" \
  --data-path "$ICESEE_RUN1_OUTPUT/_run0_precheck"
rm -rf "$ICESEE_RUN1_OUTPUT/_run0_precheck"

# --- The real calibration run ---
# num_years/obs_start_time/obs_max_time come from params_calibration.yaml
# (see this file's header for why); everything else is an ordinary,
# correctly-applied --key=value CLI override. params.yaml (production) is
# never touched.
RUN1_START=$(date +%s)
srun --mpi=pmix -n "${SLURM_NTASKS}" python3 "$IDEALIZED_PIG_DIR/run_da_icepack.py" \
  --force-params="$SCRIPT_DIR/params_calibration.yaml" \
  --Nens="$NENS" \
  --model_nprocs="$MODEL_NPROCS" \
  --execution_mode=3 \
  --compact_initialization=True \
  --initFile=data/extended_beta1000yrs_compact_idx20000.h5 \
  --random_field_method=graph \
  --use_member_streaming=True \
  --use_store_streaming_analysis=True \
  --member_store_backend=hdf5_member_major \
  --member_store_root="$ICESEE_RUN1_DATA/_mode3_member_store" \
  --data_path="$ICESEE_RUN1_DATA" \
  --verbose \
  2>&1 | tee "$ICESEE_RUN1_OUTPUT/run1_stdout.log"
RUN1_WALL=$(( $(date +%s) - RUN1_START ))
echo "RUN1 DA wall time: ${RUN1_WALL} s" | tee -a "$ICESEE_RUN1_OUTPUT/run1_stdout.log"

# --- Post-run instrumentation summary (item B's report list) ---
python3 "$SCRIPT_DIR/summarize_run1.py" \
  --output-dir "$ICESEE_RUN1_DATA" \
  --summary-dir "$ICESEE_RUN1_OUTPUT" \
  --wall-seconds "$RUN1_WALL" \
  --log-file "$ICESEE_RUN1_OUTPUT/run1_stdout.log" \
  --job-id "${SLURM_JOB_ID:-unknown}" \
  --nens "$NENS" --model-nprocs "$MODEL_NPROCS" --world-size "$SLURM_NTASKS"

# sacct is only meaningful once the job's accounting record is finalized;
# report it best-effort (never fails the job over an accounting query).
echo "=== sacct (best-effort, may be incomplete until job fully exits) ==="
sacct -j "${SLURM_JOB_ID:-}" --format=JobID,Elapsed,MaxRSS,NNodes,NTasks,AllocCPUS,State 2>/dev/null || true

echo "=== RUN1 complete: SUCCESS ==="
