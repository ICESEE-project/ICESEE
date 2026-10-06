#!/bin/bash
# ==============================================================================
# PACE RUN 0 -- environment/MPI/Firedrake smoke check (Idealized PIG Mode-3).
# Cheap, fast (~minutes). NOT the scalability result -- see run1_pig_calibration.
#
# Environment/SLURM pattern reused from the repo's own known-working
# idealized_pig launcher (../run_icesee_spack.sh / ../run_icesee_spack.sbatch),
# not invented here: module purge, source ICESEE-Spack's scripts/activate.sh,
# resolve the Spack-installed Open MPI explicitly onto PATH, srun --mpi=pmix.
#
# ==== STUDENT-EDITABLE (only these) ====
#SBATCH -J icesee-run0-smoke
#SBATCH -A gts-arobel3-atlas          # ACCOUNT
#SBATCH -p inferno                    # PARTITION/QOS
#SBATCH --partition=cpu-small
#SBATCH -t 00:15:00                   # WALLTIME -- generous for a smoke check
#SBATCH -N 1
#SBATCH -n 4                          # NP -- must be even (RUN0 exercises a P_model=2 split)
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH -o run0-smoke-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=CHANGE_ME@gatech.edu
# ==== END STUDENT-EDITABLE ====

set -euo pipefail

# Anchor every path off THIS SCRIPT's own location, never off
# $SLURM_SUBMIT_DIR (which depends on where the student happened to run
# `sbatch` from and should not matter).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IDEALIZED_PIG_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$IDEALIZED_PIG_DIR"

# --- Repository/environment identity (reproducibility manifest, item C) ---
echo "=== RUN0 startup diagnostics ==="
echo "hostname: $(hostname)"
echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-unset}"
echo "working directory: $IDEALIZED_PIG_DIR"
echo "git branch: $(git -C "$IDEALIZED_PIG_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
echo "git revision: $(git -C "$IDEALIZED_PIG_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
echo "git status (should be reviewed -- uncommitted reconciliation work is expected): $(git -C "$IDEALIZED_PIG_DIR" status --porcelain 2>/dev/null | wc -l) modified/untracked files"

# --- Environment: reuse the repo's own known-working Spack activation ---
# (../run_icesee_spack.sh's exact pattern -- not a new strategy)
# This is the ICESEE-Spack activation script, which is a SEPARATE repo
# from ICESEE itself (per README.md: "For Cluster installation and runs,
# see ICESEE-Spack") -- its location on disk is independent of this
# ICESEE checkout's own root, so it is searched for separately.
REPO_SEARCH="$SCRIPT_DIR"
while [[ "$REPO_SEARCH" != "/" && ! -f "$REPO_SEARCH/scripts/activate.sh" ]]; do
  REPO_SEARCH="$(dirname "$REPO_SEARCH")"
done
# The ICESEE repo's OWN root is a fixed, known number of levels above
# applications/icepack_model/examples/idealized_pig -- computed directly,
# never assumed to equal wherever ICESEE-Spack's activate.sh was found.
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
  echo "WARNING: scripts/activate.sh not found by walking up from $SCRIPT_DIR -- assuming the environment is already active (e.g. a pre-activated Spack/module session)."
fi

echo "python3: $(which python3)"
echo "srun: $(which srun)"
echo "mpirun: $(which mpirun 2>/dev/null || echo not-found)"

# --- Unique, never-shared output directory (item D) ---
export ICESEE_RUN0_OUTPUT="$IDEALIZED_PIG_DIR/pace_runs/run0_${SLURM_JOB_ID}"
mkdir -p "$ICESEE_RUN0_OUTPUT"
echo "output directory: $ICESEE_RUN0_OUTPUT"

# --- Reproducibility manifest ---
python3 "$ICESEE_REPO_ROOT/scripts/benchmarks/pace_run_manifest.py" \
  --out "$ICESEE_RUN0_OUTPUT/manifest.json" || echo "WARNING: manifest generation failed (non-fatal for RUN0)"

# --- The actual smoke check ---
srun --mpi=pmix -n "${SLURM_NTASKS}" python3 "$SCRIPT_DIR/run0_env_smoke.py" \
  --data-path "$ICESEE_RUN0_OUTPUT"

echo "=== RUN0 complete: SUCCESS ==="
