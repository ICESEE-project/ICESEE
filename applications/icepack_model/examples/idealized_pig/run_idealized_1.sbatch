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

# ======================================================================
# ICESEE - Idealized PIG Mode-3 Data Assimilation
#
# Automatically discovers ICESEE-Spack and its Spack environment.
#
# Supported installation layouts:
#
#   ICESEE-Spack/
#       scripts/activate.sh
#       .spack-env/icesee/
#       ICESEE/
#
# The script may be submitted from within ICESEE-Spack or its
# ICESEE subtree.
#
# Optional override:
#
#   export ICESEE_SPACK_ROOT=/path/to/ICESEE-Spack
#   sbatch run_idealized_1.sbatch
#
# ======================================================================


# ----------------------------------------------------------------------
# 1. Working directory
# ----------------------------------------------------------------------

cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "============================================================"
echo "ICESEE - Idealized PIG Mode-3"
echo "============================================================"
echo "Submission directory: $PWD"
echo "Hostname            : $(hostname)"
echo "Job ID              : ${SLURM_JOB_ID:-unknown}"
echo "============================================================"


# ----------------------------------------------------------------------
# 2. Automatically discover ICESEE-Spack
# ----------------------------------------------------------------------

find_icesee_spack() {

    local dir="$1"

    [[ -d "$dir" ]] || return 1

    dir="$(cd "$dir" && pwd -P)"

    while true; do

        # Case A: current directory is ICESEE-Spack
        if [[ -f "$dir/scripts/activate.sh" &&
              -d "$dir/.spack-env/icesee" ]]; then

            printf '%s\n' "$dir"
            return 0
        fi

        # Case B: ICESEE-Spack is a child directory
        if [[ -f "$dir/ICESEE-Spack/scripts/activate.sh" &&
              -d "$dir/ICESEE-Spack/.spack-env/icesee" ]]; then

            printf '%s\n' "$dir/ICESEE-Spack"
            return 0
        fi

        # Move upward through the directory tree
        [[ "$dir" == "/" ]] && break

        dir="$(dirname "$dir")"

    done

    return 1
}


# ----------------------------------------------------------------------
# 3. Resolve ICESEE-Spack installation
# ----------------------------------------------------------------------

ICESEE_SPACK_ROOT="${ICESEE_SPACK_ROOT:-}"

if [[ -n "$ICESEE_SPACK_ROOT" ]]; then

    ICESEE_SPACK_ROOT="$(
        find_icesee_spack "$ICESEE_SPACK_ROOT"
    )" || {
        echo "ERROR: Invalid ICESEE_SPACK_ROOT." >&2
        exit 1
    }

else

    # First: search from submission directory
    ICESEE_SPACK_ROOT="$(
        find_icesee_spack "$PWD"
    )" || true

    # Second: search from script location
    if [[ -z "$ICESEE_SPACK_ROOT" ]]; then

        SCRIPT_DIR="$(
            cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P
        )"

        ICESEE_SPACK_ROOT="$(
            find_icesee_spack "$SCRIPT_DIR"
        )" || true

    fi

fi


if [[ -z "$ICESEE_SPACK_ROOT" ]]; then

    echo "============================================================" >&2
    echo "ERROR: ICESEE-Spack installation not found." >&2
    echo >&2
    echo "Expected:" >&2
    echo "  ICESEE-Spack/scripts/activate.sh" >&2
    echo "  ICESEE-Spack/.spack-env/icesee/" >&2
    echo >&2
    echo "Set ICESEE_SPACK_ROOT explicitly if necessary." >&2
    echo "============================================================" >&2

    exit 1
fi


echo
echo "[ICESEE] Installation detected:"
echo "         $ICESEE_SPACK_ROOT"
echo


# ----------------------------------------------------------------------
# 4. Activate ICESEE-Spack environment
# ----------------------------------------------------------------------

module purge || true

ACTIVATE_SCRIPT="$ICESEE_SPACK_ROOT/scripts/activate.sh"
SPACK_ENV="$ICESEE_SPACK_ROOT/.spack-env/icesee"

if [[ ! -f "$ACTIVATE_SCRIPT" ]]; then
    echo "ERROR: Activation script missing: $ACTIVATE_SCRIPT" >&2
    exit 1
fi

if [[ ! -d "$SPACK_ENV" ]]; then
    echo "ERROR: Spack environment missing: $SPACK_ENV" >&2
    exit 1
fi

echo "[ICESEE] Activating environment..."

source "$ACTIVATE_SCRIPT"

if ! command -v spack >/dev/null 2>&1; then
    echo "ERROR: Spack is unavailable after activation." >&2
    exit 1
fi


# ----------------------------------------------------------------------
# 5. Configure OpenMPI
# ----------------------------------------------------------------------

MPI_DIR="$(spack -e "$SPACK_ENV" location -i openmpi)"

if [[ ! -d "$MPI_DIR" ]]; then
    echo "ERROR: OpenMPI installation not found: $MPI_DIR" >&2
    exit 1
fi

export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:${LD_LIBRARY_PATH:-}"

echo "[ICESEE] OpenMPI: $MPI_DIR"

command -v mpirun || true
command -v python


# ----------------------------------------------------------------------
# 6. Simulation configuration
# ----------------------------------------------------------------------

NENS=20
MODEL_NPROCS=2

NP=$((NENS * MODEL_NPROCS))

NUM_YEARS=5
DT=0.05

OBS_START=1
OBS_MAX=5
FREQ_OBS=1


# ----------------------------------------------------------------------
# 7. Validate MPI resources
# ----------------------------------------------------------------------

if [[ -n "${SLURM_NTASKS:-}" ]]; then

    if (( NP > SLURM_NTASKS )); then

        echo "ERROR: Insufficient MPI ranks allocated." >&2
        echo "Required: $NP" >&2
        echo "Allocated: $SLURM_NTASKS" >&2

        exit 1
    fi

fi


# ----------------------------------------------------------------------
# 8. Output directory
# ----------------------------------------------------------------------

OUTDIR="pig_da40_Ne${NENS}_${SLURM_JOB_ID:-local}"

mkdir -p "$OUTDIR"

OUTDIR="$(cd "$OUTDIR" && pwd -P)"

echo
echo "============================================================"
echo "Idealized PIG Mode-3 DA"
echo "============================================================"
echo "Job ID             : ${SLURM_JOB_ID:-unknown}"
echo "Nodes              : ${SLURM_JOB_NUM_NODES:-unknown}"
echo "Allocated ranks    : ${SLURM_NTASKS:-unknown}"
echo "MPI ranks used     : $NP"
echo "Ensembles          : $NENS"
echo "Ranks/model        : $MODEL_NPROCS"
echo "Years              : $NUM_YEARS"
echo "dt                 : $DT"
echo "Observation start  : $OBS_START"
echo "Observation end    : $OBS_MAX"
echo "Observation freq.  : $FREQ_OBS"
echo "Output directory   : $OUTDIR"
echo "ICESEE-Spack       : $ICESEE_SPACK_ROOT"
echo "============================================================"


# ----------------------------------------------------------------------
# 9. Per-rank Firedrake / PyOP2 / Loopy cache
#
# Prevents multiple MPI ranks from sharing the same SQLite-backed
# compilation cache.
#
# Each MPI rank receives its own cache directory.
# ----------------------------------------------------------------------

cat > "${OUTDIR}/launch_rank.sh" <<'EOF'
#!/bin/bash

set -euo pipefail

RANK="${SLURM_PROCID:-0}"
JOB="${SLURM_JOB_ID:-local}"

export XDG_CACHE_HOME="${TMPDIR:-/tmp}/icesee-${JOB}-${RANK}"

mkdir -p "$XDG_CACHE_HOME"

echo "[rank ${RANK}] host=$(hostname) cache=$XDG_CACHE_HOME"

exec python "$@"
EOF

chmod +x "${OUTDIR}/launch_rank.sh"


# ----------------------------------------------------------------------
# 10. Resolve ICESEE application
# ----------------------------------------------------------------------

PIG_DIR="$ICESEE_SPACK_ROOT/ICESEE/applications/icepack_model/examples/idealized_pig"

if [[ ! -f "$PIG_DIR/run_da_icepack.py" ]]; then

    echo "ERROR: Idealized PIG application not found." >&2
    echo "Expected: $PIG_DIR/run_da_icepack.py" >&2

    exit 1
fi

echo "[ICESEE] Application: $PIG_DIR"


# ----------------------------------------------------------------------
# 11. Run Mode-3 data assimilation
# ----------------------------------------------------------------------

echo
echo "Starting ICESEE at $(date)"
echo

/usr/bin/time -v \
srun --mpi=pmix \
     -n "$NP" \
     --ntasks-per-node=20 \
     "${OUTDIR}/launch_rank.sh" \
     "$PIG_DIR/run_da_icepack.py" \
       --Nens="$NENS" \
       --execution_mode=3 \
       --model_nprocs="$MODEL_NPROCS" \
       --num_years="$NUM_YEARS" \
       --dt="$DT" \
       --obs_start_time="$OBS_START" \
       --obs_max_time="$OBS_MAX" \
       --freq_obs="$FREQ_OBS" \
       --data_path="${OUTDIR}/data" \
     > "${OUTDIR}/stdout.log" \
     2> "${OUTDIR}/stderr.log"


# ----------------------------------------------------------------------
# 12. Completion
# ----------------------------------------------------------------------

echo
echo "============================================================"
echo "ICESEE completed successfully"
echo "============================================================"
echo "Finished at : $(date)"
echo "Results     : $OUTDIR"
echo "Stdout      : $OUTDIR/stdout.log"
echo "Stderr      : $OUTDIR/stderr.log"
echo "============================================================"
