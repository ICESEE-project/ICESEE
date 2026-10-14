#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -J PIG-DA40
#SBATCH -N 3
#SBATCH -n 60
#SBATCH --ntasks-per-node=20
#SBATCH -p cpu-medium
#SBATCH --mem=256G
#SBATCH -A gts-arobel3-atlas
#SBATCH -o pig_da40_Ne30_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bkyanjo3@gatech.edu

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

# ----------------------------------------------------------------------
# ICESEE-Spack environment
# ----------------------------------------------------------------------
module purge || true

source /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ICESEE-Spack/scripts/activate.sh

MPI_DIR="$(spack -e \
  /storage/home/hcoda1/8/bkyanjo3/r-arobel3-0/ICESEE-Spack/.spack-env/icesee \
  location -i openmpi)"

export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:${LD_LIBRARY_PATH:-}"

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
echo "Output directory   : $OUTDIR"
echo "============================================================"

# ----------------------------------------------------------------------
# IMPORTANT:
# Give every MPI rank its own Firedrake / PyOP2 / Loopy cache.
#
# We previously encountered:
#   sqlite3.DatabaseError: database disk image is malformed
#
# A fresh cache allowed Firedrake to initialize normally.
# ----------------------------------------------------------------------

cat > "${OUTDIR}/launch_rank.sh" <<'EOF'
#!/bin/bash

set -euo pipefail

export XDG_CACHE_HOME="${TMPDIR:-/tmp}/icesee-${SLURM_JOB_ID}-${SLURM_PROCID}"
mkdir -p "$XDG_CACHE_HOME"

echo "[rank ${SLURM_PROCID}] host=$(hostname) cache=$XDG_CACHE_HOME"

exec python "$@"
EOF

chmod +x "${OUTDIR}/launch_rank.sh"

# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------

echo
echo "Starting ICESEE at $(date)"
echo

/usr/bin/time -v \
srun --mpi=pmix \
     -n "$NP" \
     --ntasks-per-node=20 \
     "${OUTDIR}/launch_rank.sh" \
     run_da_icepack.py \
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

echo
echo "ICESEE completed at $(date)"
echo "Results: $OUTDIR"
