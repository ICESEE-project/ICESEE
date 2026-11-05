#!/bin/bash
#SBATCH --job-name=PIG-validate
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-medium
#SBATCH --mem-per-cpu=6G
#SBATCH --output=pig-validate-%x-%j.out
# Bounded validation stages for the Idealized PIG runtime repairs. Choose the
# stage and its resources at submission (nothing here submits anything):
#
#   STAGE=reference  -N1 -n1  -t 00:30:00  serial instrumented 5-year reference
#                                          trajectory (tools/profile_reference_trajectory.py)
#   STAGE=mode2      -N1 -n8  -t 02:00:00  the reported command: Nens=32, mode 2,
#                                          model_nprocs=1, 5 years, rank traces and
#                                          periodic stack dumps
#   STAGE=mode3      -N2 -n4 --ntasks-per-node=2 -t 00:20:00
#                                          the failing two-member, model_nprocs=2,
#                                          one-year diagnostic; expected to stop on
#                                          every rank at the k=10 analysis with
#                                          "spread in observation space has collapsed"
#
# Every run writes to a new directory under pace_runs/ and refuses an existing one.
set -euo pipefail
: "${STAGE:?Set STAGE=reference|mode2|mode3}"
: "${ICESEE_SOURCE_ROOT:?Set the ICESEE source directory}"
: "${ICESEE_SPACK_ROOT:?Set the retained Spack installation directory}"
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail
RUN_DIR="$PIG_DIR/pace_runs/validate-${STAGE}-${SLURM_JOB_ID:?}"
[[ ! -e "$RUN_DIR" ]] || { echo "Refusing existing run directory: $RUN_DIR" >&2; exit 2; }
mkdir -p "$RUN_DIR"
{ git -C "$ICESEE_ROOT" rev-parse HEAD; git -C "$ICESEE_ROOT" status --short; } > "$RUN_DIR/source.txt"
SRUN=(srun --mpi=pmix --distribution=block:block --cpu-bind=cores)

case "$STAGE" in
  reference)
    export PYOP2_CACHE_DIR="/tmp/icesee-$SLURM_JOB_ID/pyop2"
    export FIREDRAKE_TSFC_KERNEL_CACHE_DIR="/tmp/icesee-$SLURM_JOB_ID/tsfc"
    /usr/bin/time -v python "$PIG_DIR/tools/profile_reference_trajectory.py" \
      --output "$RUN_DIR/profile" --steps=100 --dt=0.05 --experiment=true \
      > "$RUN_DIR/stdout.log" 2> "$RUN_DIR/stderr.log"
    ;;
  mode2)
    export ICESEE_DEBUG_INIT=1 ICESEE_HANG_DUMP_S=300
    "${SRUN[@]}" -n 8 "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" \
      --Nens=32 --execution_mode=2 --model_nprocs=1 --num_years=5 \
      --obs_start_time=1 --obs_max_time=5 --freq_obs=1 \
      --data_path="$RUN_DIR/data" \
      > "$RUN_DIR/stdout.log" 2> "$RUN_DIR/stderr.log"
    ;;
  mode3)
    export ICESEE_DEBUG_INIT=1 ICESEE_HANG_DUMP_S=300
    set +e
    "${SRUN[@]}" -N 2 -n 4 --ntasks-per-node=2 "$PIG_DIR/pace/launch_rank.sh" \
      "$PIG_DIR/run_da_icepack.py" --Nens=2 --execution_mode=3 --mode3_experimental=True \
      --model_nprocs=2 --num_years=1 --dt=0.05 --obs_start_time=0.5 --obs_max_time=1 \
      --freq_obs=0.5 '--save_steps=[0.05,0.1,1.0]' --checkpoint_every=1 \
      --checkpoint_keep_last=2 --data_path="$RUN_DIR/data" \
      > "$RUN_DIR/stdout.log" 2> "$RUN_DIR/stderr.log"
    status=$?
    set -e
    refused=$(grep -c "spread in observation space has collapsed" "$RUN_DIR/stdout.log" "$RUN_DIR/stderr.log" | awk -F: '{s+=$2} END {print s+0}')
    echo "mode3 exit status $status; ranks reporting the refusal: $refused" | tee "$RUN_DIR/result.txt"
    ;;
  *) echo "Unknown STAGE=$STAGE" >&2; exit 2 ;;
esac
printf 'Results: %s\n' "$RUN_DIR"
