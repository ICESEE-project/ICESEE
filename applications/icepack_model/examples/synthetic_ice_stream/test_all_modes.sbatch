#!/bin/bash
#SBATCH --job-name=ICESEE-all-modes
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=1
#SBATCH --ntasks=13
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=00:20:00
#SBATCH --output=synthetic-all-modes-%j.out
set -euo pipefail
: "${ICESEE_SOURCE_ROOT:?}" "${ICESEE_SPACK_ROOT:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
unset SLURM_OVERCOMMIT SLURM_OVERLAP
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail
SYN="$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/synthetic_ice_stream"
RUN="$SYN/validation_runs/all-modes-${SLURM_JOB_ID:?}"
[[ ! -e "$RUN" ]] || exit 2
mkdir -p "$RUN/output-check"
{ git -C "$ICESEE_SOURCE_ROOT" rev-parse HEAD; git -C "$ICESEE_SOURCE_ROOT" status --short; } > "$RUN/source.txt"
# Validate read-only analysis, window rollover and restart in the MPI-HDF5 environment.
mkdir -p "$RUN/io-check"
srun --mpi=pmix --exact --exclusive -N1 -n8 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "$PIG_DIR/pace/launch_rank.sh" "$ICESEE_SOURCE_ROOT/src/tests/parallel_mpi/_enkf_io_readonly_analysis_worker.py" "$RUN/io-check" > "$RUN/io-check.log" 2>&1
[[ "$(grep -c 'RESULT rank=.*errors=0 total_errors=0' "$RUN/io-check.log")" == 8 ]]
# Exercise the replacement output path at PIG's real dimensions first.
srun --mpi=pmix --exact --exclusive -N1 -n8 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "$PIG_DIR/pace/launch_rank.sh" "$ICESEE_SOURCE_ROOT/src/tests/parallel_mpi/_analysis_production_shape_worker.py" "$RUN/output-check" > "$RUN/output-check.log" 2>&1
[[ "$(grep -c 'WRITE_RESULT rank=.*errors=0 total_errors=0' "$RUN/output-check.log")" == 8 ]]
# Concurrent independent steps consume exactly 1+4+4+4 allocated CPUs.
PIDS=()
for mode in 0 1 2 3; do
 mkdir -p "$RUN/mode$mode"
 ranks=4; [[ "$mode" != 0 ]] || ranks=1
 (
  set +e
  export ICESEE_DEBUG_INIT=1 ICESEE_HANG_DUMP_S=120
  srun --mpi=pmix --exact --exclusive --mem=8G -N1 -n "$ranks" --cpus-per-task=1 --hint=nomultithread --distribution=block:block --cpu-bind=cores "$PIG_DIR/pace/launch_rank.sh" "$SYN/run_da_icepack.py" --execution_mode="$mode" --model_nprocs=1 --Nens=8 --num_years=5 --obs_start_time=2 --obs_max_time=4 --freq_obs=1 --data_path="$RUN/mode$mode/data" > "$RUN/mode$mode/stdout.log" 2> "$RUN/mode$mode/stderr.log"
  result=$?
  printf '%s\n' "$result" > "$RUN/mode$mode/exit-status"
  exit 0
 ) &
 PIDS+=("$!")
done
for pid in "${PIDS[@]}"; do wait "$pid"; done
FAILED=0
for mode in 0 1 2 3; do
 result="$(cat "$RUN/mode$mode/exit-status")"
 printf 'Mode %s exit=%s\n' "$mode" "$result"
 [[ "$result" == 0 ]] || FAILED=1
done
printf 'Results: %s\n' "$RUN"
exit "$FAILED"
