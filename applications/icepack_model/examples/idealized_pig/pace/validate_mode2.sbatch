#!/bin/bash
#SBATCH --job-name=PIG-mode2-check
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=00:10:00
#SBATCH --output=mode2-check-%j.out
# Default: small real-MPI I/O lifecycle test. STAGE=pig requires a longer
# explicit time allocation and reproduces 32 members / 5 years / observations 2,3,4.
set -euo pipefail
: "${ICESEE_SOURCE_ROOT:?Set integration scientific source}"
: "${ICESEE_SPACK_ROOT:?Set dependency installation}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE
unset SLURM_OVERCOMMIT SLURM_OVERLAP
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail
[[ "$SLURM_NTASKS" == 8 ]] || { echo 'Validation requires exactly eight ranks' >&2; exit 2; }
RUN="$PIG_DIR/pace_runs/validation-${STAGE:-io}-${SLURM_JOB_ID:?}"
[[ ! -e "$RUN" ]] || { echo "Existing output: $RUN" >&2; exit 2; }
mkdir -p "$RUN/io"
{ git -C "$ICESEE_ROOT" rev-parse HEAD; git -C "$ICESEE_ROOT" status --short; git -C "$ICESEE_ROOT" diff HEAD --stat; } > "$RUN/source.txt"
SRUN=(srun --mpi=pmix --exact --exclusive --nodes=1 --ntasks=8 --ntasks-per-node=8 --cpus-per-task=1 --hint=nomultithread --distribution=block:block --cpu-bind=cores)
case "${STAGE:-io}" in
 io)
  "${SRUN[@]}" "$PIG_DIR/pace/launch_rank.sh" "$ICESEE_ROOT/src/tests/parallel_mpi/_enkf_io_readonly_analysis_worker.py" "$RUN/io" > "$RUN/stdout.log" 2> "$RUN/stderr.log"
  [[ "$(grep -c 'RESULT rank=.*errors=0 total_errors=0' "$RUN/stdout.log")" == 8 ]]
  "${SRUN[@]}" "$PIG_DIR/pace/launch_rank.sh" "$ICESEE_ROOT/src/tests/parallel_mpi/_analysis_production_shape_worker.py" "$RUN/production-shape" > "$RUN/production-shape.stdout.log" 2> "$RUN/production-shape.stderr.log"
  [[ "$(grep -c 'WRITE_RESULT rank=.*errors=0 total_errors=0' "$RUN/production-shape.stdout.log")" == 8 ]]
  ;;
 pig)
  export ICESEE_DEBUG_INIT=1 ICESEE_HANG_DUMP_S=300
  "${SRUN[@]}" "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" --execution_mode=2 --Nens=32 --model_nprocs=1 --num_years=5 --dt=0.05 --obs_start_time=2 --obs_max_time=4 --freq_obs=1 '--save_steps=[2,3,4]' --data_path="$RUN/_modelrun_datasets" > "$RUN/stdout.log" 2> "$RUN/stderr.log"
  ;;
 *) echo 'Set STAGE=io or STAGE=pig' >&2; exit 2;;
esac
printf 'PASS: %s; output %s\n' "${STAGE:-io}" "$RUN"
