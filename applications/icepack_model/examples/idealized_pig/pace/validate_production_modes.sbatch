#!/bin/bash
#SBATCH --job-name=PIG-production-check
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=pig-production-check-%j.out
set -euo pipefail
: "${ICESEE_SOURCE_ROOT:?}" "${ICESEE_SPACK_ROOT:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
unset SLURM_OVERCOMMIT SLURM_OVERLAP
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail
[[ "$SLURM_JOB_NUM_NODES" == 2 && "$SLURM_NTASKS" == 40 ]] || exit 2
mapfile -t HOSTS < <(scontrol show hostnames "$SLURM_JOB_NODELIST")
[[ "${#HOSTS[@]}" == 2 ]] || exit 2
RUN="$PIG_DIR/pace_runs/production-check-${SLURM_JOB_ID:?}"
[[ ! -e "$RUN" ]] || exit 2
mkdir -p "$RUN"
{ git -C "$ICESEE_SOURCE_ROOT" rev-parse HEAD; git -C "$ICESEE_SOURCE_ROOT" status --short; } > "$RUN/source.txt"
export ICESEE_DEBUG_INIT=1 ICESEE_HANG_DUMP_S=300 PYTHONUNBUFFERED=1
PIDS=()
for mode in 1 2 3; do
 mkdir "$RUN/mode$mode"
 args=(--execution_mode="$mode" --Nens=32 --num_years=5 --dt=0.05 --obs_start_time=2 --obs_max_time=4 --freq_obs=1 '--save_steps=[2,3,4]' --data_path="$RUN/mode$mode/data")
 if [[ "$mode" == 3 ]]; then
  placement=(-N2 -n16 --ntasks-per-node=8)
  args+=(--model_nprocs=2)
 else
  placement=(-N1 -n8 --ntasks-per-node=8 --nodelist="${HOSTS[$((mode-1))]}")
  args+=(--model_nprocs=1)
 fi
 (
  set +e
  srun --mpi=pmix --exact --exclusive --mem=32G "${placement[@]}" --cpus-per-task=1 --hint=nomultithread --distribution=block:block --cpu-bind=cores "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" "${args[@]}" > "$RUN/mode$mode/stdout.log" 2> "$RUN/mode$mode/stderr.log"
  printf '%s\n' "$?" > "$RUN/mode$mode/exit-status"
 ) &
 PIDS+=("$!")
done
for pid in "${PIDS[@]}"; do wait "$pid"; done
FAILED=0
for mode in 1 2 3; do
 result="$(cat "$RUN/mode$mode/exit-status")"
 printf 'Mode %s exit=%s\n' "$mode" "$result"
 [[ "$result" == 0 ]] || FAILED=1
done
printf 'Results: %s\n' "$RUN"
exit "$FAILED"
