#!/bin/bash
#SBATCH --job-name=PIG-year-1-3
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=2
#SBATCH --ntasks=12
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=00:20:00
#SBATCH --output=pig-year-1-3-%j.out
set -euo pipefail
: "${ICESEE_SOURCE_ROOT:?}" "${ICESEE_SPACK_ROOT:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
unset SLURM_OVERCOMMIT SLURM_OVERLAP
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail
RUN="$PIG_DIR/pace_runs/year-1-3-${SLURM_JOB_ID:?}"
[[ ! -e "$RUN" ]] || exit 2
mkdir -p "$RUN"
{ git -C "$ICESEE_SOURCE_ROOT" rev-parse HEAD; git -C "$ICESEE_SOURCE_ROOT" status --short; } > "$RUN/source.txt"
export ICESEE_DEBUG_INIT=1 ICESEE_HANG_DUMP_S=120 PYTHONUNBUFFERED=1
srun --mpi=pmix --exact --exclusive --mem=8G -N2 -n4 --ntasks-per-node=2 --cpus-per-task=1 --hint=nomultithread --cpu-bind=cores "$PIG_DIR/pace/launch_rank.sh" "$ICESEE_SOURCE_ROOT/src/tests/parallel_mpi/_native_process_noise_worker.py" > "$RUN/noise-check.log" 2>&1
[[ "$(grep -c 'NOISE_RESULT rank=.*errors=0 total_errors=0' "$RUN/noise-check.log")" == 4 ]]
PIDS=()
for mode in 1 2 3; do
 mkdir -p "$RUN/mode$mode"
 args=(--execution_mode="$mode" --Nens=2 --num_years=1 --dt=0.05 --obs_start_time=0.5 --obs_max_time=1 --freq_obs=0.25 '--save_steps=[0.5,0.75,1.0]' --data_path="$RUN/mode$mode/data")
 placement=(-N1 -n2 --ntasks-per-node=2)
 if [[ "$mode" == 3 ]]; then
  args+=(--model_nprocs=2 --checkpoint_keep_last=0)
  placement=(-N2 -n4 --ntasks-per-node=2)
 else
  args+=(--model_nprocs=1)
 fi
 (
  set +e
  srun --mpi=pmix --exact --exclusive --mem=24G "${placement[@]}" --cpus-per-task=1 --hint=nomultithread --distribution=block:block --cpu-bind=cores "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" "${args[@]}" > "$RUN/mode$mode/stdout.log" 2> "$RUN/mode$mode/stderr.log"
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
if [[ "$FAILED" == 0 ]]; then
 mkdir "$RUN/mode3-resume"
 python - "$RUN/mode3/data" "$RUN/mode3-resume/data" <<'PY'
from pathlib import Path
import shutil,sys
source,target=map(Path,sys.argv[1:3])
assert not target.exists()
assert sum(p.stat().st_size for p in source.rglob("*") if p.is_file()) < 1024**3
shutil.copytree(source,target)
steps=target/"_mode3_state_history/steps"
assert (steps/"checkpoint_00000014/manifest.json").is_file()
for directory in steps.glob("checkpoint_*"):
    if int(directory.name.split("_")[-1])>=15:
        shutil.rmtree(directory)
PY
 srun --mpi=pmix --exact --exclusive --mem=24G -N2 -n4 --ntasks-per-node=2 --cpus-per-task=1 --hint=nomultithread --distribution=block:block --cpu-bind=cores "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" --execution_mode=3 --model_nprocs=2 --Nens=2 --num_years=1 --dt=0.05 --obs_start_time=0.5 --obs_max_time=1 --freq_obs=0.25 '--save_steps=[0.5,0.75,1.0]' --checkpoint_keep_last=0 --resume_from_checkpoint=True --data_path="$RUN/mode3-resume/data" > "$RUN/mode3-resume/stdout.log" 2> "$RUN/mode3-resume/stderr.log"
 python "$ICESEE_SOURCE_ROOT/src/tests/parallel_mpi/_validate_native_resume.py" "$RUN/mode3/data" "$RUN/mode3-resume/data" > "$RUN/resume-comparison.log"
fi
printf 'Results: %s\n' "$RUN"
exit "$FAILED"
