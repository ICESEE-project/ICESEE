#!/bin/bash
#SBATCH --job-name=PIG-runtime-check
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-medium
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=pig-runtime-check-%j.out
set -euo pipefail
: "${ICESEE_SOURCE_ROOT:?Set the ICESEE subtree source directory}"
: "${ICESEE_SPACK_ROOT:?Set the retained Spack installation directory}"
source "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace/resolve_environment.sh"
set -euo pipefail
RUN_DIR="$PIG_DIR/pace_runs/runtime-check-${SLURM_JOB_ID:?}"
[[ ! -e "$RUN_DIR" ]] || { echo "Refusing existing diagnostic run: $RUN_DIR" >&2; exit 2; }
mkdir -p "$RUN_DIR"
DATA_DIR="$RUN_DIR/data"
export ICESEE_DEBUG_INIT=1 ICESEE_HANG_DUMP_S=300
ARGS=(--Nens=2 --execution_mode=3 --model_nprocs=2 --num_years=0.1 --dt=0.05
      --obs_start_time=0.05 --obs_max_time=0.1 --freq_obs=0.05
      '--save_steps=[0.05,0.1]' --checkpoint_every=1
      --checkpoint_keep_last=2 --data_path="$DATA_DIR")
printf 'parent commit: '; git -C "$ICESEE_ROOT" rev-parse HEAD
srun --mpi=pmix --nodes=2 --ntasks=4 --ntasks-per-node=2 --distribution=block:block --cpu-bind=cores \
  "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" "${ARGS[@]}" \
  > "$RUN_DIR/fresh.stdout.log" 2> "$RUN_DIR/fresh.stderr.log"
python - "$DATA_DIR" "$RUN_DIR/profile-before.json" <<'PY'
from pathlib import Path
import hashlib, json, sys
import h5py, numpy as np
root = Path(sys.argv[1])
paths = [root / 'hs_profiles_true.h5'] + [root / f'hs_ensemble_profiles{i}.h5' for i in range(2)]
for path in paths:
    with h5py.File(path) as handle:
        np.testing.assert_array_equal(handle['written_steps'][:], [0, 1, 2])
        assert handle['h_profiles'].shape[0] > 0
        assert np.isfinite(handle['h_profiles'][:]).all()
        assert np.isfinite(handle['s_profiles'][:]).all()
        print(path, handle['h_profiles'].shape)
Path(sys.argv[2]).write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
PY
srun --mpi=pmix --nodes=2 --ntasks=4 --ntasks-per-node=2 --distribution=block:block --cpu-bind=cores \
  "$PIG_DIR/pace/launch_rank.sh" "$PIG_DIR/run_da_icepack.py" "${ARGS[@]}" --resume_from_checkpoint=true \
  > "$RUN_DIR/resume.stdout.log" 2> "$RUN_DIR/resume.stderr.log"
python - "$RUN_DIR/profile-before.json" <<'PY'
import hashlib, json, sys
from pathlib import Path
for name, expected in json.loads(Path(sys.argv[1]).read_text()).items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == expected, name
print('Fresh run and completed-checkpoint resume exited successfully; profile files preserved.')
PY
printf 'Runtime diagnostic results: %s\n' "$RUN_DIR"
