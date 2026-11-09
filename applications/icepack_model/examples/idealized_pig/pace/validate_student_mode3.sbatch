#!/bin/bash
#SBATCH --job-name=PIG-mode3-24
#SBATCH --account=gts-arobel3-atlas
#SBATCH --partition=cpu-small
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=pig-mode3-24-%j.out
set -euo pipefail
: "${ICESEE_SOURCE_ROOT:?Set the updated ICESEE source checkout}"
: "${ICESEE_SPACK_ROOT:?Set the installed dependency environment}"
[[ "${SLURM_JOB_NUM_NODES:?}" == 1 && "${SLURM_NTASKS:?}" == 8 ]] || {
  echo "This confirmation requires one node and eight MPI ranks" >&2; exit 2;
}
export ICESEE_NENS=24 ICESEE_MODEL_NPROCS=2 ICESEE_TASKS_PER_NODE=8
export ICESEE_DATA_PATH="$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/pace_runs/student-mode3-${SLURM_JOB_ID:?}/data"
exec bash "$ICESEE_SOURCE_ROOT/applications/icepack_model/examples/idealized_pig/run_idealized_1.sbatch"
