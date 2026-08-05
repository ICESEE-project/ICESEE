#!/usr/bin/env bash
set -euo pipefail

experiment_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
example_dir="$(dirname "${experiment_dir}")"
cd "${example_dir}"

mpi_ranks="${REBUTTAL_MPI_RANKS:-8}"
nens="${REBUTTAL_NENS:-40}"
model_nprocs="${REBUTTAL_MODEL_NPROCS:-1}"

for config in param_ibf.yaml param_wbf.yaml param_ebf.yaml; do
    echo "[rebuttal-suite] Starting ${config}"
    mpirun -np "${mpi_ranks}" python run_da_issm.py \
        --Nens="${nens}" \
        --model_nprocs="${model_nprocs}" \
        -F "rebutal_experiments/${config}"
    echo "[rebuttal-suite] Completed ${config}"
done
