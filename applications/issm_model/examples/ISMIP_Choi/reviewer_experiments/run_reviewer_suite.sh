#!/usr/bin/env bash
set -euo pipefail

reviewer_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
example_dir="$(dirname "${reviewer_dir}")"
cd "${example_dir}"

reviewer_mpi_ranks="${REVIEWER_MPI_RANKS:-8}"
reviewer_nens="${REVIEWER_NENS:-40}"
reviewer_model_nprocs="${REVIEWER_MODEL_NPROCS:-1}"

for reviewer_config in \
    friction_inversion_hybrid.yaml \
    friction_enkf_only.yaml \
    wrong_friction_fixed.yaml; do
    echo "[reviewer-suite] Starting ${reviewer_config}"
    mpirun -np "${reviewer_mpi_ranks}" python run_da_issm.py \
        --Nens="${reviewer_nens}" \
        --model_nprocs="${reviewer_model_nprocs}" \
        -F "reviewer_experiments/${reviewer_config}"
    echo "[reviewer-suite] Completed ${reviewer_config}"
done
