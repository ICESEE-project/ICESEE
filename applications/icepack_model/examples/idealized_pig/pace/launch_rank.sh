#!/bin/bash
set -euo pipefail
if [[ -n "${ICESEE_JIT_GCC_MODULE:-}" ]]; then
  if type module >/dev/null 2>&1; then
    module load "$ICESEE_JIT_GCC_MODULE"
  elif [[ -x "${LMOD_CMD:-}" ]]; then
    eval "$("$LMOD_CMD" bash load "$ICESEE_JIT_GCC_MODULE")"
  else
    echo 'The requested JIT compiler requires the PACE module environment' >&2
    exit 2
  fi
  export OMPI_CC="$(command -v gcc)"
  printf '#include <stddef.h>\nint main(void){return 0;}\n' | mpicc -x c -fsyntax-only -
fi
local_tmp="${TMPDIR:-/tmp}"
# PACE's measured working layout used /tmp. Reject shared scratch for caches.
fs_type="$(stat -f -c %T "$local_tmp" 2>/dev/null || true)"
case "$fs_type" in nfs*|lustre*|gpfs*|panfs*|cifs*|smb*) local_tmp=/tmp ;; esac
# Concurrent Slurm steps reuse rank IDs; keep their SQLite/JIT namespaces separate.
job_cache="$local_tmp/icesee-${SLURM_JOB_ID:?}/step-${SLURM_STEP_ID:-batch}"
export PYOP2_CACHE_DIR="$job_cache/pyop2"
export FIREDRAKE_TSFC_KERNEL_CACHE_DIR="$job_cache/tsfc"
export XDG_CACHE_HOME="$job_cache/xdg/rank${SLURM_PROCID:?}"
mkdir -p "$PYOP2_CACHE_DIR" "$FIREDRAKE_TSFC_KERNEL_CACHE_DIR" "$XDG_CACHE_HOME"
echo "[rank ${SLURM_PROCID}] host=$(hostname) local=${SLURM_LOCALID:-} cache=$job_cache TMPDIR=${TMPDIR:-unset}"
exec python "$@"
