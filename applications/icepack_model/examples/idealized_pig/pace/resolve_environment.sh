#!/bin/bash
# Resolve the live source separately from the retained Spack installation.
icesee_find_root() {
  local candidate="$1"
  [[ -d "$candidate" ]] || return 1
  candidate="$(cd "$candidate" && pwd -P)"
  while [[ "$candidate" != / ]]; do
    if [[ -f "$candidate/scripts/activate.sh" && -d "$candidate/.spack-env/icesee" ]]; then
      printf '%s\n' "$candidate"; return 0
    fi
    candidate="$(dirname "$candidate")"
  done
  return 1
}
if [[ -n "${ICESEE_SPACK_ROOT:-}" ]]; then
  SPACK_REPO="$(icesee_find_root "$ICESEE_SPACK_ROOT")" || { echo "Invalid ICESEE_SPACK_ROOT" >&2; return 2; }
else
  SPACK_REPO="$(icesee_find_root "${SLURM_SUBMIT_DIR:-$PWD}" || icesee_find_root "$(dirname "${BASH_SOURCE[0]}")")" || {
    echo "Set ICESEE_SPACK_ROOT to the installation containing scripts/activate.sh and .spack-env/icesee" >&2; return 2;
  }
fi
if [[ -n "${ICESEE_SOURCE_ROOT:-}" ]]; then
  ICESEE_ROOT="$(cd "$ICESEE_SOURCE_ROOT" && pwd -P)"
elif [[ -f "$SPACK_REPO/applications/icepack_model/examples/idealized_pig/run_da_icepack.py" ]]; then
  ICESEE_ROOT="$SPACK_REPO"
else
  ICESEE_ROOT="$SPACK_REPO/ICESEE"
fi
PIG_DIR="$ICESEE_ROOT/applications/icepack_model/examples/idealized_pig"
[[ -f "$PIG_DIR/run_da_icepack.py" && -f "$PIG_DIR/params.yaml" && -d "$PIG_DIR/data" ]] || {
  echo "Incomplete PIG source/input layout: $PIG_DIR" >&2; return 2;
}
source "$SPACK_REPO/scripts/activate.sh" || return
# Verify package resolution without importing the scientific stack.
if [[ "$ICESEE_ROOT" == "$SPACK_REPO" ]]; then
  [[ -L "$SPACK_REPO/.pace-runtime/ICESEE" ]] || { echo "Missing standalone package alias" >&2; return 2; }
  export PYTHONPATH="$SPACK_REPO/.pace-runtime:$ICESEE_ROOT:${PYTHONPATH:-}"
else
  export PYTHONPATH="$(dirname "$ICESEE_ROOT"):$ICESEE_ROOT:${PYTHONPATH:-}"
fi
MPI_DIR="$(spack -e "$SPACK_REPO/.spack-env/icesee" location -i openmpi)" || return
export PATH="$MPI_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$MPI_DIR/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1
cd "$PIG_DIR" || return
