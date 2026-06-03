#!/usr/bin/env bash
set -euo pipefail

ICESEE_WORKDIR="$(pwd)"
ICESEE_ROOT="$(dirname "${ICESEE_WORKDIR}")"

export ICESEE_WORKDIR
export ICESEE_ROOT
export PYTHONPATH="${ICESEE_ROOT}:${PYTHONPATH:-}"

echo "ICESEE_WORKDIR=${ICESEE_WORKDIR}"
echo "ICESEE_ROOT=${ICESEE_ROOT}"
echo "PYTHONPATH=${PYTHONPATH}"

{
  echo "ICESEE_WORKDIR=${ICESEE_WORKDIR}"
  echo "ICESEE_ROOT=${ICESEE_ROOT}"
  echo "PYTHONPATH=${PYTHONPATH}"
} >> "${GITHUB_ENV}"

python - <<'PY'
import os
import sys
import importlib

root = os.environ["ICESEE_ROOT"]

if root not in sys.path:
    sys.path.insert(0, root)

print("="*80)
print("IMPORT TEST")
print("="*80)

for p in sys.path:
    print(p)

m = importlib.import_module("ICESEE")
print("SUCCESS: ICESEE ->", m.__file__)
PY