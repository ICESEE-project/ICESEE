#!/bin/bash
set -euo pipefail

python -m venv icesee-env
PY="./icesee-env/bin/python"

$PY -m pip install -U pip wheel setuptools

# Install deps (if you truly need requirements.txt)
$PY -m pip install -r requirements.txt

# Install ICESEE into the venv (this is the key)
$PY -m pip install -e .

# Sanity check
$PY -c "import ICESEE; import sys; print('ICESEE import OK from', sys.executable)"
echo "Venv ready: source icesee-env/bin/activate"