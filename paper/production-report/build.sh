#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
if grep -Eq 'undefined references|Citation .* undefined|Reference .* undefined' build/main.log; then
  echo 'Unresolved manuscript references remain.' >&2
  exit 1
fi
echo 'Built build/main.pdf; visually inspect the final PDF before distribution.'
