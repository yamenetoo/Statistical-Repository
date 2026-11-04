#!/usr/bin/env bash
# Rebuild all TikZ figures: PDF (vector) in this folder, PNG in ../../man/figures.
# Requires: pdflatex with TikZ (texlive-pictures) and pdftocairo (poppler-utils).
set -euo pipefail
cd "$(dirname "$0")"
for f in workflow bootstrap-flow arl-flow architecture chart-illustration; do
  pdflatex -interaction=nonstopmode -halt-on-error "$f.tex" > /dev/null
  pdftocairo -png -r 150 -singlefile "$f.pdf" "../../man/figures/$f"
  echo "built $f"
done
rm -f ./*.aux ./*.log
