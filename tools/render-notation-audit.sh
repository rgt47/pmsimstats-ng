#!/bin/bash
# Render the notation audit with the canonical reference bound in as its
# appendix, so the two can never drift apart.
#
#   bash tools/render-notation-audit.sh
set -euo pipefail

cd "$(dirname "$0")/.."
R=analysis/report

pandoc "$R/whitepaper-notation-audit.md" "$R/NOTATION.md" \
  -o "$R/whitepaper-notation-audit.pdf" \
  --pdf-engine=xelatex -V monofont="DejaVu Sans Mono"

echo "wrote $R/whitepaper-notation-audit.pdf"
