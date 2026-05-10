#!/bin/bash
# tools/render.sh
#
# Render a manuscript and stage the resulting PDF for collaborator
# distribution. Wraps the R helper at tools/stage-render.R so that
# terminal-based rendering produces the same staged-PDF output as
# RStudio's Knit button (which routes through the same helper via
# the YAML knit: field).
#
# Usage:
#   bash tools/render.sh <paper-slug>
#   bash tools/render.sh <paper-slug> --short          # for 04-* short variant
#   bash tools/render.sh <relative-or-absolute-rmd-path>
#
# Examples:
#   bash tools/render.sh 06-component-decomposition
#   bash tools/render.sh 04-treatment-main-effect --short
#   bash tools/render.sh analysis/report/03-latent-class-mixture-application/report.Rmd

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ $# -lt 1 ]; then
  echo "usage: $0 <paper-slug> [--short]" >&2
  echo "       $0 <path-to-rmd>" >&2
  exit 1
fi

arg="$1"
shift

# Resolve the input Rmd path. Three accepted forms:
#   1. A bare slug (e.g. 06-component-decomposition) → reports/<slug>/report.Rmd
#   2. A bare slug with --short → reports/<slug>/report_short.Rmd
#   3. An explicit path to an .Rmd file (relative or absolute)
if [[ "$arg" == *.Rmd ]]; then
  if [[ "$arg" == /* ]]; then
    rmd="$arg"
  else
    rmd="$PROJ_ROOT/$arg"
  fi
else
  base="report.Rmd"
  if [ "${1:-}" = "--short" ]; then
    base="report_short.Rmd"
    shift
  fi
  rmd="$PROJ_ROOT/analysis/report/$arg/$base"
fi

if [ ! -f "$rmd" ]; then
  echo "error: input Rmd not found at $rmd" >&2
  exit 2
fi

echo "rendering: $rmd"
Rscript -e "stage <- source('$PROJ_ROOT/tools/stage-render.R')\$value; stage('$rmd')"
