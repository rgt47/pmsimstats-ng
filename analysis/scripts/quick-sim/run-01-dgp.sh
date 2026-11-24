#!/usr/bin/env bash
set -e
# Resolve the repository root from this script's location so the wrapper
# runs from any checkout (it is three levels below the root).
cd "$(dirname "$0")/../../.."
exec Rscript analysis/scripts/quick-sim/01-dgp-prototype.R
