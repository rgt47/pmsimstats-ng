#!/usr/bin/env bash
# pre-commit-parity.sh
#
# Git pre-commit hook that runs the cross-implementation parity test in
# quick mode whenever a staged change touches either the tidyverse or
# original-extended implementation. Blocks the commit on any unexpected
# parity regression.
#
# Install:
#   ln -s ../../scripts/pre-commit-parity.sh .git/hooks/pre-commit
# or
#   cp scripts/pre-commit-parity.sh .git/hooks/pre-commit
#   chmod +x .git/hooks/pre-commit
#
# Bypass (use sparingly, only when intentionally landing divergent code
# that has been annotated in implementations/parity_divergences.csv):
#   git commit --no-verify

set -euo pipefail

# Only trigger if relevant files are staged for commit.
CHANGED=$(git diff --cached --name-only --diff-filter=ACMR | \
  grep -E '^implementations/(tidyverse|original-extended)/' || true)

if [[ -z "$CHANGED" ]]; then
  exit 0
fi

echo "[parity-hook] changes detected in tracked implementations:"
echo "$CHANGED" | sed 's/^/  /'
echo "[parity-hook] running parity test (--quick)..."

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[parity-hook] Rscript not on PATH; skipping parity test." >&2
  echo "[parity-hook] install R and retry, or commit with --no-verify." >&2
  exit 1
fi

if Rscript implementations/test-parity-extended-tidyverse.R --quick; then
  echo "[parity-hook] parity preserved; proceeding with commit."
  exit 0
else
  echo "" >&2
  echo "[parity-hook] UNEXPECTED PARITY REGRESSION." >&2
  echo "[parity-hook] if this is intentional, annotate the cells in" >&2
  echo "  implementations/parity_divergences.csv" >&2
  echo "[parity-hook] or override with: git commit --no-verify" >&2
  exit 1
fi
