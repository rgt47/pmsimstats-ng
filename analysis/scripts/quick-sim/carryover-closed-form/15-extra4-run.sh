#!/bin/bash
# Recompute the 14-analysis decision on job 37 (--extra4) alone and pooled
# with jobs 15 and 21, then bootstrap the maximin shortfall. Run from the
# repo root. The large intermediate inputs go to a temporary work
# directory (EXTRA4_WORK, default a fresh mktemp directory) that is
# removed at the end; only the decision CSVs and the .txt logs are kept
# (.txt because *.out and *.log are gitignored),
# in analysis/data/quick-sim/carryover-closed-form/strawman-stress/
# extra4-decision/. Post hoc; see docs/47.
set -e
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
S=analysis/scripts/quick-sim/carryover-closed-form
OUT=analysis/data/quick-sim/carryover-closed-form/strawman-stress/extra4-decision
WORK=${EXTRA4_WORK:-$(mktemp -d)}
mkdir -p "$OUT" "$WORK"
Rscript $S/15-extra4-pooled-inputs.R "$WORK"
for s in new pooled; do
  Rscript $S/12-strawman-stress-decision.R --strawman DS --no-symmetry \
    --with-extras --tag -x --dir "$WORK/$s" > "$OUT/decision-$s.txt" 2>&1
  mkdir -p "$OUT/$s"
  cp "$WORK/$s"/decision-*.csv "$OUT/$s/"
done
Rscript $S/15-extra4-bootstrap.R "$WORK" > "$OUT/boot.txt" 2>&1
rm -rf "$WORK"
cat "$OUT/boot.txt"
