#!/bin/bash
## Fetch the Hendrickson et al. (2020) pmsimstats sources and published
## results at commit 58b32a9 ('orig') into a local cache. Requires an
## authenticated GitHub CLI. Run from the repository root:
##   bash analysis/scripts/quick-sim/hendrickson-problems/00-fetch-58b32a9.sh
set -euo pipefail
sha=58b32a92b99a8beb145de00c7bd8a0880f06f0f0
out=analysis/data/quick-sim/hendrickson-58b32a9
mkdir -p "$out/R" "$out/data" "$out/vignettes"
for f in R/utilities.R R/buildtrialdesign.R R/generateData.R \
         R/generateSimulatedResults.R R/lme_analysis.R \
         R/censordata.R R/plottingfunctions.R \
         vignettes/Produce_Publication_Results_1_generate_data.Rmd \
         vignettes/Produce_Publication_Results_2_Look_at_data.Rmd \
         data/results_core.rda; do
  gh api "repos/rchendrickson/pmsimstats/contents/$f?ref=$sha" \
    -H 'Accept: application/vnd.github.raw' > "$out/$f"
  echo "fetched $f ($(wc -c < "$out/$f") bytes)"
done
echo "$sha" > "$out/COMMIT"
