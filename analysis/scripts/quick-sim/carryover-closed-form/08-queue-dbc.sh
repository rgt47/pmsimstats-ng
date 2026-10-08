#!/bin/bash
# Start the Dbc arms of 08-covariance-study.R once the follow-up run has
# finished, so the two runs do not compete for cores. Run from anywhere:
#   nohup bash analysis/scripts/quick-sim/carryover-closed-form/08-queue-dbc.sh &
cd "$(dirname "$0")/../../../.."
out=analysis/data/quick-sim/carryover-closed-form/covariance-study
while pgrep -f "08-covariance-study.R.*--arms followup" > /dev/null; do
  sleep 120
done
echo "$(date '+%Y-%m-%d %H:%M:%S') follow-up finished; starting Dbc arms" > "$out/queue-dbc.log"
Rscript analysis/scripts/quick-sim/carryover-closed-form/08-covariance-study.R \
  --reps-alt 250 --reps-null 500 --chunk 50 --cores 8 --arms dbc \
  > "$out/run-alt250-null500-dbc.log" 2>&1
echo "$(date '+%Y-%m-%d %H:%M:%S') Dbc arms exited with status $?" >> "$out/queue-dbc.log"
