## Validation gate {#sec-gate}

**RESULTS-PENDING.** The corrected grid has not finished. Run, from
the repository root:

`bash analysis/scripts/dgp-combined/run-combined.sh --methods e2 --cells all --reps 1000 --t1half 0,0.5,1 --workers 8`

then `Rscript analysis/scripts/dgp-combined/02-summarise-combined.R`,
and render again. The run is 9000 replicate units of nine weight
cells each. Each unit costs about 1.2 CPU-seconds on an idle machine,
about 25 to 40 minutes on eight workers. It resumes from checkpoints
if interrupted.

## Test calibration under the null

RESULTS-PENDING, as above.

## Power over the weight grid

RESULTS-PENDING, as above.

## The fitted coefficient {#sec-est}

RESULTS-PENDING, as above.

## Additivity of power {#sec-additivity}

RESULTS-PENDING, as above.

## The equal-total line and carryover {#sec-eqsum}

RESULTS-PENDING, as above.

## The G9 analysis and the pharmacological half-lives {#sec-g9}

RESULTS-PENDING, as above.
