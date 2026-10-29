# 07-gompertz-evaluation

*2026-10-28 18:01 PDT*

## Working title

Response-curve family and the biomarker-treatment interaction test
in N-of-1 trial simulation: a common-random-numbers comparison of
the modified Gompertz template with three alternatives.

## Origin

This manuscript evaluates the modified-Gompertz parametric template
used throughout pmsimstats to generate the three latent response
components (BR, PB, TV). The technical companion at
`docs/legacy/25-gompertz-model-evaluation.md` provides the
historical and biological context, the parameter-range
characterization, and the alternative-family survey on which this
manuscript is built.

## Scope

1. Characterize the modified Gompertz family in the form used by
   pmsimstats, in NOTATION.md symbols (m, d, r), with the
   vertical-offset adjustment and the asymmetry stated correctly
   (short acceleration, long deceleration tail).
2. Present three alternative families (symmetric logistic,
   hyperbolic tangent, piecewise-linear breakpoint) at matched
   asymptote and week-5 anchor, and show what the calibration
   leaves at the sampled visits: the logistic equals the Gompertz
   curve to within 0.12 points, the hyperbolic tangent is the
   logistic at twice the rate.
3. Run the four families through the standard pmsimstats linear
   mixed analysis under common random numbers, in both interaction
   encodings, and report paired power contrasts with paired MCSE,
   Type I error from pooled null cells, and bias against the
   closed-form reference -c_bm sigma_BR / sigma_bm.
4. State the design conditionality (the two-path OL+BDC schedule
   samples the rising phase at two visits) and the conditions that
   would change the answer.

## Status (2026-10-28)

Rewritten in response to `referee-report-2026-10-28.md`. The
earlier headline (a 0.039 power spread with Gompertz lowest under
mean moderation; an exactly inert covariance architecture) was
Monte Carlo noise from unpaired seeds and is retracted in the
manuscript. The supported result: families that coincide at the
sampled visits give identical decisions; the hyperbolic tangent
lowers power by about 0.01 under both encodings; the covariance
encoding is nearly but not exactly inert.

## Relationship to other papers in the series

- `01-dgp-mean-moderation-vs-mvn`: shares the Gompertz
  response-curve family and the two interaction encodings.
- `02-carryover-sensitivity`: uses Gompertz for the BR component;
  the present paper shows the family contributes little
  sensitivity at the reference design.
- `10-interaction-test-calibration`: studies the size of the
  interaction test, which the present paper finds below nominal (0.034
  against 0.05).

## Driver scripts

`analysis/scripts/gompertz-evaluation/`:

- `00-ademp-protocol.md`: working protocol (exploratory).
- `02-deviations-log.md`: departures from the protocol.
- `03-common-random-numbers-sweep.R`: production driver (about 8
  minutes on 8 cores). Outputs `analysis/data/quick-sim/
  07-gompertz-crn-{replicates,paired,meta}.rds` and
  `07-gompertz-crn-summary.txt`.
- `02-faithful-trajectory-sweep.R`: superseded unpaired driver,
  kept for the comparison in the manuscript's Section 5.1.
- `01-architecture-a-trajectory-sweep.R`: prototype, superseded.

## Author

pmsimstats team
