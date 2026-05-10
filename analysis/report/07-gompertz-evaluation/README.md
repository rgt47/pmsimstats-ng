# 07-gompertz-evaluation

*2026-05-07 15:30 PDT*

## Working title

Critical evaluation of the modified Gompertz response function
as a parametric template for treatment-response trajectories in
N-of-1 trial simulation.

## Origin

This manuscript is the methodological evaluation of the
modified-Gompertz parametric template used throughout pmsimstats
to generate the three latent response components (BR, PB, TV).
The technical companion at `docs/25-gompertz-model-evaluation.md`
provides the historical and biological context, the parameter-
range characterisation, and the alternative-family survey on
which this manuscript is built.

## Scope

1. Characterise the modified Gompertz family in the form used by
   pmsimstats, including the vertical-offset adjustment and the
   parameter regions corresponding to the calibration sets used
   across published pmsimstats simulation studies.
2. Survey three alternative response-curve families (symmetric
   logistic, hyperbolic-tangent, piecewise-linear breakpoint)
   and present them at matched marginal effect size against the
   Gompertz reference.
3. Conduct a Monte Carlo sensitivity-analysis simulation study
   that crosses the four families (Gompertz plus three
   alternatives) on the DGP side against the standard pmsimstats
   linear-mixed analysis on the analysis side. Report power, type
   I error, and bias of the biomarker-treatment interaction
   estimate for each cell.
4. Identify the trajectory regimes in which Gompertz
   misspecification meaningfully affects inference (non-saturating
   growth, biphasic patterns, breakpoint behaviour).
5. Recommend a default sensitivity grid for pmsimstats users that
   reports Gompertz-vs-alternative comparison alongside primary
   simulation cells where the simulation results inform protocol-
   design decisions.

## Relationship to other papers in the series

- `01-dgp-mean-moderation-vs-mvn`: shares the Gompertz response-
  curve family. Conclusions about Gompertz robustness in the
  present paper apply directly to that paper's simulations.
- `02-carryover-sensitivity`: independent of response-curve
  family choice but uses Gompertz for the BR component; the
  present paper's robustness analysis quantifies how much of the
  carryover-paper conclusions are conditional on the Gompertz
  assumption.
- The pedagogical companion `docs/24-component-decomposition-
  pedagogy.md` describes the role of the Gompertz family in the
  decomposition framework but does not evaluate alternatives;
  the present paper supplies that evaluation.

## Driver scripts

A driver directory at `analysis/scripts/07-gompertz-evaluation/`
will be created when the sensitivity-analysis simulation is
implemented. It does not yet exist.

## Author

pmsimstats team
