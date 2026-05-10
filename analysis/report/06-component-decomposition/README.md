# 06-component-decomposition

*2026-05-07 15:30 PDT*

## Working title

Three-component decomposition of treatment response in
aggregated N-of-1 trials: pharmacological, expectation-driven,
and natural-history components.

## Origin

This manuscript is the journal-format presentation of the
three-component (BR-PB-TV) decomposition framework used
throughout the pmsimstats simulation programme. The pedagogical
companion at `docs/24-component-decomposition-pedagogy.md`
develops the same framework with worked examples and an extended
intuition narrative for readers approaching the material for the
first time. The present manuscript is the methodology paper:
formal model specification, identifiability conditions,
linear-mixed approximations, and a simulation study quantifying
the inferential cost of failing to decompose.

## Scope

1. Formalise the three-component model: $Y_{it} = \mathrm{BL}_i -
   (BR_{it} + PB_{it} + TV_{it}) + \varepsilon_{it}$ with
   participant-specific Gompertz components.
2. Map identifiability requirements onto trial-design contrasts
   (open-label, blinded discontinuation, crossover).
3. Develop linear-mixed-effects approximations to the full
   nonlinear model that retain the pharmacological-component
   target estimand at trial-relevant sample sizes.
4. Quantify, through Monte Carlo simulation calibrated to the
   prazosin-PTSD reference parameters, the bias in the lumped
   one-component analysis as a function of the underlying PB and
   TV strengths.
5. Apply the framework to predictive-biomarker validation: derive
   the conditions under which a candidate biomarker can be
   validated as a pharmacological predictor (as opposed to a
   placebo-or-natural-history predictor).

## Relationship to other papers in the series

- `01-dgp-mean-moderation-vs-mvn`: shares the Gompertz
  parameterisation but addresses biomarker-moderation encoding
  rather than within-trajectory decomposition.
- `02-carryover-sensitivity`: complementary to the present
  paper's PB-component identification through blinded-phase
  contrasts; the carryover paper concerns mis-specification of
  the residual drug effect during off-drug timepoints.
- `04-treatment-main-effect`: clinical companion that uses the
  three-component DGP to compare N-of-1 and parallel-group
  designs for the main treatment effect.
- `06` (this paper) is methodologically foundational for `04`,
  `05`, and the within-component biomarker analyses in `01`,
  `02`, and `03`.

## Driver scripts

A driver directory at `analysis/scripts/06-component-
decomposition/` will be created when the simulation study
underlying §5 of the manuscript is implemented. It does not yet
exist.

## Author

pmsimstats team
