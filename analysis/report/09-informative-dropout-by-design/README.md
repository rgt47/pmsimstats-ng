# 09-informative-dropout-by-design

*2026-05-07 16:30 PDT*

## Working title

Informative dropout and trial-design choice in aggregated N-of-1
biomarker-validation trials.

## Origin

This manuscript directly addresses an unaddressed gap identified
during the deep read of @hendrickson2020. Hendrickson's Figure 4A
heatmap and §3.4 bias quantification analysis demonstrate that
the four candidate trial designs (open-label, open-label plus
blinded discontinuation, traditional crossover, and the proposed
hybrid N-of-1 design) differ substantially in their robustness
to informative dropout. The published methodological N-of-1
literature has not engaged this design-by-dropout coupling as a
central question, and none of the eight other manuscripts in the
pmsimstats-ng programme makes informative dropout its primary
focus. The present paper fills the gap.

## Scope

1. Formalise the informative-dropout model used in
   @hendrickson2020 (hazard depending on symptom worsening since
   baseline) and place it within the missing-data taxonomy
   (MCAR, MAR, MNAR).
2. Reproduce Hendrickson's Figure 4A finding that the four
   candidate designs differ in dropout-robustness, extending the
   analysis to (a) wider parameter ranges (longer carryover
   half-lives, higher effect sizes), (b) sensitivity to
   dropout-pattern misspecification at the analysis stage, and
   (c) the bias-quantification mode comparing the dropout-
   conditioned coefficient estimate to a gold-standard
   uncensored estimate.
3. Quantify the 'happy accident' randomization-path selection
   effect that Hendrickson identifies in OL+BDC and N-of-1
   designs: biased dropout preferentially preserves the most-
   informative crossover blocks because patients with
   inconclusive open-label response are more likely to continue
   to the later blocks.
4. Provide a design-by-dropout robustness recommendation table
   indexed by the trial-design family, the anticipated dropout
   rate, and the suspected dropout-pattern mechanism (balanced
   vs symptom-worsening-driven vs symptom-improvement-driven).

## Relationship to other papers in the series

- `01-dgp-mean-moderation-vs-mvn`: shares the MVN-covariance DGP.
  The present paper uses the same generative substrate and adds
  the informative-dropout layer.
- `02-carryover-sensitivity`: complementary; the present paper
  holds carryover at the prazosin-PTSD reference value and
  varies the dropout layer instead.
- `04-treatment-main-effect`: clinical companion that compares
  N-of-1 to parallel-group RCT for the main effect; mentions
  dropout but does not centre it.
- `05-nof1-design-sensitivity`: covers twelve parameter sweeps
  but does not include a dropout sweep; the present paper
  supplies that complementary axis.
- `08-test-procedure-design-sensitivity`: covers cycle-by-period
  design parameters within one design family. The present paper
  is design-family-comparative rather than within-family.

## Simulation substrate

The simulation programme uses the existing pmsimstats `R/censordata.R`
infrastructure (which implements the hazard-based dropout model
of @hendrickson2020) extended with explicit dropout-pattern
sweeps. The four design families (OL, OL+BDC, CO, hybrid N-of-1)
are taken directly from the @hendrickson2020 Figure 2
specification.

## Driver scripts

A driver directory at `analysis/scripts/09-informative-dropout-by-
design/` will be created when the simulation study is
implemented. It does not yet exist.

## Author

pmsimstats team
