# 12-nof1-design-efficiency

*2026-08-20*

## Working title

A closed-form approximation to the biomarker-treatment interaction
estimator in aggregated N-of-1 trials, and what it implies for
design efficiency under a mandatory open-label lead-in.

## Origin

This manuscript began as a design-consulting question: given that
the first phase of a trial must be open-label (for recruitment
reasons, a fixed constraint rather than a free design parameter),
what should the analyst prioritize among the remaining design
choices (design topology after the open-label phase, sample size
and its allocation across randomization paths, carryover half-life
and decay-shape robustness, assessment timing relative to
discontinuation, working-covariance specification) to maximize
power for the biomarker-treatment interaction test?

The first pass at an answer was a qualitative synthesis of the
eleven-paper compendium (see
`analysis/report/compendium-mathematics.pdf` and
`analysis/report/whitepaper-compendium-summary.md`). This
manuscript's contribution beyond that synthesis is a genuine
closed-form approximation to the sampling variance and power of the
interaction estimator, derived from first principles (a two-group
within-subject-contrast reduction of the analysis model, combined
with the exact variance of a mean of AR(1)-correlated, unequally
spaced observations) and validated against real Monte Carlo output
already on disk, not fabricated or asserted.

## Scope

1. Derive a closed-form approximation to $\mathrm{Var}(\hat\beta_{bm:D})$
   under the mean-moderation architecture, for the two-group
   (on-drug mean minus off-drug mean) reduction of the fixed analysis
   model, as a function of trial design (via the exact visit-time AR(1)
   sample-mean variance), total sample size $N$, and an AR(1)
   correlation parameter $\rho$.
2. Derive a companion closed-form carryover **attenuation factor**
   $A(t_{1/2})$, the fraction of the true interaction signal not
   washed out by residual on-drug exposure at nominally off-drug
   assessment times, as a function of the design's visit schedule and
   the carryover half-life.
3. Validate both against the existing Monte Carlo output of the
   dual-channel-architecture driver
   (`analysis/scripts/dgp-combined/01-run-combined-factorial.R`),
   whose $(c_{bm,A}=0.45, c_{bm,B}=0)$ boundary is exactly the pure
   mean-moderation architecture, at analysis specification E1
   (binary drug indicator, matching the two-group reduction). This
   validation is honest about where the approximation succeeds
   (the $N$- and design-dependence of the sampling variance, $R^2
   \approx 0.95$ with a single fitted scale parameter; power
   prediction MAE $\approx 0.05$, correlation $\approx 0.94$ against
   simulated power) and where it does not (the attenuation factor
   correctly *orders* the three designs by carryover vulnerability
   but substantially overstates the *magnitude* of attenuation in
   the actual fitted mixed-model coefficient, because the fitted
   model's explicit time covariate recovers information the crude
   two-group contrast discards).
4. Translate the validated parts of the closed form into concrete,
   quantitative design guidance for the open-label-first-phase
   scenario that motivated the question, superseding the purely
   qualitative version of that guidance.

## Relationship to the compendium

This is the twelfth manuscript in the pmsimstats-ng series and the
first to attempt a closed-form (rather than purely simulation-based)
treatment of the biomarker-treatment interaction estimand. Papers 04
and 08 established that closed forms exist for adjacent problems (the
treatment main effect, via Senn/Zucker; the dichotomized-biomarker
RM-ANOVA interaction test, Paper 08's own derivation) but neither
covers the continuous-biomarker, continuous-exposure ($D_{bc}$),
mixed-model estimand that is the compendium's actual estimand
throughout Papers 01, 02, 05, 06, 09, and 11. This manuscript is
explicitly scoped as a first, partial closed form for that estimand,
not a replacement for the compendium's simulation-based results.

## Status

Draft. `rgt` paragraphs are placeholders pending the author's own
prose, per the project's three-part paragraph convention.
