# ADEMP pre-registration: 06-component-decomposition

*2026-05-07 16:50 PDT*

This document pre-registers the simulation programme that
underlies `analysis/report/06-component-decomposition/`. ADEMP
structure follows Morris, White and Crowther (2019).

## Common infrastructure

**Trial design.** Hendrickson hybrid OL+BDC design with three
phases: 8-week open-label active drug, 4-week blinded
discontinuation, 4-week crossover. Within-participant assessment
at the end of each phase plus weekly during the blinded
phase.

**Sample sizes.** $N \in \{35, 70, 100, 150\}$, covering the
trial-relevant range plus one larger anchor for asymptotic
trend extrapolation.

**Drug context.** Prazosin for PTSD nightmares, calibrated to
the Hendrickson 2020 reference parameters.

**DGP architecture.** Three-component additive
generative model $Y_{it} = \mathrm{BL}_i - (BR_{it} + PB_{it} +
TV_{it}) + \varepsilon_{it}$ with each component a
participant-specific modGompertz function and AR(1)
within-participant residual correlation. The PB component
includes the phase-dependent expectancy multiplier
$\eta(\text{phase}) \in \{1.0, 0.5, 0\}$ for open-label-on-drug,
blinded, and open-label-placebo phases respectively.

**Replicates.** $n_{\text{reps}} = 1000$ per cell for power and
bias estimands; $n_{\text{reps}} = 5000$ for type I error and
component-recovery cells where rare-event resolution is needed.

**Random-number control.** Per-replicate seed via cell-descriptor
hash; base seed in `01-base-seed.txt`.

**Reporting standard.** Morris-White-Crowther ADEMP with MCSE
columns for every reported performance measure.

## Study A. Bias of one-component analysis under three-component DGP

**A.** Quantify the bias and Monte Carlo standard error of the
biomarker-treatment interaction estimate when a one-component
linear-mixed analysis is fitted to data generated under the
three-component (BR-PB-TV) decomposition, as a function of the
underlying $PB$ and $TV$ strengths.

**D.** Three-component DGP at the prazosin-PTSD reference
$m_{BR}$ values, with $m_{PB}$ varied across $\{0, 1, 3, 6, 10\}$
and $m_{TV}$ varied across $\{-1, 0, 1, 2\}$, generating a 5x4
grid of population-mean component compositions.

**E.** Primary estimand: $\hat{\beta}_{bm:D}$ from the analysis-
model fixed-effects table. The bias is computed against the
known biomarker-by-$BR$ effect specified in the DGP.

**M.** Three analysis variants compared:

1. *One-component linear-mixed*: `Sx ~ bm + t + Dbc + bm:Dbc +
   (1 | ptID)` with `corCAR1(form = ~t|ptID)`.
2. *Phase-augmented linear-mixed*: `Sx ~ bm + t + Dbc + phase +
   Dbc:phase + bm:Dbc + bm:Dbc:phase + (1 | ptID)` with the
   same residual correlation.
3. *Full nonlinear three-component*: `nlme::nlme` fit of the
   modGompertz BR + PB + TV components with participant-specific
   random effects on each $m$. Fitted only at $N \geq 100$ where
   convergence rates exceed 80% in pilot.

**P.** Power for $H_0: \beta_{bm:D} = 0$, type I error under
$\beta_{bm:D} = 0$, bias of $\hat{\beta}_{bm:D}$, empirical SE
across replicates, model SE, nominal-95% CI coverage,
convergence rate. MCSE attached to every estimate.

## Study B. Identifiability of the full decomposition

**A.** Characterise the identifiability of the participant-
specific BR, PB, TV component estimates from the full nonlinear
analysis at trial-relevant $N$.

**D.** As Study A, restricted to the $m_{PB} \in \{1, 6\}$ and
$m_{TV} \in \{0, 1\}$ combinations.

**E.** Per-participant component estimates $(\hat{m}_{BR},
\hat{m}_{PB}, \hat{m}_{TV})$ and their joint posterior
covariance.

**M.** Full nonlinear three-component fit only.

**P.** Convergence rate, mean fit time, condition number of the
parameter-covariance matrix at convergence, and mean-squared
error of each component estimate against the participant's
true $m$ value. Reported as a function of $N$ to characterise
the sample-size threshold for usable identifiability.

## Study C. Subadditivity sensitivity (mechanistic vs emergent)

**A.** Quantify the bias in $\hat{\beta}_{bm:D}$ from the
phase-augmented linear-mixed analysis under (i) a true
mechanistic $BR \times PB$ subadditive interaction with
coefficient $\gamma$, and (ii) a latent-class DGP with
class-correlated $BR$ and $PB$ components but no direct
interaction.

**D.** Two DGP families:

1. *Mechanistic subadditivity DGP*: $Y_{it} = \mathrm{BL}_i -
   [BR + PB + TV + \gamma \cdot BR \cdot PB] + \varepsilon$
   with $\gamma \in \{-0.05, 0, 0.02, 0.05, 0.10\}$.
2. *Latent-class with class-correlated PB DGP*: as in the
   Study 1 mixture cell of paper 03 with `pb_class_correlation`
   axis varied.

**E.** $\hat{\beta}_{bm:D}$ from the phase-augmented analysis,
plus $\hat{\gamma}$ from a saturated mean-structure analysis
that includes the $BR \times PB$ interaction term explicitly.

**M.** Phase-augmented linear-mixed analysis (primary), saturated
mean-structure analysis (secondary), class-aware `lcmm`
analysis (tertiary, where applicable).

**P.** Power, bias, and MCSE under each DGP family. The
falsifiability test: under the latent-class DGP with no true
$\gamma$, the saturated mean-structure analysis should still
recover non-zero $\hat{\gamma}$, demonstrating that
$\hat{\gamma}$ alone is not evidence of mechanistic subadditivity.

## Sequencing

1. **Phase 1: Study A.** Establishes the bias of the simple
   analysis as a function of population PB and TV strength.
2. **Phase 2: Study B.** Establishes the sample-size threshold
   for usable identifiability of the full decomposition.
3. **Phase 3: Study C.** Resolves the mechanistic-vs-emergent
   subadditivity question using the latent-class infrastructure
   from paper 03.

## Deliverables

- Driver scripts at `analysis/scripts/component-decomposition/`.
- ADEMP results tables per study with MCSE columns alongside
  every estimate.
- Figures: bias-by-(PB strength, TV strength) heatmap; sample-
  size-by-identifiability curve; subadditivity sensitivity plot.

## Deviations log

`02-deviations-log.md` records departures from the cells,
estimands, methods, or performance measures pre-registered
above.

---

*Source:* `~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/component-decomposition/00-ademp-pre-registration.md`
