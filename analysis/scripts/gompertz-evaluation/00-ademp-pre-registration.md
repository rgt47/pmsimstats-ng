# ADEMP pre-registration: 07-gompertz-evaluation

*2026-05-07 16:55 PDT*

This document pre-registers the simulation programme that
underlies `analysis/report/07-gompertz-evaluation/`. ADEMP
structure follows Morris, White and Crowther (2019).

## Common infrastructure

**Trial design.** Hendrickson hybrid OL+BDC design, 16-week
total trial duration, 8-week open-label active drug + 4-week
blinded discontinuation + 4-week crossover.

**Sample sizes.** $N \in \{35, 70, 100\}$.

**Drug context.** Prazosin for PTSD nightmares, calibrated to
the Hendrickson 2020 reference parameters at the saturating
asymptote.

**Replicates.** $n_{\text{reps}} = 1000$ per cell for power and
bias estimands; $n_{\text{reps}} = 5000$ for type I error
cells.

**Random-number control.** Per-replicate seed via cell-descriptor
hash; base seed in `01-base-seed.txt`.

**Reporting standard.** Morris-White-Crowther ADEMP with MCSE
columns for every reported performance measure.

**Calibration to common effect size.** All four DGP families
share the same marginal Cohen's-style standardised effect size
at week 8 in the open-label phase, computed as the marginal
on-drug Cor(B, BR) in the saturated regime. Each family's
parameters are fitted to match this effect size before the
production grid runs, via the calibration sub-study described
in Study 0 below.

## Study 0. Effect-size calibration sub-study

**A.** Fit the parameters of each of the four DGP families
(Gompertz, symmetric logistic, hyperbolic-tangent, piecewise-
linear breakpoint) so that all four produce the same marginal
on-drug Cor(B, BR) at week 8.

**D.** Single-component BR-only generative model under each of
the four families, with population-level participant
heterogeneity but no PB or TV component (those are added in
Study 1).

**E.** The set of family-specific parameters that produces
matched marginal effect size.

**M.** Numerical optimisation via `optim`, with an analytic
Jacobian where available.

**P.** Convergence to the calibration target within $|d| < 0.01$
on Cohen's-style effect size; reproducibility across random
seeds.

## Study 1. Power, bias, and type I error across families

**A.** Quantify how the standard linear-mixed analysis's
inferential properties (power, type I error, bias) depend on
the DGP family used to generate the data.

**D.** Four DGP families on the BR factor only:

1. *Modified Gompertz* (the pmsimstats default): $BR(t) =
   m \exp(-d \exp(-rt))$ with vertical-offset adjustment.
2. *Symmetric logistic*: $BR(t) = m / (1 + \exp(-r(t - t_0)))$.
3. *Hyperbolic-tangent*: $BR(t) = (m/2)(1 + \tanh(r(t - t_0)))$.
4. *Piecewise-linear breakpoint*: $BR(t) = m \cdot \min(1, t/t_0)$.

PB and TV components are held at their Hendrickson reference
specifications (modGompertz with calibrated parameters), so the
family axis varies BR alone. Each family is calibrated to the
matched marginal effect size from Study 0.

**E.** Primary estimand: $\hat{\beta}_{bm:D}$ from the linear-
mixed analysis. Secondary estimand: posterior class entropy is
not applicable here (single-component DGPs), but the empirical
SD of $\hat{\beta}_{bm:D}$ across replicates is reported as a
finite-sample efficiency measure.

**M.** Single analysis: linear mixed-effects model
`Sx ~ bm + t + Dbc + bm:Dbc + (1 | ptID)` with
`corCAR1(form = ~t|ptID)`.

**P.** Power for $H_0: \beta_{bm:D} = 0$, type I error under
$\beta_{bm:D} = 0$, bias of $\hat{\beta}_{bm:D}$, empirical SE
across replicates, model SE, nominal-95% CI coverage,
convergence rate. MCSE attached to every estimate.

**Pre-registered prediction.** Power varies less than $\pm 5\%$
across the four families at fixed sample size and effect size,
provided each is calibrated to the matched marginal effect
size. Type I error remains within $\pm 1\%$ of nominal across
families. Bias is small at all four families when the analysis
model's continuous-decay carryover specification matches the
DGP's saturating-shape regime.

**Total cells.** $4 \text{ families} \times 3 N \times 3
\text{ effect sizes} \times 4 t_{1/2,\text{DGP}} = 144$ cells.

## Study 2. Sloppiness and parameter identifiability

**A.** Characterise the parameter-redundancy structure
(sloppiness) of the modified Gompertz form fitted to short
N-of-1 trajectories.

**D.** Single Gompertz BR generative model at the Hendrickson
reference; eigenvalue analysis of the parameter Hessian at
the maximum-likelihood point, computed analytically.

**E.** The set of eigenvalues of the parameter Hessian, ordered
by magnitude. Following @jagadeesan2023sloppiness, the
sloppiness is characterised by the ratio of the largest to
the smallest eigenvalue.

**M.** Numerical Hessian computation at the MLE for each
replicate; analytical comparison to the unified-Richards
reparameterisation of @tjorve2010richards.

**P.** Mean and variance of each eigenvalue across replicates;
the sloppiness ratio; the sample-size scaling of the smallest
eigenvalue.

## Deliverables

- Driver scripts at `analysis/scripts/gompertz-evaluation/`.
- ADEMP results tables per study with MCSE columns.
- Figures: power/bias/type-I plots faceted by DGP family;
  sloppiness eigenvalue spectrum.

## Deviations log

`02-deviations-log.md` records departures from the cells,
estimands, methods, or performance measures pre-registered
above.

---

*Source:* `~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/gompertz-evaluation/00-ademp-pre-registration.md`
