# ADEMP pre-registration: 04-treatment-main-effect

*2026-05-07 16:45 PDT*

This document pre-registers the simulation programme that
underlies `analysis/report/04-treatment-main-effect/` (full and
short variants). ADEMP structure follows Morris, White and
Crowther (2019). Each study fixes Aims (A), Data-generating
mechanisms (D), Estimands (E), Methods (M), and Performance
measures (P) before any production runs are executed.

## Common infrastructure

**Sample sizes.** Aggregated N-of-1 trial: $N = 35$ participants
(per-design count, matching the Hendrickson 2020 minimum and
the per-participant convention adopted in
Blackston et al. 2019). Parallel-group RCT: $N = 35$ at 1:1
randomisation as the matched-N comparator, plus a sensitivity
row at $N = 70$ (1:2 RCT-to-N-of-1 patient ratio).

**Trial design (N-of-1 arm).** Hendrickson hybrid design:
8 weeks open-label active drug + 4 weeks blinded
discontinuation + 4 weeks crossover-1 + 4 weeks crossover-2,
with assessment timepoints at the end of each period.

**Trial design (RCT arm).** 16-week parallel-group double-blind
RCT with assessment at weeks 0, 4, 8, 12, 16.

**Drug context.** Prazosin for PTSD nightmares, calibrated to
the Raskind 2013 PG-RCT and the @hendrickson2020 N-of-1 reference.

**Carryover.** Exponential decay with half-life $t_{1/2} = 3$ days
at the baseline cell, varied across $t_{1/2} \in \{0.5, 1, 3, 5,
7\}$ days in the carryover sweep.

**Replicates.** $n_{\text{reps}} = 2000$ per cell, applied
uniformly to power, bias, Type I error, and all sensitivity
sweep cells. The pre-specified target is MCSE $\leq 1.1\%$ on
a power estimate near $0.75$ and MCSE $\leq 0.5\%$ on a Type I
estimate near $0.05$; $n_{\text{reps}} = 2000$ satisfies both.

**Random-number control.** `RNGkind("L'Ecuyer-CMRG")` set once
at program entry; `set.seed(2024L)` immediately after. Per-
replicate `.Random.seed` state is captured before each replicate
so that any anomalous result can be reproduced exactly. No
`set.seed()` calls appear inside simulation functions.

**Reporting standard.** Morris, White and Crowther (2019) ADEMP
with Monte Carlo standard errors for every reported performance
measure.

## Study M1 (primary). Power comparison at the prazosin baseline

**A.** Quantify the difference in statistical power for detecting
the treatment main effect between the aggregated N-of-1 hybrid
design and the matched parallel-group RCT, at the
prazosin-PTSD-calibrated baseline parameters.

**D.** Mean-moderation Architecture A DGP with three response
components (BR, PB, TV); true effect size $-2.0$
nightmares/week; carryover half-life 3 days; AR(1) within-
participant residual correlation $\rho = 0.3$;
between-participant random-intercept SD $\tau$ calibrated to
yield $I^2 \approx 50\%$ in the parallel-group analysis.

**E.** Primary estimand: $\hat{\beta}_D$ in the analysis-model
fixed-effects table (the average treatment effect coefficient).
Secondary estimand: power $\pi = \Pr(p < 0.05)$ for the
two-sided $H_0: \beta_D = 0$.

**M.** N-of-1 arm: linear mixed-effects model
`Sx ~ Dbc + t + (1 | ptID)` with `corCAR1(form = ~t|ptID)`.
RCT arm: ANCOVA on endpoint with baseline severity covariate.

**P.** Power, type I error under $\beta_D = 0$, bias of
$\hat{\beta}_D$, empirical SE of $\hat{\beta}_D$ across
replicates, mean within-replicate model SE, nominal-95% CI
coverage, convergence rate. Each reported with the binomial-
proportion MCSE
$\sqrt{\hat{\pi}(1-\hat{\pi})/n_{\text{reps}}}$ for proportions
and the Monte Carlo SE of the within-replicate mean estimator
for continuous performance measures.

## Studies S1-S12 (sensitivity sweeps)

Twelve sensitivity sweeps perturb one axis at a time around the
M1 baseline. Each sweep fixes its own Aims, perturbed Data-
generating mechanism axis, inherited Estimands, inherited
Methods, and inherited Performance measures.

| Sweep | Axis | Range | Cells |
|---|---|---|---|
| S1 | Effect size $\beta_D$ | $\{-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0\}$ | 7 |
| S2 | Carryover $t_{1/2}$ | $\{0.5, 1, 3, 5, 7\}$ days | 5 |
| S3 | Decay-form mis-specification | linear / exponential / Weibull | 3 |
| S4 | Between-patient SD $\tau$ | $\tau \in \{0, 0.5, 1, 1.5, 2\}$ | 5 |
| S5 | Within-patient SD $\sigma$ | $\sigma \in \{0.5, 1, 1.5, 2\}$ | 4 |
| S6 | Cycles per patient | $k \in \{1, 2, 3, 4, 6\}$ | 5 |
| S7 | Total participant count | $N \in \{15, 25, 35, 50, 70, 100\}$ | 6 |
| S8 | AR(1) $\rho$ | $\rho \in \{0, 0.15, 0.3, 0.45, 0.6\}$ | 5 |
| S9 | Period length | $L \in \{1, 2, 3, 4, 6\}$ weeks | 5 |
| S10 | Biomarker interaction $c_{bm}$ × DGP architecture | $\{0, 0.3, 0.6\}$ × $\{A, B\}$ | 6 |
| S11 | Carryover-by-AR(1) misspecification factorial | DGP $t_{1/2}$ × analysis $t_{1/2}$ | 16 |
| S12 | High-replicate null calibration | $\beta_D = 0$, $n_{\text{reps}} = 5000$ | 1 |

**P (sensitivity-sweep specific).** Each sweep reports power,
bias, MCSE, and convergence rate at every cell. Results are
plotted as a function of the perturbed axis with MCSE error
bars.

## Deliverables

- `analysis/scripts/treatment-main-effect/` driver scripts
  (already in place; documented in
  `analysis/scripts/treatment-main-effect/README.md`).
- ADEMP results table per sweep, with MCSE columns alongside
  every estimate.
- Figures: power-by-axis plot per sweep, plus the joint
  power-by-(N, effect-size) heatmap from M1.

## Deviations from this pre-registration

Any deviation from the cells, estimands, methods, or
performance measures pre-registered above is recorded in
`02-deviations-log.md` with date, reason, and reviewer signoff
before the affected production runs are committed.

---

*Source:* `~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/treatment-main-effect/00-ademp-pre-registration.md`
