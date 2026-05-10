# ADEMP pre-registration: 05-nof1-design-sensitivity

*2026-05-07 16:45 PDT*

This document pre-registers the simulation programme that
underlies `analysis/report/05-nof1-design-sensitivity/`. ADEMP
structure follows Morris, White and Crowther (2019).

## Common infrastructure

**Trial design.** Hendrickson hybrid OL+BDC N-of-1 design as
the primary aggregated variant, with the alternating A/B/A/B
N-of-1 variant retained as a sensitivity comparator at sweeps
that vary cycle count or period length. Two-period crossover
and parallel-group RCT included as cross-design comparators in
sweeps S6, S7, S8.

**Sample size.** $N = 35$ participants per design (Hendrickson
2020 minimum; per-participant parity with the parallel RCT
comparator), unless explicitly varied by the sweep.

**Drug context.** Prazosin for PTSD nightmares; effect size
$-2.0$ nightmares/week at the baseline cell.

**Carryover.** Exponential decay with half-life $t_{1/2} = 3$
days at the baseline cell.

**Replicates.** $n_{\text{reps}} = 500$ per cell for sweeps
S1-S11; $n_{\text{reps}} = 2000$ per cell for S12 (null
calibration).

**Random-number control.** Per-replicate seed derived from cell
descriptor hash and base seed; recorded in `01-base-seed.txt`.

**Primary estimand.** The treatment main-effect coefficient
$\beta_D$ in the linear mixed-effects analysis-model formula
$\mathrm{Sx} \sim \mathrm{Dbc} + t + (1\mid\mathrm{ptID})$ with
`corCAR1(form = ~t|ptID)`. The biomarker-by-treatment
interaction $\beta_{bm:D}$ is a secondary estimand in sweep S10
only.

**Performance measures.** Power $\pi = \Pr(p < 0.05)$ for the
two-sided $H_0: \beta_D = 0$, type I error rate under $\beta_D
= 0$, bias of $\hat{\beta}_D$, empirical Monte Carlo SE of
$\hat{\beta}_D$ across replicates, mean within-replicate model
SE, nominal-95% CI coverage, and convergence rate. Each
reported with the binomial-proportion MCSE
$\sqrt{\hat{\pi}(1-\hat{\pi})/n_{\text{reps}}}$ for proportions
and the Monte Carlo SE of the within-replicate mean estimator
for continuous measures.

## Sweeps S1-S12

Per the manuscript's Table 1 (Section 3 'Simulation design'):

| Sweep | Axis | Cells | $n_{\text{reps}}$ |
|---|---|---|---|
| S1 | Effect size | $\beta_D \in \{-3, -2.5, -2, -1.5, -1, -0.5, 0\}$ | 500 |
| S2 | Carryover half-life | $t_{1/2} \in \{0.5, 1, 3, 5, 7\}$ days | 500 |
| S3 | Decay-form misspecification | linear / exponential / Weibull | 500 |
| S4 | Between-patient SD | $\tau \in \{0, 0.5, 1, 1.5, 2\}$ | 500 |
| S5 | Within-patient SD | $\sigma \in \{0.5, 1, 1.5, 2\}$ | 500 |
| S6 | Cycles per patient | $k \in \{1, 2, 3, 4, 6\}$ | 500 |
| S7 | Total participant count | $N \in \{15, 25, 35, 50, 70, 100\}$ | 500 |
| S8 | AR(1) serial correlation | $\rho \in \{0, 0.15, 0.3, 0.45, 0.6\}$ | 500 |
| S9 | Period length | $L \in \{1, 2, 3, 4, 6\}$ weeks | 500 |
| S10 | Biomarker $c_{bm}$ × DGP architecture | $\{0, 0.3, 0.6\}$ × $\{A, B\}$ | 500 |
| S11 | Carryover-by-AR(1) misspecification | DGP $t_{1/2}$ × analysis $t_{1/2}$, 4×4 | 500 |
| S12 | High-replicate null calibration | $\beta_D = 0$, single cell | 2000 |

For each sweep, the **Aims** are to characterise the
sensitivity of the primary estimand to the perturbed axis;
**Data-generating mechanism** is the baseline DGP with the
named axis varied; **Estimands** and **Methods** are inherited
from the common infrastructure above; **Performance measures**
are inherited and supplemented by the sweep-specific
visualisation (power-by-axis plot with MCSE error bars).

## Deliverables

- Driver scripts at `analysis/scripts/nof1-design-sensitivity/`.
- `02-grid-summary.rds` and `02-sensitivity-summary.rds`
  checkpointed under `analysis/data/`.
- ADEMP results table per sweep with MCSE columns alongside
  every estimate.
- Figures plotting power and bias as a function of each axis
  with MCSE error bars.

## Deviations log

`02-deviations-log.md` records any departure from the cells,
estimands, methods, or performance measures pre-registered
above.

---

*Source:* `~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/nof1-design-sensitivity/00-ademp-pre-registration.md`
