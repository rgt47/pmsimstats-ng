# ADEMP protocol: 05-nof1-design-sensitivity
*2026-10-28 17:36 PDT*

**Author.** pmsimstats team

This is the working protocol for the exploratory simulation study
behind `analysis/report/05-nof1-design-sensitivity/`. It replaces the
protocol of 2026-05-07, whose grids matched none of the drivers. The
study is exploratory. No analysis plan was fixed in advance, and the
findings need confirmation in an independent study. The ADEMP
structure follows Morris, White and Crowther (2019).

## Aims

1. Describe the power of the Hendrickson et al. (2020) Hybrid design
   for the treatment main effect as one factor at a time departs from
   a baseline.
2. Find the regimes in which the Hybrid has more or less power than
   an alternating design, a two-period crossover and a parallel group
   on a common estimand.
3. Report Type I error, bias and standard-error calibration for every
   design.

## Data-generating mechanisms

- Outcome $Y_{it} = \mathrm{BL}_i - (TV_{it} + PB_{it} + BR_{it})$.
  Modified Gompertz means (displacement 1, rate 0.3 per week), with
  asymptotes 0.5 (TV), 1 (PB) and a per-design value for BR.
  Component SDs 0.5, 0.5, 1. Baseline score mean 8, SD 1.5.
  Biomarker mean 120, SD 15.
- AR(1) correlation $\rho = 0.3$ per week within components, 0.1
  between components at the same visit, 0.05 times the AR(1) factor
  at other visits. $c_{bm} = 0.3$.
- Carryover on the BR mean, exponential, half-life 3 days
  (3/7 week), with the original (compounding) rule of `docs/53`.
- Participant-level heterogeneity of the drug effect, SD $\tau$
  (baseline 0).
- $N = 35$ in total, split over the paths by the engine rule.
- Designs, all with eight follow-up visits: Hybrid (published,
  weeks 4, 8, 9, 10, 11, 12, 16, 20, four paths), Alt (alternating,
  two paths, same weeks), XO (two periods, weeks 10 and 20), RCT
  (two arms, same weeks, expectancy 0.5), RCT-OL (expectancy 1) and
  Compact2 (the two-path weekly schedule of the 2026-05 runs, S1
  only). Expectancy is 0.5 at blinded visits and 1 at open-label
  visits.
- **Common estimand.** The drug-response asymptote is set per design
  so that $\Delta_{\mathrm{on}}$, the mean drug response over the
  on-drug follow-up visits, equals $-0.5$ outcome units (S1 sweeps
  it). Every artifact stores `delta_on`.
- A cell whose covariance matrix is not positive definite is dropped,
  never repaired (`stop_on_nonpd = TRUE`). The artifact stores the
  dropped cells with the smallest eigenvalue.

## Estimands

The treatment main effect $\beta_D$: the coefficient of the exposure
indicator for the within-subject designs, and the arm coefficient of
the ANCOVA for the parallel group. The target of bias and coverage
is its large-sample value (`ETbeta`, the coefficient of a fit to all
simulated participants of the cell). The difference from
$\Delta_{\mathrm{on}}$ is reported as well.

## Methods

- Within-subject designs: `nlme::lme`, random intercept per
  participant, `corCAR1(~t | ptID)` residuals, fixed part
  `Sx ~ bm + t + Dbc` plus `De` where expectancy varies.
- Parallel groups: follow-up mean regressed on arm and baseline
  score by least squares (`parallel_analysis = 'ancova'`). The
  mixed-model analysis of 2026-05 is kept as `'legacy'` and checked
  in script 14.
- Assumed half-life of the exposure indicator: 0.2 week in all
  sweeps except S2 (0, 0.2), S3 (0, 0.2, 3/7), S11 (0, 0.2, 3/7) and
  S12 (0, 0.2).

## Performance measures

Rejection rate at 0.05, bias against `ETbeta`, bias against
$\Delta_{\mathrm{on}}$, empirical SD, mean model SE, $\kappa$,
coverage of the 95% interval, share of failed fits, each with its
Monte Carlo standard error. A Type I error is nominal when it lies
within two MCSE of 0.05 (MCSE at 0.05).

## Sweeps

| Sweep | Axis | Grid | Designs | Reps |
|---|---|---|---|---|
| S1 | $\Delta_{\mathrm{on}}$ | 0, -0.25, -0.5, -0.75, -1, -1.5, -2 | Hybrid, Compact2, Alt, XO, RCT, RCT-OL | 1000 |
| S2 | DGP half-life (days) by assumed half-life | 0, 0.5, 1, 2, 3, 5, 7, 10, 14; 0, 0.2 wk | Hybrid, Alt, XO, RCT | 500 |
| S3 | decay form by assumed half-life | exponential, linear, Weibull 0.8, Weibull 1.5, power; 0, 0.2, 3/7 wk | Hybrid, Alt, RCT | 500 |
| S4 | $\tau$; baseline SD | 0, 0.25, 0.5, 0.75, 1; 0.5, 1, 1.5, 2, 3 | Hybrid, Alt, RCT | 500 |
| S5 | multiplier on the three component SDs | 0.5, 0.75, 1, 1.25, 1.5 | Hybrid, Alt, RCT | 500 |
| S6 | periods in 16 weekly visits | 2, 4, 8, 16 (RCT: 16 visits) | Alt, RCT | 500 |
| S7 | total $N$ | 16, 24, 35, 60, 80, 120 | Hybrid, Alt, RCT | 500 |
| S8 | AR(1) $\rho$ | 0, 0.3, 0.5, 0.7, 0.9 | Hybrid, Alt, RCT | 500 |
| S9 | period (days) by half-life (days) | 3, 5, 7, 10, 14; 0.5, 3, 7 | Alt | 500 |
| S10 | $c_{bm}$ by encoding by half-life | 0, 0.2, 0.4, 0.6; mvn, mean_moderation; 0, 3 days | Hybrid, Alt, RCT | 500 |
| S11 | DGP carryover by DGP AR(1) by effect, assumed half-life | off/3 d; off/0.7; null/-0.5; 0, 0.2, 3/7 wk | Hybrid, Alt, RCT | 500 |
| S12 | null, four configurations | BL SD 1.5/1.5/0.75/2, BR SD 1/1/1.5/0.5, half-life 0/3/3/3 d; assumed 0, 0.2 wk | Hybrid, Alt, XO, RCT | 2000 |
| S12b | refit without corCAR1 | null, baseline configuration | Alt, XO | 2000 |
| S14 | parallel-group check | rho 0, 0.3, 0.7; effect 0, -0.5; legacy, ancova | RCT, RCT-OL | 1000 |

## Deliverables

- Drivers at `analysis/scripts/nof1-design-sensitivity/`. Every
  driver asserts that the parameter values stored in the artifact
  equal its grid and that run plus dropped cells equal the grid size.
  A failed assertion saves the artifact as `*-UNCHECKED.rds`.
- Artifacts `sensitivity-NN-name.rds` with the run metadata
  (`run_meta`: schema, replicates, seed, engine version, revision).
- `cell-summaries.csv` from `16-summarize-all.R`.
- Figures from `plot-all.R`.

## Seeds

Each driver calls `set.seed(42 + id)` with the sweep number as `id`.
Cells run in parallel with `furrr` and `seed = TRUE`. There is no
separate seed file.

## Deviations log

- 2026-10-28. Protocol rewritten after the referee report of
  2026-10-28 (`analysis/report/05-nof1-design-sensitivity/
  referee-report-2026-10-28.md`). Changes: published Hybrid schedule;
  common estimand; ANCOVA for the parallel group; analysis half-life
  crossed with the DGP half-life in S2, S3, S11, S12; feasibility
  check with dropped cells; S4 now varies $\tau$ (added to the DGP)
  and, separately, the baseline SD; S5 scales all component SDs; S6
  fixes the trial length; S9 passes the whole grid; S12b and S14
  added; baseline effect changed from a response asymptote of 2 to
  $\Delta_{\mathrm{on}} = -0.5$.
- The earlier references to `01-base-seed.txt`,
  `02-deviations-log.md`, `02-grid-summary.rds` and
  `docs/sensitivity-sweep-plan-2026-04-17.md` are withdrawn. None of
  these files exists.
