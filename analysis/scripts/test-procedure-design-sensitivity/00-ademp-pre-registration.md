# ADEMP pre-registration: 08-test-procedure-design-sensitivity

*2026-05-07 16:55 PDT*

This document pre-registers the simulation programme that
underlies `analysis/report/08-test-procedure-design-sensitivity/`.
ADEMP structure follows Morris, White and Crowther (2019).

## Common infrastructure

**Trial designs.** Hendrickson hybrid OL+BDC as the reference,
plus three trial-design variants on the cycle-by-period grid
(see Study 2 below).

**Sample sizes.** $N = 35$ per design at the baseline, varied
across $\{35, 70, 100\}$ in the $N$-axis sweep.

**Drug context.** Prazosin for PTSD nightmares; baseline
parameters at the Hendrickson 2020 reference values.

**Carryover.** Exponential decay with $t_{1/2} = 3$ days at
the baseline cell.

**Replicates.** $n_{\text{reps}} = 1000$ per cell for power and
bias estimands; $n_{\text{reps}} = 5000$ for type I error
cells.

**Random-number control.** Per-replicate seed via cell-descriptor
hash; base seed in `01-base-seed.txt`.

**Reporting standard.** Morris-White-Crowther ADEMP with MCSE
columns for every reported performance measure.

**Simulation substrate.** `implementations/simple/simulation.R`
extended with the cycle-by-period axes (rather than the full
pmsimstats package stack), per the substrate-choice argument
in §1 of the manuscript.

## Study 1. Test-procedure comparison at fixed design

**A.** Quantify the relative power, type I error, and bias of
three test procedures (strict classical RM-ANOVA, linear mixed-
effects with corCAR1 residual correlation, GEE with bias-
corrected sandwich variance) for the biomarker-treatment
interaction, holding the trial-design parameters at the
Hendrickson hybrid reference.

**D.** Three-component additive DGP at the prazosin baseline,
biomarker effect $c_{bm} \in \{0, 0.3, 0.6\}$, with the design
fixed at the Hendrickson hybrid (8-week OL + 4-week BDC +
4-week + 4-week crossover blocks).

**E.** Primary estimand: $\hat{\beta}_{bm:D}$ from the analysis-
model fixed-effects table for the linear-mixed and GEE
analyses; the corresponding $F$-statistic from the strict
RM-ANOVA. The $p$-value associated with each is the
inferential output that drives power and type I error.

**M.** Three analysis variants compared:

1. *Strict classical RM-ANOVA* per Greenhouse-Geisser-corrected
   F-statistic on the biomarker-stratified treatment-by-time
   interaction in a balanced split-plot design.
2. *Linear mixed-effects with corCAR1*: `Sx ~ bm + t + Dbc +
   bm:Dbc + (1 | ptID)` with `corCAR1(form = ~t|ptID)` (the
   pmsimstats default analysis).
3. *GEE with bias-corrected sandwich variance* per
   @wang2016geesmv: `geeglm(Sx ~ bm + t + Dbc + bm:Dbc, id =
   ptID, corstr = 'ar1')` with the Mancl-DeRouen correction.

**P.** Power, type I error, bias, empirical SE, model SE,
nominal-95% CI coverage, convergence rate. MCSE attached to
every estimate. Coverage is computed for the linear-mixed and
GEE analyses (where intervals are produced); for RM-ANOVA only
the rejection-rate quantities are reported.

**Pre-registered prediction.** Linear-mixed dominates GEE on
power at trial-relevant N due to the small-sample bias of the
sandwich estimator, even with correction; both dominate strict
RM-ANOVA when within-participant residuals depart from
compound symmetry, which is the empirically common case.

## Study 2. Cycle-by-period design grid at fixed test procedure

**A.** Quantify the dependence of the biomarker-treatment
interaction test on the cycle-by-period design grid, holding
the test procedure at the linear-mixed-with-corCAR1 baseline.

**D.** Three-component DGP at the prazosin baseline, with the
design parameters varied across:

- Number of treatment cycles $k \in \{1, 2, 3, 4, 6\}$.
- On-drug period length $L_{\mathrm{on}} \in \{1, 2, 3, 4\}$
  weeks.
- Off-drug period length $L_{\mathrm{off}} \in \{1, 2, 3, 4\}$
  weeks.
- On/off symmetry indicator (factorialised as $L_{\mathrm{on}} =
  L_{\mathrm{off}}$ vs $\neq$).

**E.** $\hat{\beta}_{bm:D}$ from the linear-mixed analysis;
power and bias as a function of the design grid.

**M.** Linear-mixed-with-corCAR1 analysis only.

**P.** Power, type I error, bias, MCSE, convergence rate, mean
fit time. The fit-time measure characterises the computational
cost of long-cycle designs.

**Pre-registered prediction.** Per-participant power is
monotone in cycle count $k$ with diminishing returns above
$k = 4$. On-drug period length matters more for power than
off-drug period length when carryover is non-trivial.
Asymmetric $L_{\mathrm{on}} > L_{\mathrm{off}}$ splits are
preferable when carryover is non-trivial.

## Study 3. Joint design-by-test heatmap

**A.** Synthesise Studies 1 and 2 into a joint design-by-test
recommendation table.

**D.** Same DGP grid as Study 2, with the three test procedures
from Study 1 each fitted to every cell.

**E.** Power and type I error for each (design, test) combination.

**M.** As Studies 1 and 2.

**P.** A heatmap of power as a function of design parameters
crossed with test procedure, with MCSE annotation on every
cell.

## Total cell counts

- Study 1: $3 \text{ tests} \times 3 N \times 3 c_{bm} = 27$
  cells.
- Study 2: $5 k \times 4 L_{\mathrm{on}} \times 4 L_{\mathrm{off}}
  \times 3 c_{bm} = 240$ cells.
- Study 3: subset of Studies 1-2 with all three tests applied:
  $3 \text{ tests} \times 5 k \times 4 L_{\mathrm{on}} \times 4
  L_{\mathrm{off}} = 240$ cells.

## Deliverables

- Driver scripts at `analysis/scripts/test-procedure-design-
  sensitivity/`.
- ADEMP results tables per study with **MCSE columns alongside
  every estimate**, per Morris recommendation.
- Figures: test-procedure comparison plots faceted by N and
  $c_{bm}$; design-grid heatmaps; joint design-by-test heatmap.

## Deviations log

`02-deviations-log.md` records departures from the cells,
estimands, methods, or performance measures pre-registered
above.

---

*Source:* `~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/test-procedure-design-sensitivity/00-ademp-pre-registration.md`
