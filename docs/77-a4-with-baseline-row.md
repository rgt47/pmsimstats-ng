---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
  - \usepackage{float}
  - \floatplacement{figure}{H}
---

# The Recommended Analysis A4 With and Without the Baseline Row {.unlisted .unnumbered}
*2026-11-20 15:05 PST*

**Author.** pmsimstats team

**Purpose.** The published analysis of Hendrickson et al. adds one
outcome row per patient at baseline ($t = 0$, off drug). None of paper
02's ladder analyses has that row. The refit of the published analysis
(script 54) showed that, for a binary drug indicator, the row is worth
13 to 17 points of power. Here we ask the same question of A4, the
analysis paper 02 recommends: what happens to its size, power and
estimate when the baseline row is added and nothing else changes?

**Status.** Exploratory. The comparison was designed after the script
54 result was seen. Paper 02's text and recommendation are unchanged;
this note informs that decision.

## Summary

1. Size does not change. Pooled over the six null cells, A4 rejects
   at 0.051 and A4 with the baseline row at 0.052 (paired difference
   +0.001, SE 0.004). Both are inside the Monte Carlo band in every
   cell.
2. Power rises by 17 points at the six reference alternatives (0.372
   to 0.544; paired difference +0.171, SE 0.008), and by 12 points on
   average over all 24 alternative cells. The gain is positive in
   every cell, ranges from 4 to 29 points, and is largest in the
   OL+BDC design.
3. The estimate improves as well. With the baseline row, A4's bias is
   smaller, its empirical SD is 18 to 36 percent smaller, and the
   ratio of mean SE to empirical SD moves closer to 1 (0.90 to 0.96
   under the null).
4. With the row, A4 nearly matches the oracle A3 with the row (0.544
   against 0.548) and leads the binary analysis A2 with the row by 5.9
   points (SE 0.6). Without the row, A4 led A2 by about 2 points.
5. The gain rests on a property of the data-generating process. Every
   visit outcome is $Y_t = BL - (TV_t + PB_t + BR_t)$, so the baseline
   row records each patient's intercept without any visit-level
   noise. A real baseline measurement is noisier. The size of the
   gain in practice is therefore not established here.

## 1. Methods

**Trials.** The 30 cells of job 27 (`tier1-g-vs-f-e0`) under parameter
set 1, Hendrickson et al.'s values: three designs (CO, Hybrid,
OL+BDC), $N \in \{35, 70\}$, $c_{bm} \in \{0, 0.30, 0.45\}$ at
$t_{1/2} = 1$ week, and $t_{1/2} \in \{0, 0.5\}$ weeks at
$c_{bm} = 0.45$; exponential carryover, separable AR(1) with
$\rho = 0.7$, complete data, predictive-only biomarker ($c_{bl} = 0$),
500 replicates per cell. Each trial was regenerated with job 27's seed
($24102600 + r$), exactly as script 54 does.

**Analyses.** Three ladder analyses were fitted to every trial twice,
once as stored and once with the baseline row:

| Analysis | Stored code | Fixed effects | Inference |
|---|---|---|---|
| A2 | E1cr2 | $Y \sim bm + t + Db + bm{:}Db$ | CR2, Satterthwaite df |
| A3 (oracle) | E2cr2 | $Y \sim bm + t + D_{bc} + bm{:}D_{bc}$, true $t_{1/2}$ | CR2, Satterthwaite df |
| A4 | E7null_cr2 | as A3; $t_{1/2}^a$ chosen by ML-AIC over $\{0.25, 0.5, 1, 2\}$ weeks on the model without $bm{:}D_{bc}$, then REML refit | CR2, Satterthwaite df |

All have a patient random intercept and `corCAR1(form = ~ t | ptID)`
residuals. The baseline row sets $Y = BL$, $t = 0$, and
$Db = D_{bc} = L = t_{sd} = D_e = 0$. $D_{bc}$ rebuilt at any
candidate half-life is therefore 0 on that row, and the row enters the
half-life selection as well as the final fit. The row is the only
change. Every patient is kept: the published analysis's rule that
drops never-treated patients is not used. The data are sorted by
patient and time before fitting, because the CR2 step clusters by row
position; a check on one trial per design confirmed that the CR2
$p$-value does not depend on row order.

**Guard.** The plain refits of A2, A3 and A4 must reproduce the stored
job 27 estimates within $10^{-4}$, and A4 must select the same
half-life. All 15,000 replicates passed (largest difference
$3.7 \times 10^{-5}$; no half-life differed).

**Comparator.** The published analysis with `corCAR1` added (E0pnc in
script 54: binary drug indicator, baseline row, no expectancy term,
model-based SE) was fitted to the same trials and is shown for
reference.

**Monte Carlo error.** Replicate $r$ of every cell shares seed $r$, so
the cells are not independent. Pooled paired differences are averaged
first within each replicate across cells, and the SE is that of the
500 replicate means. Paired comparisons use replicates in which both
arms converged. Non-convergence was rare: 25 of 15,000 fits for A4 and
1 for A4 with the baseline row.

## 2. Size

The six null cells all have $t_{1/2} = 1$ week. The Monte Carlo band
for 500 replicates at $\alpha = 0.05$ is about 0.031 to 0.069.

| Design | $N$ | A2 | A2 + BL | A3 | A3 + BL | A4 | A4 + BL | E0pnc |
|:-------|----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| CO     | 35 | 0.050 | 0.054 | 0.052 | 0.058 | 0.050 | 0.054 | 0.060 |
| CO     | 70 | 0.052 | 0.040 | 0.046 | 0.036 | 0.040 | 0.034 | 0.056 |
| Hybrid | 35 | 0.060 | 0.062 | 0.050 | 0.052 | 0.060 | 0.054 | 0.048 |
| Hybrid | 70 | 0.052 | 0.060 | 0.062 | 0.058 | 0.062 | 0.064 | 0.048 |
| OL+BDC | 35 | 0.058 | 0.062 | 0.044 | 0.052 | 0.044 | 0.056 | 0.036 |
| OL+BDC | 70 | 0.050 | 0.046 | 0.052 | 0.050 | 0.050 | 0.052 | 0.020 |
| Pooled |    | 0.054 | 0.054 | 0.051 | 0.051 | 0.051 | 0.052 | 0.045 |

: Null rejection rate at $\alpha = 0.05$. BL: with the baseline row.

The paired differences in pooled size are +0.001 for A4 (SE 0.004),
+0.000 for A2 and +0.000 for A3. The baseline row does not inflate the
size of any CR2 analysis. Only E0pnc leaves the band, in one cell
(OL+BDC, $N = 70$, 0.020).

The ratio of mean SE to empirical SD under the null is 0.90 for A4 and
0.96 for A4 with the row. The CR2 standard error is slightly small for
plain A4, and the row brings it close to the truth.

## 3. Power

| Design | $N$ | A2 | A2 + BL | A3 | A3 + BL | A4 | A4 + BL | E0pnc |
|:-------|----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| CO     | 35 | 0.269 | 0.344 | 0.261 | 0.350 | 0.261 | 0.360 | 0.412 |
| CO     | 70 | 0.600 | 0.712 | 0.584 | 0.718 | 0.580 | 0.718 | 0.744 |
| Hybrid | 35 | 0.236 | 0.338 | 0.294 | 0.422 | 0.277 | 0.424 | 0.339 |
| Hybrid | 70 | 0.519 | 0.658 | 0.588 | 0.746 | 0.553 | 0.738 | 0.652 |
| OL+BDC | 35 | 0.151 | 0.302 | 0.182 | 0.368 | 0.198 | 0.366 | 0.246 |
| OL+BDC | 70 | 0.329 | 0.558 | 0.379 | 0.682 | 0.365 | 0.658 | 0.495 |
| Pooled |    | 0.351 | 0.485 | 0.381 | 0.548 | 0.372 | 0.544 | 0.481 |

: Power at the six reference alternatives ($c_{bm} = 0.45$,
$t_{1/2} = 1$ week).

The pooled paired gains from the baseline row are +0.171 for A4 (SE
0.008), +0.166 for A3 and +0.135 for A2. A2 with the row (0.485)
reproduces E0pnc (0.481). Both have a binary drug indicator and the
row, and they differ in the SE (CR2 against model-based) and the
expectancy and exclusion details. That agreement is the check that the
row is coded as the published analysis codes it.

A4 with the row exceeds E0pnc by 6.2 points (SE 0.7), at similar
size (0.052 against 0.045).
In the CO design, E0pnc is still ahead by 3 to 5 points.
In the Hybrid and OL+BDC designs, A4 with the row is ahead by 8 to 16 points.

Over all 24 alternative cells, the gain for A4 ranges from 0.040 to
0.293, with a mean of 0.124 (paired SE 0.005). It is positive in every
cell. By design, the mean gain is 0.100 (CO), 0.104 (Hybrid) and
0.168 (OL+BDC). By true half-life it is 0.106 at 0, 0.128 at 0.5 and
0.131 at 1 week. The largest gain, 29 points, is in OL+BDC with
$N = 70$ at $t_{1/2} = 1$ week (0.365 to 0.658).

## 4. Estimation

At $c_{bm} = 0.45$ the target is $\theta = -c_{bm}\sigma_{BR} = -3.6$.

| Analysis | $t_{1/2}$ | Mean estimate | Bias | Empirical SD | Mean SE | SE/SD |
|:---------|----:|------:|------:|-----:|-----:|-----:|
| A3       | 0   | -3.599 |  0.001 | 1.715 | 1.615 | 0.94 |
| A3 + BL  | 0   | -3.597 |  0.003 | 1.471 | 1.396 | 0.95 |
| A4       | 0   | -3.828 | -0.228 | 1.867 | 1.741 | 0.93 |
| A4 + BL  | 0   | -3.708 | -0.108 | 1.523 | 1.445 | 0.95 |
| A3       | 0.5 | -3.603 | -0.003 | 1.928 | 1.826 | 0.95 |
| A3 + BL  | 0.5 | -3.601 | -0.001 | 1.587 | 1.517 | 0.96 |
| A4       | 0.5 | -3.802 | -0.202 | 2.281 | 1.984 | 0.87 |
| A4 + BL  | 0.5 | -3.604 | -0.004 | 1.635 | 1.531 | 0.94 |
| A3       | 1   | -3.611 | -0.011 | 2.294 | 2.172 | 0.95 |
| A3 + BL  | 1   | -3.607 | -0.007 | 1.714 | 1.648 | 0.96 |
| A4       | 1   | -3.804 | -0.204 | 2.746 | 2.350 | 0.86 |
| A4 + BL  | 1   | -3.558 |  0.042 | 1.756 | 1.636 | 0.93 |

: Estimation at $c_{bm} = 0.45$, pooled over designs and $N$.

The row leaves the target unchanged: the oracle A3 is unbiased with
and without it. The row makes the estimate more precise. For A4 it
also roughly halves or removes the bias left by half-life selection,
and it brings the SE ratio from 0.86 to 0.93 at $t_{1/2} = 1$ week.
The improvement is largest in OL+BDC. There, plain A4 has a mean
estimate of $-4.31$ and an SE ratio of 0.86 at the reference
alternatives. With the row, these are $-3.31$ and 0.91, and the
empirical SD falls from 3.63 to 1.76.

## 5. Half-life selection

| True $t_{1/2}$ | Analysis | 0.25 | 0.5 | 1 | 2 |
|----:|:---------|-----:|-----:|-----:|-----:|
| 0   | A4       | 0.70 | 0.18 | 0.07 | 0.05 |
| 0   | A4 + BL  | 0.79 | 0.10 | 0.06 | 0.05 |
| 0.5 | A4       | 0.43 | 0.24 | 0.22 | 0.11 |
| 0.5 | A4 + BL  | 0.41 | 0.33 | 0.16 | 0.10 |
| 1   | A4       | 0.22 | 0.14 | 0.27 | 0.37 |
| 1   | A4 + BL  | 0.18 | 0.20 | 0.30 | 0.31 |

: Fraction of replicates selecting each analysis half-life (weeks),
all designs, $N$ and $c_{bm}$.

The two versions select the same half-life in 62 percent of
replicates. With the row, the selection is a little more concentrated
near the truth at $t_{1/2} = 0$ and 0.5 weeks, and about the same at 1
week. The change in selection is not the main source of the gain: the
oracle A3, which does no selection, gains nearly as much (16.6 against
17.1 points).

## 6. Why the row helps, and the main caveat

In the data-generating process, each visit outcome is
$Y_{it} = BL_i - (TV_{it} + PB_{it} + BR_{it})$. Under parameter set 1
the components have SDs of 10 ($TV$), $10\eta_t$ ($PB$) and 8 ($BR$),
and $BL_i$ has SD 18.5. The baseline row records $BL_i$ itself. It is
an off-drug, zero-expectancy observation of the patient intercept
with no visit-level noise at all. Without it, the intercept must be
estimated from visits that each carry the full component noise.

A precise intercept sharpens the on-drug against off-drug contrast
within each patient. The gain is largest where off-drug visits are
fewest. Of the eight visits, OL+BDC has two or three off drug, Hybrid
three or four and CO four, and the mean gains follow that order
(0.168, 0.104, 0.100).

A real baseline assessment is a measurement with its own error,
typically of the same order as the error at any other visit. The
gains in Sections 3 and 4 are therefore an upper bound for a trial
whose baseline is as noisy as its other visits. This is the published
DGP's construction, so the published analysis has had the same
advantage in every simulation of Hendrickson et al. and of paper 02.

## 7. Limitations

- **Noise-free baseline.** As above. A sensitivity run that adds
  visit-level noise to the baseline row would measure how much of the
  gain survives. It has not been run.
- **Prognostic biomarker.** All cells have $c_{bl} = 0$. With a
  biomarker correlated with baseline severity, the row adds a direct
  $bm$-to-outcome association at $t = 0$, and its effect on size is
  untested.
- **Null cells.** There are six null cells, all at $t_{1/2} = 1$ week,
  with 500 replicates each. Size at $t_{1/2} = 0$ and 0.5 weeks was not
  checked.
- **Scope.** Parameter set 1 only, complete data, exponential decay,
  $\rho = 0.7$. Families 2 to 5 of paper 02 (misspecified decay, other
  $\rho$, dropout, patient slopes) were not rerun.
- **Coverage.** The Satterthwaite df are not stored, so interval
  coverage was not computed. The SE ratios in Section 4 are the
  available proxy.

## 8. Implications for paper 02

These results bear on two open decisions. First, the case for A4 over
A2 is stronger with the row (5.9 points) than without it (about 2
points). Second, the baseline row itself is the largest single source
of power found in the ladder, larger than any rung. If the row is
adopted, A2 to A4 would be refitted with it across all families, and
the result would need the noisy-baseline sensitivity run before it is
stated as a recommendation.

## Reproducibility

- `analysis/scripts/carryover-sensitivity/55-run-a4-baseline.R`
  (`--check` verifies the row construction and row-order invariance;
  `--time` gives the cost). The run of 2026-11-20 used six local cores
  for 61 minutes. Output:
  `analysis/data/carryover-sensitivity/published/a4-baseline/`
  (`all.rds`, `guard.rds`, `summary.rds`, `meta.rds`); the per-trial
  checkpoints under `trials/` are not tracked.
- `analysis/scripts/carryover-sensitivity/56-a4-baseline-summary.R`
  prints every number in this note from those files and from
  `published/a0-published/` (E0pnc).
