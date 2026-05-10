# Deviations from 06-component-decomposition pre-registration

*2026-05-07 PDT*

This log records departures from the cells, estimands, methods, or
performance measures specified in
`00-ademp-pre-registration.md`.

## 2026-05-07: prototype quick-sim model formula

**Pre-registered phase-augmented analysis (Section M, item 2):**

```
Sx ~ bm + t + Dbc + phase + Dbc:phase + bm:Dbc + bm:Dbc:phase
     + (1 | ptID)  with corCAR1(form = ~t|ptID)
```

**Issue.** Under the prototype trial design used in
`01-prototype-quick-sim.R` (OL+BDC hybrid; single path with on-drug
indicator `c(1,1,1,1,1, 1,0,0)`), the design matrix for the above
formula has rank 7 in 8 columns. Specifically, `Dbc = 1` whenever
`phase = OL`, so the column `Dbc:phaseBDC` is a linear combination of
the intercept, `Dbc`, and `phaseBDC`. `nlme::lme` fails with
`Singularity in backsolve at level 0, block 1` on every replicate
(observed convergence rate 0/100 across all three PB strengths).

**Deviation.** The prototype drops the redundant `Dbc:phase` term:

```
Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase
     + (1 | ptID)  with corCAR1(form = ~t|ptID)
```

The primary estimand `bm:Dbc` remains identifiable, as does the
phase-moderation term `bm:Dbc:phase` that distinguishes a true
biomarker-treatment interaction from a phase-correlated artefact.

**Scope.** Limited to `01-prototype-quick-sim.R`. Phase 1 of the full
study should resolve the identifiability issue at the design level by
either (a) using a multi-path design (OL+BDC+CO or parallel-group
OL+BDC where some participants stay on drug through BDC and others do
not), which generates `Dbc=0` rows in the OL phase and breaks the
collinearity, or (b) replacing the `phase` factor with the continuous
expectancy variable `e`, which varies within each phase only as
defined by the trial but does not collapse to `1 - Dbc` under the
prototype path.

## 2026-05-07: bias reference

**Pre-registered estimand (Section E, Study A):** bias of
`hat{beta}_{bm:D}` against the DGP-known biomarker-by-BR coefficient.

**Deviation (initial 5-min prototype, superseded).** The first
prototype reported `bias_relative_to_zero_PB`: the within-analysis
difference between each cell's mean estimate and the cell mean at
`m_PB = 0`.

**Resolution (30-min prototype, 2026-05-07).** The 30-min prototype
reports bias against the DGP-implied population coefficient under
Architecture B. With `c.bm = 0.45`, `sigma_br = 5`, `sigma_bm = 1`,
and the symptom convention `Sx = BL - (BR + PB + TV)`, the
population conditional moderation of BR by bm is
`c.bm * sigma_br / sigma_bm = 2.25`, so the population coefficient
on `bm:Dbc` for `Sx` is `TRUE_BETA = -2.25`. Bias is reported as
`mean_beta - TRUE_BETA` with MCSE = `sd_beta / sqrt(n)`.

## 2026-05-07: cell grid expanded to N x pb_strength

**Pre-registered cells (Section D, Study A):** $m_{PB} \in
\{0, 1, 3, 6, 10\}$ x $m_{TV} \in \{-1, 0, 1, 2\}$ at the
prazosin-PTSD reference $m_{BR}$.

**Deviation.** The 30-min prototype uses the 3 x 3 cell grid
$m_{PB} \in \{0, 6.5, 13\}$ x $N \in \{35, 70, 100\}$ at fixed
$m_{TV}$ defaults. This was driven by two considerations:
(i) the previous 5-min run at fixed $N = 70$ hit a rank-deficiency
that masked the framework's bias prediction, motivating an $N$
sweep; (ii) wall-time budget of 1700 s at ~0.16 s/rep limits the
total cell count when targeting 800 reps/cell. The Phase 1 study
reverts to the pre-registered $m_{PB} \times m_{TV}$ grid and
adds the $N$ sweep as a second sub-grid.

## 2026-05-07: phase-augmented formula at N >= 70

**Issue carried forward.** The structural rank deficiency in
`Dbc:phase` (`Dbc = 1` whenever `phase = OL` under the single-path
OL+BDC ondrug pattern `c(1,1,1,1,1, 1,0,0)`) is a property of the
design matrix, not of the sample size. Increasing $N$ does not
restore identifiability of `Dbc:phase` because the collinearity
holds at the row level for every replicate.

**30-min prototype behaviour.** At $N \in \{70, 100\}$ the driver
attempts the full pre-registered phase-augmented formula
including `Dbc:phase`. On every replicate the call to
`nlme::lme` fails with the singularity error and the reduced
formula (without `Dbc:phase`) is used; the dropped term is
recorded in the `formula_dropped` column. At $N = 35$ the driver
goes directly to the reduced formula.

**Implication.** Resolving `Dbc:phase` identifiability requires a
design-level change (multi-path or parallel-group OL+BDC, or
replacing `phase` with continuous `e`); it cannot be resolved by
increasing $N$ within the prototype's single-path design.

---

*Source:* `~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/component-decomposition/02-deviations-log.md`
