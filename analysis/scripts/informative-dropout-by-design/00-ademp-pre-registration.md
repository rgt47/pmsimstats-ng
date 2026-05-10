# ADEMP pre-registration: 09-informative-dropout-by-design

*2026-05-07 17:00 PDT*

This document pre-registers the simulation programme that
underlies `analysis/report/09-informative-dropout-by-design/`.
ADEMP structure follows Morris, White and Crowther (2019).

## Common infrastructure

**Trial designs (four-design factorial).**

1. *Open-label* (OL): 20 weeks of open-label active drug.
2. *Open-label plus blinded discontinuation* (OL+BDC): 16 weeks
   of OL active + 4 weeks of blinded discontinuation.
3. *Traditional crossover* (CO): 10 weeks period 1 + 10 weeks
   period 2, with within-participant randomisation.
4. *Hybrid N-of-1*: 8 weeks OL + 4 weeks BDC + 4 weeks
   crossover-1 + 4 weeks crossover-2 (the Hendrickson 2020
   primary design).

All four designs at 20 weeks total to match Hendrickson 2020
Figure 2.

**Sample sizes.** $N \in \{35, 70, 100\}$ per design.

**Drug context.** Prazosin for PTSD nightmares, calibrated to
the Hendrickson 2020 reference parameters.

**Carryover.** Exponential decay with $t_{1/2} = 3$ days at
the baseline cell.

**Replicates.** $n_{\text{reps}} = 1000$ per cell for power and
bias estimands; $n_{\text{reps}} = 5000$ for type I error
cells.

**Random-number control.** Per-replicate seed via cell-descriptor
hash; base seed in `01-base-seed.txt`.

**Reporting standard.** Morris-White-Crowther ADEMP with **MCSE
columns alongside every reported performance measure**.

## Study 1. Power × design × dropout pattern

**A.** Reproduce and extend Hendrickson 2020 Figure 4A: power
as a function of (design family) × (dropout pattern) × (sample
size) × (biomarker effect size).

**D.** Three-component DGP from the pmsimstats simulation engine
with the @hendrickson2020 hazard-based dropout layer. Five
dropout patterns:

1. *No dropout* (reference): $\beta_0 = 0$, $\beta_1 = 0$.
2. *Balanced*: $\beta_0 = 0.05$, $\beta_1 = 0.5$ (Hendrickson
   default).
3. *More-of-flat*: $\beta_0 = 0.05$, $\beta_1 = 0.2$ (treatment
   non-responders drop more).
4. *More-of-biased*: $\beta_0 = 0.05$, $\beta_1 = 0.8$
   (treatment-responders drop more on discontinuation).
5. *High-dropout*: $\beta_0 = 0.15$, $\beta_1 = 0.5$ (overall
   higher rate, balanced direction).

The hazard at each timepoint is
$\Pr(\text{drop at } t) = \beta_0 + \beta_1 \cdot
\Delta_{Sx}(t)^2$ where $\Delta_{Sx}$ is the change in
symptom score since baseline scaled to $[0, 100]$.

**E.** Primary estimand: power $\pi = \Pr(p < 0.05)$ for
$H_0: \beta_{bm:D} = 0$ in the linear-mixed analysis. Secondary
estimand: bias of $\hat{\beta}_{bm:D}$ relative to the
known-true effect at each cell.

**M.** Linear-mixed analysis `Sx ~ bm + t + Dbc + bm:Dbc +
(1 | ptID)` with `corCAR1(form = ~t|ptID)`, fitted to the
post-dropout retained data. Mixed-effects ML missing-data
inference applies under the assumption of MAR; we explicitly
note that under our DGP dropout depends on observed symptom
trajectory, which gives MNAR semantics relative to the
unobserved future trajectory but MAR semantics relative to the
observed past.

**P.** Power, type I error under $\beta_{bm:D} = 0$, bias,
empirical Monte Carlo SE of $\hat{\beta}_{bm:D}$ across
replicates, model SE, nominal-95% CI coverage, convergence
rate, mean dropout fraction. **MCSE attached to every reported
estimate**, computed as
$\sqrt{\hat{\pi}(1-\hat{\pi})/n_{\text{reps}}}$ for proportions
and the standard Monte Carlo SE of the within-replicate mean
for continuous measures.

**Pre-registered prediction.** Power loss from dropout exceeds
50% in the most-vulnerable design-by-dropout cells (high
dropout in OL+BDC and N-of-1 designs) and remains negligible in
the OL design across all dropout patterns. Design rankings on
power reverse depending on the dropout pattern: OL+BDC and
N-of-1 dominate under no dropout but fall below CO under
symptom-improvement-driven dropout.

## Study 2. Bias quantification (two-mode)

**A.** Quantify the bias of the biomarker-treatment interaction
estimate under each dropout pattern, in two reference modes:
(a) bias relative to the known true value at the cell, (b)
bias relative to the gold-standard estimate from the full
uncensored dataset.

**D.** As Study 1, but with bias as the primary estimand and the
$c_{bm} \in \{0.3, 0.6\}$ cells the focus.

**E.** $\Delta\beta = \hat{\beta}_{bm:D}^{\text{dropped}} -
\beta_{bm:D}^{\text{true}}$ (mode a) and
$\Delta\beta = \hat{\beta}_{bm:D}^{\text{dropped}} -
\hat{\beta}_{bm:D}^{\text{gold}}$ (mode b), where
$\hat{\beta}_{bm:D}^{\text{gold}}$ is the analysis-model
estimate from the full uncensored data.

**M.** As Study 1.

**P.** Mean and Monte Carlo SE of $\Delta\beta$ across
replicates; the proportion of replicates where $\Delta\beta$
crosses zero (effect-direction reversal under dropout).

## Study 3. Randomization-path 'happy accident' decomposition

**A.** Quantify the @hendrickson2020 'happy accident' selection
effect: in OL+BDC and hybrid N-of-1 designs, the multiple
randomization paths interact with biased dropout to
preferentially preserve the most-informative crossover blocks.

**D.** Three-component DGP at the $c_{bm} = 0.3$ effect size,
with the OL+BDC and hybrid N-of-1 designs only, under the
balanced and more-of-biased dropout patterns.

**E.** Power conditional on randomization path. Within OL+BDC,
the path indicator records which 1-week sub-block of the
4-week BDC phase the participant was discontinued; within
hybrid N-of-1, the indicator records the on/off pattern across
the two crossover blocks.

**M.** As Study 1, with path-conditional analysis.

**P.** Path-conditional power; the difference between the
weighted average power (weighted by the empirical path
distribution under dropout) and the unweighted average power
(weighted by the design's nominal path distribution). The
difference quantifies the 'happy accident' selection effect.

## Study 4. Sensitivity to dropout-pattern misspecification

**A.** Quantify the inferential cost of analysing the data
under an incorrect missing-data assumption (MCAR when truth is
MNAR; MAR with the wrong covariate structure).

**D.** As Studies 1-2, with the analysis-side missing-data
treatment varied:

1. Naive complete-case (drop participants on first missing).
2. MAR-conditional MMRM (the Study 1 analysis).
3. Multiple-imputation under the correct MNAR mechanism.

**E.** Bias and power as in Studies 1 and 2.

**M.** Three analysis variants per the D specification.

**P.** Power, bias, MCSE, coverage. The differences between
the three analysis variants quantify the cost of mis-handling
the missing-data mechanism.

## Total cell counts

- Study 1: $4 \text{ designs} \times 5 \text{ patterns} \times
  3 N \times 3 c_{bm} = 180$ cells.
- Study 2: subset of Study 1 at $c_{bm} \in \{0.3, 0.6\}$ =
  120 cells.
- Study 3: $2 \text{ designs} \times 2 \text{ patterns} \times
  3 N$ = 12 cells, with path-conditional decomposition adding
  per-cell complexity.
- Study 4: $4 \text{ designs} \times 5 \text{ patterns} \times 3
  \text{ analyses} \times N = 35$ = 60 cells.

## Deliverables

- Driver scripts at `analysis/scripts/informative-dropout-by-
  design/`.
- ADEMP results tables per study **with MCSE columns alongside
  every estimate** (per Morris recommendation, addressing the
  manuscript-stage gap flagged in the project ADEMP audit).
- Figures: power-by-(design, dropout-pattern) heatmap (Study 1
  reproduces Hendrickson Figure 4A and extends); bias plots
  (Study 2); path-conditional power decomposition (Study 3);
  analysis-method sensitivity (Study 4).

## Deviations log

`02-deviations-log.md` records departures from the cells,
estimands, methods, or performance measures pre-registered
above.

---

*Source:* `~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/informative-dropout-by-design/00-ademp-pre-registration.md`
