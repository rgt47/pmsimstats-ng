---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# Carryover-Representation and Half-Life-Selection Methods for Paper 02: A Catalog of the G- and F-Families {.unlisted .unnumbered}
*2026-10-21 11:05 PDT*

**Author.** pmsimstats team

**Purpose.** Two separate, only partly reconciled lines of work bear
on which analysis specification paper 02 should recommend. The first
is paper 02's own nine specifications, `G1`-`G9`, compared on a
simple model (`Sx ~ bm + t + X + bm:X`, no expectancy term) that
differ only in how carryover enters as `X`. The second is a
project-wide research thread on choosing the analysis-side carryover
half-life, built on the richer "F1" model
(`Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De`, with expectancy terms),
documented in `docs/57-data-chosen-carryover-half-life.md` and
`docs/series/ch23-data-chosen-half-life.md`. This note catalogs both
families in one place: what each method is, what has been determined
about it so far (evidence-labeled), and its approximate computational
cost where available. It does not propose a recommendation; it is a
reference document for deciding what paper 02's own remediation
(`docs/60-trials-submission-gap-analysis.md`, Part 2) should test
next. It is written for a reader with statistical training (mixed
models, cluster-robust inference, model selection by AIC) but no
prior exposure to this project's own method names; the "Notation and
evidence labels" section below defines the recurring statistical
devices (CR2, AIC selection, bootstrap calibration) in that light
before the catalog itself begins.

**What "best" would mean.** No method cataloged here has been checked
against all of the following, and that gap is itself the main finding
of this note. Reading the "Determined" entries against this list is
the intended way to use the catalog — a method's entry states which of
these it has been shown to satisfy, which it has been shown to fail,
and which remain untested, rather than offering a single verdict.

1. **Type I error (size).** Does the test hold its nominal level under
   the null, across the cells actually checked — not just at the one
   reference cell most specifications have been tested at?
2. **Power**, both in absolute terms and relative to the other
   candidates, at the reference cell and, ideally, across the grid.
3. **Robustness to carryover half-life misspecification** — does
   performance hold when the analyst's assumed half-life (or the
   method's own data-driven estimate of it) diverges from the true
   one?
4. **Robustness to the decay-shape assumption** (exponential vs.
   Weibull, light- vs. heavy-tailed) — does the method depend on
   getting the *shape* of carryover decay right, not just its rate?
5. **Robustness to trial design** — does a specification validated at
   the Hybrid reference cell still perform well under CO or OL+BDC, or
   is its apparent advantage a Hybrid-specific artifact (as `G3`/`G4`
   turned out to be)?
6. **Robustness to the DGP architecture** (covariance vs. mean
   moderation) — does the ranking among specifications survive a
   change in how the biomarker-treatment interaction is encoded?
7. **Robustness to a prognostic biomarker** — does the method still
   hold nominal size when the biomarker is also associated with
   baseline severity or the natural course, not only predictive of
   treatment response? Paper 01's Study 3 found this the single
   sharpest discriminator among otherwise-similar analyses, and no
   `G`-family specification has been tested against it at all.
8. **Computational cost**, both per-replicate (for simulation studies)
   and per-trial (for analyzing a single real dataset), since methods
   that scale by a bootstrap multiplier $M$ are impractical for some
   uses and routine for others.
9. **Validation completeness itself** — how many of 1-8 have actually
   been checked for a given method, as distinct from how well it did
   on the ones that have been. `G8` and `F1` are furthest along this
   axis; `G9` is the best performer on power alone but the least
   validated against the rest of this list.

```{=latex}
\clearpage
\tableofcontents
\clearpage
```

## Notation and evidence labels

`Db`: binary on-drug indicator. `Dbc`: continuous exposure-decayed
drug indicator, `1` on drug, decaying with time since discontinuation
off drug, `0` before first exposure. `L`: lagged-treatment nuisance
indicator (Jones and Kenward). `bm`/`bmc`: biomarker, raw or centered.
`De`: design-expectancy covariate ($\eta_t$); absent from the
`G`-family, present in the `F`-family. $t_{1/2}^a$: the analysis-side
carryover half-life used to build `Dbc`. $\mathcal{T}_a$: the
candidate half-life grid the data-driven methods search over. Paper
02 uses $\mathcal{T}_a = \{0.25, 0.5, 1.0, 2.0\}$ weeks; `docs/57`/
`ch23` use $\mathcal{T}_a = \{0, 0.25, 0.5, 1, 2\}$ weeks (one more
point, including 0). All models below share a random intercept per
patient plus `corCAR1` residuals — the continuous-time AR(1) residual
correlation `nlme::lme()` fits, matched to the AR(1) structure
assumed to generate the data (`docs/63`, Section 2).

**CR2** (Bell-McCaffrey bias-reduced cluster-robust standard error,
with Satterthwaite degrees of freedom) addresses a specific small-
sample problem: a mixed model's *model-based* standard error for a
fixed effect is derived assuming the random-effects and residual
covariance structure is correctly specified; when it is not (or when
the number of clusters, here patients, is modest), that standard
error can be anti-conservative, giving a test that rejects more often
than its nominal level even under the null. CR2 instead estimates the
variance of the fixed-effect estimator directly from the
patient-level residuals, with a small-sample bias correction, and
pairs it with a data-dependent (Satterthwaite) degrees-of-freedom
approximation for the reference $t$-distribution — the same logic as
a sandwich/robust variance estimator in a GEE, adapted for mixed
models via `clubSandwich`.

**AIC selection**, in this catalog, always means: fit the candidate
model at every half-life in $\mathcal{T}_a$, compute Akaike's
Information Criterion for each fit, and take the half-life whose fit
minimizes it — the usual AIC-based model-selection rule, applied here
to a single nuisance parameter ($t_{1/2}^a$) rather than to choosing
among different covariate sets.

**Bootstrap calibration** (`F1AB`/`F1A0B`, Section "Method 5" below)
answers a specific objection to naive AIC selection: if the same data
both choose $t_{1/2}^a$ and are then used to test the interaction at
that chosen value, the nominal $p$-value no longer has its usual
interpretation, because the selection step is not accounted for (a
post-selection, or "look-elsewhere," inference problem familiar from
subgroup-selection and stepwise-regression settings). The fix
simulated here is a parametric bootstrap: refit a null model to the
observed data, repeatedly simulate new datasets from it, and on *each*
simulated dataset repeat the *entire* selection-then-test procedure
(not just the test). The resulting distribution of test statistics,
under a model where the null is exactly true and the same selection
step has been applied, gives a calibrated critical value or
$p$-value that does account for the selection.

Each claim below carries one of four evidence labels, the project
convention: *quoted* (text reproduced from a source), *inspected*
(read from source code or a manuscript table), *derived* (follows
algebraically from a quoted or inspected fact) or *verified*
(computed by running code in this or an earlier session). Status
entries not otherwise labeled are *inspected* summaries of the cited
section's reported numbers.

## Part 1: the G-family (paper 02's own specifications)

All nine share the base model `Sx ~ bm + t + X + bm:X`, fit by
`nlme::lme` with a random intercept and continuous-time AR(1)
residuals (`corCAR1`); none include an expectancy term. Source:
`analysis/report/02-carryover-sensitivity/report.Rmd`, Table
`tab-spec-summary` and the sections cited per row below.

| Code | Name | Formula for $X$ |
|---|---|---|
| `G1` | Unadjusted (`E1`) | $X = D_{it}$, binary |
| `G2` | Lag-adjusted (`E3`) | $X = D_{it}$ + nuisance $L_{it} = D_{i,t-1}(1-D_{it})$ |
| `G3` | Exposure-weighted (`E2`) | $X = D_{bc,it}$, half-life assumed/fixed |
| `G4` | AIC-selected $t_{1/2}^a$ (`E7`) | $X = D_{bc,it}$, $t_{1/2}^a$ by AIC, interaction already in the selection model |
| `G5` | Paired-difference (`E9`) | per-patient on/off mean difference, OLS on biomarker |
| `G6`-`G8` | `G1`/`G2`/`G3` + CR2 | same $X$ as base spec |
| `G9` | `G4` + CR2 (`E7cr2`) | same $X$ as `G4` |
| retired | `E4`/`E5`/`E6` | biomarker-interacted carryover terms |

Table: The `G`-family, quick index

### `G1` Unadjusted (`E1`)

**Determined.** Baseline benchmark; no carryover representation.
Statistically indistinguishable from `G2` everywhere tested (largest
gap 0.003-0.004 across every Tier 1 and Tier 2 cell; identical to two
decimal places in the higher-precision 24-cell paired rerun) — the
manuscript's clearest negative result. Trails `G3`-`G5` at the
Hybrid/OL+BDC reference cell (0.760 vs. 0.836-0.850 at $t_{1/2}=1.0$)
but is the best-performing specification under CO (0.834, vs.
`G3`/`G4` collapsing to 0.48-0.55). CR2 (`G6`) lifts it modestly
(0.748 $\to$ 0.776) without changing its ranking relative to `G3`-`G5`.

**Cost.** Cheapest of the nine: one `lme` fit per replicate.

### `G2` Lag-adjusted (`E3`)

**Determined.** Adds the classical Jones-and-Kenward lagged-treatment
term. Does no measurable work anywhere in the grid: indistinguishable
from `G1` at every half-life, design, architecture and Tier 2 block
checked. CR2 (`G7`) gives the same modest, ranking-preserving gain as
`G1`$\to$`G6` (0.748 $\to$ 0.764).

**Cost.** Same order as `G1`: one extra fixed-effect column,
negligible extra compute.

### `G3` Exposure-weighted (`E2`)

**Determined.** The current validated default and the most thoroughly
tested of the nine. Beats `G1`/`G2` at Hybrid/OL+BDC; holds under
half-life mis-specification short of a large over-assumption (Block
S2) and under light-tailed decay-shape mis-specification (narrowing
to non-significance only at the two most extreme heavy-tailed shapes,
Section 3.5). CR2 (`G8`) restores nominal size without disturbing the
ranking (Block S6; validated at $N\in\{35,70\}$ in Hybrid
specifically, not itself stress-tested in OL+BDC). Fails markedly
under CO (0.482-0.546 vs. `G1`/`G2`'s 0.834-0.866) — this is inherited
by `G4` almost entirely (Section 3.9). Reverses marginally under mean
moderation (Block S8: trails `G1`/`G2` by 2-3 points instead of
leading).

**Cost.** One `lme` fit at a single pre-specified half-life — same
order as `G1`/`G2`, no grid search.

### `G4` AIC-selected $t_{1/2}^a$ (`E7`)

Selects $t_{1/2}^a = \arg\min_{\mathcal{T}_a} \text{AIC}$ of the model
**with** `bm:Dbc` already in it.

**Determined.** Confirmed 2026-10-21: this is the liberal "F1A"-style
selection (same pattern as the package's own
`characterize_carryover()`, `R/carryover_analysis.R`). Beats `G3` at
the Hybrid/OL+BDC reference cell (0.850 vs. 0.836 at $t_{1/2}=1.0$).
Fails catastrophically under CO: AIC drifts to the longest grid
candidate in 86-97% of replicates regardless of the true half-life,
but this mis-selection costs almost nothing beyond `G3`'s own CO
collapse with the true half-life in hand (0.474 vs. 0.482) — the
dominant failure is Exposure-weighted's CO weakness, not the AIC step
itself. Type I error (Block S11, one cell, $n=2000$): 0.038, "not
inflated" — but this session's quick replication (500 reps, corrected
DGP) reproduced 0.038 at the same cell, while `docs/57`/`ch23`'s
pooled check of the same naive-selection pattern over 11 null cells
found 0.067-0.074; whether `G4`'s single tested cell represents the
grid is open. Not yet tested with an expectancy term or a prognostic
biomarker (paper 01 Study 3 found exactly that combination breaks
simpler models).

**Cost.** Five `lme` ML fits (one per $\mathcal{T}_a$ candidate) plus
one REML refit per replicate — roughly 5-6$\times$ a single `G3` fit.
Matches `docs/57` Method 2's reported cost for the same fit count
(~0.5 sec/trial, *quoted*, Section 4.2). This session's quick check:
~0.6 sec/replicate running an `E7`-equivalent alongside an `F1A0`-style
comparator together (500 reps, ~303 sec total).

### `G5` Paired-difference (`E9`)

Per-patient on-drug minus off-drug mean, OLS on biomarker; no
carryover term, no repeated-measures structure.

**Determined.** Beats `G1` at the reference cell (0.842) but trails
`G4` narrowly. The most architecture-robust of the five Tier 1
candidates — never reverses to a large disadvantage either way (Block
S8: 0.964 under mean moderation vs. `G1`/`G2`'s 0.972/0.976, a smaller
gap than `G3`/`G4`'s). Underperforms `G1`/`G2` under CO but mildly
(0.750-0.778 vs. 0.834-0.866), reflecting that CO's own within-subject
contrast already captures most of the signal a paired-difference
regression targets. No CR2 variant applies (collapses to one row per
patient before fitting). Tested at only one design/cell pairing
outside the Hybrid comparison (Block S9, CO only).

**Cost.** Cheapest of the data-driven alternatives: a single OLS fit
after a deterministic per-patient aggregation — no `nlme::lme` call at
all, materially faster than `G3`/`G4`.

### `G6`-`G8` (`G1`/`G2`/`G3` + CR2)

**Determined.** CR2 lifts every base specification 2-5 points without
disturbing the ranking among `G1`-`G5` (Block S6). `G8` (`G3`+CR2) is
the validated extension of the current default, checked against the
same robustness battery as `G3` itself.

**Cost.** One extra `clubSandwich`/Satterthwaite computation on top of
the base fit; small fixed overhead.

### `G9` (`E7cr2`, `G4`+CR2)

**Determined.** The single best power result anywhere in the
manuscript (0.870 at the reference cell, a 3.8-point CR2 gain over
`G4`'s 0.832 — the largest CR2 correction for any specification),
overtaking `G8` (0.863). Also the least validated: rests on one cell,
with none of `G4`'s own robustness gaps (CO failure, mean-moderation
reversal, untested prognostic-biomarker behavior) closed, and its
Type I error claim (0.052 at that one cell) carries the same
single-cell-generalization caveat as `G4`'s.

**Cost.** Same grid-search cost as `G4`, plus CR2's small fixed
overhead.

### Retired `E4`/`E5`/`E6` (biomarker-interacted carryover terms)

A lagged term, and a continuous elapsed-time term, each crossed with
`bm`.

**Determined.** Underperformed `G1` at every cell examined across
three blocks, two designs, two architectures and three half-lives.
Traced to a structural cause: the DGP has a single true decay
parameter, so crossing a second carryover term with the biomarker
adds a free parameter with no counterpart in the truth being
estimated. Replaced by `G4` and `G5`.

**Cost.** Not reported separately; superseded before production-scale
timing was recorded.

**Batch-level run times** (quoted from `report.Rmd`'s Reproducibility
section; 8 cores unless noted), for context on the relative cost of
running the full nine-specification grid versus a subset:

| Run | Cells | Specs | $n_{sim}$ | Time |
|---|---|---|---|---|
| Tier 1 grid, `G1`-`G9` extension | 30 | all 9 | 500 | ~3300 sec (~55 min) |
| Decay-shape sensitivity, `G1`-`G3` only | 30 | 3 | 500 | ~770 sec |
| Decay-shape sensitivity, `G1`-`G9` extension | 30 | all 9 | 500 | ~2500 sec |
| S1-S4 `G1`-`G9` extension | 25 | all 9 | 500 | ~13800 sec (well above linear scaling from the 100-rep run; not further diagnosed, *quoted*) |
| S6 `G1`-`G9` extension | 10 | all 9 | 500 | ~1210 sec |
| S7 | 1 | 3 | 500 | ~40 sec |
| S8 (architecture sensitivity) | 2 architectures $\times$ 3 half-lives | 5 | 500 each | ~5.5 min |
| S9 (CO) | 2 | 5 (`G1`-`G5`) | 500 | ~8 min |
| S10 (CR2 comparison) | 1 | 8 | 500 | ~2 min |
| S11 (Type I error) | 1 | 6 (`G1`-`G5`,`G9`) | 2000 | ~7 min |
| S12 (mean-moderation counterpart, Tier 1) | 30 | all 9 | 500 | ~55 min on 8 cores |

Table: Batch run times by block, as reported in `report.Rmd`'s
Reproducibility section (*quoted*)

## Part 2: the F-family (half-life-selection methods, `docs/57`/`ch23`)

All built on the richer model
`Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De`; the methods differ only
in how $t_{1/2}^a$ (feeding `Dbc`) is chosen. Source:
`docs/57-data-chosen-carryover-half-life.md`, section numbers as
cited.

| Method | What it does |
|---|---|
| F1 | $t_{1/2}^a$ fixed in advance (0.2 week) |
| F1A | AIC over $\mathcal{T}_a$, interaction already in the selection model |
| F1A0 (+ F1G) | AIC over $\mathcal{T}_a$ on the null model; F1G is a guarded fixed variant |
| Method 3 (DS/DA/ES/EA) | sample-splitting and cross-fitting |
| Method 4 | profile-likelihood interval for $t_{1/2}^a$, not a test |
| Method 6 | joint/continuous estimation of $t_{1/2}^a$ |
| F1AB / F1A0B | bootstrap calibration of F1A / F1A0 |
| Method 7 (MaxT) | $\max_h \lvert T(h) \rvert$ over the grid |
| Method 8 (MA0/MA1/MAB/MA0B) | marginalize over $t_{1/2}^a$ (Bayesian / AIC-weight) |
| Permutation | residual permutation, Freedman-Lane style |

Table: The `F`-family, quick index

### F1 (Method 1, fixed)

$t_{1/2}^a$ fixed in advance (0.2 week, current policy), with a
non-inferential sensitivity grid reported alongside.

**Determined.** Size 0.051 pooled (at most 0.062 in one cell), stored
setting; 0.053 pooled over 12 cells under the refined-process
extension; at the current 0.2-week policy value, size 0.050 with no
power gain over the historical 0.5-week default. "The validated
simple fallback." A fixed value can lose real power when far from the
best-fitting one (size-adjusted power 0.496 at true half-life 0 vs.
0.272 at true half-life 1 week, fixed at 0.5, *quoted* `docs/47`
Section 6).

**Cost.** One fit per trial; five more for the (non-inferential)
sensitivity grid.

### F1A (naive)

Narrative baseline, not separately numbered in `docs/57` — the
failure case motivating Methods 2 and 5. AIC over $\mathcal{T}_a$,
selection model **already includes** the interaction; test as if
$t_{1/2}^a$ were fixed.

**Determined.** Liberal: the same data drive selection and test.
Pooled Type I error 0.074 across 11 non-prognostic null cells (vs.
F1's 0.061 on the same trials); 0.079 specifically where it selected
the longest candidate (2 weeks). This is the pattern `G4`/`G9`
currently implement.

**Cost.** Same fit count as F1A0: five ML fits across the grid plus
one REML refit.

### F1A0 (Method 2) + F1G guard variant

AIC over $\mathcal{T}_a$ on the model **without** the interaction,
refit with it at the selected value; F1G keeps $t_{1/2}^a=0.5$ fixed
unless a null-model profile-likelihood test rejects it.

**Determined.** Removes "about half" of F1A's excess size, not all of
it: 0.067 pooled (illustration) vs. F1A's 0.074; 0.053 pooled over 12
cells in the refined-process extension, matching F1 almost exactly —
"as nominal tests neither gains power over F1." The orthogonality
argument motivating it is not borne out cleanly at $N=70$. F1G
performs almost identically (0.053 pooled), switching away from the
fixed value in only a minority of trials. With bootstrap calibration
($\to$ F1A0B), "the primary data-driven candidate."

**Cost.** Five ML fits and one REML fit with CR2 per trial, about 0.5
seconds (*quoted*).

### Method 3: sample-splitting / cross-fitting (`ch23`'s DS/DA/ES/EA)

Single split: select on half, test on the other half. Cross-fit
average: average the two swapped single-split estimates. Bonferroni
and Cauchy $p$-value combination across the split.

**Determined.** Single split: 0.044 pooled null, power 0.208
(conservative; loses power testing on $N/2$). Cross-fit average:
0.061 pooled, power 0.372 (closer to nominal, more power, but the two
halves' statistics are not independent, so this calibration is
approximate). Cauchy combination: 0.050 and 0.292 (closest to nominal
of the three). Bonferroni is exact and conservative by construction.

**Cost.** About the cost of F1A0 (ten ML fits on half the data plus
two REML fits with CR2); the cost that matters is the power lost
testing on half the sample (noncentrality falls by about $\sqrt{2}$
for a single split).

### Method 4: report the half-life, not a test

A profile-likelihood interval for $t_{1/2}^a$, reported alongside
whichever test is used for inference.

**Determined.** Not itself a test of the interaction. Illustrated on
one trial; coverage not evaluated. "Recommended for reporting in
every case" regardless of which test is used.

**Cost.** 11 seconds per trial (52-point grid, two ML fits per point,
plus a REML refit and CR2 for the profile) — once per trial, not per
bootstrap resample.

### Method 6: continuous/joint estimation

Estimate $t_{1/2}^a$ as a continuous nonlinear parameter rather than
selecting from a candidate grid.

**Determined.** Proposed, not adopted as a working method here. Its
main contribution is a theoretical orthogonality argument supporting
F1A0's validity case. No CR2 analog exists yet for this nonlinear
parameterization.

**Cost.** One `nlme` fit per trial, typically several times an `lme`
fit's cost, with convergence failures reported near parameter
boundaries; or a ~25-point profile plus a sandwich correction.

### F1AB / F1A0B (Method 5, parametric bootstrap calibration)

Calibrate the select-then-test statistic's critical value via a
parametric null-model bootstrap that repeats the entire selection
step in every resample.

**Determined.** This is what actually fixes the size problem.
Illustrated (warp-speed bootstrap): F1AB 0.062 pooled (vs. F1A's own
0.061 on the same trials), power 0.364. Evaluated in the
refined-process extension: F1A0B 0.050 and F1AB 0.050 pooled, mean
size-adjusted power 0.395 and 0.400, maximin shortfall 0.076 and
0.078 — both hold size at this scale. **F1A0B is the project's
current primary data-driven candidate, still described as
"provisional"; it has not been tested in paper 02 or cited by it**
(see `docs/60-trials-submission-gap-analysis.md`, Part 2, items 7 and
12).

**Cost.** Each resample costs one run of the selection procedure
(~1.1 sec for the five-fit grid). For one trial, $M=199$ resamples
took 128 sec on 7 workers; $M=999$ would take ~11 minutes. In a
simulation study this multiplies per-replicate cost by $M$ unless the
"warp-speed" shortcut (one resample per Monte Carlo replicate, pooled
across the cell) is used, bringing it to about twice the cost of the
plain, uncalibrated procedure.

### Method 7: maximum statistic over the grid (MaxT)

Test $\max_{h \in \mathcal{T}_a} \lvert T(h) \rvert$ against its joint
null distribution, rather than picking one $h$; Bonferroni /
multivariate-normal / bootstrap-calibrated variants.

**Determined.** Bonferroni: exact and conservative (0.024 pooled in
the illustration; 0.018, "too conservative," in the extension).
Multivariate-normal reference: approximate, liberal in the extension
(0.072). Warp-speed bootstrap calibration: 0.048 pooled, power 0.384
in the illustration — similar performance to F1A0B. Reports a
$p$-value and the half-life of largest $\lvert T \rvert$, not a
single effect-size estimand, which "weakens its face validity for a
reader who wants an effect size."

**Cost.** Five REML fits with CR2 (about 1 second); bootstrap
calibration adds $M$ times that per trial, the same multiplier
structure as F1AB/F1A0B.

### Method 8: marginalizing over the half-life

Average over $\mathcal{T}_a$ weighted by (posterior or AIC-derived)
model weight, rather than selecting a single $h$ (Bayesian model
averaging; AIC-weight stand-in).

**Determined.** Full Bayesian MCMC version proposed but not run (cost
prohibitive at simulation scale). AIC-weight stand-in illustrated:
0.053 pooled with interaction-included weights, 0.049 with null-model
weights; evaluated in the extension: 0.0455 and 0.040 (both slightly
conservative); with bootstrap calibration (MAB/MA0B): 0.048 and
0.049, mean size-adjusted power 0.401 and 0.399 — essentially
matching F1AB/F1A0B.

**Cost.** Full Bayesian: an assumed (not measured) 1 to 5 minutes per
trial with a custom Stan program, making large-scale calibration
expensive. AIC-weight stand-in: no extra fitting cost beyond the grid
fits already computed for F1A/F1A0.

### Permutation tests

Residual permutation (Freedman-Lane style, within path), as a
nonparametric alternative to the bootstrap.

**Determined.** Naive permutation of the biomarker is invalid when it
is prognostic (destroys an association in the resamples that remains
in the observed data). A valid within-path residual-permutation
version is "proposed with conditions; not recommended over the
parametric bootstrap, which needs fewer exchangeability assumptions
and handles unequal path sizes and dropout directly."

**Cost.** Approximately that of the parametric bootstrap, per
resample, without a clear cost advantage to offset its added
assumptions.

## Part 3: how the two families connect

**The basic difference.** The two families differ first in the base
model, not in the carryover question. Every `G`-family specification
is fit on the simple model `Sx ~ bm + t + X + bm:X`, with no
expectancy term; the families differ among themselves only in how
carryover enters as `X` (binary `G1`, lag-adjusted `G2`,
exposure-weighted `G3`, AIC-selected `G4`, paired-difference `G5`,
each with a CR2 variant `G6`-`G9`). Every `F`-family method is fit on
the richer model
`Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De`, which adds the design-
expectancy covariate `De` and its interaction with the biomarker; the
`F`-family methods differ among themselves only in how the
analysis-side half-life $t_{1/2}^a$ feeding `Dbc` is chosen (fixed in
`F1`; naive AIC in `F1A`; null-model AIC in `F1A0`; either with
bootstrap calibration in `F1AB`/`F1A0B`; and the further grid-based,
profile-likelihood, sample-splitting, max-statistic and
marginalization methods of Part 2). `De` matters because paper 01's
Study 3 found it necessary for robustness to a *prognostic*
biomarker, one associated with baseline severity or natural course
rather than treatment response; no `G`-family specification includes
it or has been tested against that failure mode.

`F1`'s `Dbc` term is the identical construct to `G3`/`G8`'s
Exposure-weighted predictor (both: 1 on drug, exponential decay off
drug after first exposure, 0 before first exposure). `F1A`/`F1A0`/
`F1AB`/`F1A0B` are therefore asking the same question `G4`/`G9` ask
— how to choose $t_{1/2}^a$ — on top of the richer, expectancy-term
model, and (for the `B`-suffixed variants) with the bootstrap
calibration `G4`/`G9` currently lack. Three gaps follow directly from
the two tables above:

1. **`G4`/`G9` implement F1A, not F1A0 or F1A0B.** The naive,
   interaction-included AIC selection is the variant `docs/57`/`ch23`
   found most liberal. Paper 02 has not tested the cheaper partial
   fix (F1A0) or the calibrated fix (F1A0B) on its own grid.
2. **The `G`-family has never been tested with an expectancy term or
   a prognostic biomarker**; the `F`-family's methods are all built
   on a model that already includes one, because paper 01's Study 3
   found it necessary. `G4`/`G9`'s single-cell Type I error check
   (Section S11) cannot speak to this, since its model lacks the
   term that was found to matter.
3. **No single document currently compares the two families
   directly** — `G4`/`G9`'s power against `F1`'s, or either against
   F1A0B, on a shared grid. `docs/60-trials-submission-gap-analysis.md`
   Part 2, items 11-13 track closing this gap.

## References

- `analysis/report/02-carryover-sensitivity/report.Rmd`, Table
  `tab-spec-summary` and Sections 2.4, 2.7, 3.4-3.9.
- `docs/57-data-chosen-carryover-half-life.md`, Section 4 (Methods
  1-9) and Section 8 (the refined-process extension).
- `docs/series/ch23-data-chosen-half-life.md`, construct box and
  candidate grid.
- `docs/60-trials-submission-gap-analysis.md`, Part 2, for the
  remediation items this catalog supports.
- `CLAUDE.md` (pmsimstats-ng project file), "Recommended analysis for
  the interaction (F1)" and "Known open item" entries.
