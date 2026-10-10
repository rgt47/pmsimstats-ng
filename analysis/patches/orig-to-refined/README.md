# From the published process to the refined process: patches to pmsimstats 58b32a9
*2026-10-08 15:32 PDT*

Author: pmsimstats team

## Purpose

These patches turn Hendrickson's pmsimstats at commit 58b32a9
(github.com/rchendrickson/pmsimstats, "orig") into the refined process of
Paper 01 (`analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`,
Tables `tab:refinements` and `tab:spec`, Appendix B), and its analysis
into the analyses of the stress test, with F1 as the default. They are
written against orig's own code and style, so that the changes can be
read line by line against the published package. The names follow
`analysis/report/NOTATION.md`, "Named data-generating processes": the
published process is orig's `generateData()` and `buildtrialdesign()`;
the refined process is the recommended data-generating process (DGP).

## Files

| File | Content |
|---|---|
| `orig-to-refined.patch` | Core refinements: `R/generateData.R`, `R/generateSimulatedResults.R` |
| `orig-to-refined-options.patch` | Stress-test options, on top of the core patch: `R/generateData.R`, `R/censordata.R` |
| `orig-to-refined-analysis.patch` | The stress-test analyses: `R/lme_analysis.R`, `R/generateSimulatedResults.R`, `DESCRIPTION`; applies to orig alone or on top of the other two |
| `verify.sh` | Fetches 58b32a9, applies the patches, pins the reference driver, runs `verify.R` |
| `verify.R` | The checks, one PASS/FAIL line each |

## How to apply and verify

On a checkout or tarball of 58b32a9, from its root:

```bash
patch -p1 < orig-to-refined.patch
patch -p1 < orig-to-refined-options.patch    # optional
patch -p1 < orig-to-refined-analysis.patch   # optional
```

`git apply` works the same way. All three apply without fuzz, in this
order or with the analysis patch alone on orig.

To verify, from the pmsimstats-ng repository root (needs an
authenticated GitHub CLI; installs nothing):

```bash
bash analysis/patches/orig-to-refined/verify.sh [workdir [--fast]]
```

`verify.sh` fetches the 58b32a9 tarball, builds four trees (`orig`;
`core`; `options` = core + options; `full` = core + options + analysis),
checks each patch with `git apply --check` and `patch` (no fuzz, no
rejects), and extracts the reference driver
`analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R`
and its helper `hendrickson-problems/04-power-simulation.R` as committed
in 3238916. That commit is pinned because the drivers are being changed
to make anchored carryover and randomized allocation their defaults; the
pinned version is the one that produced the stress-test results.
`--fast` uses smaller Monte Carlo samples.

## The reference implementation and where it departs from the refined process

The stress-test driver is the reference for what it shares with the
refined process: the separable AR(1) correlation, graded coupling and
the stop on a matrix that is not positive definite. It does not
implement two elements of `tab:spec` (inspected):

- **Carryover.** Driver line 100 is the published compounding recursion,
  `br[p] + br[p - 1] * w[p]` with the cumulative time since
  discontinuation, not anchored carryover.
- **Allocation.** Driver line 303 splits N over the paths in fixed
  numbers, `N %/% nP` with the remainder to the first paths (18, 18, 17,
  17 for N = 70 in the Hybrid design), exactly as orig's
  `generateSimulatedResults()`. It is not randomized.

The patches implement both as `tab:spec` specifies (anchored carryover,
Appendix B.3; randomized allocation), and keep the driver's behavior as
options, so the stress-test DGP is reproduced by
`carryover = "recursive"` with `allocation = "fixed"`. Accordingly,
`verify.R` checks anchored carryover against the closed form and
randomized allocation against its sampling distribution, not against
the driver, and checks the `recursive` option against the driver's
means.

## Core patch (`orig-to-refined.patch`)

New arguments have defaults that select the refined process; every
refinement can be switched off to recover the published process.
Signatures are otherwise unchanged.

| Hunk (orig lines) | File | Change | Why | Report 01 |
|---|---|---|---|---|
| `@@ -15` | generateData.R | Documents the optional `modelparam` columns `rho` and `c1` | | `tab:spec` |
| `@@ -40` | generateData.R | Documents `coupling`, `corstructure`, `carryover`, `nonPD` | | `tab:spec` |
| `@@ -47` | generateData.R | New arguments `coupling = "graded"`, `corstructure = "separable"`, `carryover = "anchored"`, `nonPD = "stop"`, checked with `match.arg` | Defaults select the refined process; the alternatives restore orig | `tab:refinements` |
| `@@ -81` | generateData.R | Anchored carryover: off drug after exposure the br mean is the last on-drug mean times `(1/2)^(tsd/t_half)` (B.3 code); `"recursive"` keeps orig lines 84 to 92 | The recursion decays an already decayed mean, so the decay compounds at the second and later off-drug visits | Step 3; B.3 |
| `@@ -98` | generateData.R | Separable AR(1), within factors: one kernel `A = rho^abs(t - s)` in calendar weeks for all three factors (`rho` from `modelparam$rho`, default 0.7), replacing `c.tv`, `c.pb`, `c.br` | Compound symmetry overstates long-lag correlation; one rho is required for the Kronecker form | Steps 4 and 5; `tab:spec` |
| `@@ -116` | generateData.R | Separable AR(1), across factors: `c1` at the same visit and `c1 * rho^abs(t - s)` across visits (`c1` from `modelparam$c1`, default 0.2), replacing `c.cf1t` and `c.cfct` | The factor block is K x A, positive definite for any schedule; applied to the correlation because the pb SD scales with expectancy | Step 5; "Why separability" |
| `@@ -132` | generateData.R | Graded coupling: `c.bm * g`, with `g = 1` on drug and `(1/2)^(tsd/t_half)` off drug after exposure; 0 before exposure | Step coupling removes the interaction at any positive half-life | Step 2; B.1 |
| `@@ -142` | generateData.R | `nonPD = "stop"`: stop, with c.bm and the smallest eigenvalue in the message, when the correlation matrix is not positive definite; `"repair"` keeps orig's `make.positive.definite()` path | Silent repair changes the correlations, and so the effect | Step 1; B.2 |
| `@@ -43` | generateSimulatedResults.R | Documents `dgpoptions` and `allocation` | | |
| `@@ -65` | generateSimulatedResults.R | New arguments `dgpoptions = list()` (passed to `generateData()`) and `allocation = "randomized"` | Lets a simulation run the published process | |
| `@@ -131` | generateSimulatedResults.R | Randomized allocation: within each replicate each participant is assigned to a path with probability 1/nP (simple randomization), so the counts per path vary between replicates; an allocation with an empty path is redrawn. `"fixed"` keeps orig's split and its replicate labeling | `tab:spec`: randomized to paths (RH: fixed numbers) | `tab:spec` |

The published process is

```r
generateSimulatedResults(..., allocation = "fixed",
  dgpoptions = list(coupling = "step", corstructure = "published",
                    carryover = "recursive", nonPD = "repair"))
```

and with these settings the patched code reproduces orig bit for bit
(verified, Section 4 below). The stress-test DGP is
`dgpoptions = list(carryover = "recursive")` with
`allocation = "fixed"` and `modelparam` columns `rho = 0.7`,
`c1 = 0.2`; its compound-symmetry curve is
`corstructure = "published"` with graded coupling.

**Re-exposure is unchanged.** Time on drug restarting at zero when the
drug is resumed comes from orig's `buildtrialdesign()`, which the patch
does not touch (`tab:spec`: retained because re-titration justifies it
for prazosin). Combined with anchored carryover this gives the downward
step in the mean at the first re-exposure visit described in B.3.

### Design choices forced by orig's structure

- **Coupling gate.** Orig assigns the coupling only where the br mean is
  nonzero. That set is exactly the set where the graded weight is
  positive (on drug, or off drug after exposure with a positive
  half-life), so the patch keeps orig's test and multiplies by the
  weight, instead of rewriting the block as in B.1. The two agree except
  in the degenerate case of a drug response with maximum 0, where B.1
  would still couple and the patch, like orig, does not. Keeping the
  test also keeps the hunks for separability and coupling apart.
- **Separability inside orig's loops.** The within- and cross-factor
  loops are kept and only the values change, so that
  `corstructure = "published"` is orig's code path.
- **What the stop tests.** The correlation matrix, not the covariance:
  the pb SD is the expectancy times its SD, so the covariance is
  singular whenever an expectancy is 0, without anything being wrong.
  With `nonPD = "stop"`, MASS::mvrnorm's own tolerance (which silently
  truncates small negative eigenvalues) is never reached.
- **Where allocation happens.** Orig allocates in
  `generateSimulatedResults()` by generating N times Nreps participants
  per path at once and labeling them by replicate; the randomized
  version draws the per-replicate counts first and labels accordingly.
  The redraw on an empty path is needed because orig's `lme_analysis()`
  offsets participant IDs by the maximum of the previous path, which
  fails for an empty path. It conditions the randomization on every path
  being used; at N = 70 over four paths the probability of an empty path
  is about 1e-8.
- **Repair paths in orig.** The only repair is
  `make.positive.definite(sigma, tol = 1e-3)` in `generateData()` (orig
  lines 145 to 150; searched in all of `R/`). It is kept under
  `nonPD = "repair"`.

## Options patch (`orig-to-refined-options.patch`)

All five options of the stress test are included; none was tangled with
driver code that has no counterpart in orig. All are off by default
and are set by optional `modelparam` columns (or a `censorparam`
column for dropout).

| Hunk | File | Option | Stress-test cell |
|---|---|---|---|
| `@@ -21`, `@@ -74` | generateData.R | Documentation; reads the option columns | |
| `@@ -116` | generateData.R | Weibull decay: `decay.k` (1 = exponential, the default) in the br mean, both anchored and recursive, through `carryoverdecay()` | weibull0.5, weibull2 |
| `@@ -185` | generateData.R | Weibull decay in graded coupling; prognostic couplings `c.bm.tv`, `c.bm.pb`, `c.bm.bl` (correlation of the biomarker with tv and pb at every visit, and with BL) | prog-TV-const, prog-PB, prog-BL |
| `@@ -220` | generateData.R | Mean moderation of tv, `b.bm.tv`: tv shifted by `b.bm.tv * sd_tv * z * G(t) / max G` | prog-TV-grow |
| `@@ -245` | generateData.R | Measurement error `me.ol`, `me.bl`: independent N(0, SD^2) on the observed total, SD by phase (expectancy 1 = open label) | meas-error (3 and 6) |
| `@@ -252` | generateData.R | `carryoverdecay()`; with `k = 1` it keeps orig's `(1/2)^(tsd/t_half)` expression, so the default is bit-identical | |
| `@@ -22`, `@@ -43` | censordata.R | Optional `type` column: `"MCAR"` or `"MAR"`, rate `beta0` | dropout-MCAR, dropout-MAR |
| `@@ -85` | censordata.R | `censordropout()`: monotone dropout from the third visit; MCAR uniform visit; MAR hazard `plogis(a + z)` on the previous standardized score, standardized and calibrated within each path of each replicate, as the driver does per path and trial | |

Notes. Measurement error is added to the observed totals only, not to
the components or to orig's `D_` columns, so orig's own (published)
censoring mechanism, which reads the `D_` columns, ignores it. MAR
dropout needs the `path` and `replicate` columns that
`generateSimulatedResults()` adds; without them the whole data set is
one group.

## Analysis patch (`orig-to-refined-analysis.patch`)

`lme_analysis()` keeps its signature. A new element of `op`,
`op$analysis` (a column of `analysisparams`), selects the analysis;
the default is `"F1"`, and `"orig"` runs orig's code unchanged. The
other codes call the new internal function `lme_analysis_refined()`.
Two further options: `op$t_half_a` (assumed half-life, default 0.5
week) and `op$t_half_grid` (AIC grid, default 0, 0.25, 0.5, 1, 2).

| Hunk | File | Change |
|---|---|---|
| `@@ -18` | lme_analysis.R | Documents `op$analysis`, `op$t_half_a`, `op$t_half_grid` |
| `@@ -49` | lme_analysis.R | Dispatch: any analysis other than `"orig"` goes to `lme_analysis_refined()` |
| `@@ -167` | lme_analysis.R | New `exposuredecay()` (the exposure-decayed indicator with the pre-exposure guard of B.4) and `lme_analysis_refined()` |
| `@@ -38` | generateSimulatedResults.R | Documents the `analysis` and `population_fit` columns |
| `@@ -153` | generateSimulatedResults.R | The pooled fit (`ETbeta`, `ETbetaSE`) by `population_fit` (below) |
| `@@ -180`, `@@ -201` | generateSimulatedResults.R | Records the analysis code, the test's df and the kind of pooled fit in the results (`analysis`, `df`, `ETbetaSE_type`) |
| `@@ -20` | DESCRIPTION | `Imports: nlme`; `Suggests: clubSandwich, mmrm` |

All stress-test analyses use the long data with the baseline as a
t = 0 row (t in weeks), the biomarker centered at its trial mean
`bmc = bm - mean(bm at t = 0)` (centered, not scaled), and the
exposure-decayed indicator `Dbc`: 1 on drug,
`(1/2)^(tsd/t_half_a)` off drug after first exposure, 0 at baseline
and before first exposure (`Dbc = Db` when `t_half_a = 0`). The core
mean model is `Sx ~ bmc + t + Dbc + bmc:Dbc`; the tested coefficient
is `bmc:Dbc`, two-sided at 0.05. Specifications are those of the
driver (`fit_all`, `fit_extra`, `fit_extra2`) and of the compendium
table `tab:invd-analyses` (`analysis/report/compendium-mathematics.tex`).

| Code | Mean model beyond the core | t_half_a | Covariance | Inference | df | Role |
|---|---|---|---|---|---|---|
| orig | orig: `Sx ~ bm + Db + t + bm*Db` (+ `De` if `useDE` and the expectancy varies), uncentered bm, binary Db, `lmer` | none | random intercept (or intercept and slope) | lmerTest t | Satterthwaite | the published analysis |
| S | none | 0.5 | RC | CR2 | Satterthwaite (CR2) | provisional reference |
| DS | none | AIC | RC | CR2 | Satterthwaite (CR2) | reference, AIC half-life |
| C1 | none | 0.5 | RCH | CR2 | Satterthwaite (CR2) | alternative |
| C2 | none | 0.5 | RCH | model-based | containment | alternative |
| F1 | `De + bmc:De` | 0.5 | RC | CR2 | Satterthwaite (CR2) | selected; default |
| F1A | `De + bmc:De` | AIC | RC | CR2 | Satterthwaite (CR2) | post hoc |
| G1 | `De + bmc:De` | 0.5 | RCH | CR2 | Satterthwaite (CR2) | post hoc |
| G2 | `De + bmc:De` | 0.5 | RCH | model-based | containment | post hoc |
| F2 | `bmc:t` | 0.5 | RC | CR2 | Satterthwaite (CR2) | alternative |
| ES | `bmc:bl0` | 0.5 | RC | CR2 | Satterthwaite (CR2) | alternative |
| A | none | 0.5 | US | Kenward-Roger (linear) | Kenward-Roger | alternative |
| DA | none | AIC | US | Kenward-Roger (linear) | Kenward-Roger | alternative, AIC half-life |
| EA | `bmc:bl0` | 0.5 | US | Kenward-Roger (linear) | Kenward-Roger | alternative |
| R | none, `Dbc = Db` | 0 | CS | model-based | Satterthwaite (mmrm) | published-type reference |

RC: `nlme::lme`, random intercept, `corCAR1(form = ~t | ptID)`, REML,
`lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)`. RCH: RC
with `varIdent(form = ~1 | phase)` over baseline, open label and
blinded. US and CS: `mmrm` `us(visit | ptID)` and `cs(visit | ptID)`,
no random intercept. CR2: `clubSandwich::coef_test(vcov = "CR2",
test = "Satterthwaite")`, clustered on participant. Model-based: the
`nlme` Wald t on its containment df. AIC: the half-life minimizing the
ML AIC over the grid, refitted by REML (Kenward-Roger for DA). De is
the expectancy (0 at baseline), bl0 the baseline indicator.

G1 and G2 (phase-specific residual variances) were reported in the
stress-test documents as effectively equivalent to F1; they are
included because they were easy to add, not re-examined here.

Dependencies: nlme (recommended, always installed) is imported;
clubSandwich (CR2: S, DS, C1, F1, F1A, F2, ES, G1) and mmrm (A, DA, EA,
R) are suggested, and an analysis that needs a missing one stops with
an error naming the package. A fit that fails returns `NA` estimates
with the error message in `warning`, as the driver records a failed fit.

The long data are assembled in the driver's row order (by path, then
visit with the baseline first, then participant), so that the fits are
numerically identical to the driver's, not merely equal in theory.

### The pooled fit (`ETbeta`) and `population_fit`

Orig's `generateSimulatedResults()` fits the analysis once to all
N times Nreps simulated participants before censoring and stores the
estimate and SE as `ETbeta` and `ETbetaSE` in every result row. The
article uses this point estimate as the reference for the effect of
censoring (Hendrickson et al. 2020, Results, "Variability and Bias in
Effect Size Estimates"): "the β from each replicate was compared to
the β obtained when the model was applied to all simulated
participants across all replicates in the absence of censoring
(Figure 7B)". The Figure 7 caption defines the "gold standard"
coefficient as "the coefficient calculated across all simulated
participants from all replicates with no censoring". Only the point
estimate is used; the SE is not.

The analysis patch keeps the pooled fit and adds a column
`population_fit` to `analysisparams`:

| Value | Pooled fit | `ETbetaSE_type` |
|---|---|---|
| `"full"` | the analysis exactly as applied to each replicate; default for `"orig"`, which is orig's behavior (same call, same `ETbeta`, `ETbetaSE`, so Figure 7B reproduces) | `"full"` |
| `"model"` | the analysis's own mean and covariance model (and, for DS, F1A and DA, its AIC choice of half-life) with model-based inference only: the `nlme` Wald SE instead of CR2, `mmrm` `method = "Residual"` instead of Kenward-Roger or Satterthwaite; default for the nlme-based analyses (S, DS, C1, C2, F1, F1A, G1, G2, F2, ES) | `"model"` |
| `FALSE` | no pooled fit; `ETbeta`, `ETbetaSE` are `NA`; default for the mmrm-based analyses (A, DA, EA, R) | `NA` |

`"model"` changes only the inference, so its point estimate is the
estimate the full analysis gives on the same data. One exception was
found: on the pooled data of 7,200 participants the stress test's
`lmeControl(opt = "optim")` fails numerically ("NA/NaN/Inf in foreign
function call"), and so would the `"full"` fit. The `"model"` fit then
refits with nlme's default optimizer (nlminb), which reaches the same
REML optimum to optimizer tolerance: on 1,440 participants, where both
converge, the two optimizers agree to 1e-7 in the estimate. The
fallback is applied only to the pooled `"model"` fit, never to the
per-replicate analyses, which keep the driver's behavior (a failed fit
is `NA`).

**Timing and the defaults** (host session, single core, Hybrid design,
refined DGP, c_bm 0.3, t_half 0.5 week). The `"model"` pooled fit of F1
took 6.7 s at 1,440 participants (20 replicates of N = 72) and 11 to
27 s at 7,200 (100 replicates; 27 s when the optim attempt failed and
the fallback ran), well under
the 5-minute threshold, so `"model"` is the default for the
nlme-based analyses. The mmrm analyses do not scale: the `"model"` fit
of R took 18.6 s at 1,440 participants, 219 s and 3.7 GB of memory at
2,880, and the process was killed (out of memory) at 7,200 for R and
for A. They therefore default to `FALSE`. The DS and F1A pooled fits
consist of six fits (five ML fits for the AIC choice and one REML fit)
and were not timed at 7,200 participants; by the F1 figure, about
2 to 3 minutes is expected (inferred, not measured).

## Verification results

Full run of `verify.sh` on 2026-10-08, repeated on 2026-10-09 after the
`population_fit` change (R 4.6.1, host session, renv not active): all
checks PASS, 0 failed, in both runs. The numbers:

- **Patch application.** Each patch passes `git apply --check` and
  applies with `patch -p1` without fuzz or rejects, in the trees core,
  options (core + options), analysis (orig + analysis) and full (all
  three).
- **Refined DGP against the driver** (base cell: Hybrid, N = 70,
  t_half 0.5 week, rho 0.7, c1 0.2, c_bm 0.3; also the CO design).
  Correlation matrix: max abs difference 1.1e-16; covariance 2.2e-16
  (entries up to 342); graded coupling vector against the closed form
  c_bm g_t: 0, against the driver 1.7e-18; coupling before first
  exposure (CO placebo-first path): exactly 0; tv and pb means 1.8e-15.
  The driver hard-codes the published parameters rounded to five
  decimals (differences up to 1e-4); the comparisons use the same rounded values.
- **Anchored carryover** equals the closed form
  `mu_last * 2^(-t_sd/t_half)` exactly (difference 0). With
  `carryover = "recursive"` the means equal the driver's (1.8e-15).
  Anchored and recursive means differ by up to 0.478 points (Hybrid)
  and 0.0096 (CO) at this half-life.
- **CS curve.** `corstructure = "published"` with graded coupling
  equals the driver's compound-symmetry construct (rho 0.8): 1.8e-15.
- **Ceiling and stop.** The Hybrid ceiling of the base cell is 0.5182.
  At c_bm 0.5172 the generator runs (driver smallest eigenvalue
  +0.0016); at 0.5192 it stops (driver -0.0016). `nonPD = "repair"`
  restores one repair at c_bm 0.9.
- **Sampler** (20,000 participants per path): empirical
  Cor(bm, BR_t) against c_bm g_t, max |z| 2.30 over 32 visits;
  BR means against the anchored closed form, max |z| 1.72.
- **Randomized allocation** (1,000 replicates, N = 70): every
  replicate has 70 participants and no empty path; path shares 0.2494,
  0.2481, 0.2507, 0.2517 (SE 0.0016); per-replicate count variances
  13.06, 13.74, 13.43, 12.63 against the binomial 13.125; 769 distinct
  count vectors. Fixed allocation gives 18, 18, 17, 17 every time.
- **Published process reproduced.** With the published options,
  `generateData()` output is identical to orig's (5 cells by 4 paths,
  in all three patched trees); orig repairs at c_bm 0.3, t_half 0 (4 of
  4 paths) and not at 0.25; `generateSimulatedResults()` with the
  published DGP, fixed allocation and `analysis = "orig"` gives
  identical ETbeta, beta, SE and p in all 15 result rows (3 replicates
  by 5 censoring sets). A smoke test with the refined defaults (F1, S,
  orig; N = 35, published censoring) runs with no failed fits.
- **Pooled fit** (`population_fit`, run of 2026-10-09). With the
  default for `"orig"` (`"full"`), ETbeta and ETbetaSE are identical to
  unpatched orig's. In an end-to-end run the defaults give type
  `"model"` for F1 and S (ETbeta populated) and `"full"` for orig;
  `FALSE`, and the default of the mmrm analysis R, give `NA`. On a
  pooled data set of 1,440 participants (20 replicates of N = 72) the
  `"model"` and `"full"` point estimates are identical (absolute
  difference 0) for F1, C1, DS (same AIC half-life, 0.5), A and R; the
  SEs differ as expected (F1: model 0.0201, full CR2 0.0197). The
  `"model"` fit of F1 on 7,200 participants took 11.1 s in this run
  (27 s in a separate run on other data in which the optim control
  failed and the nlminb fallback was used).
- **Analyses.** On 5 data sets from the refined generator (3 complete,
  one with MAR and one with MCAR dropout), every one of the 14
  stress-test analyses gives estimates, SEs, p-values and selected
  half-lives identical to the driver's code (max abs difference 0).
  Mean df over the data sets: S 21.2, F1 20.5, C2 528, A 128, R 530.
  `analysis = "orig"` is identical to unpatched `lme_analysis()`
  (useDE FALSE and TRUE). Dbc on the CO placebo-first path is
  0 0 0 0 1 1 1 1; without the guard it would be 1 at every visit.
- **Options.** Matrices and recursive means equal the driver's in the
  cells prog-TV-const, prog-PB, prog-BL, weibull0.5, weibull2, thalf0,
  thalf1 (R 1.1e-16, Sigma at most 8.9e-16, means 1.8e-15); the tv
  shift of `b.bm.tv` equals the driver's formula (5.3e-15);
  measurement-error SDs 3.03, 3.00 (open label) and 5.98 to 6.03
  (blinded), max |z| 2.25, components unchanged; MCAR and MAR dropout
  visits identical to the driver's `dropout_visit()` under the same
  seed (5 seeds each); dropout monotone from the third visit; rates
  0.2019 (MCAR) and 0.1994 (MAR), SE 0.0033.

The analysis equivalence was obtained by handing the driver's
`make_trial()` the same participants (its `draw()` and
`dropout_visit()` replaced); its data assembly and fitting code ran
unchanged.

## What is not included or not verified

- **buildtrialdesign() and re-exposure.** Unchanged by design.
- **Rd files and NAMESPACE.** `man/` and `NAMESPACE` are not
  regenerated; the new functions (`exposuredecay`,
  `lme_analysis_refined`, `carryoverdecay`, `censordropout`) are
  internal. A package build would need `roxygen2` to be rerun; no
  `R CMD check` of the patched tree was run.
- **Vignettes and plotting functions.** Unchanged. With the new
  defaults, orig's publication vignettes would run the refined process
  and F1; they reproduce the publication only with the published
  options above and `analysis = "orig"`.
- **Designs without within-participant variation in Dbc** (the
  open-label design). Orig's analysis switches to a biomarker-by-time
  model there; the stress-test analyses have no such fallback and are
  not identifiable (F1's `De` is then collinear with `Dbc`). Not tested.
- **The pooled fit with the mmrm analyses.** No pooled fit by default
  for A, DA, EA and R (see "The pooled fit" above); a pooled `"model"`
  or `"full"` fit for them is possible on request but was not timed at
  7,200 participants (it was killed for memory).
- **Draws are not identical to the driver's.** The random streams
  differ (orig draws per path in one block), so the sampler is checked
  distributionally; matrices, means, couplings, dropout visits and
  analyses are checked exactly.
- **Designs checked.** Matrix checks cover the Hybrid and crossover
  designs, analysis checks the Hybrid design; OL and OL+BDC were not
  run.
- **orig's own defects** that do not bear on the refinements (for
  example, `lme_analysis()` failing when `op` is missing, and the
  one-sided assignment of the same-visit cross-factor correlation,
  which is symmetric only because the loop visits both orders) are left
  as they are.
- **The 2024 upstream commits are not part of the base.** The patches
  apply to `58b32a9`, the code behind the published results, not to
  current `master` (`06dac83`). The 2024 changes are superseded rather
  than merged:
  - `8609f12` (scale factor on the compounding recursion; the ratio
    coupling rule marked `RON THOMAS VERSION`) is replaced by anchored
    carryover and direct graded coupling. The ratio rule is correct
    only while the means compound (`docs/53`, Section 9.5; `docs/55`).
  - `325314f` (a continuous `Dbc` in `lme_analysis()`) and its fixes in
    PR #2 (`1cbdb39`, the `op$carryover_scalefactor` typo) and PR #3
    (`78c8635`, the coefficient name `bm:Dbc`) are replaced by the
    analysis patch's own `exposuredecay()`, which also codes
    pre-exposure visits as unexposed (report 01, Appendix B.4).
  - `f6ee86d` (cross-factor index typo) does not arise under the
    separable structure, which builds that block directly.

  A version of these patches against current `master` would need to be
  rebased and verified separately; it has not been produced.
