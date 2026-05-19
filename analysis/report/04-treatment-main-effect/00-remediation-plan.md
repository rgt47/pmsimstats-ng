# Remediation Plan: Paper 04 -- Treatment Main Effect
*2026-05-19 07:39 PDT*

This document records the findings of a deep review of `report.Rmd` against
the Morris et al. (2019) ADEMP framework and establishes a prioritized
remediation plan. The review identified eight distinct problems, which we
organize below by category. Each item carries a priority label: **Critical**
(the paper cannot be submitted without it), **Major** (it affects the
validity or interpretability of results), or **Minor** (it affects code
quality and reproducibility but not the conclusions).

We note at the outset that the report currently contains two independent
simulation implementations that have not been reconciled: an older section
(approximately lines 394-476 of the current Rmd, using DerSimonian-Laird
meta-analysis and n_sim = 1000) and a newer "Prototype Monte Carlo
confirmation" section (approximately lines 1189-1237, using `nlme::lme`
with `corCAR1` correlation and n_sim = 1500). The abstract draws from the
prototype section; the Methods narrative describes the older implementation.
This inconsistency is the root cause of several of the problems that follow
and must be addressed first.

---

## Problem 1. Dual simulation implementations (Critical)

The paper presents two simulation frameworks that share a parameter space
but differ in analysis model, sample size, and performance reporting. The
abstract and ADEMP estimand block describe `nlme::lme` with `corCAR1`; the
Methods narrative and the main simulation table describe DerSimonian-Laird
meta-analysis on per-patient period means (Chen's "M4" estimator, known to
inflate Type I error at small n). A reader cannot determine which results
are authoritative.

**Remediation.** Retire the DerSimonian-Laird section entirely. Promote the
`lme`/`corCAR1` prototype as the single canonical implementation. The
prototype section should be expanded into a full Methods and Results
narrative; the old section should be removed. Once this consolidation is
complete, the abstract, Methods, and results tables will all describe the
same estimator.

---

## Problem 2. Seed and RNG discipline (Critical)

Two failures exist in the current implementation.

First, `set.seed(42)` appears to be called inside the simulation function
rather than once at program entry. Per-function seeding destroys the ability
to reproduce any individual replicate in isolation, and it confounds the
RNG state across parallel workers.

Second, the code does not specify `RNGkind('L'Ecuyer-CMRG')` before the
seed is set. This is required for correct parallel stream generation (Morris
§4.6) and for compatibility with `parallel::mclapply` or any future move to
parallel execution.

**Remediation.** At the top of the simulation driver script (or, once the
code is moved to `R/`, in the simulation entry-point function), add:

```r
RNGkind("L'Ecuyer-CMRG")
set.seed(2024L)
```

Remove any `set.seed()` calls from inside simulation functions. Capture the
per-replicate `.Random.seed` state at the start of each replicate so that
any anomalous result can be reproduced exactly. A minimal pattern is:

```r
seed_log[[i]] <- .Random.seed
result[[i]]   <- run_one_replicate(params)
```

---

## Problem 3. Non-uniform and insufficiently justified n_sim (Critical)

The current report uses at least three different values for n_sim: 1000 in
the primary simulation, 100 in the scenario sweep, and 500 in the Type I
error evaluation. Morris §4.4 recommends choosing a single n_sim before
examining any results, justified by the worst-case Monte Carlo SE across all
performance measures of interest.

The n_sim = 100 used for the scenario sweep is particularly problematic. At
that n, the Monte Carlo SE on a power estimate near 0.80 is approximately
4.0%; at power near 0.50 it is approximately 5.0%. Neither is adequate for
the comparative claims made about relative efficiency across conditions. The
report acknowledges the limitation in a footnote, but the sweep results are
nonetheless presented as evidence for design guidance.

**Remediation.** Establish a single pre-specified n_sim for all conditions
(primary, scenario sweep, and Type I) before production runs begin. We
recommend n_sim = 2000, which gives a worst-case MCSE of approximately
1.1% on a 75% power estimate and approximately 0.5% on a 5% Type I
estimate. If the scenario sweep remains computationally prohibitive at
n_sim = 2000 across all 540 cells, reduce the number of cells by collapsing
less informative dimensions rather than reducing n_sim. Document the
pre-specified n_sim in the ADEMP pre-registration file described in
Problem 8 below.

---

## Problem 4. Analysis model inconsistency with declared estimand (Major)

The ADEMP estimand block declares the primary estimand as $\hat{\beta}_D$
from `nlme::lme` with `corCAR1`. The abstract states that the analysis model
for the RCT arm is "a baseline-adjusted linear model." The Methods narrative
(approximately line 429 of the current Rmd) instead refers to independent
samples t-tests. These three descriptions cannot simultaneously be correct.

Furthermore, the prototype section notes that the N-of-1 arm exhibits a
Type I error rate of 0.061 (approximately 1.8 Monte Carlo SEs above the
nominal 0.05) and coverage in the range 0.935 to 0.949. Both observations
are consistent with small-sample Wald test inflation, which has a
well-understood remedy.

**Remediation.** For the RCT arm, replace any t-test reference with ANCOVA
(baseline-adjusted linear model), as the abstract already states. For the
N-of-1 arm, apply the Kenward-Roger or Satterthwaite small-sample degree-of-
freedom correction to the `lme` Wald test. This is available via
`pbkrtest::KRmodcomp()` or via the `lmerTest` package. Re-run the prototype
at the full production n_sim after this correction and verify that Type I
coverage returns to the nominal level before finalizing results. If coverage
remains below 0.95 after the correction, this should be reported as a
finding rather than suppressed.

---

## Problem 5. Estimand alignment with Gompertz DGP (Major)

The prototype section notes that the observed bias in $\hat{\beta}_D$ grows
with the magnitude of the true effect: from approximately -0.20 at the null
to approximately -2.12 at a true effect of -2.0. The current text flags this
as a DGP artifact (the Gompertz trajectory is not linear at the indicator
value) but does not formally declare it as part of the estimand definition.

This matters because a bias that is a deterministic function of the true
parameter value is an estimand alignment problem, not a stochastic variance
problem. The estimator is not estimating the declared quantity in the
presence of the Gompertz trajectory. Presenting the estimator as
"essentially unbiased" under these conditions is, at minimum, misleading.

**Remediation.** Revise the estimand block to explicitly state the target
quantity under the Gompertz DGP. If the intended estimand is the linear
approximation to the Gompertz slope at a specified indicator value, state
this, and explain why the linear approximation is the quantity of clinical
interest. Report the bias at each effect-size cell as a substantive finding
rather than a nuisance. If the bias is large enough to affect the practical
interpretation of the results, this warrants discussion as a limitation.

---

## Problem 6. Performance measures table incomplete (Major)

The Morris (2019) summary table currently in the report is missing the
following recommended performance measures:

- 95% Wald coverage of $\hat{\beta}_D$ (mentioned in the prototype notes
  but absent from the table)
- Monte Carlo SE of the coverage estimate
- Model-based SE (average of the estimated standard errors across
  replicates, to assess the calibration of the analytical SE)
- Relative error of the model SE (ratio of model SE to empirical SE,
  minus one)

Without coverage and model SE, the table cannot confirm that the
confidence interval procedure has the correct operating characteristics.
The prototype notes suggest coverage is already slightly below 0.95 in
at least some conditions; this should be prominently reported.

**Remediation.** Add a column for empirical coverage, its Monte Carlo SE,
average model SE, and relative error to the performance table. Compute
these from the production simulation output. If the table becomes unwieldy
for print, restrict the displayed rows to primary design cells and provide
full results in a supplementary table.

---

## Problem 7. Simulation code infrastructure (Major)

The simulation code currently lives in a workspace file (`sim_workspace.RData`)
loaded at the top of the Rmd. This approach has several practical problems:
the workspace is not version-controlled in a readable form, it cannot be
tested with `testthat`, and it cannot be called from other papers in the
pmsimstats-ng package without duplication.

Additionally, stub placeholder variables at the top of the Rmd (approximately
lines 47-64) indicate that the integration between the Rmd and the simulation
output is not complete. A chunk that assigns `n_sim <- NULL` and `results <- NULL`
and then uses those variables in a results table will silently produce empty
output rather than an error.

The report also uses `%>%` in at least one chunk (approximately line 940),
inconsistent with the `|>` standard established in the project's CLAUDE.md.

**Remediation.** Move all simulation driver and helper functions into
properly documented `R/` package functions following the pmsimstats-ng
package infrastructure. Replace workspace loading with a call to a package
function (or, at minimum, `source()` of a tracked `.R` file). Replace all
stub placeholder assignments with `stopifnot()` or `assertthat::assert_that()`
guards that fail loudly if the expected objects are not present. Replace all
`%>%` with `|>`.

---

## Problem 8. Absent ADEMP pre-registration file (Minor)

Morris et al. (2019) recommend, and the pmsimstats-ng project's own ADEMP
audit notes, that the simulation plan be fixed in a pre-registration document
before production runs begin. No such file exists for paper 04.

**Remediation.** Create `00-ademp-pre-registration.md` in the
`04-treatment-main-effect/` directory. The file should specify, in their
final form before any production run: the Aims; the Data-generating
mechanisms and their full parameter grid; the Estimands (formal definition
under the Gompertz DGP); the Methods (lme/corCAR1, ANCOVA for RCT, KR
correction); the Performance measures (power, Type I, coverage, bias,
model SE, relative error); and the pre-specified n_sim with its MCSE
justification. Once production runs begin, this file should not be modified;
any post-hoc sensitivity analyses must be clearly labeled as exploratory.

---

## Summary and sequencing

We may organize the remediation work into three phases:

**Phase 1 (prerequisite for any further work).** Problems 1, 2, and 3.
Consolidate to a single simulation implementation, fix RNG discipline, and
establish a pre-specified n_sim. Until these are resolved, subsequent work
on performance reporting (Problem 6) will produce results that may change.

**Phase 2 (correctness).** Problems 4 and 5. Verify estimand alignment,
apply the KR correction, and re-run the production simulation with the
corrected model and full n_sim. Confirm that Type I and coverage are at
nominal levels before reporting final results.

**Phase 3 (completeness).** Problems 6, 7, and 8. Expand the performance
table, move the code into `R/`, and write the pre-registration file.
Although the pre-registration file is formally a prerequisite for Phase 1,
it is listed here because it can be drafted after the consolidation (Problem 1)
clarifies what the canonical plan actually is.

In conclusion, three points are to be emphasized. First, the consolidation
of the two simulation frameworks (Problem 1) is the most consequential
single action: it resolves the contradiction between the abstract and
Methods, removes an entire estimator that is known to have poor small-sample
properties, and provides the foundation for all subsequent work. Second, the
KR correction (Problem 4) is the most likely change to alter a reported
numerical result, as the current prototype data suggest non-nominal Type I
and coverage in the N-of-1 arm. Third, the pre-specified n_sim and ADEMP
pre-registration (Problems 3 and 8) are the most important changes for
long-term credibility of the simulation findings, since they prevent
post-hoc selection of conditions that favor a desired conclusion.

---

*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/04-treatment-main-effect/00-remediation-plan.md*
