# pmsimstats-ng Project Guidelines

Project-specific conventions, architectural patterns, and workflow
guidance for the pmsimstats-ng R package. General R coding style
is inherited from the user-level CLAUDE.md; this file records only
what is specific to this project.

Canonical notation and vocabulary: `analysis/report/NOTATION.md`
(revised 2026-10-09). Where this file and NOTATION.md differ,
NOTATION.md governs.

## Project overview

**Purpose.** Monte Carlo simulation-based power analysis for
aggregated N-of-1 clinical trial designs, targeting the detection
of biomarker-treatment interactions in the presence of carryover.

**Key architectural feature.** Two interaction architectures
(formerly labeled "Architecture A/B"; the letters are retired,
NOTATION.md Part 3):

- **The covariance-moderation architecture** (`mvn`). The
  biomarker-response interaction is encoded through the biomarker's
  correlation with the drug response in the multivariate-normal
  covariance (Hendrickson's formulation).
- **The mean-moderation architecture** (`mean_moderation`). The
  interaction is encoded by shifting the drug-response mean by the
  centered biomarker.

At matched total N = 70 in the Hybrid design, the relative power loss
under carryover is 15.4% for covariance moderation against 6.9% for
mean moderation. This figure was computed for a pre-restatement
version of paper 01 (2026-08-08 N-matching fix; verified then,
NOTATION.md Part 2, "In paper 01 only the Hybrid column moved...").
Paper 01's scope was narrowed on 2026-10-07 to the covariance-
moderation DGP specification only, and this cross-architecture
comparison was removed from its text along with the rest of the
cross-architecture material; a 2026-10-07 referee pass
(`revisions/33-paper01-dgp-architecture-review.md`) had also flagged
the original comparison as using non-comparable absolute-power
baselines across architectures (Section 3.3 there), a problem paper
01's current Discussion resolves qualitatively via a conditioning/
rank-one argument but without rerunning the matched-N numbers under
it. The 15.4%/6.9% figure is therefore **not currently narrated in
any live manuscript**: it survives only in NOTATION.md's historical
record above, the archived draft
(`analysis/report/01-dgp-mean-moderation-vs-mvn/archive/report-2026-10-07-pre-restatement.Rmd`),
and `docs/series/ch07-architectures.md` (which cites it onward from
"paper 01" rather than from a live source, and runs no simulation of
its own). Cite NOTATION.md Part 2 directly for this figure, not
"paper 01," until it is rerun and re-narrated somewhere live.
Conditional on the biomarker, covariance binding is mean moderation
with a smaller residual variance; the divergence reflects the
exposure profile of the interaction (graded coupling fades off
drug), not a difference in kind (compendium, sec:architectures).
The earlier "40-60% vs 8-13%" contrast is superseded.

**Named processes.** "The published process" is Hendrickson's code at
commit 58b32a9; "the refined process" is the recommended DGP (graded
coupling, separable AR(1) with rho = 0.7 and c_1 = 0.2, anchored
carryover, a stop on non-positive-definite matrices). Path allocation
is RH's original fixed-numbers approach, unchanged; randomization of
path allocation is not part of the refined process's definition
(2026-10-21). See NOTATION.md "Named data-generating processes" and
paper 01 Table tab:spec.

## Repository layout

```
pmsimstats-ng/
|- R/                   Package source (exported API)
|- man/                 Roxygen-generated documentation
|- tests/               tinytest suite
|- vignettes/           Package vignettes
|- inst/                Installed files
|- data/                Package datasets (*.rda)
|- NAMESPACE, DESCRIPTION
|
|- analysis/
|  |- scripts/          Analysis drivers, one subdir per topic
|  |  |- quick-sim/              Closed-form, stress-test and covariance-study drivers
|  |  |- carryover-sensitivity/   Paper 02 driver
|  |  |- treatment-main-effect/   Paper 04 driver and ADEMP working protocol
|  |  |- compendium-mathematics/  Check scripts for the mathematics compendium
|  |  |- figure4/                 Hendrickson (2020) Figure 4 reproduction
|  |  \- parity/                  Cross-implementation parity tests + baselines
|  |- patches/          orig-to-refined patches for commit 58b32a9
|  |- data/             Shared simulation outputs
|  |- remote_runs/      simrun job results (gitignored)
|  \- report/           Manuscripts, one numbered subdir per paper
|     |- 01-dgp-mean-moderation-vs-mvn/ ... 14-two-component-drop-tv/
|     |- compendium-mathematics.tex     Mathematics compendium
|     |- NOTATION.md                    Canonical notation
|     \- references.bib                 Shared candidate bibliography
|
|- docs/
|  |- NN-*.md            Methodological whitepapers (numbered)
|  |- *.tex             Older whitepapers (00-20; to be consolidated)
|  |- figures/          Figures for the whitepapers
|  |- references/       Literature source PDFs and bib
|  \- reviews/          Annotated review PDFs (non-regenerable)
|
|- implementations/
|  |- original/         Hendrickson reference (covariance moderation only, data.table)
|  |- original-extended/ Data.table collection, both architectures
|  |- tidyverse/        Tidyverse collection, both architectures
|  |- simple/           Reduced-complexity base-R sandbox
|  \- nof1power/        Vendored R package for paper 04 main-effect sims
|
|- CLAUDE.md, DESCRIPTION, Dockerfile, Makefile, renv.lock
```

## Documentation map

The whitepaper series lives at `docs/series/`. Start with
`docs/series/fm0-reading-guide.md` (reading paths and audience
guide) and `docs/series/fm1-executive-summary.md` (headline
results). The index of legacy documents is `docs/legacy/INDEX.md`.

**External archive.** Historical and exploratory material is kept
outside the repo at `/Users/zenn/Dropbox/prj/alz/10-pmsimstats-ng/
archive/`, including `figure5/`, `architecture_comparison/`,
`carryover_factorial/`, `analysis-archive/`, `figure4-history/`,
`gompertz-clinician-guide/`, `2025/`, and
`pmsimstats2025-2026-04-23/` (superseded pmsimstats2025 package
and its own archive, migrated here on 2026-04-23). These
directories are preserved for reference but are not reachable
from any current manuscript or driver in the repo.

## Core R package pipeline

Seven files in `R/` implement the simulation pipeline:

1. `buildtrialdesign.R` — trial design specification (paths,
   timepoints, treatment phases)
2. `generateData.R` — MVN data generation with both architectures
3. `censordata.R` — missing-data mechanisms (dropout patterns)
4. `lme_analysis.R` — LME analysis with `corCAR1` residual
   correlation
5. `generateSimulatedResults.R` — simulation orchestration and
   parallel execution
6. `carryover_analysis.R` — extended analysis (half-life
   characterization, clinical thresholds)
7. `plottingfunctions.R` — heatmap and trajectory visualization

Utility modules: `utilities.R` (cumulative effects, modGompertz),
`packagedocumentation.R`, `datadocumentation.R`.

Resolved 2026-10-21: `buildSigma()`/`generateData()`/
`generateSimulatedResults()` in `R/generateData.R` now accept
`stopOnNonPD` (default `FALSE`, preserving all existing behavior).
When `TRUE`, a non-positive-definite covariance matrix raises an
error instead of being silently repaired, both at the explicit PD
check and in the Cholesky-factorization fallback (the refined
process's recommendation, paper 01 B.2 and Table tab:spec). Callers
that want the refined process's validity guarantee must opt in
explicitly; the package default remains silent repair for backward
compatibility. Paper 02's driver core
(`analysis/scripts/carryover-sensitivity/simulation-core.R`) is the
exception: since 2026-10-26 it defaults to `stop_on_nonpd = TRUE`,
after two jobs silently repaired infeasible prognostic-confound
cells (docs/60, 2026-10-26 update); a non-PD cell there shows as a
non-converged replicate. Pass `stop_on_nonpd = FALSE` to reproduce
earlier paper 02 runs.

## Implementations

Five collections live under `implementations/`. Four are
biomarker-treatment-interaction simulators (one historical,
two production, one reduced-complexity sandbox); the fifth is
a vendored R package supporting paper 04's treatment-main-
effect comparison.

| Collection | Style | Architecture | Role |
|---|---|---|---|
| `original` | data.table | covariance moderation only | Historical reference; back-compat testing |
| `original-extended` | data.table | both | Production data.table collection |
| `tidyverse` | tidyverse | both | Production tidyverse collection |
| `simple` | base R | mean moderation only | Pedagogical sandbox |
| `nof1power` | R package | mean moderation only | Vendored simulation engine for paper 04; install with `devtools::install_local('implementations/nof1power')` |

Common parameters accepted by all collections that support them:

- `c.bm`: under covariance moderation, the biomarker-BR correlation
  ($c_{bm}$); under mean moderation, the coefficient scaling BR by
  the centered biomarker ($\beta_{bm}$).
- `carryover_t1half`: half-life of the carryover effect (weeks).
- `lambda_cor`: correlation decay rate, auto-derived from
  `carryover_t1half` if `NA`.
- `dgp_architecture`: `'mvn'` (covariance moderation) or
  `'mean_moderation'`. Absent in `original` (fixed to MVN); present
  in the others with default `'mvn'`.

Parity artifacts (`analysis/scripts/parity/parity_baseline.rds`,
`parity_report.rds`, `parity_divergences.csv`,
`test-parity-*.R`) validate numerical agreement between
collections under matching seeds.

## Naming conventions

- **Functions:** camelCase (`generateData`, `lme_analysis`,
  `buildtrialdesign`).
- **Variables:** snake_case (`c.bm`, `carryover_t1half`,
  `t_random_slope`), plus domain abbreviations:
  - `tod` — time on drug (cumulative; Gompertz input)
  - `tsd` — time since discontinuation (carryover input)
  - `tpb` — time x placebo belief accumulation
  - `Dbc` — continuous drug indicator in the analysis model
- **Data objects:** lowercase with underscores (`trial_design`,
  `resp_params`, `bl_params`, `model_params`).
- **Directories under analysis/:** hyphen-separated lowercase
  (`carryover-sensitivity`, not `02-carryover_sensitivity`).

## Statistical parameterization

- **AR(1) correlation:** $\rho^{|w_t - w_s|}$ in calendar weeks
  (the refined process); compound symmetry is written $\rho_{cs}$.
- **Carryover half-life.** Distinguish the data-generating half-life
  $t_{1/2}$ (pharmacological range 0.1-0.2 week for prazosin; the
  stress test's base setting is 0.5 week; longer values are labeled
  behavioral scenarios) from the analysis half-life $t_{1/2}^a$ used
  to build $D_{bc}$. Policy (2026-10-11): when fixed, $t_{1/2}^a =
  0.2$ week for both the interaction and the main effect. Grid for
  data-chosen half-lives in new runs: $\{0, 0.2, 0.5, 1, 2\}$ weeks.
- **Lambda decay:** `lambda_cor = ln(2) / carryover_t1half`,
  auto-derived unless overridden for sensitivity.
- **Modified Gompertz:** $y = m \exp(-d \exp(-r t))$ (code `maxr`,
  `disp`, `rate`) with vertical offset adjustment.
- **Exponential decay:** `exp(-lambda_cor * tsd)` for all
  time-dependent effects (biomarker correlation, carryover).
- **Carryover on the mean:** anchored at the last on-drug mean,
  $\mu_{BR,\text{last}} 2^{-t_{sd}/t_{1/2}}$. The published code
  compounds the decay (`docs/53`); the drivers default to anchored
  carryover since 2026-10-08 (`--carryover recursive` reproduces
  stored runs).

**Biomarker parameterization.**

- **Centered biomarker:** `(bm - mean_bm) / sd_bm`, yielding a
  population-level moderation effect of zero.
- **BR scaling for architecture equivalence:** effect size
  calibrated to `c.bm * sigma_br / sigma_bm`.

## Data generation pipeline

**Four parameter sets** are required by `generateData`:

1. `modelparam`: sample size `N`, biomarker moderation `c.bm`.
2. `respparam`: Gompertz parameters (`maxr`, `rate`, `disp`) for
   each factor (br, tv, pb).
3. `blparam`: means and SDs for biomarker and nuisance factors.
4. `trialdesign`: output of `buildtrialdesign` (paths, timepoints,
   treatment phases).

Optional: `dgp_architecture`, `lambda_cor`, `carryover_t1half`,
`cached_sigma` (pre-built covariance for vectorized speed).

**Sigma matrix construction** encodes:

- AR(1) within-factor autocorrelation.
- Cross-factor correlations (BR-TV, BR-PB, TV-PB).
- Under covariance moderation only: biomarker-BR graded coupling,
  `Cor(BM, BR) = c.bm` on drug, `c.bm * exp(-lambda_cor * tsd)` off
  drug after first exposure, 0 before first exposure. (The published
  code uses a step rule that keeps full coupling while the mean is
  nonzero, which removes the interaction signal at any positive
  half-life; `docs/55`.)
- Bivariate structure: each timepoint replicates the 4-variate
  block (BM, BR, TV, PB).

**Coupling ceiling.** The largest feasible $c_{bm}$ is 0.256 under the
published compound-symmetry construct without carryover; 0.440 under
separable AR(1) at rho = 0.8; 0.474 (Hybrid 0.480) at rho = 0.7
(`docs/34`, `docs/35`, compendium).

Sigma caching across parallel workers is performed in
`generateSimulatedResults` to avoid redundant matrix construction.

**Covariance moderation (default).** Interaction encoded in the
covariance matrix via `dgp_architecture = 'mvn'`.

**Mean moderation.** With `dgp_architecture = 'mean_moderation'`,
after the MVN draw, BR is additively shifted at on-drug timepoints
by `c.bm * bm_z * br_sd` where `bm_z` is the standardized
biomarker.

## Analysis model

**Recommended analysis for the interaction (F1).**
`Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De` (bmc the biomarker
centered at its baseline mean; De the design expectancy $\eta_t$),
random intercept plus `corCAR1(form = ~t|ptID)`, CR2 cluster-robust
standard errors with Satterthwaite df (clubSandwich). It is robust to
a prognostic biomarker, which the plain model is not
(`docs/47`-`docs/49`, `docs/56`). Data-chosen half-life (provisional):
F1A0B, the half-life chosen by AIC without the interaction term and
the test calibrated by a parametric bootstrap repeating the selection;
F1 at a fixed $t_{1/2}^a = 0.2$ week is the validated fallback
(`docs/57` Section 8). Choosing the half-life by AIC and testing as
if it were fixed (F1A) is liberal.

**Historical analysis.** The published and earlier package analysis
was `Sx ~ bm + t + Dbc + bm:Dbc` with `nlme::lme`, residual
`corCAR1(form = ~t|ptID)`, random intercept only (`~1|ptID`).

`Dbc` (continuous drug indicator):

- On-drug: `Dbc = 1`.
- Off-drug after first exposure: `Dbc = exp(-lambda_cor * tsd)`.
- Before first exposure: `Dbc = 0` (pre-exposure guard; without it a
  placebo-first participant is coded exposed; paper 01 B.4).

**Main effect (paper 04).** See the working protocol
`analysis/scripts/treatment-main-effect/00-ademp-protocol.md` and
`docs/58`; main-effect analysis IDs carry an `M-` prefix.

## Simulation workflow

### Single trial

```r
td <- buildtrialdesign(
  n_reps = 1000,
  design = 'Hybrid',          # OL, CO, Hybrid, OL+BDC
  carryover_t1half = 1.0
)

dat <- generateData(
  modelparam = list(N = 50, c.bm = 0.45),
  respparam = list(...),
  blparam = list(...),
  trialdesign = td,
  dgp_architecture = 'mvn'
)

result <- lme_analysis(td, dat, op, full_model_out = TRUE)
beta_interaction <- result$summary['bm:Dbc', 'Value']
p_interaction    <- result$summary['bm:Dbc', 'p-value']
```

### Batch sweep

```r
results <- generateSimulatedResults(
  designparam     = list(carryover_t1half = c(0, 0.5, 1.0)),
  respparamsets   = expand.grid(...),
  blparamsets     = ...,
  modelparamsets  = expand.grid(N = c(30, 50, 70),
                                 c.bm = c(0.3, 0.45, 0.6)),
  n_reps          = 1000,
  n_cores         = -1,
  dgp_architecture = 'mvn',
  save_chunks     = TRUE
)
```

### Simulation server

Long runs go to the simulation host with `simrun` (see
`analysis/scripts/quick-sim/carryover-closed-form/RERUN-anchored.md`).
The host has 2 physical cores (4 threads); jobs run one at a time.
Results from the container and from a local session are not
draw-for-draw identical (different BLAS/LAPACK), so arms that must be
paired trial by trial have to run on the same platform.

## Testing

Framework: `tinytest`. Tests live in `inst/tinytest/` and are
invoked via `tests/tinytest.R`.

Tested components: `buildtrialdesign`, `utilities`
(modGompertz, cumulative effects), `generateData` (structure and
outcome computation), `lme_analysis` (model fit and coefficient
extraction), sigma construction (positive definiteness, structure),
`analyze_trial_extended` (threshold on the fitted biomarker scale).

Run tests:

```bash
make test
# or from R:
tinytest::test_all()
```

Cross-implementation parity tests
(`analysis/scripts/parity/test-parity-*.R`) validate numerical agreement
between `original-extended` and `tidyverse` under matching seeds.

## Docker workflow

Container-first development (zzcollab framework):

```bash
make r                 # enter container shell
make rstudio           # RStudio Server at localhost:8787
make check-renv        # strict renv validation
make check-system-deps # system dependency check
```

## Key technical decisions

**AR(1) vs. compound symmetry.** AR(1) within-factor temporal
correlation (`nlme::lme` with `corCAR1(form = ~t|ptID)`). Chosen
for clinical realism (nearby measurements more correlated) and
numerical tractability: the coupling ceiling rises from 0.256 under
the published compound-symmetry construct to 0.440 (rho = 0.8) or
0.474 (rho = 0.7) under the separable AR(1) form.

**Exponential carryover decay.** Biomarker-BR correlation and
drug indicator both decay as `exp(-lambda_cor * tsd)` during
off-drug periods. Pharmacokinetically plausible; tied to drug
half-life via `lambda_cor = ln(2) / t_half`.

**Two interaction architectures.** Support for both covariance
moderation (Hendrickson's formulation) and mean moderation (the
standard psychometric form); see "Key architectural feature" above
for their relation and the matched-N power comparison. See
`docs/legacy/19-mean-moderation-implementation-notes.md` for the
implementation history and `analysis/report/01-dgp-mean-moderation-
vs-mvn/` for the manuscript treatment.

**Continuous `Dbc` in analysis model.** Uses exponential-decay
continuous drug indicator rather than binary on/off. Matches the
probabilistic carryover structure of the DGP and improves power
estimation by leveraging off-drug information.

## Common workflows

### Reproduce Hendrickson Figure 4

```r
source('analysis/scripts/figure4/04-generate-power-head.R')
```

### Run the carryover-sensitivity study

```bash
Rscript analysis/scripts/carryover-sensitivity/01-run-factorial.R
Rscript analysis/scripts/carryover-sensitivity/04-run-sensitivity-blocks.R
Rscript analysis/scripts/carryover-sensitivity/02-summarise-grid.R
Rscript analysis/scripts/carryover-sensitivity/03-render-figures.R
Rscript analysis/scripts/carryover-sensitivity/04-sensitivity-figures.R
bash tools/render.sh analysis/report/02-carryover-sensitivity/report.Rmd
```

Data summaries are checkpointed to
`analysis/data/`; figures to
`analysis/figures/`.

### Render a manuscript or whitepaper

```bash
bash tools/render.sh analysis/report/<NN-slug>/report.Rmd
```

Always render through `tools/render.sh` (it stamps provenance and
stages a dated copy in `share/`); do not call `rmarkdown::render()`
directly. Render whitepapers from inside `docs/` so figure paths
resolve.

### Analyze real trial data

See vignette `vignettes/03-analyze-clinical-trial-data.Rmd`:

```r
result <- analyze_trial_extended(td, dat, op, threshold = 5)
print_trial_summary(result)
```

### Characterize carryover in new data

```r
cc <- characterize_carryover(
  td, dat,
  half_lives = c(0.25, 0.5, 1.0, 2.0, 4.0)
)
print_carryover_summary(cc)
# cc$best_t_half gives the AIC-preferred half-life
```

The AIC-preferred half-life is descriptive: testing at it as if it
had been fixed is liberal (`docs/57`, `docs/58`).

## Performance

- Use data.table vectorized operations (`:=`) over loops in
  `original` and `original-extended` collections.
- Pre-allocate matrices where possible.
- Leverage Cholesky factorization for MVN sampling efficiency.
- `furrr::future_map` for parameter grid sweeps; auto-detect
  cores with `n_cores = -1`.
- Sigma caching: build all unique matrices once and reuse across
  parallel workers.
- Progressive save: checkpoint replicate results every 500
  iterations for long runs (`save_chunks = TRUE`).
- Set `VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  OMP_NUM_THREADS=1` for parallel simulation runs, and keep the
  total worker count within memory (16 GB here).

## Debugging

**Positive definiteness failures.** Covariance matrices can fail
PD at extreme `c.bm` values (above the coupling ceiling).
Pre-flight validation in `generateSimulatedResults` flags PD
failures before the simulation loop; PD diagnostics at
`analysis/scripts/figure4/05-pd-diagnostics.R`.

**Type I error validation.** Run simulations under the null
(`c.bm = 0`) and check rejection rate at `alpha = 0.05` against a
Monte Carlo band (2 MCSE). The figure of 2.5-5.8% attributed to
Hendrickson et al. comes from `docs/03` and has not been verified
against the publication.

**Singularity in random effects.** If `lme` fails to converge,
check (a) within-subject variation in the trial design,
(b) residual degrees of freedom for nested designs, and (c)
whether a random slope is needed in addition to the random
intercept.

## Documentation conventions

- Academic, scholarly tone. No em dashes, no hyperbole. US English.
- The studies are exploratory: never describe them as
  pre-registered, pre-specified or confirmatory.
- Author all whitepapers and audit reports as 'pmsimstats team'.
- File naming: hyphen-separated lowercase.
- Markdown line wrapping: 78 characters.
- Notation per `analysis/report/NOTATION.md`; run
  `perl tools/notation-lint.pl` on edited documents.
- PDF rendering via `tools/render.sh`.
- LaTeX build artifacts and rendered output are gitignored;
  rendering produces fresh output from source on demand.

## Related documents

- `analysis/report/NOTATION.md` — canonical notation and vocabulary.
- `analysis/report/compendium-mathematics.tex` — mathematics
  compendium (the series' mathematical reference).
- `docs/59-whitepaper-refactoring-plan.md` — plan for the whitepaper
  series.
- `docs/52-whitepaper-review-2026-10-07.md` — last whitepaper review
  record.
- `docs/60-trials-submission-gap-analysis.md` — gap analysis and
  remediation plan for submitting papers 01 and 02 to *Trials*
  (2026-10-20).
- `docs/61-cross-factor-correlation-bug-paper02-driver.md` — the
  `c.cf1t`/`c.cfct` mismatch in paper 02's simulation driver
  (`simulation-core.R`), its exact size, and the fix (2026-10-21).
- `docs/62-carryover-specification-and-half-life-method-catalog.md`
  — catalog of paper 02's `G1`-`G9` specifications and the project's
  `F`-family half-life-selection methods (`docs/57`/`ch23`), what's
  determined about each and their relative cost (2026-10-21).
- `docs/63-why-correlation-fades-with-time-a-clinicians-guide.md` —
  clinician-facing explanation of the compound-symmetry-to-separable-
  AR(1) change (`docs/35`, `ch10`), in the style of `docs/53`/`docs/55`
  (2026-10-21).
- `docs/64-f1a0b-first-pass-evaluation.md` — first-pass evaluation of
  F1A0B against `docs/62`'s criteria: no Type I error inflation,
  including under a prognostic biomarker (prog-BL), but no power
  advantage over plain F1 (2026-10-22). Adds the F-family model to
  `simulation-core.R`; documents that paper 02's `CO` design preset
  is unidentifiable for this model (constant expectancy) and that
  `OLBDC` was substituted for the design-robustness arm.
- `docs/legacy/INDEX.md` — maps every old `docs/NN` number to its
  new home (archived files live in `docs/legacy/`).
- `docs/legacy/reorg-plan-2026-04-22.md` — the April 2026 repository
  reorganization plan (historical).
- `docs/legacy/23-documentation-compendium.md` — April 2026 documentation
  hub (historical; superseded by the documents above).
- `analysis/report/11-combined-dgp-architecture/planning/plan-architecture-c.md` — original implementation plan for
  the dual-channel architecture (formerly "Architecture C"), written
  when it was to be a third approach in paper 01. Superseded as to
  manuscript placement; see below.

## Open work items

### The dual-channel architecture -- paper 11

The dual-channel architecture (formerly "Architecture C") activates
both channels simultaneously with independent parameters `c.bm_a`
(mean channel) and `c.bm_b` (covariance channel), and recovers the
mean- and covariance-moderation architectures as its single-channel
limits. It was moved out of paper 01 into its own companion paper,
paper 11 (`analysis/report/11-combined-dgp-architecture/`). The move
is documented in
`analysis/report/01-dgp-mean-moderation-vs-mvn/architecture-c-extraction-plan.md`.

- **Paper 01** is restricted to the mean-versus-covariance contrast.
  As of 2026-09-30 it contains no reference to the dual-channel
  architecture or `@pmsimstats-paper11`. Do not reintroduce them.
- **Package**: the `dgp_architecture = 'combined'` branch and the
  `c.bm_b` lookup in `buildSigma` are implemented in
  `R/generateData.R`.
- **Driver and data**: `analysis/scripts/dgp-combined/` and
  `analysis/data/dgp-combined/`.
- **Outstanding for paper 11** (from the extraction plan; status
  not re-checked): re-run the boundary cells
  `(c.bm_a=0.45, c.bm_b=0)` and `(c.bm_a=0, c.bm_b=0.45)` under a
  driver whose parameters match the paper 01 production driver, so
  they reproduce the pure single-channel cases; derive the
  super-additivity result analytically; add MCSE columns and report
  all three $t_{1/2}$ panels.

### Paragraph structure -- all papers

The three-part `bullets / rgt / orig` paragraph scaffolding was
dropped on 2026-10-07. Manuscripts are plain prose: do not
reintroduce the `bullets`, `rgt` or `orig` environments, and do not
treat any prose as reserved for the author. No `report.Rmd` uses the
environments. The `claudecode.tex` files that define them, and the
`in_header: claudecode.tex` lines in some manuscripts, are inert
leftovers.

## Project lexicon and controlled vocabulary

A detailed glossary of N-of-1 trial and precision-medicine
terminology, with plain-language definitions and cross-field
synonyms (biostatistics, single-case psychometrics, precision
medicine), is maintained at
`docs/29-nof1-precision-medicine-lexicon.md` (241 terms across 13
thematic sections, plus a same-concept / different-name map). It
was mined from the project Zotero `nof1` collection. Consult it
when drafting manuscripts or naming code constructs so usage stays
consistent.

**Cross-field synonyms.** The same concept is often named
differently by field. Prefer the canonical term on the left in
project prose; recognize the aliases when reading sources.

- Patient-by-treatment interaction = subject-by-treatment
  interaction (psychometrics) = heterogeneity of treatment effects
  / HTE (precision medicine); its variance is $\tau^2$.
- Biomarker-treatment interaction ($\beta_{bm:D}$) is the project's
  core estimand.
- N-of-1 trial = single-patient / single-subject trial = single-
  case experimental design (SCED) = idiographic clinical trial.
- Predictive biomarker = effect modifier / moderator (psychology)
  = treatment-by-covariate interaction term (biostatistics).
- Random effects / variance components (biostatistics) = between-
  case variance (single-case) ; shrinkage / empirical Bayes / BLUP
  = borrowing strength / partial pooling.
- Serial correlation = autocorrelation = AR(1); the intraclass
  correlation (ICC) summarizes within-subject correlation only under
  compound symmetry.
- Aggregated N-of-1 analysis = mixed-model / IPD meta-analysis =
  multilevel meta-analysis of single cases.

**Core terms used throughout this codebase and manuscripts.**

- Design: crossover, period, cycle, sequence, washout (incl.
  active / analytical washout), run-in, open-label (OL), blinded
  discontinuation (BDC), Hybrid design, aggregated N-of-1.
- Carryover: simple vs complex carryover; behavioral /
  psychometric carryover; carryover half-life (`carryover_t1half`);
  exponential decay `exp(-lambda_cor * tsd)`; anchored vs
  compounding carryover; treatment-by-period interaction.
- Architecture: the mean-moderation and covariance-moderation
  architectures, and the dual-channel architecture (paper 11 only);
  step vs graded coupling; the biomarker-treatment interaction term
  (`bm:Dbc`).
- Analysis: linear mixed-effects model (LME / MMRM), `corCAR1`
  AR(1) residual correlation, random intercept, REML, CR2
  cluster-robust standard errors, expectancy terms, Type I error
  validation, statistical power, design matrix.
- Precision medicine: predictive vs prognostic biomarker;
  qualitative vs quantitative interaction; individual treatment
  effect (ITE); responder analysis; biomarker-stratified and
  enrichment designs.
- Reporting: CENT 2015 (CONSORT extension for N-of-1), PATH
  statement, RoBiNT scale, operating characteristics.

---

*Last updated: 2026-10-11*
