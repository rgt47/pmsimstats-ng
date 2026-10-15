# Software and parity workflow

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** One-page pointer to the parity workflow and software
documentation.
**Prerequisites:** None

---

## Contents

1. [Notation](#notation)
2. [Summary](#summary)
3. [Construct box](#construct-box)
4. [Setting and argument](#setting-and-argument)
5. [Derivations](#derivations)
6. [Numerical evidence](#numerical-evidence)
7. [Relation to the manuscripts](#relation-to-the-manuscripts)
8. [Limitations and open items](#limitations-and-open-items)
9. [Reproduce](#reproduce)
10. [References](#references)
11. [Change log](#change-log)

---

## Notation

*Local symbols used in this chapter only. For project-wide notation
see `analysis/report/NOTATION.md`.*

Code identifiers appear in backticks only (`c_bm`, `t_half`,
`n_cores`, `Dbc`).

---

## Summary

This appendix has two parts:

1. **Parity workflow** -- a verbatim copy of
   `analysis/scripts/parity/README.md` (2026-04-15, with a statement
   check against source 2026-10-07). The parity workflow is the
   safety net for tidyverse-side development; this document describes
   how to integrate it into day-to-day development rather than treat
   it as an after-the-fact audit.
2. **Software pointer** -- a short description of the five
   implementation collections under `implementations/`.

---

## Construct box

*No Monte Carlo simulation results in this chapter.*

---

## Setting and argument

<!-- From analysis/scripts/parity/README.md (2026-04-15) -->

### Part I: Parity workflow for tidyverse-side development

#### Notation

- **Parity.** Numerical agreement between `original-extended`
  (data.table) and `tidyverse` under matching seeds, to the script's
  tolerances: $10^{-10}$ on the covariance matrix and the generated
  data, $10^{-6}$ on the interaction estimate and its standard error,
  $10^{-8}$ on its p-value.
- **Cell.** One combination of design, data-generating architecture
  (`mvn` or `mean_moderation`), carryover half-life, `c_bm`, and seed
  in the parity grid.
- **Layer.** One of the compared stages: covariance matrix (Sigma),
  generated data, analysis, and a batch smoke test.
- **Baseline.** A saved parity report (`parity_baseline.rds`) against
  which later reports are compared.
- **Expected divergence.** A cell-layer failure annotated as
  intentional in `analysis/scripts/parity/parity_divergences.csv`.

#### 1. Baseline before you start

Run the full parity sweep once on the clean branch and save the
report as the frozen "this is what parity looks like before I touched
anything" artifact:

```bash
git checkout -b tidyverse-refactor
Rscript analysis/scripts/parity/test-parity-extended-tidyverse.R
cp analysis/scripts/parity/parity_report.rds \
   analysis/scripts/parity/parity_baseline.rds
```

Commit `parity_baseline.rds` with the branch. At the end of the
development cycle, diff the current `parity_report.rds` against it to
confirm the only differences (if any) are the ones you intended.

#### 2. Pick a tight inner loop

The default mode runs in approximately 10 seconds on current hardware.
Use `--quick` (approximately 2 seconds, 16 cells) while iterating on
a specific function and escalate to the full sweep before committing.
The default mode and `--full` run the same grid: 144 cells (four
designs, two architectures, three half-lives, three values of `c_bm`,
two seeds); `--full` differs only in the mode label stored in the
report (inspected, 2026-10-07).

An optional git pre-commit hook can enforce this:

```bash
Rscript analysis/scripts/parity/test-parity-extended-tidyverse.R \
  --quick || \
  { echo 'parity broken; run the full sweep for details'; exit 1; }
```

No such hook is installed in the repository as of 2026-10-07
(`.git/hooks/` holds only the sample hooks; inspected).

#### 3. Edit semantics versus edit shape

Parity is preserved if a change is "same output, different code." It
is broken if the change is "different output, possibly better." Be
deliberate about which is intended:

- **Pure refactors** (rename, split function, swap `dplyr` for
  `purrr`, extract helper): parity must hold. If it breaks, a bug has
  been introduced.
- **Feature additions** that extend behavior without changing the
  existing path (new optional argument with default value, new
  analysis option): parity must still hold for the default-argument
  call paths the test exercises.
- **Behavior changes** (different DGP parameterization, different
  MVN-draw ordering for speed, different Sigma caching scheme): parity
  will break by design. Before making such a change, decide whether
  this also means breaking parity with `original-extended` (in which
  case the parity test becomes the wrong oracle) or whether the change
  should be applied to both implementations in lockstep.

#### 4. When parity breaks mid-refactor

The test reports which layer broke. Diagnose in this order:

1. **Sigma failure:** the change likely touched `build_sigma_matrix`
   or `build_correlation_matrix`. The BM-BR correlation branch is the
   most common culprit.
2. **Data failure with Sigma passing:** the change touched the MVN
   sampling path, the post-draw mean-moderation shift, or the
   participant/path assembly logic.
3. **Analysis failure with data passing:** the change touched formula
   construction, `Dbc` computation, or the `nlme::lme` invocation in
   `lme_analysis`.

The test prints a pass count per layer and, on failure, the table of
unexpected failing cells for each layer, then exits with status 1.

#### 5. Scope edits to one concern at a time

The `--quick` run covers 16 cells in approximately 2 seconds.
Re-running after every approximately 20 lines of changes localizes
bugs to the edit that introduced them.

#### 6. Guardrails around parallelization and RNG

The parity test calls `set.seed()` at the per-cell level. If a
refactor consumes random draws in a different order (for example,
changing `purrr::map` to `furrr::future_map` or reordering per-path
MVN calls), parity will break even though the code is semantically
correct. Two defenses apply:

- Do not change RNG-consuming loops unless `original-extended` is
  updated in the same commit.
- Keep `n_cores = 1` for any parity-relevant comparison. Introduce
  parallelism only after functional parity is locked in.

#### 7. Expand the parity test as coverage gaps appear

When a feature is added that the current grid does not exercise, add a
row to the `edge_cells` tibble in the script's edge-cell section
(currently five cells: zero `c_bm` with carryover under each
architecture, two small-N cells, one empirical draw; inspected
2026-10-07). A one-line edit per new feature is inexpensive, and the
habit compounds.

#### 8. When divergence is intentional

If the intent is to change tidyverse behavior in a way that
`original-extended` should not follow:

- Record which cells of the parity test are allowed to fail and why:
  add a row (design, architecture, half-life, `c_bm`, seed, layer,
  reason) to `analysis/scripts/parity/parity_divergences.csv`, and
  say so in the commit message.
- The script reads that file and counts matching failures as expected
  divergences, so the summary reports them separately and they do not
  set a failing exit status. (Implemented as the CSV annotation file;
  inspected 2026-10-07.)
- Re-baseline after the intentional divergence lands:
  `cp analysis/scripts/parity/parity_report.rds analysis/scripts/parity/parity_baseline.rds`.

#### 9. Before merging

Run the full sweep (not `--quick`) one final time. If the diff between
`parity_report.rds` and `parity_baseline.rds` is empty (or matches
documented expected divergences), the branch is safe to merge.

The diff script does this and exits 0 when the reports agree:

```bash
Rscript analysis/scripts/parity/parity-diff.R
```

A one-liner for the same comparison of the main grid:

```r
a <- readRDS('analysis/scripts/parity/parity_baseline.rds')
b <- readRDS('analysis/scripts/parity/parity_report.rds')
all.equal(a$grid, b$grid)  # TRUE, or a character vector of diffs
```

Compare reports produced in the same mode: a `--quick` report has 16
grid rows and a default or `--full` report 144.

#### Checklist

- [ ] Copy `parity_report.rds` to `parity_baseline.rds` at the start
      of the branch.
- [ ] `--quick` after every focused edit (~2 s).
- [ ] Full sweep before every commit (~10 s).
- [ ] If parity breaks, read the failing-cell table; fix before
      proceeding.
- [ ] Never change RNG consumption order in a pure refactor.
- [ ] Add edge cells to the script when new features are added.
- [ ] Record any intentional divergence in `parity_divergences.csv`
      and the commit message.

#### Optional infrastructure status (2026-10-07, inspected)

1. Pre-commit hook that runs `--quick`: not installed.
2. Expected-divergence annotation: implemented as
   `parity_divergences.csv`.
3. Baseline-diff wrapper: implemented as
   `analysis/scripts/parity/parity-diff.R`.

---

### Part II: Software collections

<!-- New to write: one-page pointer to implementations/ -->

Five collections live under `implementations/`. Four are
biomarker-treatment-interaction simulators (one historical, two
production, one reduced-complexity sandbox); the fifth is a vendored
R package supporting paper 04's treatment-main-effect comparison.

| Collection | Style | Architecture | Role |
|---|---|---|---|
| `original` | data.table | covariance moderation only | Historical reference (Hendrickson commit 58b32a9); back-compat testing |
| `original-extended` | data.table | both | Production data.table collection; oracle for parity testing |
| `tidyverse` | tidyverse | both | Production tidyverse collection; subject of the parity workflow |
| `simple` | base R | mean moderation only | Pedagogical sandbox; reduced-complexity |
| `nof1power` | R package | mean moderation only | Vendored simulation engine for paper 04; install with `devtools::install_local('implementations/nof1power')` |

Table: Implementation collections under `implementations/`

Common parameters accepted by all collections that support them:

- `c.bm`: under covariance moderation, the biomarker-BR correlation
  $c_{bm}$; under mean moderation, the coefficient scaling BR by the
  centered biomarker $\beta_{bm}$.
- `carryover_t1half`: half-life of the carryover effect (weeks).
- `lambda_cor`: correlation decay rate, auto-derived from
  `carryover_t1half` if `NA`.
- `dgp_architecture`: `'mvn'` (covariance moderation) or
  `'mean_moderation'`. Absent in `original` (fixed to MVN); present
  in the others with default `'mvn'`.

Parity artifacts in `analysis/scripts/parity/`:

- `parity_baseline.rds`: frozen baseline report.
- `parity_report.rds`: current sweep report.
- `parity_divergences.csv`: annotated expected divergences.
- `test-parity-extended-tidyverse.R`: parity test script.
- `parity-diff.R`: baseline-diff wrapper.

---

## Derivations

<!-- one row per result: statement | CM label | verification
script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

The parity workflow ensures that results reported from the tidyverse
collection agree with `original-extended` under matching seeds. Any
manuscript result that uses the tidyverse collection should be
verified against the parity baseline before submission.

---

## Limitations and open items

- The run times quoted in the parity README (approximately 2 s quick,
  approximately 10 s full) are from 2026-04-15 and have not been
  re-measured on current hardware.
- The five-cell edge-cell set does not cover the `dgp_architecture =
  'combined'` path (dual-channel, paper 11) or the `nof1power`
  collection.
- The pre-commit hook is not installed; the checklist must be followed
  manually.

---

## Reproduce

Parity test scripts are in `analysis/scripts/parity/`. Source:

```bash
cat analysis/scripts/parity/README.md
```

```bash
Rscript analysis/scripts/parity/test-parity-extended-tidyverse.R
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Verbatim parity workflow from `analysis/scripts/parity/README.md`; software collections pointer from `pmsimstats-ng/CLAUDE.md` | Phase 2 content move |
