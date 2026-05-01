# pmsimstats Consolidation Plan

## Current State

Three repositories developing the same simulation framework
in parallel, with increasing overlap:

| Repo | Package | Role | Language | Status |
|------|---------|------|----------|--------|
| orig | pmsimstats v0.2.0 | Revised publication code | data.table | Active, 20 commits ahead |
| 2025 | pmsimstats2025 | Tidyverse rewrite | dplyr/tibble | Active, audit changes applied |
| simple | pmsimstats-simple | Pedagogical version | tidyverse | Stable, minor fixes applied |

### What Converged

After the audit, `orig` and `2025` now share:

- AR(1) autocorrelation (rho = 0.7)
- Exponential BM-BR correlation decay (lambda_cor = ln(2)/t_half)
- nlme::lme with corCAR1 analysis model
- No scale factor on carryover
- Timepoint-1 correlation fix
- Hendrickson-aligned Gompertz function
- Default parameters: c.bm <= 0.45, carryover 0-1.0 weeks

### What Diverges

| Feature | orig | 2025 | simple |
|---------|------|------|--------|
| Data manipulation | data.table | dplyr/tibble | tidyverse |
| Sigma caching | buildSigma() | build_sigma_matrix() | N/A |
| Parallel processing | furrr/future | furrr/future (infra, unused) | None |
| Carryover models | Exponential | Exponential + linear + Weibull | Exponential |
| Biomarker interaction | MVN correlation | MVN correlation | Direct mean moderation |
| Analysis model | nlme + corCAR1 | nlme + corCAR1 | lm on phase means |
| Trial design builder | buildtrialdesign() | buildtrialdesign() (copy) | Inline functions |
| Censoring | censordata() | censordata() (copy) | Simple MCAR |
| Plotting | PlotModelingResults() | Custom ggplot scripts | Custom ggplot |
| Docker/renv | None | Yes | Yes (zzcollab) |
| Tests | None | 14 unit tests | Placeholder |
| White papers | 7 documents | 20+ docs | simplification-plan.md |
| Commit history | v0.1.0 + v0.2.0 | 39 commits | Clean |

## Proposed Structure

One repository with three modules:

```
pmsimstats-ng/
|-- R/                        # Core package (shared)
|   |-- generateData.R        # DGP: buildSigma, generateData
|   |-- lme_analysis.R        # Analysis: nlme + corCAR1
|   |-- carryover_analysis.R  # Extended analysis functions
|   |-- buildtrialdesign.R    # Trial design specification
|   |-- censordata.R          # Dropout simulation
|   |-- utilities.R           # Gompertz, cumulative, helpers
|   |-- plottingfunctions.R   # Heatmap and trajectory plots
|   +-- generateSimulatedResults.R  # Simulation loop + parallel
|
|-- analysis/
|   |-- figure4/              # Figure 4 reproduction
|   |   |-- generate.R        # Simulation script
|   |   |-- plot.R            # Plotting script
|   |   +-- output/           # Results and heatmaps
|   |
|   |-- figure5/              # Figure 5 reproduction
|   |   |-- generate.R
|   |   |-- plot.R
|   |   +-- output/
|   |
|   |-- archive/              # Historical commit code
|   |   |-- commit_42ac030/   # v0.1.0 (publication)
|   |   |-- commit_8609f12/   # Ron Thomas edits
|   |   +-- ...
|   |
|   |-- carryover/            # Carryover characterization
|   |   |-- characterize.R
|   |   +-- output/
|   |
|   +-- gompertz/             # Gompertz fitting
|       |-- gompertz_fitting.Rmd
|       +-- output/
|
|-- simple/                   # Pedagogical module
|   |-- simulation.R          # Self-contained simplified sim
|   |-- docs/
|   |   +-- simplification-plan.md
|   +-- output/
|
|-- docs/                     # White papers and documentation
|   |-- pmsimstats_revision_summary.pdf
|   |-- figure5_comparison.pdf
|   |-- correlation_decay_derivation.pdf
|   |-- analysis_model_ar1.pdf
|   |-- pd_failure_white_paper.pdf
|   |-- hendrickson.pdf       # Published paper
|   |-- murray.pdf            # Raskind RCT paper
|   +-- CODEBASE-OVERVIEW.md
|
|-- data/                     # Bundled datasets
|   |-- CTdata.rda
|   |-- extracted_bp.rda
|   |-- extracted_rp.rda
|   +-- results_core.rda      # v0.1.0 publication results
|
|-- tests/
|   |-- testthat/
|   |   |-- test-gompertz.R
|   |   |-- test-sigma.R
|   |   |-- test-lme-analysis.R
|   |   |-- test-type1-error.R
|   |   +-- test-carryover.R
|   +-- testthat.R
|
|-- vignettes/                # Original vignettes (preserved)
|
|-- DESCRIPTION               # Single package
|-- NAMESPACE
|-- CLAUDE.md
|-- Makefile                  # Docker workflow
|-- Dockerfile
|-- renv.lock
+-- .github/workflows/       # CI/CD
```

## Migration Steps

### Phase 1: Establish the consolidated repo

1. Fork `orig` as the base (it has the audited code,
   version tags, and commit history).
2. Create `pmsimstats-ng` as the new repo name (drop the
   version suffixes).
3. Reorganize `analysis/` into topic subdirectories
   (figure4, figure5, carryover, gompertz, archive).
4. Move white papers from `analysis/output/` to `docs/`.

### Phase 2: Initialize zzc and import infrastructure

5. Run `zzc analysis` to initialize the workspace with
   Dockerfile, Makefile, .Rprofile, zzcollab.yaml,
   and renv scaffolding.
6. Import renv.lock from `2025` (has the right
   dependency set), then `renv::snapshot()` to add
   any missing packages (nlme, corpcor, furrr, etc.).
7. Import GitHub Actions CI from `simple`.
8. Import unit tests from `2025`
   (`test-hendrickson-dgp.R`) and expand.
9. Update DESCRIPTION Imports to match renv.lock.
   Remove unused dependencies.

### Phase 3: Import simple module

9. Copy `simple/simulation.R` as a self-contained module
   in `simple/`. It uses a different architecture (direct
   mean moderation, OLS) and serves a different purpose
   (pedagogy).
10. Keep `simple/` independent -- it does not share code
    with `R/`. Document the relationship in
    `simple/docs/simplification-plan.md`.

### Phase 4: Retire source repos

11. Archive `pmsimstats2025` (the audited functions are
    now in the consolidated `R/`).
12. Archive `pmsimstats-simple` (the simulation is now in
    `simple/`).
13. Keep `orig` as a read-only archive with the v0.1.0
    and v0.2.0 tags for reproducibility.

## What Gets Dropped

| Item | Reason |
|------|--------|
| 2025 `mod_gompertz()` | Does not pass through origin; use `modgompertz()` |
| 2025 `calculate_carryover_adjusted_correlations()` | Deprecated |
| 2025 30+ exploratory scripts | Archive only |
| 2025 17 unused DESCRIPTION dependencies | Never used |
| orig `comb.R` | Dead code (duplicate functions) |
| orig `scratch_WorkingScript.R` | Development scratch pad |
| Timestamped intermediate output files | Keep only final versions |

## What Gets Preserved

| Item | Location |
|------|----------|
| v0.1.0 publication code | `analysis/archive/commit_42ac030/` |
| All white papers | `docs/` |
| Original vignettes | `vignettes/` |
| Bundled .rda data | `data/` |
| Hendrickson-aligned functions | `R/utilities.R` (modgompertz) |
| 2025 unit tests | `tests/testthat/` |
| simple simulation | `simple/simulation.R` |
| Git history | Preserved via fork from orig |

## Naming Convention

- Package: `pmsimstats-ng` (no year suffix)
- Version: Semantic versioning from v0.2.0
- Tags: v0.1.0 (publication), v0.2.0 (revised DGP),
  v0.3.0 (consolidated)
- Functions: snake_case for new code, preserve existing
  camelCase for backward compatibility

## Estimated Effort

| Phase | Tasks | Effort |
|-------|------:|-------:|
| Phase 1: Reorganize orig | 4 | 1 day |
| Phase 2: Import from 2025 | 4 | 1 day |
| Phase 3: Import simple | 2 | 2 hours |
| Phase 4: Archive repos | 3 | 1 hour |
| Testing and validation | -- | 1 day |
| **Total** | | **~3 days** |

## zzcollab Workspace

The consolidated repo will be a zzcollab (zzc) workspace,
providing the full reproducible research environment:

```
pmsimstats-ng/
|-- zzcollab.yaml             # zzc project config
|-- .zzcollab/
|   +-- manifest.json
|-- Dockerfile                # From zzc analysis profile
|-- Makefile                  # zzc targets (make r, make rstudio, etc.)
|-- .Rprofile                 # renv + container detection
|-- renv.lock                 # Dependency snapshot
+-- ...
```

**Setup:** `zzc analysis` to initialize, then `make r` to
enter the container. All dependencies (nlme, data.table,
furrr, corpcor, MASS, ggplot2) locked via renv.

**Key zzc targets:**

- `make r` -- Enter container terminal
- `make rstudio` -- Start RStudio Server
- `make check-renv` -- Validate dependencies
- `make docker-build` -- Build from renv.lock
- `make docker-rebuild` -- Force fresh build

**Profile:** `analysis` (tidyverse + nlme + simulation
dependencies). The zzcollab.yaml:

```yaml
docker:
  default_profile: "analysis"
```

## Decision Points

1. **data.table vs tidyverse:** The core `R/` functions use
   data.table (from orig). The 2025 code uses tidyverse.
   Recommendation: keep data.table for the core (it is
   faster and the code works), add tidyverse as a Suggests
   for analysis scripts.

2. **Package vs project:** Should this be an installable R
   package or a research compendium? Recommendation:
   maintain as a package (with DESCRIPTION, NAMESPACE,
   exports) for `devtools::load_all()` workflow, but do
   not submit to CRAN.

3. **Docker base image:** Use `rocker/tidyverse:4.5` (from
   zzc analysis profile) with renv for dependency management.

4. **Where to host:** GitHub under the same organization
   as the original. The original `rchendrickson/pmsimstats`
   URL could point to the consolidated repo.

5. **zzc profile:** `analysis` provides tidyverse, renv,
   and the base R environment. Additional packages (nlme,
   corpcor, furrr, data.table, MASS) are installed via
   renv.lock.
