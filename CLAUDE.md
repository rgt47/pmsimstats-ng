# pmsimstats-ng Project Guidelines

This document captures project-specific conventions, architectural patterns, and workflow guidance for the pmsimstats-ng R package simulation framework.

## Project Overview

**Purpose:** Monte Carlo simulation-based power analysis for aggregated N-of-1 clinical trial designs, with focus on detecting biomarker-treatment interactions under carryover effects.

**Key architectural feature:** Dual DGP (data-generating process) architectures:
- **Architecture B (MVN):** Biomarker-response interaction via differential correlation in covariance structure
- **Architecture A (Mean Moderation):** Interaction via direct mean scaling by biomarker value

Both architectures produce qualitatively different power losses under carryover (40-60% vs. 8-13%).

## Code Organization

### R/ Package Directory

**Core data pipeline:**
1. `buildtrialdesign.R` - Trial design specification (paths, timepoints, treatment phases)
2. `generateData.R` - MVN data generation with dual DGP architecture support
3. `censordata.R` - Missing data mechanisms (dropout patterns)
4. `lme_analysis.R` - LME analysis model with corCAR1 residual structure
5. `generateSimulatedResults.R` - Simulation orchestration and parallel execution
6. `carryover_analysis.R` - Extended analysis (characterize half-lives, clinical thresholds)
7. `plottingfunctions.R` - Visualization (heatmaps, trajectory plots)

**Utilities:**
- `utilities.R` - Cumulative effects, modGompertz function
- `packagedocumentation.R`, `datadocumentation.R` - Package-level documentation

### analysis/ Directory

Organized by research question:
- `figure4/`, `figure5/` - Publication figure reproduction and sensitivity analysis
- `architecture_comparison/` - DGP architecture contrast studies
- `carryover_factorial/` - Carryover effect decomposition
- `2025/` - Alignment work between original and tidyverse implementations
- `archive/` - Historical versions and snapshots

### docs/ Documentation

**Strategic documentation** (43 Markdown and LaTeX files):
- `00-documentation-index.tex` - Master index and navigation
- `01-codebase-overview.md` - Architecture and function reference
- `02-dgp-mean-moderation-vs-mvn.md` - Dual architecture white paper
- `03-audit-and-revision-report.md` - Code audit findings and corrections
- `04-revised-power-analysis.tex` - Publication revision summary
- `06-09` - Technical deep dives (AR(1), PD failures, correlation decay)
- `10-20` - Implementation guides and case studies
- Individual PDFs generated from .md and .tex for distribution

All documents authored as "pmsimstats team" (not individual names).

## Code Conventions

### R Naming

- **Functions:** camelCase (`generateData`, `lme_analysis`, `buildtrialdesign`)
- **Variables:** snake_case (`c.bm`, `carryover_t1half`, `t_random_slope`) or descriptive abbreviations:
  - `tod` = time on drug (cumulative, Gompertz input)
  - `tsd` = time since discontinuation (carryover input)
  - `tpb` = time × placebo belief accumulation
  - `Dbc` = continuous drug indicator (analysis model)
- **Data objects:** lowercase underscore (`trial_design`, `resp_params`, `bl_params`, `model_params`)

### Documentation Style

All exported functions use roxygen2 documentation:
```r
#' Function title
#'
#' Description (separate paragraph on purpose)
#'
#' @param name Type and meaning
#' @param name Type and meaning
#'
#' @returns Type and structure of return value
#'
#' @examples
#' # Show usage with realistic parameters
#'
#' @export
```

Examples in docstrings should use vignette-style complete workflows, not minimal toy examples.

### Mathematical/Statistical Patterns

**Parameterization conventions:**
- **AR(1) correlation:** `rho^|t_i - t_j|` (uses abs difference of timepoint indices)
- **Carryover half-life:** `carryover_t1half` (weeks, default 0.5-1.0 for PTSD/prazosin context)
- **Lambda decay:** `lambda_cor = ln(2) / carryover_t1half` (auto-derived unless overridden for sensitivity)
- **Modified Gompertz:** `y = maxr * exp(-disp * exp(-rate * t))` with vertical offset adjustment
- **Exponential decay:** `exp(-lambda_cor * tsd)` for all time-dependent effects (biomarker correlation, carryover)

**Biomarker parameterization:**
- **c.bm:** Under MVN architecture, correlation between biomarker and BR. Under mean moderation, regression coefficient scaling BR by centered biomarker.
- **Centered biomarker:** `(bm - mean_bm) / sd_bm` for population-level effect = 0
- **BR scaling:** Effect size calibrated to `c.bm * sigma_br / sigma_bm` for architecture equivalence

## Data Generation Pipeline

### Parameter Structure

**Four parameter types** (all required for generateData):
1. `modelparam`: N (sample size), c.bm (biomarker moderation/correlation)
2. `respparam`: Gompertz response parameters (maxr, rate, disp for each factor: br, tv, pb)
3. `blparam`: Biomarker and nuisance factor means and SDs
4. `trialdesign`: Paths, timepoints, treatment phases (output of buildtrialdesign)

**Optional parameters:**
- `dgp_architecture` ("mvn" default, or "mean_moderation")
- `lambda_cor` (auto-derived from carryover_t1half if NA)
- `carryover_t1half` (0 for no carryover; realistically 0.5-1.0 weeks)
- `cached_sigma` (pre-built covariance for vectorized speed)

### Sigma Matrix Construction

The covariance matrix encodes:
- AR(1) within-factor autocorrelation
- Cross-factor correlations (BR-TV, BR-PB, TV-PB)
- Biomarker-BR differential correlation (Architecture B only):
  - On-drug: `Cor(BM, BR) = c.bm`
  - Off-drug: `Cor(BM, BR) = c.bm * exp(-lambda_cor * tsd)`
- Bivariate structure: each timepoint replicates the 4-variate structure (BM, BR, TV, PB)

**Performance optimization:** Sigma caching in generateSimulatedResults avoids redundant matrix construction.

### DGP Architectures (Switchable)

**Architecture B (MVN, default):**
```r
dgp_architecture = "mvn"
# Interaction encoded as differential correlation in sigma matrix
# Off-drug correlation decays: reduces power (40-60% loss under carryover)
```

**Architecture A (Mean Moderation):**
```r
dgp_architecture = "mean_moderation"
# After MVN draw, additively shift BR at on-drug timepoints:
# dat[[br_col]] += c.bm * bm_z * br_sd  (where bm_z = standardized biomarker)
# Proportional relationship preserved under carryover (8-13% loss)
```

### Analysis Model

**Fixed formula:**
```
Sx ~ bm + t + Dbc + bm:Dbc
```

**Dbc (continuous drug indicator):**
- On-drug: 1
- Off-drug: exponential decay `exp(-lambda_cor * tsd)`
- Carryover half-life parameterized in analysis options (can differ from DGP)

**Error structure:**
- `nlme::lme` with `corCAR1(form = ~t|ptID)` (AR(1) continuous-time autocorrelation)
- Random intercept only (`~1|ptID`)

## Simulation Workflow

### Single Trial Simulation

```r
# Define trial design
td <- buildtrialdesign(
  n_reps = 1000,        # observations per participant
  design = "Hybrid",    # OL, CO, Hybrid, or OL+BDC
  carryover_t1half = 1.0
)

# Generate data
dat <- generateData(
  modelparam = list(N = 50, c.bm = 0.45),
  respparam = list(...),  # Gompertz parameters
  blparam = list(...),    # biomarker/nuisance distributions
  trialdesign = td,
  dgp_architecture = "mvn"
)

# Analyze
result <- lme_analysis(td, dat, op, full_model_out = TRUE)

# Extract
beta_interaction <- result$summary["bm:Dbc", "Value"]
p_interaction <- result$summary["bm:Dbc", "p-value"]
```

### Batch Simulation

```r
# Parameter grid expansion
results <- generateSimulatedResults(
  designparam = list(carryover_t1half = c(0, 0.5, 1.0)),
  respparamsets = expand.grid(...),  # Gompertz sensitivity grid
  blparamsets = ...,
  modelparamsets = expand.grid(N = c(30, 50, 70), c.bm = c(0.3, 0.45, 0.6)),
  n_reps = 1000,
  n_cores = -1,  # Auto-detect available cores
  dgp_architecture = "mvn",
  save_chunks = TRUE  # Progressive saving for long runs
)
```

## Testing Framework

**Framework:** testthat (3rd edition)

**Test organization:**
- `tests/testthat.R` - Test configuration
- `tests/testthat/test-*.R` - Individual test files

**Tested components:**
- `buildtrialdesign.R` - Design construction, timing variables
- `utilities.R` - Cumulative effects, Gompertz curve
- `generateData.R` - Data structure, outcome computation
- `lme_analysis.R` - Model fitting, coefficient extraction
- `sigma.R` - Positive definiteness, covariance structure

**Test discipline:**
- Input-output validation (what you put in, what comes out)
- Functional correctness (does the math work?)
- Avoid testing implementation details
- Vignettes are live documentation; tests validate core functionality

**Run tests:**
```bash
make test          # or
devtools::test()   # in R
```

## Docker Workflow

**Container first:** All development and testing via Docker (zzcollab framework).

```bash
# Enter container
make r
# Equivalent to: docker-compose run --rm r bash

# Run RStudio Server
make rstudio
# Access at http://localhost:8787

# Validate environment
make check-renv         # Strict mode: all packages match lock file
make check-system-deps  # System dependency check
```

## Documentation Standards

### Markdown Files (.md)

- **Line wrapping:** 78 characters max (code and comments)
- **Lists:** Blank line before first item
- **Author block:** "pmsimstats team" (not individual names)
- **Rendering:** Pandoc with xelatex, DejaVu Sans Mono for code
- **Timestamp footer:** "Rendered on YYYY-MM-DD at HH:MM TZ. Source: ~/path/to/file.md"

### LaTeX Files (.tex)

- Use `pmsimstats-preamble.tex` for common formatting
- Author block in preamble section (changed to "pmsimstats team")
- Date field should match document content date (e.g., last major revision)

### White Papers and Audit Reports

These documents establish the scientific and technical record:
- 02-dgp-mean-moderation-vs-mvn.md: Formalizes dual architecture distinction
- 03-audit-and-revision-report.md: Lists seven specific DGP corrections with justification
- 04-revised-power-analysis.tex: Before/after power estimates
- 08-biomarker-correlation-decay.tex: Mathematical derivation of lambda_cor rule
- Others: Deep dives on specific technical issues

**Key principle:** Every claim has a source (code commit, derivation, or citation).

## Key Technical Decisions

### 1. AR(1) vs. Compound Symmetry

**Decision:** AR(1) autocorrelation for within-factor temporal correlation.

**Rationale:** More clinically realistic (nearby measurements more correlated) and numerically beneficial (expands feasible parameter space for c.bm from max 0.25 to 0.49+).

**Implementation:** nlme::lme with corCAR1(form = ~t|ptID).

### 2. Exponential Carryover Decay

**Decision:** Biomarker-BR correlation decays as `c.bm * exp(-lambda_cor * tsd)` during off-drug periods.

**Rationale:** Pharmacokinetically plausible; tied to drug half-life via lambda_cor = ln(2) / t_half.

**Sensitivity:** carryover_t1half is a free parameter (0, 0.5, 1.0 weeks for PTSD/prazosin).

### 3. Dual DGP Architectures

**Decision:** Support both MVN differential correlation (Architecture B, original) and direct mean moderation (Architecture A, standard in literature).

**Rationale:** Discovered they produce 40-60% divergent power losses under carryover; both are scientifically defensible depending on biomarker biology.

**Impact:** switchable via dgp_architecture parameter; enables sensitivity analysis and architectural clarity.

### 4. Continuous Dbc in Analysis Model

**Decision:** Use exponential-decay continuous drug indicator in analysis, not binary on/off.

**Rationale:** Matches the probabilistic carryover structure of the DGP; improves power estimation by leveraging off-drug information.

**Alternative considered:** Separate analysis with carryover as a nuisance parameter (less efficient).

## Common Workflows

### Reproduce Publication Figure

```r
# Figure 4 (power heatmaps under carryover)
source("analysis/figure4/04-generate-power-head.R")

# Figure 5 (response parameter sensitivity)
source("analysis/figure5/01-generate-parameter-sensitivity.R")
```

### Compare DGP Architectures

```r
# See analysis/architecture_comparison/01-compare-dgp-architectures.R
# Runs dual simulations, generates comparative heatmaps and RDS outputs
```

### Analyze Real Trial Data

```r
# See vignette 03-analyze-clinical-trial-data.Rmd
# Load trial data, apply lme_analysis(), extract interaction coefficient
result <- analyze_trial_extended(td, dat, op, threshold = 5)
print_trial_summary(result)
```

### Characterize Carryover in New Data

```r
cc <- characterize_carryover(
  td, dat,
  half_lives = c(0.25, 0.5, 1.0, 2.0, 4.0)
)
print_carryover_summary(cc)
# cc$best_t_half gives optimal half-life via AIC comparison
# cc$best_model gives full fitted model
```

## Performance Considerations

### Vectorization

- Use data.table vectorized operations (`:=`) over loops
- Pre-allocate matrices when possible
- Leverage Cholesky factorization for MVN sampling efficiency

### Parallelization

- furrr::future_map() for parameter grid sweeps
- Auto-detect cores with n_cores = -1
- Sigma caching: build all unique matrices once, reuse in parallel workers
- Progressive save: checkpoint results every 500 replicates

### Memory

- Large parameter grids (>100K parameter sets): Enable save_chunks = TRUE
- Monitor RAM with complex designs (>1000 timepoints per participant)

## Debugging and Diagnostics

### Positive Definiteness Failures

Covariance matrices can fail positive definiteness at extreme c.bm values. Diagnostics:
```r
# Pre-flight validation (in generateSimulatedResults)
validateParameterGrid()  # Flags PD failures before simulation

# Empirical diagnostics
analysis/figure4/05-pd-diagnostics.R  # Compute eigenvalues, condition numbers
```

### Type I Error Validation

Run simulations under null (c.bm = 0) to check inflation:
```r
# Should see ~5% rejection at alpha=0.05
# Hendrickson et al. reported 2.5-5.8% (nominal) across all designs
```

### Singularity Issues

If random effects model fails to converge:
- Check trial design (must have within-subject variation)
- Verify residual degrees of freedom (nested design)
- Consider random slope models if random intercept alone insufficient

## Tone and Communication

**Documentation tone:** Academic and scholarly, third person, no em dashes or hyperbole.

**Examples:**
- ✓ "The revised code implements exponential decay tied to pharmacokinetic half-life."
- ✗ "The original code had a nasty step-function bug that needed fixing."

**File naming:** Use hyphens, not underscores (02-dgp-mean-moderation-vs-mvn.md, not 02_dgp_...).

**Commit messages:** Link to specific issues; reference code lines when fixing bugs.

---

*Last updated: 2026-04-08*
