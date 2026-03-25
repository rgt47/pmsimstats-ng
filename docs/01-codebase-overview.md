# pmsimstats Codebase Overview

The `pmsimstats` R package implements Monte Carlo
simulation-based power analysis for aggregated N-of-1
clinical trial designs, comparing four trial designs on their
ability to detect a biomarker-treatment interaction. Based on
Hendrickson et al. (2020), *Frontiers in Digital Health* 2:13.

## Core R Package Files (`R/`)

### generateData.R

The data-generating process (DGP). Two exported functions:

- **`generateData()`** -- Draws simulated participant data
  from a multivariate normal distribution. Accepts a
  pre-built covariance matrix (`cached_sigma`) or builds
  one via `buildSigma()`. Uses pre-computed Cholesky factor
  for efficient MVN draws (`L %*% Z + mu`). Computes
  outcome columns via vectorized `data.table` operations.

- **`buildSigma()`** -- Constructs the 26x26 covariance
  matrix for a single trial path. Implements:
  - AR(1) within-factor autocorrelation: `rho^|t_i - t_j|`
  - Cross-factor correlations (same-time and different-time
    with AR(1) decay)
  - Biomarker-BR correlation with exponential decay during
    off-drug periods: `c.bm * exp(-lambda_cor * tsd)`
  - Auto-derives `lambda_cor = ln(2) / carryover_t1half`
  - Positive definiteness check with verbose warning
  - Pre-computes Cholesky factor for caching

- **`validateParameterGrid()`** -- Tests all parameter
  combinations for PD before simulation. Reports failures
  upfront.

### generateSimulatedResults.R

The simulation loop. Single exported function:

- **`generateSimulatedResults()`** -- Orchestrates the full
  parameter sweep:
  1. PD pre-validation of all parameter combinations
  2. Sigma matrix caching (builds all unique matrices once)
  3. Parameter grid expansion across designs, response
     params, baseline params, and model params
  4. Per-parameter-set worker function
     (`run_one_paramset`) that generates data, fits LME,
     applies censoring, and collects results
  5. Parallel execution via `furrr::future_map()` with
     closure-based function passing to workers
  6. Progressive save support for long runs

  Key parameters: `lambda_cor` (correlation decay, default
  NA = auto), `n_cores` (parallelism, default auto-detect).

### lme_analysis.R

The analysis model. Single exported function:

- **`lme_analysis()`** -- Fits a linear mixed-effects model
  to simulated or actual trial data:
  - Uses `nlme::lme` with `corCAR1(form = ~t|ptID)` for
    continuous-time AR(1) residual correlation
  - Dynamic formula selection: `Sx ~ bm + t + Dbc +
    bm:Dbc` when within-subject drug variation exists;
    `Sx ~ bm + t + bm:t` for open-label
  - `Dbc` is a continuous drug indicator (1 when on drug,
    exponential decay when off drug)
  - Graceful handling of rank-deficient models (returns NA)
  - Fallback to `lme` without `corCAR1` if optimization
    fails
  - Extracts interaction coefficient from `nlme` tTable

### carryover_analysis.R

Extended analysis functions for actual trial data:

- **`characterize_carryover()`** -- Sweeps candidate
  carryover half-lives (default 0.1 to 8 weeks), fits
  `nlme::lme` at each, compares AIC/BIC to identify the
  best-fitting half-life. Returns drug effect, interaction,
  AR(1) phi at each half-life.

- **`analyze_trial_extended()`** -- Extracts both the main
  drug effect (`Dbc`) and biomarker interaction (`bm:Dbc`)
  with 95% CIs. Computes variance components, predicted
  drug effect as a function of biomarker value, and
  clinical decision threshold.

- **`print_carryover_summary()`** / **`print_trial_summary()`**
  -- Clinical reporting summaries.

### buildtrialdesign.R

Trial design specification:

- **`buildtrialdesign()`** -- Constructs trial paths from
  intuitive inputs (timepoints, expectancies, on-drug
  vectors). Computes derived timing variables: `tod`
  (time on drug), `tsd` (time since discontinuation),
  `tpb` (time in positive belief).

### utilities.R

Helper functions:

- **`modgompertz()`** -- Modified Gompertz function for
  response trajectories. Passes through origin, asymptotes
  at `maxr`.
- **`cumulative()`** -- Converts interval durations to
  cumulative time.
- **`reknitsimresults()`** -- Recombines progressive save
  chunks.

### censordata.R

Dropout simulation:

- **`censordata()`** -- Two-component dropout: time-dependent
  (flat rate) and outcome-dependent (biased by symptom
  change). Implements last-observation-carried-forward
  censoring.

### plottingfunctions.R

Visualization:

- **`PlotModelingResults()`** -- Power heatmaps across
  4-dimensional parameter grids with `ggplot2 geom_tile`
  and `facet_grid`.
- **`plotfactortrajectories()`** -- Factor trajectory plots
  by trial path.

## Analysis Scripts (`analysis/scripts/`)

### Simulation Scripts

| Script | Purpose |
|--------|---------|
| `01d_generate_head.R` | Figure 4: current code, parallel, c.bm={0,0.25,0.45} |
| `01b_generate_results_core_original.R` | Figure 4: original commit code |
| `01c_generate_by_commit.R` | Figure 4: parameterized by commit |
| `01_generate_results_core.R` | Figure 4: current code, sequential |
| `01e_pd_diagnostics.R` | PD failure rate instrumentation |
| `03_generate_figure5.R` | Figure 5: revised code |
| `03b_generate_figure5_original.R` | Figure 5: original code |
| `02_plot_figure4.R` | Heatmap plotting with provenance |

### Extracted Code by Commit

Code from each historical commit extracted for comparison:

```
analysis/scripts/commit_42ac030/   # Initial (publication)
analysis/scripts/commit_8609f12/   # Ron Thomas edits
analysis/scripts/commit_f6ee86d/   # Autocorrelation typo fix
analysis/scripts/commit_f70f86d/   # Carryover added
analysis/scripts/commit_325314f/   # Modified Db (broken)
analysis/scripts/original_R/      # Copy of 42ac030
```

## White Papers (`docs/pdf/`)

| Document | Pages | Content |
|----------|------:|---------|
| `04-revised-power-analysis.pdf` | 13 | Main paper: problems, fixes, validated results |
| `18-response-parameter-sensitivity.pdf` | 6 | Response parameter sensitivity before/after |
| `08-biomarker-correlation-decay.pdf` | 10 | Pharmacokinetic derivation of lambda_cor |
| `06-ar1-residual-correlation.pdf` | 8 | lmer to nlme transition with Type I error validation |
| `07-positive-definiteness-failures.pdf` | 10 | PD failure prevalence and causes |
| `09-carryover-correlation-artifact.pdf` | 10 | BM-BR correlation artifact analysis |
| `05-commit-impact-report.pdf` | 13 | Commit-by-commit Figure 4 comparison |

## Bundled Data (`data/`)

| File | Content |
|------|---------|
| `results_core.rda` | Published Figure 4 results (initial commit, 500 reps) |
| `extracted_bp.rda` | Baseline parameters (BP mean/SD, CAPS mean/SD) |
| `extracted_rp.rda` | Response parameters (Gompertz max/disp/rate/sd) |
| `CTdata.rda` | Actual clinical trial data |
| `results_maxes.rda` | Figure 5A results (response maxima) |
| `results_rates.rda` | Figure 5B results (response rates) |
| `results_trajectories.rda` | Factor trajectory data |

## Key Parameters

### Current Configuration (Revised Code)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Autocorrelation | AR(1), rho = 0.7 | PD-feasible, clinically realistic |
| c.bm | {0, 0.25, 0.45} | Max PD-feasible = 0.45 (OL+BDC) |
| Carryover t_half | {0, 0.5, 1.0} weeks | Clinically realistic |
| lambda_cor | ln(2) / t_half (auto) | Pharmacokinetic derivation |
| Cross-factor (same time) | c.cf1t = 0.2 | Unchanged from publication |
| Cross-factor (diff time) | c.cfct * rho^|t_i-t_j| | AR(1) decay applied |
| Analysis model | nlme::lme + corCAR1 | Required for AR(1) Type I error control |
| PD corrections | 0/162 (0%) | All matrices valid |

### Publication Configuration (Original Code)

| Parameter | Value | Issue |
|-----------|-------|-------|
| Autocorrelation | CS, rho = 0.8 | PD failures at c.bm > 0.25 |
| c.bm | {0, 0.3, 0.6} | 0.6 exceeds PD-feasible range |
| Carryover t_half | {0, 0.1, 0.2} weeks | Too short for genuine effects |
| BM-BR correlation | Step function | Power drop artifact |
| Analysis model | lmer (no corCAR1) | Correct under CS |
| PD corrections | 40/162 (24.7%) | Silent, uncontrolled |

## Dependencies

Core: `data.table`, `nlme`, `lme4`, `lmerTest`, `corpcor`,
`MASS`, `ggplot2`

Performance: `furrr`, `future` (parallel processing)

Legacy (in DESCRIPTION but used minimally): `svMisc`,
`tictoc`, `ggpubr`, `gridExtra`
