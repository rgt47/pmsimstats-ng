# Carryover sensitivity simulation

Orchestration for manuscript 02
(`manuscripts/02-carryover-sensitivity/`).

## Scope

Factorial simulation of biomarker-treatment interaction power across:

- **3 DGP carryover decay forms:** linear, exponential, Weibull
- **3 analysis-model carryover specifications:** binary-treatment
  (A1), exposure-weighted (A2), lagged-treatment (A3)
- **2 DGP architectures:** mean moderation (A), MVN differential
  correlation (B)
- **3 carryover half-lives:** 0, 0.5, 1.0 weeks
- **3 trial designs:** CO, N-of-1 (Hybrid), OL+BDC
- **2 sample sizes:** 35, 70
- **3 biomarker moderation strengths:** 0.0, 0.30, 0.45

Total cells: $3 \times 3 \times 2 \times 3 \times 3 \times 2 \times
3 = 972$. Replicates per cell: 50 (development), 500 (production).

## Scripts

| Script | Purpose |
|---|---|
| `simulation-core.R` | Shared helpers: design presets (OL, CO, Hybrid, OLBDC), multi-path data generation wrapper, long-form preparation with `Db`/`Dbc`/`L` columns, and the three analysis-spec fitters (A1/A2/A3). Sourced by `01-run-factorial.R`. |
| `01-run-factorial.R` | Primary simulation driver; writes `output/01-factorial.rds` |
| `02-summarise-grid.R` | Aggregate replicate-level results into power and type-I-error grids; writes `output/02-grid-summary.rds` and `manuscripts/data/02-grid-summary.rds` |
| `03-render-figures.R` | Generate PDF figures into `manuscripts/figures/` |

## Prerequisites

The linear and Weibull decay forms require the `carryover_decay()`
helper and `carryover_form` / `weibull_shape` plumbing in
`implementations/tidyverse/R/functions.R` (added in the same
working session that created this directory). The
`implementations/original-extended/` collection does **not** yet
support the alternative forms and is therefore not used here;
cross-validation against that collection is limited to the
`carryover_form = "exponential"` cells.

## Usage

```r
# Smoke test (2 reps x 4 designs, one parameter setting; ~2 seconds)
# Confirms the full pipeline runs. Writes output/01-smoke.rds.
Rscript analysis/02-carryover-sensitivity/01-run-factorial.R --smoke

# Development run (50 replicates per cell, minutes)
Rscript analysis/02-carryover-sensitivity/01-run-factorial.R --dev

# Production run (500 replicates per cell, hours)
Rscript analysis/02-carryover-sensitivity/01-run-factorial.R

# Summarise + write manuscript-side RDS
Rscript analysis/02-carryover-sensitivity/02-summarise-grid.R

# Figures
Rscript analysis/02-carryover-sensitivity/03-render-figures.R
```

## Pipeline status

- Full pipeline end-to-end validated on the smoke configuration:
  OL, CO, Hybrid, OLBDC designs generate data, prepare long-form
  with `Db`/`Dbc`/`L` columns, and fit all three analysis
  specifications (A1 binary, A2 exposure-weighted, A3
  binary + lagged) with non-NA p-values and estimates.
- OL designs correctly report `no off-drug observations` for A1
  and `no lagged-on timepoints` for A3 in the `reason` column;
  A2 on OL fails silently because `Dbc` is identically 1 (aliased
  with `bm`).
- The full factorial grid (162 cells x 500 reps = 81000 cell-
  replicate evaluations) has not been run. Production timing on
  the current hardware should be estimated from a dev run before
  committing to a full production run.

## Output conventions

- Raw replicate-level results: `output/01-factorial.rds` (large;
  gitignored).
- Summarised grid: `output/02-grid-summary.rds` and mirrored into
  `manuscripts/data/02-grid-summary.rds` (small; committed).
- Figures: `manuscripts/figures/02-*.pdf` (committed).

## Reproducibility

Each output RDS includes a `meta` list capturing R version,
`renv.lock` hash, `implementations/tidyverse/` commit SHA,
simulation seed, and parameter grid.
