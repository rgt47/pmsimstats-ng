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
# Development run (50 replicates, ~few minutes)
Rscript analysis/02-carryover-sensitivity/01-run-factorial.R --dev

# Production run (500 replicates, hours)
Rscript analysis/02-carryover-sensitivity/01-run-factorial.R

# Summarise + write manuscript-side RDS
Rscript analysis/02-carryover-sensitivity/02-summarise-grid.R

# Figures
Rscript analysis/02-carryover-sensitivity/03-render-figures.R
```

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
