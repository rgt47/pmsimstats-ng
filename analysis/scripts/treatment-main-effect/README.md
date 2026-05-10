# treatment-main-effect: driver scripts

*2026-05-06 10:17 PDT*

Driver scripts for the manuscript at
`analysis/report/04-treatment-main-effect/`.

## Provenance

Migrated on 2026-05-06 from
`~/prj/res/06-nof1-power/nof1_power/analysis/scripts/`. On
2026-05-08 the simulation engine (the `nof1power` R package) was
vendored into `pmsimstats-ng` at `implementations/nof1power/`.
Drivers source the vendored package; the external compendium is
no longer required for re-execution.

## Directory contents

The migration preserved file names from the source compendium.
Files break into four loose groups.

### Top-level Rmd notebooks

| File | Role |
|---|---|
| `monte_carlo_design_comparison.Rmd` | Primary design-comparison notebook |
| `carryover_only.Rmd` | Carryover-isolation notebook |
| `carryover_sensitivity_analysis.Rmd` | Carryover sensitivity analysis |
| `carryover_sensitivity_analysis_fixed.Rmd` | Revised version of the above |
| `just_before_simulation.Rmd` | Pre-simulation diagnostics |
| `just_after_simulation.Rmd` | Post-simulation diagnostics |
| `pmsim_development.qmd` | Development-history Quarto notebook |

### Top-level R drivers

| File | Role |
|---|---|
| `simulation.R` | Main simulation orchestration |
| `sim.R`, `sim1.R` | Auxiliary simulation drivers |
| `vig_carry.R`, `vig_carryover.R`, `vig5.R` | Vignette-style drivers |

### Sensitivity factorial

`sensitivity/` holds the factorial sensitivity analysis with one
file per axis (effect size, half-life, decay form, tau2, sigma2,
cycles, patients, AR1, period, biomarker, misspecification, Morris
null). Run via `sensitivity/run-all.R`; plot via
`sensitivity/plot-all.R`.

### Planning and partitioning notes

| File | Role |
|---|---|
| `partitioning.md` | Notes on partitioning the simulation grid |
| `plan1.md` | Early planning document |
| `README_VIGNETTES.md` | Notes on the vignette-style drivers |
| `README_CARRYOVER.md` | Notes on the carryover analyses |

## Dependencies

All drivers assume `nof1power` is loaded:

```r
library(nof1power)
```

Install instructions are in
`analysis/report/04-treatment-main-effect/README.md`.

## Path adjustments after migration

The source compendium rooted paths at
`~/prj/res/06-nof1-power/nof1_power/`. After migration, `here()`
resolves to the `pmsimstats-ng` root. Each driver should be
reviewed for hard-coded paths to `analysis/data/`,
`analysis/figures/`, or `analysis/tables/` before being run, and
those paths reconciled with the `pmsimstats-ng` layout. The
vendored `nof1power` package source at
`implementations/nof1power/` does not include the original
compendium's `analysis/` tree; only the package skin
(DESCRIPTION, NAMESPACE, R/, man/, inst/, tests/, vignettes/)
was vendored.

## Author

pmsimstats team
