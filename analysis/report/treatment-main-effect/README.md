# treatment-main-effect

*2026-05-06 10:17 PDT*

## Working title

Power comparison of N-of-1 and parallel-group designs for detecting
the treatment main effect in the presence of carryover.

## Provenance

This manuscript was migrated into `pmsimstats-ng` on 2026-05-06 from
the standalone research compendium

```
~/prj/res/06-nof1-power/nof1_power/
```

The simulation engine that produced the results in this manuscript
remains in that compendium and is consumed here as an installed R
package, `nof1power`. The manuscript and its driver scripts live
inside `pmsimstats-ng` so that the full set of papers in development
share a common build, bibliography, and review surface.

## Directory contents

| File | Role |
|---|---|
| `report.Rmd` | Canonical manuscript source |
| `report_short.Rmd` | Short-format companion |
| `report.tex`, `report.pdf` | Last rendered LaTeX and PDF |
| `report_short.tex`, `report_short.pdf` | Short-format outputs |
| `references.bib` | Manuscript-specific bibliography |
| `nof1-pgt.bib` | Auxiliary bibliography (PGT) |
| `statistics-in-medicine.csl` | Citation style |
| `PLAN_power_methods_section.md` | Methods-section planning note |

Driver scripts are at
`analysis/scripts/treatment-main-effect/`.

## Dependency on `nof1power`

The manuscript and its driver scripts assume `nof1power` is
installed and load it with

```r
library(nof1power)
```

Install with

```r
devtools::install_local('~/prj/res/06-nof1-power/nof1_power')
```

The `pmsimstats-ng` `DESCRIPTION` lists `nof1power` under `Suggests`.

### Function-name collisions to be aware of

`nof1power` and `pmsimstats` both export `lme_analysis`,
`modgompertz`, and `cumulative`. Whichever is loaded second masks
the first three. Manuscript scripts attached to this paper should
call only `library(nof1power)` and qualify any `pmsimstats` symbols
explicitly with `pmsimstats::` if needed.

## Path adjustments after migration

The original compendium used `here::here()` rooted at
`~/prj/res/06-nof1-power/nof1_power/`. After migration, `here()`
will resolve to the `pmsimstats-ng` root. Driver scripts that read
or write data, figures, or tables under the prior layout
(`analysis/data/`, `analysis/figures/`, `analysis/tables/`) need
their paths reviewed before re-running. Inspect each driver in
`analysis/scripts/treatment-main-effect/` before invoking.

## Build

Per project convention, this `report.Rmd` should load pre-computed
summaries rather than re-running simulations inline. The driver
scripts that produce those summaries live in
`analysis/scripts/treatment-main-effect/`.

## Author

pmsimstats team
