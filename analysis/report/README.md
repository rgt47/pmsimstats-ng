# analysis/report/

Shared infrastructure consumed by the manuscripts under
`analysis/report/`. Per-manuscript sources (`report.Rmd`,
paper-specific `references.bib`, `statistics-in-medicine.csl`)
live inside each manuscript's own subdirectory.

## Contents

- `data/` holds numerical outputs checkpointed from simulation
  drivers in `analysis/` (e.g., `analysis/scripts/carryover-sensitivity/`
  writes `02-grid-summary.rds` and `02-sensitivity-summary.rds`
  here).
- `figures/` holds rendered publication figures produced by the
  same drivers.
- `tables/` holds publication tables.
- `references.bib` is a master candidate bibliography; individual
  manuscripts may use a paper-specific `references.bib` inside
  their own subdirectory.

## Current manuscripts

Each lives under `analysis/report/<slug>/` with `report.Rmd` as
the canonical source.

| Slug                                    | Title (working)                                                                                  |
|----------------------------------------|--------------------------------------------------------------------------------------------------|
| `dgp-mean-moderation-vs-mvn`           | Two data-generating architectures for biomarker-treatment interaction in N-of-1 trials          |
| `carryover-sensitivity`                 | Robustness of carryover-mitigation analysis strategies for biomarker-treatment interaction       |
| `carryover-analysis-model-assessment`   | Assessing carryover in the analysis model: standard methodology and application to pmsimstats    |

## Reproducibility

Simulation runs that feed a manuscript are orchestrated from
`analysis/` (per project convention) and checkpointed here under
`data/`. Each manuscript `.Rmd` loads finalised results rather
than re-running simulations inline.

## Relationship to `docs/`

- `docs/` contains the internal technical reports, audit notes,
  and white papers used to develop the arguments.
- `analysis/report/` contains the journal-format distillations
  intended for external submission.

A manuscript may cite specific `docs/` reports as supplementary
material when a journal permits linking to a companion technical
report; in particular, the
`analysis/report/dgp-mean-moderation-vs-mvn/longform.tex`
document is the extended technical companion to that
manuscript's `report.Rmd`.

## Splitting out at submission time

If a target journal requires a standalone submission bundle,
use `git subtree split` on the manuscript subdirectory to create
a focused repository preserving the relevant history.
