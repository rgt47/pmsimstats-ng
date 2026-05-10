# nof1-design-sensitivity

*2026-05-06 18:28 PDT*

## Working title

Power sensitivity of a hybrid open-label-plus-blinded-discontinuation
N-of-1 design under carryover, serial correlation, and
variance-component variation.

## Provenance

This directory holds the methodological companion to the clinical
manuscript at `analysis/report/04-treatment-main-effect/`. The two
manuscripts were originally drafted as a single document at
`treatment-main-effect/report.Rmd`; the methodological sensitivity
content was split off on 2026-05-06 because the combined draft
exceeded the typical word budget of biostatistics journals
(*Statistics in Medicine*, *Statistical Methods in Medical Research*,
*Pharmaceutical Statistics*).

The simulation engine consumed by this manuscript is the `nof1power`
R package, originally developed in the standalone compendium at
`~/prj/res/06-nof1-power/nof1_power/` and migrated alongside the
clinical manuscript on 2026-05-06.

## Scope

Twelve parameter-sensitivity sweeps (S1 through S12) executed against
the pmsimstats-ng tidyverse pipeline with the Hendrickson 2020 hybrid
N-of-1 design as the primary aggregated variant and the alternating
A/B/A/B variant retained as a sensitivity comparator at sweeps that
vary cycle count or period length. The target estimand is the main
treatment effect (the `Dbc` coefficient in the `nlme::lme`
fixed-effects formula `Sx ~ bm + t + Dbc + bm:Dbc`).

## Baseline parameters

- N = 35 participants per design (Hendrickson 2020 minimum,
  per-participant parity with the parallel RCT comparator)
- 8 weekly timepoints
- Effect size approximately -2 nightmares/week (prazosin-calibrated)
- Carryover half-life 3 days (exponential decay)
- AR(1) within-factor autocorrelation 0.3
- 500 replicates per cell for S1-S11; 2,000 replicates per cell for
  S12 (null calibration)

## Sweep index

| Sweep | Axis | Precedent |
|-------|------|-----------|
| S1 | Effect size | Blackston 2019 |
| S2 | Carryover half-life | Blackston 2019, Hendrickson 2020 |
| S3 | Decay-form misspecification | Gärtner 2023 |
| S4 | Between-patient heterogeneity | Senn 2018, Zucker 2010 |
| S5 | Within-patient SD | Wang-Schork 2019 |
| S6 | Cycles per patient k | Senn 2018 |
| S7 | Total participant count N | Senn 2018, Blackston 2019 |
| S8 | AR(1) correlation rho | Wang-Schork 2019 |
| S9 | Period length crossed with carryover | Wang-Schork 2019 |
| S10 | Biomarker interaction crossed with DGP architecture | Hendrickson 2020 |
| S11 | Carryover-by-AR(1) misspecification factorial | Gärtner 2023 |
| S12 | Null-effect calibration (Type I) | Morris 2019 |

## Directory contents

| File | Role |
|------|------|
| `report.Rmd` | Canonical manuscript source |
| `report.tex`, `report.pdf` | Last rendered LaTeX and PDF |
| `references.bib` | Manuscript-specific bibliography |
| `statistics-in-medicine.csl` | Citation style |
| `README.md` | This file |

Driver scripts live at
`analysis/scripts/nof1-design-sensitivity/`, with figures and
RDS artifacts written to `analysis/figures/sensitivity/` and
`analysis/data/derived_data/sensitivity/` respectively.

## Build

```bash
Rscript analysis/scripts/nof1-design-sensitivity/run-all.R
Rscript analysis/scripts/nof1-design-sensitivity/plot-all.R
Rscript -e "rmarkdown::render('analysis/report/05-nof1-design-sensitivity/report.Rmd')"
```

`run-all.R` takes approximately 32 minutes of wall-clock time on an
eight-core machine using the furrr multisession backend.

## Author

pmsimstats team
