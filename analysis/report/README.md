# analysis/report/

Shared infrastructure consumed by the manuscripts under
`analysis/report/`. Per-manuscript sources (`report.Rmd`,
paper-specific `references.bib`, `statistics-in-medicine.csl`)
live inside each manuscript's own subdirectory.

## Contents

- `sim-preamble.tex` is the shared LaTeX preamble used by all
  manuscripts (line numbers, math-line patches, page footer with
  source path and render timestamp). Each `report.Rmd`
  references it via `includes: in_header: ../sim-preamble.tex`
  in the YAML.
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

## Canonical YAML for all manuscripts

All manuscripts share an identical YAML targeting Statistics in
Medicine submission style (continuous line numbers, 1.5 line
spacing, 12pt, asymmetric geometry to leave wide right margin
for editor and reviewer comments, numbered sections, TOC depth 2,
SIM CSL bibliography). The canonical YAML is:

```yaml
---
title: "[paper title]"
author:
  - pmsimstats team
date: "`r format(Sys.Date(), '%Y-%m-%d')`"
output:
  bookdown::pdf_document2:
    latex_engine: xelatex
    keep_tex: true
    toc: true
    toc_depth: 2
    number_sections: true
    includes:
      in_header: ../sim-preamble.tex
bibliography: references.bib
csl: statistics-in-medicine.csl
link-citations: true
linestretch: 1.2
fontsize: 12pt
geometry: "left=2.5cm,right=6cm,top=2.5cm,bottom=2.5cm"
---
```

The `treatment-main-effect/report.Rmd` and
`treatment-main-effect/report_short.Rmd` use
`bibliography: "nof1-pgt.bib"` instead of `references.bib` (the
two are kept distinct in that paper's directory). Otherwise the
YAML is identical across manuscripts.

After the YAML and the `setup` chunk, every manuscript also
includes a `footer-setup` chunk that injects the source-file path
and render timestamp into the page footer:

```r
{r footer-setup, echo = FALSE, results = 'asis'}
if (knitr::is_latex_output()) {
  src <- tryCatch(normalizePath(knitr::current_input()),
                  error = function(e) 'unknown')
  src_disp <- gsub('^/Users/zenn/Dropbox', '~/Dropbox', src)
  src_disp <- gsub('^/Users/zenn', '~', src_disp)
  src_tex <- gsub('([_~^&%#$])', '\\\\\\1', src_disp)
  ts <- format(Sys.time(), '%Y-%m-%d %H:%M %Z')
  cat(sprintf('\\renewcommand{\\rmdsource}{Source: \\texttt{%s}}\n',
              src_tex))
  cat(sprintf('\\renewcommand{\\rendertimestamp}{Rendered: %s}\n',
              ts))
}
```

The `\rmdsource` and `\rendertimestamp` macros are defined as
placeholders in `sim-preamble.tex` and overridden at knit time.

## Current manuscripts

Each lives under `analysis/report/<NN-slug>/` with `report.Rmd`
as the canonical source. Manuscripts are numbered in dependency
and build order: foundational papers first, extensions and
companions after. The numerical prefix is part of the directory
name and is used in cross-references throughout the repo.

| Slug                                       | Title (working)                                                                                  |
|--------------------------------------------|--------------------------------------------------------------------------------------------------|
| `01-dgp-mean-moderation-vs-mvn`            | Two data-generating architectures for biomarker-treatment interaction in N-of-1 trials (with `longform.tex` companion: 'Two Architectures for Simulating ...'). Foundational DGP comparison. |
| `02-carryover-sensitivity`                 | Robustness of carryover-mitigation analysis strategies for biomarker-treatment interaction. Extends `01-` with three carryover decay forms and three analysis-side specifications. |
| `03-latent-class-mixture-application`      | Latent-class and mixture-model formulations for biomarker-treatment interaction in N-of-1 trials. Extends `01-` with mixture and class-aware DGPs; companion simulation plan committed. |
| `04-treatment-main-effect`                 | Comparative statistical power of N-of-1 and parallel-group designs for detecting the treatment main effect (migrated from `~/prj/res/06-nof1-power/` on 2026-05-06; the `nof1power` simulation engine was vendored to `implementations/nof1power/` on 2026-05-08). Clinical companion. |
| `05-nof1-design-sensitivity`               | Power sensitivity of a hybrid OL+BDC N-of-1 design under carryover, serial correlation, and variance-component variation (split off from `04-` on 2026-05-06). Methodological companion. |
| `06-component-decomposition`               | Three-component decomposition of treatment response in aggregated N-of-1 trials: pharmacological (BR), expectation-driven (PB), and natural-history (TV) components. Methodological foundation for biomarker-validation analyses across `01-`, `02-`, and `03-`. |
| `07-gompertz-evaluation`                   | Critical evaluation of the modified Gompertz response function as a parametric template for treatment-response trajectories in N-of-1 simulation, with logistic / hyperbolic-tangent / piecewise-linear breakpoint sensitivity comparators. |
| `08-test-procedure-design-sensitivity`     | Test-procedure (strict RM-ANOVA, linear-mixed, GEE) and trial-design (cycle count, period length, on/off symmetry) choices for the biomarker-treatment interaction in aggregated N-of-1 trials. Simulation substrate: `implementations/simple/`. |
| `09-informative-dropout-by-design`         | Informative dropout (hazard depending on cumulative symptom worsening) and its coupling with the four canonical N-of-1 design families (OL, OL+BDC, traditional crossover, hybrid). Reproduces and extends Hendrickson 2020 Figure 4A; quantifies the 'happy accident' randomization-path selection effect. |

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
`analysis/report/01-dgp-mean-moderation-vs-mvn/longform.tex`
document is the extended technical companion to that
manuscript's `report.Rmd`.

## Sharing manuscript drafts with collaborators

Every render automatically copies the produced PDF into a
project-local staging directory at `analysis/report/share/`
with a reproducible filename
`<NN-slug>[-short]-<YYYY-MM-DD>-d<N>.pdf`, where `<N>` is the
per-day draft counter (incremented for every render on the same
calendar day). The original `report.pdf` continues to land in
the manuscript directory unchanged; the staged copy is added
alongside.

Three invocation paths all trigger staging:

- **RStudio Knit button**: respects each Rmd's YAML `knit:` hook
  and routes the render through `tools/stage-render.R`.
- **Terminal wrapper**: `bash tools/render.sh <paper-slug>` (or
  `--short` for paper 04's short variant).
- **R session**: source the helper and call it directly with
  `source('tools/stage-render.R')$value('path/to/report.Rmd')`.

The staging directory is gitignored except for `MANIFEST.md`
(auto-generated by the render hook with one row per
render-and-stage event listing the staged filename, source path,
and render timestamp). To share a specific draft with a
collaborator, copy the file from `analysis/report/share/` and
attach it to email or upload to wherever they prefer.

Plain `Rscript -e "rmarkdown::render('report.Rmd')"` calls do
**not** trigger staging because they bypass the YAML `knit:`
hook. Use the terminal wrapper or the R-session helper for
distribution-ready renders.

## Splitting out at submission time

If a target journal requires a standalone submission bundle,
use `git subtree split` on the manuscript subdirectory to create
a focused repository preserving the relevant history.
