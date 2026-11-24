# 10-interaction-test-calibration
*2026-11-23 13:04 PST*

Factual file and build index, added by the 2026-11-23 deferred-item
pass (whitepaper of 2026-08-16, checklist item 16 and Minor Issue
3.2). It records what is in the directory and how the manuscript is
built. It does not summarize the findings; the Abstract of
`report.Rmd` is the authoritative statement of those.

## Manuscript

Title in `report.Rmd`: 'Type I calibration of a mixed-model
interaction test in an aggregated N-of-1 simulation: mean-model
misfit absorbed by REML and a cluster-robust remedy'.

Last render in this directory: `report.pdf`, 23 pages, created
2026-10-29 22:40 PDT (`pdfinfo`).

## Files

- `report.Rmd`: the manuscript. It runs no simulation; every table,
  figure and quoted number is computed inline from stored `.rds`
  files.
- `references.bib`, `statistics-in-medicine.csl`: bibliography and
  citation style named in the YAML header.
- `title-page.md`: submission front matter.
- `referee-report-2026-10-28.md`: the internal referee report.
- `type1-examination-plan.md` (and `.pdf`): the examination plan that
  preceded the manuscript.
- `bullets.Rmd` (and `.pdf`, `.tex`): a structured summary under an
  earlier title ('Type I calibration of fixed-effect interaction
  tests in linear mixed models: a staged diagnosis of
  working-covariance misspecification and a cluster-robust remedy').
  It is not maintained with the manuscript.
- `report.tex`, `report_files/`, `tmp-pdfcrop-50901.tex`: build
  by-products.
- `claudecode.tex`: role not documented here.

There is no `cover-letter.md`.

## Inputs

All inputs are under `analysis/scripts/carryover-sensitivity/output/`
(paths relative to the repository root):

- `paper10-refined/`: runs under the extended DGP, written by
  `run-paper10-refined.sh`. A table uses a file here when it exists
  with its full replicate count.
- `diag-*.rds`: the pre-correction run of 2026-05-20, used as the
  fallback; each table caption states which source it used.
- `paper10-refined/diag-kappa-analytic.rds`: the analytic file,
  always used.

## Rebuilding

From the repository root:

```
bash analysis/scripts/carryover-sensitivity/run-paper10-refined.sh
bash tools/render.sh \
  analysis/report/10-interaction-test-calibration/report.Rmd
```

Render only through `tools/render.sh`, never `rmarkdown::render()`.
