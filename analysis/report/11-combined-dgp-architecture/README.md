# 11-combined-dgp-architecture
*2026-11-23 13:04 PST*

Factual file and build index, added by the 2026-11-23 deferred-item
pass (whitepaper of 2026-08-16, checklist item 16 and Minor Issue
3.2). It records what is in the directory and how the manuscript is
built. It does not summarize the findings; the Abstract of
`report.Rmd` is the authoritative statement of those.

## Manuscript

Title in `report.Rmd`: 'A dual-channel data-generating process for
biomarker-treatment interactions in aggregated N-of-1 trials: two
parameterizations of one conditional-mean slope'.

Last render in this directory: `report.pdf`, 22 pages, created
2026-10-29 13:27 PDT (`pdfinfo`).

## Files

- `report.Rmd`: the manuscript. It includes `results-main.Rmd` as a
  child when the stored main run is present, and `results-main.Rmd`
  in turn includes `results-extra.Rmd`. If the stored run is absent,
  the text of `results-pending.md` is inserted instead.
- `references.bib`, `statistics-in-medicine.csl`: bibliography and
  citation style named in the YAML header.
- `referee-report-2026-10-28.md`: the internal referee report.
- `planning/plan-architecture-c.md`: the 2026-06-10 plan for the
  combined DGP, written when this material was part of paper 01.
- `bullets.Rmd` (and `.pdf`, `.tex`): a structured summary under an
  earlier title that names 'channel super-additivity'. The current
  manuscript withdraws that claim (commit `6ea37805e`), so the
  summary is superseded and is not maintained with the manuscript.
- `report.tex`: build by-product.
- `claudecode.tex`: role not documented here.

There is no `cover-letter.md`.

## Inputs

All inputs are in `analysis/data/dgp-combined/` (relative to the
repository root). The main grid has untagged file names
(`01-combined-*.rds`, `02-benchmarks.rds`, `03-analytic.rds`); the
0.2-week half-life pass and the G9 pass carry the suffixes `-hl02` and
`-g9`. The run of 2026-07-31 reported by the earlier draft is kept in
`archive-2026-07-31/`.

The stored main run passed its validation gate: all 1,800 paired
boundary comparisons in `01-combined-gate.rds` show zero difference in
data, estimate and p-value against the single-encoding DGPs, and
`01-combined-meta.rds` records `gate_pass = TRUE` (checked
2026-11-23).

## Rebuilding

The drivers are in `analysis/scripts/dgp-combined/`:
`01-run-combined-factorial.R` (run through `run-combined.sh`, which
sets the BLAS thread variables to 1; its header lists the exact
invocations), `02-summarise-combined.R` and
`03-analytic-projection.R`. The stored main run took about 169
minutes on eight workers (`elapsed_secs` in the meta file).

Render from the repository root, only through the wrapper:

```
bash tools/render.sh \
  analysis/report/11-combined-dgp-architecture/report.Rmd
```
