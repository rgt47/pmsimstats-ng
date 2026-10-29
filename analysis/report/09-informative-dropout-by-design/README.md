# 09-informative-dropout-by-design

*2026-05-07 16:30 PDT; revised 2026-10-29 09:49 PDT*

## Working title

Informative dropout and trial-design choice in aggregated N-of-1
biomarker-validation trials: a paired simulation at one parameter
point.

## Status (2026-10-29)

The manuscript was revised in response to the internal referee report
`referee-report-2026-10-28.md`. It now reports a paired
(common-random-numbers) simulation in which every dropout condition is
applied to the same complete datasets. The earlier unpaired prototype
is kept only as a point of comparison.

What the manuscript reports:

- One parameter point: total $N = 70$ per design, mean-shift encoding
  with $\beta_{bm} = 0.45$, carryover half-life 0.5 week.
- Three designs: OL+BDC, CO and Hybrid. The OL design is excluded
  because it has no off-drug contrast and the package analysis tests
  `bm:t` for it, not the interaction.
- Seven dropout conditions: none, MCAR 25 percent,
  magnitude-of-change 25 and 40 percent (the former 'biased'
  hazard), worsening-driven, improvement-driven and symptom-return
  25 percent.
- Two analysis half-lives (0.5 week, the historical choice, and 0,
  which matches the generating interaction) with model-based and CR2
  inference.
- A null regime ($\beta_{bm} = 0$) for test size.

What was not run: the four-study protocol grid (sample size,
effect size and carryover sweeps), a path-stratified decomposition
of the 'happy accident', analysis-stage sensitivity analyses
(pattern-mixture, delta-adjusted or reference-based multiple
imputation), and shared-parameter dropout.

## Files

- `report.Rmd`: the manuscript. Numbers in it are read from the
  summary files below by inline R.
- `references.bib`, `statistics-in-medicine.csl`: bibliography and
  style.
- `cover-letter.md`, `title-page.md`: submission front matter.
- `referee-report-2026-10-28.md`: the internal referee report.
- `bullets.Rmd`: an outline file, not maintained with the
  manuscript.

## Drivers and data

- Paired driver:
  `analysis/scripts/informative-dropout-by-design/20-paired-dropout.R`,
  run through `20-run-paired.sh` (single-threaded BLAS).
- Summaries:
  `analysis/scripts/informative-dropout-by-design/21-summarise-paired.R`.
- Output: `analysis/data/informative-dropout-by-design/paired-*`.
- Superseded unpaired prototype:
  `analysis/scripts/quick-sim/09-dropout-prototype-10min.R`, output
  `analysis/data/quick-sim/09-dropout-*`.
- Working protocol (not executed as planned):
  `analysis/scripts/informative-dropout-by-design/00-ademp-protocol.md`
  and the unrun scaffolds `00-common.R`, `01-study1.R` to
  `04-study4.R`, `run-all.R`.

## Rendering

```bash
bash tools/render.sh \
  analysis/report/09-informative-dropout-by-design/report.Rmd
```

## Author

pmsimstats team
