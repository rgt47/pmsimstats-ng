# Paper 02: carryover-mitigation analysis strategies
*2026-07-30 17:12 PDT*

Directory conventions and draft lineage for the carryover-sensitivity
manuscript.

## Master source

`report.Rmd` is the single maintained source. Its scope is:

- Architecture B (MVN differential correlation) only.
- Exponential and Weibull DGP decay forms.
- Analysis specifications E2 (exposure-weighted) versus E3
  (lagged-treatment).
- Sensitivity blocks S1 through S4, plus the S6 cluster-robust
  recalibration.
- 216 cells at 500 replicates, filtered from the 540-cell production
  grid.

It consumes the `02xs-*` filtered figures under `analysis/figures/`
and the shared cell-level summaries under `analysis/data/`. It runs
no simulation code during knit. The build pipeline is documented in
the HTML comment at the top of the file.

Render with:

```bash
Rscript -e "rmarkdown::render('analysis/report/02-carryover-sensitivity/report.Rmd')"
```

## Other live files

- `report-devresults.Rmd` is not a manuscript variant. It is a
  development-run review document: 150 replicates per cell with the
  `--robust` flag, five analysis methods including lme+CR2 and
  GEE+MD, all figures computed inline from
  `output/01-dev.rds` and `output/04-sensitivity-dev.rds`. Use it for
  ranking and direction checks, not for power estimates.
- `cover-letter.md`, `title-page.md` are submission materials.
- `references.bib`, `statistics-in-medicine.csl` are the master's
  bibliography and citation style.
- `whitepaper-carryover-specification-summary.md` is the standalone
  summary whitepaper.

## Directory layout

```
report.Rmd                 master manuscript
report-devresults.Rmd      development-run review (not a variant)
archive/                   superseded manuscript drafts
revisions/                 referee reports and revision planning
```

## Archived drafts

Four drafts preceded the master, in this order. Each is retained with
its final rendered PDF. None is maintained; the master supersedes all
of them.

| File | Scope | Superseded because |
|---|---|---|
| `archive/report.Rmd` | Full: both architectures, three decay forms (linear, exponential, Weibull), E1+E2+E3, S1-S5, full front and back matter, TOC | Scope narrowed for a self-contained evaluation |
| `archive/report-slim.Rmd` | Same scope; condensed Discussion, adds an N-of-1 primer, no TOC | Same |
| `archive/report-extra-slim.Rmd` | Narrowed to the master's scope; carries the `bullets`/`rgt`/`orig` three-part scaffolding with `rgt` placeholders | Scaffolding stripped and prose rewritten |
| (promoted to `report.Rmd`) | Narrative rewrite of the above | Current master |

**Content held only by the archived full-scope drafts.** The
following is not in the master and must be taken from `archive/` if
needed:

- Architecture A (mean moderation) results and the Discussion
  subsection on architecture dependence.
- The linear decay form.
- The E1 (binary on-drug) specification results.
- Sensitivity block S5 (autocorrelation by carryover, exploratory).
- Discussion subsections on the latent-subtype question and on
  comparison with prior simulation studies, and the long-form
  Limitations section.

The journal back matter (conflict of interest, funding, data
availability, reproducibility) was restored into the master during
consolidation and does not need to be recovered from `archive/`.

**Rendering an archived draft.** Their relative paths were repointed
one level deeper on archiving (`../../../figures/`,
`../references.bib`, `../statistics-in-medicine.csl`), and
`archive/claudecode.tex` supplies the scaffolding environments they
require. They should therefore still render in place, though this has
not been verified since the move. They read the full-scope `02-*`
figures, which are produced by the unsuffixed scripts
`03-render-figures.R`, `04-sensitivity-figures.R`, and
`07-type1-cr2-figure.R` under
`analysis/scripts/carryover-sensitivity/`. Those scripts are retained
solely for the archived drafts and are not part of the master's
build.

**Caveat on archived claims.** The full-scope drafts state the
exposure-weighted-dominance result more broadly than the data
support. The master qualifies it: the exposure-weighted
specification leads only under the Hybrid and OL+BDC designs and is
markedly inferior under the classical crossover design. Read the
archived drafts with that correction in mind.

**Specification labels.** The simulation output stores the three
analysis specifications as `A1`, `A2`, `A3`; every document here
displays them as E1, E2, E3. The mapping is applied at display time
in the `spec_display()` helper in each document's setup chunk, and
the stored values are never rewritten.

## Revision history

`revisions/` holds the referee reports of 2026-06-16 (full and slim
variants), the revision plan of the same date, the sister-paper
consistency review against paper 01, and the readability and
consistency review with its supporting figures. These are historical
records, not manuscript inputs.
