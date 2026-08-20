# Paper 01: DGP architectures for biomarker-treatment interaction
*2026-07-31 11:51 PDT*

Directory conventions and draft lineage for the
mean-moderation-versus-MVN manuscript.

## Master source

`report.Rmd` is the single maintained source. Its scope is:

- Two DGP architectures: mean moderation (the biomarker scales the
  treatment effect in the conditional mean) and covariance
  moderation (differential correlation, the multivariate normal
  approach).
- Three within-subject designs: crossover (CO), Hybrid, and
  open-label plus blinded discontinuation (OL+BDC).
- Carryover half-lives of 0, 0.5, and 1.0 weeks, with $N = 35$ per
  randomization path and an $N = 70$ robustness check.
- Tutorial framing aimed at applied statisticians new to the N-of-1
  setting, including an N-of-1 primer in the Introduction.
- The three-part `bullets` / `rgt` / `orig` paragraph scaffolding.
  As of 2026-08-16, `report.Rmd` carries zero remaining `rgt`
  placeholder markers (verified by grep); this manuscript has
  reached zero-placeholder state, unlike most of the other series
  in this compendium.

All numerical results are typeset inline; the document loads no
data and runs no simulation code during knit. It requires
`claudecode.tex` for the scaffolding environments, `references.bib`,
and `statistics-in-medicine.csl`, all in this directory.

Render with:

```bash
Rscript -e "rmarkdown::render('analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd')"
```

The current render is 21 pages.

## Other live files

- `title-page.md` is submission material.
- `references.bib`, `statistics-in-medicine.csl` are the master's
  bibliography and citation style, shared with the archived drafts.
- `claudecode.tex` supplies the `bullets` / `rgt` / `orig`
  environments. A copy lives in `archive/` for the archived
  comprehensive draft.
- `whitepaper-dgp-architectures-summary.md` and
  `whitepaper-architecture-c-latent-class-assessment.md` are
  standalone summary whitepapers.
- `architecture-c-extraction-plan.md` is the plan that moved the
  dual-channel architecture to paper 11. Its file names and line
  numbers refer to the pre-consolidation layout, in which the
  master was named `report-slim.Rmd`.
- `referee-report-2026-06-13.md` and `revisions/` are historical
  records, not manuscript inputs.

## Directory layout

```
report.Rmd                 master manuscript
archive/                   superseded manuscript drafts
revisions/                 referee reports and revision planning
```

## Archived drafts

Two drafts preceded the master. Each is retained with its rendered
PDF. Neither is maintained; the master supersedes both.

| File | Scope | Superseded because |
|---|---|---|
| `archive/report.Rmd` | Comprehensive two-architecture manuscript: numbered subsections throughout Sections 4, 6, and 7, a TOC, line numbers, and full journal back matter | Tutorial rewrite adopted; its long-form subsections were condensed into narrative sections |
| `archive/report-with-responses.Rmd` | Response draft of 2026-05-06: the older two-architecture body with tracked changes (`\add`, `\del`, `\repl`) plus a point-by-point reply to eleven reviewer comments | Superseded by the revisions it proposed |

**Content held only by the archived comprehensive draft.** The
following is not in the master and must be taken from
`archive/report.Rmd` if needed:

- The numbered biology subsections 4.1 (when mean moderation is
  appropriate), 4.2 (when covariance moderation is appropriate),
  and 4.3 (the prazosin-PTSD case), which the master compresses
  into the single narrative Section 4.
- The seven numbered alternative-analysis subsections 6.1 through
  6.7 (on-drug-only restriction, weighted analysis, within-subject
  contrast, two-stage random slopes, exclusion of contaminated
  observations, design-level solutions, qualitative comparison),
  which the master compresses into Section 6 plus a single
  comparison subsection 6.1.
- The four numbered literature subsections 7.1 through 7.4
  (treatment effect heterogeneity, prevalence of each architecture,
  crossover and N-of-1 methodology, precision medicine trial
  design), which the master compresses into Section 7.
- Separate journal back matter under its own headings: conflict of
  interest, funding, data availability, reproducibility. The master
  carries this material condensed into a single
  "Notes and reproducibility" section.

**Content held only by the archived response draft.** The editorial
preface and the point-by-point responses C1 through C11 to the
reviewer (Rebecca C. Hendrickson) exist nowhere else in this
directory. The tracked-changes markup itself is of historical
interest only.

**Rendering an archived draft.** Their `bibliography` and `csl`
paths were repointed one level deeper on archiving
(`../references.bib`, `../statistics-in-medicine.csl`), and
`archive/claudecode.tex` supplies the scaffolding environments the
comprehensive draft requires. Both were re-rendered in place on
2026-07-31 and produced their PDFs without error.

## The dual-channel architecture (paper 11)

A third architecture, activating the mean and covariance channels
simultaneously with independent parameters $c_{bm,A}$ and
$c_{bm,B}$, was developed in this manuscript through 2026-07-01 and
then extracted to paper 11
(`analysis/report/11-combined-dgp-architecture/`) under
`architecture-c-extraction-plan.md`. The extraction was applied to
the comprehensive draft on 2026-07-01 but not to the draft that
became the master, which continued to carry the material until the
consolidation of 2026-07-31.

The master now retains only forward-pointing references to paper
11, in the Introduction, at the close of Section 2.2, in Section 4,
and in Section 5.1. Material removed from the master during
consolidation, and available from paper 11 or from the
pre-consolidation history:

- Subsection 2.3 defining the combined architecture (sequential
  MVN-draw-then-mean-shift construction, algebraic boundary
  reductions, the $\alpha$ mixing parameterization) and its
  paragraph in the Mathematical Comparison.
- Subsection 3.4 and its $3 \times 3$ power table
  (`tab:arch-c-power`), covering the $c_{bm,A} \times c_{bm,B} \in
  \{0, 0.22, 0.45\}^2$ grid at $t_{1/2} = 1.0$ week, together with
  the super-additivity commentary.
- The Architecture C column of the Section 5.1 decision table and
  the mixing-weight guidance that followed it.

**Caveat on the archived boundary-validation claim.** The
subsection 3.4 that was removed asserted that the combined
parameterization recovers each single-channel limit, while also
noting that its boundary cells are not numerically comparable to
the Section 3.1 tables because that grid used $N = 35$ as a total
rather than per randomization path. The extraction plan records
this as referee finding MC-1 and assigns the corrective re-run to
paper 11. Read the archived grid with that correction in mind.

## Known defects in the master

- One citation renders as a raw citekey: `[@pmsimstats-paper08]` in
  the first Introduction bullet block. Citations placed inside the
  raw-LaTeX `bullets` environments bypass pandoc citation
  processing. Pre-existing, not introduced by consolidation.
- The `rgt` blocks throughout are placeholders reading
  `rgt to complete.` and await the author's prose. The
  `~/bin/strip-claudecode` script finalizes the paper once that
  writing is done.
