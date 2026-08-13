# latent-class-mixture-application

*2026-05-06 10:17 PDT*

## Working title

Latent-class and mixture-model formulations for biomarker-treatment
interaction in N-of-1 trials.

## Origin

This manuscript was scaffolded on 2026-05-06 in response to a
reviewer comment on the companion paper
`analysis/report/01-dgp-mean-moderation-vs-mvn/`. The reviewer asked:
if Architecture B is functionally an imperfect responder indicator,
why model the heterogeneity through covariance structure rather
than through an explicit participant-level latent-variable
construction? This paper takes that suggestion as its scientific
question.

The full exposition of the latent-variable framing that motivates
the paper is preserved verbatim at
`notes/source-revision-latent-class.md`. That note develops ten
threads, the first five concerning statistical substance and the
last five concerning the psychometric and econometric literatures
from which the relevant machinery is drawn. `report.Rmd` has since
developed past the scaffold stage into a full manuscript reporting a
240-replicate Study 5 pilot; the extended literature review and
taxonomy material that note motivated now live in `supplement.Rmd`.

## Directory contents

| File | Role |
|---|---|
| `report.Rmd` | Main manuscript (abstract, DGP, MVN approximation, identifiability, Study 5 pilot results) |
| `supplement.Rmd` | Companion supplement: extended literature review, full finite-mixture taxonomy, full biomarker-moderation spectrum, secondary theoretical notes, and pre-flight pilot detail. Split out from `report.Rmd` on 2026-08-12 to bring the main manuscript to a working page length; section numbers and citation keys match the main text |
| `bullets.Rmd` | Structured, section-by-section bullet summary of the manuscript (reading aid) |
| `whitepaper-latent-class-mixture-summary.md` | Two-page white-paper summary |
| `references.bib` | Seeded from `analysis/report/references.bib`; trim and extend as the paper develops |
| `statistics-in-medicine.csl` | Citation style |
| `notes/source-revision-latent-class.md` | Source exposition that motivates this paper |

A `_cache/` and rendered PDFs appear once the Rmds are knit; `report.Rmd`
and `supplement.Rmd` render independently via `tools/render.sh`.

## Scope

The paper has four developmental tasks:

1. Formalise the latent-class and finite-mixture formulations of
   the biomarker-treatment interaction problem and document
   identifiability conditions specific to the N-of-1 design family.
2. Characterise estimation: EM and Bayesian samplers, candidate R
   packages (`mclust`, `flexmix`, `OpenMx`, `Stan`).
3. A simulation study mirroring the carryover-sensitivity factorial
   with an added latent-class-separability axis.
4. Re-analysis of the prazosin-PTSD data alongside the
   Architecture-A and Architecture-B fits.

## Relationship to other papers in this repository

- `dgp-mean-moderation-vs-mvn/`: scientific origin; the reviewer
  question that this paper addresses lives in the revisions
  subdirectory of that paper.
- `carryover-sensitivity/`: provides the simulation infrastructure
  this paper will extend.
- `carryover-analysis-model-assessment/`: methodological cousin;
  shares the analysis-model-misspecification framing.
- `treatment-main-effect/`: parallel-vs-N-of-1 power comparison;
  related but addresses the main-effect rather than the
  biomarker-interaction question.

## Driver scripts

A directory at `analysis/scripts/latent-class-mixture-application/`
will hold simulation drivers once the simulation study is designed.
That directory is not yet created.

## Author

pmsimstats team
