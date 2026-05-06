# latent-class-mixture-application

*2026-05-06 10:17 PDT*

## Working title

Latent-class and mixture-model formulations for biomarker-treatment
interaction in N-of-1 trials.

## Origin

This manuscript was scaffolded on 2026-05-06 in response to a
reviewer comment on the companion paper
`analysis/report/dgp-mean-moderation-vs-mvn/`. The reviewer asked:
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
from which the relevant machinery is drawn. The current `report.Rmd`
is a section-stubbed scaffold that the manuscript will fill in.

## Directory contents

| File | Role |
|---|---|
| `report.Rmd` | Manuscript skeleton (sections stubbed) |
| `references.bib` | Seeded from `analysis/report/references.bib`; trim and extend as the paper develops |
| `statistics-in-medicine.csl` | Citation style |
| `notes/source-revision-latent-class.md` | Source exposition that motivates this paper |

A `_cache/` and rendered PDFs will appear once the Rmd is knit.

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
