# Reading guide and index

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Three reading paths, chapter status, and old-to-new map
for every docs/NN.
**Prerequisites:** None

---

## About this series

The pmsimstats-ng whitepaper series documents the design and
evaluation of a Monte Carlo simulation tool for detecting
biomarker-treatment interactions in aggregated N-of-1 clinical
trials. The series covers: the data-generating process (DGP),
including its published predecessor and the defects corrected in
the refined process; the coupling ceiling and its dependence on
covariance structure; carryover and its effect on power; analysis
model selection; and sample-size planning.

The series does not report clinical trial results and does not
provide pre-specified or confirmatory analysis plans. The simulation
studies are exploratory. The series is not a user manual for the
pmsimstats-ng R package; package documentation lives in `man/` and
`vignettes/`.

For all project-wide notation, see `fm2` (Notation, vocabulary and
conventions). The canonical source is
`analysis/report/NOTATION.md`; where `fm2` and `NOTATION.md`
differ, `NOTATION.md` governs.

---

## Contents

1. [Reading paths](#reading-paths)
2. [Chapter status](#chapter-status)
3. [Old-to-new map](#old-to-new-map)
4. [Change log](#change-log)

---

## Reading paths

<!-- From docs/59 Section 3.1 -->

Three paths through the series. Each path is a linear reading order;
chapters not on a path are accessible once their stated prerequisites
are met. See the chapter status table for the current state of each
chapter (all are stubs in Phase 2).

### Path A. New to N-of-1 trials

Readers who are new to aggregated N-of-1 trials and want to
understand the program's motivation and key findings before the
technical details.

`fm2` → `ch01` → `ch02` → `ch03` → `ch10` → `ch22`

- **fm2** -- notation and vocabulary; read before any chapter.
- **ch01** -- the question: what is a biomarker-treatment
  interaction, why N-of-1 trials can detect it, and where carryover
  complicates that goal.
- **ch02** -- trial designs, randomization paths and exposure
  clocks: the mechanics of how carryover arises and how the exposure
  indicator is constructed.
- **ch03** -- the published process (Hendrickson 2020) at commit
  58b32a9, documented before any modification.
- **ch10** -- the refined process: the four corrections made, their
  rationale, and the resulting specification.
- **ch22** -- the recommended analysis F1 term by term: what each
  part of the model does, what goes wrong without it, and what it
  costs.

### Path B. Collaborator reviewing the published process defects

Readers who want to understand what is wrong with the published
process and why the refined process addresses it.

`ch03` → `ch04` → `ch10`

- **ch03** -- the published code at commit 58b32a9 in detail: data
  structure, parameter choices, and the covariance construction.
- **ch04** -- the three defects: compound-symmetry covariance,
  step coupling profile, and compounding carryover decay; evidence
  and consequence for each.
- **ch10** -- the refined process: how each defect is addressed,
  the configuration selected, and its validation.

### Path C. Statistical reviewer checking the analysis recommendation

Readers who want to evaluate the analysis F1 and its evidence base.

`ch17` → `ch18` → `ch19` → `ch20` → `ch21` → `ch22`

- **ch17** -- the analysis model: exposure codings, time
  adjustment, and the role of baseline as an outcome row.
- **ch18** -- inference and calibration: CR2 cluster-robust
  standard errors, Satterthwaite degrees of freedom, and the size
  criterion.
- **ch19** -- the prognostic biomarker and the shared biomarker
  effect: why expectancy terms are necessary when the biomarker is
  not purely predictive.
- **ch20** -- the covariance study: 14 candidate analyses compared
  on 378,000 simulated trials.
- **ch21** -- the stress test: 41 cells, 28,000 trials, 11
  protocol analyses plus post hoc arms; protocol and results.
- **ch22** -- the recommended analysis F1, the maximin rule, and
  the basis for the recommendation.

---

### Additional reading paths

#### Methodologist

New to the project, or reviewing the main ideas before reading the
companion papers.

`fm2` → `ch01` → `ch07` → `ch10` → `ch21` → `ch22` → `ch23`

#### Statistician evaluating the analysis

`ch17` → `ch18` → `ch19` → `ch21` → `ch22` → `ch23`

#### Implementation reader

Readers who need to work with the code, patch the published
process, or understand the parity workflow.

`ch03` → `ch04` → `ch10` → `ch28` → `appD` → `appJ`

#### Manuscripts only

Readers consulting the companion papers who want only the key DGP
and analysis context.

`ch22` → `ch10` → `ch15` → `ch16`

---

## Chapter status

<!-- From docs/series/registry/chapters.csv -->

Status as of 2026-10-15 after Phases 2-4. **substantial** = most
sections have real content, construct box filled where applicable,
evidence labels present; **partial** = real prose in some sections
but material TODO items or incomplete extractions remain; **stub** =
section headers and source-pointer comments only.
Summary: 24 substantial, 16 partial, 1 stub (appB).

| ID | Title | Part | Needs | Est. pp. | Status |
|----|-------|------|-------|----------|--------|
| `fm0` | Reading guide and index | Front matter | | 6 | substantial |
| `fm1` | Executive summary | Front matter | | 4 | substantial |
| `fm2` | Notation, vocabulary and conventions | Front matter | | 8 | substantial |
| `ch01` | The question: a predictive biomarker in aggregated N-of-1 trials | I | `fm2` | 12 | partial |
| `ch02` | Trial designs, randomization paths and exposure clocks | I | `ch01` | 14 | substantial |
| `ch03` | The published process at commit `58b32a9` | I | `ch02` | 22 | substantial |
| `ch04` | Three defects of the published process | I | `ch03` | 18 | partial |
| `ch05` | Response components and the modified Gompertz trajectory | II | `ch02` | 16 | partial |
| `ch06` | Carryover on the mean | II | `ch05` | 10 | substantial |
| `ch07` | Interaction architectures: mean, covariance and dual-channel | II | `ch06` | 20 | partial |
| `ch08` | Biomarker coupling under carryover | II | `ch07` | 12 | partial |
| `ch09` | The response covariance and the coupling ceiling | II | `ch08` | 30 | substantial |
| `ch10` | The refined process: specification, the adjustments evaluated, and the roadmap | II | `ch04`, `ch09` | 18 | substantial |
| `ch11` | Features of an N-of-1 data-generating process | II | `ch10` | 26 | substantial |
| `ch12` | Component decomposition, identification and reduction | II | `ch05`, `ch07` | 22 | partial |
| `ch13` | Informative dropout and missingness | II | `ch02`, `ch10` | 10 | partial |
| `ch14` | Beyond the refined process: a pharmacodynamic sensitivity model | II | `ch09`, `ch10` | 10 | substantial |
| `ch15` | Closed-form theory: E9, the coupling gap and power | III | `ch08`, `ch09` | 30 | substantial |
| `ch16` | Design efficiency and the choice of schedule | III | `ch15` | 20 | substantial |
| `ch17` | The analysis model: exposure codings, time and baseline | IV | `ch15` | 22 | partial |
| `ch18` | Inference and calibration | IV | `ch17` | 16 | substantial |
| `ch19` | The prognostic biomarker and the shared biomarker effect | IV | `ch17`, `ch18` | 16 | partial |
| `ch20` | The covariance study: size of candidate analyses | IV | `ch18`, `ch19` | 16 | substantial |
| `ch21` | The stress test: protocol, results and the refined-process rerun | IV | `ch19`, `ch20` | 26 | substantial |
| `ch22` | The recommended analysis | IV | `ch21` | 22 | substantial |
| `ch23` | Choosing the analysis half-life from the data | IV | `ch22` | 30 | substantial |
| `ch24` | Testing the drug main effect when the half-life is unknown | V | `ch15`, `ch23` | 30 | substantial |
| `ch25` | Analyzing a real trial | V | `ch22`, `ch23`, `ch24` | 12 | partial |
| `ch26` | Note: the carryover recursion in the published code | VI | | 20 | substantial |
| `ch27` | Note: step coupling in the published code | VI | | 24 | substantial |
| `ch28` | Note: patching `58b32a9` | VI | `ch26`, `ch27` | 8 | partial |
| `appA` | Catalog of analyses and the ID registry | Appendix | | 20 | partial |
| `appB` | Decision history | Appendix | | 8 | stub |
| `appC` | Errata to earlier documents and manuscripts | Appendix | | 8 | partial |
| `appD` | Implementation status and open-items register | Appendix | | 8 | partial |
| `appE` | Lexicon | Appendix | | 16 | substantial |
| `appF` | Style guide | Appendix | | 6 | substantial |
| `appG` | Bibliography | Appendix | | 14 | partial |
| `appH` | Classical test equivalents: repeated-measures ANOVA | Appendix | | 14 | substantial |
| `appI` | The legacy record (`docs/00` to `docs/20`) | Appendix | | 40 | substantial |
| `appJ` | Software and parity workflow | Appendix | | 1 | substantial |

---

## Notation cross-reference

All project-wide notation is defined in `fm2` (Notation, vocabulary
and conventions). The canonical source is
`analysis/report/NOTATION.md`; where any chapter or appendix
departs from that source, `NOTATION.md` governs. Run
`perl tools/notation-lint.pl` on edited documents to catch
inconsistencies.

---

## Old-to-new map

<!-- From docs/59 Section 4 -->

Disposition of every file in `docs/`. "Archive" means `docs/legacy/`
with the file name unchanged and a one-line banner naming the
replacing chapter. "Remove from series" means the file moves out of
`docs/` (destination named) and is not part of the read-through.

| Old doc | New home | Notes |
|---------|----------|-------|
| `00-documentation-index.tex` | archive | superseded; replaced by the generated index `fm0` |
| `01-codebase-overview.tex` | archive | repository description belongs in `CLAUDE.md` and package docs; nothing unique |
| `02-dgp-mean-moderation-vs-mvn-rchcoments.pdf` | `docs/reviews/`, renamed `review-comments-2026-04-29.pdf`; Appendix I Section 8 | record of the co-author's comments; the questions enter Appendix I |
| `03-audit-and-revision-report.tex` | archive; Appendix I Sections 1, 9 | claims superseded; corrected commit attribution in `ch03` |
| `04-revised-power-analysis.tex` | archive; Appendix I Sections 1, 3 | commit history and Figure 4 reproduction are unique; mechanism and ceilings superseded |
| `06-ar1-residual-correlation.tex` | archive; Appendix I Section 4 | `lmer` Type I evidence unique; derivation superseded by CM |
| `07-positive-definiteness-failures.tex` | archive; Appendix I Section 1 | counts and condition numbers kept; root cause refuted by `docs/36` Section 8 |
| `07a-pd-constraints-theory.tex`, `.pdf` | archive | superseded and partly wrong; nothing carried |
| `08-biomarker-correlation-decay.tex` | archive; `ch08` credit; Appendix I Section 7 | origin of graded coupling; trait-against-state question unique |
| `09-carryover-correlation-artifact.tex` | archive; Appendix I Section 7 | superseded by `docs/55` and `docs/36`; clinical-plausibility point kept |
| `09a-carryover-bias-variance.md`, `.pdf`, `.tex` | archive all three; Appendix I Section 6 | unique content placed in Appendix I; the `.tex` contradicts the corrected `.md` |
| `13-figure4-walkthrough.tex` | archive; Appendix I Section 2 | grid and censoring specification unique; DGP sections mislabel the 2024 code |
| `14-trial-analysis-implementation-guide.tex` | archive; Appendix I Section 9; `ch25` | goals record; formulas superseded by the 2026-10-08 code fix |
| `15-tidyverse-alignment-plan.tex` | archive; Appendix I Section 10 | planning record |
| `16-tidyverse-code-review.tex` | archive; Appendix I Section 10 | planning record |
| `17-simplification-plan.tex` | archive; Appendix I Sections 10, 11 | origin credit for E9 and RM-ANOVA; a suspect reference |
| `18-response-parameter-sensitivity.tex` | archive; Appendix I Section 3 | figures missing; result kept as historical |
| `19-mean-moderation-implementation-notes.md`, `.pdf`, `.tex` | archive all three; Appendix I Section 5; `ch07` pointer | unique content placed in Appendix I |
| `20-annotated-bibliography.tex` | archive; Appendix G; Appendix I Section 11 | misattributions; corrected list in `docs/23` Section 17 |
| `21-architecture-comparison-summary.md`, `.pdf` | archive; `ch07` (rationale) | superseded; 50-replicate tables unreproducible |
| `22-journal-target-recommendation.md`, `.pdf` | remove from series: the paper 01 folder, `planning/` | planning, not methods |
| `22-tidyverse-development-parity-workflow.md`, `.pdf` | remove from series: a README in `analysis/scripts/parity/`; Appendix J pointer | current developer guide, not methods; frees the duplicate number |
| `23-documentation-compendium.md`, `.pdf` | archive; salvage 3.1 to `ch03`, 8.1 to 8.2 to `ch07`, 14 to Appendix D, 17 to Appendix G | every section has a fuller current home |
| `24-component-decomposition-pedagogy.md`, `.pdf` | `ch05` (components) and `ch12` (identification) | current pedagogy; analysis advice rewritten to point to `ch22` |
| `25-gompertz-model-evaluation.md`, `.pdf` | `ch05` | literature, alternatives and ranked impact kept; numerics from CM |
| `26-rm-anova-interaction-test.md`, `.pdf` | Appendix H | current derivation, reference character |
| `27-cycle-period-design-sweep-plan.md`, `.pdf` | remove from series: the paper 05 folder, `planning/`; Section 4 to `ch16` | planning record; corrected information result kept |
| `28-hendrickson-hybrid-design-and-decomposition.md`, `.pdf` | archive after porting its 2026-10-08 corrections to `ch12` | superseded by `docs/31` |
| `29-nof1-precision-medicine-lexicon.md`, `.pdf` | Appendix E | reference; fix $\tau^2$ |
| `30-terminology-consistency-plan.md`, `.pdf` | Appendix F (extended); execution record to Appendix B | prose authority beside NOTATION |
| `31-component-decomposition-and-reductions.md`, `.pdf` | split: Parts II to V to `ch12`; Parts I and VI removed (pointer to `ch05`); then archive | executes `plan-doc31-revision` Phase 0 |
| `32-why-covar-departs-from-orig.md`, `.pdf` | `ch10` (framing), `ch09` (numbers labeled configuration F); then archive | survey b recommended merging into `docs/54`; `docs/54` largely superseded |
| `33-paper01-dgp-architecture-review.md`, `.pdf` | remove from series: the paper 01 folder, `revisions/`; findings 2.3, 2.4 to `ch07`, `ch18`; Section 2.5 to `ch02` | manuscript review, not a methods chapter |
| `34-cbm-feasible-range.md`, `.pdf` | `ch09` (formula, schedule dependence), `ch16` (Section 5 levers); Sections 6.1, 7 archive; 6.2 to Appendix D | partly superseded at $\rho_{cs} = 0.7$ |
| `35-separable-response-covariance.md`, `.pdf` | `ch09` spine; Sections 7, 9 to Appendix D | current |
| `36-hendrickson-pd-and-carryover.md`, `.pdf` | split: Sections 1 to 5 to `ch04`; 6.1 to 6.6, 6.8 to `ch10`; 6.7 and Appendix A.7 to `ch09`; 8 to Appendix C; 9 to Appendix D | current; drop "Working draft" |
| `37-carryover-power-closed-form.md`, `.pdf` | `ch15` spine; Appendix A to `ch28` | current |
| `38-pharmacodynamic-carryover-model.md`, `.pdf` | `ch14` spine | proposal; symbols renamed to `docs/44`'s |
| `39-nof1-design-optimization.md`, `.pdf` | `ch16` spine | current within its construct |
| `40-hendrickson-hybrid-design-primer.md`, `.pdf` | split: Sections 2 to 3 to `ch02`; 3 to 6 and 8 to `ch03`, `ch05`; 9 to `ch17`; 11 to Appendix C | primer; duplicated by `docs/42` |
| `41-baseline-as-covariate.md`, `.pdf` | `ch17` (Propositions 1 to 5); Proposition 4 to `ch15`; recommendation to Appendix B | partly superseded |
| `42-flexible-time-adjustment.md`, `.pdf` | `ch17` (time adjustment, merged with `docs/40` Section 9) | duplicate tables removed |
| `43-analysis-approaches-catalog.md`, `.pdf` | Appendix A (regenerated) | reference, incomplete since `docs/56` |
| `44-dgp-roadmap.md`, `.pdf` | `ch10` (step table); Steps 7 to 8 to `ch14`; status to Appendix D | planning, current as status board |
| `45-correctly-sized-hybrid-analysis.md`, `.pdf` | `ch20` spine; Section 6 to `ch19`; 7 to `ch17`; 8 to Appendix C; 9 to Appendix B | covariance-study record |
| `46-strawman-stress-test-protocol.md`, `.pdf` | `ch21` annex (protocol of record) | historical protocol; Amendment 2 item 4 cut to a pointer |
| `47-stress-test-results.md`, `.pdf` | `ch21` spine | current record of the stored run |
| `48-shared-biomarker-effect-and-stress-test-issues.md`, `.pdf` | Sections 1 to 3 to `ch19`; 4 to `ch21`; 5 dropped | partly superseded; Section 4 table corrected (CR-20) |
| `49-predictive-prognostic-biomarker-expectancy-term.md`, `.pdf` | Sections 2 to 4 to `ch19`; 6 to `ch22`; 5 and 7 to Appendix B | recommendation superseded |
| `50-dgp-and-analysis-status-brief.md`, `.pdf` | archive; regenerated as `fm1` | superseded brief |
| `51-ten-features-of-an-nof1-dgp.md`, `.pdf` | `ch11` (renamed "seventeen features"); Section 20 to `ch10`; 21 to paper 01 | current framework |
| `52-whitepaper-review-2026-10-07.md`, `.pdf` | archive; Section 2 to Appendix B; Section 4 to Appendix D | review log |
| `53-carryover-recursion-in-published-code.md`, `.pdf` | `ch26` (kept standalone, DQ-3) | current note for the original authors |
| `54-covar-vs-orig-dgp.md`, `.pdf` | archive; Section 3.6.2 to `ch27`; Section 5 to `ch09` if CM lacks it | largely superseded; assessment of the ratio rule overtaken by `docs/55` |
| `55-step-coupling-in-published-code.md`, `.pdf` | `ch27` (kept standalone, DQ-3) | current note for the original authors |
| `56-recommended-analysis-expectancy-cr2.md`, `.pdf` | `ch22` spine; Section 6 to `ch23` | current, internally stale in Section 9 |
| `57-data-chosen-carryover-half-life.md`, `.pdf` | `ch23` spine; Sections 6 and 7.3 to an annex | current |
| `58-main-effect-testing-and-carryover-half-life.md`, `.pdf` | `ch24` spine | current |
| `59-whitepaper-refactoring-plan.md`, `.pdf` | archive at release; Section 4 becomes `fm0`'s old-to-new map | planning |
| `plan-architecture-c.md` | remove from series: the paper 11 folder, `planning/` | superseded as to placement |
| `plan-coherent-analysis-comparison.md` | remove from series: the paper 02 folder, `revisions/`; Q1 to Q3 to `ch17` | planning; the forbidden word on line 20 left in the moved copy only |
| `plan-doc31-revision.md`, `.pdf` | archive | executed by this refactor |
| `notes.md` | `docs/reviews/` with a header naming its source and date; Appendix I Section 8 | reviewer comments |
| `reorg-plan-2026-04-22.md` | archive | executed in a different form |
| `ZZCOLLAB_USER_GUIDE.md` | remove from series: `docs/tooling/` | framework guide, no project content |
| `pub_review_whitepaper_2026-08-16.md`, `pub_review_remediation_2026-08-20.md` | `docs/reviews/` | manuscript reviews of the August campaign |
| `Architecting_Precision_Medicine.pdf` | `docs/archive/`, labeled unreviewed (DQ-4) | 12.2 MB AI-generated deck, no text layer, cited nowhere, unreviewed |
| `pmsimstats-preamble.tex` | archive, beside the archived `.tex` set | needed to compile the archived `.tex` |
| `figures/` (44 files) | stays; new figures named `chNN-figK-*.png`; `40-fig2` and `40-fig3` renamed to print order | archived documents keep their paths |
| `archive/` (4 files) | stays, for non-whitepaper historical material; receives the slide deck (DQ-4) | already historical |
| `references/` (13 files) | stays; `nof1-simulation.bib` becomes the series bibliography source | `pmsimstats-nof1-whitepaper.md` there was not surveyed |
| `reviews/` (1 file) | stays; receives `docs/02` PDF, `notes.md`, the `pub_review` files | review records |
| `.claude/` | unchanged | settings, not documentation |

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | all | Phase 2 content move: reading paths, chapter status table, old-to-new map | Phase 2 skeleton fill |
| 2026-10-15 | About this series, Reading paths, Notation cross-reference | Phase 3 GP-01: purpose section, path annotations, additional paths, notation cross-reference | GP-01 fill |
| 2026-10-15 | Chapter status | Update status column from "stub" for all to actual state after Phases 3-4: 24 substantial, 16 partial, 1 stub | Post-Phase-4 survey |
