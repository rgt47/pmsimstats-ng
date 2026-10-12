# Index of Legacy Whitepapers: Old Number to New Home
*2026-10-11 18:04 PDT*

**Author.** pmsimstats team

**Purpose.** This index keeps every old document number resolvable
while the whitepapers are rebuilt as one series under `docs/series/`
(plan: `docs/59-whitepaper-refactoring-plan.md`, Sections 4, 9.2 and
9.3). A citation `docs/NN` keeps meaning "legacy document NN". Each
row gives the new home of the document, or of each section where the
plan splits it, using the identifiers of `docs/59` Section 3 (`fm0`
to `fm2`, `ch01` to `ch28`, Appendices A to J, Volume II for the
mathematics compendium); the file's location on 2026-10-11 after
Phase 0; and its status. The series chapters do not exist yet; a new
home names where the content is to go. The repository state before
Phase 0 is tagged `whitepapers-pre-series`.

Status values: **archived** (moved to `docs/legacy/` in Phase 0 with a
one-line banner; frozen); **removed from series** (moved out of
`docs/` to the named place in Phase 0; not part of the read-through);
**still in docs/** (stays at its old path until its content is moved
in Phase 2 and the file is retired in Phase 5); **deleted** (removed
in Phase 0 as a stale duplicate; recoverable from the tag).

## 1. Numbered documents

| No. | Old file(s) | Section(s) | New home | Current location | Status |
|----|----------------------------|----------|--------------------------|----------------------------|----------------|
| 00 | `00-documentation-index.tex` | all | `fm0` (generated index) | `docs/legacy/` | archived |
| 01 | `01-codebase-overview.tex` | all | archive only (`CLAUDE.md` and package documentation) | `docs/legacy/` | archived |
| 02 | `02-dgp-mean-moderation-vs-mvn-rchcoments.pdf` | annotations (pages 2, 3, 5, 7, 13) | Appendix I Section 8 | `docs/reviews/review-comments-2026-04-29.pdf` | removed from series |
| 03 | `03-audit-and-revision-report.tex` | all | Appendix I Sections 1, 9; commit attribution in `ch03` | `docs/legacy/` | archived |
| 04 | `04-revised-power-analysis.tex` | all | Appendix I Sections 1, 3 | `docs/legacy/` | archived |
| 05 | (no document) | | | | number unused |
| 06 | `06-ar1-residual-correlation.tex` | all | Appendix I Section 4 | `docs/legacy/` | archived |
| 07 | `07-positive-definiteness-failures.tex` | all | Appendix I Section 1 (counts and condition numbers) | `docs/legacy/` | archived |
| 07a | `07a-pd-constraints-theory.tex`, `.pdf` | all | archive only (nothing carried) | `docs/legacy/` | archived |
| 08 | `08-biomarker-correlation-decay.tex` | all | `ch08` (credit); Appendix I Section 7 | `docs/legacy/` | archived |
| 09 | `09-carryover-correlation-artifact.tex` | all | Appendix I Section 7 | `docs/legacy/` | archived |
| 09a | `09a-carryover-bias-variance.md`, `.pdf`, `.tex` | all | Appendix I Section 6 | `docs/legacy/` | archived |
| 10 to 12 | (no documents) | | | | numbers unused |
| 13 | `13-figure4-walkthrough.tex` | all | Appendix I Section 2 | `docs/legacy/` | archived |
| 14 | `14-trial-analysis-implementation-guide.tex` | all | Appendix I Section 9; `ch25` | `docs/legacy/` | archived |
| 15 | `15-tidyverse-alignment-plan.tex` | all | Appendix I Section 10 | `docs/legacy/` | archived |
| 16 | `16-tidyverse-code-review.tex` | all | Appendix I Section 10 | `docs/legacy/` | archived |
| 17 | `17-simplification-plan.tex` | all | Appendix I Sections 10, 11 | `docs/legacy/` | archived |
| 18 | `18-response-parameter-sensitivity.tex` | all | Appendix I Section 3 | `docs/legacy/` | archived |
| 19 | `19-mean-moderation-implementation-notes.md`, `.pdf`, `.tex` | all | Appendix I Section 5; pointer from `ch07` | `docs/legacy/` | archived |
| 20 | `20-annotated-bibliography.tex` | all | Appendix G; Appendix I Section 11 | `docs/legacy/` | archived |
| 21 | `21-architecture-comparison-summary.md`, `.pdf` | all | `ch07` (biological rationale) | `docs/legacy/` | archived |
| 22 (journal targets) | `22-journal-target-recommendation.md`, `.pdf` | all | none (planning) | `analysis/report/01-dgp-mean-moderation-vs-mvn/planning/` | removed from series |
| 22 (parity workflow) | `22-tidyverse-development-parity-workflow.md`, `.pdf` | all | Appendix J (pointer) | `analysis/scripts/parity/README.md` (PDF beside it under its old name) | removed from series |
| 23 | `23-documentation-compendium.md`, `.pdf` | 3.1 | `ch03` | `docs/legacy/` | archived |
| 23 | | 8.1 to 8.2 | `ch07` | `docs/legacy/` | archived |
| 23 | | 14 | Appendix D | `docs/legacy/` | archived |
| 23 | | 17 | Appendix G | `docs/legacy/` | archived |
| 23 | | other sections | archive only (fuller current homes elsewhere) | `docs/legacy/` | archived |
| 24 | `24-component-decomposition-pedagogy.md`, `.pdf` | components | `ch05` | `docs/` | still in docs/ until Phase 2/5 |
| 24 | | identification; analysis advice | `ch12`; advice points to `ch22` | `docs/` | still in docs/ until Phase 2/5 |
| 25 | `25-gompertz-model-evaluation.md`, `.pdf` | all | `ch05` | `docs/` | still in docs/ until Phase 2/5 |
| 26 | `26-rm-anova-interaction-test.md`, `.pdf` | all | Appendix H | `docs/` | still in docs/ until Phase 2/5 |
| 27 | `27-cycle-period-design-sweep-plan.md`, `.pdf` | all | planning record; Section 4 to `ch16` | `analysis/report/05-nof1-design-sensitivity/planning/` | removed from series |
| 28 | `28-hendrickson-hybrid-design-and-decomposition.md`, `.pdf` | 2026-10-08 corrections | `ch12`; then archive | `docs/` | still in docs/ until Phase 2/5 |
| 29 | `29-nof1-precision-medicine-lexicon.md`, `.pdf` | all | Appendix E | `docs/` | still in docs/ until Phase 2/5 |
| 29 | `29-nof1-precision-medicine-lexicon_2026-06-18_1719.pdf` | dated duplicate | none | (removed) | deleted |
| 30 | `30-terminology-consistency-plan.md`, `.pdf` | all; execution record | Appendix F; execution record to Appendix B | `docs/` | still in docs/ until Phase 2/5 |
| 30 | `30-terminology-consistency-plan_2026-06-18_0834.pdf` | dated duplicate | none | (removed) | deleted |
| 31 | `31-component-decomposition-and-reductions.md`, `.pdf` | Parts II to V | `ch12` | `docs/` | still in docs/ until Phase 2/5 |
| 31 | | Parts I, VI | removed (pointer to `ch05`); then archive | `docs/` | still in docs/ until Phase 2/5 |
| 32 | `32-why-covar-departs-from-orig.md`, `.pdf` | framing | `ch10` | `docs/` | still in docs/ until Phase 2/5 |
| 32 | | numbers | `ch09` (labeled configuration F); then archive | `docs/` | still in docs/ until Phase 2/5 |
| 33 | `33-paper01-dgp-architecture-review.md`, `.pdf` | 2.3, 2.4, 2.5 | `ch07`, `ch18`, `ch02` | `analysis/report/01-dgp-mean-moderation-vs-mvn/revisions/` | removed from series |
| 34 | `34-cbm-feasible-range.md`, `.pdf` | formula, schedule dependence | `ch09` | `docs/` | still in docs/ until Phase 2/5 |
| 34 | | 5 (levers) | `ch16` | `docs/` | still in docs/ until Phase 2/5 |
| 34 | | 6.1, 7 | archive only | `docs/` | still in docs/ until Phase 2/5 |
| 34 | | 6.2 | Appendix D | `docs/` | still in docs/ until Phase 2/5 |
| 35 | `35-separable-response-covariance.md`, `.pdf` | body | `ch09` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 35 | | 7, 9 | Appendix D | `docs/` | still in docs/ until Phase 2/5 |
| 36 | `36-hendrickson-pd-and-carryover.md`, `.pdf` | 1 to 5 | `ch04` | `docs/` | still in docs/ until Phase 2/5 |
| 36 | | 6.1 to 6.6, 6.8 | `ch10` | `docs/` | still in docs/ until Phase 2/5 |
| 36 | | 6.7, Appendix A.7 | `ch09` | `docs/` | still in docs/ until Phase 2/5 |
| 36 | | 8 | Appendix C | `docs/` | still in docs/ until Phase 2/5 |
| 36 | | 9 | Appendix D | `docs/` | still in docs/ until Phase 2/5 |
| 37 | `37-carryover-power-closed-form.md`, `.pdf` | body | `ch15` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 37 | | Appendix A | `ch28` | `docs/` | still in docs/ until Phase 2/5 |
| 38 | `38-pharmacodynamic-carryover-model.md`, `.pdf` | all | `ch14` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 39 | `39-nof1-design-optimization.md`, `.pdf` | all | `ch16` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 40 | `40-hendrickson-hybrid-design-primer.md`, `.pdf` | 2 to 3 | `ch02` | `docs/` | still in docs/ until Phase 2/5 |
| 40 | | 3 to 6, 8 | `ch03`, `ch05` | `docs/` | still in docs/ until Phase 2/5 |
| 40 | | 9 | `ch17` | `docs/` | still in docs/ until Phase 2/5 |
| 40 | | 11 | Appendix C | `docs/` | still in docs/ until Phase 2/5 |
| 41 | `41-baseline-as-covariate.md`, `.pdf` | Propositions 1 to 5 | `ch17` | `docs/` | still in docs/ until Phase 2/5 |
| 41 | | Proposition 4 | `ch15` | `docs/` | still in docs/ until Phase 2/5 |
| 41 | | recommendation | Appendix B | `docs/` | still in docs/ until Phase 2/5 |
| 42 | `42-flexible-time-adjustment.md`, `.pdf` | all | `ch17` (merged with 40 Section 9) | `docs/` | still in docs/ until Phase 2/5 |
| 43 | `43-analysis-approaches-catalog.md`, `.pdf` | all | Appendix A (regenerated) | `docs/` | still in docs/ until Phase 2/5 |
| 44 | `44-dgp-roadmap.md`, `.pdf` | step table | `ch10` | `docs/` | still in docs/ until Phase 2/5 |
| 44 | | Steps 7 to 8 | `ch14` | `docs/` | still in docs/ until Phase 2/5 |
| 44 | | status | Appendix D | `docs/` | still in docs/ until Phase 2/5 |
| 45 | `45-correctly-sized-hybrid-analysis.md`, `.pdf` | body | `ch20` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 45 | | 6 | `ch19` | `docs/` | still in docs/ until Phase 2/5 |
| 45 | | 7 | `ch17` | `docs/` | still in docs/ until Phase 2/5 |
| 45 | | 8 | Appendix C | `docs/` | still in docs/ until Phase 2/5 |
| 45 | | 9 | Appendix B | `docs/` | still in docs/ until Phase 2/5 |
| 46 | `46-strawman-stress-test-protocol.md`, `.pdf` | all | `ch21` annex (protocol of record) | `docs/` | still in docs/ until Phase 2/5 |
| 47 | `47-stress-test-results.md`, `.pdf` | all | `ch21` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 48 | `48-shared-biomarker-effect-and-stress-test-issues.md`, `.pdf` | 1 to 3 | `ch19` | `docs/` | still in docs/ until Phase 2/5 |
| 48 | | 4 | `ch21` | `docs/` | still in docs/ until Phase 2/5 |
| 48 | | 5 | dropped | `docs/` | still in docs/ until Phase 2/5 |
| 49 | `49-predictive-prognostic-biomarker-expectancy-term.md`, `.pdf` | 2 to 4 | `ch19` | `docs/` | still in docs/ until Phase 2/5 |
| 49 | | 6 | `ch22` | `docs/` | still in docs/ until Phase 2/5 |
| 49 | | 5, 7 | Appendix B | `docs/` | still in docs/ until Phase 2/5 |
| 50 | `50-dgp-and-analysis-status-brief.md`, `.pdf` | all | `fm1` (regenerated) | `docs/legacy/` | archived |
| 51 | `51-ten-features-of-an-nof1-dgp.md`, `.pdf` | body | `ch11` ("seventeen features") | `docs/` | still in docs/ until Phase 2/5 |
| 51 | | 20 | `ch10` | `docs/` | still in docs/ until Phase 2/5 |
| 51 | | 21 | paper 01 | `docs/` | still in docs/ until Phase 2/5 |
| 52 | `52-whitepaper-review-2026-10-07.md`, `.pdf` | 2 | Appendix B | `docs/legacy/` | archived |
| 52 | | 4 | Appendix D | `docs/legacy/` | archived |
| 53 | `53-carryover-recursion-in-published-code.md`, `.pdf` | all | `ch26` (standalone, DQ-3) | `docs/` | still in docs/ until Phase 2/5 |
| 54 | `54-covar-vs-orig-dgp.md`, `.pdf` | 3.6.2 | `ch27` | `docs/legacy/` | archived |
| 54 | | 5 | `ch09` if the compendium lacks it | `docs/legacy/` | archived |
| 54 | | other sections | archive only (largely superseded; ratio-rule assessment overtaken by 55) | `docs/legacy/` | archived |
| 55 | `55-step-coupling-in-published-code.md`, `.pdf` | all | `ch27` (standalone, DQ-3) | `docs/` | still in docs/ until Phase 2/5 |
| 56 | `56-recommended-analysis-expectancy-cr2.md`, `.pdf` | body | `ch22` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 56 | | 6 | `ch23` | `docs/` | still in docs/ until Phase 2/5 |
| 57 | `57-data-chosen-carryover-half-life.md`, `.pdf` | body | `ch23` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 57 | | 6, 7.3 | `ch23` annex | `docs/` | still in docs/ until Phase 2/5 |
| 58 | `58-main-effect-testing-and-carryover-half-life.md`, `.pdf` | all | `ch24` (spine) | `docs/` | still in docs/ until Phase 2/5 |
| 59 | `59-whitepaper-refactoring-plan.md`, `.pdf` | all; Section 4 | archive at release (Phase 5); Section 4 becomes the old-to-new map of `fm0` | `docs/` | still in docs/ until Phase 5 |

Table: Old document numbers and their new homes (from `docs/59`
Section 4)

## 2. Unnumbered files

| Old file(s) | New home | Current location | Status |
|----------------------------|--------------------------|----------------------------|----------------|
| `plan-doc31-revision.md`, `.pdf` | archive only (executed by the refactor) | `docs/legacy/` | archived |
| `reorg-plan-2026-04-22.md` | archive only (executed in a different form) | `docs/legacy/` | archived |
| `pmsimstats-preamble.tex` | archive only (needed to compile the archived `.tex` set) | `docs/legacy/` | archived |
| `plan-architecture-c.md` | none (planning; paper 11) | `analysis/report/11-combined-dgp-architecture/planning/` | removed from series |
| `plan-coherent-analysis-comparison.md` | Q1 to Q3 to `ch17` | `analysis/report/02-carryover-sensitivity/revisions/` | removed from series |
| `notes.md` | Appendix I Section 8 | `docs/reviews/notes.md` (header added) | removed from series |
| `ZZCOLLAB_USER_GUIDE.md` | none (framework guide) | `docs/tooling/` | removed from series |
| `pub_review_whitepaper_2026-08-16.md`, `pub_review_whitepaper_2026-08-20.md`, `pub_review_remediation_2026-08-20.md` | none (manuscript reviews of the August campaign) | `docs/reviews/` | removed from series |
| `Architecting_Precision_Medicine.pdf` | none (unreviewed AI-generated deck, DQ-4) | `docs/archive/` | removed from series |
| `figures/`, `references/`, `reviews/`, `archive/`, `.claude/` | unchanged in Phase 0 | `docs/` | unchanged |

Table: Unnumbered files of `docs/` and their Phase 0 disposition

## 3. Notes on use

- **Archived files are frozen.** Apart from the banner line (after
  the title and timestamp of a `.md` file; a `%` comment on the first
  line of a `.tex` file), their content is as at tag
  `whitepapers-pre-series`. PDFs were moved without modification.
- **Compiling an archived `.tex`.** The preamble moved with the `.tex`
  set, so `\input{pmsimstats-preamble}` still resolves from
  `docs/legacy/`. Its `\graphicspath{{figures/}}` is relative, so
  `docs/legacy/04-revised-power-analysis.tex` would need the figures
  path adjusted to compile from its new location; `docs/18` already
  could not render (its figures under `analysis/figure5/` are gone).
- **Citations elsewhere.** Scripts, manuscripts, the compendium,
  `NOTATION.md` and `CLAUDE.md` still cite `docs/NN` by old path; per
  `docs/59` Section 9.4 they are updated in Phase 5, and until then
  this index resolves them.
