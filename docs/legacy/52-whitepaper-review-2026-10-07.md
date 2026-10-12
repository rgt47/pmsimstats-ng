# Review of the pmsimstats-ng Whitepapers: Notation, Accuracy and Conventions {.unlisted .unnumbered}
*2026-10-07 19:03 PDT*

*Archived 2026-10-11; current treatment: Section 2 in Appendix B (decision history), Section 4 in Appendix D (open-items register); see `docs/legacy/INDEX.md`.*

**Author.** pmsimstats team

**Purpose.** This document records a deep review of the 34 numbered
whitepapers in `docs/` (09a, 19 to 51, and
`whitepaper-covar-vs-orig-dgp.md`, renumbered `docs/54` on 2026-10-09) for three things: consistency with
`analysis/report/NOTATION.md`, accuracy, and the project's documentation
conventions. Problems found were fixed in place; superseded statements
received dated notes rather than silent rewrites. It lists what changed
in each document and the items left open.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Beyond symbol names quoted from the reviewed documents, this document
uses only the following symbols of `analysis/report/NOTATION.md`
(index added 2026-10-09).

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $\eta_t$, $\eta$ | design expectancy at visit $t$ (1 open label, 0.5 blinded, 0 at baseline); $\eta$ as the 2026-10-07 `NOTATION.md` wrote it (Section 2) | canonical |
| $e_t$ | former symbol for $\eta_t$, quoted only to record the rename | retired |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $\rho_{cs}$ | common correlation of the compound-symmetry construct (row added 2026-10-09) | canonical |
| $N$ | total sample size across randomization paths | canonical |

Table: Notation index of symbols used in this document

**Glossary.**

- **Dated note.** A short paragraph, dated 2026-10-07, that records a
  correction or a superseding result and points to the newer document.
- **Linter.** `tools/notation-lint.pl`, which checks a few mechanical
  rules of `NOTATION.md` (bare `Dbc` and `Sx`, legacy specification
  labels, day-denominated half-lives, per-path values of $N$).

## 1. Method

Each document was read end to end against one rubric:

- **Notation.** Every symbol checked against `NOTATION.md`; colliding or
  non-canonical symbols renamed (not merely declared as departures); a
  notation index and glossary placed at the front; the linter run
  before and after.
- **Accuracy.** Derivations re-derived; numbers in prose, tables and
  captions cross-checked; claims about code checked against the current
  source and, for the published process, commit `58b32a9`; numerical
  claims checked against the cited data files and scripts; references
  checked against PubMed or Crossref.
- **Conventions.** Header block (title, timestamp, Author, Purpose),
  contents block and captions, US English, no em dashes, evidence
  labels, no hand-written render stamps, single titles, figure captions
  not duplicated, resolvable citations.

Every document was rendered through `tools/render.sh` from `docs/` and
checked with `pdftotext` and `pdfimages`.

**Final state (verified).** The linter is clean on all 34 documents.
No em dashes, no hand-written "Rendered on" lines and no use of the word
"epistemic" remain. Every table carries a caption, and every figure is
embedded and listed (docs/50, the two-page brief, has no contents block
by design).

## 2. Changes that apply across documents

- **Expectancy symbol.** `NOTATION.md` defines $\eta$ as the canonical
  expectancy covariate. The design expectancy written $e_t$ (1 open
  label, 0.5 blinded, 0 at baseline; code `e`, regressor `De`) is now
  $\eta_t$ in every document. In `docs/46` the rename is recorded in
  Amendment 2 as typographical. (Notation, 2026-10-09: `NOTATION.md`
  now defines $\eta_t$ itself as the design expectancy and reserves
  $\tilde\eta_{it}$ for a measured placebo-belief covariate.)
- **Ceiling values.** Statements no longer mix serial correlations: the
  compound-symmetry ceiling is 0.256 without carryover; the separable
  AR(1) ceiling is 0.440 at $\rho = 0.8$ and 0.474 (Hybrid 0.480) at
  $\rho = 0.7$. The published process uses compound symmetry with
  $\rho_{cs} = 0.8$ (`docs/34`). (Notation, 2026-10-09: the
  compound-symmetry correlation was written $\rho$.)
- **Stress-test outcome.** Documents that describe the stress test now
  give 41 cells, 28,000 trials and the 14-analysis decision rerun
  (`docs/46`, Amendment 2): F1 selected, the phase-variance versions G1
  and G2 effectively equivalent, F1A failing size.
- **Superseded recommendations.** The unstructured Kenward-Roger
  recommendation of `docs/45` and its repetitions are marked superseded
  by `docs/46` and `docs/49`. The flexible-time size excess is stated as
  real under model-based and Kenward-Roger inference and absent under
  CR2 (`docs/42`).
- **References.** Every reference list was checked. Several documents
  had copied invented author lists from `docs/references/nof1-simulation.bib`
  (corrected separately the same day); these are now attributed to the
  real papers (for example Shi and Ye 2024, Chen and Chen 2014,
  Hawksworth et al. 2024, Chatters et al. 2024, Senn 2019, Lim and In
  2021, Schmid and Yang 2022). DOIs were added where they exist.

## 3. Per-document summary

Counts are approximate; each count is a group of related edits.

Table: Review outcome by document

| Document | Notation | Accuracy | Main accuracy findings |
|------------------|--------|--------|------------------------------------------------|
| 09a | 6 | 13 | Headline "drop then recovery" not supported (within 2 MCSE); model misdescribed; three of six references wrong; conclusion 4 withdrawn |
| 19 | 6 | 2 | Noise SD 5 corrected to 8; shift 650 corrected to 558 |
| 21 | 10 | 7 | Hybrid loss 68% corrected to 39%; range "40-60%" corrected to 20-61%; tables not reproducible |
| 22 (journal) | 3 | 9 | Six misattributed citations corrected |
| 22 (parity) | 1 | 6 | Duplicate title removed; script behavior corrected against the code |
| 23 | 25 | 17 | Multiplicative/additive equivalence, AR(1) residual formula, step-rule row label, seven references |
| 24 | 8 | 13 | Gompertz parameter meaning reversed; Hybrid design misdescribed; variance identity |
| 25 | 8 | 10 | Parameter table (TV row) wrong; heterogeneity claim withdrawn; superseded by paper 07 |
| 26 | 12 | 9 | Sphericity conflated with compound symmetry; misattributed Type I figure |
| 27 | 8 | 9 | Tier 1 cell count 264 corrected to 324; E2 decision premise false; duplicate title |
| 28 | 10 | 3 | Superseded by docs/31 (light review); crossover row and power pair corrected |
| 29 | 0 | 2 | US English; two definitions corrected |
| 30 | 3 | 4 | Stale process statements updated |
| 31 | 15 | 22 | Gompertz section; covariance dimension; paper 06 sample sizes; useDE claim; references |
| 32 | 4 | 8 | Commit dating; "one-week" collapse corrected to any positive half-life; ceiling decline with visits |
| 33 | 6 | 3 | One value 20.0% corrected to 19.9%; CS versus AR(1) reversal attributed to step coupling |
| 34 | 2 | 8 | Published-process numbers used $\rho = 0.7$, not 0.8; binding path A, not C; second silent repair |
| 35 | 9 | 11 | Ceiling formula term; stress-test description (26,500 to 28,000; not nominal everywhere) |
| 36 | 8 | 7 | Five configurations and 261,000 fits corrected to six and 309,000 |
| 37 | 4 | 7 | Six small numerical misquotes corrected against the data files; derivations verified |
| 38 | 12 | 4 | Sign of a correlation corrected; Hobbs 1978 subjects |
| 39 | 8 | 8 | 48 conditions corrected to 96; noise ratio a fifth, not a quarter |
| 40 | 1 | 5 | Figure references swapped; time on drug restarts at zero |
| 41 | 6 | 6 | Power floor range; recommendation superseded |
| 42 | 3 | 8 | MCSE multiples; flexible-time correction applied |
| 43 | 2 | 20 | 14-analysis stress-test table added; status entries updated |
| 44 | 3 | 9 | "Carryover raises every ceiling" false at $\rho = 0.5$; step statuses dated |
| 45 | 6 | 7 | Section 4.4 correction; superseded primary recommendation |
| 46 | 5 | 4 | Amendment 2 added; original protocol text unchanged in substance |
| 48 | 3 | 7 | Pending results filled in; maximin table scope stated |
| 49 | 9 | 8 | Pending results filled in; cancellation claim qualified |
| 50 | 1 | 5 | Analysis table and recommendation updated; still two pages |
| 51 | 18 | 7 | Package carryover is anchored, not compounding; Feature 6 results; dropout size wording |
| 54 (covar vs orig) | 100 | 13 | "orig" means the 2024 copy, not `58b32a9`; AR(1)-alone result depends on an uncoupled first visit |

## 4. Open items

**Code.**

- `R/carryover_analysis.R`: `threshold_bm` and `predicted_effect` added
  `bm_mean`, but `lme_analysis()` fits the uncentered biomarker, so the
  threshold was offset. Fixed 2026-10-08: both are now computed on the
  fitted scale, and `inst/tinytest/test_carryover_analysis.R` checks the
  predicted effect against the model's own contrast (verified; with a
  biomarker of mean 98 the old code was off by 7.0 points and put the
  threshold at 172.6 instead of 74.5).
- `R/generateData.R:391-393`: the Cholesky fallback still repairs a
  non-positive-definite matrix silently.
- `analysis/scripts/carryover-sensitivity/simulation-core.R`: a comment
  says the TV parameters follow the published defaults; they do not.

**Data and reproducibility.**

- `docs/21`: the tables cannot be reproduced; the driver's output file
  is missing. Re-run the driver at a realistic replicate count, or
  retire the tables in favor of paper 01.
- `analysis/scripts/parity/`: the committed `parity_report.rds` is a
  quick-mode report (16 rows) while `parity_baseline.rds` is default
  mode (144 rows), so `parity-diff.R` compares unlike reports.
- `docs/34`: $\rho = 0.8$ is missing from `rho.csv` and Figure 2; the
  added table row was computed independently. Figure 10's subtitle
  ("binding paths") is wrong for Hybrid and needs the figure script.
- `docs/40`: the figure files `40-fig2` and `40-fig3` are numbered
  opposite to the order in which they print.
- Manuscript 02: whether its crossover cells were rerun with the
  corrected exposure coding is not established.

**Documents.**

- `docs/47`, the stress-test results report referred to by `docs/46`
  and `docs/48`, did not exist; it was written on 2026-10-08 from the
  stored decision outputs.
- The 14-analysis rerun ranks G1 (0.308) and G2 (0.312) ahead of F1
  (0.358) over all 26 alternative cells, prognostic cells included. The
  working protocol does not say whether those cells count; F1 remains
  selected.
- `docs/09a`: the source of its Section 1.1 table is unknown, and its
  robustness recommendation cites simulations it does not report.
- `docs/24`, `docs/31`: three claims have no source and are flagged in
  the text (residual expectancy of 30-60%, a chronic-fatigue example,
  convergence below 90%). Corrections made in `docs/31` to Parts I-III
  should be carried into `docs/24`.
- `docs/28`, `docs/31` and paper 06 may still quote stale per-path
  sample sizes for Studies A and C.
- `docs/20-annotated-bibliography.tex` probably repeats the reference
  errors corrected in `docs/23`; not checked.
- The AHRQ user-guide chapter authors rest on web-search reports; the
  AHRQ site blocks automated access.

**Not done.** No Monte Carlo cell was re-simulated; numerical checks
used stored outputs. The `.tex` companions of the documents were not
reviewed.

## 5. Reproducibility

The per-document review logs and check scripts were written to the
session scratchpad and are not part of the repository. The final checks
can be repeated with:

```
perl tools/notation-lint.pl docs/[0-9]*.md
```
