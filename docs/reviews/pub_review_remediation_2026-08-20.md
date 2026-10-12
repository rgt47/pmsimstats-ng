# pub_review remediation log
*2026-08-20 13:16 PDT*

Remediation pass against
`docs/pub_review_whitepaper_2026-08-16.md`, as instructed. Prioritized
series 04, 05, and 09 (the correctness-tier defects in headline
results per that whitepaper), plus cheap, high-value correctness and
acceptance items elsewhere. Time-boxed aggressively; did not attempt
all eleven series.

## 0. Important context discovered before starting

Before making any edits, `git log` showed that an earlier remediation
pass already occurred, bundled into commit `8817308` ("Auto-backup:
2026-08-16 17:37:00"), roughly 90 minutes after the 2026-08-16 16:10
PDT whitepaper review. That commit touched `DESCRIPTION`,
`vignettes/01-simulate-trial-data.Rmd`,
`vignettes/02-visualize-power-results.Rmd`, and the `report.Rmd`/
`cover-letter.md` files for series 04, 05, 06, 07, and 09, and already
resolved most of Section 4(a)'s correctness checklist. That prior
pass's own remediation log
(`docs/pub_review_remediation_2026-08-16.md`, referenced by the
revised `07/cover-letter.md`) was never actually written or committed
- a genuine gap, noted under "New issues found" below. This session's
work is therefore (1) verification of that prior pass's fixes, (2) a
handful of new fixes the prior pass missed or left incomplete, and
(3) this log, which now exists where the prior one did not.

A newer whitepaper, `docs/pub_review_whitepaper_2026-08-20.md`, also
exists in this repository. Per the remediation instructions, this
session worked from the 2026-08-16 whitepaper as directed. The 08-20
whitepaper's own method note states it did not re-verify series
01/03-11, the vignettes, or `DESCRIPTION` against the current file
state for its update pass (it only re-read series 02 in full), so its
carried-forward text for Major Issues 2.2, 2.3, 2.4, and 2.9 is stale
relative to the actual repository state confirmed directly in this
session. The 08-20 whitepaper does contain one genuinely new item,
Major Issue 2.14 (series 02's `cover-letter.md` contradicts the
current manuscript), which was out of this session's assigned scope
(2026-08-16 whitepaper, series 04/05/09 priority) and is listed under
Deferred below.

## 1. Fixed

- **Whitepaper Major Issue 2.9 / Checklist item 5 (series 04's two
  irreconcilable pipelines).** `[verified]` Already resolved by the
  prior pass: `analysis/report/04-treatment-main-effect/
  report_short.Rmd` carries a `[SUPERSEDED, NOT AUTHORITATIVE]` title
  and an in-document warning block naming `report.Rmd` as sole
  authoritative source; `README.md` updated to match. Confirmed by
  direct reading this session; no further change made.
- **Whitepaper Major Issue 2.10 / Checklist item 6 (series 05's stale
  N=40 baseline).** `[verified]` Already resolved by the prior pass:
  the boxed Results note was rewritten from a vague "in progress" to
  a specific diagnosis (the driver's `baseline_model_param()` already
  defaults to N=35; the twelve `.rds` sweep artifacts have simply not
  been regenerated) plus the exact rerun command. This session
  independently re-verified both claims: `grep` confirms
  `baseline_model_param <- function(N = 35, ...)` in
  `analysis/scripts/nof1-design-sensitivity/00-common.R`, and running
  `analysis/scripts/nof1-design-sensitivity/00-smoke-test.R` against
  the loaded package this session (via `pkgload::load_all()`) printed
  "Smoke test PASSED" with valid non-NA dimensions and a valid beta/p
  from `generate_simulated_results`. The full 12-sweep production
  rerun (~32 min per the manuscript's own estimate) was correctly
  deferred by the prior pass and remains deferred here (see Deferred).
- **Whitepaper Major Issue 2.11 / Checklist item 7 (series 06 Study B
  numbers not traceable).** `[verified]` Already resolved by the prior
  pass: a `sim_study_b <- readRDS(...)` call and a `study_b_cell()`
  helper were added to the setup chunk, and every Study B bias/
  coverage/power number in prose was replaced with inline `r`
  expressions reading from `sim_study_b`. This session independently
  loaded
  `analysis/data/derived_data/component-decomposition/study-b-balanced/
  study-b-balanced-summary.rds` directly and confirmed it is a valid
  16-row data.table with `bias`, `coverage`, `power`, and `conv_rate`
  columns matching what the inline expressions read (verified via
  `Rscript`, this session).
- **Whitepaper Major Issue 2.12 / Checklist item 8 (series 07 cover
  letter contradicts current finding; sloppiness diagnostic
  overclaimed).** `[verified]` Already resolved by the prior pass:
  `cover-letter.md` rewritten to report the architecture-conditional
  0.039 power spread and explicitly supersede the pre-fix null
  finding; `report.Rmd`'s Introduction narrowed the sloppiness/
  identifiability framing to "motivates... but is not itself computed
  in the present... submission... remains future work." Confirmed by
  direct reading this session.
- **Whitepaper Major Issue 2.6 / Checklist item 11 (series 04
  comparison-fairness disclosure).** `[verified]` Already resolved by
  the prior pass: a new `.orig` paragraph was added to the Limitations
  section disclosing the participant-matched-but-not-observation-
  matched confound explicitly, and stating no measurement-matched
  sensitivity arm has been run. Confirmed by direct reading this
  session.
- **Whitepaper Major Issues 2.2, 2.3, 2.4 (vignette `browser()` call,
  `simrsults`/`simresults` typo, undeclared vignette dependencies).**
  `[verified]` All three were already resolved by the prior pass, at
  the file level and confirmed structurally by this session: no
  `browser()` call exists anywhere in `vignettes/*.Rmd` (grep, this
  session); `rerun_simulations <- FALSE` at line 362 of vignette 01;
  vignette 02's `results_maxes` chunk (lines 140-146) consistently
  uses `simresults` throughout, no stray `simrsults`; `DESCRIPTION`
  declares `data.table`, `corpcor`, `ggplot2`, `MASS` under `Depends`
  (a stronger declaration than the `Suggests` the whitepaper
  recommended, which still satisfies `R CMD check`) and `merTools`
  under `Suggests`. All packages loaded by all three vignettes
  (including vignette 03's `library(lme4)`, `lmerTest`, `svMisc`,
  `tictoc`, `ggpubr`, `gridExtra`) are covered.
- **Whitepaper Major Issue 2.13 / Checklist item 9 (series 09
  placeholder and overclaimed mechanism), residual overclaim.**
  `[applied, unverified]` The prior pass already replaced the literal
  "(Forthcoming detail...)" placeholder with real Simulation Design
  prose. This session found and fixed a remaining overclaim the prior
  pass missed: the Introduction's "What this paper contributes" and
  the Abstract's Methods paragraph both still stated Sections 4-5
  deliver "the randomization-path-conditional power decomposition
  that quantifies the 'happy accident' selection effect" / "an
  explicit randomization-path decomposition," contradicting the
  Discussion's own honest statement that this is "consistent with,"
  not formally tested. Edited both passages in
  `analysis/report/09-informative-dropout-by-design/report.Rmd` to
  state the decomposition is scoped as future work and not delivered
  as a quantified result in the present submission, matching what
  Results/Discussion actually reports. Not re-rendered to PDF (see
  Deferred).
- **Whitepaper Minor Issue 3.9 (notation-linter violations in series
  01 and 02).** `[verified]` `tools/notation-lint.pl` flagged a
  `bare-Dbc` in `01/report.Rmd` line 352 and a `bare-Sx` in
  `02/report.Rmd` line 573; both were false negatives of the linter's
  protected-span stripping caused by an inline code/`\texttt{}` span
  being split across a markdown line wrap, which the linter's
  single-line regex does not stitch back together. Reflowed both
  spans onto single lines (`01/report.Rmd`,
  `02/report.Rmd`). Re-ran `perl tools/notation-lint.pl
  analysis/report/*/report.Rmd`: all eleven series now report
  "clean" (this session, verified).
- **Whitepaper Minor Issue 3.1 (series 01 README oversells
  incompleteness).** `[verified]` `README.md` and its "Known defects"
  section both claimed `rgt` blocks in `report.Rmd` still carry
  placeholders; `grep -c` confirms zero placeholder markers remain in
  the current `report.Rmd`. Updated the top-level bullet in
  `analysis/report/01-dgp-mean-moderation-vs-mvn/README.md` to state
  the zero-placeholder status; left the separate "Known defects"
  bullet about the same claim and the `[@pmsimstats-paper08]`
  raw-citekey concern untouched (see Deferred - could not verify the
  citation-rendering claim without a PDF render).
- **Whitepaper Major Issue 2.5 (empty notation-audit stub cited as
  load-bearing).** `[applied, unverified]` Writing the actual audit
  (cross-referencing `NOTATION.md`'s symbol table against all eleven
  manuscripts' usage) was out of budget. Per the remediation
  instructions' stated fallback, removed the false claim instead:
  edited `analysis/report/NOTATION.md` and
  `analysis/report/README.md` to state plainly that
  `whitepaper-notation-audit.md` is a front-matter-only stub as of
  2026-08-20 and that consistency currently rests on
  `tools/notation-lint.pl`, not a written audit trail.
- **Trivial placeholder test (`inst/tinytest/test-basic.R`).**
  `[verified]` Not directly named in the 2026-08-16 whitepaper, but
  matches its general "placeholder/trivial test suite" defect class
  and was cheap to fix. Replaced the single `expect_true(TRUE)` with
  three `expect_true(exists(..., mode = "function"), info = ...)`
  checks for the package's three core exported functions
  (`buildtrialdesign`, `generateData`, `generateSimulatedResults`).
  Ran the full suite: `tinytest::run_test_dir("inst/tinytest")`
  reports "All ok, 70 results" (this session, verified).
- **Whitepaper Minor Issue 3.5 (`[___]` placeholder citation in
  vignettes 01 and 02).** `[applied, unverified]` No target venue has
  been chosen (Section 5's partitioning decision is explicitly
  unresolved and requires author judgment), so per the remediation
  instructions' fallback, rephrased both vignettes' opening sentences
  to describe what the vignette demonstrates without promising a
  specific, not-yet-existing publication.

## 2. Deferred

- **Whitepaper Major Issue 2.1 / Checklist item 1 (series 11's
  boundary-validation gate).** Not fixed. Requires modifying
  `analysis/scripts/dgp-combined/01-run-combined-factorial.R` to use
  a per-randomization-path sample-size convention matching series
  01's production driver (currently hardcoded to `N == 35L` as a
  total, confirmed by `grep` this session), then re-running the full
  combined-architecture factorial (crossover/Hybrid/OL+BDC designs x
  3 carryover half-lives x 9 channel-weight cells x 1,000 replicates)
  and confirming the boundary rows/columns of Table 1 match series
  01's published power values within Monte Carlo tolerance. This is a
  code change plus a simulation run whose scope and runtime were not
  assessed as "a few minutes"; genuinely out of this session's budget
  and outside the assigned series 04/05/09 priority. The manuscript
  already honestly discloses this as "Status: pending re-run" in
  Section 4.2 and in its own abstract, so no false claim is currently
  being made. To complete: fix the sample-size convention in
  `01-run-combined-factorial.R`, run
  `Rscript analysis/scripts/dgp-combined/01-run-combined-factorial.R`
  followed by
  `Rscript analysis/scripts/dgp-combined/02-summarise-combined.R`,
  then re-stamp Table 1 and the abstract.
- **Whitepaper Checklist item 6 / Major Issue 2.10, full-scale rerun
  (series 05's twelve sensitivity sweeps at N=35).** Not run. The
  underlying code bug does not exist (the driver already defaults to
  N=35, confirmed by this session's smoke test); only the twelve
  `.rds` artifacts need regenerating. Manuscript's own estimate is
  ~32 minutes wall-clock, above the "a few minutes" threshold for
  this session. To complete: `Rscript
  analysis/scripts/nof1-design-sensitivity/run-all.R` from the
  repository root, then re-render `05/report.Rmd`.
- **Whitepaper Checklist item 9, full Type I error cells for series
  09 (OL+BDC, crossover, hybrid designs).** Not run. This is a
  simulation-scale expansion (three designs x additional cells,
  5,000 replicates/cell per the pre-registration), well beyond a
  "few minutes." The manuscript already discloses this gap honestly
  in the Discussion ("no formal Type I error cells... were run for
  the non-OL designs... A canonical production run should add
  explicit beta_bm = 0 cells"), so this is a completeness gap, not a
  false claim. No exact single command exists yet for this expansion;
  it requires extending the existing `09-*` driver scripts under
  `analysis/scripts/informative-dropout-by-design/` with explicit
  null cells before it can be run.
- **Whitepaper Checklist item 17 (series 09 through an internal
  referee cycle).** Not done; requires a human/author-level review
  process, not a code or prose fix this session can perform.
- **Whitepaper Checklist item 10 (writing `rgt` prose for the eight
  series still carrying placeholder blocks).** Not attempted. This is
  the largest single item in the whitepaper by volume (110+101+59+
  54+51+51+37+13 = 476 placeholder markers across 8 series) and is
  explicitly authorial/narrative work, not a defect fix; genuinely out
  of scope for a time-boxed remediation pass.
- **Whitepaper Checklist item 13 (partitioning decision - fold/promote
  series 11 relative to 01/03; settle 05's relationship to 04; settle
  07's relationship to 06).** Not done; the whitepaper's Section 5
  itself frames this as requiring authorial judgment about target
  venues and paper boundaries, not a mechanical fix.
- **Whitepaper Checklist item 14 (production run for series 03, or
  continued explicit pilot disclosure).** Not assessed this session;
  outside the series 04/05/09 priority and the whitepaper reports the
  pilot is already explicitly disclosed as such.
- **Whitepaper Checklist item 15 (second worked application for
  series 10) and item 16 (README/cover-letter for series 10 and
  11).** Not attempted; outside the series 04/05/09 priority, and
  item 15 in particular is a substantive new-analysis addition, not a
  defect fix.
- **Whitepaper Checklist item 18-25 (desirable polish: series 06
  word-count reconciliation, series 07 README/scope reconciliation,
  `@unpublished` cross-citation audit, etc.).** Not attempted beyond
  items 3.1, 3.5, and 3.9 above, per the instructions' guidance not to
  let (c)-tier polish crowd out (a)/(b) work.
- **2026-08-20 whitepaper's new Major Issue 2.14 (series 02
  `cover-letter.md` contradicts the current decay-shape-robustness
  finding).** Out of this session's assigned scope (2026-08-16
  whitepaper, series 04/05/09 priority); flagged here for a future
  pass. `analysis/report/02-carryover-sensitivity/cover-letter.md`
  needs the same treatment as series 07's cover letter received in
  the prior remediation pass.
- **No PDF re-render performed.** Per the task's hard constraints,
  rendering was optional and `rmarkdown::render()` was not to be
  called directly. `bash tools/render.sh` exists in this repository
  but was not invoked this session to keep the pass source-level and
  time-boxed; all edited `.Rmd`/`.md` files should be re-rendered
  before any of them is treated as submission-ready.

## 3. New issues found while fixing

- The prior 2026-08-16 remediation pass's own cover-letter edit for
  series 07 cites `docs/pub_review_remediation_2026-08-16.md` as the
  record of that pass, but that file does not exist anywhere in the
  repository (`ls docs/pub_review_remediation*.md` before this
  session's write found nothing). The prior pass's work is real and
  verified in this session, but its own documentation trail was never
  written or was lost. This session's log (this file) is the first
  actual remediation log in this repository.
- `DESCRIPTION` declares `data.table`, `corpcor`, `ggplot2`, and
  `MASS` under `Depends` rather than `Suggests`. This satisfies `R CMD
  check` (the whitepaper's actual concern) but is a stronger coupling
  than the whitepaper's stated remediation suggested; `Depends`
  forces these packages into every user's search path on
  `library(pmsimstats)`, which is worth revisiting as a design choice
  independent of the CRAN-compliance question, though it is not a
  correctness defect.
- Series 01's README "Known defects" section separately claims a
  citation `[@pmsimstats-paper08]` renders as a raw, unprocessed
  citekey because it sits inside a raw-LaTeX `bullets` environment
  that bypasses pandoc citation processing. Direct reading of
  `report.Rmd` this session found the citation in ordinary Markdown
  prose in Section 1 (Introduction), not inside a `bullets`/`rgt`/
  `orig` div, so the stated mechanism for the defect does not match
  its stated location; this claim could not be confirmed or refuted
  without a PDF render, which was out of scope this session. Worth
  a render-and-check pass to confirm whether this citation actually
  renders correctly or the README's underlying concern is valid but
  mislocated.
