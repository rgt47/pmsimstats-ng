# pub_review remediation log: series 01
*2026-11-23 12:24 PST*

Remediation pass against
`analysis/report/01-dgp-mean-moderation-vs-mvn/docs/pub_review_whitepaper_2026-08-26.md`.
Scope: series 01 only, per the task assignment. No commits, staging
or pushes were made; all changes are in the working tree.

## 0. Principal finding: the whitepaper is stale against the current manuscript

The 2026-08-26 whitepaper reviewed a manuscript titled 'Mean
moderation and covariance moderation as data-generating architectures
for biomarker-treatment interactions in aggregated N-of-1 trials, and
their divergent power under carryover' (17 pages). The current
`report.Rmd` is a different paper:

- Title: 'Toward an optimal data-generating process for simulating
  predictive-biomarker interactions in aggregated N-of-1 trials'.
- Render of 2026-11-21: 72 pages (`pdfinfo report.pdf`).
- None of the whitepaper's flagged numbers or sections are present.
  `grep` for `0.978`, `0.945`, `0.885`, `30.6`, `7.64`, `Analytical
  comparison`, `Robustness`, `9,600`, `0.042`, `0.074`, `mean
  moderation`, `None of these strategies`, `wider literature` returns
  nothing.

How the supersession was confirmed (inspected, from git history):

- `git log -S0.885` on `report.Rmd` shows the old Section 3.3 (the
  N = 140 robustness check, table 0.978/0.945/0.885) was removed in
  commit `4b02b5ba0` (auto-backup, 2026-08-26 10:42), about two
  hours after the review. That diff deletes
  `\subsection{3.3 Robustness check at \emph{N} = 140}` and its
  table.
- Commit `610688bb4` (2026-10-27, 'Papers 01 and 02: remove
  mean-moderation mentions and tabled values') removed the remaining
  mean-moderation material. Later commits through `a8c7d0a40`
  (2026-11-21) reframed the paper around an optimal DGP.

The live review record for the current manuscript is therefore the
later referee reports in the series directory
(`referee-report-2026-11-01.md`, `referee-report-2026-11-05.md`), not
this whitepaper. This pass did not address those reports; they are
outside the assigned whitepaper.

## 0.1 What the 2026-08-20 workspace-level pass already fixed for series 01

From `docs/reviews/pub_review_remediation_2026-08-20.md`:

- Notation-linter `bare-Dbc` false negative in `01/report.Rmd` (then
  line 352), fixed by reflowing a split inline span. Re-checked this
  session: `perl tools/notation-lint.pl
  analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd` reports
  'clean' [verified].
- README top-level bullet updated to the zero-placeholder status.
  The 08-20 pass left the 'Known defects' raw-citekey bullet
  untouched; by the start of this pass it read 'None recorded'
  (presumably an author edit after 08-20; the commit was not
  traced). This pass corrected a residual
  inaccuracy in that bullet (see 3.1 below).

## 1. Fixed

- **Major 2.1 / checklist (a)1 (N = 140 robustness check with no
  traceable source).** In `report.Rmd`: superseded, not edited; the
  section was removed in `4b02b5ba0` (see section 0). The same
  unsourced content survived in the companion file `bullets.Rmd`
  (Section 3.3: 12 cells at N = 140, 800 replicates, 9,600 fits,
  9.3-point loss, z = 7.4). Applied the whitepaper's option (c)
  there: added a visible '[Unverified legacy result.]' note at the
  head of Section 3.3 stating that no script or data in the
  repository generates it and that it must not be cited, plus a
  boxed 'Superseded' notice under the opening paragraph stating that
  the summarized manuscript no longer exists in this form. File:
  `analysis/report/01-dgp-mean-moderation-vs-mvn/bullets.Rmd`.
  `[applied, unverified]` (source edit only; not rendered).
- **Current N = 140 claims in `report.Rmd` (successor to 2.1).** The
  current manuscript reports a different N = 140 result: the code of
  Hendrickson et al. in their Figure 4B layout, extended to N = 140
  (lines about 1935-1964). Checked traceability: the driver
  `analysis/scripts/quick-sim/hendrickson-problems/15-fig4b-graded-58b32a9.R`
  is git-tracked and parses `--N` (comma-separated list, default
  35,70); its output
  `analysis/data/quick-sim/hendrickson-problems/fig4b-graded/fig4b-graded-cbm0.25-reps250-N140.csv`
  is git-tracked. Every typed number in the N = 70 bullets (Hybrid
  0.66 to 0.08/0.08 versus 0.62/0.68; OL+BDC 0.57 to 0.12/0.13
  versus 0.57/0.60; crossover 0.34/0.33 versus 0.71/0.78) and the
  N = 140 sentence (Hybrid decaying 0.94/0.95 versus constant
  0.12/0.10) reproduces from the stored CSVs. The claim 'the same
  pattern holds at N = 35 and N = 140 and under every censoring
  pattern ... At N = 140 the gap widens' also holds: across all 90
  non-open-label cells with positive half-life (3 sample sizes x 3
  designs x 2 half-lives x 5 censoring patterns), decaying power
  exceeds constant power in every cell (minimum gap 0.052), mean
  gap rises 0.19, 0.36, 0.53 at N = 35, 70, 140, and the two
  versions are identical at zero half-life. No repairs in any cell.
  No edit needed. `[verified]` (recomputed with R from the CSVs this
  session).
- **Checklist (b)3 / Minor 3.1 and 3.2 (README stale claims and page
  count).** File:
  `analysis/report/01-dgp-mean-moderation-vs-mvn/README.md`.
  `[verified]` for the facts stated (each checked by grep or
  `pdfinfo` this session):
    - page count changed from '21 pages' to '72 pages' with the render
      date;
    - the false statement 'the document loads no data' replaced with
      an accurate one (it calls `read.csv()` on stored results under
      `analysis/data/`, about 19 calls, and includes stored figures);
    - required CSL corrected from `statistics-in-medicine.csl` to
      `vancouver.csl` (the YAML header's value);
    - the Known-defects bullet corrected: `pmsimstats-paper08` is not
      cited in the current `report.Rmd` at all, rather than 'now in
      plain prose';
    - added a note that the scope bullets describe the 2026-07-31
      consolidation and are stale, with the confirmed facts (current
      title, mean moderation removed 2026-10-27, Hybrid at N = 70 as
      primary case), leaving the full rewrite to the pending outline;
    - the `bullets.Rmd` description corrected (it summarizes the
      pre-2026-08-26 manuscript and is no longer cited by the
      master), and a Known-defects bullet added recording its
      superseded status.
  The 'rgt placeholder' self-contradiction flagged in 3.1 was already
  gone before this pass.
- **Checklist (c)7 (`run-01-dgp.sh` hard-codes a Dropbox path).**
  File: `analysis/scripts/quick-sim/run-01-dgp.sh`. Replaced the
  absolute `cd /Users/zenn/Dropbox/...` with a script-relative
  `cd "$(dirname "$0")/../../.."`. Confirmed `01-dgp-prototype.R`
  still exists and is tracked, and that no other file references the
  wrapper. `bash -n` passes. `[applied, unverified]`: `bash -n` is a
  syntax check only; the simulation was not launched (it is the
  driver for the superseded 36-cell mean-versus-covariance factorial
  and is not cited by the current manuscript).

## 2. Not applicable (superseded) or not acted on

- **Major 2.2 (Section 3.1/3.4 numbers reproduce from
  `01-dgp-summary.txt`).** Positive finding; no action. The numbers
  it verified are no longer in the manuscript.
- **Major 2.3 / checklist (a)2 (add a direction-versus-magnitude
  sentence to Section 2.2.3).** Not applicable. Section 2.2.3
  'Analytical comparison' no longer exists. The current manuscript
  derives closed-form power for the paired-difference statistic E9
  (Appendix A.1-A.2) and checks it against Monte Carlo, so it now
  states magnitudes analytically for that statistic; the concern
  does not carry over. No edit.
- **Minor 3.3 (Type I error median 0.035 versus 0.0345).** Not
  applicable; the sentence is gone.
- **Minor 3.4 / checklist (c)5 (Section 4.3 mitigation strategies
  untested).** Not applicable; that section is gone.
- **Checklist (b)4 (re-read the extended Section 4.1 against Section
  4.4 for redundancy).** Not applicable; neither section exists in
  the current Discussion (headings: 'The correlation encoding of the
  interaction', 'What remains uncalibrated', 'Recommendations',
  'Limitations'). Whether the boundary-regime material added during
  the review survives anywhere was not traced.
- **Checklist (c)6 (data-availability pointer to
  `01-dgp-summary.txt` and `01-dgp-prototype.R`).** Not applied as
  written: pointing the current paper to data it no longer reports
  would mislead. The intent is already met. The current 'Availability
  of data and materials' section lists each driver script, its flags
  and its output directory.

## 3. Deferred

- **Render.** No `.Rmd` that feeds the master was changed. Only
  `bullets.Rmd` (a superseded companion) received source edits. Not
  rendered, to avoid adding share/ and MANIFEST rows for a historical
  file. To refresh its PDF:
  `bash tools/render.sh analysis/report/01-dgp-mean-moderation-vs-mvn/bullets.Rmd`
  from the repository root.
- **Rewriting the README scope bullets in full.** Left as a stale
  flagged note pending the author's restated outline (README
  'Known defects', 2026-10-07 entry). This needs author judgment
  about the paper's current framing.
- **Retiring or rewriting `bullets.Rmd`.** Labeled rather than
  deleted or rewritten. Deciding whether to delete it, move it to
  `archive/`, or rewrite it for the current paper is the author's
  call.
- **Referee reports of 2026-11-01 and 2026-11-05.** Out of scope for
  this pass (the assignment names the 2026-08-26 whitepaper only).
  They are the live review record and should be checked for open
  items separately.

## 4. New issues found while fixing

- `bullets.Rmd` still carried the unsourced N = 140 robustness
  numbers that the whitepaper flagged in `report.Rmd`, under the old
  title and framing. It also still uses `statistics-in-medicine.csl`
  and cites `pmsimstats-paper08`. Labeled as above.
- Two scripts in
  `analysis/scripts/quick-sim/hendrickson-original-comparison/`
  (`01-hendrickson-orig-driver.R`, `02-pd-sweep.R`) also contain the
  old Dropbox path. The current manuscript does not cite them, so
  they were not edited.
- The N = 70 and N = 140 power values in the 'Differences from the
  DGP of Hendrickson et al.' bullets (about lines 1947-1963 of
  `report.Rmd`) are typed by hand rather than computed inline from
  the CSVs. They match the data today (verified above), but unlike
  the neighboring tables they could drift if the run is redone.
  Consider replacing them with inline `r` expressions.
