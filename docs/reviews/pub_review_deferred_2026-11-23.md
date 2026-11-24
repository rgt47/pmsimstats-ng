# pub_review deferred-item pass
*2026-11-23 13:06 PST*

Scope: every Deferred item in

- `docs/reviews/pub_review_remediation_2026-08-20.md` (workspace
  level, answering `docs/reviews/pub_review_whitepaper_2026-08-16.md`;
  'A' below), and
- `analysis/report/01-dgp-mean-moderation-vs-mvn/docs/pub_review_remediation_2026-11-23.md`
  (series 01, answering that series'
  `pub_review_whitepaper_2026-08-26.md`; 'B' below).

Every item was checked against the files as they stand today and
against `git log --since=2026-08-20` (272 commits,
`git rev-list --count`). No commits,
staging or pushes were made. The working tree already held
uncommitted author changes before this pass (01 and 02 `report.pdf`
and `report.tex`, 01 `README.md` and `bullets.Rmd` from this morning's
series-01 pass, `share/MANIFEST.md`, many untracked data folders and
`share/` PDFs); none of those were touched except as listed in
section 2.

Counting convention: bundled bullets in the source logs were split
into one row per whitepaper item. A's 'items 15 and 16' bullet became
three rows (second application; READMEs; cover letters). A's 'items
18-25' bullet became five rows (18, 20, 22, 23, 25); items 19, 21 and
24 were already fixed by the 08-20 pass and are not counted. B's two
'new issues' that call for action (Dropbox paths; hand-typed numbers)
are counted as items B5 and B6.

Evidence labels: **verified** (ran code and observed the result),
**inspected** (read the file or diff), **inferred**.

## 1. Inventory

| ID | Source | Item | Class |
|---|---|---|---|
| A1 | 08-20 | Series 11 boundary-validation gate and factorial rerun | RESOLVED-SINCE |
| A2 | 08-20 | Series 05 twelve sweeps rerun at N = 35 | RESOLVED-SINCE |
| A3 | 08-20 | Series 09 Type I cells for OL+BDC, CO, Hybrid | RESOLVED-SINCE |
| A4 | 08-20 | Series 09 internal referee cycle (item 17) | RESOLVED-SINCE |
| A5 | 08-20 | `rgt` prose for eight series (item 10) | RESOLVED-SINCE |
| A6 | 08-20 | Partitioning decision (item 13) | AUTHOR-DECISION |
| A7 | 08-20 | Series 03 production run or pilot disclosure (item 14) | AUTHOR-DECISION |
| A8a | 08-20 | Series 10 second application (item 15) | AUTHOR-DECISION |
| A8b | 08-20 | README for series 10 and 11 (item 16) | DO-NOW |
| A8c | 08-20 | Cover letter for series 10 and 11 (item 16) | AUTHOR-DECISION |
| A9-18 | 08-20 | Series 06 word-count reconciliation | DO-NOW |
| A9-20 | 08-20 | Vignette 02 proofreading | DO-NOW |
| A9-22 | 08-20 | Series 09 robustness to the other encoding and wider effect size | AUTHOR-DECISION |
| A9-23 | 08-20 | Series 07 README promises unsimulated regimes | RESOLVED-SINCE |
| A9-25 | 08-20 | `@unpublished{pmsimstats-paperNN}` cross-citation audit | DO-NOW |
| A10 | 08-20 | Series 02 `cover-letter.md` contradicts manuscript (2.14) | DO-NOW (label); rewrite AUTHOR-DECISION |
| A11 | 08-20 | No PDF re-render after the 08-20 edits | RESOLVED-SINCE |
| B1 | 11-23 | Render `bullets.Rmd` | DO-NOW |
| B2 | 11-23 | Rewrite series 01 README scope bullets | AUTHOR-DECISION |
| B3 | 11-23 | Retire or rewrite `bullets.Rmd` | AUTHOR-DECISION |
| B4 | 11-23 | Referee reports of 2026-11-01 and 2026-11-05 | DO-NOW (status table) |
| B5 | 11-23 | Dropbox paths in two Hendrickson comparison scripts | DO-NOW |
| B6 | 11-23 | Hand-typed N = 70 / N = 140 power values in series 01 | AUTHOR-DECISION |

Totals: 23 items; 7 resolved-since; 0 no-longer-applies; 8 done now;
0 long-run; 8 author-decision.

### Evidence for the RESOLVED-SINCE rows

- **A1 (series 11).** Commit `6ea37805e` (2026-10-29) rewrote the
  driver with common random numbers and a draw-for-draw validation
  gate, and reran at N = 70. Verified: `01-combined-gate.rds` holds
  1,800 paired boundary comparisons and the maximum absolute
  difference in data, estimate and p-value is 0 in every row;
  `01-combined-meta.rds` has `gate_pass = TRUE`, 1,000 replicates,
  three designs, half-lives 0, 0.5, 1. `report.tex` contains no
  'pending' text (the Rmd's pending strings are fallback branches).
  Note: the original target of the gate (series 01's published power
  values) no longer exists; the current gate compares each boundary
  cell with the single-encoding DGP of the same driver, and the
  manuscript says explicitly that it does not compare with paper 01's
  numbers. The run was made with uncommitted driver changes, which
  the Reproducibility section discloses.
- **A2 (series 05).** Commit `e238025cf` ('report the full rerun,
  zeus job 51'). Verified: every `sensitivity-01` to `-12` `.rds` is
  dated 2026-10-29, and `cell-summaries.csv` has 414 cells at N = 35
  and 3 each at the swept N values 16, 24, 60, 80, 120. The
  manuscript (line 94 and the Abstract) states N = 35 as the
  baseline.
- **A3 (series 09).** Commit `8918e5c6e` rebuilt the paper on a
  paired design. Verified: `paired-cells.csv` in
  `analysis/data/informative-dropout-by-design-published/` has 28
  null-regime cells for each of CO, Hybrid and OLBDC at 1,000
  replicates, and the manuscript reports size per design ('Size MB',
  'Size CR2'). The pre-registration's 5,000 replicates per null cell
  were not used; the MCSE is reported.
- **A4.** `09/referee-report-2026-10-28.md` exists and `8918e5c6e`
  addresses it (inspected).
- **A5.** `git log -S'To be completed by rgt'` places the removal in
  `3ce116216` (2026-09-07). Verified by grep: zero `Lorem ipsum`,
  `To be completed by rgt`, `To be finalized` or `.rgt` markers in
  any of the fifteen `report.Rmd` files or their child Rmds (a
  recursive grep over all `*.Rmd` under `analysis/report/`). Outside
  the canonical files they remain only in `archive/` copies (01, 02)
  and in 06 `report-slim.Rmd` (section 5, item 3).
- **A9-23.** `07/README.md` was revised in `0c0f1e260` (2026-10-28)
  and now lists the families actually simulated (logistic,
  hyperbolic tangent, piecewise-linear breakpoint), which
  `report.Rmd` computes (inspected).
- **A11.** The rendered files the 08-20 pass edited have since been
  re-rendered: `share/` holds 09 renders of 2026-11-05 and 01 and 02
  renders up to 2026-11-21/22; `NOTATION.pdf` is dated 2026-10-09,
  after the 08-20 edit (inspected). `NOTATION.md` was edited again by
  the author on 2026-11-18 (`531deb672`) and its PDF predates that
  edit; that is new author work, not the 08-20 item. The vignettes
  were not rebuilt (not checked).

## 2. Done now

- **A8b. READMEs for series 10 and 11.** New files
  `analysis/report/10-interaction-test-calibration/README.md` and
  `analysis/report/11-combined-dgp-architecture/README.md`. Strictly
  factual: current title, last render (`pdfinfo`), file roles,
  inputs, driver and render commands taken from the manuscripts'
  own build notes; no findings summarized. `[verified]` for each
  stated fact (the role of `claudecode.tex` is left undocumented):
  page counts and dates by `pdfinfo`; series 11 gate
  result and run time (10,109 s, about 169 minutes, eight workers)
  read from the stored `.rds` files; child-document structure by
  grep.
- **A9-18. Series 06 word count.** Re-measured by one method
  (`wc -w`): `report.Rmd` 21,602 words; `report-slim.Rmd` 12,653
  words, of which 67 lines are `Lorem ipsum` placeholders (the slim
  file was last committed 2026-10-09). The slim variant is therefore
  not a drop-in canonical version. `[verified]`. Whether to promote,
  complete or drop it is the author's call.
- **A9-20. Vignette 02 proofreading.** The 'comapered' typo cited by
  the whitepaper is already gone. Fixed the remaining misspelled axis
  label 'Signifiance' (two occurrences) in
  `vignettes/02-visualize-power-results.Rmd`. Greps for common
  misspellings (all three vignettes) and for British spellings (all
  three vignettes) found nothing else; this is a pattern check, not a
  full proofread. `[applied, unverified]`: the vignette was not rebuilt
  (label text only).
- **A9-25. Companion-paper citation audit.** Script compared every
  `@unpublished{pmsimstats-paperNN}` entry in the fifteen series
  `references.bib` files with the current `title:` of the cited
  series' `report.Rmd`. 24 of 38 entries carried superseded titles
  (for example, paper 02 as 'Carryover representation and inference
  ... under covariance moderation: a factorial simulation comparison
  of nine analysis specifications'; paper 01 as 'Simulating
  predictive-biomarker interactions ...' or 'Two architectures ...').
  Rewrote the `title` field of those 24 entries to the current titles
  (papers 01, 02, 05, 09, 10, 11), keeping keys, authors, notes and
  years, with `{N}`, `{I}`, `{REML}`, `{Hybrid}` brace-protected.
  Files: the `references.bib` of series 01, 03, 08, 10, 11, 12, 13,
  14, 15. `[verified]`: re-running the audit reports 0 stale titles
  of 38; `pandoc -f bibtex -t csljson` parses all nine files;
  `git diff -U0` shows only `title` lines and their continuation
  lines changed, with no end-of-file newline changes. Originals are backed up at
  `/private/tmp/claude-501/-Users-zenn-prj-res-01-fisher-exact-rx2-fisherexacttestrx2-analysis-report/be4b228f-9df5-41a4-9039-7499d4f23926/scratchpad/deferred/bib_backup/bibs.tar`
  (scratch space; may be cleared). Not re-rendered, so the nine PDFs
  still print the old titles. The whitepaper's broader point (whether
  each cited manuscript will be available when the citing paper is
  submitted) remains an author decision.
- **A10. Series 02 cover letter.** Verified that the letter's title
  and 'Tier 1/Tier 2' design occur nowhere in the current
  `report.Rmd`, `report.tex` or `supplement.tex`; the current title is
  'Testing a predictive biomarker under drug carryover in aggregated
  N-of-1 trials: a staged simulation study'. (0.830 does occur once in
  `report.tex`, in an unrelated table row for A1, so the notice makes
  no claim about the letter's figures.) Prepended a dated 'Superseded,
  do not send' notice to `analysis/report/02-carryover-sensitivity/cover-letter.md`
  (the same device used for 04 `report_short.Rmd` and 01
  `bullets.Rmd`); the letter body is unchanged. `[verified]` by grep.
  Rewriting it is left to the author (section 4).
- **B1. Render `bullets.Rmd`.** Ran
  `bash tools/render.sh analysis/report/01-dgp-mean-moderation-vs-mvn/bullets.Rmd`.
  It succeeded in seconds. `[verified]`: `pdftotext` shows the
  'Superseded (noted 2026-11-23)' notice and the '[Unverified legacy
  result.]' label in the new `bullets.pdf`. Side effects: `bullets.pdf`
  and `bullets.tex` rewritten; new
  `analysis/report/share/bullets-2026-11-23-1305-whitepapers-pre-series-226-ga832d566c-wip.pdf`;
  one row appended to `analysis/report/share/MANIFEST.md` (confirmed by
  `diff` against a pre-render copy; the author's existing uncommitted
  changes to that file are otherwise untouched).
- **B4. Referee reports 2026-11-01 and 2026-11-05 (series 01).**
  Bounded status check, not a remediation; all verdicts are
  **inspected** (grep and reading of the current `report.Rmd`, 3,084
  lines). The report's line numbers are stale, so text was located by
  phrase or label. G10 is now labeled A4.

  | 11-05 item | Status | Evidence in current `report.Rmd` |
  |---|---|---|
  | M1 rho dependence | Partly addressed | 'no PTSD symptom series' removed; sensitivity at rho 0.5 and 0.9 reported; 'largest feasibility bound' still given as the reason for 0.7 (about line 1238), which remedy (c) asked to drop or restate |
  | M2 A4 under prognostic biomarker | Partly addressed | Limitation states A4 loses size with a prognostic biomarker (about line 2221); prognostic cells moved to the companion paper; whether A4 was run in the prognostic null cells not determined |
  | M3 validity table under DGP as defined | Addressed | Third run with fixed allocation is the reported one; the three runs are distinguished (lines 692-700, 1385, 1570) |
  | M4 overstated case for the correlation encoding | Partly addressed | 'property of the design' removed; implementing the random-coefficient sensitivity model not verified |
  | M5 literature | Largely addressed | Kraemer, Ding, Higham 2016, Ghosh, Joe, Siepe, Pawel, Senn 1992/1998/2004/2011, Jones-Kenward, Haller, Kent, Tueller, Naik, Holford all cited; Web of Science and Scopus not searched, stated as a limitation (line 214) |
  | M6 internal documents, Declarations | Open | `docs/69` and `docs/47` no longer cited; CSP 563 data-use confirmation is still a visible `[Placeholder ...]` (about line 2485); competing interests says 'no competing interests' without the overlap statement the referee asked for |
  | Minor 9 `tab:nonpd` typed | Addressed | Table now read from `nonpd-by-rho.csv` |
  | Minor 14 Siepe, Pawel | Addressed | Both cited |
  | Minor 15 companion title | Addressed now | `pmsimstats-paper02` title updated in A9-25 |
  | Other minors (1-8, 10-13, 16) | Not checked | Out of time box |

  For 11-01, the existing `revisions/referee-2026-11-01-response-audit.md`
  listed M2, M7, Minor 2, Minor 9 and the Section 6 references as
  partly or not addressed. M2 (cause of changes separated) and M7
  (database coverage stated) are now addressed as in M3 and M5 above;
  the Section 6 references are cited except `tang2020`; Minor 9
  (Declarations) remains open as in M6.
- **B5. Dropbox paths.** In
  `analysis/scripts/quick-sim/hendrickson-original-comparison/01-hendrickson-orig-driver.R`
  and `02-pd-sweep.R`, replaced
  `load_all('/Users/zenn/Library/CloudStorage/Dropbox/...')` with
  `load_all('.')`, the convention of every sibling quick-sim script
  (run from the repository root). `[applied, unverified]`: both files
  parse; the scripts were not run (the current manuscript does not
  cite them).

## 3. Long runs

None of the remaining items is a pure compute job that can be
started as is:

- A7 (series 03): the corrected study has already run (see section
  4); what remains is integration.
- A9-22 (series 09): `20-paired-dropout.R` has no option for the
  encoding, N or effect size (its environment variables are
  `P09_REPS`, `P09_NULL_REPS`, `P09_CORES`, `P09_OUT`, `P09_TAG`), so a
  code change must come first. For scale, the stored paired run took
  3,552 s (about 59 minutes) for one parameter point (`paired-meta.rds`).

For context, `docs/75-manuscript-readiness-assessment-2026-11-09.md`
lists the queued reruns the author already tracks (paper 10 at about
10 CPU-hours; the merged 13/14 grid; 05 after an engine fix).

## 4. Author decisions

- **A6. Partitioning.** Options: keep fifteen manuscripts; adopt the
  five merges proposed in `docs/66` (13 + 14; 08 into 10; 07 + 08; 12
  into 01; 11 into 01 or 03), which `docs/75` refines ('11 into 01, 12
  retired, 07's placement'). Evidence: both documents are the
  author's own and mark the merges as unagreed. Recommendation: decide
  scope for 04, 07, 08, 11 and 12 before further rewriting, as
  `docs/66` itself orders.
- **A7. Series 03 corrected study.** Both blocks ran under the
  published parameter set (verified:
  `analysis/data/latent-class-mixture-application/published/referee-rerun-null/`
  with 2,000 replicates per null cell, 21 summary rows;
  `.../referee-rerun-alt/` with 1,000 per alternative cell, 105 rows;
  `smoke = FALSE`), but these folders are untracked and `report.Rmd`
  reads `analysis/data/latent-class-mixture-application/referee-rerun/summary.rds`,
  which does not exist, so the rendered Results print 'has not
  completed'. `docs/75` already records this and that the results do
  not support the draft's anticipated conclusion. Options: point the
  Rmd at the published folders and rewrite the Results, Discussion and
  Abstract; or keep the pilot-only paper. Recommendation: integrate;
  the run is complete and the current text is out of date. Not done
  here because the conclusions change.
- **A8a. Series 10 second application.** A new analysis; scope
  decision. Recommendation: defer until the paper 10 published-set
  rerun is integrated (section 5, item 2).
- **A8c. Cover letters for 10 and 11; rewrite of 02's.** Need a
  target venue and the author's framing. Recommendation: write them
  from the current Abstracts at submission time; 11's depends on the
  A6 decision (fold or keep).
- **A9-22. Series 09 robustness.** The manuscript discloses 'one
  parameter point' and 'one encoding' as limitations. Options: extend
  the driver and run the correlation encoding and at least one other
  effect size (about an hour per point at the stored settings), or
  keep the disclosed scope. Recommendation: keep the scope for a
  first submission; the limitation is stated plainly.
- **B2. Series 01 README scope bullets.** Needs the author's
  restated outline. Recommendation: rewrite after the current
  revision of paper 01 settles.
- **B3. Series 01 `bullets.Rmd`.** Options: delete, move to
  `archive/`, or rewrite for the current paper. Precedent:
  `5dfc1eafa` archived paper 02's bullets outline. Recommendation:
  archive it the same way; the same applies to the bullets files of
  10 and 11 (section 5, item 4).
- **B6. Hand-typed N = 70 / N = 140 values in series 01.** They match
  the stored CSVs today (verified by the 11-23 series-01 pass).
  Converting them to inline `r` is mechanical, but `report.Rmd` is
  under active uncommitted revision (its PDF and TeX are modified in
  the working tree), so it was not edited here. Recommendation:
  convert them in the author's next commit to that file.

## 5. New issues found

1. **Series 03 results not wired** (A7 above). Already recorded in
   `docs/75`; listed here because it decides A7.
2. **Series 10 has an unused published-set run.** An untracked
   `analysis/scripts/carryover-sensitivity/output/paper10-published/`
   exists, but `report.Rmd` reads only `paper10-refined/` and the
   pre-correction files. Other series have moved to the published
   parameter set. Whether 10 should is an author decision (inspected;
   contents of that folder not examined).
3. **Series 06 `report-slim.Rmd`** carries 67 `Lorem ipsum`
   placeholder lines (verified by grep). Harmless while it is not the
   canonical file; it must not be promoted as is.
4. **Stale side files beyond 02.** The 04 cover letter names a title
   ('Comparative Statistical Power of Aggregated N-of-1 and
   Parallel-Group Trial Designs ...') that differs from the current
   04 title; the 10 and 11 `bullets.Rmd` carry superseded titles, and
   11's names 'channel super-additivity', which the paper withdrew in
   `6ea37805e` (inspected). `docs/66` notes the general pattern.
5. **Nine PDFs need re-rendering** for the bib title updates (A9-25):
   series 01, 03, 08, 10, 11, 12, 13, 14, 15. Not rendered here:
   01 is under active revision and several others are long builds.
6. **Host cannot load the package.** `pkgload::load_all()` stops with
   'The packages "furrr", "svMisc", and "tictoc" are required', so
   the tinytest suite could not be run on the host this session
   (verified). The 08-20 pass ran it successfully; the container
   (`make r`) presumably still has them (not checked).

## Tests and renders run

- `Rscript` evidence checks (series 05, 09, 11 stored outputs; series
  03 published outputs; series 09 run metadata).
- Companion-title audit script before and after the bib edits; pandoc
  parse of the nine edited bib files.
- `parse()` of the two edited R scripts.
- `bash tools/render.sh .../01-dgp-mean-moderation-vs-mvn/bullets.Rmd`
  (succeeded).
- tinytest suite: attempted, could not run (section 5, item 6). No
  package code was changed.
