# Decision history

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Chronological record of every methodological decision:
date, recommendation, replacement, evidence.
**Prerequisites:** None

---

## Contents

1. [Notation](#notation)
2. [Summary](#summary)
3. [Construct box](#construct-box)
4. [Setting and argument](#setting-and-argument)
5. [Derivations](#derivations)
6. [Numerical evidence](#numerical-evidence)
7. [Relation to the manuscripts](#relation-to-the-manuscripts)
8. [Limitations and open items](#limitations-and-open-items)
9. [Reproduce](#reproduce)
10. [References](#references)
11. [Change log](#change-log)

---

## Notation

*Local symbols used in this chapter only. For project-wide notation
see `analysis/report/NOTATION.md`.*

<!-- notation entries -->

---

## Summary

This appendix assembles every dated methodological decision from the
numbered whitepapers in chronological order. The format is: Date |
Decision | What replaced it | Evidence | Source.

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: the dated notes of 35, 41 to 45, 49, 50, 56, 57; 52
  Section 2; 30 execution record; 45 Section 9; 57 superseded plans.
- New to write: all.
-->

<!-- TODO: [Build the full decision history table from the dated
notes in each source document. Each row should record:
  Date | Decision | What it replaced | Evidence label | Source] -->

The table below collects dated decisions extracted from the source
documents on 2026-10-15 (inspected). The decision content is
inferred from the dated notes in each source; the full evidence
behind each item is in the cited document.

| Date | Decision | Replaced by | Evidence | Source |
|---|---|---|---|---|
| 2026-04-15 | Parity workflow established for tidyverse development; `--quick` (16 cells, ~2 s) and full (144 cells, ~10 s) modes defined | Still current; updated to include divergence CSV mechanism 2026-10-07 | inspected | docs/22 (now `analysis/scripts/parity/README.md`) |
| 2026-04-16 | Code audit: T1 (unqualified dplyr calls), C1 (silent PD repair), C2 (singularity flag), C5 (negative half-life), C6 (stale imports), C7 (missing test coverage) identified as P0/P1 blockers | T1 resolved; C6, C7 still open (2026-10-07) | inspected | docs/legacy/23 Section 14 |
| 2026-07-30 | Terminology decisions D1 to D7 locked and Tier 1 applied across 21 manuscript files; US English direction set | Current | verified | docs/30 Section 8 |
| 2026-07-31 | Analysis specifications relabeled from M1/M2/M3 to E1/E2/E3 in stored data, drivers and manuscripts | Current | inspected | docs/30; NOTATION Part 3 |
| 2026-09-30 | docs/35 proposed the separable AR(1) covariance structure as the refined process covariance-moderation DGP | Steps 1-3 (package integration) still pending 2026-10-09 | verified | docs/35 Section 9; docs/44 Section 4 |
| 2026-10-07 | Whitepaper review: expectancy symbol $e_t$ renamed $\eta_t$ across all 34 documents; ceiling values harmonized (0.256 CS, 0.440 at $\rho=0.8$, 0.474 Hybrid at $\rho=0.7$) | Current | verified | docs/legacy/52 Section 2 |
| 2026-10-07 | Stress-test outcome reported: F1 selected over 14 analyses; unstructured KR recommendation of docs/45 superseded | F1 is the recommended analysis | verified | docs/legacy/52 Section 2; docs/46 Amendment 2; docs/49 |
| 2026-10-07 | docs/45 Section 9 primary recommendation (MMRM unstructured KR2 at $t_{1/2}^a = 0.5$ week) superseded by F1 | F1 (`docs/49` Section 7; `docs/46` Amendment 2) | verified | docs/45 Section 9.1; docs/49 |
| 2026-10-07 | Three-part `bullets`/`rgt`/`orig` paragraph scaffolding dropped from all manuscripts; plain prose | Current; no `rgt` or `bullets` environments remain | inspected | CLAUDE.md; docs/30 Section 8 |
| 2026-10-07 | `R/carryover_analysis.R` threshold and predicted-effect formulas corrected (centered biomarker; with a mean-98 biomarker, old code was off by 7 points) | Fixed in code; test added | verified | docs/legacy/52 Section 4 |
| 2026-10-08 | Anchored carryover set as driver default in `04-power-simulation.R`, `08-covariance-study.R`, and `11-strawman-stress-test.R`; recursive carryover behind `--carryover recursive` | Current | inspected | docs/44 Section 4; `RERUN-anchored.md` |
| 2026-10-08 | Pre-exposure guard for `Dbc` corrected in `exposure_weight()`: pre-exposure off-drug visits now coded 0 instead of 1 | Fixed | inspected | docs/44 Section 5; docs/45 Section 8 |
| 2026-10-09 | Stress-test randomized path allocation set as driver default | Current | inspected | docs/44 Section 4 |
| 2026-10-09 | NOTATION.md revised: $\eta_t$ for design expectancy; $\rho_{cs}$ for compound-symmetry correlation; $\beta_T$ for linear-time coefficient | Current | inspected | docs/legacy/52; NOTATION.md |
| 2026-10-10 | Analysis half-life policy: $t_{1/2}^a = 0.2$ week for both interaction and main-effect fixed-half-life analyses; grid for data-chosen: {0, 0.2, 0.5, 1, 2} weeks | Current | decided | CLAUDE.md; docs/57 Section 8.2 |
| 2026-10-11 | DQ-1: Mathematics compendium becomes Volume II of the series; chapters cite CM labels for proofs | Current | decided | docs/59 Section 2.3 |
| 2026-10-11 | DQ-5: Analysis ID scheme finalized: `M-` prefix for main-effect track, `CV-` for covariance study, `PUB-` for published, paper-local prefix for collision-prone IDs | Current | decided | docs/59 Section 2.7 |
| 2026-10-11 | DQ-6: $\sigma_v^2$ renamed $\tau^2$ (between-patient variance) in docs/38, docs/51 and the lexicon | Current; tau^2 fix needed in appE | decided | docs/59 Section 7.2 |
| 2026-10-11 | DQ-8: Follow-up and $D_{bc,it}$ passes under anchored carryover cancelled for the covariance study | Stored passes use recursive carryover; anchored rerun of main pass only | decided | docs/45 Rerun section |
| 2026-10-11 | DQ-10: Half-life grid unified across docs/57 and docs/58: {0, 0.2, 0.5, 1, 2} weeks | Current | decided | docs/59 Section 7.2 |
| 2026-10-11 | Phase 1 conflict-resolution (CR-31, CR-43, CR-44, CR-51, GP-29): closed-form anchored figures verified; `characterize_carryover()` labeled F1A-type (liberal); published censordata form documented; failed-fit counts verified; four misattributed references corrected | Current | verified | docs/series/notes/phase1/phase1-cr-log.md |

Table: Chronological decision history (inferred from dated notes in
source documents; evidence label applies to the extraction, not
necessarily the original decision)

<!-- TODO: [add remaining dated decisions from docs/35 Sections 7
and 9 (recommendations timeline for separable AR(1)); docs/41
(baseline-covariate decision); docs/42 (flexible-time size-excess
decision); docs/43 (analysis-catalog status updates); docs/49
(F1 recommendation timeline); docs/56 (expectancy-term and CR2
decision); docs/57 (superseded plans for data-chosen half-life);
docs/legacy/50 (DGP and analysis status brief decisions)] -->

---

## Derivations

<!-- one row per result: statement | CM label | verification
script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the
numbers agree? -->

---

## Limitations and open items

- This table is constructed from inspected dated notes; not all
  decisions have been extracted. The sources docs/35, docs/41,
  docs/42, docs/43, docs/49, docs/56, docs/57, and docs/legacy/50
  contain additional dated notes that were not extracted in this
  pass.
- The "evidence" column describes how the decision was extracted
  here, not the original verification status.

---

## Reproduce

<!-- one command per output -->

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Initial decision table from dated notes in source documents | Phase 2 content move |
