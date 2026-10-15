# Implementation status and open-items register

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Implementation status and open-items register with one
maintained inspection date.
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

This appendix is the single maintained register of all open
implementation items. **Last inspection date: 2026-10-15.** Items
are tagged by status: Open, Partial, Fixed, or Deferred.

---

## Construct box

*No Monte Carlo simulation results in this chapter.*

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 35 Sections 7, 9; 44 Section 4; 52 Section 4; 23
  Section 14; 34 Sections 4.7, 6.2; 36 Section 9.
- New to write: one maintained register with one inspection date.
-->

### Open-items register (inspection date: 2026-10-15)

<!-- From docs/35 Sections 7, 9 (separable response covariance)
     and docs/44 Section 4 (DGP roadmap status) -->

#### Package: DGP and simulation core

| ID | Item | Status | Source | Notes |
|---|---|---|---|---|
| PKG-01 | Add separable AR(1) form to `buildSigma()` as a new cross-component option | Open | docs/35 Sec 9; docs/44 Sec 4 | Step 1 of the refined-process adoption; no implementation as of 2026-10-09 (inspected) |
| PKG-02 | Add factored generator to `generateSimulatedResults()` as the sampler for the separable form; stops on infeasible $c_{bm}$ | Open | docs/35 Sec 9; docs/44 Sec 4 | Step 2; removes the silent Cholesky repair |
| PKG-03 | Remove (or make into an error) the Cholesky fallback in `R/generateData.R` lines 391-392 | Open | docs/35 Sec 10; docs/44 Sec 4; docs/legacy/52 Sec 4 | The fallback still repairs a non-PD matrix silently (inspected 2026-10-07) |
| PKG-04 | Add `cbmCeiling()` implementing the closed-form coupling ceiling; call from `generateSimulatedResults()` and `validateParameterGrid()` | Open | docs/35 Sec 9 | |
| PKG-05 | Apply separable AR(1) change to `implementations/tidyverse` and `implementations/nof1power`, or retire them in favor of the package | Open | docs/35 Sec 9 | |
| PKG-06 | Parity: regenerate parity baselines and rerun affected papers after PKG-01 to PKG-03 | Open | docs/35 Sec 9 | Depends on PKG-01 to PKG-05 |

#### Package: analysis and inference

| ID | Item | Status | Source | Notes |
|---|---|---|---|---|
| ANA-01 | `R/carryover_analysis.R`: fix `characterize_carryover()` to use ML (not REML) likelihoods when comparing models with different fixed effects, and F1A0 selection rule (select without interaction); label output as descriptive | Open | docs/series/notes/phase1/phase1-cr-log.md, CR-43 | Current code is F1A-type (liberal); no driver calls it |
| ANA-02 | `R/carryover_analysis.R`: add 0 to the half-life grid, align with the DQ-10 grid {0, 0.2, 0.5, 1, 2}; guard against all-failed fits before the AIC refit | Open | docs/series/notes/phase1/phase1-cr-log.md, CR-43 | |
| ANA-03 | `R/carryover_analysis.R`: update `print_carryover_summary()` wording to not present selected-half-life p-value as confirmatory | Open | docs/series/notes/phase1/phase1-cr-log.md, CR-43 | |
| ANA-04 | `R/lme_analysis.R`: singularity flag (`issingular`) is hardcoded FALSE; detect from `fit$apVar` | Partial | docs/legacy/23 Sec 14 (C2); docs/legacy/52 Sec 4 | VarCorr check added in root package (2026-10-07, inspected); original audit finding C2 |

#### Package: code quality

| ID | Item | Status | Source | Notes |
|---|---|---|---|---|
| CQ-01 | `R/censordata.R` and `R/plottingfunctions.R`: replace `eval(parse(text=...))` with `get()` or `dat[[name]]` | Open | docs/legacy/23 Sec 14 (C3) | No behavior change |
| CQ-02 | `R/lme_analysis.R`: positional coefficient extraction; harmonize around `carryover_analysis` pattern | Open | docs/legacy/23 Sec 14 (C4) | |
| CQ-03 | DESCRIPTION: remove stale Imports (`lme4`, `lmerTest`, `ggpubr`, `gridExtra`) | Open | docs/legacy/23 Sec 14 (C6) | All four still listed under Imports (inspected 2026-10-07) |
| CQ-04 | Tidyverse `functions.R` lines 421-423: qualify `filter()` and `pull()` calls with `dplyr::` | Open | docs/legacy/23 Sec 14 (T1) | CRAN blocker; parity test masks issue because it runs in a tidyverse-loaded session |

#### Testing

| ID | Item | Status | Source | Notes |
|---|---|---|---|---|
| TST-01 | Add tinytest case for `dgp_architecture = 'mean_moderation'` without tidyverse loaded | Open | docs/legacy/23 Sec 14 (C7) | |
| TST-02 | Add tinytest case for null-effect Type I error at nominal alpha | Open | docs/legacy/23 Sec 14 (C7) | |
| TST-03 | Add tinytest case for `corCAR1` correctness under known `phi` ground truth | Open | docs/legacy/23 Sec 14 (C7) | |
| TST-04 | Fix mismatched parity reports: `parity_report.rds` is quick mode (16 rows), `parity_baseline.rds` is default mode (144 rows) | Open | docs/legacy/52 Sec 4 | `parity-diff.R` compares unlike reports |

#### Analysis/scripts

| ID | Item | Status | Source | Notes |
|---|---|---|---|---|
| SCR-01 | Manuscript 02: rerun CO cells with corrected pre-exposure coding | Open | docs/45 Sec 8; docs/legacy/52 Sec 4 | Not yet confirmed |
| SCR-02 | `analysis/scripts/carryover-sensitivity/simulation-core.R`: comment says TV parameters follow published defaults; they do not | Open | docs/legacy/52 Sec 4 | |
| SCR-03 | `docs/21`: tables cannot be reproduced; driver output file missing; rerun or retire tables in favor of paper 01 | Open | docs/legacy/52 Sec 4 | |
| SCR-04 | `docs/34`: $\rho=0.8$ missing from `rho.csv` and Figure 2; Figure 10 subtitle "binding paths" wrong for Hybrid | Open | docs/legacy/52 Sec 4 | |
| SCR-05 | `docs/40`: figure files `40-fig2` and `40-fig3` numbered opposite to print order | Open | docs/legacy/52 Sec 4 | |

<!-- From docs/36 Section 9 -->

#### DGP open items (from docs/36 Section 9)

| ID | Item | Status | Source | Notes |
|---|---|---|---|---|
| DGP-01 | Continuous exposure in analysis: test `Dbc` (vs binary) under configuration C at varying half-lives; under mean moderation, test that it introduces paper 01's mismatch | Open | docs/36 Sec 9 | Prediction; not yet tested |
| DGP-02 | Add $\rho = 0.7$ to configuration C in docs/36 power simulation | Open | docs/36 Sec 9 | Current simulation used published $\rho = 0.8$ |
| DGP-03 | Run smaller samples ($N = 35$) and multi-cycle designs for PD ceiling comparisons | Open | docs/36 Sec 9 | Only $N = 70$ and three published designs simulated |
| DGP-04 | Explain elevated Type I error in OL+BDC under compound-symmetry configurations (0.06-0.11 range) | Open | docs/36 Sec 9 | Neither analysis model accounts for this |
| DGP-05 | Integrate the factored generator (Section 6.7 of docs/36) into `buildSigma()`; confirm it reproduces configuration C's power | Open | docs/36 Sec 9 | Overlaps with PKG-01 to PKG-02 |

<!-- From docs/44 Section 4 (longer-term) -->

#### Longer-term: docs/38 model

| ID | Item | Status | Source | Notes |
|---|---|---|---|---|
| LT-01 | Implement `docs/38` PK/PD model as a new `dgp_architecture` value after step-8 checks; run alongside the refined process | Deferred | docs/44 Sec 4 | Longer-term target |

<!-- From docs/34 Sections 4.7 and 6.2 (cbm feasible range) -->

#### Coupling ceiling: outstanding items

<!-- TODO: [Read docs/34 Sections 4.7 and 6.2 and extract their
specific open items] -->

<!-- TODO: [docs/34 Section 4.7: behavior of the PD repair at
varying rho and carryover; Section 6.2: project-specific action
items for the coupling ceiling] -->

---

## Package implementation status

Last inspection date: 2026-10-15. Status values: Implemented, Open,
Partial, Deferred.

| Feature | Status | Notes |
|---|---|---|
| Separable AR(1) covariance matrix | Open | `buildSigma()` uses the published compound-symmetry construct by default; PKG-01 |
| Graded coupling profile ($c_{bm} \cdot e^{-\lambda t_{sd}}$ off drug) | Implemented | Encoded in Sigma matrix for on/off-drug timepoints in `generateData.R` |
| Anchored carryover | Implemented | Default in drivers since 2026-10-08; `--carryover recursive` reproduces stored runs |
| Pre-exposure guard on $D_{bc}$ ($D_{bc} = 0$ before first dose) | Implemented | Guarded in `lme_analysis.R` |
| Randomized path allocation | Implemented | `randomize_paths` option in `buildtrialdesign()` |
| Stop on non-positive-definite matrix | Open | Package still silently repairs non-PD matrix via Cholesky fallback; PKG-03 (`R/generateData.R` known open item) |
| Factored generator (step 5 of refined process) | Open | Efficiency improvement, not blocking correctness; PKG-02 |
| `cbmCeiling()` coupling-ceiling function | Open | Not yet exported from package; PKG-04 |
| Prognostic biomarker couplings | Implemented | `blparam` list accepts `c.bm.pb` for contaminated-biomarker DGP |
| Mean-moderation architecture | Implemented | `dgp_architecture = 'mean_moderation'` in `generateData()` |
| Dual-channel architecture | Implemented | `dgp_architecture = 'combined'` with `c.bm_a` and `c.bm_b` |

---

## Open items register

Computational gaps from `docs/59` Section 8 and implementation open
items. Status: open, in progress, deferred, or blocked.

| Item | Chapter | Gap | Status |
|---|---|---|---|
| F1A0B full bootstrap check ($M=999$, $\approx$170 worker-hours) | ch23 | GP-15 | deferred |
| Sample-size planning for F1 and F1A0B | ch16 | GP-16 | deferred |
| CO and OL+BDC interaction analysis | ch02, ch17 | GP-17 | deferred |
| Prognostic coupling strengths other than 0.3 | ch19 | GP-18 | deferred |
| Features 11-17 of docs/51 simulation | ch11 | GP-19 | deferred |
| Calibration to prazosin trial data | ch10 | GP-20 | blocked (data access) |
| Main-effect stress test under rebuilt P04 | ch24 | GP-21 | in progress (per P04 Section 7.3) |
| Coverage of profile-likelihood intervals | ch18 | GP-22 | deferred |
| Design ranking under refined process at matched $N$ | ch16 | GP-25 | deferred |
| Paper 02 CO rerun with corrected pre-exposure coding | appC | GP-31 | open (computation) |
| OL+BDC Type I error under compound symmetry | ch04 | GP-32 | open (analysis) |
| Paper 11 boundary cells and super-additivity proof | ch07 | GP-33 | open (computation) |
| F1 at 0.2 week in 41-cell stress test | ch21 | GP-39 | deferred (optional) |
| Silent Cholesky repair in `generateData.R` | appD | PKG-03 | open |
| Factored generator in `generateSimulatedResults()` | appD | PKG-02 | open |
| `cbmCeiling()` exported function | appD | PKG-04 | open |
| Parity baselines after refined-process adoption | appD | PKG-06 | open (depends on PKG-01 to PKG-05) |

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

- docs/34 Sections 4.7 and 6.2 were not extracted in this pass
  and require a separate read.
- Status labels are as of 2026-10-15; items flagged Fixed in the
  source documents may have been re-opened by subsequent edits.
- "Deferred" means the item is acknowledged but not scheduled.

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
| 2026-10-15 | Setting and argument | Initial register from docs/35 Secs 7,9; docs/44 Sec 4; docs/legacy/52 Sec 4; docs/legacy/23 Sec 14; docs/36 Sec 9 | Phase 2 content move |
| 2026-10-15 | Package implementation status, Open items register | GP-27 and GP-38: implementation status table and computational gap register | Phase 3 gap fill |
