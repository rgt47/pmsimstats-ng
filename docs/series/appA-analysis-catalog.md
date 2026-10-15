# Catalog of analyses and the ID registry

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Complete registry of all analysis IDs, with formula,
track, source and stored-data key.
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

| Symbol | Meaning |
|--------|---------|
| `id` | short analysis label (e.g. S, F1, M-S) |
| `track` | one of: interaction, main_effect, covariance_study, exposure_coding, published |

---

## Summary

This appendix is the machine-readable registry of all named analysis
IDs used in the series. The registry lives at
`docs/series/registry/analyses.csv`; this appendix renders it as a
human-readable table and adds the global-ID-scheme rules from
`docs/59` Section 2.7.

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- From docs/59 Section 2.7 -->

### Global ID scheme

The stored interaction-track labels (S, C1, C2, A, ES, EA, F1, F2,
R, DS, DA, G1, G2, F1A, F1A0, F1A0B, F1AB, MA, MaxT and their
suffixes) are keys in the stored stress-test files, and NOTATION
Part 2 forbids a display-time mapping layer. The scheme keeps them
and changes what is new or local. A document names the track at
first use ("interaction analysis F1",
"main-effect analysis `M-F1`").

| Class | Form | Examples | Change |
|---|---|---|---|
| Interaction track (stress test, covariance study selection, data-chosen half-life) | bare label | S, F1, F1A0B, G2 | none; the default track |
| Main-effect track (P04, `docs/58`) | `M-` prefix | `M-S`, `M-F1`, `M-VF`, `M-RV`, `M-F1A` | rename in P04, `docs/58` and `ch24` |
| Covariance study analyses (`docs/45`, `docs/43` Sections 3.7–3.8) | `CV-` prefix | `CV-V1` to `CV-V7`, `CV-F1` to `CV-F7` | removes collision with stress-test F1 |
| Paper-local catalog IDs | paper prefix | `02-G1` to `02-G9`, `06-C1` to `06-C4`, `08-T1`, `10-M0` | removes collisions with G1, G2, C1, C2 |
| Published analyses (`docs/43` Section 3.1) | `PUB-` prefix | `PUB-O1` to `PUB-O5` | none in substance |
| Summary measures and time adjustments | `SM-`, `TA-` prefix | `SM-1`, `TA-3` | removes collisions with S1, S2 and $D_{it}$ |
| Exposure codings | E1 to E9, E9b | unchanged | NOTATION extended E series |
| Constructs | always spelled out with $\rho$ | "configuration C at $\rho = 0.7$" | never a bare letter |

Table: The analysis-ID scheme for the series (decided 2026-10-11,
DQ-5; source: `docs/59` Section 2.7)

<!-- From docs/series/registry/analyses.csv -->

### Analysis ID registry

The table below is generated from
`docs/series/registry/analyses.csv`. Each row gives the analysis
label, its track, the model formula, the inference method, the
source documents, and any notes.

#### Interaction track

| ID | Formula | Inference | Source | Notes |
|---|---|---|---|---|
| S | Y ~ bmc + t + Dbc + bmc:Dbc; corCAR1; RI; E9 at fixed $t_{1/2}^a$ | CR2/Satterthwaite | docs/46,47; plan Sec 2.7 | Baseline interaction model; stress-test reference arm |
| C1 | Y ~ bmc + t + Dbc + bmc:Dbc; corCAR1; RI; compound-symmetry Dbc coding | CR2/Satterthwaite | docs/46,47 | Compound-symmetry exposure coding variant of S |
| C2 | Y ~ bmc + t + Dbc + bmc:Dbc; corCAR1; RI; compound-symmetry Dbc coding | Model-based (REML) | docs/46,47 | Model-based inference variant of C1 |
| A | Y ~ bmc + t + Dbc + bmc:Dbc; MMRM; E9 at fixed $t_{1/2}^a$ | Kenward-Roger (KR2) | docs/45,46,47 | MMRM with KR2 inference; from covariance study |
| ES | Y ~ bmc + t + Dbc + bmc:Dbc; corCAR1; RI; E9 exposure coding | CR2/Satterthwaite | docs/46,47 | E9 exposure coding, CR2 inference |
| EA | Y ~ bmc + t + Dbc + bmc:Dbc; MMRM; E9 exposure coding | Kenward-Roger (KR2) | docs/46,47 | E9 exposure coding, KR2 inference |
| F1 | Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De; corCAR1; RI; E9 at fixed $t_{1/2}^a$ = 0.2 week | CR2/Satterthwaite | docs/56,57; plan Sec 2.11 | Recommended analysis; expectancy terms |
| F2 | Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De; corCAR1; RI; E9 at fixed $t_{1/2}^a$ | Model-based (REML) | docs/46,47 | Model-based variant of F1 |
| R | Y ~ bmc + t + Dbc + bmc:Dbc; corCAR1; RI; ratio-graded coupling rule | CR2/Satterthwaite | docs/46,47 | Ratio-coupling exposure coding |
| DS | Y ~ bmc + t + Dbc + bmc:Dbc; corCAR1; RI; grid search over $t_{1/2}^a$ in {0, 0.25, 0.5, 1, 2} | CR2/Satterthwaite | docs/46,47 | Data-selected half-life; S formula |
| DA | Y ~ bmc + t + Dbc + bmc:Dbc; MMRM; grid search over $t_{1/2}^a$ in {0, 0.25, 0.5, 1, 2} | Kenward-Roger (KR2) | docs/46,47 | Data-selected half-life; A formula |
| G1 | Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De; corCAR1; $t_{1/2}^a$ chosen by AIC over null model | CR2/Satterthwaite | docs/47; plan Sec 2.7 | AIC-selected half-life without interaction in selection model |
| G2 | Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De; corCAR1; $t_{1/2}^a$ chosen by AIC over null model | Model-based (REML) | docs/47; plan Sec 2.7 | Model-based variant of G1 |
| F1A | F1 formula; $t_{1/2}^a$ chosen by AIC with interaction term; tested as if fixed | CR2/Satterthwaite | docs/57 Sec 1; plan Sec 2.7 | Liberal: AIC selection with interaction inflates size |
| F1A0 | F1 formula; $t_{1/2}^a$ chosen by AIC over model without interaction | CR2/Satterthwaite | docs/57 Secs 1,8; plan Sec 2.7 | AIC selection excludes interaction from model |
| F1A0B | F1 formula; $t_{1/2}^a$ chosen by AIC over null model; inference calibrated by parametric bootstrap | Bootstrap/CR2 | docs/57 Secs 1,8; plan Sec 2.11 | Provisional primary; bootstrap recalibrates size |
| F1AB | F1 formula; $t_{1/2}^a$ chosen by AIC with interaction; bootstrap calibration | Bootstrap/CR2 | docs/57 Sec 1; plan Sec 2.7 | Alternative data-driven candidate |
| MA | Max-t over interaction term across half-life grid | Permutation or asymptotic | docs/57; plan Sec 2.7 | Max-|t| aggregate over grid half-lives |
| MaxT | Max-t statistic over analysis half-life grid | Permutation | docs/57,58; plan Sec 2.7 | Maximum-t permutation test across half-life grid |
| F1h02 | F1 formula at fixed $t_{1/2}^a$ = 0.2 week | CR2/Satterthwaite | docs/57 Sec 8.2; plan Sec 2.11 | Pooled null rate 0.050 over 12 non-prognostic cells |
| F1h05 | F1 formula at fixed $t_{1/2}^a$ = 0.5 week | CR2/Satterthwaite | docs/47; stress test stored run | Historical stress-test base setting |

Table: Interaction-track analysis IDs with formula and inference
(source: `docs/series/registry/analyses.csv`)

#### Main-effect track

| ID | Formula | Inference | Source | Notes |
|---|---|---|---|---|
| M-S | Y ~ t + Dbc + De; corCAR1; RI; fixed $t_{1/2}^a$ = 0.2 week | CR2/Satterthwaite | docs/58; P04 protocol; plan Sec 2.7 | Standard main-effect test without biomarker terms |
| M-F1 | Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De; corCAR1; fixed $t_{1/2}^a$ = 0.2 week; test bmc:Dbc | CR2/Satterthwaite | docs/58; plan Sec 2.7 | Full F1 model; main-effect estimand extracted |
| M-VF | Vertically-flipped test; row-group weighted estimand | CR2/Satterthwaite | docs/58; plan Sec 2.7 | Addresses null bias from unknown half-life |
| M-RV | Randomization-test variant of main-effect analysis | Permutation | docs/58; plan Sec 2.7 | Permutation test for main effect |
| M-F1A | M-F1 formula; $t_{1/2}^a$ chosen by AIC over null model | CR2/Satterthwaite | docs/58; plan Sec 2.7 | Data-chosen half-life for main-effect track |

Table: Main-effect track analysis IDs (source: `analyses.csv`)

#### Covariance study track

| ID | Description | Inference | Source |
|---|---|---|---|
| CV-V1 through CV-V7 | Candidate CV analyses V1 to V7; see `docs/43` Section 3.7 | CR2/Satterthwaite | docs/45, docs/43 Sec 3.7 |
| CV-F1 through CV-F7 | Candidate CV analyses F1 to F7; see `docs/43` Section 3.8; the `CV-` prefix avoids collision with stress-test F1 | CR2/Satterthwaite | docs/45, docs/43 Sec 3.8 |

Table: Covariance study analysis IDs (source: `analyses.csv`)

#### Exposure-coding series

| ID | Description | Source |
|---|---|---|
| E1 | Binary Dbc: 1 on-drug, 0 off-drug and pre-exposure | CM E-series; docs/43 Sec 2; NOTATION |
| E2 | Binary Dbc with pre-exposure guard: 0 before first dose | CM E-series; docs/37 App A; NOTATION |
| E3 | Exponential decay off-drug; step to 1 on-drug | CM E-series; NOTATION |
| E4 | Exponential decay off-drug; step to 1 on-drug; pre-exposure guard | CM E-series; NOTATION |
| E5 | Compound-symmetry Dbc: constant off-drug value | CM E-series; NOTATION |
| E6 | AR(1) Dbc: correlated off-drug decay | CM E-series; NOTATION |
| E7 | Cumulative-exposure Dbc: integral of on-drug indicator | CM E-series; NOTATION |
| E8 | Phase-indicator Dbc: per-phase average | CM E-series; NOTATION |
| E9 | Exponential decay with pre-exposure guard and restart on re-exposure: exp(-lambda_cor * tsd) | CM E-series; docs/37; NOTATION |
| E9b | E9 with Dbc_it variant (time-varying biomarker-interaction term) | docs/45; CM E-series; NOTATION |

Table: Exposure-coding IDs (source: `analyses.csv`)

#### Published analyses

| ID | Description | Inference | Source |
|---|---|---|---|
| PUB-O1 | Hendrickson lmer: Y ~ bm + t + D + bm:D; RI; binary D_it | Model-based (lmer REML) | docs/43 Sec 3.1; docs/03,04 |
| PUB-O2 through PUB-O5 | Published analysis variants O2 to O5; see `docs/43` Section 3.1 | Model-based | docs/43 Sec 3.1 |

Table: Published analysis IDs (source: `analyses.csv`)

<!-- TODO: [extend coverage to docs/56-58 analyses not yet in the
CSV: confirm F1h02, F1h05, M-series are complete; add any
covariance-study analyses from docs/43 Secs 3.7-3.8 not yet in
the CSV as CV- entries] -->

<!-- TODO: [add paper-local catalog IDs from docs/43 (02-G1 to
02-G9, 06-C1 to 06-C4, 08-T1, 10-M0) with their original
meanings] -->

---

## Derivations

<!-- one row per result: statement | CM label | verification
script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

This appendix is the authoritative ID source for all manuscripts.
Each manuscript that names an analysis ID must use the ID from this
registry without a separate display-time mapping.

<!-- TODO: [table of which manuscript sections depend on which IDs;
confirm the numbers agree with the stored CSV] -->

---

## Limitations and open items

- The CSV covers the interaction, main-effect, covariance-study,
  exposure-coding and published tracks. Paper-local IDs (02-G1 to
  02-G9, etc.) are specified in the ID scheme but not yet rows in
  the CSV.
- Coverage of docs/56-58 analyses requires verification.

---

## Reproduce

The CSV is the machine-readable source:

```bash
cat docs/series/registry/analyses.csv
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Initial content move: analyses.csv expanded to tables; global ID scheme from docs/59 Sec 2.7 | Phase 2 content move |
