# Errata to earlier documents and manuscripts

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Errata to earlier documents and manuscripts with status
of each correction.
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

This appendix collects corrections to earlier project documents and
manuscripts. Each entry gives the document, the section, what was
wrong, the correction, and its current status.

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 36 Section 8; 40 Section 11; 45 Sections 8 and 10; CM
  main-effect section.
- New to write: the status of each item.
-->

### Corrections from docs/36: Hendrickson PD and carryover

<!-- From docs/36 Section 8 -->

The following corrections apply to the documents listed in
`docs/36` Section 8 ("Relation to earlier project documents"):

| Document | Statement corrected | Correction | Status |
|---|---|---|---|
| `docs/07` | 40 PD failures in 162 checks; breakdown 32/40 in OL+BDC, 2/24 in Hybrid; explanation via expectancy scaling of placebo variance | The breakdown and explanation do not reproduce. Positive diagonal scaling cannot change positive definiteness; the published code's own matrices fail as in docs/36 Section 4.1 | Superseded by docs/36 Section 4 |
| `docs/07a` | Empirical ceiling of ~0.6 under an AR(1) structure the published code does not use | Superseded | Superseded |
| `docs/09` | Correctly identifies the step rule but attributes power loss to standard-error inflation with an unbiased estimate; arithmetic slips in residual fractions | The published estimates attenuate by ~three quarters; arithmetic corrected | Superseded by docs/36 Section 5.3 |
| `docs/13` | Labels code as the initial commit but describes the 2024 revision's coupling, scale factor and analysis variable | Corrected: the published process is commit 58b32a9 | Fixed in docs/13 banner |
| Manuscript 02 | Exposure-weighted analysis found inferior in CO; attributed to carryover | Likely an artifact of the pre-exposure coding: CO's placebo-first path got $D_{bc}=1$ at every visit | Pre-exposure guard fixed in `exposure_weight()`; CO cells not yet rerun (open) |

Table: Corrections from docs/36 Section 8 (inspected)

### Corrections from docs/40: Hendrickson hybrid design primer

<!-- From docs/40 Section 11 -->

`docs/40` Section 11 ("Corrections to earlier documents") corrects
earlier documents as follows:

| Document | Statement corrected | Correction | Status |
|---|---|---|---|
| `docs/24` | Gompertz parameter meaning reversed | Corrected in docs/40 | Fixed |
| `docs/24` | Hybrid design misdescribed | Corrected | Fixed |
| `docs/24` | Variance identity | Corrected | Fixed |
| `docs/25` | Parameter table (TV row) wrong | Corrected | Fixed |
| `docs/25` | Heterogeneity claim withdrawn | Withdrawn | Fixed |
| `docs/28` | Crossover row and power pair corrected | Corrected | Fixed |
| Multiple | Figure 2 figures 40-fig2 and 40-fig3 numbered opposite to print order | Open | Open (`docs/legacy/52` Section 4) |

Table: Corrections from docs/40 Section 11 (inspected)

`docs/40` also records an update (2026-10-09): results at $t_{1/2}>0$
now use anchored carryover; values changed by at most 0.008 in
rejection rates and 0.53 points in the drug coefficient. Conclusions
unchanged.

### Corrections from docs/45: correctly sized hybrid analysis

<!-- From docs/45 Sections 8 and 10 -->

`docs/45` Section 8 ("The pre-exposure coding of $D_{bc,it}$")
identified an error in manuscript 02's CO cells:

| Document | Statement corrected | Correction | Status |
|---|---|---|---|
| Manuscript 02, Section 3.4 and blocks S7-S10 | Exposure-weighted analysis collapses under CO | Likely an artifact of the pre-exposure coding; CO power was 0.484 vs 0.744 under the corrected coding (CS) | Fixed in driver; manuscript CO cells not yet rerun (open) |

`docs/45` Section 10 ("Corrections to earlier documents"):

| Document | Statement corrected | Correction | Status |
|---|---|---|---|
| `docs/41`, Section 7 table | cLDA and ANCOVA coincide for the unstructured row model | True only with $\beta_B$ free at baseline; the published row coding shares it and differs (docs/45 Section 6) | Fixed |
| `docs/41`, Section 8 | Recommends baseline as a covariate for the trial analysis | Superseded by docs/45 Section 9: the covariate form is the sensitivity analysis, or primary if constancy cannot be defended | Superseded |
| `docs/42` and `docs/43`, Section 4 | Flexible-time size excess persists under every covariance | Confined to one cell, absent under CR2 (docs/45 Section 4.4); real under model-based and KR inference | Fixed in docs/42 and docs/43 |
| `docs/43`, Sections 1 and 5 | `us-KR` sized at 0.046 to 0.072 | Slightly liberal, mean 0.058; CR2 analyses with parsimonious covariances are nominal | Fixed |
| `docs/43`, Section 6 | Follow-up questions open | Answered in docs/45 Sections 4 to 7 | Fixed |
| Manuscript 02, Section 3.4 | CO pre-exposure coding error | See above | Open |

Table: Corrections from docs/45 Section 10 (inspected)

### Corrections from the mathematics compendium (main-effect section)

<!-- TODO: [search compendium for errata/erratum/correction in the
main-effect section; extract any corrections to earlier documents
stated there] -->

### Phase 1 corrections (GP-29): misattributed references

<!-- From docs/series/notes/phase1/phase1-cr-log.md, GP-29 -->

Four references in `docs/legacy/20` and `docs/legacy/23` were
misattributed. Corrected in `docs/legacy/23` Section 17 (inspected
2026-10-07) and in the bibliography:

| Cited as | Correct citation | Status |
|---|---|---|
| Chen Y, Li P, Yang K. *PLoS One* 2020;15(1):e0216949 | Blackston JW, Chapple AG, McGree JM, McDonald S, Nikles J. *Healthcare* 2019;7(4):137, doi:10.3390/healthcare7040137 | Corrected in docs/legacy/23 Section 17 and bib |
| Eng KH, Konings P. *Trials* 2018;19:72, doi:10.1186/s13063-018-2491-0 | Haller B, Ulm K. *Trials* 2018;19:128, doi:10.1186/s13063-018-2491-0 | Corrected |
| Pequignat F, Balaam B, et al. *Behav Res Methods* 2023;55:3824 | Tueller S, Ramirez D, et al. *Behav Res Methods* 2023;55(8):4175-4199, doi:10.3758/s13428-022-02012-1; `analysis/report/references.bib` key `pequignatIdiographicPower2023` already corrected | Corrected in bib |
| De Bock M, Van den Noortgate W, Gielen M. *Behav Res Methods* 2021;53:1782 | Natesan Batley P, Hedges LV. *Behav Res Methods* 2021;53:1782-1798, doi:10.3758/s13428-020-01522-0 | Corrected |

Table: GP-29 misattributed references corrected in Phase 1
(verified, phase1-cr-log.md)

Note: the annotations in `docs/legacy/20` describe content
attributed to invented authors; they are not carried into appG
without re-reading the real papers.

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

- Status column entries are as of 2026-10-15 (inspected from source
  documents); some "Fixed" items may have been re-introduced by
  subsequent edits and have not been re-verified here.
- The mathematics compendium main-effect section has not been
  searched for errata in this pass.
- Manuscript 02's CO cells have not been rerun with the corrected
  pre-exposure coding (open item).

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
| 2026-10-15 | Setting and argument | Initial errata tables from docs/36 Sec 8, docs/40 Sec 11, docs/45 Secs 8 and 10, GP-29 | Phase 2 content move |
