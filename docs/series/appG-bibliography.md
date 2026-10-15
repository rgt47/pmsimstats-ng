# Bibliography

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Curated bibliography with verification status for
suspect entries.
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

This appendix is the curated bibliography for the series. It contains
three things:

1. **Literature synthesis** from `docs/legacy/23` Section 17 (50
   references organized in 7 subsections; 7 entries corrected from
   misattributed originals on 2026-10-07).
2. **Phase 1 corrections table** (GP-29): four misattributed
   references in `docs/legacy/20` and `docs/legacy/23` that were
   corrected during the Phase 1 review.
3. **Pointer to the BibTeX file** `docs/references/nof1-simulation.bib`
   (the machine-readable source).

Notes in `docs/legacy/20` describe content attributed to invented
authors; those annotations are not reproduced here. The machine-
readable bibliography is `docs/references/nof1-simulation.bib`.

---

## Construct box

*No Monte Carlo simulation results in this chapter.*

---

## Setting and argument

<!-- From docs/legacy/23-documentation-compendium.md Section 17 -->

*Note (2026-10-07).* Each entry in the literature synthesis below was
checked against Crossref and PubMed. Seven were misattributed or had
wrong metadata and have been corrected in place, with the earlier
citation recorded. The same errors probably occur in
`docs/20-annotated-bibliography.tex`, which was not checked. DOIs are
given for verified entries.

### 17.1 Core N-of-1 Simulation Methodology

- **Hendrickson, Thomas, Schork, and Raskind (2020)**, Frontiers in
  Digital Health 2:13, doi:10.3389/fdgth.2020.00013. Simulation study
  (four designs, LME analysis, 1,000 repetitions) of aggregated
  N-of-1 designs for detecting predictive biomarker effects. Finding:
  the N-of-1 (Hybrid) design outperformed OL and OL+BDC but had
  slightly lower power than the traditional crossover. This is the
  framework pmsimstats-ng revises. (Corrected 2026-10-07: previously
  cited as "Thomas et al.", described as deriving analytic power
  formulae, and summarized as needing fewer participants than
  parallel-group designs, none of which matches the abstract.)
- **Blackston et al. (2019)**, Healthcare 7(4):137,
  doi:10.3390/healthcare7040137. 5,000 simulated trials per scenario
  across four variance structures comparing aggregated N-of-1,
  parallel RCT, crossover. Finding: aggregated N-of-1 outperforms
  traditional designs but has elevated Type I error under moderate-to-
  strong unmodeled carryover.
- **Wang and Schork (2019)**, Healthcare 7(3):84,
  doi:10.3390/healthcare7030084. Analytical power formulae under AR(1)
  serial correlation. Finding: serial correlation substantially reduces
  power; increasing evaluation periods partially compensates.

### 17.2 Carryover and Temporal Correlation

- **Shi and Ye (2024)**, "Behavioral carry-over effect and power
  consideration in crossover trials," Biometrics 80(2):ujae023,
  doi:10.1093/biomtc/ujae023 (corrected 2026-10-07; previously
  attributed to "Yap, Wang, and Harhay"). Formal framework for
  behavioral carryover in crossover trials (psychological mechanisms,
  not pharmacokinetic). Closest comparator for the pmsimstats
  manuscript. Finding: behavioral carryover biases the basic estimator
  toward underestimating the treatment effect, costing power.
  Directly relevant to the prazosin-PTSD application.
- **Sturdevant and Lumley (2021)**, Contemporary Clinical Trials
  Communications 22:100711, doi:10.1016/j.conctc.2021.100711.
  Mixed-effects model approach for testing carryover. Finding:
  'carryover has been largely ignored except as a nuisance parameter'.
- **Schork (2022)**, Harvard Data Science Review Special Issue 3,
  doi:10.1162/99608f92.f1eef6f4. Serial correlation in personalized
  and aggregated studies. Finding: ignoring serial correlation inflates
  false positives.
- **Tang and Landes (2020)**, PLoS One 15(2):e0228077,
  doi:10.1371/journal.pone.0228077. Modified t-tests accounting for
  serial correlation in N-of-1 data. Finding: proposed tests maintain
  nominal Type I error where standard t-tests fail.

### 17.3 Biomarker-Treatment Interactions

- **Haller and Ulm (2018)**, Trials 19:128,
  doi:10.1186/s13063-018-2491-0. Simulation study of biomarker-
  treatment interaction estimation in randomized trials with time-to-
  event outcomes.
- **Grenet et al. (2020)**, Basic and Clinical Pharmacology and
  Toxicology 126(1):59-64, doi:10.1111/bcpt.13293. Simulation
  comparing crossover versus parallel-group for binary predictive
  biomarker identification. Finding: crossover more efficient when
  within-subject correlation is high.

### 17.4 Analysis Methods

- **Chen and Chen (2014)**, "A comparison of four methods for the
  analysis of N-of-1 trials," PLoS One 9(2):e87752,
  doi:10.1371/journal.pone.0087752 (corrected 2026-10-07; previously
  attributed to "Senn, Julious, and Araujo"; the Araujo, Julious, and
  Senn paper is a different one, PLoS One 2016, 11(12):e0167167).
  Four methods for N-of-1 analysis (fixed-effects, random-effects,
  Bayesian hierarchical, paired). Finding: all similar for population
  inference; differ substantially for individual treatment effect
  estimation.
- **Peterlin, Kejzar, and Blagus (2023)**, TEST 32(1):184-210,
  doi:10.1007/s11749-022-00830-1. Diagnostic tests for correct
  specification of fixed- and random-effects design matrices in LMM.

### 17.5 Sample Size

- **Senn (2019)**, Statistical Methods in Medical Research
  28(2):372-383, doi:10.1177/0962280217726801. Five distinct criteria
  for planning sample size in N-of-1 trials. Finding: criterion
  depends on goal (individual-level, population-level, heterogeneity
  detection).
- **Schmid and Duan (2014)**. Statistical design and analytic
  considerations for N-of-1 trials. Chapter 4 in Kravitz RL, Duan N
  (eds.), *Design and Implementation of N-of-1 Trials: A User's
  Guide*, AHRQ Publication No. 13(14)-EHC122-EF, Rockville, MD:
  AHRQ, 2014. (Chapter authorship confirmed from AHRQ product page;
  automated access blocked at AHRQ 2026-10-07.)
- **Hawksworth et al. (2024)**, "A methodological review of randomised
  n-of-1 trials," Trials 25:263, doi:10.1186/s13063-024-08100-1
  (corrected 2026-10-07; previously cited as "Senn, Julious, and Chen
  (2024), Trials 25:264"). Methodological review of published
  randomized N-of-1 trials. Finding: widespread inconsistencies in
  reporting and analytic practice.

### 17.6 Bayesian Approaches

- **Schmid and Yang (2022)**, "Bayesian models for N-of-1 trials,"
  Harvard Data Science Review, doi:10.1162/99608f92.3f1772ce
  (corrected 2026-10-07; previously "Schmid and Staudenmayer
  (2021)"). Review of Bayesian modeling for N-of-1 trials.
- **Daza (2018)**, "Causal analysis of self-tracked time series data
  using a counterfactual framework for N-of-1 trials," Methods of
  Information in Medicine 57(1):e10-e21,
  doi:10.3414/ME16-02-0044 (corrected 2026-10-07; previously cited
  as Daza, Wac, and Oppezzo, AMIA Annual Symposium Proceedings
  2018:532-541, which was not found). Counterfactual framework for
  causal inference in N-of-1 trials via self-tracked time series.

### 17.7 Prazosin-PTSD Clinical Context

- **Raskind et al. (2003)**, American Journal of Psychiatry
  160(2):371-373, doi:10.1176/appi.ajp.160.2.371. Initial
  proof-of-concept (N=10), double-blind crossover; significant
  reduction in nightmares.
- **Raskind et al. (2013)**, American Journal of Psychiatry
  170(9):1003-1010, doi:10.1176/appi.ajp.2013.12081133.
  Parallel-group RCT (N=67) in active-duty soldiers; significant
  effects; observed effect sizes and heterogeneity calibrate
  simulations.
- **Raskind et al. (2018)**, New England Journal of Medicine
  378(6):507-517, doi:10.1056/NEJMoa1507598. PACT trial (N=304); no
  significant difference versus placebo on primary outcome. The null
  result motivates the pmsimstats framework search for a predictive
  biomarker.
- **Mayer et al. (2023)**, "Randomized controlled pilot trial of
  prazosin for prophylaxis of posttraumatic headaches in active-duty
  service members and veterans," Headache 63(6):751-762,
  doi:10.1111/head.14529 (corrected 2026-10-07; previously attributed
  to Raskind et al.). Pilot for posttraumatic headache prophylaxis
  (N=48; sample size not re-checked).

---

### Phase 1 corrections (GP-29): misattributed references

<!-- From docs/series/notes/phase1/phase1-cr-log.md, GP-29;
     also in Appendix C -->

Four references in `docs/legacy/20` and `docs/legacy/23` were
misattributed. Corrected in `docs/legacy/23` Section 17 (inspected
2026-10-07) and in the bibliography:

| Cited as | Correct citation | Status |
|---|---|---|
| Chen Y, Li P, Yang K. *PLoS One* 2020;15(1):e0216949 | Blackston JW, Chapple AG, McGree JM, McDonald S, Nikles J. *Healthcare* 2019;7(4):137, doi:10.3390/healthcare7040137 | Corrected in `docs/legacy/23` Section 17 and bib |
| Eng KH, Konings P. *Trials* 2018;19:72 | Haller B, Ulm K. *Trials* 2018;19:128, doi:10.1186/s13063-018-2491-0 | Corrected |
| Pequignat F, Balaam B, et al. *Behav Res Methods* 2023;55:3824 | Tueller S, Ramirez D, et al. *Behav Res Methods* 2023;55(8):4175-4199, doi:10.3758/s13428-022-02012-1; `analysis/report/references.bib` key `pequignatIdiographicPower2023` already corrected | Corrected in bib |
| De Bock M, Van den Noortgate W, Gielen M. *Behav Res Methods* 2021;53:1782 | Natesan Batley P, Hedges LV. *Behav Res Methods* 2021;53:1782-1798, doi:10.3758/s13428-020-01522-0 | Corrected |

Table: GP-29 misattributed references corrected in Phase 1 (verified,
`phase1-cr-log.md`)

---

### BibTeX file pointer

The machine-readable bibliography is:

```
docs/references/nof1-simulation.bib
```

This file contains the full reference list for the series. Per-chapter
reference lists are generated from it at render time. The shared
candidate bibliography for manuscripts is:

```
analysis/report/references.bib
```

<!-- TODO: [verify the suspect entries listed in Appendix I Section 11
against `docs/references/nof1-simulation.bib`; compare to the
corrected versions in docs/legacy/23 Section 17; confirm that the
GP-29 corrections have been applied in the .bib files and not only
in the prose] -->

---

## Derivations

<!-- one row per result: statement | CM label | verification
script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

All manuscript reference lists draw from `analysis/report/references.bib`
and `docs/references/nof1-simulation.bib`. The GP-29 corrections apply
to both files.

---

## Limitations and open items

- The annotations in `docs/legacy/20` describe content attributed to
  invented authors; they are not reproduced here and require reading
  the real papers before any summary can be written.
- Verification of suspect entries listed in Appendix I Section 11
  against the `.bib` files has not been done in this pass (see TODO
  above).
- The coverage of the literature synthesis is the 50 references mined
  from the Zotero `nof1` collection; more recent publications are not
  included.

---

## Reproduce

```bash
cat docs/references/nof1-simulation.bib
```

```bash
cat analysis/report/references.bib
```

---

## References

See `docs/references/nof1-simulation.bib`.

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Content move from docs/legacy/23 Section 17 and GP-29 corrections table; pointer to .bib file | Phase 2 content move |
