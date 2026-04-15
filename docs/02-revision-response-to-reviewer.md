---
title: "Revision response to reviewer comments on doc 02"
author: "pmsimstats team"
date: "2026-04-15"
---

# Revision response: reviewer comments on `02-dgp-mean-moderation-vs-mvn.tex`

This document summarises the edits applied to
`02-dgp-mean-moderation-vs-mvn.tex` in response to reviewer RH's
annotated PDF (`02-dgp-mean-moderation-vs-mvn rchcomments_rh.pdf`)
and the accompanying email summary. Edits were applied in two
passes: an initial surgical pass addressing the line-level sticky
notes, followed by a second pass addressing the two higher-level
concerns raised in the email summary.

## Pass 1: edits from sticky notes

| # | Comment                                                                       | Location           | Edit                                                                                                                                                                                                       |
|---|-------------------------------------------------------------------------------|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1 | "whether" -> "the extent to which"                                            | Abstract           | Reworded                                                                                                                                                                                                   |
| 2 | Why the "MVN" abbreviation?                                                   | Abstract           | Expanded to "multivariate normal (MVN)" on first use                                                                                                                                                       |
| 3 | "i.e., when on drug" is inaccurate; moderation also persists during decay     | Section 2.2        | Replaced with the attenuated off-drug correlation formula and noted the moderation persists with decay                                                                                                     |
| 4 | Does carryover asymmetry really apply only to B? (raised three times)         | Section 3.2        | Explicitly labelled mean-blurring as a shared channel across both architectures; identified correlation erosion as the B-unique channel; added explanation of why Architecture A ratios survive amplitude loss |
| 5 | "Deterministic" too strong given biological variability                       | Section 4.1        | Clarified that determinism refers to the population mean structure, not individual outcomes                                                                                                                |
| 6 | "(in preparation)" placeholder citation                                       | Section 4.3        | Replaced with "ongoing reanalysis by the pmsimstats team"                                                                                                                                                  |
| 7 | Section 6.1 lists only effective-sample-size cost; confounding cost missing   | Section 6.1        | Added placebo and time-varying-response confounding as the dominant cost; explicitly linked to the observation that OL does not dominate in Section 3                                                      |
| 8 | Section 6.6 "Front-loaded" rationale unclear                                  | Section 6.6        | Explained the trade-off versus the reverse ordering (off-drug baseline first)                                                                                                                              |

## Pass 2: edits from email summary

| # | Comment                                                                                                                                                 | Location    | Edit                                                                                                                                                                                                                                                                                                                                                |
|---|---------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 9 | Is the divergent carryover impact an artifact of Architecture A being implemented without the non-biologic components?                                  | Section 3.1 | Added an explicit methodological paragraph: TV, placebo, AR(1) autocorrelation, and cross-factor correlations come from the same MVN draw in both architectures. The DGPs differ only in BM-BR linkage (covariance entry in B, post-draw additive shift in A). Cross-referenced the scaling derivation in document 19.                              |
| 10| If Architecture B is really "imperfect responder indicator," why not model that explicitly at the participant level rather than via covariance structure? | Section 4.2 | Added a paragraph acknowledging the tension: if the biomarker is truly a noisy latent-class indicator, an explicit finite-mixture DGP is the biologically faithful representation. MVN differential correlation is best understood as a tractable second-moment approximation. Framed mixture modelling versus analysis-strategy mitigation as an open question. |

## Skipped with rationale

| Comment                                                                             | Rationale                                                                                                                   |
|-------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| Section 2.2 BR vocabulary needs upfront introduction                                | Already handled in the Section 2 "Notation" paragraph                                                                       |
| Section 2.2 properties list should be parallel with Section 2.1, ideally a table    | Reviewer later retracted at Section 5.1: "Ah, just what I was suggesting above! Great!"                                     |
| Section 6.5 regrouping with Section 6.1 / Section 6.2                               | Structural reorganisation beyond the scope of targeted revision                                                             |
| Section 7.1 predictive and prognostic framing feels vacuous                         | Section already substantially expanded with PATH framing in the current draft                                               |
| Section 7.2 redundant with introduction                                             | Partial redundancy acknowledged; implication paragraph carries unique content                                               |
| Section 8 Conclusions: rename A and B by underlying model                           | Terminology is load-bearing across documents 02, 19, and 21; renaming is out of scope for this revision                     |

## Output

- Source: `02-dgp-mean-moderation-vs-mvn.tex` (rendered timestamp updated to 2026-04-15 at 10:48 PDT).
- PDF: 19 pages (was 18 pages before the Pass 2 additions).

---

*Rendered on 2026-04-15 at 10:55 PDT.*
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/02-revision-response-to-reviewer.md*
