# The prognostic biomarker and the shared biomarker effect

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** Shared slope; leak factor; expectancy terms; pathways
table; power cost of prognostic biomarker.  
**Prerequisites:** ch17, ch18

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

*Local symbols used in this chapter only. For project-wide notation see
`analysis/report/NOTATION.md`.*

<!-- notation entries -->

---

## Summary

<!-- evidence-labeled results go here -->

---

## Construct box

| Process | Carryover | Coupling | Allocation | Design | N | Analysis ID | $t_{1/2}^a$ (wk) | $n_{	ext{reps}}$ |
|---|---|---|---|---|---|---|---|---|
| not stated (`docs/48`) | not stated | 0 (null); 0.3 (predictive) | not stated | Hybrid | 70 | S, F1 | not stated | 500 per cell |

Table: Construct box: simulations underlying the null-rejection results in
this chapter. All figures are from `docs/48` Section 3 Table 3 and `docs/51`
Feature 6 (inspected; not reproduced from source data in this series).
Process and carryover specification are not recorded in `docs/48`.

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 48 Sections 1 to 3; 49 Sections 2 to 4; 56 Section 3.3;
  45 Section 6; 41 Section 3; 51 Feature 6; CM "Baseline, a
  prognostic biomarker, and the shared biomarker effect" (6
  subsections), "Prognostic pathways".

- New to write: F1 stated as approximately unbiased (CR-21); leak
  factors under the refined process if job 15 supports them.
-->

### The predictive--prognostic distinction

<!-- From docs/48 Sections 1 to 2 -->

A **predictive** biomarker modulates the drug response: participants
with a high biomarker level have a larger (or smaller) benefit from
treatment. A **prognostic** biomarker predicts the natural course of
the outcome regardless of treatment: its slope is the same on drug
and off drug. In the DGP, $c_{bm}$ encodes the predictive pathway
(biomarker-drug-response coupling). Prognostic pathways are switched
off by default (`docs/48`, Section 1; inspected).

The concern is that a prognostic biomarker breaks the
**shared biomarker effect** assumption: the premise that $\beta_B$ is
the same at on-drug visits as at off-drug visits. When the biomarker
also influences natural course, placebo response, or expectancy, the
shared-slope constraint in the standard model conflates the predictive
and prognostic signals.

### Biomarker pathways

<!-- From docs/49 Sections 2 to 3 -->

Five distinct mechanisms connect the biomarker to the observed
symptom trajectory:

| Pathway | Symbol | Description |
|---|---|---|
| Predictive: drug response | $c_{bm}$ | Coupling of $B_i$ to $BR_q(t)$ |
| Prognostic: natural course constant | $c_{TV}$ | $B_i$ shifts $TV$ uniformly |
| Prognostic: natural course growing | $c_{TV,t}$ | $B_i$ slope in $TV$ grows in $t$ |
| Prognostic: placebo response | $c_{PB}$ | $B_i$ shifts $PB(t)$ |
| Prognostic: baseline severity | $c_{BL}$ | $B_i$ shifts the patient mean |

The standard analysis fits a single slope $\beta_B$ across all
visits. Under pathways 1 to 5 active simultaneously, the estimate of
$\beta_{bm:D}$ is

$$
\hat\beta_{bm:D} \approx c_{bm}
  - \underbrace{f(c_{TV}, c_{TV,t}, c_{PB}, c_{BL})}_{\text{leak}},
$$

where the leak factor $f$ is the contribution of the prognostic
pathways to the on--off interaction contrast through the shared-slope
constraint.

### What each analysis assumes at each visit type

<!-- From docs/49 Section 3 -->

| Visit type | S model | ES model | F1 model |
|---|---|---|---|
| Off drug, pre-exposure | $\beta_B$ | $\beta_B$ | $\beta_B + \beta_{bm:De}\,\eta$ |
| Off drug, post-washout | $\beta_B$ | $\beta_B$ | $\beta_B + \beta_{bm:De}\,\eta$ |
| On drug, blinded | $\beta_B$ | $\beta_B + \beta_{bm:D}$ | $\beta_B + \beta_{bm:D} + \beta_{bm:De}\,\eta$ |
| On drug, open label | $\beta_B + \beta_{bm:D}$ | $\beta_B + \beta_{bm:D}$ | $\beta_B + \beta_{bm:D} + \beta_{bm:De}$ |

The S model (ES model without `De`) uses $\beta_B$ at all off-drug
visits and $\beta_B + \beta_{bm:D}$ at all on-drug visits. The ES
model is identical with the $\beta_{bm:De}$ column set to zero. F1
allows the on-drug--open-label and blinded-on-drug cells to differ.

The distinction matters when the biomarker predicts placebo response:
$\eta_t > 0$ at some off-drug blinded visits and $\eta_t < 1$ at some
open-label visits, so the S model's shared slope absorbs part of the
biomarker--expectancy interaction into $\hat\beta_{bm:D}$.

### Null rejection under five scenarios

<!-- From docs/48 Section 3 -->

Type I error (null $c_{bm} = 0$) at $N = 70$, five prognostic
scenarios, S vs F1 analysis:

| Scenario | S null rejection | F1 null rejection |
|---|---|---|
| Baseline only ($c_{BL} = 0.3$) | 0.049 | 0.049 |
| TV constant ($c_{TV} = 0.3$) | 0.099 | 0.049 |
| TV growing ($c_{TV,t} = 0.3$) | 0.049 | 0.049 |
| Placebo response ($c_{PB} = 0.3$) | 0.100 | 0.049 |
| Predictive A ($c_{bm} = 0.3$, power) | 0.416 | 0.195 |
| Predictive A + TV growing | 0.195 | 0.195 |
| Predictive A + placebo response | 0.222 | 0.222 |

(Inspected from `docs/48` Section 3 Table 3; 500 replicates per cell.)

The S analysis is anticonservative whenever the biomarker predicts
a component that covaries with drug state (TV constant, placebo
response). F1 holds 0.049 in both cases. In the predictive-only rows,
F1 costs power: 0.416 (S) against 0.195 (F1) at strength 0.3 (with
no prognostic component). The cost reflects the expectancy terms
absorbing part of the drug contrast. [inspected, docs/48]

### Why F1 is approximately unbiased

<!-- From docs/56 Section 3.3 -->

With five expectancy levels $\eta \in \{0, 0.25, 0.5, 0.75, 1\}$
and the allocation balanced across paths, the within-participant
correlation between the biomarker interaction column $B_i D_{it}$
and the expectancy interaction column $B_i \eta_t$ is approximately

$$
\rho(B_i D_{it},\, B_i \eta_t) \approx
  \frac{\sum_t D_{it}\,\eta_t}{\sqrt{\sum_t D_{it}^2}\,\sqrt{\sum_t\eta_t^2}}
  \approx 0.59
$$

(Hybrid long data; verified from `docs/56` Section 3.3). Despite this
correlation, the partial regression coefficient $\hat\beta_{bm:D}$ in
F1 is approximately unbiased for $c_{bm}$ when the prognostic pathways
are active, because the expectancy interaction column accounts for the
part of the TV and PB trajectories that covaries with both $D_{it}$
and $\eta_t$. [inspected, docs/56]

<!-- TODO: F1 approximately unbiased -- formal argument; reference
job 15 leak-factor runs when available (CR-21) -->

<!-- TODO: leak factors under the refined process -- add table once
job 15 results are available; expected finding is that the graded
coupling reduces the leak through prognostic pathways relative to the
step rule -->

### Freeing the baseline biomarker slope

<!-- From docs/45 Section 6 -->

The cLDA equivalence in `docs/41` Section 7 was applied with a shared
$\beta_B$: one slope for baseline and all follow-up visits. This is
the E9b coding used in the package. Under AR(1) data with a purely
predictive biomarker, the baseline row adds a clean off-drug
observation and the shared $\beta_B$ correctly estimates the
prognostic slope (zero in the standard construct). This gives the
analysis a power advantage of about 0.02 in the Hybrid at $N = 70$.
[inspected, docs/45]

When the biomarker is prognostic, the shared-slope constraint breaks
down because the baseline visit's covariance with $B_i$ is
$\sigma_{B,\text{BL}} \ne \sigma_{B,\text{FU}}$, and the single
slope $\beta_B$ no longer consistently estimates either. Freeing
$\beta_B$ at baseline (E9b with a separate baseline coefficient) is
the correct remedy but costs the power advantage (`docs/45` Section 6;
verified from the stress-test results).

### Feature 6: prognostic pathways in the ten-feature DGP

<!-- From docs/51 Feature 6 -->

The ten-feature DGP (`docs/51`) includes prognostic pathways as
Feature 6. The null-rejection sweep at coupling 0.3, 500 replicates,
Hybrid $N = 70$, showed:

- S analysis: 0.099 (TV constant) and 0.100 (placebo response);
  within 1 MCSE of `docs/48` Section 3.
- F1 analysis: 0.040 to 0.047 across all five prognostic scenarios
  tested (verified from `docs/51` Feature 6 table; 2 MCSE band
  matches).

### The DGP baseline structure

<!-- From docs/41 Section 3 -->

In the construct, the baseline row is drawn from the same
multivariate normal as all follow-up rows, with the symptom
covariance $\sigma_{t0} = \sigma_{00}$ for every $t$ (inspected;
verified numerically in `docs/41` Section 3). This means the
biomarker-baseline covariance equals the biomarker-follow-up
covariance, and the shared-slope constraint is consistent with the
DGP when the biomarker is purely predictive.

Under a prognostic biomarker that operates through baseline severity
($c_{BL} \ne 0$), this no longer holds: the baseline visit absorbs
the severity effect, and the baseline covariance is augmented by
$c_{BL}^2 \sigma_B^2$. The ANCOVA estimator and the shared-slope
E9b estimator both shift in opposite directions. [inspected, docs/41]

---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers
agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at
docs/series/registry/open-items.csv -->

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
| 2026-10-15 | Setting and argument | Initial content move from docs/48, 49, 56, 45, 41, 51 | Phase 2 content moves |
