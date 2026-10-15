# The stress test: protocol, results and the refined-process rerun

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** 14 analyses; 41 cells; stored run; job 15 refined-process rerun; jobs 21,22.  
**Prerequisites:** ch19, ch20

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

| Process | Carryover | Coupling | Allocation | Design | N | Analysis ID | $t_{1/2}^a$ (wk) | $n_{\text{reps}}$ |
|---|---|---|---|---|---|---|---|---|
| config C at $\rho = 0.7$ | compounding (recursive) | graded | fixed | Hybrid | 70$^*$ | S, DS, C1, C2, F1, F1A, G1, G2, F2, ES, EA, A, DA, R | 0.5 | 1000 (null) / 500 (alt); 41 cells |
| refined | anchored | graded | randomized | Hybrid | 70$^*$ | S, DS, C1, C2, F1, F1A, G1, G2, F2, ES, EA, A, DA, R | 0.5 | 1000 (null) / 500 (alt); 41 cells |
| refined | anchored | graded | randomized | Hybrid | 70$^*$ | S, C1, C2, ES, F1, F2, G1, G2 | 0.5 | 1500; 41 cells |

Table: Construct box: distinct DGP configurations producing numeric results in
this chapter. Row 1 is the stored run (Sections 1-11, 2026-10-07). Row 2 is the
refined-process rerun: job 15 (11 analyses) plus jobs 21 and 22 (G1, G2 and F1A
fitted to the same job 15 trials). Row 3 is job 37 (1,500 fresh trials per
cell, post hoc). $^*$Stress cells also vary $N \in \{35, 150\}$.

*Note: stored runs used $t_{1/2}^a = 0.5$ week; the 0.2 week policy applies to
new runs.*

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 47 (spine); 46 (annex, with Amendment 2 item 4 cut to a pointer);
  48 Section 4; 49 Section 5; 50 Section 2; CM "The stress-test construct and
  its fourteen analyses" (5 subsections), "Size, power screens and the
  selection criteria".

- New to write: the job 15 results section; the job 21 and 22 results; the
  comparison of the two runs as independent proportions (`RERUN-anchored.md`
  Section 7).
-->

<!-- Spine: docs/47-stress-test-results.md (copied verbatim; Phase 2) -->

# Stress Test of Candidate Analyses for the Biomarker-by-Drug Interaction in the Hybrid Design: Exploratory Results {.unlisted .unnumbered}
*2026-10-08 09:49 PDT*

**Author.** pmsimstats team

**Purpose.** This document reports the results of the stress test set
out in the working protocol `docs/46`. Fourteen analyses of the
biomarker-by-drug interaction $\beta_{bm:D}$ in the Hybrid design were
fitted to the same 28,000 simulated trials in 41 cells, and judged for
size, robustness to a prognostic biomarker, sensitivity and power. The
study is exploratory: its criteria and its set of analyses were revised
as results came in (`docs/46`, Amendments 1 and 2), so the selected
analysis is a recommendation to be confirmed in an independent study.
The reasoning behind several revisions is in `docs/48`; the
recommendation for practice is in `docs/49`.

**Rerun note (2026-10-11).** Sections 1 to 11 are the stored run
(recursive carryover, fixed allocation; 2026-10-07), and are left as
written. The main run (eleven interaction analyses) has been repeated
under the refined process, with anchored carryover and randomized path
allocation; see the section "Rerun under the refined process
(2026-10-11)" before the references. In brief: size and the prognostic
flags are unchanged within Monte Carlo error; DS and EA now narrowly
fail the size screen; among the eleven analyses F1 remains the
unflagged survivor with the smallest maximin shortfall (0.240, stored
0.242); and at the AR(1) reference cell F1 has 11% less size-adjusted
power than S (0.448 against 0.504), not the 29% of the stored run.

**Decision note (2026-10-15).** G1, G2 and F1A have since been rerun on
the trials of the main rerun (jobs 21 and 22), and 1,500 fresh trials
per cell were added for the eight surviving analyses (job 37, post
hoc). F1 fails nothing; F1A fails the size screen. On the 28,000 rerun
trials G2 led F1 by 0.026 in maximin shortfall, beyond the tie margin
of 0.02; on the fresh trials alone and on all trials pooled, F1, G1 and
G2 lie within the margin (pooled 0.216, 0.221 and 0.221), and the
tie-break selects **F1**. F1 is retained as the recommended analysis
under the refined process; G1 and G2 are post hoc alternatives that tie
with it, and the earlier G2 lead is attributed to Monte Carlo error.
See the subsections "Fourteen-analysis decision" and "Fresh trials for
the F1, G1, G2 choice" of the rerun section.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here or in `docs/46`.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $Y_{it}$ | symptom score of participant $i$ at visit $t$ (code `Sx`); $t = 0$ is baseline | canonical |
| $B_i$, $b_i$ | biomarker; standardized biomarker | canonical |
| $B^c_i$ | biomarker centered at the trial's baseline sample mean, $B_i - \bar B$, not scaled (code `bmc`, the form the analyses use; see the note below the table) | canonical |
| $D_{bc,it}$ | exposure-decayed drug indicator, 0 before any exposure (code `Dbc`) | canonical |
| $\eta_t$ | expectancy: 1 open label, 0.5 blinded, 0 at baseline (code `De`) | canonical |
| $t_{1/2}$, $t_{1/2}^{a}$ | true carryover half-life; half-life assumed by the analysis | canonical; local |
| $c_{bm}$ | biomarker-BR coupling, the predictive effect | canonical |
| $c_{bm,TV}$, $c_{bm,PB}$, $c_{bm,BL}$ | biomarker correlations with the natural course, the placebo response and baseline severity (prognostic couplings) | canonical |
| $\beta_{bm}^{TV}$ | mean moderation of the natural course by $b_i$ (growing prognostic pathway) | canonical |
| $\beta_{bm:D}$ | biomarker-by-drug interaction, the estimand | canonical |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $N$ | total sample size across paths | canonical |
| $\alpha$ | nominal level, 0.05 | canonical |
| $\pi$, $\hat\pi$ | power; an estimated rejection rate (size or power) | canonical |
| $n_{\text{reps}}$ | Monte Carlo replicates per cell | canonical |

Table: Notation index and glossary of symbols used in this document

(Corrected 2026-10-08: `bmc` was described here as the standardized
biomarker. In `11-strawman-stress-test.R` it is $B_i - \bar B$, centered
at the trial's baseline sample mean and not divided by an SD. Test
statistics and $p$-values are invariant to rescaling the biomarker, so
no reported rate changes; the fitted $\beta_{bm:D}$ is per unit of the
raw biomarker. The terms `bmc:De` are therefore written
$B^c_i\eta_t = (B_i - \bar B)\eta_t$ below (renamed from $b_i\eta_t$ in
the notation pass of 2026-10-09). See
`analysis/report/compendium-mathematics.tex`, Section "The stress-test
construct and its fourteen analyses".)

The fourteen analyses (`docs/46`, Section 4 and Amendment 2):

| ID | Analysis |
|----|------------------------------------------------------------|
| S | random intercept plus continuous-time AR(1) residuals (`corCAR1`), CR2 standard errors with Satterthwaite degrees of freedom, baseline as a response row, linear time, $D_{bc,it}$ at $t_{1/2}^{a} = 0.5$ week |
| DS | S with $t_{1/2}^{a}$ chosen by AIC over {0, 0.25, 0.5, 1, 2} |
| C1, C2 | S with residual variances by phase (baseline, open label, blinded); CR2 (C1) or model-based (C2) inference |
| F1 | S plus the expectancy terms $\eta_t$ and $B^c_i \eta_t$ (`De + bmc:De`) |
| F1A | F1 with $t_{1/2}^{a}$ chosen by AIC |
| G1, G2 | F1 with residual variances by phase; CR2 (G1) or model-based (G2) inference |
| F2 | S plus `bmc:t` |
| ES, EA | S and A with the biomarker's baseline association freed (`+ bmc:bl0`) |
| A | unstructured covariance with Kenward-Roger inference, otherwise as S |
| DA | A with $t_{1/2}^{a}$ chosen by AIC |
| R | the published reference arm: compound-symmetry covariance (the `cs` structure of `mmrm`, by visit within participant, equivalent to a random intercept with independent errors), model-based standard errors with Satterthwaite degrees of freedom, binary drug indicator (wording corrected 2026-10-08; `docs/46` described it as a random intercept with model-based inference) |

Table: The fourteen analyses compared in the stress test

**Glossary.**

- **Null and alternative cells.** A null cell has $c_{bm} = 0$ (1,000
  replicates); its alternative twin has $c_{bm} = 0.30$, or 0.25 under
  compound symmetry (500 replicates).
- **Prognostic cell.** A cell in which the biomarker is associated with
  the outcome regardless of treatment: through the natural course
  (constant, $c_{bm,TV} = 0.3$, or growing, $\beta_{bm}^{TV} = 0.3$),
  the placebo response ($c_{bm,PB} = 0.3$) or baseline severity
  ($c_{bm,BL} = 0.3$).
- **Size-adjusted power.** Power when each analysis rejects at the 5%
  quantile of its own p-values in the matching null cell.
- **Shortfall.** In an alternative cell, the gap between an analysis's
  size-adjusted power and that of the best nominal analysis; the
  *maximin shortfall* of an analysis is its largest shortfall over the
  cells considered. (Clarified 2026-10-08: the term is kept because the
  decision outputs use it; selecting the analysis with the smallest
  largest shortfall is a minimax-regret rule. The best analysis in each
  cell is taken over all survivors of the size and sensitivity screens,
  including the flagged C1 and C2, so shortfalls are measured against
  the phase-variance models. See
  `analysis/report/compendium-mathematics.tex`, Section "Size, power
  screens and the selection criteria of the stress test".)
- **MCSE.** Monte Carlo standard error, $\sqrt{\hat\pi(1-\hat\pi)/n_{\text{reps}}}$ for a
  rejection rate $\hat\pi$ over $n_{\text{reps}}$ replicates.

## 1. Summary

All results are from stored simulation outputs (verified; Section 11).

1. **Size.** Ten of the fourteen analyses hold their size over the 11
   non-prognostic null cells (pooled 0.047 to 0.053, within the 2-MCSE
   band 0.046 to 0.054, and no cell above 0.075)
   [inspected, summary-full.csv]. Four fail: the
   unstructured model A (0.057 pooled, 0.087 at $N = 35$), its AIC
   version DA (0.065), the published reference arm R (0.069), and F1A
   (0.057, and 0.073 when $t_{1/2} = 1$ week)
   [inspected, summary-full.csv].
2. **A prognostic biomarker is the decisive stress.** Every analysis
   that shares the biomarker's main effect between the baseline row and
   the later visits rejects too often in at least one prognostic null
   cell: S 0.100, C1 0.382, A 0.416 [inspected, summary-full.csv].
   Five analyses stay within the cap
   in all four prognostic null cells: ES, EA, F1, G1 and G2
   [inspected, summary-full.csv].
3. **Sensitivity.** Every analysis has monotone size-adjusted power
   curves in both constructs [inspected, summary-full.csv]. The
   symmetry screen was dropped, because power is symmetric in the sign
   of $c_{bm}$ by construction.
4. **Power.** Robustness costs power. At $c_{bm} = 0.30$ under AR(1),
   size-adjusted power is 0.514 for S and 0.364 for F1, 29% less
   [inspected, summary-full.csv]; under compound symmetry at
   $c_{bm} = 0.25$ the difference is small (0.554 against 0.522)
   [inspected, summary-full.csv]. F1 has significantly lower power
   than DS in 15 of the 22 non-prognostic alternative cells (McNemar
   test at 0.05) and higher power in none [inspected,
   decision-power-ds-nosym-all14.csv].
5. **Selection.** Among the five robust analyses the maximin shortfall
   over the non-prognostic alternative cells is G2 0.232, F1 0.242, G1
   0.252, ES 0.278 and EA 0.318 [inspected,
   decision-maximin-ds-nosym-all14.csv]. G2 and F1 lie within the tie
   margin of 0.02, neither reaches power 0.80 on the grid
   [inspected, summary-full.csv], and the rule's last tie-break (fewer
   covariance parameters) selects **F1**. G1 and G2 are effectively
   equivalent alternatives.
6. **A second reading.** Over all 26 alternative cells, prognostic cells
   included, G1 (0.308) and G2 (0.312) lead F1 (0.358) by about 0.05
   [inspected, decision-maximin-ds-nosym-all14.csv]. The protocol does
   not say which set of cells governs the selection, so both readings
   are reported (Section 7).

## 2. Design in brief

The design is set out in `docs/46`; this section records what was run.

- **Design.** Hybrid, visits at weeks 4, 8, 9, 10, 11, 12, 16 and 20
  plus baseline, four paths with the published allocation.
- **Data-generating process.** Configuration C of `docs/36`: separable
  AR(1) correlation at $\rho = 0.7$, graded coupling of the biomarker
  to the drug response, the published recursive carryover on the
  drug-response mean, the published Gompertz means and SDs, the fixed
  published allocation, and a stop on any matrix that is not positive
  definite; base cell $N = 70$, $t_{1/2} = 0.5$ week, no dropout,
  exponential decay. (Corrected 2026-10-08: this bullet said "anchored
  carryover". `11-strawman-stress-test.R` decays the previous visit's
  mean with the cumulative time since discontinuation, the published
  recursion (inspected). The construct is therefore the covariance
  construct of the refined process of paper 01, not the refined
  process in full.)
- **Cells.** 41: two power curves (AR(1), $c_{bm}$ from $-0.30$ to
  0.45; compound symmetry, $-0.20$ to 0.25) and stress cells varying
  one factor each (measurement error, $N = 35$ and 150, $t_{1/2} = 0$
  and 1 week, Weibull decay with shape 0.5 and 2, 20% MCAR and MAR
  dropout, and the four prognostic pathways). There are 15 null and 26
  alternative cells and 28,000 simulated trials.
- **Analyses.** All fourteen are fitted to the same trials, so every
  comparison is paired.
- **How the protocol changed** (`docs/46`, Amendments 1 and 2):
  1. the reference analysis became DS, and the symmetry screen was
     dropped;
  2. the two baseline-severity cells (40 and 41) were added and judged
     with the other prognostic cells;
  3. G1, G2 and F1A were added and judged by the same rule;
  4. failed fits are excluded from rejection rates, rather than counted
     as non-rejections as first written; the two conventions differ by
     at most 0.0025 in any rate.

## 3. Size

| Analysis | Pooled | Largest cell | Nominal | TV, constant | TV, growing | Placebo | Baseline | Prognostic flag |
|------|------|--------------------|------|------|------|------|------|------|
| S | 0.048 | 0.059 (MCAR dropout) | yes | 0.099 | 0.049 | 0.100 | 0.050 | yes |
| DS | 0.053 | 0.068 (MCAR dropout) | yes | 0.115 | 0.058 | 0.112 | 0.055 | yes |
| C1 | 0.053 | 0.069 ($N = 35$) | yes | 0.382 | 0.173 | 0.205 | 0.054 | yes |
| C2 | 0.052 | 0.066 ($N = 35$) | yes | 0.393 | 0.172 | 0.193 | 0.050 | yes |
| F1 | 0.051 | 0.062 (AR(1) curve) | yes | 0.047 | 0.046 | 0.040 | 0.040 | no |
| F1A | 0.057 | 0.073 ($t_{1/2} = 1$) | no | 0.059 | 0.052 | 0.047 | 0.045 | no |
| G1 | 0.049 | 0.060 (AR(1) curve) | yes | 0.055 | 0.044 | 0.042 | 0.042 | no |
| G2 | 0.047 | 0.066 (AR(1) curve) | yes | 0.062 | 0.043 | 0.045 | 0.046 | no |
| F2 | 0.049 | 0.059 (MCAR dropout) | yes | 0.090 | 0.048 | 0.098 | 0.048 | yes |
| ES | 0.048 | 0.059 (AR(1) curve) | yes | 0.049 | 0.063 | 0.049 | 0.050 | no |
| EA | 0.052 | 0.070 ($N = 35$) | yes | 0.049 | 0.057 | 0.050 | 0.066 | no |
| A | 0.057 | 0.087 ($N = 35$) | no | 0.416 | 0.195 | 0.222 | 0.051 | yes |
| DA | 0.065 | 0.103 ($N = 35$) | no | 0.509 | 0.245 | 0.256 | 0.064 | yes |
| R | 0.069 | 0.081 ($N = 35$) | no | 0.107 | 0.053 | 0.159 | 0.075 | yes |

Table: Null rejection rates by analysis: pooled over the 11 non-prognostic
null cells, largest single cell, and the four prognostic null cells.
*Source: summary-full.csv [inspected]*

The pooled rate is over 11,000 replicates, so its 2-MCSE band is 0.046
to 0.054; a single cell has 1,000 replicates (MCSE 0.0069), and the
criterion caps each cell at 0.075. An analysis is flagged when any
prognostic null cell exceeds 0.075.

**Reading.**

- **Inference matters most at small $N$.** The unstructured model with
  Kenward-Roger inference [5] rejects 0.087 at $N = 35$, and with an
  AIC-chosen half-life 0.103 [inspected, summary-full.csv]. The
  random-intercept models with CR2 standard errors [3, 4] stay at or
  below 0.069 in every cell, except F1A (0.073), whose half-life is
  chosen by AIC [inspected, summary-full.csv].
- **Choosing the half-life by AIC inflates size slightly.** DS (0.053)
  passes; F1A (0.057) and DA (0.065) do not [inspected,
  summary-full.csv]. The selection step adds variability that the
  standard errors do not reflect.
- **The published reference arm is liberal** throughout (0.069 pooled)
  [inspected, summary-full.csv], as earlier documents found (`docs/45`).
- **Baseline severity alone does little harm.** A prognostic association
  with baseline severity enters every visit equally, so every analysis
  stays within the cap in that cell (0.040 to 0.075; R reaches the cap
  exactly, and EA reaches 0.066) [inspected, summary-full.csv].

## 4. Robustness to a prognostic biomarker

Analyses that keep baseline as a response row and use one biomarker
main effect for it and for the later visits attribute part of a
prognostic association to the interaction (derived, `docs/48`, Section
2). The size table shows this directly: S, DS, C1, C2, F2, A, DA and R
are flagged, with null rejection rates up to 0.509 when the natural
course is associated with the biomarker [inspected,
summary-full.csv]. F2's `bmc:t` term remedies the growing
natural-course pathway (0.048) but not the constant one (0.090) or the
placebo pathway (0.098) [inspected, summary-full.csv].

Two families stay within the cap in every prognostic cell:

- **the expectancy-term models** F1, G1 and G2, which add $\eta_t$ and
  $B^c_i \eta_t$ and so let the biomarker's association differ between
  baseline, open-label and blinded visits;
- **the free-baseline models** ES and EA, which give the baseline row
  its own biomarker slope.

Staying within the cap is an empirical result, not exact invariance. In
the Hybrid design F1 leaves a small residual bias under a constant
natural-course pathway, because the blinded on-drug and off-drug visits
carry unequal weight, and ES leaks a little under the placebo pathway
through the open-label visits (derived; `docs/49`, Section 3, and
`docs/48`, Section 2).

The protocol's rule recommends a flagged analysis only together with a
remedy (ES, EA, F1 or F2) that is nominal in the failing cell. The
remedies that pass every prognostic cell are themselves candidates, so
the selection is made among the unflagged analyses.

## 5. Sensitivity

| Analysis | Monotone (AR(1) / CS) | Power at $c_{bm} = -0.30$ / $+0.30$ | Steepness | Detectable coupling |
|---|---|---|---|---|
| S | yes / yes | 0.450 / 0.514 | 2.32 | 0.446 |
| DS | yes / yes | 0.422 / 0.486 | 2.21 | not reached |
| C1 | yes / yes | 0.524 / 0.578 | 2.49 | 0.403 |
| C2 | yes / yes | 0.554 / 0.606 | 2.65 | 0.403 |
| F1 | yes / yes | 0.352 / 0.364 | 1.80 | not reached |
| F1A | yes / yes | 0.326 / 0.338 | 1.67 | not reached |
| G1 | yes / yes | 0.340 / 0.354 | 1.72 | not reached |
| G2 | yes / yes | 0.360 / 0.380 | 1.75 | not reached |
| F2 | yes / yes | 0.446 / 0.516 | 2.36 | not reached |
| ES | yes / yes | 0.336 / 0.364 | 1.71 | not reached |
| EA | yes / yes | 0.300 / 0.306 | 1.48 | not reached |
| A | yes / yes | 0.510 / 0.534 | 2.48 | 0.430 |
| DA | yes / yes | 0.472 / 0.498 | 2.25 | 0.434 |
| R | yes / yes | 0.394 / 0.398 | 1.76 | not reached |

Table: Size-adjusted power curves under AR(1): monotonicity in both
constructs, power at plus and minus 0.30, steepness between 0.15 and
0.30, and the smallest coupling reaching power 0.80.
*Source: summary-full.csv [inspected]*

**Reading.**

- **Every curve is monotone**, so the sensitivity screen removes no
  analysis [inspected, summary-full.csv].
- **The symmetry screen was defective.** The data-generating process
  treats $+c_{bm}$ and $-c_{bm}$ alike, so power is symmetric in theory;
  the observed differences between $-0.30$ and $+0.30$ (at most 0.070
  for F2 [inspected, summary-full.csv], against an MCSE of about 0.03
  for a difference of two independent cells of 500 [inferred,
  Monte Carlo theory]) are consistent with Monte Carlo error, and the
  screen would have eliminated analyses for noise. It was dropped in
  Amendment 1 (derived, `docs/48`, Section 4).
- **The robust analyses are flatter.** Their steepness is 1.5 to 1.8,
  against 2.2 to 2.7 for the flagged analyses other than R (1.76)
  [inspected, summary-full.csv], and none reaches power 0.80 by
  $c_{bm} = 0.45$ at $N = 70$ [inspected, summary-full.csv]. Steepness
  and detectable coupling therefore cannot break ties among them.

## 6. Power

| Analysis | AR(1), $c_{bm} = 0.30$ | Compound symmetry, $c_{bm} = 0.25$ |
|---|---|---|
| S | 0.514 (0.508) | 0.554 (0.554) |
| DS | 0.486 (0.516) | 0.480 (0.476) |
| C1 | 0.578 (0.602) | 0.660 (0.658) |
| C2 | 0.606 (0.606) | 0.656 (0.630) |
| F1 | 0.364 (0.410) | 0.522 (0.544) |
| F1A | 0.338 (0.436) | 0.446 (0.508) |
| G1 | 0.354 (0.408) | 0.546 (0.552) |
| G2 | 0.380 (0.426) | 0.528 (0.504) |
| F2 | 0.516 (0.496) | 0.552 (0.552) |
| ES | 0.364 (0.412) | 0.662 (0.634) |
| EA | 0.306 (0.392) | 0.602 (0.596) |
| A | 0.534 (0.566) | 0.618 (0.618) |
| DA | 0.498 (0.592) | 0.572 (0.596) |
| R | 0.398 (0.476) | 0.596 (0.594) |

Table: Size-adjusted power (unadjusted in parentheses) in the two
reference alternative cells, $N = 70$.
*Source: summary-full.csv [inspected]*

| Cell | S | C1 | F1 | G1 | G2 | ES | EA |
|--------------------------|-----|-----|-----|-----|-----|-----|-----|
| $N = 150$ | 0.828 | 0.918 | 0.748 | 0.752 | 0.754 | 0.704 | 0.698 |
| $N = 35$ | 0.222 | 0.272 | 0.192 | 0.214 | 0.218 | 0.206 | 0.196 |
| $t_{1/2} = 0$ | 0.594 | 0.666 | 0.496 | 0.510 | 0.532 | 0.504 | 0.490 |
| $t_{1/2} = 1$ week | 0.388 | 0.502 | 0.272 | 0.282 | 0.306 | 0.288 | 0.294 |
| Weibull decay, shape 0.5 | 0.404 | 0.514 | 0.308 | 0.296 | 0.378 | 0.270 | 0.230 |
| Weibull decay, shape 2 | 0.588 | 0.682 | 0.504 | 0.526 | 0.526 | 0.458 | 0.464 |
| Measurement error | 0.400 | 0.444 | 0.318 | 0.320 | 0.330 | 0.336 | 0.324 |
| 20% MCAR dropout | 0.426 | 0.542 | 0.378 | 0.394 | 0.406 | 0.360 | 0.346 |
| 20% MAR dropout | 0.418 | 0.518 | 0.382 | 0.404 | 0.378 | 0.356 | 0.322 |
| Prognostic, TV constant | 0.606 | 0.774 | 0.482 | 0.520 | 0.524 | 0.370 | 0.378 |
| Prognostic, TV growing | 0.516 | 0.672 | 0.378 | 0.484 | 0.490 | 0.192 | 0.208 |
| Prognostic, placebo | 0.652 | 0.736 | 0.446 | 0.428 | 0.424 | 0.498 | 0.438 |
| Prognostic, baseline | 0.444 | 0.552 | 0.396 | 0.394 | 0.406 | 0.336 | 0.290 |

Table: Size-adjusted power at $c_{bm} = 0.30$ in the stress cells: the
reference analysis S, the most powerful analysis C1, and the five
robust analyses. *Source: summary-full.csv [inspected]*

**Reading.**

- **The phase-variance models are the most powerful** (C1 has the
  highest power in every stress cell shown [inspected,
  summary-full.csv]), but they fail the prognostic cells badly, so
  their power is not available to a trial whose biomarker may be
  prognostic.
- **The cost of robustness depends on the correlation structure.** Under
  AR(1) the robust analyses have 26% (G2) to 40% (EA) less
  size-adjusted power than S at $c_{bm} = 0.30$ [inferred,
  summary-full.csv]; under compound symmetry ES even exceeds S (0.662
  against 0.554), and F1 is close (0.522) [inspected,
  summary-full.csv].
- **The power reported for flagged analyses in prognostic cells is
  size-adjusted**, so it is not the power a trial would see: unadjusted,
  those analyses reject far too often under the null.
- **The free-baseline models lose most when the natural course grows
  with the biomarker** (0.192 and 0.208) [inspected, summary-full.csv],
  which is why their all-cells shortfall is large (Section 7).

## 7. Selection

The rule (`docs/46`, Section 7, as revised by Amendment 1) is applied to
the analyses that pass the size screen and the sensitivity screen and
are not flagged in a prognostic cell: ES, EA, F1, G1 and G2.

| Analysis | Largest shortfall, non-prognostic cells | Worst cell | Largest shortfall, all cells | Worst cell |
|---|---|---|---|---|
| C1 (flagged) | 0.034 | Weibull 0.5 | 0.064 | TV growing |
| C2 (flagged) | 0.040 | Weibull 2 | 0.040 | Weibull 2 |
| S (flagged) | 0.144 | Weibull 0.5 | 0.220 | TV growing |
| F2 (flagged) | 0.160 | Weibull 0.5 | 0.244 | TV growing |
| DS (flagged) | 0.182 | CS curve, 0.25 | 0.214 | TV growing |
| G2 | 0.232 | AR(1) curve, 0.375 | 0.312 | placebo |
| F1 | 0.242 | AR(1) curve, 0.30 | 0.358 | TV growing |
| G1 | 0.252 | Weibull 0.5 | 0.308 | placebo |
| ES | 0.278 | Weibull 0.5 | 0.544 | TV growing |
| EA | 0.318 | Weibull 0.5 | 0.528 | TV growing |

Table: Maximin shortfall in size-adjusted power among the analyses that
pass the size and sensitivity screens, over the non-prognostic
alternative cells and over all alternative cells.
*Source: decision-maximin-ds-nosym-all14.csv [inspected]*

**Outcome.**

1. **Smallest shortfall.** Over the non-prognostic alternative cells,
   G2 (0.232) and F1 (0.242) are the two best robust analyses
   [inspected, decision-maximin-ds-nosym-all14.csv].
2. **Tie.** They differ by 0.010 [inferred,
   decision-maximin-ds-nosym-all14.csv], inside the rule's tie margin
   of 0.02.
3. **First tie-break.** The smaller detectable coupling in the base
   construct is undefined: neither reaches power 0.80 on the grid.
4. **Last tie-break.** The analysis with fewer covariance parameters,
   **F1**, is selected. F1 is also the only one of the two with CR2
   inference; G2 relies on model-based standard errors.
5. **Strawman rule.** The reference analysis DS is flagged in the
   prognostic cells, so it cannot be retained; F1 replaces it. The
   paired comparison makes the price explicit: F1 has significantly
   lower size-adjusted power than DS in 15 of the 22 non-prognostic
   alternative cells, and higher power in none [inspected,
   decision-power-ds-nosym-all14.csv].

**The all-cells reading.** If the prognostic alternative cells are
included in the maximin set, the order among the robust analyses is G1
0.308, G2 0.312, F1 0.358, EA 0.528, ES 0.544, and the phase-variance
versions lead F1 by about 0.05 [inspected,
decision-maximin-ds-nosym-all14.csv]. The working protocol does not say
whether those cells count. The non-prognostic reading is reported as
the outcome because the prognostic cells are judged by their own
screen. Under the all-cells reading G1 and G2 would tie (0.004 apart)
ahead of F1, and the tie-breaks could not separate them, since they
share one model and differ only in inference. The difference from F1
arises mainly when the natural course is associated with the biomarker
(G1 0.484 against F1 0.378 in that cell) [inspected, summary-full.csv].

## 8. Choosing the half-life by AIC

DS, DA and F1A choose $t_{1/2}^{a}$ by AIC from {0, 0.25, 0.5, 1, 2}
weeks (0 meaning a binary drug indicator), with maximum likelihood for
the choice and a REML refit.

| Analysis | 0 | 0.25 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|
| DS | 0.170 | 0.144 | 0.299 | 0.208 | 0.179 |
| DA | 0.188 | 0.144 | 0.294 | 0.181 | 0.193 |
| F1A | 0.236 | 0.163 | 0.292 | 0.190 | 0.120 |

Table: Share of fits choosing each analysis half-life, over all cells.
*Source: replicates-full.rds [inspected]*

The choice is spread across the grid even though the true half-life is
0.5 week in most cells, and the selection step costs size: F1A fails
the size criterion where F1 passes, and DS loses about 0.03 of
size-adjusted power against S at the AR(1) reference cell (0.486
against 0.514) [inspected, summary-full.csv]. A fixed
$t_{1/2}^{a} = 0.5$ week is therefore preferred with F1.

## 9. Failed fits

A fit that returned an error was recorded as missing. Failures number
0 to 12 per analysis over 28,000 trials (G1 and G2 12; S, C1 and C2 6;
F1, F1A and ES 3; DS and F2 2; A, DA, EA and R 0), and at most 0.6% in
any cell (G1 and G2), below the 2% at which an analysis would be
flagged [inspected, replicates-full.rds].

## 10. Limitations

- **Exploratory.** The criteria, the reference analysis and the set of
  analyses were revised after results were seen. The selection is
  therefore a recommendation from an exploratory comparison, and
  should be confirmed in an independent simulation study before it is
  relied on.
- **One design and one effect size.** Only the Hybrid design was
  studied, with the alternative cells at $c_{bm} = 0.30$ (0.25 under
  compound symmetry) and prognostic couplings of 0.3.
- **The set of cells governs the maximin.** The outcome depends on
  whether the prognostic alternative cells count (Section 7).
- **Monte Carlo error.** Alternative cells have 500 replicates (MCSE of
  power at most 0.022), so differences in shortfall below about 0.03
  between two analyses are within noise for a single cell.
- **The prognostic mechanisms exist only in the stress-test driver,**
  not yet in the package (`docs/44`).

## 11. Reproducibility

All numbers in this document were read from stored outputs on
2026-10-08:

- `analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R`
  produced the 28,000 trials and the fits (main run, `--extra` for G1
  and G2, `--extra2` for F1A); replicate and summary files are in
  `analysis/data/quick-sim/carryover-closed-form/strawman-stress/`
  (`replicates-full*.rds`, `summary-full*.csv`).
- `12-strawman-stress-decision.R --strawman DS --no-symmetry
  --with-extras` applies the rule to all fourteen analyses and writes
  `decision-screen-ds-nosym-all14.csv`,
  `decision-sensitivity-ds-nosym-all14.csv`,
  `decision-power-ds-nosym-all14.csv` and
  `decision-maximin-ds-nosym-all14.csv`.
- Seeds are set per chunk as $20261007 + 1000 \times \text{cell} +
  \text{chunk}$, so results do not depend on the number of cores.

## Rerun under the refined process (2026-10-11)

Sections 1 to 11 above report the stored run (recursive carryover,
fixed allocation; 2026-10-07). This section reports the main run of the
stress test repeated under the refined process. It covers the eleven
interaction analyses of the main run (S, DS, C1, C2, F1, F2, ES, EA, A,
DA and R). The three analyses of the `--extra` and `--extra2` runs (G1,
G2 and F1A) were rerun on the same trials afterwards (jobs 21 and 22),
and the fourteen-analysis decision of Section 7 is repeated in the
subsection "Fourteen-analysis decision" below. The subsection "Fresh
trials for the F1, G1, G2 choice" reports a further, post hoc run that
was made because that decision was within Monte Carlo error of a tie.
(Updated 2026-10-15; the first paragraphs of this section were written
on 2026-10-11, before jobs 21, 22 and 37.)

### What changed in the data-generating process

Two features of the stored run departed from the refined process
(`analysis/report/NOTATION.md`, "Named data-generating processes";
Section 2 of this document). Both were changed in
`11-strawman-stress-test.R` (`RERUN-anchored.md`, Section 1):

1. **Anchored carryover.** Off drug, the mean drug response is
   $\mu_{\text{last}} \, 2^{-t_{sd}/t_{1/2}}$, with $\mu_{\text{last}}$
   the mean at the last on-drug visit, weighted in the Weibull cells by
   the cell's decay weight. The published recursion decayed the previous
   visit's already-decayed mean by the cumulative time since
   discontinuation, so it counted elapsed time twice from the second
   off-drug visit of a run onward (at $t_{1/2} = 0.5$ week the
   recursive mean at the second off-drug visit is a quarter of the
   anchored one). The means of the $t_{1/2} = 0$ cell are unchanged.
2. **Randomized path allocation.** In every trial each participant is
   assigned to one of the four Hybrid paths uniformly at random,
   redrawn when a path receives fewer than 2 participants, instead of
   the fixed published allocation.

The 41 cells, the 28,000 trials (1,000 per null cell, 500 per
alternative cell), the eleven analyses and the per-chunk seeds are
those of the stored run. The trials are nevertheless not the same:
the means differ, and the allocation draws consume random numbers. The
stored and rerun rates are therefore independent estimates, and a
difference is judged by $z = (\hat\pi_{\text{new}} -
\hat\pi_{\text{old}}) / \sqrt{\text{MCSE}_{\text{old}}^2 +
\text{MCSE}_{\text{new}}^2}$. Within the rerun, all eleven analyses
are fitted to the same trials, so comparisons between analyses remain
paired.

The run is simulation-server job 15 at commit `e077526`
(`replicates-full-anchored.rds`, `summary-full-anchored.csv`). The
decision rule was applied with
`12-strawman-stress-decision.R --strawman DS --no-symmetry --tag
-anchored` (Amendment 1 reading) and with `--tag -anchored` alone (the
reading as first specified), writing the screen, sensitivity, power
and maximin files `decision-*-ds-nosym-anchored.csv` and
`decision-*-anchored.csv`. The stored-run figures below were recomputed
from `replicates-full.rds` with the same script, for the same eleven
analyses and all 41 cells (verified: the recomputation with
`--with-extras` reproduces the four stored `-all14` files exactly). As
in the stored decision files, a failed fit counts as a non-rejection.
There were 26 failed fits in 308,000 (stored run 28; 0 to 5 per
analysis) [inspected, replicates-full-anchored.rds].

### Size

| Analysis | Pooled, stored | Pooled, rerun | Largest cell, stored | Largest cell, rerun | Nominal, stored / rerun |
|----------|--------|--------|--------------------|--------------------|------------|
| S | 0.0483 | 0.0516 | 0.059 (MCAR dropout) | 0.061 ($N = 35$) | yes / yes |
| DS | 0.0528 | 0.0569 | 0.068 (MCAR dropout) | 0.067 (Weibull 0.5) | yes / **no** |
| C1 | 0.0531 | 0.0529 | 0.069 ($N = 35$) | 0.065 (meas. error) | yes / yes |
| C2 | 0.0519 | 0.0492 | 0.066 ($N = 35$) | 0.055 (meas. error) | yes / yes |
| F1 | 0.0511 | 0.0514 | 0.062 (AR(1) curve) | 0.064 ($t_{1/2} = 0$) | yes / yes |
| F2 | 0.0486 | 0.0519 | 0.059 (MCAR dropout) | 0.060 ($N = 35$) | yes / yes |
| ES | 0.0478 | 0.0490 | 0.059 (AR(1) curve) | 0.060 (Weibull 0.5) | yes / yes |
| EA | 0.0517 | 0.0575 | 0.070 ($N = 35$) | 0.076 ($N = 35$) | yes / **no** |
| A | 0.0573 | 0.0571 | 0.087 ($N = 35$) | 0.084 ($N = 35$) | no / no |
| DA | 0.0653 | 0.0656 | 0.103 ($N = 35$) | 0.093 ($N = 35$) | no / no |
| R | 0.0687 | 0.0686 | 0.081 ($N = 35$) | 0.084 ($N = 35$) | no / no |

Table: Rerun under the refined process: null rejection rates of the
eleven main-run analyses over the 11 non-prognostic null cells (pooled,
2-MCSE band 0.0458 to 0.0542; largest single cell, cap 0.075), stored
run and rerun. *Source: summary-full-anchored.csv,
summary-full.csv [inspected]*

| Analysis | TV, constant | TV, growing | Placebo | Baseline | Prognostic flag |
|----------|-------------|-------------|-------------|-------------|----------|
| S | 0.099 (0.099) | 0.055 (0.049) | 0.114 (0.100) | 0.056 (0.050) | yes |
| DS | 0.127 (0.115) | 0.068 (0.058) | 0.130 (0.112) | 0.057 (0.055) | yes |
| C1 | 0.367 (0.382) | 0.196 (0.173) | 0.191 (0.205) | 0.056 (0.054) | yes |
| C2 | 0.382 (0.393) | 0.190 (0.172) | 0.178 (0.193) | 0.051 (0.050) | yes |
| F1 | 0.050 (0.047) | 0.050 (0.046) | 0.048 (0.040) | 0.052 (0.040) | no |
| F2 | 0.094 (0.090) | 0.052 (0.048) | 0.112 (0.098) | 0.054 (0.048) | yes |
| ES | 0.045 (0.049) | 0.073 (0.063) | 0.063 (0.049) | 0.057 (0.050) | no |
| EA | 0.049 (0.049) | 0.057 (0.057) | 0.063 (0.050) | 0.061 (0.066) | no |
| A | 0.411 (0.416) | 0.210 (0.195) | 0.193 (0.222) | 0.063 (0.051) | yes |
| DA | 0.515 (0.509) | 0.276 (0.245) | 0.229 (0.256) | 0.065 (0.064) | yes |
| R | 0.119 (0.107) | 0.077 (0.053) | 0.167 (0.159) | 0.069 (0.075) | yes |

Table: Rerun under the refined process: null rejection rates in the four
prognostic null cells (rerun, stored run in parentheses; cap 0.075) and
the prognostic flag, which is the same in both runs.
*Source: summary-full-anchored.csv [inspected]*

**Reading.**

- **Size is unchanged within Monte Carlo error.** Over the 165 null
  rates (15 null cells by 11 analyses), 7 differ from the stored run by
  more than 2 in $|z|$, about the 8 expected by chance; the largest is
  EA at $t_{1/2} = 0$ (0.040 stored, 0.069 rerun, $z = 2.9$)
  [inspected, summary-full-anchored.csv].
- **Two analyses cross a boundary.** DS (pooled 0.0569, outside the
  band) and EA (pooled 0.0575, and 0.076 at $N = 35$ against the cap of
  0.075) now fail the size screen [inspected,
  summary-full-anchored.csv]. Neither change is large against Monte
  Carlo error ($z = 1.4$ and 2.0 for the pooled rates [inferred,
  summary-full-anchored.csv]); both were close to the boundary in the
  stored run. A, DA and R fail as before.
- **The prognostic flags are unchanged** for every analysis. F1 stays
  at or below 0.052 in all four prognostic null cells [inspected,
  summary-full-anchored.csv]. ES stays within the cap with a narrow
  margin in the growing natural-course cell (0.073, MCSE 0.008)
  [inspected, summary-full-anchored.csv].

### Survivors

Every size-adjusted power curve is again monotone in both constructs,
so the sensitivity screen (Amendment 1, no symmetry) removes no
analysis.

| Reading | Stored run (11 analyses) | Rerun |
|---|---|---|
| Size and sensitivity screens (Amendment 1) | S, DS, C1, C2, F1, F2, ES, EA | S, C1, C2, F1, F2, ES |
| Of these, not flagged in a prognostic cell | F1, ES, EA | F1, ES |
| As first specified (strawman S, symmetry screen on) | C2, F1 | C1, F1, ES |

Table: Rerun under the refined process: analyses passing the screens,
stored run and rerun.
*Source: decision-screen-ds-nosym-anchored.csv [inspected]*

Under Amendment 1 the strawman DS fails the size screen, so the
strawman rule (Section 7, item 5) is not reached; DS was in any case
flagged in the prognostic cells. Under the reading as first specified,
S fails the symmetry screen on the compound-symmetry curve in the rerun
(it passed in the stored run) [inspected,
decision-screen-ds-nosym-anchored.csv], and the survivors change from
C2 and F1 to C1, F1 and ES. The symmetry screen thus changes its verdict between
two runs of the same design, as expected of a screen that tests
Monte Carlo noise (Section 5; `docs/48`, Section 4).

### Maximin shortfall

| Analysis | Non-prognostic cells, rerun (stored) | Worst cell, rerun | All cells, rerun (stored) | Worst cell, rerun |
|----------|------------------|----------------|------------------|----------------|
| C2 (flagged) | 0.022 (0.040) | AR(1) curve, 0.45 | 0.058 (0.040) | TV constant |
| C1 (flagged) | 0.032 (0.034) | measurement error | 0.032 (0.064) | measurement error |
| F2 (flagged) | 0.164 (0.160) | Weibull 2 | 0.208 (0.244) | TV growing |
| S (flagged) | 0.172 (0.144) | Weibull 2 | 0.172 (0.220) | Weibull 2 |
| DS (flagged) | fails size (0.182) | | fails size (0.214) | |
| F1 | 0.240 (0.242) | Weibull 0.5 | 0.348 (0.358) | placebo |
| ES | 0.266 (0.278) | Weibull 0.5 | 0.498 (0.544) | TV growing |
| EA | fails size (0.318) | | fails size (0.528) | |

Table: Rerun under the refined process: maximin shortfall in
size-adjusted power among the survivors of the size and sensitivity
screens (eleven analyses; the best analysis in each cell is taken over
all survivors, flagged ones included).
*Source: decision-maximin-ds-nosym-anchored.csv [inspected]*

Among the eleven analyses, F1 has the smallest maximin shortfall of
the unflagged survivors in both readings: 0.240 against 0.266 for ES
over the non-prognostic cells, a lead of 0.026, just outside the tie
margin of 0.02, and 0.348 against 0.498 over all cells [inspected,
decision-maximin-ds-nosym-anchored.csv]. F1's own shortfall is almost
the same as in the stored run (0.240 against 0.242) [inspected,
decision-maximin-ds-nosym-all14.csv and
decision-maximin-ds-nosym-anchored.csv]. Whether G1 or G2 again come within the tie margin of F1, as G2
did in the stored run (0.232 against 0.242), is taken up in the
subsections "Fourteen-analysis decision" and "Fresh trials for the F1,
G1, G2 choice" below.

### Power

| Analysis | AR(1), $c_{bm} = 0.30$, rerun | stored | Compound symmetry, $c_{bm} = 0.25$, rerun | stored |
|---|---|---|---|---|
| S | 0.504 (0.506) | 0.514 (0.508) | 0.490 (0.498) | 0.554 (0.554) |
| DS | 0.502 (0.522) | 0.486 (0.516) | 0.372 (0.404) | 0.480 (0.476) |
| C1 | 0.612 (0.594) | 0.578 (0.602) | 0.624 (0.632) | 0.660 (0.658) |
| C2 | 0.606 (0.602) | 0.606 (0.606) | 0.616 (0.604) | 0.656 (0.630) |
| F1 | 0.448 (0.424) | 0.364 (0.410) | 0.510 (0.534) | 0.522 (0.544) |
| F2 | 0.482 (0.490) | 0.516 (0.496) | 0.488 (0.494) | 0.552 (0.552) |
| ES | 0.376 (0.392) | 0.364 (0.412) | 0.602 (0.572) | 0.662 (0.634) |
| EA | 0.390 (0.392) | 0.306 (0.392) | 0.564 (0.548) | 0.602 (0.596) |
| A | 0.556 (0.562) | 0.534 (0.566) | 0.566 (0.556) | 0.618 (0.618) |
| DA | 0.540 (0.568) | 0.498 (0.592) | 0.524 (0.542) | 0.572 (0.596) |
| R | 0.448 (0.482) | 0.398 (0.476) | 0.488 (0.546) | 0.596 (0.594) |

Table: Rerun under the refined process: size-adjusted power (unadjusted
in parentheses) in the two reference alternative cells, $N = 70$,
rerun and stored run.
*Source: summary-full-anchored.csv; summary-full.csv [inspected]*

| Cell | S | C1 | F1 | ES |
|----------------|-------------|-------------|-------------|-------------|
| $N = 150$ | 0.840 (0.828) | 0.906 (0.918) | 0.758 (0.748) | 0.746 (0.704) |
| $N = 35$ | 0.258 (0.222) | 0.360 (0.272) | 0.208 (0.192) | 0.188 (0.206) |
| $t_{1/2} = 0$ | 0.498 (0.594) | 0.556 (0.666) | 0.398 (0.496) | 0.428 (0.504) |
| $t_{1/2} = 1$ week | 0.374 (0.388) | 0.454 (0.502) | 0.274 (0.272) | 0.276 (0.288) |
| Weibull decay, shape 0.5 | 0.392 (0.404) | 0.522 (0.514) | 0.288 (0.308) | 0.262 (0.270) |
| Weibull decay, shape 2 | 0.494 (0.588) | 0.666 (0.682) | 0.498 (0.504) | 0.466 (0.458) |
| Measurement error | 0.396 (0.400) | 0.462 (0.444) | 0.282 (0.318) | 0.286 (0.336) |
| 20% MCAR dropout | 0.446 (0.426) | 0.556 (0.542) | 0.356 (0.378) | 0.352 (0.360) |
| 20% MAR dropout | 0.458 (0.418) | 0.508 (0.518) | 0.344 (0.382) | 0.356 (0.356) |
| Prognostic, TV constant | 0.598 (0.606) | 0.752 (0.774) | 0.456 (0.482) | 0.402 (0.370) |
| Prognostic, TV growing | 0.544 (0.516) | 0.684 (0.672) | 0.430 (0.378) | 0.206 (0.192) |
| Prognostic, placebo | 0.586 (0.652) | 0.746 (0.736) | 0.398 (0.446) | 0.448 (0.498) |
| Prognostic, baseline | 0.490 (0.444) | 0.582 (0.552) | 0.416 (0.396) | 0.406 (0.336) |

Table: Rerun under the refined process: size-adjusted power at
$c_{bm} = 0.30$ in the stress cells for S, C1 and the two unflagged
survivors F1 and ES (stored run in parentheses).
*Source: summary-full-anchored.csv [inspected]*

**F1 against S.** At the AR(1) reference cell F1's size-adjusted power
is 0.448 against 0.504 for S, 11% less, where the stored run gave
0.364 against 0.514, 29% less (Section 1, item 4)
[inspected, summary-full-anchored.csv]. F1 rose by 0.084 ($z = 2.7$)
while S changed by $-0.010$ [inferred, summary-full-anchored.csv].
Over the 20 non-prognostic alternative cells with $c_{bm} > 0$, the
ratio of F1's to S's size-adjusted power averages 0.88 in the rerun
and 0.83 in the stored run [inferred, summary-full-anchored.csv]. The
stored 29% was therefore larger than F1's typical deficit (17% on
average in the stored run, 12% in the rerun) [inferred,
summary-full-anchored.csv], and part of it is plausibly Monte Carlo
error in a single cell of 500 trials [inferred: the two runs are
independent, and the rerun cannot separate the effect of the DGP change
from chance]. Under compound symmetry F1 is now slightly ahead of S
(0.510 against 0.490; stored 0.522 against 0.554)
[inspected, summary-full-anchored.csv].

**Paired comparisons within the rerun** (McNemar test at 0.05 on paired
size-adjusted rejections, 22 non-prognostic alternative cells):

- F1 against DS: F1 lower in 12 cells and higher in 3 (stored run: 15
  and 0) [inspected, decision-power-ds-nosym-anchored.csv];
- F1 against S: F1 lower in 12 and higher in none (stored: 13 and 0)
  [inspected, decision-power-ds-nosym-anchored.csv];
- ES against S: ES lower in 14 and higher in 5 (stored: 15 and 3)
  [inspected, decision-power-ds-nosym-anchored.csv];
- C1 against S: C1 higher in 20 and lower in none (stored: 20 and 0)
  [inspected, decision-power-ds-nosym-anchored.csv].

**Reading.**

- **The cost of robustness is smaller than the stored run suggested at
  the AR(1) reference cell,** but it remains: F1 has significantly
  lower power than S in more than half of the non-prognostic cells
  [inspected, decision-power-ds-nosym-anchored.csv], and the largest
  deficits are in the measurement-error, dropout, $t_{1/2} = 1$ and
  Weibull 0.5 cells (F1/S ratios 0.71 to 0.80) [inferred,
  summary-full-anchored.csv].
- **Cell-level shifts between the runs are of the size expected from
  Monte Carlo error.** The mean change in unadjusted power over the 11
  analyses ranges from $-0.048$ (compound-symmetry curve at 0.25) to
  $+0.054$ (prognostic baseline cell) over the 26 alternative cells
  [inferred, summary-full-anchored.csv], against an MCSE of about 0.03
  for the difference of two independent cells of 500 trials [inferred,
  Monte Carlo theory]. The changes come in blocks (for example every
  analysis lower at $t_{1/2} = 0$, a cell whose means anchoring does
  not change), because all analyses in a cell share the same trials. The
  $t_{1/2} = 0$ shift can only come from the allocation change or from
  chance; which of the two was not separated (a run with
  `--allocation fixed` would separate them).
- **Choosing the half-life by AIC spreads further.** The share of DS
  fits choosing $t_{1/2}^{a} = 0.5$ falls from 0.299 to 0.250, and the
  shares at 0 and 2 weeks rise (0.170 to 0.216, and 0.179 to 0.213)
  [inspected, summary-full-anchored.csv].

### Fourteen-analysis decision (jobs 21 and 22; added 2026-10-15)

G1 and G2 (simulation-server job 21, `--extra`) and F1A (job 22,
`--extra2`) were rerun on the trials of job 15, so the fourteen
analyses are paired (`replicates-full-extra-anchored.rds`,
`replicates-full-extra2-anchored.rds`). The pairing was checked on the
S arm, which every run refits: 27,997 of the 28,000 S estimates are
identical between job 15 and job 21, and the other 3 are failed fits in
both [verified, replicates-full-extra-anchored.rds]. The decision was
applied with
`12-strawman-stress-decision.R --strawman DS --no-symmetry --with-extras
--tag -anchored`, which writes
`decision-*-ds-nosym-all14-anchored.csv` (committed `a4216f52`); the
post hoc arm summaries are `posthoc-extra-anchored-*.csv` and
`posthoc-extra2-anchored-*.csv` (`14-posthoc-arm-summary.R`).

| Analysis | Pooled size, rerun (stored) | Largest cell, rerun | Prognostic flag | Maximin, non-prognostic, rerun (stored) | Maximin, all cells, rerun (stored) |
|---|---|---|---|---|---|
| G1 | 0.0503 (0.049) | 0.062 (meas. error) | no | 0.236 (0.252) | 0.350 (0.308) |
| G2 | 0.0499 (0.047) | 0.062 (Weibull 0.5) | no | 0.214 (0.232) | 0.342 (0.312) |
| F1A | 0.0592 (0.057) | 0.078 (Weibull 0.5) | no | fails size | fails size |
| F1 | 0.0514 (0.051) | 0.064 ($t_{1/2} = 0$) | no | 0.240 (0.242) | 0.348 (0.358) |
| ES | 0.0490 (0.048) | 0.060 (Weibull 0.5) | no | 0.266 (0.278) | 0.498 (0.544) |

Table: Rerun under the refined process, fourteen analyses: size,
prognostic flag and maximin shortfall of the robust analyses (stored run
in parentheses; the shortfall is measured against the best survivor in
each cell, flagged survivors included).
*Source: decision-maximin-ds-nosym-all14-anchored.csv,
posthoc-extra-anchored-size.csv [inspected]*

**Reading.**

- **Survivors.** The size and sensitivity screens (Amendment 1) leave
  S, C1, C2, ES, F1, F2, G1 and G2. F1A fails size (pooled 0.059,
  outside the band 0.0458 to 0.0542, and 0.078 in the Weibull 0.5 cell
  against the cap of 0.075), as in the stored run; DS and EA fail as in
  the eleven-analysis rerun. Of the survivors, ES, F1, G1 and G2 are
  not flagged in any prognostic null cell (G1 and G2 at most 0.063).
- **Maximin shortfall.** Over the non-prognostic alternative cells the
  order is G2 0.214, G1 0.236, F1 0.240, ES 0.266 [inspected,
  decision-maximin-ds-nosym-all14-anchored.csv]. Under the rule as
  written, G2 leads F1 by 0.026, beyond the tie margin of 0.02, and
  would be selected; G1 ties with F1. Over all cells G2 (0.342), F1
  (0.348) and G1 (0.350) lie within 0.008 of one another [inspected,
  decision-maximin-ds-nosym-all14-anchored.csv]. In the stored run G2
  led F1 by 0.010, within the margin [inspected,
  decision-maximin-ds-nosym-all14.csv].
- **Power at the AR(1) reference cell** ($c_{bm} = 0.30$, size
  adjusted): S 0.504, F1 0.448, G2 0.436, G1 0.426, F1A 0.404
  [inspected, summary-full-anchored.csv and
  posthoc-extra-anchored-power.csv].

A lead of 0.026 is of the order of the Monte Carlo error of a maximin
shortfall at 500 trials per alternative cell [inferred: the MCSE of
one size-adjusted power near 0.4 is about 0.022, and the maximin takes
the worst of 22 such cells]. Rather than select G2 on a margin of that
size, fresh trials were run.

### Fresh trials for the F1, G1, G2 choice (job 37; added 2026-10-15)

Job 37 (`11-strawman-stress-test.R --extra4`, driver committed
`e24d71a6`, data `a3ea8f5f`) fitted the eight survivors S, C1, C2, ES,
F1, F2, G1 and G2 to 1,500 fresh trials in each of the 41 cells (seed
base 47261012, distinct from jobs 15, 21 and 22), under the refined
process (`replicates-full-extra4-anchored.rds`,
`summary-full-extra4-anchored.csv`). There were 104 failed fits in
492,000 [inspected, replicates-full-extra4-anchored.rds]. The run is post hoc: it was made after the decision above had
been seen, to settle a choice that lay within Monte Carlo error.

The unchanged decision script was applied to job 37 alone and to the
trials of job 37 pooled with those of jobs 15 and 21 (job 37 replicate
ids offset by 10,000; 2,500 trials per null cell and 2,000 per
alternative cell). The shortfalls were bootstrapped by resampling trials
within cell, paired across analyses (300 resamples). The scripts are
`15-extra4-pooled-inputs.R`, `15-extra4-bootstrap.R` and
`15-extra4-run.sh` in `analysis/scripts/quick-sim/carryover-closed-form/`;
the decision files and logs are in
`analysis/data/quick-sim/carryover-closed-form/strawman-stress/extra4-decision/`
(`decision-new.txt`, `decision-pooled.txt`, `boot.txt`, and the
decision CSVs under `new/` and `pooled/`). A rerun of `15-extra4-run.sh`
reproduced every output file of the first, scratch computation byte
for byte [verified, extra4-decision/].

| Analysis | Pooled size, job 37 | Pooled size, all trials (largest cell) | Prognostic flag, all trials | Maximin, non-prognostic, job 37 | Maximin, non-prognostic, all trials (bootstrap SE) | Maximin, all cells, all trials |
|---|---|---|---|---|---|---|
| F1 | 0.0487 | 0.0497 (0.057) | no | 0.257 | 0.2155 (0.014) | 0.373 |
| G1 | 0.0482 | 0.0491 (0.056) | no | 0.239 | 0.2205 (0.015) | 0.3585 |
| G2 | 0.0474 | 0.0484 (0.059) | no | 0.245 | 0.2205 (0.017) | 0.359 |
| ES | 0.0495 | 0.0493 (0.057) | no | 0.243 | 0.235 | 0.4385 |
| S | 0.0496 | 0.0504 (0.059) | yes | 0.151 | 0.122 | 0.133 |
| F2 | 0.0504 | 0.0510 (0.059) | yes | 0.157 | 0.1415 | 0.1455 |
| C1 | 0.0487 | 0.0504 (0.060) | yes | 0.000 | 0.0025 | 0.0025 |
| C2 | 0.0458 | 0.0471 (0.051) | yes | fails size | fails size | fails size |

Table: Job 37 and job 37 pooled with jobs 15 and 21: size over the 11
non-prognostic null cells, prognostic flag and maximin shortfall in
size-adjusted power (the shortfall is measured against the best
survivor in each cell; the bootstrap SE is reported for the three
analyses in contention).
*Source: extra4-decision/decision-pooled.txt,
extra4-decision/decision-new.txt [inspected]*

**Reading.**

- **Size.** F1, G1 and G2 hold their size: over the 27,500 pooled null
  trials of the 11 non-prognostic null cells their rates are 0.0497,
  0.0491 and 0.0484 (band 0.0474 to 0.0526), with no cell above 0.059,
  and they stay at or below 0.062 in every prognostic null cell (ES at
  most 0.063) [inspected, extra4-decision/decision-pooled.txt]. C2
  now fails the size screen as conservative (pooled 0.0471, below the
  band; 0.0458 in job 37 alone) [inspected,
  extra4-decision/decision-pooled.txt], so the pooled survivors are S,
  C1, ES, F1, F2, G1 and G2.
- **Maximin shortfall.** On the fresh trials alone the order is G1
  0.239, ES 0.243, G2 0.245, F1 0.257: G1, G2 and F1 lie within 0.018
  [inspected, extra4-decision/decision-new.txt]. On all trials pooled
  it is F1 0.2155, G1 0.2205, G2 0.2205, ES 0.235, with F1 now first
  by 0.005 [inspected, extra4-decision/decision-pooled.txt]. The three
  worst cells are all the measurement-error cell. Over all cells G1
  (0.3585) and G2 (0.359) lead F1 (0.373) by 0.014, also within the
  margin [inspected, extra4-decision/decision-pooled.txt].
- **Monte Carlo error.** The paired bootstrap SE of the F1 to G1, F1 to
  G2 and G1 to G2 differences in maximin shortfall is 0.012 to 0.016
  on the pooled trials (0.011 to 0.017 on job 37 alone)
  [inspected, extra4-decision/boot.txt]. The tie margin of 0.02 is
  therefore not much larger than one SE of the difference at 2,000
  trials per cell, and smaller than its SE at 500 trials per cell
  (about twice as large). This is why the first rerun appeared to favor
  G2: its lead of 0.026 does not recur on fresh trials [inferred from
  the two independent sets of trials and the bootstrap SE]. (The
  bootstrap takes the best analysis in each cell over all eight,
  C2 included, so its ES shortfall, 0.2455 pooled, differs from the
  decision script's 0.235 [inspected, extra4-decision/boot.txt]; the
  F1, G1 and G2 values agree.)
- **Power at the AR(1) reference cell** ($c_{bm} = 0.30$, all trials
  pooled, 2,000 per cell): S 0.496, G2 0.444, F1 0.431, G1 0.426
  [inspected, extra4-decision/decision-pooled.txt]; F1 has 13% less
  size-adjusted power than S (11% on the 500 trials of job 15 alone)
  [inferred, extra4-decision/decision-pooled.txt and
  summary-full-anchored.csv].

### Decision (2026-10-15)

Under the rule of Section 7 on the pooled trials, F1, G1 and G2 lie
within the tie margin of 0.02 and the tie-break (fewer covariance
parameters) selects **F1**. F1 is retained as the recommended analysis
for the biomarker-treatment interaction under the refined process. G1
and G2 are post hoc alternatives that tie with it; the lead of G2 on
the 28,000 rerun trials is attributed to Monte Carlo error. Because job
37 was run after the first decision had been seen, and the margin is
small against Monte Carlo error, the order of F1, G1 and G2 is not
established by this study; an independent study with more trials per
cell would be needed to separate them.

Not done: paired McNemar comparisons of F1 against G1 and G2 in each
cell (the decision script runs them only against a surviving
strawman); and the 1,500 fresh trials did not include F1A, DS, EA, A,
DA or R, which had already failed a screen.

## 12. References

1. Morris TP, White IR, Crowther MJ. Using simulation studies to
   evaluate statistical methods. *Statistics in Medicine*
   2019;38(11):2074-2102. doi:10.1002/sim.8086.
2. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020;2:13. doi:10.3389/fdgth.2020.00013.
3. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology*
   2002;28(2):169-181.
4. Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust
   variance estimation and hypothesis testing in fixed effects models.
   *Journal of Business and Economic Statistics* 2018;36(4):672-683.
   doi:10.1080/07350015.2016.1247004.
5. Kenward MG, Roger JH. Small sample inference for fixed effects from
   restricted maximum likelihood. *Biometrics* 1997;53(3):983-997.
   doi:10.2307/2533558.
6. pmsimstats team. `docs/44`, `docs/45`, `docs/46`, `docs/48`,
   `docs/49`.


---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at docs/series/registry/open-items.csv -->

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
