---
geometry: margin=1.4cm
fontsize: 10pt
---

# Simulating and Testing the Biomarker-by-Drug Interaction in the Hybrid Design: Where We Are
*2026-10-07 11:06 PDT*

*Archived 2026-10-11; current treatment: the executive summary `series/fm1`; see `docs/legacy/INDEX.md`.*

**Author.** pmsimstats team. **Scope.** The recommended data-generating
process (DGP), called the refined process, the recommended analysis,
and the evidence for each.
Details: `docs/36` to `docs/38`, `docs/44` to `docs/49`, and paper 01.
**Symbols.** $c_{bm}$ biomarker coupling to the drug response; $t_{1/2}$
carryover half-life (weeks); $t_{sd}$ time since discontinuation;
$\rho$ AR(1) correlation per week; $c_1$ cross-component correlation;
$N$ total sample size.

## 1. The data-generating process

**Starting point.** Hendrickson et al. (2020; "RH") simulate the score
as baseline minus three Gompertz responses: drug (BR), expectancy (PB)
and natural course (TV). The biomarker is bound to the drug response as
a *correlation* $c_{bm}$ in a joint multivariate normal, a
parameterization new to the N-of-1 literature to our knowledge. Two
defects follow from the published implementation.

| Defect in RH's code | Consequence | Evidence | Fix |
|---|---|---|---|
| Biomarker coupled at full strength wherever any drug effect remains (step rule) | interaction vanishes at any positive half-life | Hybrid, $N = 70$, $c_{bm} = 0.25$, no censoring: power 0.66 without carryover to 0.08 at both $t_{1/2} = 0.1$ and 0.2 (RH's code, Figure 4B reproduction, 250 replicates per cell) | **graded coupling**: $c_{bm}$ on drug, $c_{bm}2^{-t_{sd}/t_{1/2}}$ off drug, 0 before exposure; power 0.62 at 0.1 and 0.68 at 0.2 in the same run |
| Compound symmetry allows $c_{bm}$ only below 0.256; larger values silently repaired | RH's $c_{bm} = 0.6$ results describe a different process | 20 of 54 matrices repaired | **separable AR(1)** ($\rho = 0.7$): ceiling 0.48; **stop, do not repair** |
| Carryover recursion compounds the decay | off-drug means too small | code inspection | **anchored carryover** from the last on-drug value |

*(Clarified 2026-10-09: the first row quotes the 250-replicate
Figure 4B reproduction of paper 01, `fig4b-graded-cbm0.25-reps250.csv`.
`docs/51` quotes a separate 500-replicate run with other seeds
(`docs/37`, Appendix A), in which published power at 0.1 week is 0.13
and graded power 0.64, against 0.67 without carryover. The two runs
describe the same comparison and differ within Monte Carlo error
except in the published cell at 0.1 week, 0.08 against 0.13.)*

**How we got there.** First, closed-form power for a paired-difference
statistic showed that power depends on carryover only through the
*coupling gap* (on minus off coupling). The step rule sets that gap to
0; graded coupling lets it fall smoothly, and provably monotonically
(docs/37). Next, the eight-line fix restored power in RH's own code.
Then the separable form was shown to be valid for every visit schedule,
to give an exact ceiling as a sum of per-switch costs, and to match the
random-intercept-plus-AR(1) analysis model (docs/35).

**Key insight.** Conditional on the biomarker, covariance binding is
mean moderation with a smaller residual variance. How much carryover
costs therefore depends on whether the interaction fades with the drug
effect, not on where it is bound. Mean binding with the same fading
profile loses power almost identically: 0.454 against 0.473 at one week.

**The refined process (recommended DGP).**

- **Core:** graded coupling, separable AR(1) with $\rho = 0.7$ and
  $c_1 = 0.2$, anchored carryover, a stop on infeasible matrices.
- **Re-exposure (retained from the published process):** the
  drug-response curve restarts at zero, as in Hendrickson's code; this
  is consistent with re-titration of prazosin.
- **Carryover range:** pharmacological half-life 0.1 to 0.2 week.
- **Options:** prognostic biomarker couplings (natural course, placebo,
  baseline severity), phase-specific measurement error, Weibull decay,
  dropout.

Mean binding, or the drug-sensitivity model of `docs/38` (not yet
implemented), covers effects above the ceiling.

## 2. The analysis

**How we got there.** Three rounds of simulation, each on identical
trials, with size judged before power:

1. **Covariance study** (`docs/45`, 378,000 fits). Kenward-Roger
   corrects an unstructured covariance; CR2 corrects a parsimonious one.
   Neither substitutes for the other.
2. **Exploratory stress test** (`docs/46`, 41 cells, 28,000 trials,
   11 protocol analyses). It varied correlation structure, $N$,
   carryover, decay form, dropout and measurement error, and added a
   *prognostic* biomarker.
3. **Post hoc arms** on the same trials, judged by the same rule
   (`docs/46`, Amendment 2): phase-specific variances with the
   expectancy terms (G1, G2), AIC selection with them (F1A), and two
   baseline-severity cells.

The decisive finding concerns analyses that keep baseline as an outcome
row. They share one biomarker slope between baseline and the off-drug
visits. If the biomarker also predicts the natural course or the
placebo response, part of that association is reported as a drug
interaction (`docs/48`, `docs/49`).

| Analysis (random intercept + CAR(1) residuals unless stated) | Size, standard cells | Size, prognostic cells | Power, AR(1) $c_{bm} = 0.3$ | Power, CS $c_{bm} = 0.25$ |
|------------------------------------|---------------|-----------|---------|---------|
| Published (`mmrm` compound symmetry, Satterthwaite df) | 0.069 (0.081 at $N = 35$) | 0.053 to 0.159 | 0.398 | 0.596 |
| Unstructured, Kenward-Roger | 0.057 (0.087 at $N = 35$) | 0.195 to 0.416 | 0.534 | 0.618 |
| S: fixed half-life 0.5, CR2 | 0.048 | 0.049 to 0.100 | 0.514 | 0.554 |
| DS: S with AIC-selected half-life | 0.053 | 0.058 to 0.115 | 0.486 | 0.480 |
| C1: S with phase-specific variances | 0.053 | 0.173 to 0.382 | 0.578 | 0.660 |
| **F1: S plus expectancy terms `De + bmc:De`** | **0.051** | **0.040 to 0.047** | **0.364** | **0.522** |
| G1: F1, phase variances, CR2 (post hoc) | 0.049 | 0.042 to 0.055 | 0.354 | 0.546 |
| G2: G1, model-based (post hoc) | 0.047 | 0.043 to 0.062 | 0.380 | 0.528 |
| F1A: F1, AIC half-life (post hoc) | 0.057 (0.073 at $t_{1/2} = 1$) | 0.047 to 0.059 | 0.338 | 0.446 |

*Size: null rejection, nominal 0.05 (standard cells pooled over 11,000
trials; band 0.046 to 0.054). Power is size-adjusted, 500 trials per
cell. "Prognostic": the biomarker correlated 0.3 with the natural
course (constant or growing) or with the placebo response. In a
baseline-severity null cell every analysis is within the 0.075 cell
limit (S 0.050, F1 0.040).*

**Reading the table.**

- **Phase-specific variances (C1)** are the most powerful, but only
  because they weight the shared baseline slope heavily. Once the
  baseline is freed (G1, G2), the power advantage disappears and the
  robustness is kept. On maximin shortfall G2 (0.232) and F1 (0.242)
  tie within 0.02, and the rule gives the tie to F1, which has fewer
  covariance parameters.
- **AIC selection (DS)** holds its size and estimates the carryover
  half-life sensibly: it picks 0 82% of the time when there is none, and
  1 or 2 weeks 86% of the time when it is 1. It costs about 0.03 power
  against a fixed 0.5 week. It does not protect against a
  prognostic biomarker by itself, and with the expectancy terms (F1A)
  it fails the size criterion.
- **The expectancy terms are what make an analysis robust.** Their cost
  is about 30% of power under strong serial correlation, and small under
  compound symmetry.

**Recommended analysis.**

- **Primary:** F1. Random intercept with `corCAR1` residuals, CR2
  standard errors, linear time, baseline as an outcome row, the
  exposure-decayed drug indicator at a fixed half-life of 0.5
  week (not AIC-selected, since F1A fails size), and `De + bmc:De`.
- **Sensitivity analyses:** S (biomarker assumed purely predictive),
  the freed-baseline version, and the binary drug indicator.
  Disagreement between F1 and S is evidence of a prognostic biomarker.
- **Sample size:** plan it for F1, stating the serial-correlation
  assumption.

**Open items.** The refinements still have to move into the package,
and its sampler's silent repair has to go. The covariance study and
the stress test ran configuration C of `docs/36` at $\rho = 0.7$, the
covariance construct of the refined process, with the published
recursive carryover and fixed allocation, so the refined process with
anchored carryover has not itself been simulated (inspected in
`08-covariance-study.R` and `11-strawman-stress-test.R`). Effect sizes are not yet
calibrated to prazosin data: blood pressure's placebo-arm association
was null in only 35 participants (Raskind et al. 2016). Any post hoc
choice (G1, G2) needs an independent run.
