# Three defects of the published process

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** Three defects: ceiling violation, step coupling, compounding
recursion.  
**Prerequisites:** ch03

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

*No Monte Carlo simulation results in this chapter.*

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 36 Sections 1, 4, 5 (spine); 55 Sections 1, 5 to 8
  (summary only); 53 Sections 1, 8 (summary only); 37 Section 5.1;
  34 Section 4.6; 50 Section 1; 07 counts via Appendix I.

- New to write: one table reconciling the published power at
  $t_{1/2} = 0.1$ week (0.126, 0.078, 0.100, 0.076: cell, code,
  seeds), from 55 Section 7.3; one reference cell quoted throughout
  the series.
-->

<!-- TODO: one table reconciling the published power at $t_{1/2}
     = 0.1$ week: the four values 0.126, 0.078, 0.100, 0.076 as
     read from cell, code, and different seed sets; from docs/55
     Section 7.3; this is the reference cell quoted throughout the
     series -->

### The three defects: a summary

<!-- From docs/legacy/50 Section 1 "The data-generating process" -->

Hendrickson et al. (2020; "RH") simulate the score as baseline
minus three Gompertz responses: drug (BR), expectancy (PB) and
natural course (TV). The biomarker is bound to the drug response
as a correlation $c_{bm}$ in a joint multivariate normal, a
parameterization new to the N-of-1 literature. Two defects follow
from the published implementation; a third is smaller.

| Defect | Consequence | Evidence | Fix |
|---|---|---|---|
| Biomarker coupled at full strength wherever any drug effect remains (step rule) | interaction vanishes at any positive half-life | Hybrid, $N = 70$, $c_{bm} = 0.25$, no censoring: power 0.66 without carryover to 0.08 at $t_{1/2} = 0.1$ (RH's code, 250 replicates per cell) | **graded coupling**: $c_{bm}$ on drug, $c_{bm}2^{-t_{sd}/t_{1/2}}$ off drug, 0 before exposure; power 0.62 at 0.1 week in the same run |
| Compound symmetry allows $c_{bm}$ only below 0.256; larger values silently repaired | RH's $c_{bm} = 0.6$ results describe a different process | 20 of 54 matrices repaired | **separable AR(1)** ($\rho = 0.7$): ceiling 0.48; **stop, do not repair** |
| Carryover recursion compounds the decay | off-drug means too small at half-lives of 0.5 weeks or longer | code inspection and paired simulation | **anchored carryover** from the last on-drug value |

Table: Three defects of the published process with consequences, evidence and fixes (docs/legacy/50, Section 1)

**Key insight.** Conditional on the biomarker, covariance binding
is mean moderation with a smaller residual variance. How much
carryover costs therefore depends on whether the interaction fades
with the drug effect, not on where it is bound. Mean binding with
the same fading profile loses power almost identically: 0.454
against 0.473 at one week (docs/legacy/50, Section 1).

### Defect 1: correlation matrices that are not positive definite

<!-- From docs/36 Sections 1 and 4 -->

**What fails.** Running the published `generateData()` at the
published parameters, 40 of the 162 design-path-cell matrices fail
the code's own positive definiteness check (verified, docs/36).
The failures are systematic: every matrix with $c_{bm} > 0$ and
no carryover in the CO, OL+BDC and Hybrid designs, plus the
placebo-first CO path at every half-life. The response block itself
is always valid.

**Why.** The failures come because the requested
biomarker-response correlation exceeds the largest value the
compound-symmetry (CS) structure can carry, 0.256 to 0.262 in
these designs, well below the published values of 0.3 and 0.6. The
exact ceiling is

$$
c_{bm}^{\ast} = \bigl(\tilde g^{\top} M^{-1} \tilde g\bigr)^{-1/2},
$$

where $M$ is the response block and $\tilde g$ the coupling pattern
over the $BR$ visits, padded with zeros at the $TV$ and $PB$
entries. For the published construct (verified):

| Design | $t_{1/2} = 0$ | $t_{1/2} = 0.1$ | $t_{1/2} = 0.2$ |
|---|---|---|---|
| OL | 0.893 | 0.893 | 0.893 |
| OL+BDC | 0.262 | 0.893 | 0.893 |
| CO | 0.256 | 0.256 | 0.256 |
| Hybrid | 0.256 | 0.893 | 0.893 |

Table: Ceiling $c_{bm}^{\ast}$ for the published construct by
design and half-life (docs/36, Section 4.2)

The ceiling is high when the coupling pattern is constant across
visits (OL always; OL+BDC and Hybrid with any carryover under the
step rule) and low when it has an on/off contrast (no carryover, or
CO's placebo-first path). Both published effect sizes (0.3 and 0.6)
exceed the 0.256 ceiling wherever it applies.

**What the repair does.** The code replaces the covariance matrix
by a nearby positive definite one when the check fails. The effect
on the correlations actually simulated (verified, docs/36,
Section 4.3):

| $c_{bm}$ specified | On-drug correlation used | Off-drug correlation used | Largest change elsewhere |
|---|---|---|---|
| 0.3 | 0.298-0.300 | 0.000-0.001 | 0.007 |
| 0.6 | 0.540-0.568 | 0.029-0.042 | 0.097 |

Table: Correlations simulated after the published positive
definiteness repair, by specified $c_{bm}$ (docs/36, Section 4.3)

At $c_{bm} = 0.3$ the repair is cosmetic. At $c_{bm} = 0.6$ the
simulated effect is materially smaller than specified, some
interaction signal is moved into off-drug visits, and the within-
and cross-factor correlations are altered. The published
no-carryover results at $c_{bm} = 0.6$ describe a DGP other than
the one stated.

### Defect 2: step coupling and the power collapse

<!-- From docs/55 Sections 1, 5 to 8 and docs/37 Section 5.1 -->

**What the code does.** Lines 131 to 138 of `R/generateData.R` at
`58b32a9` set the correlation between the biomarker and the drug
response at visit $t$ to $c_{bm}$ if the expected drug response at
that visit is nonzero, and leave it at 0 otherwise (inspected).
The expected drug response used in the test already includes
carryover. For any $t_{1/2} > 0$ it is positive at every off-drug
visit after the first exposure, so those visits receive the full
on-drug correlation (derived; verified in all 40 path and half-life
cases).

**The consequence.** With full coupling on and off drug, the
coupling gap ($\bar g_{\text{on}} - \bar g_{\text{off}}$) is 0 and
the slope of the on-minus-off contrast on the biomarker is exactly
0. Power is a step function of the half-life: its no-carryover
value at $t_{1/2} = 0$, and a constant at every positive half-life,
whether 0.01 weeks or 10 (derived, docs/37, Section 5.1). The
constant for E9 is the nominal size, 0.050, because the slope is
exactly zero; for E9b it is 0.081, carried by the baseline row
alone. The magnitude of the half-life plays no role.

**The closed-form evidence.** Power for the paired-difference
statistic E9 (and E9b, the proxy for the published analysis) by
coupling rule (verified, docs/37, Section 5):

| Rule | Variant | $t_{1/2}$ = 0 | 0.1 | 0.2 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|---|---|
| Step (published) | E9 | 0.721 | 0.050 | 0.050 | 0.050 | 0.050 | 0.050 |
| Step (published) | E9b | 0.654 | 0.081 | 0.081 | 0.081 | 0.081 | 0.081 |
| Graded | E9 | 0.721 | 0.721 | 0.712 | 0.631 | 0.473 | 0.278 |
| Graded | E9b | 0.654 | 0.653 | 0.647 | 0.586 | 0.472 | 0.327 |
| Mean moderation | E9 | 0.681 | 0.681 | 0.682 | 0.682 | 0.682 | 0.682 |
| Mean moderation | E9b | 0.619 | 0.619 | 0.620 | 0.620 | 0.620 | 0.620 |

Table: Closed-form power by coupling rule, variant, and
$t_{1/2}$, Hybrid design, $N = 70$, $c_{bm} = 0.25$, CS construct
(docs/37, Section 5.1)

**The published results confirm the pattern.** Hybrid power at
$N = 70$, $c_{bm} = 0.3$ falls from 0.74 to 0.12 at a half-life
of 0.1 weeks (published saved results, verified, docs/36,
Section 5.1). A half-life of 0.1 weeks leaves a residual drug
effect of 0.1% at the first off-drug visit (one week after
stopping), yet the power falls by 84 percent. The carryover is
negligible; the coupling rule is not.

**The mean estimate shows the same attenuation.** Mean published
Hybrid interaction estimate at $N = 70$ (verified, docs/36):

| $c_{bm}$ | $t_{1/2} = 0$ | 0.1 | 0.2 |
|---|---|---|---|
| 0.3 | $-0.151$ | $-0.036$ | $-0.036$ |
| 0.6 | $-0.279$ | $-0.065$ | $-0.081$ |

Table: Mean published Hybrid interaction estimate at $N = 70$ by
$c_{bm}$ and half-life (docs/36, Section 5.3)

The estimate shrinks toward zero by roughly three quarters. The
loss of power is primarily attenuation of the effect being
estimated, not only inflation of its standard error.

**Graded coupling restores power.** In a paired simulation with 500
replicates per cell (verified, docs/53, Section 8.3):

| $t_{1/2}$ (weeks) | Published | Anchored only | Graded | Graded and anchored |
|---|---|---|---|---|
| 0 | 0.628 | 0.628 | 0.628 | 0.628 |
| 0.1 | 0.078 | 0.078 | 0.658 | 0.658 |
| 0.2 | 0.082 | 0.082 | 0.670 | 0.670 |
| 1 | 0.088 | 0.088 | 0.486 | 0.484 |
| 2 | 0.078 | 0.072 | 0.344 | 0.344 |

Table: Power of the published analysis for the interaction, Hybrid
design, $N = 70$, $c_{bm} = 0.25$, 500 replicates per cell, by
version of the 58b32a9 data generator (docs/53, Section 8.3;
verified; MCSE at most 0.022)

Correcting the mean alone (anchored only) leaves the published
analysis essentially unchanged. Graded coupling restores power
to near its no-carryover value at short half-lives and produces a
gradual monotone decline at longer ones.

### The common cause

<!-- From docs/36 Section 5.4 "The common cause" -->

| | No carryover ($t_{1/2} = 0$) | Any carryover ($t_{1/2} > 0$) |
|---|---|---|
| Coupling pattern (Hybrid, OL+BDC) | on/off contrast | constant |
| Matrix at $c_{bm} = 0.3, 0.6$ | invalid, repaired | valid |
| Interaction identified | yes | no |

Table: Coupling pattern, matrix validity and identification with
and without carryover (docs/36, Section 5.4)

The step rule makes the coupling a discontinuous function of the
half-life: any positive value, however small, flips the pattern
from contrast to constant. Problem 1 (invalid matrices) and
problem 2 (power collapse) are both consequences of the step rule.
Problem 1 is the price of the contrast at $t_{1/2} = 0$, and
problem 2 is the loss of the contrast at $t_{1/2} > 0$.

CO is affected less because its placebo-first path retains a zero
coupling before first exposure, which preserves some contrast (and
keeps that path's matrix invalid at every half-life). OL has no
off-drug visits, so its test does not rest on a drug contrast in
the first place.

### Defect 3: the compounding carryover recursion

<!-- From docs/53 Sections 1 and 8 -->

**What the code does.** At `58b32a9` line 88, the expected drug
response at an off-drug visit is set to the previous visit's
expected response times $2^{-t_{sd}/t_{1/2}}$, where $t_{sd}$ is
the *cumulative* time since the drug was stopped (inspected,
docs/53, Section 1). Because the previous visit's value has
already been decayed, elapsed time is counted again at every visit.

**What it should do.** The article states that "the value of the
BR at the last timepoint on active drug was exponentially decayed,
with the half-life of this decay being maintained as a model
parameter." That is $\mu_{\text{last}}\,2^{-t_{sd,m}/t_{1/2}}$
(the anchored form).

**How large is the discrepancy.** At the half-lives of the
article's power results (0.1 and 0.2 weeks), the entire carryover
is nearly gone after one week under either form, and the largest
shortfall in any off-drug mean is 0.0104 points. The discrepancy
is therefore numerically negligible for the published power
analyses (verified, docs/53, Section 8.1). It grows with the
half-life and with the length of the off-drug run. On Hybrid
path 1 at $t_{1/2} = 1$ week, the mean two weeks after stopping is
1.273 points where the intended value is 2.547 (docs/53,
Section 1).

**What it does not touch.** The variances and the correlation
matrix are unchanged (inspected, verified). A paired simulation
shows no detectable change in power from correcting the mean alone
at short half-lives (Table 8, docs/53, Section 8.3). This is the
third defect in severity: it is a real implementation error, but
its effects on the power results the article reports are
negligible.

### The ceiling reversal under step coupling

<!-- From docs/34 Section 4.6 "The ceiling reversal under step
     coupling" -->

Under step coupling, the ceiling reverses between no-carryover and
any-carryover cells. For Hybrid and OL+BDC, at $\rho = 0.7$:

- At $t_{1/2} = 0$: ceiling 0.346 (an on/off contrast is present).
- At any $t_{1/2} > 0$: ceiling 0.841 (the coupling is constant).

The mechanism is that once any carryover is present, step coupling
assigns the full $c_{bm}$ at every post-exposure occasion, making
$g$ constant on paths that begin on drug. A constant vector aligns
with the leading eigenvector of a compound-symmetric block, raising
the feasible ceiling (docs/34, Section 4.6, verified).

Two qualifications:

- **The reversal does not occur in CO.** CO's placebo-first path
  has pre-exposure zeros in $g$ at every half-life, so $g$ never
  becomes constant.
- **The large apparent gain depends on compound symmetry.** With
  step coupling but an AR(1) block, the ceiling at positive
  half-lives is 0.519 rather than 0.841. The large gain in the
  published configuration requires both the step rule and the
  compound-symmetric block.

As docs/36 argues, the higher ceiling at positive half-lives is not
an advantage: it is available precisely because the coupling vector
has stopped carrying any on-drug versus off-drug contrast. Both
defects share one cause. The step rule: without carryover the
coupling has a contrast and the matrix is invalid; with any
carryover the coupling is constant and the interaction disappears.

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
| 2026-10-15 | Setting and argument | Initial content moves from docs/36, docs/55, docs/53, docs/37, docs/34, docs/legacy/50 | Phase 2 content move |
