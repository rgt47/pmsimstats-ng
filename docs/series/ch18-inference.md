# Inference and calibration

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** Random-intercept absorption; Type I error decomposition;
KR vs CR2; Satterthwaite df.  
**Prerequisites:** ch17

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

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: CM "Inference" (6 subsections); 45 Section 4; 33 Section
  2.4; 56 Section 3.5; 36 Section 6.8; Appendix I Section 4 (the
  `lmer` evidence); Appendix H (RM-ANOVA).

- New to write: a kappa section (33 Section 2.4 has the only
  analysis); the 0.14 to 0.16 against 0.14 to 0.23 ranges labeled
  by construct (CR-28).
-->

### Type I error decomposition

<!-- From analysis/report/compendium-mathematics.tex lines 7625-7722
(Section "Type I error decomposition") -->

The rejection rate under the null $\beta_{bm:D} = 0$ decomposes as

$$
\alpha_{\text{actual}} = \alpha_{\text{nominal}}
  + \delta_{\text{size}},
$$

where $\delta_{\text{size}}$ is the miscoverage attributable to
misspecified inference. Write $\kappa$ for the ratio of the
model-based standard error to the empirical standard deviation of the
estimate across replications:

$$
\kappa = \frac{\overline{\widehat{\mathrm{SE}}(\hat\beta)}}
              {\widehat{\mathrm{SD}}(\hat\beta)}.
$$

A value below 1 means the model-based standard error is
anticonservative; a value at 1 is correctly sized. The Wald test
statistic is $\hat\beta / \widehat{\mathrm{SE}}$, so the actual size is
approximately $2\Phi(-z_{1-\alpha/2}/\kappa)$ for a two-sided test
at nominal $\alpha$. At $\kappa = 0.78$, $\alpha = 0.05$, the actual
size is about 0.10 (`docs/45`, Section 4; inspected).

<!-- TODO: kappa section -- docs/33 Section 2.4 is the only source
for the kappa analysis; that document does not exist in the project.
Add the derivation and numerical values when docs/33 is located or
reconstructed. -->

### What the random intercept absorbs under AR(1) data

<!-- From analysis/report/compendium-mathematics.tex lines 7724-7808
(Section "What the random intercept absorbs") -->

Suppose the true residual structure is AR(1): the within-participant
covariance is $\sigma^2 \rho^{|t_1 - t_2|}$. The random-intercept
model posits $\mathrm{Cov}(Y_{it}, Y_{it'}) = \sigma_u^2 + \sigma_e^2
\mathbf 1_{t=t'}$. The two structures agree at

$$
\sigma_u^2 = \sigma^2 \rho^{|t - t'|} \text{ only when } |t-t'| = 0,
$$

which forces $\sigma_u^2 = 0$ and maps the AR(1) covariance to white
noise. For the Hybrid at $\rho = 0.7$, $T = 9$ timepoints, the best
compound-symmetry fit to the AR(1) Gram matrix yields REML estimates
$\hat\sigma_u^2 \approx 0$ and inflated $\hat\sigma_e^2$ (inspected;
`docs/45` Section 4.3). The Kenward-Roger correction is designed for
balanced data with a correctly specified covariance; applied here it
does not fix the mismatch between compound symmetry and AR(1).

### Kenward-Roger is not interchangeable with CR2

<!-- From docs/45 Section 4 -->

`docs/45` compared five inference procedures on the same simulated
trials (500 to 1000 replicates; Hybrid design; CS and AR(1)
constructs; $N = 70$):

| Inference | CS data | AR(1) data |
|---|---|---|
| `us-KR` | 0.054 (Hybrid) | 0.062 (Hybrid) |
| `ri-CAR1-KR` | 0.056 | 0.062 |
| `ri-CAR1-CR2` | 0.051 | 0.054 |
| `ri-CAR1-Satt` | 0.053 | 0.061 |
| `us-CR2` | 0.052 | 0.052 |

(Rounded; verified from `docs/45` Table 4, CO column omitted here.)

KR with either the unstructured or the random-intercept-plus-CAR(1)
covariance holds size at 0.054 to 0.062 under CS data but inflates to
0.062 in Hybrid under AR(1). CR2 holds 0.051 to 0.054 in both
constructs. Under the crossover, `ri-CAR1-KR` fails badly: 0.14 to
0.23 (inspected; construct-specific ranges in `docs/36` Section 6.8).
KR and CR2 are not interchangeable once the data are AR(1).

### CR2 intuition, plain sandwich, and degrees of freedom

<!-- From docs/56 Section 3.5 -->

The plain sandwich estimator is consistent under any working
covariance but converges slowly: with only $N$ independent clusters
(participants) the denominator degrees of freedom are at most $N - p$,
and at small $N$ the plain sandwich is liberal even relative to
KR. CR2 (Bell-McCaffrey bias-reduced sandwich) inflates each
cluster's score by $(I - H_i)^{-1/2}$ before squaring, where $H_i$
is the hat matrix block for participant $i$. The effect is to correct
for the finite-sample downward bias in the plain sandwich.

Satterthwaite degrees of freedom for the CR2 test statistic are
derived by matching the first two moments of the $F$ distribution
to the empirical distribution of the Wald ratio. This is a
moment-matching approximation and does not require balance or equal
cluster sizes.

**CR2 in the crossover.** Under the crossover at 1 to 5 MCSE, `us-KR`
runs 0.14 to 0.23 under the published compound-symmetry construct with
AR(1) data (`docs/36`, Section 6.8). `ri-CAR1-CR2` under the same
runs 0.14 to 0.16 (verified from `docs/56`, Section 3.5). The gap
narrows but does not close; using AR(1) working covariance is the
correct fix.

<!-- TODO: table distinguishing 0.14 to 0.16 (ri-CAR1-CR2, CO,
AR(1) data) from 0.14 to 0.23 (us-KR, CO, AR(1) data), with
construct labels, N, and MCSE bands (CR-28) -->

### Evidence from `lmer`-based inference

<!-- From analysis/report/compendium-mathematics.tex Appendix I
Section 4 -->

The parallel analysis using `lmer` (R package `lme4`) with
`lmerTest::summary()` replicates the Satterthwaite df for the
KR-adjacent inference. The relevant finding for this chapter is that
`lme4`'s REML-KR approximation agrees closely with `nlme::lme`
KR at sample sizes $N \ge 50$. At $N = 30$ the two can diverge by
up to 0.005 in rejection rate (inspected; Appendix I Section 4).

### Why the published analysis fails in the crossover under AR(1)

<!-- From docs/45 Section 4.2 -->

The crossover has no between-path comparison: every participant sees
both drug and placebo, so $D_{it}$ varies only within participant.
The random intercept absorbs all between-participant variation, leaving
the drug effect and the interaction identified entirely by
within-participant contrasts. When the data are AR(1), the within-
participant residuals are serially correlated, and a white-noise
working covariance underestimates the variance of within-participant
contrasts. The result is the 0.14 to 0.23 inflation.

The Hybrid design mixes between- and within-participant information,
so the inflation is smaller (0.062 for `ri-CAR1-KR` against 0.054
for `ri-CAR1-CR2`; verified from `docs/45` Table 4).

### RM-ANOVA as special case

See Appendix H (`docs/series/appH-rm-anova.md`) for the derivation
of RM-ANOVA as a restricted case of the MMRM; the size inflation
under the RM-ANOVA's sphericity assumption when the true structure is
AR(1) follows the same mechanism as above.

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
| 2026-10-15 | Setting and argument | Initial content move from CM, docs/45, docs/56, docs/36 | Phase 2 content moves |
