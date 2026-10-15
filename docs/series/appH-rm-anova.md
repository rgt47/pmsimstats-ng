# Classical test equivalents: repeated-measures ANOVA

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** Classical repeated-measures ANOVA equivalents for the interaction test.  
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

*Local symbols used in this chapter only. For project-wide notation see
`analysis/report/NOTATION.md`.*

<!-- notation entries -->

---

## Summary

<!-- evidence-labeled results go here -->

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N, analysis ID,
t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 26 (whole; reserved letters renamed).

- New to write: none beyond renaming.
-->

<!-- Spine: docs/26-rm-anova-interaction-test.md (copied verbatim; Phase 2) -->

# Strict classical repeated-measures ANOVA for the biomarker-treatment interaction: closed-form derivation and limitations {.unlisted .unnumbered}

*2026-05-07 06:30 PDT*

**Author.** pmsimstats team

**Purpose.** This note derives, in closed form, the strict classical
repeated-measures ANOVA (RM-ANOVA) test of the biomarker-treatment
interaction for a dichotomized biomarker, states its power formula,
and sets out how it relates to the mixed-model Wald test that
pmsimstats-ng uses. The derivations are analytic; the few numerical
statements are checked against stored simulation output or flagged.
The simulation study this note proposed was later carried out as
paper 08 (`analysis/report/08-test-procedure-design-sensitivity/`),
and dated notes in Sections 6 and 8 point to its results.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. The classical
split-plot symbols are *local*; they were renamed on 2026-10-07 so
that none collides with a canonical symbol (formerly $\alpha_i$,
$\beta_j$, $\pi_{l(i)}$, $e_{ijlk}$, $\delta$, $\lambda$ and
$1 - \beta$).

| Symbol | Meaning | Status |
|---|---|---|
| $g \in \{1, \ldots, a\}$ | biomarker stratum (factor A); $a = 2$ in the main case | local |
| $j \in \{1, \ldots, b\}$ | treatment phase (factor B); $b = 2$ in the main case | canonical index |
| $i \in \{1, \ldots, n\}$ | participant within stratum | canonical index |
| $t \in \{1, \ldots, T\}$ | within-phase timepoint | canonical index |
| $n$, $N$ | participants per stratum; total, $N = an$ | local; canonical |
| $Y_{gjit}$ | symptom score (code `Sx`) | canonical |
| $\mu$, $\xi_g$, $\phi_j$, $\tau_t$ | grand mean, stratum, phase and time effects | local |
| $(\xi\phi)_{gj}$, $(\xi\tau)_{gt}$, $(\phi\tau)_{jt}$, $(\xi\phi\tau)_{gjt}$ | interaction effects | local |
| $u_{i(g)}$, $\sigma_u^2$ | participant random intercept, nested in stratum, and its variance | canonical |
| $\varepsilon_{gjit}$, $\sigma_\varepsilon^2$ | within-participant residual and its variance | canonical |
| $\bar Y$ with dots | means over the dotted subscripts | local |
| $SS$, $MS$, $df$ | sums of squares, mean squares, degrees of freedom | local |
| $L$, $\hat L$ | population and sample interaction contrast, eq. (9) | local |
| $L_1$ | value of $L$ under the alternative | local |
| $\Delta_{gi}$, $\bar\Delta_g$ | within-participant phase difference and its stratum mean | local |
| $F_{AB}$, $t_{AB}$ | interaction $F$- and $t$-statistics | local |
| $\epsilon_{\mathrm{Box}}$, $\hat\epsilon_{\mathrm{Box}}$ | Box's sphericity measure and its estimate | local |
| $\Lambda$ | non-centrality parameter of the $F$ distribution | local |
| $\alpha$, $\pi$ | significance level; power | canonical |
| $B_i$, $D_{bc,it}$, $\beta_{bm:D}$ | biomarker, exposure-decayed drug indicator, interaction coefficient | canonical |
| $t_{1/2}$ | carryover half-life | canonical |
| $m$ | number of levels of the factor whose sphericity is at issue | local |

Table: Notation index and glossary of symbols used in this document

- **Strict RM-ANOVA.** Balanced split-plot analysis of variance with
  categorical factors, method-of-moments estimation and the
  sphericity assumption.
- **Sphericity.** All pairwise differences of the within-participant
  levels have equal variance. Compound symmetry implies it; the
  converse does not hold. Any factor with two levels is spherical.
- **E1, E2.** The binary on-drug indicator and the exposure-decayed
  continuous indicator $D_{bc,it}$ (NOTATION.md, Part 3; formerly
  labeled `A1` and `A2`).
- **Box's measure.** $m$ denotes the number of levels of the
  within-participant factor whose sphericity is at issue.

## 1. Motivation

The founding pmsimstats question is the test of the
biomarker-treatment interaction in an aggregated N-of-1 trial
[15]. The pmsimstats-ng analysis layer fits a continuous-time
mixed-effects model [8, 9], in R
`nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc)` with a random intercept
`~ 1 | ptID` and residual correlation `corCAR1(form = ~ t | ptID)`,
and tests the
interaction coefficient $\beta_{bm:D}$ via a Wald $t$-statistic.
A reviewer trained in
classical analysis-of-variance methodology may reasonably ask
whether a strict repeated-measures ANOVA (RM-ANOVA) test of the
same interaction is available, what its closed-form expression
looks like, and how it relates to the mixed-model test the
framework currently produces. This white paper supplies that
answer for the strict classical case.

The presentation is deliberately confined to **strict** RM-ANOVA:
balanced split-plot design, categorical predictors, sphericity
assumption, method-of-moments estimation. Modern extensions (the
mixed-effects model for repeated measures, MMRM [10];
generalized estimating equations, GEE [11]; ANCOVA with random
participant effects) are mentioned for context but are not strict
RM-ANOVA in the textbook sense.

## 2. Notation and design

Consider the simplest design that supports a biomarker-treatment
interaction test in a within-participant setting:

- Two biomarker strata $g \in \{1, 2\}$, formed by dichotomizing
  the continuous biomarker $B_i$ at a fixed cut-point.
  Strict RM-ANOVA requires categorical predictors; the
  dichotomization is the price of admission to the framework
  [12, 13, 14].
- $n$ participants per stratum, balanced, so $N = 2n$ in total.
  Index $i \in \{1, \ldots, n\}$ within stratum.
- Two treatment phases $j \in \{1, 2\}$ (active drug, placebo).
  In an N-of-1 design these are two phases through which every
  participant passes; in a parallel-group design they are
  between-participant and the design is not a split-plot.
- $T$ within-phase timepoints, indexed $t \in \{1, \ldots, T\}$.
  In the simplest case $T = 1$, the analyst pre-averages
  within-phase observations to a single mean per participant per
  phase.

The biomarker stratum is **between-participants** (each
participant sits in exactly one stratum); the treatment phase and
the timepoints are **within-participants** (each participant sees
both phases at every timepoint). This is the canonical
*split-plot* repeated-measures design [1].

## 3. The strict RM-ANOVA model

Write $Y_{gjit}$ for the symptom outcome at biomarker stratum $g$,
treatment phase $j$, participant $i$ within stratum $g$, and
timepoint $t$. The strict classical split-plot model is

$$
\begin{aligned}
Y_{gjit} = {} & \mu + \xi_g + \phi_j + (\xi\phi)_{gj}
+ u_{i(g)} \\
& + \tau_t + (\xi\tau)_{gt} + (\phi\tau)_{jt}
+ (\xi\phi\tau)_{gjt}
+ \varepsilon_{gjit} ,
\end{aligned} \qquad (1)
$$

with sum-to-zero constraints on each fixed effect, the random
participant effect $u_{i(g)} \sim N(0, \sigma_u^2)$ nested in
stratum, and within-participant error
$\varepsilon_{gjit} \sim N(0, \sigma_\varepsilon^2)$, independent.
Model (1) therefore implies compound symmetry of the
within-participant covariance: every pair of within-participant
observations has the same covariance $\sigma_u^2$, and every
observation has the same total variance
$\sigma_u^2 + \sigma_\varepsilon^2$. The classical $F$-tests need
only the weaker sphericity condition on each within-participant
factor (glossary), which compound symmetry implies.

The biomarker-treatment interaction is the
$(\xi\phi)_{gj}$ term. The time-resolved version of the
interaction is the three-way $(\xi\phi\tau)_{gjt}$ term.

## 4. Sum-of-squares decomposition

The total sum of squares partitions into between-participants and
within-participants components:

$$
SS_{\text{tot}} = SS_{\text{between}} + SS_{\text{within}} . \qquad (2)
$$

The between-participants component decomposes further into
$SS_A$ (biomarker stratum) and $SS_{S(A)}$ (participant within
stratum):

$$
SS_A = bnT \sum_{g=1}^a (\bar Y_{g\cdots} - \bar Y_{\cdots\cdot})^2 ,
\qquad (3)
$$

$$
SS_{S(A)} = bT \sum_{g=1}^a \sum_{i=1}^n
(\bar Y_{g\cdot i\cdot} - \bar Y_{g\cdots})^2 . \qquad (4)
$$

The within-participants component contains $SS_B$ (treatment
phase), $SS_{AB}$ (the **biomarker $\times$ treatment
interaction**, the target of the test), and the within-participants
error term for both, $SS_{B \times S(A)}$, together with the time
terms of Section 5d:

$$
SS_B = anT \sum_{j=1}^b
(\bar Y_{\cdot j\cdot\cdot} - \bar Y_{\cdots\cdot})^2 , \qquad (5)
$$

$$
SS_{AB} = nT \sum_{g,j}
(\bar Y_{gj\cdot\cdot} - \bar Y_{g\cdots} -
\bar Y_{\cdot j\cdot\cdot} + \bar Y_{\cdots\cdot})^2 , \qquad (6)
$$

$$
SS_{B \times S(A)} = T \sum_{g,j,i}
(\bar Y_{gji\cdot} - \bar Y_{g\cdot i\cdot} -
\bar Y_{gj\cdot\cdot} + \bar Y_{g\cdots})^2 . \qquad (7)
$$

The corresponding degrees of freedom are
$df_A = a - 1$, $df_{S(A)} = a(n-1)$, $df_B = b - 1$,
$df_{AB} = (a-1)(b-1)$, and
$df_{B \times S(A)} = a(n-1)(b-1)$ (the latter being the
within-participants error degrees of freedom appropriate for the
treatment effect and its interaction with the between-participants
factor).

## 5. The test statistic: closed form

The biomarker-treatment interaction is tested by the $F$-ratio

$$
F_{AB} = \frac{MS_{AB}}{MS_{B \times S(A)}}
= \frac{SS_{AB}/df_{AB}}{SS_{B \times S(A)}/df_{B \times S(A)}} .
\qquad (8)
$$

Under the null $H_0: (\xi\phi)_{gj} \equiv 0$ and the
sphericity assumption on the phase factor,
$F_{AB} \sim F(df_{AB}, df_{B \times S(A)})$.

### 5a. The $2 \times 2$ split-plot case

Specialize to $a = b = 2$ (two biomarker strata, two phases). The
formulas below keep $T$ general; $T = 1$ is the pre-averaged case.
Writing $\bar Y_{gj}$ for the cell mean $\bar Y_{gj\cdot\cdot}$,
the sample interaction contrast is

$$
\hat L \;=\; (\bar Y_{11} - \bar Y_{12}) -
            (\bar Y_{21} - \bar Y_{22}) , \qquad (9)
$$

i.e.\ the difference of within-participant treatment differences
across biomarker strata, an estimate of the population contrast
$L$. This is the natural biomarker-treatment interaction
statistic.

The interaction parameter under the sum-to-zero ANOVA
parameterization satisfies
$\widehat{(\xi\phi)}_{11} = \hat L / 4$ (and the other three cells
take values $\pm \hat L / 4$; derived), so the interaction sum of
squares is

$$
SS_{AB} \;=\; nT \cdot 4 \cdot (\hat L / 4)^2 \;=\; \frac{nT\,\hat L^2}{4} .
\qquad (10)
$$

With $df_{AB} = 1$, the interaction mean square equals the
interaction sum of squares.

The within-participants error mean square is the pooled estimator
of $\sigma_\varepsilon^2$ from the
treatment-by-participant-within-stratum component of variance.
Under model (1) and balanced data, $MS_{B \times S(A)}$ is an
unbiased estimator $\hat\sigma_\varepsilon^2$ of
$\sigma_\varepsilon^2$.

Substituting into equation (8):

$$
\boxed{\,F_{AB} \;=\; \frac{nT\,\hat L^2}{4\,\hat\sigma_\varepsilon^2}\,} ,
\qquad df_1 = 1, \quad df_2 = 2(n-1) . \qquad (11)
$$

This is the closed-form expression for the strict classical
RM-ANOVA test statistic for the biomarker-treatment interaction
in a $2 \times 2 \times T$ split-plot design.

### 5b. Within-participant contrast representation

Equation (11) is equivalent to a two-sample $t$-test on the
within-participant difference. Define
$\Delta_{gi} = \bar Y_{g1i\cdot} - \bar Y_{g2i\cdot}$ as the
within-participant treatment difference for participant $i$ in
biomarker stratum $g$, averaged over $T$ within-phase timepoints.
Then

$$
\hat L \;=\; \bar\Delta_1 - \bar\Delta_2 , \qquad
\bar\Delta_g \;=\; \frac{1}{n} \sum_{i=1}^n \Delta_{gi} . \qquad (12)
$$

Under model (1), $\Delta_{gi}$ are iid within stratum with
$\text{Var}(\Delta_{gi}) = 2\sigma_\varepsilon^2 / T$ (the random
participant effect cancels in the within-participant difference).
The variance of $\hat L$ is therefore

$$
\text{Var}(\hat L) \;=\; \frac{2 \cdot 2\sigma_\varepsilon^2}{nT}
\;=\; \frac{4\sigma_\varepsilon^2}{nT} . \qquad (13)
$$

The pooled-variance two-sample $t$-statistic for the contrast is

$$
t_{AB} \;=\; \frac{\hat L}{\sqrt{4\hat\sigma_\varepsilon^2/(nT)}}
\;=\; \sqrt{\frac{nT}{4}} \cdot \frac{\hat L}{\hat\sigma_\varepsilon} ,
\qquad (14)
$$

and $t_{AB}^2 = F_{AB}$, recovering equation (11) (derived). The
pooled variance of the $\Delta_{gi}$ equals
$2\hat\sigma_\varepsilon^2/T$, so this is the ordinary pooled $t$,
not a Welch statistic. The strict classical RM-ANOVA $F$-test for
the biomarker-treatment interaction is therefore identical to a
two-sample $t$-test on the within-participant treatment
differences, with degrees of freedom $2(n-1)$. Because that
$t$-test needs only that the $\Delta_{gi}$ be iid normal with a
common variance in the two strata, the identity, and the exactness
of the test, hold for any within-participant covariance, AR(1)
included, once each phase is reduced to its mean (derived). This
identity is well known in the split-plot literature [1, ch.\ 12;
2] and is the analytical anchor for the
statement that the strict-RM-ANOVA test of a biomarker-treatment
interaction reduces to a simple two-sample comparison of
treatment-effect estimates between biomarker strata.

### 5c. Generalization: $a$ strata, $b$ phases

For $a$ biomarker strata, $b$ treatment phases, $n$ participants per
stratum, and $T$ within-phase timepoints, equation (8) generalizes
to

$$
F_{AB} = \frac{SS_{AB}/[(a-1)(b-1)]}{SS_{B \times S(A)}/[a(n-1)(b-1)]} .
\qquad (15)
$$

Under sphericity of the phase factor (automatic when $b = 2$),
$F_{AB} \sim F((a-1)(b-1), a(n-1)(b-1))$ under $H_0$. The closed
form of the interaction sum of squares follows directly from
equation (6).

### 5d. The time-resolved three-way case

Including the time factor $\tau_t$ explicitly, the test of
whether the biomarker modifies the treatment effect across time is
the three-way interaction $(\xi\phi\tau)_{gjt}$. The
$F$-statistic is

$$
F_{AB\tau} = \frac{MS_{AB\tau}}{MS_{B\tau \times S(A)}} ,
\qquad (16)
$$

with $df_1 = (a-1)(b-1)(T-1)$ and $df_2 = a(n-1)(b-1)(T-1)$. The
sums of squares follow the standard split-plot pattern [1, ch.\ 14].

## 6. Sphericity, Mauchly, and the Greenhouse-Geisser correction

The $F$-tests of Section 5 presuppose sphericity of the
within-participant factor being tested: all pairwise differences
of its levels have equal variance. Compound symmetry is sufficient
but not necessary. A factor with two levels is always spherical,
so equations (11) and (15) with $b = 2$ are exact whatever the
covariance, as shown at the end of Section 5b; sphericity bites on
the time-resolved tests of Section 5d, whose within-participant
factor has $T > 2$ levels, and on any test with $b > 2$. When
sphericity fails, the nominal $F$-distribution is no longer exact
under $H_0$, and the uncorrected test inflates Type I error. One
diagnostic, one summary measure and two corrections are standard:

- **Mauchly's test** [3]. Likelihood-ratio test of sphericity
  against an unstructured alternative. Often rejects in moderate
  samples; not informative about the direction of departure.
- **Box's $\epsilon_{\mathrm{Box}}$** [4, 5]. A scalar measure of
  departure from sphericity, with $\epsilon_{\mathrm{Box}} = 1$
  under sphericity and $\epsilon_{\mathrm{Box}} \geq 1/(m-1)$ at
  maximum departure for a factor with $m$ levels ($m = T$ for the
  time factor).
- **Greenhouse-Geisser correction** [6]. Estimate
  $\hat\epsilon_{\mathrm{Box}}$ from the observed
  within-participant covariance and multiply both degrees of
  freedom of the test by it. Conservative.
- **Huynh-Feldt correction** [7]. Less conservative variant,
  recommended when $\hat\epsilon_{\mathrm{Box}} > 0.75$.

The pmsimstats DGP imposes AR(1) within-factor correlation, which
is **not** compound-symmetric. Under AR(1), Box's
$\epsilon_{\mathrm{Box}}$ for a factor with more than two levels
is strictly below 1, and the Greenhouse-Geisser correction
attenuates the degrees of freedom of the $F$-tests on that
factor. The project audit `docs/06-ar1-residual-correlation.tex`
records Type I error of 0.13 (OL) and 0.17 (CO) at $N = 35$ (200
replicates) under the AR(1) DGP for a random-intercept mixed model
with conditionally independent residuals (`lme4::lmer`), which
imposes the same compound-symmetric working covariance as model
(1) (inspected). That figure is indicative only: it is not a
simulation of the RM-ANOVA $F$-test, and the Greenhouse-Geisser
corrected Type I error and power were not simulated for this
note. The mixed-model analysis with `corCAR1` avoids the issue
because it models the AR(1) structure directly.

*Correction, 2026-10-07.* An earlier version of this section said
that equations (11) and (15) presuppose compound symmetry and
attributed the 13 to 17 percent figure to strict RM-ANOVA. Paper 08
(`analysis/report/08-test-procedure-design-sensitivity/`, Section
2) makes the same correction and reports that its pre-averaged
strict RM-ANOVA arm holds nominal Type I error (inspected).

## 7. Power calculation under the closed form

For the $2 \times 2$ case, equation (11) yields a clean
closed-form power formula. Under the alternative
$H_1: L = L_1$, the statistic
$F_{AB} = nT\,\hat L^2 / (4\hat\sigma_\varepsilon^2)$ follows a
non-central $F$ distribution with non-centrality parameter

$$
\Lambda \;=\; \frac{nT\,L_1^2}{4\sigma_\varepsilon^2} . \qquad (17)
$$

Power at level $\alpha$ is

$$
\pi \;=\; \Pr\!\left[ F_{AB} > F_{1, 2(n-1), 1-\alpha}
\;\Big|\; \Lambda \right] , \qquad (18)
$$

which is computable via standard non-central $F$ functions
(`stats::pf` in R with the `ncp` argument; `PROC POWER` in SAS).
For sample-size calculation, equation (17) inverts to

$$
n \;=\; \frac{4\sigma_\varepsilon^2 \Lambda}{T L_1^2} , \qquad (19)
$$

with $\Lambda$ chosen to deliver the target power at the design
$\alpha$. Because the critical value in (18) depends on $n$
through $df_2 = 2(n-1)$, the required $\Lambda$ itself depends on
$n$, and (19) is solved by iteration (derived).

This formula assumes (i) the variance $2\sigma_\varepsilon^2/T$
of each phase difference implied by model (1), (ii) categorical
biomarker at the dichotomized cut-point with equal stratum sizes,
(iii) no carryover, (iv) no within-treatment relapse, (v) iid
within-phase observations after within-participant differencing.
All five assumptions are violated to varying degrees in the
pmsimstats DGP; see Section 8. Under AR(1) errors the test stays
exact (Section 5b), but the variance of $\Delta_{gi}$ is no longer
$2\sigma_\varepsilon^2/T$, so (17) must use the actual variance of
the phase difference.

## 8. Limitations of the strict classical formulation

The closed-form simplicity of equation (11) is bought with five
strong assumptions, each of which is empirically problematic in
the pmsimstats setting.

1. **Categorical biomarker.** The pmsimstats biomarker (resting
   SBP in the prazosin-PTSD application) is continuous and has
   no defensible cut-point. Dichotomizing at the median or at a
   clinical threshold loses the within-stratum biomarker
   variation, attenuates the interaction estimator, and inflates
   Type II error [12, 13, 14]. The information loss can be
   substantial. An earlier version of this point quoted a median
   split standard error 1.4 to 1.7 times the continuous one and a
   30 to 50 percent power penalty "in the carryover-sensitivity
   production data"; no script or stored output producing those
   figures was found (2026-10-07), so they are withdrawn. Paper 08
   has since measured the cost directly: replacing the
   dichotomized biomarker by the continuous one raises power by
   0.08 to 0.21 across its sample sizes (0.180 at $N = 25$),
   while the further step to the full mixed model adds at most
   0.04 (paper 08 abstract, inspected).

2. **Sphericity.** Compound-symmetric within-participant
   covariance is empirically rejected for time-series symptom
   data: nearby timepoints are more correlated than distant ones.
   The AR(1) structure that pmsimstats implements violates
   sphericity for any within-participant factor with more than two
   levels. The pre-averaged two-phase test of equation (11) is
   unaffected (Section 6); the time-resolved test of Section 5d
   needs a Greenhouse-Geisser or Huynh-Feldt correction, which
   restores approximately nominal Type I error at a loss of
   degrees of freedom and power.

3. **No carryover model.** Strict RM-ANOVA has no provision for
   residual drug effect during washout; the model treats off-drug
   observations as if no drug had been administered. This is the
   binary-indicator specification E1 of the carryover-sensitivity
   manuscript (paper 02). In its stored grid
   (`analysis/data/02-grid-summary.rds`, $N = 70$, 500 replicates
   per cell, exponential decay), E1 at $t_{1/2} = 1.0$ week
   attenuates $\hat\beta_{bm:D}$ toward zero by 30 to 34 percent
   in the Hybrid design and by 40 to 42 percent in OL+BDC under
   covariance moderation, but by at most 4 percent in CO and by
   at most 1.4 percent in any design under mean moderation
   (verified, 2026-10-07). The attenuation of $\hat L$ is
   expected to follow the same pattern (argued, not simulated).

4. **Balanced design required.** Equations (11) and (15) assume
   equal stratum sizes and equal numbers of within-phase
   timepoints. Real trials have missing data; classical RM-ANOVA
   handles missingness only through complete-case deletion. The
   MMRM framework's likelihood-based handling of MAR data [10] is
   better behaved.

5. **Equal time spacing.** Strict RM-ANOVA treats time as a
   categorical factor with implicitly equal-spaced levels. Real
   N-of-1 trials have unequally-spaced visits; pmsimstats
   implements continuous-time AR(1) (`corCAR1`), which the
   strict-RM-ANOVA framework cannot do.

## 9. Relationship to the pmsimstats MMRM analysis

The strict-RM-ANOVA $F$-statistic in equation (11) and the
pmsimstats `nlme::lme` Wald $t$-statistic on `bm:Dbc` are
**asymptotically equivalent** under three conditions: (i) the
biomarker enters the mixed model as the stratum indicator, (ii)
the within-participant covariance is exactly compound-symmetric,
and (iii) the drug indicator is binary (specification E1 rather
than the exposure-decayed $D_{bc,it}$ of E2). When all three hold,
the two tests reach the same coefficient and the same standard
error up to small-sample correction factors (argued). When any
fails, the mixed-model test extracts more information than the
strict RM-ANOVA: the continuous biomarker contributes
within-stratum variation, the AR(1) covariance is modeled
directly, and the continuous $D_{bc,it}$ exposure-decay predictor
(E2) uses the off-drug trajectory.

For the founding pmsimstats interaction question [15], the
strict-RM-ANOVA test is therefore a **lower-power reduction** of
the available analysis rather than a competing alternative. Where
a reviewer requests 'RM-ANOVA results', the natural response is
one of (a) point at the existing mixed-model analysis as a
continuous-time member of the repeated-measures family, (b) report
the strict-RM-ANOVA test from equation (11) as a sensitivity
comparator under dichotomized biomarker and binary phase, or (c)
report both with the explicit power-loss tabulation.

## 10. Recommendation for pmsimstats

We recommend retaining the MMRM-with-continuous-time analysis as
primary and adding a strict-RM-ANOVA comparator as a planned
sensitivity:

1. **Categorize the biomarker at the median.** Conventional
   choice; the alternative of a fixed clinical cut-point is
   acceptable when one exists.
2. **Pre-average within-phase timepoints to a single mean per
   participant per phase.** Reduces the $T$-axis to $T = 1$ for
   the sensitivity test and makes the test exact under any
   within-participant covariance (Section 5b).
3. **Compute $\hat L$ per equation (9) and $F_{AB}$ per equation
   (11).** This is computable by hand from the four cell means
   and the within-participant error mean square.
4. **Report the strict-RM-ANOVA $p$-value alongside the mixed-model
   $p$-value in the manuscript table.** The two should be in the
   same order of magnitude. Substantial divergence between them
   indicates that the dichotomization, the pre-averaging, or the
   binary phase indicator is doing meaningful work, and
   the discrepancy is informative for interpretation.

The carryover-sensitivity manuscript at
`analysis/report/02-carryover-sensitivity/` already evaluates the
E1 binary-phase specification, which is the analysis-model
counterpart to the strict-RM-ANOVA closed form derived here. The
manuscript reports the E1 specification's bias and coverage under
carryover. The strict-RM-ANOVA test is expected to inherit the E1
specification's attenuation properties (argued; see Section 8,
point 3). Paper 08 implements this comparator as its `rmanova`
arm (`analysis/scripts/test-procedure-design-sensitivity/
01-study1-test-procedure.R`).

## References

Metadata were checked on 2026-10-07 against PubMed (entries 8, 12 to
15) and Crossref (entries 4 to 7, 10); entries 1, 2 and 9 are books
and entries 3 and 11 are cited from their standard records. An
uncited Gompertz (1825) entry was removed.

1. Maxwell SE, Delaney HD, Kelley K. *Designing Experiments and
   Analyzing Data: A Model Comparison Perspective*. 3rd ed.
   Routledge; 2017.
2. Kirk RE. *Experimental Design: Procedures for the Behavioral
   Sciences*. 4th ed. SAGE; 2013.
3. Mauchly JW. Significance test for sphericity of a normal
   $n$-variate distribution. *Annals of Mathematical Statistics*.
   1940;11(2):204-209.
4. Box GEP. Some theorems on quadratic forms applied in the study
   of analysis of variance problems, I. Effect of inequality of
   variance in the one-way classification. *Annals of Mathematical
   Statistics*. 1954;25(2):290-302. doi:10.1214/aoms/1177728786.
5. Box GEP. Some theorems on quadratic forms applied in the study
   of analysis of variance problems, II. Effects of inequality of
   variance and of correlation between errors in the two-way
   classification. *Annals of Mathematical Statistics*.
   1954;25(3):484-498. doi:10.1214/aoms/1177728717.
6. Greenhouse SW, Geisser S. On methods in the analysis of profile
   data. *Psychometrika*. 1959;24(2):95-112. doi:10.1007/BF02289823.
7. Huynh H, Feldt LS. Estimation of the Box correction for degrees
   of freedom from sample data in randomized block and split-plot
   designs. *Journal of Educational Statistics*. 1976;1(1):69-82.
   doi:10.3102/10769986001001069.
8. Laird NM, Ware JH. Random-effects models for longitudinal data.
   *Biometrics*. 1982;38(4):963-974.
9. Verbeke G, Molenberghs G. *Linear Mixed Models for Longitudinal
   Data*. Springer; 2000.
10. Mallinckrodt CH, Lane PW, Schnell D, Peng Y, Mancuso JP.
    Recommendations for the primary analysis of continuous
    endpoints in longitudinal clinical trials. *Drug Information
    Journal*. 2008;42(4):303-319. doi:10.1177/009286150804200402.
11. Liang K-Y, Zeger SL. Longitudinal data analysis using
    generalized linear models. *Biometrika*. 1986;73(1):13-22.
12. Royston P, Altman DG, Sauerbrei W. Dichotomizing continuous
    predictors in multiple regression: a bad idea. *Statistics in
    Medicine*. 2006;25(1):127-141. doi:10.1002/sim.2331.
13. Altman DG, Royston P. The cost of dichotomising continuous
    variables. *BMJ*. 2006;332(7549):1080.
    doi:10.1136/bmj.332.7549.1080.
14. Senn S. Mastering variation: variance components and
    personalised medicine. *Statistics in Medicine*.
    2016;35(7):966-977. doi:10.1002/sim.6739.
15. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
    aggregated N-of-1 trial designs for predictive biomarker
    validation: statistical methods and theoretical findings.
    *Frontiers in Digital Health*. 2020;2:13.
    doi:10.3389/fdgth.2020.00013.

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
