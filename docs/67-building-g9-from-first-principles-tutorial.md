---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# From a Paired Difference to G9: Building the Analysis of an N-of-1 Biomarker Trial from First Principles {.unlisted .unnumbered}
*2026-10-30 09:20 PDT*

**Author.** pmsimstats team

**Audience.** Undergraduate statistics students who have taken a
course in linear regression. You should know what a regression
coefficient, a standard error, a $p$-value and a $t$-test are. You do
not need to know mixed models, sandwich estimators or model selection.
This tutorial introduces each of them at the point where the problem
requires it.

**Purpose.** The pmsimstats project recommends an analysis called G9
for one kind of clinical trial (paper 02,
`analysis/report/02-carryover-sensitivity/`). G9 is a linear
mixed-effects model with a random intercept, a continuous-time
autoregressive residual correlation, an exposure-weighted drug
indicator whose half-life is chosen by AIC, and cluster-robust CR2
standard errors with Satterthwaite degrees of freedom. Stated in one
sentence, that model looks arbitrary. It is not. Each component
answers a specific question that the data raise. This tutorial starts
with the simplest analysis a student could propose and adds one
component at a time. At each step it states the problem, the remedy,
and the evidence that the remedy works.

**How to read it.** A notation table and a glossary come first.
Consult them whenever a symbol or term is unfamiliar. Sections 1 to 3
set up the scientific question
and the tools. Sections 4 to 10 are the ladder, one rung per section.
Section 11 assembles G9 in R. Section 12 lists what G9 does not
solve. Exercises close the document. A companion
script, `analysis/scripts/tutorial-g9/01-toy-ladder.R`, runs a small
simulation that reproduces the direction of each step. Its numbers
appear in Section 3.4 and are referred to throughout.

```{=latex}
\clearpage
\tableofcontents
\clearpage
```

## Notation

Symbols follow the project's canonical notation,
`analysis/report/NOTATION.md`. The section column gives where each
symbol first appears.

| Symbol | Meaning | Code | Section |
|------------|----------------------------------------------------------------|----------------|---------|
| $N$ | total number of participants, across all paths | | 1.3 |
| $i$ | participant, $i = 1, \ldots, N$ | `ptID` | 1.3 |
| $t$ | visit, $t = 1, \ldots, 8$ | | 1.3 |
| $w_t$ | calendar week of visit $t$ (4, 8, 9, 10, 11, 12, 16, 20) | `t` | 1.3 |
| $Y_{it}$ | symptom score of participant $i$ at visit $t$; lower is better | `Sx` | 2 |
| $B_i$ | participant's biomarker, measured once before treatment | `bm` | 1.1 |
| $\bar B$ | sample mean of the biomarker; $B_i - \bar B$ is the centered biomarker | `bmc` | 3.1 |
| $b_i$ | standardized biomarker, mean 0 and SD 1 (toy DGP) | `bm_z` | 3.4 |
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug | `on` | 1.3 |
| $D_{bc,it}$ | continuous exposure-decayed drug indicator | `Dbc` | 8.2 |
| $D^{*}_{it}$ | true drug exposure in the toy DGP, $D_{bc,it}$ at the true half-life | | 3.4 |
| $t_{sd}$ | time since the drug was stopped, in weeks | `tsd` | 1.3 |
| $t_{1/2}$ | true carryover half-life | | 8.1 |
| $t_{1/2}^a$, $\hat t_{1/2}^a$ | half-life assumed by the analysis; its AIC choice | `h_hat` | 8.2, 9.2 |
| $\mathcal T$ | candidate half-life grid, $\{0.25, 0.5, 1, 2\}$ weeks | `grid` | 9.2 |
| $\eta_t$ | design expectancy: 1 open label, 0.5 blinded | `De` | 12 |
| $\Delta_i$ | participant's on-drug mean minus off-drug mean | | 2 |
| $\gamma_0$, $\gamma_1$ | intercept and slope of $\Delta_i$ on $B_i$ | | 2 |
| $\beta_0$ | intercept of the visit-level model | `(Intercept)` | 3.1 |
| $\beta_B$ | biomarker main effect (prognostic slope, off drug) | `bm` | 3.1 |
| $\beta_T$ | time trend, points per week | `t` | 3.1 |
| $\beta_D$ | drug main effect at the mean biomarker | `D`, `Dbc` | 3.1 |
| $\beta_{bm:D}$ | biomarker-by-treatment interaction, the target | `bm:D`, `bm:Dbc` | 3.1 |
| $\varepsilon_{it}$, $e_{it}$ | residual error; $e_{it}$ once a random intercept is split off | | 3.1, 5.1 |
| $u_i$ | participant random intercept, $\mathcal N(0, \sigma_u^2)$ | `random` | 5.1 |
| $\sigma_u^2$, $\sigma_e^2$ | random-intercept and residual variances | | 5.1 |
| $\rho$ | AR(1) correlation per week, $\rho^{\lvert w_t - w_s \rvert}$ | `corCAR1` | 6.2 |
| $X_i$, $V_i$, $r_i$ | participant $i$'s design matrix, working covariance and residual vector | | 7.2 |
| $\kappa$ | calibration ratio: mean model SE over the empirical SD of the estimates | | 7.1 |
| $L$, $k$ | maximized likelihood and number of parameters in AIC | | 9.2 |
| $\alpha$ | significance level, 0.05 throughout | | 3.2 |
| $\hat\pi$ | estimated rejection rate (size or power) | | 3.3 |
| $R$ | number of simulation replicates | `n_reps` | 3.3 |
| MCSE | Monte Carlo standard error, $\sqrt{\hat\pi(1-\hat\pi)/R}$ | | 3.3 |

Table: Notation

The symptom score falls as patients improve. Treatment effects are
therefore negative numbers. A negative interaction means that
participants with a higher biomarker improve more on drug.

## Glossary

Terms are listed alphabetically, with analysis codes grouped at the
end. The number in brackets is the section where the term is
introduced.

- **Aggregated N-of-1 trial.** Many single-patient on/off trials
  analyzed together as one data set [1.2].
- **AIC.** Akaike information criterion, $-2\log L + 2k$. Among
  candidate models fitted to the same data, the smallest AIC is
  preferred [9.2].
- **AR(1), corCAR1.** A residual correlation $\rho^{\lvert w_t - w_s
  \rvert}$ that fades with the time between visits. `corCAR1` is the
  continuous-time version in the R package `nlme`, which handles
  unequally spaced visits [6.2].
- **Attenuation.** Bias of a coefficient toward zero caused by
  measurement error in its predictor [8.1].
- **Blinded discontinuation.** A phase in which the patient does not
  know whether the drug has been switched to placebo [1.3].
- **Calibration ratio $\kappa$.** Mean model-based standard error
  divided by the empirical standard deviation of the estimates across
  simulated trials. Above 1 means a conservative test [7.1].
- **Carryover.** Drug effect that persists after the drug is stopped
  [8.1].
- **Centering.** Subtracting the sample mean from a predictor, so that
  the other main effect in an interaction is read at the average
  value [3.1].
- **Cluster.** A group of observations that are correlated with one
  another; here the eight visits of one patient [4].
- **Cluster-robust (sandwich) variance.** A variance estimate for
  regression coefficients built from per-cluster residuals rather than
  from the assumed covariance [7.2].
- **Compound symmetry.** Equal correlation between every pair of
  visits, whatever the time between them [5.2].
- **Conservative, liberal.** A test whose size is below, or above,
  its nominal level [3.2].
- **CR2.** The Bell-McCaffrey small-sample correction of the
  cluster-robust variance [7.3].
- **Crossover.** A design in which each patient receives drug and
  placebo in randomized order [1.3].
- **Data-generating process (DGP).** The set of rules a simulation
  uses to produce artificial data [3.3].
- **Exposure-weighted predictor.** The drug indicator $D_{bc,it}$,
  equal to 1 on drug and decaying by half every $t_{1/2}^a$ weeks
  after the drug stops [8.2].
- **Extended DGP.** The project's main simulation model, built on that
  of Hendrickson et al. (2020), with drug-response, placebo-response
  and natural-history components [3.3].
- **Fixed effects.** In a mixed model, the population coefficients
  $\beta$. In econometrics, one indicator variable per unit. The two
  uses differ [5.1].
- **Half-life.** The time for the remaining drug effect to fall by
  half [8.1].
- **Hierarchy principle.** Keep both main effects in a model that
  contains their interaction [3.1].
- **Hybrid design.** Open-label titration, blinded discontinuation,
  then a brief crossover [1.3].
- **Interaction.** A coefficient on a product of two predictors. Here,
  the change in drug effect per unit of biomarker [3.1].
- **Intraclass correlation.** $\sigma_u^2 / (\sigma_u^2 + \sigma_e^2)$,
  the correlation between two visits of one patient implied by a
  random intercept alone [5.2].
- **Main effect.** The coefficient of a single predictor. In a model
  with an interaction, it is the effect when the other predictor is
  zero, or at its mean after centering [3.1].
- **MCSE.** Monte Carlo standard error of a simulated rejection rate
  [3.3].
- **Mean-shift and correlation encodings.** Two ways a simulation can
  build in the interaction: by shifting the mean drug effect with the
  biomarker (the toy DGP), or by a correlation between the biomarker
  and the drug response (the extended DGP) [3.4].
- **Mixed-effects model.** A regression with both population
  coefficients and patient-specific random effects [5.1].
- **ML, REML.** Maximum likelihood and restricted maximum likelihood.
  REML gives less biased variance components. Only ML likelihoods can
  be compared across different fixed-effects designs [5.3, 9.2].
- **Model-based standard error.** A standard error computed from the
  fitted model's assumed covariance [7.1].
- **N-of-1 trial.** A trial in one patient, who is observed on and off
  drug several times [1.2].
- **Open label.** A phase in which the patient knows they are taking
  the drug [1.3].
- **Paired difference.** The difference between two analyses applied
  to the same simulated trials. Its standard error is smaller than
  that of either rate alone [3.4].
- **Path.** One sequence of drug states through the trial. The Hybrid
  design has four [1.3].
- **Post-selection inference.** Inference after the same data chose
  the model. Nominal $p$-values can then be too small [9.3].
- **Power.** The probability of rejecting the null hypothesis when it
  is false [3.2].
- **Pre-exposure rule.** $D_{bc,it} = 0$ at off-drug visits before the
  first dose [8.2].
- **Predictive, prognostic.** A biomarker that modifies the treatment
  effect, and one that predicts the outcome regardless of treatment
  [1.1].
- **Random effects.** Patient-specific deviations modeled as random
  draws from a distribution [5.1].
- **Satterthwaite degrees of freedom.** Estimated degrees of freedom
  for the $t$ reference distribution of a robust test [7.3].
- **Size (Type I error).** The probability of rejecting the null
  hypothesis when it is true [3.2].
- **Working covariance.** The covariance the analysis model assumes,
  which may differ from the true one [7.1].

**Analysis codes.** These are the project's labels. Paper 02, Table
`tab-spec-summary`, is the full definition.

- **E0.** The original analysis of Hendrickson et al. (2020): random
  intercept, binary drug indicator, no residual correlation, `lme4`
  with Satterthwaite $p$-values, and an expectancy term [5.4].
- **G1.** Random intercept, `corCAR1` residuals, binary drug
  indicator, model-based inference [6.2].
- **G3.** G1 with the exposure-weighted predictor at a given
  half-life [8.2].
- **G4.** G3 with the half-life chosen by ML AIC over $\mathcal T$,
  then refitted by REML [9.2].
- **G5 (stored as E9).** The paired-difference statistic: each
  patient's on-minus-off mean regressed on the biomarker [2].
- **G6, G8, G9.** G1, G3 and G4 with CR2 standard errors and
  Satterthwaite degrees of freedom [7.3, 8.2, 9.2].
- **F1, F1A0, F1A0B.** Analyses that add the design expectancy
  $\eta_t$ and its interaction with the biomarker. F1 fixes the
  half-life. F1A0 chooses it on a model without the interaction.
  F1A0B calibrates F1A0 by a parametric bootstrap [9.3, 12].

## 1. The scientific question

### 1.1 Who benefits?

A drug rarely helps every patient equally. Suppose a blood test, the
biomarker $B_i$, could tell in advance who will benefit. A doctor
could then prescribe the drug to those patients and spare the others
its side effects. Before anyone relies on such a test, a trial must
show that the drug effect really changes with the biomarker. In
regression language the trial must show a **biomarker-by-treatment
interaction**. This is the estimand of the whole project.

A biomarker that predicts the size of the drug effect is called
**predictive**. A biomarker that predicts how sick the patient is,
regardless of treatment, is called **prognostic**. A biomarker can be
both. Keep the distinction in mind. It returns in Section 12.

### 1.2 Why an N-of-1 design?

In an ordinary parallel-group trial each patient receives either the
drug or a placebo. We never see the same patient both ways. The drug
effect for one person is not observed. We can only compare groups.

In an **N-of-1 trial** a single patient is observed on drug and off
drug several times. The patient's own on-minus-off difference
estimates that patient's drug effect. An **aggregated N-of-1 trial**
runs this for many patients and analyzes them together. Because each
patient serves as their own control, the person-to-person
variation in baseline severity cancels out of the comparison. The
interaction question becomes concrete: do the patients' own drug
effects line up with their biomarker values?

### 1.3 The Hybrid design

The project studies the **Hybrid design** of Hendrickson et al.
(2020): eight weeks of open-label drug, then a blinded discontinuation
block, then a short blinded crossover. There are eight visits, at
weeks 4, 8, 9, 10, 11, 12, 16 and 20. Patients are allocated to four
**paths**, which differ in when the drug stops and in the order of
the crossover (Table 1, from `docs/40`).

| Path | Drug state at weeks 4, 8, 9, 10, 11, 12, 16, 20 | Weeks since stopping |
|---|---|---|
| 1 | 1 1 1 1 0 0 1 0 | 0 0 0 0 1 2 0 4 |
| 2 | 1 1 1 1 0 0 0 1 | 0 0 0 0 1 2 6 0 |
| 3 | 1 1 1 0 0 0 1 0 | 0 0 0 1 2 3 0 4 |
| 4 | 1 1 1 0 0 0 0 1 | 0 0 0 1 2 3 7 0 |

Table: The four Hybrid paths. A 1 means on drug.

Three features of this table drive the rest of the tutorial.

1. Each patient contributes eight correlated measurements, not eight
   independent ones (Sections 5 to 7).
2. The visits are unequally spaced: four weeks apart early and late,
   one week apart in the middle (Section 6).
3. Several off-drug visits come one or two weeks after the drug
   stopped. The drug may still be acting (Section 8).

## 2. The first analysis: one number per patient

The simplest analysis follows directly from Section 1.2. For each
patient compute

$$
\Delta_i = \overline{Y}_{i,\text{on}} - \overline{Y}_{i,\text{off}},
$$

the mean score on drug minus the mean score off drug. Then regress
$\Delta_i$ on $B_i$ by ordinary least squares:

$$
\Delta_i = \gamma_0 + \gamma_1 B_i + \varepsilon_i .
$$

The slope $\gamma_1$ is the interaction. If the drug effect does not
depend on the biomarker, $\gamma_1 = 0$. The test of $\gamma_1 = 0$
is a plain regression $t$-test with $N - 2$ degrees of freedom.

This is a good analysis. Paper 02 calls it the paired-difference
statistic (G5, stored as E9) and it performs respectably there. It
has three weaknesses.

- It treats every off-drug visit alike. A visit one week after
  stopping is averaged with a visit seven weeks after stopping, even
  though the first may still contain drug effect.
- It discards the time order. Symptoms drift over twenty weeks, and
  in the Hybrid design every patient is on drug early and off drug
  late. Drift is then confused with the drug effect.
- It gives every patient equal weight, regardless of how many on and
  off visits that patient has.

The remaining sections build an analysis that keeps every visit as
its own row and handles these problems inside one model.

## 3. Tools for the rest of the tutorial

### 3.1 Main effects and interactions in one regression

Write the score at every visit as

$$
Y_{it} = \beta_0 + \beta_B B_i + \beta_T w_t + \beta_D D_{it}
       + \beta_{bm:D}\, B_i D_{it} + \varepsilon_{it}. \tag{1}
$$

Read it one term at a time.

- $\beta_B$ is the **biomarker main effect**. It is the difference in
  score per unit of biomarker *when off drug*. It describes a
  prognostic biomarker.
- $\beta_T$ is the time trend: the change in score per week that
  happens regardless of treatment.
- $\beta_D$ is the **drug main effect**. It is the drug effect for a
  patient whose biomarker equals zero.
- $\beta_{bm:D}$ is the **interaction**. The drug effect for a
  patient with biomarker $B_i$ is $\beta_D + \beta_{bm:D} B_i$. The
  interaction is the slope of the drug effect on the biomarker. It
  plays the role of $\gamma_1$ in Section 2.

Two habits follow from this reading.

**Keep the main effects when you fit an interaction.** If $\beta_B$
is dropped, the product $B_i D_{it}$ must absorb any prognostic
association between $B_i$ and $Y$. The interaction estimate is then
biased. The rule is sometimes called the hierarchy principle.

**Center the biomarker.** "The drug effect when $B_i = 0$" is
meaningless if no patient has a biomarker near zero. Replace $B_i$ by
$B_i - \bar B$. Then $\beta_D$ is the drug effect for an average
patient. The interaction coefficient and its test do not change. The
project's code calls the centered biomarker `bmc`.

In R, model (1) is written `Sx ~ bm + t + D + bm:D`, or equivalently
`Sx ~ bm * D + t`.

### 3.2 Size and power

Every analysis in this tutorial ends in a test of
$H_0: \beta_{bm:D} = 0$ at level $\alpha = 0.05$. Two numbers judge a
test.

- **Size**, or Type I error, is the probability of rejecting $H_0$
  when it is true. A correct test has size 0.05. A test with size
  above 0.05 is **liberal**: it finds biomarkers that do not work. A
  test with size below 0.05 is **conservative**: it is safe but wastes
  power.
- **Power** is the probability of rejecting $H_0$ when it is false.
  Higher is better, *but only among tests with correct size*. A test
  that rejects at random half the time has power 0.5 and is useless.

### 3.3 Simulation as a laboratory

Real trials cannot tell us the size of a test, because we never know
whether $H_0$ is true. A **simulation study** can. We write a
data-generating process (DGP), generate many artificial trials in
which we choose the truth, apply the analysis to each, and count.

- With the interaction set to zero, the fraction of trials in which
  the test rejects estimates the size.
- With the interaction set to a nonzero value, the same fraction
  estimates the power.

A fraction estimated from $R$ replicates has a **Monte Carlo standard
error**

$$
\mathrm{MCSE}(\hat\pi) = \sqrt{\hat\pi (1 - \hat\pi) / R}.
$$

With $R = 2000$ and a true size of 0.05, the MCSE is 0.005. An
estimate between about 0.04 and 0.06 is then consistent with a
correct test.

The project's main DGP, called the extended DGP in the manuscripts,
generates each patient's score as the sum of a drug-response
component, a placebo-response component and a natural-history
component. It builds a large correlation matrix across components and
visits. The interaction enters as a correlation between the
biomarker and the drug response that is present on drug and fades
after the drug stops. That DGP is described in paper 01 and is not
needed here.

### 3.4 A toy laboratory for this tutorial

To let you run every step yourself, the companion script uses a much
simpler DGP:

$$
Y_{it} = 50 + u_i - 0.2\, w_t + (-5 + \beta_{bm:D} b_i)\, D^{*}_{it}
       + e_{it},
$$

with $u_i \sim \mathcal N(0, 4^2)$, a standard normal biomarker
$b_i$, residuals $e_{it}$ with standard deviation 3 and correlation
$0.7^{|w_t - w_s|}$ between visits, and $D^{*}_{it}$ the true drug
exposure with a carryover half-life of one week (Section 8). There
are $N = 40$ patients, ten per path. The null cell sets
$\beta_{bm:D} = 0$. The alternative sets $\beta_{bm:D} = -0.8$. Each
cell has 2,000 replicates, and both cells use the same seeds. Here
the interaction shifts the mean of the drug effect, which is simpler
than the project's correlation encoding. Given the biomarker, the two
encodings behave alike, so the toy shows the same qualitative
behavior.

Table 2 gives the results. Each later section refers to one row.

| Rung | Analysis | Size (MCSE) | Power (MCSE) |
|---|---|---|---|
| A | OLS, every visit independent | 0.002 (0.001) | 0.095 (0.007) |
| B | random intercept, binary indicator | 0.082 (0.006) | 0.501 (0.011) |
| C | G1: random intercept + AR(1) | 0.050 (0.005) | 0.426 (0.011) |
| D | G6: G1 + CR2 | 0.058 (0.005) | 0.400 (0.011) |
| E | G3: exposure-weighted, true half-life | 0.056 (0.005) | 0.508 (0.011) |
| E | G8: G3 + CR2 | 0.056 (0.005) | 0.448 (0.011) |
| F | G9: AIC half-life + CR2 | 0.056 (0.005) | 0.450 (0.011) |

Table: Toy ladder, $N = 40$, 2,000 replicates per cell, the same
seeds in both cells. Size is the rejection rate at
$\beta_{bm:D} = 0$, power at $\beta_{bm:D} = -0.8$. Script
`analysis/scripts/tutorial-g9/01-toy-ladder.R`, output
`analysis/data/tutorial-g9/`.

Because every analysis was applied to the same simulated trials, two
rows can be compared trial by trial. The paired standard error of
such a difference is much smaller than either MCSE alone. The
sections below quote paired differences where they matter.

The rest of this tutorial explains why the rows change as they do.

## 4. Rung A: ordinary least squares, and why it fails

The obvious first fit of model (1) is ordinary least squares (OLS),
with every row of the long data set treated as an independent
observation. A trial with 40 patients and eight visits has 320 rows.
OLS believes it has 320 independent pieces of information.

It does not. A patient who is sicker than average at week 4 is
usually sicker than average at week 20 as well. The eight rows of one
patient are positively correlated. OLS still gives an unbiased
estimate of $\beta_{bm:D}$, but its standard error is computed from
the wrong variance formula. The test statistic is then not $t$
distributed, and the reported size is not 0.05.

The direction of the error depends on which comparison the
coefficient uses. For a between-patient coefficient such as $\beta_B$,
OLS overstates the information and is liberal. For a within-patient
contrast the patient-level variation cancels, and the OLS standard
error, which still includes it, can be too large. The interaction
$\beta_{bm:D}$ is mostly a within-patient contrast. In the toy, OLS
rejected a true null in 0.002 of trials and had power 0.095 (Table
2). Its standard error was far too large, so the test was almost
never significant. Either way, the fault is the same: the model's
variance assumption is wrong.

**Lesson.** Repeated measurements on the same unit are clustered. The
analysis must model, or at least respect, the clustering.

## 5. Rung B: fixed and random effects

### 5.1 Two ways to give every patient an intercept

The cheapest correction lets every patient have their own baseline
level:

$$
Y_{it} = \beta_0 + u_i + \beta_B B_i + \beta_T w_t + \beta_D D_{it}
       + \beta_{bm:D} B_i D_{it} + e_{it}. \tag{2}
$$

There are two ways to treat the $u_i$.

**As fixed effects.** Add one indicator column per patient and
estimate every $u_i$ as an ordinary coefficient. This is the
"fixed-effects estimator" of econometrics. It makes no assumption
about how the $u_i$ are distributed. It has a cost that matters here.
The biomarker $B_i$ is constant within a patient, so the column for
$B_i$ is a linear combination of the patient indicators. The main
effect $\beta_B$ cannot be estimated. The interaction
$B_i D_{it}$ changes within a patient, so it remains estimable.

**As random effects.** Treat the $u_i$ as random draws,
$u_i \sim \mathcal N(0, \sigma_u^2)$, independent of everything else
in the model. Only one extra parameter, $\sigma_u^2$, is estimated.
The model can still estimate $\beta_B$, because it borrows the
between-patient comparison as well as the within-patient one. The
model requires an assumption: $u_i$ must be unrelated to $B_i$. If the
biomarker is prognostic, $\beta_B$ in the model absorbs the
association, so this is usually acceptable.

A model with both kinds of terms is a **linear mixed-effects model**.
Be careful with the vocabulary. In a mixed model, "fixed effects"
means the population coefficients $\beta$, the things we want to
report, and "random effects" means the patient-specific deviations.
In econometrics "fixed effects" means patient indicators. The two
uses are different.

### 5.2 What a random intercept implies about correlation

Model (2) with independent residuals $e_{it}$ of variance $\sigma_e^2$
implies that any two visits of the same patient have correlation

$$
\mathrm{Cor}(Y_{it}, Y_{is}) = \frac{\sigma_u^2}{\sigma_u^2 + \sigma_e^2},
$$

the **intraclass correlation**. The value is the same for visits one
week apart and for visits sixteen weeks apart. This pattern is called
**compound symmetry**.

### 5.3 Estimation by REML

Mixed models are fitted by maximum likelihood (ML) or by restricted
maximum likelihood (REML). ML estimates of variance components are
biased downward in small samples, for the same reason that dividing a
sum of squares by $n$ instead of $n-1$ underestimates a variance.
REML corrects this by maximizing the likelihood of residual contrasts
that do not depend on the fixed effects. Use REML for the final fit.
Section 9 shows when ML is needed instead.

### 5.4 The original analysis

The published analysis of Hendrickson et al. (2020), which the
project calls **E0**, is essentially model (2) with a binary drug
indicator, fitted with `lme4::lmer`. In the toy, rung B
rejected a true null in 0.082 of trials (MCSE 0.006). A random
intercept alone gave a liberal test. Its power of 0.501 cannot be
compared with the other rows, because a liberal test rejects too
often under any hypothesis. In paper 02's simulations under the extended DGP, E0
rejected a true null in 0.110 of trials, more than twice the nominal
rate. A random intercept was not enough.

## 6. Rung C: correlation that fades with time (G1)

### 6.1 The problem

Compound symmetry says that weeks 9 and 10 are as alike as weeks 4
and 20. For symptom scores that is implausible. Visits close in time
share more of the patient's current state than visits far apart. In
the Hybrid design the middle visits are only a week apart, so this
matters.

### 6.2 The AR(1) remedy

A first-order autoregressive correlation, AR(1), says that the
correlation between two visits falls geometrically with the time
between them:

$$
\mathrm{Cor}(e_{it}, e_{is}) = \rho^{|w_t - w_s|}.
$$

With $\rho = 0.7$ per week, visits one week apart correlate at 0.7,
four weeks apart at 0.24, and sixteen weeks apart at 0.003. Because
the Hybrid visits are unequally spaced, the exponent must be the
calendar gap in weeks, not the count of visits. The continuous-time
version, `corCAR1(form = ~ t | ptID)` in the R package `nlme`, does
exactly this.

**G1** keeps the random intercept and adds AR(1) residuals:

$$
Y_{it} = \beta_0 + u_i + \beta_B B_i + \beta_T w_t + \beta_D D_{it}
       + \beta_{bm:D} B_i D_{it} + e_{it}, \quad
\mathrm{Cor}(e_{it}, e_{is}) = \rho^{|w_t - w_s|}. \tag{3}
$$

The total correlation between two visits is now a mixture: a constant
part from $u_i$ and a fading part from $e_{it}$.

### 6.3 What happened

Paper 02 compared E0 and G1 on the same simulated trials. E0 rejected
at 0.110. G1 rejected at 0.038. E0 and G1 differ in three ways: the
residual correlation, an extra expectancy term, and the software used
to compute $p$-values. A separate isolation study changed one at a
time. The residual correlation accounted for the liberal size. The
software and the degrees-of-freedom method had no effect.

In the toy, G1 rejected at 0.050 (MCSE 0.005). The toy generates
residuals that are exactly a random intercept plus AR(1), so G1 is
the correctly specified model there. This difference between the toy
and the extended DGP matters in the next section.

**Lesson.** Model the correlation pattern you expect in the data.
Fit the time scale on which it decays.

## 7. Rung D: cluster-robust standard errors (G6)

### 7.1 A model can be close and still wrong

G1's size of 0.038 in paper 02 is below 0.05. The test is
conservative. The reason is that the extended DGP generates a
covariance across visits that no random intercept plus AR(1)
residual can reproduce exactly: three components, each with its own
variance and its own time course. G1's working covariance is an
approximation. When the working covariance is wrong, the model-based
standard error is wrong, even though the point estimate is still
unbiased.

The project measures the error with the **calibration ratio**

$$
\kappa = \frac{\text{mean of the model-based standard errors}}
              {\text{empirical SD of the estimates across replicates}}.
$$

A correct standard error gives $\kappa = 1$. A value above 1 means
the model overstates its own uncertainty and the test is
conservative. Paper 10 found $\kappa$ of about 1.05 to 1.07 for the
model-based test under the extended DGP. Most of the excess came from
a slightly misfitted time trend, not from the covariance.

### 7.2 The sandwich idea

Instead of trusting the model's covariance, estimate the variance of
$\hat\beta$ from the data. The patients are independent of one
another, even though the visits within a patient are not. For a
generalized least squares estimator with working covariance $V_i$
for patient $i$,

$$
\widehat{\mathrm{Var}}(\hat\beta)
 = \underbrace{\Big(\textstyle\sum_i X_i^\top V_i^{-1} X_i\Big)^{-1}}_{\text{bread}}
   \underbrace{\Big(\textstyle\sum_i X_i^\top V_i^{-1}\, r_i r_i^\top\, V_i^{-1} X_i\Big)}_{\text{meat}}
   \underbrace{\Big(\textstyle\sum_i X_i^\top V_i^{-1} X_i\Big)^{-1}}_{\text{bread}},
$$

where $X_i$ is patient $i$'s design matrix and $r_i$ is the vector
of that patient's residuals. If the working covariance is correct,
the meat equals the inverse bread on average and the sandwich reduces
to the model-based variance. If it is wrong, the meat still measures
the actual spread, because it uses the observed residuals. This is
the **cluster-robust** or **sandwich** variance.

### 7.3 Small samples: CR2 and Satterthwaite

The plain sandwich is biased downward when the number of clusters is
small. Residuals are on average smaller than the true errors, because
the fit has been pulled toward the data. With 35 or 70 patients this
bias matters. Two corrections fix it.

- **CR2** (Bell and McCaffrey, 2002) inflates each patient's
  residuals by a matrix that undoes the shrinkage the fit imposed. It
  plays the same role as dividing by $n - 1$ instead of $n$.
- **Satterthwaite degrees of freedom** replace the normal reference
  distribution by a $t$ distribution whose degrees of freedom are
  estimated from the data. With few clusters, or with a coefficient
  that only a few patients inform, the degrees of freedom can be far
  below $N$.

The R package `clubSandwich` provides both for `nlme::lme` fits.
**G6** is G1 with CR2 standard errors and Satterthwaite degrees of
freedom. The point estimate is unchanged. Only the standard error,
the reference distribution and the $p$-value change.

### 7.4 What happened

In paper 02, CR2 moved the size from 0.038 to 0.055. A test with
correct size can reject more often when the alternative is true, so
CR2 also raised power.

The toy behaved differently, and the difference is instructive. G1 is
correct in the toy, so CR2 has nothing to repair. G6 rejected at
0.058 (MCSE 0.005) and lost 0.026 of power against G1 (paired SE
0.007). A robust standard error is noisier than a correct
model-based one, and its Satterthwaite degrees of freedom are
smaller. When the working covariance is right, CR2 costs some power.
When it is wrong, as under the extended DGP, CR2 restores the size,
and for a conservative test it restores power as well. A real analyst
cannot tell which case holds. Paper 02 therefore accepts the small
cost in the first case to be protected in the second.

**Lesson.** Use the model to get efficient estimates. Do not rely on
it alone for their uncertainty when the number of patients is
moderate and the covariance model is an approximation.

## 8. Rung E: carryover and the exposure-weighted predictor (G8)

### 8.1 The problem

A drug does not stop working the moment the patient stops taking it.
Its effect decays. Suppose the decay is exponential with half-life
$t_{1/2}$. One half-life after stopping, half the effect remains.
Two half-lives after, a quarter remains.

The binary indicator $D_{it}$ codes every off-drug visit as 0. In
path 1 at a one-week half-life, the visit at week 11 still carries
half the drug effect, and week 12 carries a quarter. The analysis
treats them as fully off drug. This is **measurement error in a
predictor**. Its classical consequence is **attenuation**: the
interaction estimate is pulled toward zero. The visits also vary in a
way the model does not represent, which inflates the residual
variance. Power falls for both reasons.

Carryover misspecification costs power. It does not cost size. Under
the null there is no interaction to attenuate, so the size of the
test does not depend on how exposure is coded.

### 8.2 The remedy

Replace $D_{it}$ by a predictor that decays as the drug does:

$$
D_{bc,it} =
\begin{cases}
1 & \text{on drug,}\\
(1/2)^{t_{sd,it} / t_{1/2}^a} & \text{off drug after first exposure,}\\
0 & \text{before first exposure.}
\end{cases}
$$

Here $t_{1/2}^a$ is the half-life the analysis assumes. Table 3
shows the values for path 3 at $t_{1/2}^a = 1$ week.

| Week | 4 | 8 | 9 | 10 | 11 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|
| $D_{it}$ | 1 | 1 | 1 | 0 | 0 | 0 | 1 | 0 |
| $t_{sd}$ | 0 | 0 | 0 | 1 | 2 | 3 | 0 | 4 |
| $D_{bc,it}$ | 1 | 1 | 1 | 0.5 | 0.25 | 0.125 | 1 | 0.0625 |

Table: Binary and exposure-weighted drug indicators, path 3,
assumed half-life 1 week.

The third case matters. Without it, a patient who starts on placebo
would have $t_{sd}$ undefined or zero and could be coded as exposed
before ever taking the drug. In the Hybrid design every path starts
on drug, so the case never arises there, but it does in a crossover.

As $t_{1/2}^a \to 0$, $D_{bc,it} \to D_{it}$. The binary indicator is
the zero-half-life member of this family.

**G3** is G1 with $D_{bc,it}$ in place of $D_{it}$, and **G8** is G3
with CR2 inference. In rung E the analysis is given the true half-life.
No real analyst knows it. Rung E measures how much power correct
exposure coding can recover. Rung F removes the unrealistic
knowledge.

### 8.3 What happened

In paper 02, moving from G6 to G8 raised power by 0.053 (paired SE
0.013) on the same simulated trials. In the toy, G8
gained 0.048 of power over G6 (paired SE 0.008), and the model-based
G3 gained 0.083 over G1 (paired SE 0.008). Size did not change: the
paired difference in null rejection between G8 and G6 was $-0.001$
(SE 0.005). This is the pattern Section 8.1 predicts. Exposure coding
changes power, not size.

## 9. Rung F: choosing the half-life from the data (G9)

### 9.1 The problem

Pharmacology gives a half-life for the drug in blood. The half-life of
its effect on symptoms may be longer, and it may differ by patient.
The analyst must choose $t_{1/2}^a$ without knowing the truth.

### 9.2 Model selection by AIC

Fit the model at each value of a candidate grid,
$\mathcal T = \{0.25, 0.5, 1, 2\}$ weeks, and keep the one with the
smallest **Akaike information criterion**,

$$
\mathrm{AIC} = -2 \log L + 2k,
$$

where $L$ is the maximized likelihood and $k$ the number of
parameters. Every candidate has the same number of parameters, so
here AIC simply picks the best-fitting half-life.

Two technical points matter.

- **Compare candidates by ML, then refit by REML.** REML maximizes the
  likelihood of contrasts that depend on the fixed-effects design
  matrix. Candidates that differ in a fixed-effects column, here
  $D_{bc,it}$, have REML likelihoods of different data. They cannot
  be compared. ML likelihoods can. After the choice, refit the chosen
  model by REML for the final estimate.
- **Report one fit.** The interaction estimate, standard error and
  $p$-value come from the single REML refit at the chosen half-life.

**G4** is this procedure with model-based inference. **G9** is the
same with CR2 inference.

### 9.3 The cost: post-selection inference

The same data choose the half-life and then test the interaction at
it. The reported $p$-value pretends the half-life was fixed in
advance. Selecting among several analyses and reporting the best is
a classic route to liberal tests. The concern is real and must be
measured.

In paper 02, pooled over 7,987 null fits, G9 rejected at 0.055 and
G8 at 0.049. The paired difference was 0.006 with standard error
0.001. Selection costs about half a percentage point of size. G9
stayed within 0.04 to 0.06 in all eight scenario families tested, and
it matched the power of G8, which was given the true half-life. In
the toy, G9 rejected a true null at 0.056, the same as G8 (paired
difference $-0.001$, SE 0.002). Its power was 0.450 against 0.448
for G8 (paired difference 0.003, SE 0.003). Under the alternative,
AIC chose the true half-life of 1 week in 1,552 of 2,000 trials,
0.5 week in 232 and 2 weeks in 216. The toy has only four candidates
and a clear signal, so selection cost nothing measurable here.
Paper 02's larger grid of scenarios detected a cost of about half a
point.

Two stronger fixes exist and are discussed in `docs/57`. One selects
the half-life on a model *without* the interaction, so the selection
cannot favor the coefficient under test (the project's F1A0). The
other calibrates the whole select-then-test procedure by a
parametric bootstrap (F1A0B). They cost more computing and, in the
project's simulations, gained no power over a well-chosen fixed
half-life.

## 10. The ladder at a glance

| Rung | Code | What changes | Question answered |
|---|---|---|---|
| | G5 (E9) | per-patient difference regressed on $B_i$ | simplest valid analysis |
| A | | OLS on every visit | ignores clustering |
| B | E0 | random intercept | clustering, compound symmetry |
| C | G1 | AR(1) residuals in calendar weeks | correlation fades with time |
| D | G6 | CR2 and Satterthwaite df | working covariance is approximate |
| E | G8 | exposure-weighted $D_{bc}$ at the true half-life | carryover attenuates the interaction |
| F | G9 | half-life chosen by ML AIC, REML refit | the true half-life is unknown |

Table: The analysis ladder. Each rung changes one thing from the rung
above it, except B to C in paper 02, where E0 and G1 also differ in
software and an expectancy term.

Each rung is a response to evidence, not a decoration. Rungs B and C
fix the size. Rung D protects the size when the working covariance is
wrong. It returns the power a conservative test wasted, at a small
cost in power when the covariance happens to be right. Rung E recovers power lost to
carryover. Rung F removes an assumption no analyst can satisfy, at a
measured and small cost in size.

## 11. G9 in R

The following code fits G9 to one data set `dat` in long format, one
row per patient visit, with columns `ptID`, `t` (week), `on` (binary
drug state), `tsd` (weeks since stopping) and `bm` (centered
biomarker). It mirrors `select_half_life_aic()` and `cr2_extract()` in
`analysis/scripts/carryover-sensitivity/simulation-core.R`.

```r
library(nlme)
library(clubSandwich)

exposure_weight <- function(on, tsd, half_life) {
  ifelse(on == 1, 1, 0.5^(tsd / half_life))
}

fit_g <- function(dat, method) {
  lme(Sx ~ bm + t + Dbc + bm:Dbc,
      random = ~ 1 | ptID,
      correlation = corCAR1(form = ~ t | ptID),
      data = dat, method = method)
}

grid <- c(0.25, 0.5, 1, 2)
aics <- purrr::map_dbl(grid, function(h) {
  dat$Dbc <- exposure_weight(dat$on, dat$tsd, h)
  AIC(fit_g(dat, 'ML'))
})
h_hat <- grid[which.min(aics)]

dat$Dbc <- exposure_weight(dat$on, dat$tsd, h_hat)
fit <- fit_g(dat, 'REML')

coef_test(fit, vcov = 'CR2', cluster = dat$ptID,
          test = 'Satterthwaite')
```

The row `bm:Dbc` of the last output gives $\hat\beta_{bm:D}$, its CR2
standard error, the Satterthwaite degrees of freedom and the
$p$-value. Report $\hat t_{1/2}^a$ beside it. This sketch omits the
pre-exposure rule because every Hybrid path starts on drug. A design
with placebo-first paths must set $D_{bc,it} = 0$ before first
exposure.

## 12. What G9 does not solve

A tutorial that ended at Section 11 would overstate the case. The
following limits are documented in the project.

1. **It is a recommendation for one design under one DGP.** Paper 02
   recommends G9 for the Hybrid design under the correlation encoding
   of the interaction. It does not establish G9 for other designs or
   other data-generating processes.
2. **A prognostic biomarker can break it.** If the biomarker also
   predicts the placebo response or the natural-history trend, a
   model without terms for those pathways can attribute them to the
   interaction. The project's F1 model adds the design expectancy
   $\eta_t$ (1 open label, 0.5 blinded) and its interaction with the
   biomarker for this reason (`docs/49`, `docs/56`). Paper 02 found
   that a biomarker predicting the natural-history trend biased every
   analysis it tested.
3. **Selection is not fully accounted for.** G9's size of about 0.055
   is close to nominal but not exactly nominal. A bootstrap-calibrated
   procedure (F1A0B) holds 0.05 at higher computing cost.
4. **The working covariance is still an approximation.** CR2 repairs
   the standard error. It does not make the estimator efficient. A
   better covariance model could give more power.
5. **Simulation evidence is conditional.** Every number in this
   tutorial comes from a simulated world whose rules we wrote. A real
   trial can violate them in ways no simulation anticipated.

## 13. Exercises

1. Show that in model (1) the drug effect for a patient with
   biomarker $B_i$ is $\beta_D + \beta_{bm:D} B_i$. Show that
   centering $B_i$ changes $\beta_D$ but not $\beta_{bm:D}$.
2. Show that with patient indicator variables in the model, the
   column for $B_i$ is a linear combination of the indicators. Explain
   why $B_i D_{it}$ is not.
3. Derive the intraclass correlation of Section 5.2 from model (2).
4. With $\rho = 0.7$, compute the AR(1) correlation between every pair
   of Hybrid visits. Which pairs are most alike? Compare with a
   compound-symmetry correlation of 0.5.
5. Compute $D_{bc,it}$ for path 4 at $t_{1/2}^a = 0.5$ and 2 weeks.
   At which half-life does $D_{bc}$ differ most from $D_{it}$?
6. Run `analysis/scripts/tutorial-g9/01-toy-ladder.R` with 200
   replicates. Compute the MCSE of each rejection rate. Which
   differences in Table 2 are distinguishable from Monte Carlo noise?
7. Change the true half-life in the toy to 0.25 week and rerun. Does
   rung E still gain power over rung D? Explain using Section 8.1.
8. Explain in words why the AIC comparison in Section 9.2 must use ML
   rather than REML.
9. Modify the toy so that the biomarker also shifts the time trend,
   `beta_time = -0.2 + 0.1 * b`. Which rungs become liberal? Relate
   your answer to Section 12, item 2.

## References

- Bell RM, McCaffrey DF (2002). Bias reduction in standard errors for
  linear regression with multi-stage samples. *Survey Methodology*
  28(2):169-181.
- Hendrickson RC, Thomas RG, Schork NJ, Raskind MA (2020). Optimizing
  aggregated N-of-1 trial designs for predictive biomarker validation:
  statistical methods and theoretical findings. *Frontiers in Digital
  Health* 2:13.
- Pustejovsky JE, Tipton E (2018). Small-sample methods for
  cluster-robust variance estimation and hypothesis testing in fixed
  effects models. *Journal of Business and Economic Statistics*
  36(4):672-683.
- Satterthwaite FE (1946). An approximate distribution of estimates of
  variance components. *Biometrics Bulletin* 2(6):110-114.
- Project documents: paper 02 (`analysis/report/02-carryover-sensitivity/`,
  rendered 2026-10-29), paper 10
  (`analysis/report/10-interaction-test-calibration/`), `docs/40`
  (Hybrid design primer), `docs/49` and `docs/56` (expectancy terms),
  `docs/57` (data-chosen half-life), `docs/62` (G- and F-family
  catalog), `docs/63` (why correlation fades with time).
