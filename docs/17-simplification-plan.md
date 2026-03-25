# Simplification Plan for pmsimstats-simple

## Overview

This document outlines possible simplifications to the N-of-1 trial simulation
for biomarker-moderated treatment response. The goal is to reduce complexity
while preserving the core scientific question: **What is the statistical power
to detect a biomarker × treatment interaction across different trial designs?**

## Current Complexity in pmsimstats2025

| Component | Current Implementation | Complexity |
|-----------|----------------------|------------|
| Covariance structure | Full multivariate normal with AR(1) temporal correlation, cross-component correlations | High |
| Data generation | Two-stage conditional sampling from joint distribution | High |
| Analysis model | Mixed-effects (lme) with CAR(1) correlation | High |
| Carryover | Exponential decay via half-life parameter | Medium |
| Censoring | 5 different dropout/censoring patterns | Medium |
| Response components | BR, ER, TR with separate random effects per timepoint | High |

## Proposed Simplifications

### 1. Data Generation

#### 1a. Covariance Structure

**Current**: Joint multivariate normal for BR, ER, TR random effects across all
timepoints, plus biomarker and baseline, with AR(1) temporal structure and
cross-component correlations.

**Simplified options**:

| Option | Description | Pros | Cons |
|--------|-------------|------|------|
| Independent errors | `response = signal + rnorm(1, 0, sigma)` | Trivial to implement | Ignores repeated measures |
| Compound symmetry | Equal correlation between all timepoints | Simple, one parameter | Unrealistic for longitudinal |
| Random intercept only | `response = signal + u_i + epsilon_ij` | Standard mixed model | No temporal correlation |
| AR(1) residuals only | Single AR(1) process, no component separation | Moderate complexity | Loses BR/ER/TR decomposition |

**Recommendation**: Start with **random intercept + independent errors**. This
captures between-subject variability without temporal complexity.

```r
# Simplified generation
u_i <- rnorm(1, 0, sigma_between)  # random intercept
epsilon <- rnorm(n_timepoints, 0, sigma_within)  # independent errors
response <- baseline + signal + u_i + epsilon
```

#### 1b. Response Components

**Current**: Separate BR, ER, TR components each with own random effects.

**Simplified options**:

| Option | Description | Implementation |
|--------|-------------|----------------|
| Single treatment effect | Collapse BR into one effect | `effect = on_drug * beta * (1 + bm_mod * biomarker)` |
| No ER component | Remove expectancy response | Reduces parameters by 1/3 |
| No TR component | Remove time/regression to mean | Assumes stable baseline |
| Combined ER+TR | Single "non-drug" trajectory | `non_drug = week * rate` |

**Recommendation**: Keep BR with biomarker moderation (core hypothesis).
Consider dropping TR initially; add back if baseline drift is important.

### 2. Statistical Analysis

#### 2a. Model Choice

**Current**: `nlme::lme()` with CAR(1) correlation structure.

**Simplified options**:

| Option | Formula | Pros | Cons |
|--------|---------|------|------|
| ANCOVA | `lm(post ~ pre + drug * biomarker)` | Fast, familiar, closed-form | Requires summary measure |
| t-test on change | `t.test(delta ~ group)` | Simplest | Ignores biomarker |
| GEE | `geepack::geeglm()` | Robust to correlation misspecification | Slightly complex |
| Random intercept LMM | `lme4::lmer(y ~ x + (1|id))` | Simpler than full lme | Still mixed model |
| Simple linear regression | `lm(response ~ on_drug * biomarker)` | Ignores clustering | Biased SE but fast |
| RM-ANOVA | `aov(y ~ bm_group * drug + Error(id/drug))` | Familiar to clinicians, handles within-subject factor | Requires dichotomized biomarker (~36% information loss); sphericity assumption |

**Recommendation**: Use **ANCOVA on change scores** as primary analysis:

```r
# Compute within-subject change (or mean during treatment vs control)
summary_data <- trial_data |>
  group_by(participant_id, on_drug) |>
  summarise(mean_response = mean(response), .groups = "drop") |>
  pivot_wider(names_from = on_drug, values_from = mean_response) |>
  mutate(change = `TRUE` - `FALSE`)

# ANCOVA
model <- lm(change ~ biomarker, data = summary_data)
p_value <- coef(summary(model))["biomarker", "Pr(>|t|)"]
```

#### 2b. Handling Repeated Measures

**Current**: Full longitudinal model with all timepoints.

**Simplified options**:

| Option | Description |
|--------|-------------|
| Summary statistics | Mean/AUC during each phase |
| Last observation | Use final timepoint only |
| First vs last | Compare endpoints |
| Phase means | Average within treatment phases |

**Recommendation**: **Phase means** - average response during on-drug vs
off-drug periods, then analyze the difference.

### 2c. Mathematical Specification of the Interaction Test

This section formally defines the biomarker-treatment interaction test
under each analysis method and connects the precision of the interaction
estimate to trial design features: the number of subjects ($N$), the
number of assessment timepoints per subject ($T$), the number of
within-subject treatment transitions (crossovers), and the proportion of
time spent on vs. off drug.

#### Notation

| Symbol | Definition |
|--------|-----------|
| $Y_{it}$ | Observed outcome for subject $i$ at timepoint $t$ |
| $bm_i$ | Baseline biomarker value for subject $i$ (time-invariant) |
| $D_{it}$ | Drug exposure indicator (1 = on drug, 0 = off drug) at timepoint $t$ |
| $t$ | Time (weeks since baseline) |
| $u_i$ | Subject-specific random intercept, $u_i \sim N(0, \sigma^2_u)$ |
| $\varepsilon_{it}$ | Residual error, $\varepsilon_{it} \sim N(0, \sigma^2_\varepsilon)$ |
| $N$ | Number of subjects |
| $T_i$ | Total number of assessment timepoints for subject $i$ |
| $n_{i,\text{on}}$ | Number of on-drug assessments for subject $i$ |
| $n_{i,\text{off}}$ | Number of off-drug assessments for subject $i$ |
| $K_i$ | Number of treatment transitions (drug-to-placebo or placebo-to-drug) for subject $i$ |

Throughout, the **target hypothesis** is:

$$H_0: \text{the biomarker does not predict the magnitude of
drug-specific response}$$

The mathematical form of this hypothesis differs depending on whether
the design includes within-subject drug variation.

---

#### The fundamental distinction: within-subject drug variation

**Designs with within-subject drug variation** (OL+BDC, CO, N-of-1):
Subject $i$ has observations both on drug ($D_{it} = 1$) and off drug
($D_{it} = 0$). The interaction is $bm_i \times D_{it}$, testing whether
the biomarker predicts the within-subject drug effect.

**Designs without within-subject drug variation** (OL): All observations
for subject $i$ have $D_{it} = 1$. There is no within-subject drug
contrast. The interaction must be reformulated as $bm_i \times t$,
testing whether the biomarker predicts the rate of improvement over time
while on drug. This is a fundamentally weaker test because temporal
improvement confounds BR, ER, and TV.

This distinction is the primary reason the open-label design has low
power (Hendrickson Figure 4A): the interaction term targets a different,
noisier quantity.

---

#### Method 1: ANCOVA on phase-mean change scores

**Step 1.** Collapse repeated measures into within-subject phase means:

$$\bar{Y}_{i,\text{on}} = \frac{1}{n_{i,\text{on}}}
\sum_{t: D_{it}=1} Y_{it}, \qquad
\bar{Y}_{i,\text{off}} = \frac{1}{n_{i,\text{off}}}
\sum_{t: D_{it}=0} Y_{it}$$

**Step 2.** Compute the within-subject drug effect:

$$\Delta_i = \bar{Y}_{i,\text{on}} - \bar{Y}_{i,\text{off}}$$

**Step 3.** Regress the drug effect on the biomarker:

$$\Delta_i = \gamma_0 + \gamma_1 \cdot bm_i + e_i$$

**Interaction test:** $H_0: \gamma_1 = 0$ via $t$-test on
$\hat{\gamma}_1$, with $N - 2$ degrees of freedom.

**Precision of $\hat{\gamma}_1$:**

$$\text{Var}(\hat{\gamma}_1) =
\frac{\sigma^2_\Delta}{N \cdot \text{Var}(bm)}$$

where $\sigma^2_\Delta = \text{Var}(\Delta_i)$ is the between-subject
variance of the within-subject drug effect. This variance decomposes as:

$$\sigma^2_\Delta = \text{Var}(\bar{Y}_{i,\text{on}}) +
\text{Var}(\bar{Y}_{i,\text{off}}) -
2\,\text{Cov}(\bar{Y}_{i,\text{on}},
\bar{Y}_{i,\text{off}})$$

Each phase-mean variance is approximately:

$$\text{Var}(\bar{Y}_{i,\text{phase}}) \approx
\sigma^2_u + \frac{\sigma^2_\varepsilon}{n_{i,\text{phase}}}$$

The covariance term $\text{Cov}(\bar{Y}_{i,\text{on}},
\bar{Y}_{i,\text{off}})$ arises from the shared random intercept $u_i$
and equals $\sigma^2_u$ under compound symmetry.

**Substituting:**

$$\sigma^2_\Delta \approx
\frac{\sigma^2_\varepsilon}{n_{i,\text{on}}} +
\frac{\sigma^2_\varepsilon}{n_{i,\text{off}}}
= \sigma^2_\varepsilon \left(
\frac{1}{n_{i,\text{on}}} + \frac{1}{n_{i,\text{off}}}
\right)$$

The random intercept cancels in the difference, which is a key advantage
of within-subject contrasts. What remains is purely residual variance,
reduced by the number of observations in each phase.

**Design implications:**

| Design feature | Effect on $\text{Var}(\hat{\gamma}_1)$ |
|:---------------|:---------------------------------------|
| More subjects ($N \uparrow$) | Direct reduction: $\propto 1/N$ |
| More on-drug assessments ($n_\text{on} \uparrow$) | Reduces $\sigma^2_\Delta$ via $1/n_\text{on}$ term |
| More off-drug assessments ($n_\text{off} \uparrow$) | Reduces $\sigma^2_\Delta$ via $1/n_\text{off}$ term |
| Balanced phases ($n_\text{on} \approx n_\text{off}$) | Minimizes the harmonic mean, hence $\sigma^2_\Delta$ |
| More crossovers ($K \uparrow$) | No direct effect beyond increasing $n_\text{on}$ and $n_\text{off}$; however, distributing observations across more transitions reduces sensitivity to period effects |
| Wider biomarker range ($\text{Var}(bm) \uparrow$) | Direct reduction: $\propto 1/\text{Var}(bm)$ |

**Design-specific behavior:**

- **OL:** $n_{i,\text{off}} = 0$, so $\Delta_i$ is undefined. The ANCOVA
  must be reformulated as $\bar{Y}_i = \gamma_0 + \gamma_1 \cdot bm_i + e_i$
  (regress total response on biomarker), which confounds BR with ER and TV.
- **OL+BDC:** $n_{i,\text{on}} \approx 5$, $n_{i,\text{off}} \approx 2$
  (the blinded discontinuation contributes only 2 off-drug assessments for
  half of participants). Imbalanced phases inflate $\sigma^2_\Delta$.
- **CO:** $n_{i,\text{on}} \approx 4$, $n_{i,\text{off}} \approx 4$.
  Balanced phases. One crossover transition.
- **N-of-1:** $n_{i,\text{on}} \approx 6$, $n_{i,\text{off}} \approx 4$.
  Multiple crossover transitions. Highest total observation count and
  reasonable balance.

**OL special case.** Without drug variation, the ANCOVA tests a different
quantity. One approach is to use the slope of improvement over time as
the subject-level summary:

$$b_i = \text{slope of } Y_{it} \text{ vs. } t \text{ for subject } i$$

$$b_i = \gamma_0 + \gamma_1 \cdot bm_i + e_i$$

Here $\gamma_1$ tests whether the biomarker predicts the rate of
improvement. This conflates BR rate with ER and TV rates, explaining
the poor power of the OL design.

---

#### Method 2: Paired $t$-test on change scores

The paired $t$-test applies to the same $\Delta_i$ as ANCOVA but tests a
different hypothesis:

$$H_0: E[\Delta_i] = 0 \quad \text{(no average drug effect)}$$

$$t = \frac{\bar{\Delta}}{s_\Delta / \sqrt{N}}, \quad df = N - 1$$

**This does not test the biomarker interaction.** It tests whether the
drug works on average, collapsing over biomarker values. To recover a
biomarker interaction test from a $t$-test framework, one must
dichotomize:

$$\text{Split subjects into } bm_\text{high} \text{ and }
bm_\text{low} \text{ by median}$$

$$t = \frac{\bar{\Delta}_\text{high} -
\bar{\Delta}_\text{low}}{s_p \sqrt{1/N_\text{high} +
1/N_\text{low}}}, \quad df = N - 2$$

This is equivalent to ANCOVA with a binary biomarker and inherits the
same design dependencies. Dichotomization discards information, reducing
power relative to ANCOVA with a continuous biomarker by approximately
$2/\pi \approx 64\%$ relative efficiency (Cohen 1983).

**Design implications:** Same as ANCOVA, but uniformly less powerful due
to dichotomization. Not recommended when the biomarker is continuous.

---

#### Method 3: GEE (generalized estimating equations)

**Model (population-averaged):**

$$E[Y_{it}] = \beta_0 + \beta_1 \cdot bm_i + \beta_2 \cdot D_{it}
+ \beta_3 \cdot t + \beta_4 \cdot bm_i \cdot D_{it}$$

with working correlation matrix $R_i(\alpha)$ and sandwich (robust)
variance estimator.

**Interaction test:** $H_0: \beta_4 = 0$ via Wald test using the
sandwich-estimated variance:

$$z = \frac{\hat{\beta}_4}{\widehat{\text{SE}}_\text{robust}
(\hat{\beta}_4)}$$

**Precision of $\hat{\beta}_4$:**

The sandwich estimator is consistent regardless of the working
correlation structure, but efficiency depends on how well the working
correlation approximates the truth. Under a correct working correlation:

$$\text{Var}(\hat{\beta}_4) \approx
\frac{1}{N \cdot \text{Var}(bm) \cdot
\text{Var}_\text{within}(D)}
\cdot \sigma^2_\text{eff}$$

where $\text{Var}_\text{within}(D)$ is the within-subject variance of
drug exposure and $\sigma^2_\text{eff}$ is the effective residual
variance after accounting for the working correlation.

**Within-subject variance of drug exposure** is the critical
design-dependent quantity:

$$\text{Var}_\text{within}(D_i) =
\bar{D}_i(1 - \bar{D}_i)$$

where $\bar{D}_i = n_{i,\text{on}} / T_i$ is the proportion of time on
drug. This is maximized when $\bar{D}_i = 0.5$ (equal time on and off
drug) and equals zero when $\bar{D}_i = 1$ (all observations on drug,
i.e. the OL design).

**Design implications:**

| Design feature | Effect on $\text{Var}(\hat{\beta}_4)$ |
|:---------------|:--------------------------------------|
| More subjects ($N \uparrow$) | Direct reduction: $\propto 1/N$ |
| More timepoints ($T \uparrow$) | Reduces $\sigma^2_\text{eff}$ (more observations to estimate correlation); modest gain under strong autocorrelation |
| Balanced drug exposure ($\bar{D} \to 0.5$) | Maximizes $\text{Var}_\text{within}(D)$, directly reducing $\text{Var}(\hat{\beta}_4)$ |
| More crossovers ($K \uparrow$) | Does not change $\text{Var}_\text{within}(D)$ directly, but distributes drug variation across time, improving robustness to temporal confounds |
| Autocorrelation ($\rho \uparrow$) | Increases effective variance; observations close in time contribute less independent information |

**Design-specific behavior:**

- **OL:** $\text{Var}_\text{within}(D) = 0$. The interaction $bm \times D$
  is not estimable. Must substitute $bm \times t$ (same limitation as
  ANCOVA).
- **OL+BDC:** $\bar{D} \approx 0.75$ (most time on drug).
  $\text{Var}_\text{within}(D) \approx 0.19$. Suboptimal but nonzero.
- **CO:** $\bar{D} = 0.5$. $\text{Var}_\text{within}(D) = 0.25$.
  Optimal balance. One transition.
- **N-of-1:** $\bar{D} \approx 0.55$ (open-label phase tips the balance
  slightly toward on-drug). $\text{Var}_\text{within}(D) \approx 0.25$.
  Near-optimal, with multiple transitions.

**GEE vs. LMM:** GEE estimates population-averaged effects (marginal
model), while LMM estimates subject-specific effects (conditional model).
For the linear case with continuous outcomes, $\beta_4$ has the same
interpretation under both approaches. The practical difference is that
GEE provides valid inference under working-correlation misspecification
(via the sandwich estimator), while LMM requires correct specification
of the random effects structure but is more efficient when correctly
specified.

---

#### Method 4: Random intercept LMM (lme4::lmer)

**Model (subject-specific):**

$$Y_{it} = \beta_0 + \beta_1 \cdot bm_i + \beta_2 \cdot D_{it}
+ \beta_3 \cdot t + \beta_4 \cdot bm_i \cdot D_{it}
+ u_i + \varepsilon_{it}$$

$$u_i \sim N(0, \sigma^2_u), \qquad
\varepsilon_{it} \sim N(0, \sigma^2_\varepsilon)$$

**Interaction test:** $H_0: \beta_4 = 0$ via Wald test (or Kenward-Roger
/ Satterthwaite $F$-test for small $N$).

**Precision of $\hat{\beta}_4$:**

Under the random-intercept model with independent residuals, the
information matrix for $\beta_4$ depends on the 'effective' within-subject
information after projecting out the random intercept. The approximate
variance is:

$$\text{Var}(\hat{\beta}_4) \approx
\frac{\sigma^2_\varepsilon}{
\sum_{i=1}^{N} \widetilde{bm}_i^2 \cdot
\text{SS}_\text{within}(D_i)}$$

where $\widetilde{bm}_i = bm_i - \overline{bm}$ is the centered
biomarker value and:

$$\text{SS}_\text{within}(D_i) =
\sum_{t=1}^{T_i} (D_{it} - \bar{D}_i)^2
= T_i \cdot \bar{D}_i (1 - \bar{D}_i)$$

is the within-subject sum of squares of drug exposure for subject $i$.
This quantity has a direct interpretation: it is the total 'information'
about drug variation contained in subject $i$'s data.

**Design implications:**

| Design feature | Effect on $\text{Var}(\hat{\beta}_4)$ | Mechanism |
|:---------------|:--------------------------------------|:----------|
| More subjects ($N \uparrow$) | $\propto 1/N$ | More independent biomarker-by-drug contrasts |
| More timepoints ($T \uparrow$) | $\propto 1/T$ (given fixed $\bar{D}$) | Each additional timepoint adds to $\text{SS}_\text{within}(D_i)$ |
| Balanced drug exposure ($\bar{D} \to 0.5$) | Maximizes $\bar{D}(1 - \bar{D})$, hence $\text{SS}_\text{within}(D)$ | Equal on/off time maximizes within-subject contrast |
| More crossovers ($K \uparrow$) | No direct effect on $\text{SS}_\text{within}(D)$ under independence | Under AR(1) residuals, more transitions reduce the effective autocorrelation penalty (see below) |
| Wider biomarker range | Increases $\sum \widetilde{bm}_i^2$ | More leverage for the interaction term |
| Autocorrelation ($\rho > 0$) | Inflates effective $\sigma^2_\varepsilon$ | Correlated residuals reduce the effective sample size per subject |

**The autocorrelation-crossover interaction.** Under AR(1) residuals
(as implemented in the full nlme::lme model with corCAR1), the effective
information from each on-drug or off-drug block depends on block length
$L$ and autocorrelation $\rho$:

$$\text{effective observations per block} \approx
\frac{L}{1 + (L-1)\rho} \cdot L$$

Short blocks (high $K$, low $L$) suffer less from autocorrelation than
long blocks (low $K$, high $L$) because fewer consecutive observations
are highly correlated. This is why the N-of-1 design (four 4-week
blocks) can outperform the crossover (two 10-week blocks) when
autocorrelation is moderate, despite having the same total duration.

Conversely, when the drug has substantial carryover, each transition
introduces a contaminated observation at the block boundary. More
transitions ($K \uparrow$) means more contaminated observations. This
is the mechanism behind the 'precipitous decline in power' with
carryover that Hendrickson et al. report for the N-of-1 design (their
Figure 4B).

**Design-specific $\text{SS}_\text{within}(D)$:**

| Design | $T$ | $n_\text{on}$ | $n_\text{off}$ | $\bar{D}$ | $\text{SS}_\text{within}(D)$ | $K$ |
|--------|-----|---------------|-----------------|-----------|------------------------------|-----|
| OL | 8 | 8 | 0 | 1.0 | 0 | 0 |
| OL+BDC | 8 | ~6 | ~2 | ~0.75 | ~1.5 | 1 |
| CO | 8 | 4 | 4 | 0.5 | 2.0 | 1 |
| N-of-1 | 12 | ~7 | ~5 | ~0.58 | ~2.9 | 3 |

The N-of-1 design has the highest $\text{SS}_\text{within}(D)$ due to
both more total timepoints and near-balanced drug exposure. The OL
design has zero within-subject drug information.

---

#### Method 5: Simple linear regression (ignoring clustering)

**Model:**

$$Y_{it} = \beta_0 + \beta_1 \cdot bm_i + \beta_2 \cdot D_{it}
+ \beta_3 \cdot t + \beta_4 \cdot bm_i \cdot D_{it}
+ \varepsilon_{it}$$

This is the same fixed-effects structure as Method 4 but without the
random intercept $u_i$. All $N \times T$ observations are treated as
independent.

**Interaction test:** $H_0: \beta_4 = 0$ via OLS $t$-test.

**Precision of $\hat{\beta}_4$ (nominal):**

$$\text{Var}_\text{OLS}(\hat{\beta}_4) =
\frac{\hat{\sigma}^2}{\text{SS}(bm \cdot D \mid
\text{other predictors})}$$

where $\hat{\sigma}^2$ is the OLS residual variance. Because the model
ignores the random intercept, $\hat{\sigma}^2$ absorbs both
$\sigma^2_u$ and $\sigma^2_\varepsilon$, but the denominator treats all
$N \times T$ observations as independent. This produces **downward-biased
standard errors** because the effective sample size for the
between-subject component ($bm_i$) is $N$, not $N \times T$.

**Type I error inflation:** The degree of inflation depends on
$\sigma^2_u / \sigma^2_\varepsilon$ (the intraclass correlation) and the
number of observations per subject $T$. With $T = 8$ and ICC = 0.3
(typical for CAPS data), the nominal $\alpha = 0.05$ test has an actual
Type I error of approximately 0.15--0.25.

**Design implications:** The same as Method 4 for power ordering across
designs, but with inflated Type I error that invalidates nominal
$p$-values. Not recommended for inference but can serve as a fast
screening tool when followed by a properly specified model.

---

#### Method 6: Repeated measures ANOVA (biomarker dichotomized)

**Setup.** Because repeated measures ANOVA requires a grouping factor
(not a continuous covariate), the biomarker must be dichotomized:

$$G_i = \begin{cases} 1 & \text{if } bm_i > \text{median}(bm) \\
0 & \text{otherwise} \end{cases}$$

The data are arranged as a $N \times T$ matrix of outcomes, with two
between-subject groups ($G = 0, 1$) and a within-subject factor $D$
(on-drug vs. off-drug). If the design has multiple timepoints within
each drug condition, there is an additional within-subject factor for
time nested within drug phase.

**Model (split-plot / mixed ANOVA):**

$$Y_{it} = \mu + \alpha_{G_i} + \delta_{D_{it}} +
(\alpha\delta)_{G_i, D_{it}} + \pi_i + \varepsilon_{it}$$

where:

- $\mu$ is the grand mean
- $\alpha_{G_i}$ is the between-subject biomarker group effect
- $\delta_{D_{it}}$ is the within-subject drug effect
- $(\alpha\delta)_{G_i, D_{it}}$ is the biomarker-by-drug interaction
- $\pi_i \sim N(0, \sigma^2_\pi)$ is the subject effect (random)
- $\varepsilon_{it} \sim N(0, \sigma^2_\varepsilon)$ is the residual

**Interaction test:** $H_0: (\alpha\delta)_{G,D} = 0$ for all
combinations, tested via the $F$-statistic:

$$F = \frac{\text{MS}_{G \times D}}{\text{MS}_{\text{error}(within)}}$$

with $df_1 = (g-1)(d-1)$ and $df_2 = (N-g)(d-1)$, where $g = 2$
(biomarker groups) and $d = 2$ (on-drug, off-drug). Under the simple
two-group, two-condition case: $df_1 = 1$ and $df_2 = N - 2$.

**Sphericity.** When more than two levels of the within-subject factor
are present (e.g., multiple timepoints within each drug phase), the
$F$-test assumes sphericity: equal variances of all pairwise differences
between within-subject levels. Violation inflates the Type I error rate.
The Greenhouse-Geisser (GG) and Huynh-Feldt (HF) corrections multiply
both $df_1$ and $df_2$ by a correction factor
$\hat{\varepsilon} \in [1/(k-1), 1]$:

$$F_\text{corrected}: \quad df_1' = df_1 \cdot \hat{\varepsilon},
\quad df_2' = df_2 \cdot \hat{\varepsilon}$$

When the within-subject factor is simply on-drug vs. off-drug ($k = 2$),
sphericity is satisfied trivially (there is only one pairwise
difference), and no correction is needed. When timepoints within phases
are retained as separate levels, the correction becomes relevant.

**Precision of the interaction test:**

For the simple case ($g = 2$ groups, $d = 2$ drug conditions, collapsing
timepoints within phases into means):

$$\text{MS}_{G \times D} =
\frac{N}{4} \left(
(\bar{\Delta}_\text{high} - \bar{\Delta}_\text{low})
\right)^2$$

$$\text{MS}_\text{error(within)} =
\frac{1}{N-2} \sum_{i=1}^{N}
(\Delta_i - \bar{\Delta}_{G_i})^2$$

where $\Delta_i = \bar{Y}_{i,\text{on}} - \bar{Y}_{i,\text{off}}$ as
in the ANCOVA approach. The noncentrality parameter for the $F$-test is:

$$\lambda = \frac{N}{4} \cdot
\frac{(\mu_{\Delta,\text{high}} -
\mu_{\Delta,\text{low}})^2}{\sigma^2_\Delta}$$

**Comparison to ANCOVA (Method 1).** The repeated measures ANOVA on
phase means with dichotomized biomarker is algebraically equivalent to a
two-sample $t$-test on $\Delta_i$ between biomarker groups. The ANCOVA
on change scores with continuous biomarker tests the same conceptual
hypothesis but retains the full biomarker distribution. The efficiency
loss from dichotomization is well established: for a normally distributed
biomarker, the median split retains only $2/\pi \approx 63.7\%$ of the
information, corresponding to a relative efficiency of approximately
0.64 (Cohen 1983, Royston et al. 2006). In sample-size terms, the
repeated measures ANOVA requires approximately $N/0.64 \approx 1.57N$
subjects to achieve the same power as ANCOVA with a continuous biomarker.

**When timepoints are not collapsed.** If the full set of $T$ timepoints
is retained as within-subject levels (rather than collapsing into phase
means), the model becomes a two-way mixed ANOVA with $G$ (between) and
timepoint (within). The interaction $G \times \text{timepoint}$ tests
whether the two biomarker groups have different temporal profiles, which
is broader than the targeted $bm \times D$ interaction. This formulation:

- Gains power from using all observations (no information loss from
  averaging)
- Loses specificity by testing a composite hypothesis across all
  timepoints rather than the drug contrast specifically
- Requires the sphericity assumption across all $T$ timepoints, which
  is typically violated in longitudinal data and necessitates GG/HF
  correction
- The corrected degrees of freedom can substantially reduce power when
  $\hat{\varepsilon}$ is small (strong autocorrelation, many timepoints)

**Design implications:**

| Design feature | Effect on power | Mechanism |
|:---------------|:----------------|:----------|
| More subjects ($N \uparrow$) | $\propto N$ (via noncentrality $\lambda$) | More independent group comparisons |
| More on/off-drug assessments | Reduces $\sigma^2_\Delta$ when collapsing to phase means (same as ANCOVA) | More precise within-subject drug effect estimate |
| Balanced phases ($\bar{D} \to 0.5$) | Same benefit as ANCOVA | Minimizes variance of $\Delta_i$ |
| Balanced biomarker groups ($N_\text{high} \approx N_\text{low}$) | Maximizes $F$-test power | Median split guarantees this for even $N$ |
| More crossovers ($K \uparrow$) | Indirect benefit via $n_\text{on}$, $n_\text{off}$ | Same as ANCOVA |
| More timepoints retained as levels | Increases $df_1$ but triggers sphericity issues | GG/HF correction can offset gains |
| Autocorrelation ($\rho \uparrow$) | Reduces $\hat{\varepsilon}$ when timepoints retained | Corrected $df$ shrink, widening the rejection threshold |

**Design-specific behavior:**

- **OL:** No within-subject drug factor. The within-subject factor
  becomes time, and the interaction $G \times \text{time}$ tests whether
  biomarker groups have different slopes. Same confounding problem as all
  other methods.
- **OL+BDC:** The within-subject drug factor has two levels (on/off) but
  is heavily imbalanced. Collapsing to phase means is necessary; the
  sparse off-drug phase inflates $\sigma^2_\Delta$.
- **CO:** Clean two-level within-subject drug factor with balanced
  phases. Sphericity is trivially satisfied (two levels). The repeated
  measures ANOVA is most natural here.
- **N-of-1:** Multiple drug-on and drug-off blocks. If collapsed to
  phase means, behaves like ANCOVA with good balance. If blocks are
  retained as separate within-subject levels, sphericity becomes a
  concern.

**Summary.** Repeated measures ANOVA is conceptually accessible and
widely implemented, making it attractive for clinical audiences.
However, for the biomarker interaction test specifically, it suffers
from two compounding limitations: (1) dichotomization of the biomarker
discards approximately 36% of the information, and (2) when timepoints
are retained rather than collapsed, sphericity violations and
GG/HF corrections can further erode power. For these reasons, ANCOVA
on change scores (Method 1) or the random intercept LMM (Method 4)
are generally preferable when the biomarker is continuous.

---

#### The OL design: a fundamentally different interaction test

For the open-label design, all six methods face the same problem: there
is no within-subject drug variation. The interaction must be formulated
as biomarker-by-time:

$$Y_{it} = \beta_0 + \beta_1 \cdot bm_i + \beta_2 \cdot t
+ \beta_3 \cdot bm_i \cdot t + u_i + \varepsilon_{it}$$

Here $\beta_3$ tests whether the biomarker predicts the slope of
improvement over time. This is a weaker test than $\beta_4$ (the
biomarker-by-drug interaction) for two reasons:

1. **Confounding.** The slope $dY/dt$ includes BR, ER, and TV. The
   biomarker-by-time interaction conflates all three: $\beta_3$ detects a
   biomarker association with the total improvement rate, not just the
   drug-specific component.

2. **Precision.** The within-subject information for the time slope
   depends on $\text{SS}_\text{within}(t) = \sum(t - \bar{t})^2$, which
   is determined by the spacing and number of assessment points. This
   quantity is typically smaller than $\text{SS}_\text{within}(D)$ for
   designs with balanced crossovers, because the time variable has less
   'contrast' than a binary drug indicator that switches between 0 and 1.

---

#### Summary: interaction test precision across designs and methods

The following table summarizes the effective information for the
biomarker interaction test across all design-method combinations. The
entries are proportional to $1/\text{Var}(\hat{\beta}_\text{interaction})$
(higher = more power):

| | OL | OL+BDC | CO | N-of-1 |
|:----------|:---------|:---------|:---------|:---------|
| **Interaction term** | $bm \times t$ | $bm \times D$ | $bm \times D$ | $bm \times D$ |
| **Target** | Total improvement rate | Drug-specific effect | Drug-specific effect | Drug-specific effect |
| **Confounded with** | ER, TV slopes | -- | -- | -- |
| **ANCOVA info** | $N \cdot \text{Var}(bm) \cdot \text{Var}(b_i)^{-1}$ | $N \cdot \text{Var}(bm) \cdot \sigma^{-2}_\Delta$ | $N \cdot \text{Var}(bm) \cdot \sigma^{-2}_\Delta$ | $N \cdot \text{Var}(bm) \cdot \sigma^{-2}_\Delta$ |
| **LMM info** | $\sum \widetilde{bm}^2_i \cdot \text{SS}(t_i)$ | $\sum \widetilde{bm}^2_i \cdot 1.5$ | $\sum \widetilde{bm}^2_i \cdot 2.0$ | $\sum \widetilde{bm}^2_i \cdot 2.9$ |
| **GEE info** | ~Same as LMM | ~Same as LMM | ~Same as LMM | ~Same as LMM |
| **RM-ANOVA info** | ~$0.64 \times$ ANCOVA | ~$0.64 \times$ ANCOVA | ~$0.64 \times$ ANCOVA | ~$0.64 \times$ ANCOVA |
| **$K$ (transitions)** | 0 | 1 | 1 | 3 |
| **Autocorrelation vulnerability** | Low (no transitions) | Medium (1 long block + 1 short block) | Medium (2 long blocks) | Low per block (shorter blocks) but more boundary effects |
| **Carryover vulnerability** | None (no transitions) | Low (1 transition) | Low (1 transition) | High (3 transitions) |

The RM-ANOVA row reflects the ~$2/\pi$ efficiency penalty from
dichotomizing the biomarker at the median. If timepoints within phases
are retained as separate within-subject levels rather than collapsed to
phase means, GG/HF corrections further reduce effective degrees of
freedom under autocorrelation, widening the gap relative to the
regression-based methods.

The N-of-1 design achieves the highest information for the biomarker
interaction under the LMM and GEE approaches when carryover is absent,
owing to both its larger total $T$ and its near-balanced drug exposure.
The crossover achieves optimal balance ($\bar{D} = 0.5$) with fewer
transitions, making it more robust to carryover at the cost of fewer
total timepoints. The OL+BDC design has the lowest information among
the designs with within-subject drug variation, because its short
discontinuation block contributes few off-drug observations.

---

### 3. Trial Designs (NON-SIMPLIFIABLE)

**CRITICAL**: The four Hendrickson trial designs must be preserved exactly as
specified. The design structure (measurement weeks, randomization paths,
blinding, expectancy effects) is the core scientific comparison and cannot be
simplified.

#### The Four Designs (from Hendrickson et al. 2020)

| Design | Weeks | Structure | Blinding |
|--------|-------|-----------|----------|
| Open-Label (OL) | 1-8 | All on drug, unblinded | None |
| OL + Blinded Discontinuation (OL+BDC) | 1-4 OL, 5-8 randomized drug/placebo | Partial |
| Traditional Crossover (CO) | 4 weeks drug, 4 weeks placebo (order randomized) | Full |
| N-of-1 (Hybrid) | 1-4 OL, 5-12 two 4-week crossover cycles | Partial |

#### Implementation Requirements

Each design must specify:
1. **Week-by-week on_drug indicator** (TRUE/FALSE per measurement week)
2. **Blinding status** (affects expectancy response)
3. **Randomization path** (which weeks are randomized)

```r
# Exact design specifications from Hendrickson
designs <- list(
  ol = list(
    weeks = 1:8,
    on_drug = rep(TRUE, 8),
    blinded = rep(FALSE, 8),
    expectancy = rep(TRUE, 8)  # Unblinded: full expectancy

),
  ol_bdc = list(
    weeks = 1:8,
    on_drug = c(rep(TRUE, 4), NA, NA, NA, NA),  # NA = randomized
    blinded = c(rep(FALSE, 4), rep(TRUE, 4)),
    expectancy = c(rep(TRUE, 4), rep(FALSE, 4))  # Blinded: no expectancy
  ),
  crossover = list(
    weeks = 1:8,
    on_drug = NA,  # Fully randomized order
    blinded = rep(TRUE, 8),
    expectancy = rep(FALSE, 8)
  ),
  nof1 = list(
    weeks = 1:12,
    on_drug = c(rep(TRUE, 4), NA, NA, NA, NA, NA, NA, NA, NA),
    blinded = c(rep(FALSE, 4), rep(TRUE, 8)),
    expectancy = c(rep(TRUE, 4), rep(FALSE, 8))
  )
)
```

**Note**: The `NA` values indicate randomized assignment at simulation time.

### 4. Carryover Effects

#### 4a. Carryover Model

**Current**: Exponential decay with half-life parameter.

**Simplified options**:

| Option | Formula | Description |
|--------|---------|-------------|
| No carryover | `carryover = 0` | Assume instant washout |
| Binary carryover | `carryover = ifelse(time_since < washout, 1, 0)` | All-or-nothing |
| Linear decay | `carryover = max(0, 1 - time_since/washout)` | Linear decrease |
| Step function | `carryover = ifelse(time_since < 1, 0.5, 0)` | Partial for one period |

**Recommendation**: Start with **no carryover**, add **linear decay** as
sensitivity analysis.

```r
# Linear carryover
carryover_linear <- function(time_since_drug, washout_period) {
  pmax(0, 1 - time_since_drug / washout_period)
}
```

### 5. Censoring/Dropout

#### 5a. Dropout Model

**Current**: Five patterns (none, balanced, high dropout, more biased, more
flat) with complex probability functions.

**Simplified options**:

| Option | Description |
|--------|-------------|
| No dropout | Complete data for all |
| MCAR | Random dropout, constant probability |
| Single MAR | One dropout pattern tied to response |

**Recommendation**: Start with **no dropout**. Add MCAR (e.g., 10% per
timepoint) as sensitivity.

### 6. Parameter Space

#### 6a. Reduce Dimensions

**Current**: Multiple sensitivity analyses varying different parameter sets.

**Simplified approach**:

| Parameter | Current | Simplified |
|-----------|---------|------------|
| BR_rate | 0.05, 0.15, 0.3 | 0.1, 0.3 (2 levels) |
| ER_rate | 0.05, 0.15, 0.3 | 0.1 (fixed) |
| TR_rate | 0.05, 0.15, 0.3 | 0 (removed) |
| biomarker_mod | 0, 0.3, 0.6 | 0.3 (fixed for power), vary for sensitivity |
| carryover | 0, 0.1, 0.2 | 0 (default), 0.2 (sensitivity) |
| censoring | 5 levels | none (default) |
| N | 35, 70 | 35, 50, 70 |

**Recommendation**: Single main simulation with key parameters, separate
one-at-a-time sensitivity analyses.

## Implementation Strategy

### Phase 1: Minimal Viable Simulation

1. Single script `simulation.R` (~200 lines)
2. Random intercept model for data generation
3. ANCOVA on phase means for analysis
4. Four designs with simplified specification
5. No carryover, no dropout
6. Fixed BR_rate, ER_rate; vary biomarker_mod and N

### Phase 2: Add Complexity as Needed

1. Add carryover (linear decay)
2. Add dropout (MCAR)
3. Sensitivity on BR_rate
4. Compare ANCOVA vs LMM results

### Phase 3: Validation

1. Compare power estimates to pmsimstats2025
2. Document where simplifications matter vs don't

## Code Architecture

```
01b-pmsimstats-simple/
├── R/
│   └── functions.R          # Shared utility functions
├── analysis/
│   ├── scripts/
│   │   ├── simulation.R     # Main simulation (single file)
│   │   └── sensitivity.R    # Parameter sensitivity (optional)
│   └── output/
└── docs/
    └── simplification-plan.md
```

### Key Design Principles

1. **Single file first**: One simulation.R that does everything
2. **Explicit over clever**: Write out loops rather than complex abstractions
3. **Base R preferred**: Minimize dependencies (tidyverse for data manipulation only)
4. **No parallel by default**: Add parallelization only if needed
5. **Inline parameters**: Define parameters at top of script, not in config files

## Expected Outcomes

| Metric | pmsimstats2025 | pmsimstats-simple (target) |
|--------|----------------|---------------------------|
| Lines of code | ~700 per script | ~200 total |
| Dependencies | nlme, lme4, MASS, zoo, patchwork, future, furrr | tidyverse only |
| Runtime (20 iter) | ~160 sec | ~10-20 sec |
| Convergence issues | Occasional lme failures | None (closed-form) |
| Parameters | ~15 varying | ~5 varying |

## Literature Support for Simplifications

### Key Reference: Hendrickson et al. (2020)

The primary reference for this work is [Hendrickson et al. (2020)](https://doi.org/10.3389/fdgth.2020.00013)
"Optimizing Aggregated N-Of-1 Trial Designs for Predictive Biomarker Validation"
in *Frontiers in Digital Health*.

**Key methodological points from Hendrickson:**

1. **Three-factor response model**: BR (biologic response), ER (expectancy
   response), TR (time-dependent response). Each modeled via Gompertz function
   with max, rate, displacement parameters.

2. **Analysis approach**: Linear mixed effects model (lme4) with random intercept
   per subject. For designs with on/off drug: `response ~ on_drug + biomarker +
   on_drug:biomarker + (1|id)`. For open-label only: `response ~ time + biomarker
   + time:biomarker + (1|id)`.

3. **Carryover**: Exponential decay with half-life parameter. Critical finding:
   even short (0.1 week) carryover causes "precipitous decline in power" for
   N-of-1 and OL+BDC designs. Crossover design more robust.

4. **Censoring/dropout**: Probability = β₀ + β₁×(change in symptoms)². Biased
   dropout where responders more likely to drop out during placebo periods.

5. **Key result**: N-of-1 design has power between OL+BDC (lower) and traditional
   crossover (higher), but allows all participants to start on open-label
   treatment.

### Analysis Method Comparisons

**Chen & Chen (2014)** - [A Comparison of Four Methods for the Analysis of N-of-1
Trials](https://doi.org/10.1371/journal.pone.0087752):

Compared paired t-test, mixed effects model, mixed effects model of difference,
and meta-analysis for N-of-1 trials (3-4 cycles, N=1-30).

| Method | Type I Error | Power | Bias | MSE |
|--------|-------------|-------|------|-----|
| Paired t-test | Near nominal | Highest | Comparable | Smallest |
| Mixed effects | Similar | Lower | Smaller | Bigger |
| ME of difference | Far from nominal | Low | Large | Large |
| Meta-analysis | Far from nominal | Low | Large | Large |

**Recommendation**: "Paired t-test was recommended to use for normally
distributed data of N-of-1 trials irrespective with or without carryover
effect."

**Critique (Araujo et al. 2016)**: Chen & Chen models "do not include a treatment
by patient interaction" advocated in meta-analysis literature.

### ANCOVA for RCTs

**BMC Medical Research Methodology (2014)** - [Bias, precision and statistical
power of ANCOVA](https://doi.org/10.1186/1471-2288-14-49):

"ANCOVA should be the analysis of choice, a priori, for RCTs with a single
post-treatment outcome measure previously measured at baseline; its superiority
is particularly marked when baseline imbalance is present, but also—in terms of
precision—when groups are balanced at baseline."

### Biomarker-Treatment Interaction Power

**Trials (2018)** - [A simulation study on estimating biomarker–treatment
interaction effects](https://doi.org/10.1186/s13063-018-2491-0):

- "Tests of interactions are often lacking statistical power"
- "The probability of detecting a biomarker–treatment interaction can be
  increased by including prognostic variables"
- "If the investigation of a biomarker–treatment interaction is of major
  importance for a clinical trial, this should be considered in the design stage"

### N-of-1 Power Analysis

**Wang & Schork (2019)** - [Power and Design Issues in Crossover-Based N-Of-1
Clinical Trials](https://doi.org/10.3390/healthcare7030084):

Analytical and simulation-based power analysis for N-of-1 designs considering:
- Number of cycles
- Serial correlation between measurements
- Washout periods
- Heteroscedasticity
- Carry-over phenomena

### Implications for Simplification

| Literature Finding | Implication for pmsimstats-simple |
|--------------------|-----------------------------------|
| Paired t-test has highest power (Chen 2014) | ANCOVA on change scores is reasonable |
| ANCOVA superior for baseline-adjusted outcomes (BMC 2014) | Supports phase-mean approach |
| Carryover devastating to N-of-1 power (Hendrickson 2020) | Must include carryover sensitivity |
| LME may have lower power than simpler methods | Justifies avoiding complex mixed models |
| Interaction tests lack power (Trials 2018) | Need adequate sample sizes |

## Decision Log

| Decision | Rationale | Date |
|----------|-----------|------|
| Use ANCOVA | Chen 2014 shows paired t-test ≈ ANCOVA has highest power; BMC 2014 supports ANCOVA | 2026-02-14 |
| Keep 4 designs | Core scientific comparison per Hendrickson 2020 | 2026-02-14 |
| Keep ER component | Hendrickson shows expectancy (blinded vs open) affects power | 2026-02-14 |
| Include carryover sensitivity | Hendrickson shows carryover "precipitous" effect on power | 2026-02-14 |
| No carryover default | Simplifies initial analysis; add as sensitivity | 2026-02-14 |
| Random intercept + iid | Chen 2014 shows simpler models perform well | 2026-02-14 |

## References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA (2020). Optimizing Aggregated
   N-Of-1 Trial Designs for Predictive Biomarker Validation: Statistical Methods
   and Theoretical Findings. *Front. Digit. Health* 2:13.
   https://doi.org/10.3389/fdgth.2020.00013

2. Chen X, Chen P (2014). A Comparison of Four Methods for the Analysis of N-of-1
   Trials. *PLoS ONE* 9(2): e87752.
   https://doi.org/10.1371/journal.pone.0087752

3. Araujo A, Julious S, Senn S (2016). Understanding Variation in Sets of N-of-1
   Trials. *PLoS ONE* 11(12): e0167167.
   https://doi.org/10.1371/journal.pone.0167167

4. Egbewale BE, Lewis M, Sim J (2014). Bias, precision and statistical power of
   analysis of covariance in the analysis of randomized trials with baseline
   imbalance. *BMC Med Res Methodol* 14:49.
   https://doi.org/10.1186/1471-2288-14-49

5. Eng KH, Konings P (2018). A simulation study on estimating biomarker–treatment
   interaction effects in randomized trials with prognostic variables. *Trials*
   19:72.
   https://doi.org/10.1186/s13063-018-2491-0

6. Wang Y, Schork NJ (2019). Power and Design Issues in Crossover-Based N-Of-1
   Clinical Trials with Fixed Data Collection Periods. *Healthcare* 7(3):84.
   https://doi.org/10.3390/healthcare7030084

7. Kravitz R, Duan N, editors (2014). Design and Implementation of N-of-1 Trials:
   A User's Guide. AHRQ Publication No. 13(14)-EHC122-EF.
   https://effectivehealthcare.ahrq.gov/products/n-1-trials/research-2014-5
