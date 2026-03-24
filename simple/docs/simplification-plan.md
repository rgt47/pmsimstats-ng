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
