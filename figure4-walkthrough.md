# Figure 4 Walkthrough: Power Heatmaps from Hendrickson et al. (2020)

This document provides a complete technical walkthrough of the power
heatmaps in Figure 4 of Hendrickson RC, Thomas RG, Schork NJ, and
Raskind MA (2020), 'Optimizing Aggregated N-Of-1 Trial Designs for
Predictive Biomarker Validation,' *Frontiers in Digital Health*
2:13. It covers (1) the methodology behind the simulations, (2) the
code that implements them in the `pmsimstats` package, and (3) how
to reproduce the figures from the bundled data.

---

## 1. Scientific Question

Figure 4 addresses a central question in precision medicine trial
design: *What is the statistical power to detect a
biomarker-treatment interaction across four candidate trial
designs, and how does that power depend on sample size, the
strength of the biomarker-response correlation, carryover effects,
and participant dropout patterns?*

The biomarker is standing systolic blood pressure; the outcome is
PTSD symptom severity (CAPS score); the treatment is the alpha-1
adrenoceptor antagonist prazosin. The interaction of interest is
whether baseline blood pressure predicts individual treatment
response -- not whether the drug works on average, but whether the
biomarker identifies who will respond.

---

## 2. The Four Trial Designs

Each design spans 20 weeks with 8 assessment timepoints (excluding
baseline). They differ in blinding, randomization structure, and
the number of distinct paths a participant can follow through the
trial.

### 2a. Open-Label (OL) -- 1 path

All participants receive active drug for 20 weeks, unblinded.
Expectancy is 1.0 at all timepoints. There is no within-subject
drug variation, so the analysis tests whether the biomarker
predicts the *rate* of improvement over time (`bm:t`
interaction).

```r
tdOL <- buildtrialdesign(
  name_longform = "open label",
  name_shortform = "OL",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste("OL", 1:8, sep = ""),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)
```

### 2b. Open-Label + Blinded Discontinuation (OL+BDC) -- 2 paths

16 weeks open-label on drug, then a 4-week blinded
discontinuation block. In the blinded block, participants are
randomized to receive active drug for the first 1 or 2 weeks
before switching to placebo. Expectancy drops to 0.5 during the
blinded phase.

```r
tdOLBDC <- buildtrialdesign(
  name_longform = "open label+blinded discontinuation",
  name_shortform = "OL+BDC",
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptname = c("OL1","OL2","OL3","OL4",
                 "BD1","BD2","BD3","BD4"),
  expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)
```

### 2c. Traditional Crossover (CO) -- 2 paths

Two 10-week blocks: one on drug, one on placebo, order
randomized. All blinded (expectancy = 0.5). No washout period.

```r
tdCO <- buildtrialdesign(
  name_longform = "traditional crossover",
  name_shortform = "CO",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = c(paste("COa", 1:4, sep = ""),
                 paste("COb", 1:4, sep = "")),
  expectancies = rep(.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)
```

### 2d. N-of-1 (Hybrid) -- 4 paths

8-week open-label phase followed by a 4-week blinded
discontinuation and two 4-week crossover blocks. Two independent
randomization points create 4 possible paths. The timepoints are
non-uniform: c(4, 8, 9, 10, 11, 12, 16, 20), with weekly
assessments clustered around the blinded-discontinuation
transition.

```r
tdNof1 <- buildtrialdesign(
  name_longform = "primary N-of-1 design",
  name_shortform = "N-of-1",
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptname = c("OL1","OL2","BD1","BD2",
                 "BD3","BD4","COd","COp"),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)
```

The `buildtrialdesign()` function computes three derived timing
variables from the `ondrug` vector:

- **tod** (time on drug): cumulative weeks on active treatment,
  resetting to 0 when off drug
- **tsd** (time since discontinuation): cumulative weeks off drug
  since last on-drug period, zeroed before the participant is
  ever on drug
- **tpb** (time in positive belief): cumulative weeks where
  expectancy > 0

These drive the Gompertz trajectory calculations in the DGP.

---

## 3. The Data-Generating Process

### 3a. Response Decomposition

The observed outcome (symptom score) at each timepoint t for
participant i is:

```
Symptom_{i,t} = Baseline_i - [TV_{i,t} + ER_{i,t} + BR_{i,t}]
```

Higher factor values represent greater improvement, so they are
subtracted from baseline. The three factors are:

- **TV (time-variant response, TR in the paper):** Natural history
  unrelated to treatment. Includes regression to the mean, the
  effect of study participation, and time-dependent fluctuations.
- **ER (expectancy-related response):** Response attributable to
  the participant's belief that they are receiving active
  treatment. Scaled by the expectancy value (1.0 for open-label,
  0.5 for blinded, 0 when expectancy is absent).
- **BR (biologic response):** Direct pharmacologic effect of the
  drug. This is the only factor that depends on actual drug
  exposure.

### 3b. Modified Gompertz Function

Each factor's expected trajectory is modeled by a modified
Gompertz curve:

```r
modgompertz <- function(t, maxr, disp, rate) {
  y <- maxr * exp(-disp * exp(-rate * t))
  vert_offset <- maxr * exp(-disp * exp(-rate * 0))
  y <- y - vert_offset
  y <- y * (maxr / (maxr - vert_offset))
  y
}
```

This function guarantees f(0) = 0 and f(infinity) = maxr. The
three parameters control:

- `maxr`: asymptotic maximum response
- `disp`: displacement (controls the inflection point)
- `rate`: growth rate

The 'tabula rasa' (TR) parameters extracted from pilot data:

| Factor | max | disp | rate | sd |
|--------|------:|-----:|-----:|---:|
| TV | 6.506 | 5 | 0.35 | 10 |
| ER | 6.506 | 5 | 0.35 | 10 |
| BR | 10.986 | 5 | 0.42 | 8 |

TV and ER share parameters because pilot data did not allow
separating them (only the combined placebo response was
observed).

### 3c. Mean Vector Construction

For a single path with nP = 8 timepoints, the mean vector has
2 + 3*8 = 26 elements:

```
means = [bm_mean, BL_mean,
         TV_1 ... TV_8,
         ER_1 ... ER_8,
         BR_1 ... BR_8]
```

Where:

```r
# Baseline (extracted from pilot data):
bm_mean = 124.328   # systolic BP
BL_mean = 83.069    # CAPS score

# Factor means at each timepoint:
TV_t = modgompertz(t_cumulative, 6.506, 5, 0.35)
ER_t = modgompertz(tpb, 6.506, 5, 0.35) * expectancy_t
BR_t = modgompertz(tod, 10.986, 5, 0.42) + carryover_t
```

### 3d. Carryover Effect

When a participant transitions from on-drug to off-drug, the BR
mean does not immediately drop to zero. Instead, the carryover
from the previous timepoint decays exponentially:

```r
if (!onDrug[t] && tsd[t] > 0) {
  BR_t = BR_t + BR_{t-1} *
         (1/2)^(scalefactor * tsd[t] / carryover_t1half)
}
```

Key features:

- **Sequential propagation:** Carryover at t depends on the
  carryover-adjusted value at t-1, not the raw value. This
  allows residual effects to chain across multiple off-drug
  periods.
- **Scale factor = 2** (default): The drug effect decays faster
  than the raw half-life would suggest, reflecting the combined
  effect of pharmacokinetic clearance and clinical re-emergence
  of symptoms.
- **Half-life values tested:** 0, 0.1, 0.2 weeks. Even 0.1 weeks
  (~17 hours) produces substantial power loss.
- **BR only:** Carryover is not applied to TV or ER.

### 3e. Covariance Matrix

The 26x26 covariance matrix encodes correlations between all
variables. It is constructed from a correlation matrix and a
vector of standard deviations:

```r
sigma <- outer(sds, sds) * correlations
```

**Standard deviations:**

```
sds = [bm_sd, BL_sd,
       rep(TV_sd, 8),
       rep(ER_sd, 8) * expectancy,
       rep(BR_sd, 8)]
```

Where `bm_sd = 15.362`, `BL_sd = 18.483`, `TV_sd = 10`,
`ER_sd = 10`, `BR_sd = 8`. The ER standard deviations are
scaled by the expectancy value at each timepoint (0 expectancy
means 0 variance in ER).

**Correlation structure (fixed across all simulations):**

| Correlation | Value | Meaning |
|-------------|------:|---------|
| c.tv, c.pb, c.br | 0.8 | Within-factor autocorrelation (compound symmetry) |
| c.cf1t | 0.2 | Cross-factor correlation at the same timepoint |
| c.cfct | 0.1 | Cross-factor correlation at different timepoints |

**Biomarker-BR correlation (the critical interaction mechanism):**

The correlation between the biomarker (`bm`) and the bio-response
(`BR_t`) at each timepoint determines the biomarker-treatment
interaction. This is set according to Ron Thomas's logic:

- **Timepoint 1:** Always 0 (no drug effect at first measurement)
- **On-drug timepoints (raw BR > 0):** `Cor(bm, BR_t) = c.bm`
- **Off-drug, no carryover residual (BR = 0):** `Cor(bm, BR_t) = 0`
- **Off-drug, with carryover residual (raw BR was 0 but adjusted
  BR > 0):** `Cor(bm, BR_t) = (BR_t / BR_{t-1}) * c.bm`
  (scaled proportionally to the carryover fraction)

This differential correlation IS the biomarker-treatment
interaction. There is no population-level mean shift; the
interaction emerges from the fact that participants with higher
biomarker values have BR draws that are more correlated (when on
drug) or less correlated (when off drug) with the biomarker.

**Positive definiteness:** After construction, the covariance
matrix is checked via `corpcor::is.positive.definite()`. If it
fails, it is forced positive definite via
`corpcor::make.positive.definite(sigma, tol = 1e-3)`.

### 3f. Multivariate Normal Draw

With the mean vector and covariance matrix in hand, N
participants are drawn from a multivariate normal distribution:

```r
dat <- MASS::mvrnorm(n = N, mu = means, Sigma = sigma)
```

This produces a data.table where each row is a participant and
the columns are: biomarker, baseline, and 24 factor values
(3 factors x 8 timepoints). The outcome at each timepoint is
then computed as:

```
Symptom_t = BL - (TV_t + ER_t + BR_t)
```

This is done via `eval(parse(text = ...))` in `generateData.R`
lines 200-215, constructing column expressions dynamically.

---

## 4. Censoring (Dropout) Model

The censoring model simulates realistic participant dropout
with two components:

1. **Time-dependent (flat) dropout:** Probability increases
   linearly with the duration of each interval. Controlled
   by `beta0 * (1 - beta1)`.

2. **Outcome-dependent (biased) dropout:** Probability is
   higher for participants whose symptoms have not improved
   (or have worsened). Computed as `(symptom_change + 100)^eb1`,
   where the shift by 100 ensures all values are positive.
   Controlled by `beta0 * beta1`.

Once a timepoint is censored, all subsequent timepoints for that
participant are also set to NA (last-observation-carried-forward
censoring).

**Censoring parameter sets used in Figure 4:**

| Name | beta0 | beta1 | eb1 | Description |
|------|------:|------:|----:|-------------|
| No censoring | -- | -- | -- | Complete data |
| balanced | 0.05 | 0.5 | 2 | Equal flat and biased dropout |
| more of flat | 0.05 | 0.2 | 2 | Mostly time-dependent dropout |
| more of biased | 0.05 | 0.8 | 2 | Mostly outcome-dependent dropout |
| high dropout | 0.15 | 0.5 | 2 | Triple the overall dropout rate |

---

## 5. The Analysis Model

Each simulated trial is analyzed using a linear mixed-effects
model (`lme4::lmer`). The model formula is selected dynamically
based on the trial design structure.

### 5a. Designs with within-subject drug variation (CO, N-of-1, OL+BDC)

The formula is:

```
Sx ~ bm + t + Dbc + bm:Dbc + (1|ptID)
```

Where:

- `Sx`: symptom score at each timepoint
- `bm`: baseline biomarker value (continuous)
- `t`: time since baseline (continuous)
- `Dbc`: continuous drug indicator. When on drug: `Dbc = 1`.
  When off drug with carryover modeling:
  `Dbc = (1/2)^(scalefactor * tsd / carryover_t1half)`.
  When off drug without carryover: `Dbc = 0`.
- `bm:Dbc`: the biomarker-treatment interaction (primary
  outcome of interest)
- `(1|ptID)`: random intercept per participant

The coefficient of `bm:Dbc` represents the extent to which the
biomarker predicts differential treatment response. A significant
negative coefficient indicates that higher biomarker values
predict greater symptom improvement when on drug.

### 5b. Designs without within-subject drug variation (OL)

When all participants are always on drug, there is no variation
in `Dbc`, so the model cannot estimate `bm:Dbc`. Instead:

```
Sx ~ bm + t + bm:t + (1|ptID)
```

The `bm:t` coefficient tests whether the biomarker predicts the
rate of improvement over time.

### 5c. Model construction in code

The formula is built as a string and evaluated:

```r
# lme_analysis.R lines 133-151
modelbase <- "Sx~bm+t"
if (varInDb) {
  modelbase <- paste0(modelbase, "+Dbc+bm*Dbc")
} else {
  modelbase <- paste0(modelbase, "+bm*t")
}
# Optional terms for expectancy, carryover, random slopes...
modelbase <- paste0(modelbase, "+(1|ptID)")
eval(parse(text = paste0("form<-", modelbase)))
fit <- lmer(form, data = datamerged)
```

### 5d. Output extraction

The p-value for the interaction term is extracted from the
`lmerTest` summary (which uses Satterthwaite degrees of freedom):

```r
c <- summary(fit)$coefficients
if (varInDb) {
  p     <- c["bm:Dbc", "Pr(>|t|)"]
  beta  <- c["bm:Dbc", "Estimate"]
  betaSE <- c["bm:Dbc", "Std. Error"]
} else {
  p     <- c["bm:t", "Pr(>|t|)"]
  beta  <- c["bm:t", "Estimate"]
  betaSE <- c["bm:t", "Std. Error"]
}
```

---

## 6. The Simulation Loop

`generateSimulatedResults()` orchestrates the parameter sweep:

1. **Expand the parameter grid:**
   `expand.grid(trialdesign, respparamset, blparamset,
   modelparamset)` -- note that censoring is applied after
   data generation, not included in the grid.

2. **For each parameter combination:**
   - Generate `N * Nreps` participants at once (for efficiency)
     across all paths
   - Split into `Nreps` replicates
   - Fit the LME model to the full population to get the
     'empirical true beta' (ETbeta)
   - Fit the LME model to each replicate to get per-replicate
     beta, betaSE, p-value
   - Optionally apply each censoring pattern and re-analyze

3. **Output:** A data.table with one row per
   (parameter-combination x replicate x censoring-pattern),
   containing: beta, betaSE, p, issingular, warning, ETbeta,
   ETbetaSE, frac_NA.

### 6a. Model parameters for Figure 4

```r
coremodelparams <- expand.grid(
  N = c(35, 70),
  c.bm = c(0, 0.3, 0.6),
  carryover_t1half = c(0, 0.1, 0.2),
  c.tv = 0.8, c.pb = 0.8, c.br = 0.8,
  c.cf1t = 0.2, c.cfct = 0.1
)
```

This gives 2 x 3 x 3 = 18 model parameter combinations, times
4 trial designs, times 1 response parameter set, times 1
baseline parameter set = 72 cells in the design grid. Each cell
is replicated 1,000 times. Each replicate is analyzed with and
without each of 4 censoring patterns, yielding 5 analysis
conditions per replicate.

Analysis parameters: `useDE = FALSE` (expectancy not included
as a fixed effect), `t_random_slope = FALSE` (random intercept
only), `full_model_out = FALSE`.

---

## 7. Figure 4 Heatmap Construction

### 7a. Figure 4A: Effect of trial design, censoring, and
biomarker effect on power

```r
param2vary <- c("trialdesign", "N", "c.bm",
                "censorparamset")
param2hold <- data.table(blparamset = 1,
                         respparamset = 1,
                         carryover_t1half = 0)
param2nothold <- c("c.cfct", "modelparamset")
param2plot <- "power"
```

| Mapping | Parameter | Values |
|---------|-----------|--------|
| Facet columns | trialdesign | OL, OL+BDC, CO, N-of-1 |
| Facet rows (right axis) | N | 35, 70 |
| X-axis | c.bm | 0, 0.3, 0.6 |
| Y-axis | censorparamset | No censoring, balanced, more of flat, more of biased, high dropout |
| Fill color | power | 0.0 to 1.0 |
| Held constant | carryover_t1half | 0 |

### 7b. Figure 4B: Effect of nonzero carryover term on power

```r
param2vary <- c("trialdesign", "N",
                "carryover_t1half",
                "censorparamset")
param2hold <- data.table(blparamset = 1,
                         respparamset = 1,
                         c.bm = 0.6)
param2nothold <- c("c.cfct", "modelparamset")
param2plot <- "power"
```

| Mapping | Parameter | Values |
|---------|-----------|--------|
| Facet columns | trialdesign | OL, OL+BDC, CO, N-of-1 |
| Facet rows (right axis) | N | 35, 70 |
| X-axis | carryover_t1half | 0, 0.1, 0.2 |
| Y-axis | censorparamset | No censoring, balanced, more of flat, more of biased, high dropout |
| Fill color | power | 0.0 to 1.0 |
| Held constant | c.bm | 0.6 |

### 7c. How `PlotModelingResults()` computes power

The function (plottingfunctions.R:26-191) operates as follows:

1. **Filter** `data$results` to the held-constant parameters.
   For parameters not explicitly held or varied, the function
   takes the first available value (line 92-94).

2. **Convert facet parameters** (first two in `param2vary`) to
   factors with human-readable labels from the trial design
   metadata and censoring parameter names.

3. **Compute summary statistics** grouped by the four varying
   parameters:

   ```r
   d[, sig_p := (p < 0.05), by = tID]
   d[, power := mean(sig_p), by = op]
   ```

   Power is the fraction of replicates where the interaction
   term is significant at alpha = 0.05.

4. **Plot** using `ggplot2`:

   ```r
   ggplot(data = ds, aes_string(op[3], op[4])) +
     geom_tile(aes(fill = toplot), colour = "white") +
     scale_fill_distiller(type = "seq", palette = "RdYlBu") +
     geom_text(aes(label = round(toplot, 2)), size = 3) +
     facet_grid(reformulate(op[1], op[2]))
   ```

   The RdYlBu palette maps low power (red) through yellow to
   high power (blue).

---

## 8. Reproduction Script

The bundled `results_core.rda` dataset contains the pre-computed
simulation results (1,000 replicates). To reproduce Figure 4:

```r
library(pmsimstats)
library(data.table)
library(ggplot2)
library(gridExtra)

data(results_core)
simresults <- results_core

# --- Figure 4A ---
p1 <- PlotModelingResults(
  simresults,
  param2plot = "power",
  param2vary = c("trialdesign", "N", "c.bm",
                 "censorparamset"),
  param2hold = data.table(blparamset = 1,
                          respparamset = 1,
                          carryover_t1half = 0),
  param2nothold = c("c.cfct", "modelparamset")
)

# --- Figure 4B ---
p2 <- PlotModelingResults(
  simresults,
  param2plot = "power",
  param2vary = c("trialdesign", "N",
                 "carryover_t1half",
                 "censorparamset"),
  param2hold = data.table(blparamset = 1,
                          respparamset = 1,
                          c.bm = 0.6),
  param2nothold = c("c.cfct", "modelparamset")
)

# --- Compose ---
p1 <- p1 +
  ggtitle(paste("A. Effect of trial design, censoring",
                "and biomarker effect on power")) +
  theme(plot.title = element_text(hjust = 0))

p2 <- p2 +
  ggtitle(paste("B. Effect of nonzero carryover term",
                "on power")) +
  theme(plot.title = element_text(hjust = 0))

grid.arrange(p1, p2, nrow = 2)
```

To re-run the simulations from scratch (rather than using the
bundled data), set `Nreps <- 1000` in the first vignette and
run the `coremodelparams` block. On a modern machine, expect
approximately 8-12 hours for the full parameter space.

---

## 9. Code Map

| Step | File | Lines | Function |
|------|------|------:|----------|
| Define trial designs | buildtrialdesign.R | 75-119 | `buildtrialdesign()` |
| Compute tod/tsd/tpb | buildtrialdesign.R | 88-106 | (within `buildtrialdesign()`) |
| Gompertz trajectories | utilities.R | 34-40 | `modgompertz()` |
| Cumulative time helper | utilities.R | 8-18 | `cumulative()` |
| Build mean vector | generateData.R | 63-103 | (within `generateData()`) |
| Carryover on BR means | generateData.R | 86-100 | (within `generateData()`) |
| Build correlation matrix | generateData.R | 105-163 | (within `generateData()`) |
| BM-BR correlation logic | generateData.R | 146-162 | Ron Thomas scaling |
| Covariance from correlation | generateData.R | 172 | `outer(sds, sds) * correlations` |
| PD enforcement | generateData.R | 177-181 | `corpcor::make.positive.definite()` |
| MVN draw | generateData.R | 188 | `MASS::mvrnorm()` |
| Compute outcomes | generateData.R | 200-215 | `eval(parse())` column construction |
| Censoring/dropout | censordata.R | 39-89 | `censordata()` |
| Parameter sweep loop | generateSimulatedResults.R | 66-255 | `generateSimulatedResults()` |
| LME model construction | lme_analysis.R | 133-151 | Dynamic formula via `eval(parse())` |
| Interaction extraction | lme_analysis.R | 196-205 | `bm:Dbc` or `bm:t` coefficient |
| Heatmap visualization | plottingfunctions.R | 26-191 | `PlotModelingResults()` |
| Trajectory visualization | plottingfunctions.R | 220-338 | `plotfactortrajectories()` |

---

## 10. Key Findings from Figure 4

### Figure 4A (no carryover)

- **N-of-1 has the highest power** across all conditions, reaching
  0.95 at N=35, c.bm=0.6 without censoring.
- **OL has the lowest power** because it lacks within-subject drug
  variation and can only detect the biomarker-time interaction.
- **CO has slightly lower power than N-of-1** despite full blinding,
  because the N-of-1 design includes both open-label and
  crossover information.
- **Censoring reduces power** across all designs, with the strongest
  effect on OL+BDC (which has the shortest blinded phase,
  making it most vulnerable to pre-blinding dropout).
- **c.bm = 0 produces power near 0.05** across all designs,
  confirming appropriate Type I error control.

### Figure 4B (carryover at c.bm = 0.6)

- **Even short carryover (0.1 weeks) precipitously reduces power**
  for N-of-1 and OL+BDC designs. N-of-1 drops from 0.95 to
  ~0.18 at N=35.
- **CO is least affected** because it has only one transition from
  on-drug to off-drug (or vice versa), while N-of-1 has multiple
  transitions amplifying carryover contamination.
- **OL is unaffected** by carryover because participants are always
  on drug.
- The vulnerability of the N-of-1 design to carryover is the
  paper's most clinically important finding: the design that
  offers the best power under ideal conditions is also the most
  sensitive to even minimal carryover effects.
