# Response to Collaborator Goals

Date: 2026-03-23 (Updated 2026-04-08)
Context: Following the pmsimstats code audit and revision

**Note**: The code implements **Architecture A (direct mean moderation)** for
biomarker-treatment interactions. See `docs/02-dgp-mean-moderation-vs-mvn.md`
for detailed comparison with differential correlation approaches and implications
for power under carryover.

## Code Locations

All new functions are in the `orig` repo at
`~/Dropbox/prj/alz/01c-pmsimstats-orig/pmsimstats/`:

| Function | File | Purpose |
|----------|------|---------|
| `analyze_trial_extended()` | `R/carryover_analysis.R` | Main drug effect + interaction + predicted curves |
| `characterize_carryover()` | `R/carryover_analysis.R` | Sweep half-lives, compare AIC, find best fit |
| `print_trial_summary()` | `R/carryover_analysis.R` | Clinical summary of trial analysis |
| `print_carryover_summary()` | `R/carryover_analysis.R` | Summary of carryover characterization |
| `lme_analysis()` | `R/lme_analysis.R` | Core LME analysis (nlme + corCAR1) |
| `generateData()` | `R/generateData.R` | Data generation with cached sigma |
| `buildSigma()` | `R/generateData.R` | Covariance matrix construction |
| `validateParameterGrid()` | `R/generateData.R` | PD pre-validation |
| `generateSimulatedResults()` | `R/generateSimulatedResults.R` | Full simulation loop with parallel processing |

---

## Goal 1: Did prazosin improve PTSD symptoms?

*How to analyze this in this trial design? -- RCH implement
with RGT guidance*

### What we have

The revised `lme_analysis()` uses `nlme::lme` with `corCAR1`
(continuous-time AR(1) residual correlation), which is
directly applicable to actual trial data. Vignette 3 in the
original package demonstrates applying `lme_analysis()` to
real clinical trial data (`CTdata.rda`). The model formula
for the N-of-1 design is:

```
Sx ~ bm + t + Dbc + bm:Dbc, random = ~1|ptID,
  correlation = corCAR1(form = ~t|ptID)
```

### What is needed

The current model tests the biomarker-treatment
*interaction* (`bm:Dbc`): does blood pressure predict who
responds? To test whether prazosin improves symptoms
*overall*, the coefficient of interest is `Dbc` (the main
drug effect), not `bm:Dbc`. This is a simpler test -- check
whether the drug indicator is significant.

The existing infrastructure supports this; it only requires
extracting a different coefficient from the model output.

### Action -- DONE

Implemented in `R/carryover_analysis.R` as
`analyze_trial_extended()`. Extracts both the main drug
effect (`Dbc`) and the biomarker interaction (`bm:Dbc`)
with 95% CIs, variance components, predicted response
curves, and clinical decision thresholds.

Usage:
```r
result <- analyze_trial_extended(td, dat, op, threshold = 5)
print_trial_summary(result)
```

---

## Goal 2: Characterize carryover in actual data

*Need to figure out how to structure the carryover
quantification and identify the best sequence of response
measures. -- RCH implement with RGT guidance*

### What we have

- The `lambda_cor = ln(2)/t_half` derivation provides a
  principled framework for relating carryover half-life to
  the correlation structure (see
  `pdf/08-biomarker-correlation-decay.pdf`).
- The `Dbc` continuous drug indicator with exponential decay
  is already implemented and parameterized by
  `carryover_t1half`.
- The `nlme::lme` with `corCAR1` analysis model handles the
  AR(1) residual structure.

### What is needed

#### 2a. Structure the carryover quantification

Fit the actual trial data with multiple carryover models
(different half-lives) and compare fit via AIC/BIC or
likelihood ratio test:

```r
results <- list()
for (t_half in c(0.25, 0.5, 1.0, 2.0, 4.0)) {
  # Compute Dbc with this half-life
  # Fit nlme::lme with corCAR1
  # Record AIC, BIC, log-likelihood
  results[[as.character(t_half)]] <- list(
    t_half = t_half,
    aic = AIC(fit),
    bic = BIC(fit),
    loglik = logLik(fit)
  )
}
# Select best-fitting half-life
```

This is a straightforward extension of the existing
`lme_analysis()`: loop over candidate half-lives and
compare model fit. The `Dbc` computation is already
parameterized by `carryover_t1half`.

A new function `characterize_carryover()` would:

1. Accept the trial data and design specification
2. Sweep a range of candidate half-lives
3. Fit the LME model at each half-life
4. Return a comparison table (half-life, AIC, BIC,
   log-likelihood, Dbc coefficient, p-value)
5. Identify the best-fitting half-life

#### 2b. Identify best sequence of response measures

This likely means: at which timepoints does the carryover
effect manifest most strongly? The existing trial design
structure (with `tod`, `tsd`, `tpb` at each timepoint)
provides the timing information. The analysis would:

1. Fit the model with the best carryover half-life
2. Extract per-timepoint residuals
3. Examine whether residuals at off-drug timepoints
   systematically differ from expectations (indicating
   carryover signal)
4. Plot the residual pattern across timepoints to
   visualize where carryover is strongest

### Action -- DONE

Implemented in `R/carryover_analysis.R` as
`characterize_carryover()`. Sweeps candidate half-lives,
fits nlme::lme at each, compares AIC/BIC, and returns the
best-fitting half-life with full model output.

Usage:
```r
cc <- characterize_carryover(td, dat,
  half_lives = c(0.25, 0.5, 1.0, 2.0, 4.0))
print_carryover_summary(cc)
# cc$best_t_half gives the optimal half-life
# cc$best_model gives the full model fit
```

Validated: correctly recovers the true t_half = 1.0 from
simulated data generated with that half-life.

---

## Goal 3: Analysis methods to add carryover

### 3a. Simulation power calculations

*Done, partially; perhaps needs to be more modular/flexible.
-- RGT implement*

#### Status: Essentially complete

The revised simulation code provides:

- `lambda_cor` parameter for correlation decay
  (auto-derived from half-life or manual override for
  sensitivity analysis)
- Carryover half-life as a simulation parameter:
  `c(0, 0.5, 1.0)` weeks (clinically realistic values;
  updated from the publication's `c(0, 0.1, 0.2)` which
  were too short to capture realistic pharmacokinetics)
- `Dbc` continuous drug indicator with exponential decay
  in the analysis model
- 1,000-replicate validated power estimates with proper
  Type I error control (2.5-5.8% across all designs)
- 8-minute runtime (200 reps) enabling rapid sensitivity
  analysis; 47 minutes for 1,000 reps
- Zero positive definiteness failures across all 162
  parameter combinations

#### Remaining modularity/flexibility items

1. **Carryover decay model options.** Currently exponential
   only. The `pmsimstats2025` repo has linear and Weibull
   decay in `calculate_carryover()`. These could be
   backported to allow sensitivity analysis over decay
   shape.

2. **Separate correlation decay rate.** The `lambda_cor`
   parameter currently defaults to `ln(2)/t_half`
   (matching the drug elimination rate). Making it an
   explicit sensitivity parameter allows exploring
   scenarios where biomarker predictive power decays
   faster or slower than the drug effect.

3. **With/without carryover in analysis model.** The
   simulation could compare power under two analysis
   scenarios: (a) carryover modeled in the analysis
   (continuous Dbc), and (b) carryover ignored (binary
   Db). This addresses the question: how much power is
   gained by correctly modeling carryover in the analysis?

4. **Modular design specification.** The trial designs are
   currently hardcoded in the simulation scripts. A
   configuration file or function that generates design
   specifications from high-level parameters (number of
   phases, block durations, randomization points) would
   make it easier to explore alternative designs.

### 3b. Actual analysis

*In progress, reviewed 2/14 how to get the information
needed in the output structure. -- RGT implement*

#### What we have

The `lme_analysis()` with `full_model_out = TRUE` returns:

- `form`: the model formula
- `fit`: the full `nlme::lme` model object
- `data`: the analysis-ready data with all derived
  variables (`Dbc`, `tsd`, `De`, etc.)
- `summary`: the standard output (beta, betaSE, p-value)

The full model object provides access to:

- `summary(fit)$tTable`: all coefficients with SEs and
  p-values
- `fitted(fit)`: predicted values per observation
- `residuals(fit)`: residuals per observation
- `ranef(fit)`: random effects (per-participant intercepts)
- `intervals(fit)`: confidence intervals for all parameters
- `VarCorr(fit)`: variance components (random intercept
  variance, residual variance, AR(1) parameter)

#### What is needed for the output structure

1. **Estimated carryover half-life** from the model
   comparison in Goal 2 (the best-fitting `t_half`).

2. **Predicted response curves per participant** using the
   fitted model: `predict(fit, newdata = ...)` with
   participant-specific random effects.

3. **Biomarker-specific treatment effect as a continuous
   function of blood pressure.** The `bm:Dbc` coefficient
   gives the slope: for each unit increase in blood
   pressure, the drug effect changes by `beta_4` CAPS
   points. Combined with the main `Dbc` effect, the
   predicted drug effect for a participant with blood
   pressure `bm` is:
   `effect(bm) = beta_Dbc + beta_bm:Dbc * (bm - mean_bm)`

4. **Clinical decision threshold.** At what blood pressure
   value does the predicted benefit exceed a clinically
   meaningful threshold (e.g., 5-point CAPS improvement)?
   This is: `bm_threshold = mean_bm +
   (threshold - beta_Dbc) / beta_bm:Dbc`

### Action

**DONE.** Implemented as `analyze_trial_extended()` in
`R/carryover_analysis.R`. Items 1-4 are all included:
estimated carryover half-life (via `characterize_carryover`),
predicted response curves (`result$predicted_effect(bm)`),
biomarker-specific treatment effect with CIs, and clinical
decision threshold (`result$threshold_bm`).

---

## Summary of Actions

| Goal | Status | Action | Location |
|------|--------|--------|----------|
| 1. Prazosin efficacy | **DONE** | `analyze_trial_extended()` | `R/carryover_analysis.R` |
| 2a. Carryover quantification | **DONE** | `characterize_carryover()` | `R/carryover_analysis.R` |
| 2b. Best response sequence | Conceptual | Residual analysis at each timepoint | Needs implementation |
| 3a. Simulation power | **DONE** | Full simulation with optimizations | `R/generateSimulatedResults.R` |
| 3b. Actual analysis output | **DONE** | Extended output with predictions | `R/carryover_analysis.R` |

---

## Key Deliverables from the Audit

The following documents and code are available:

### White Papers

- `pdf/04-revised-power-analysis.pdf` (13 pages): Complete
  revision summary for the authors, with before/after
  Figure 4 heatmaps and all code changes
- `pdf/18-response-parameter-sensitivity.pdf` (6 pages):
  Response parameter sensitivity analysis
- `pdf/08-biomarker-correlation-decay.pdf` (10 pages):
  Pharmacokinetic derivation of the correlation decay rule
- `pdf/06-ar1-residual-correlation.pdf` (8 pages): Analysis
  model transition from lmer to nlme with corCAR1
- `pdf/07-positive-definiteness-failures.pdf` (10 pages):
  Positive definiteness failure analysis
- `13-figure4-walkthrough.md`: Complete technical walkthrough
  of Figure 4 methodology and code

### Simulation Scripts

- `analysis/figure4/04-generate-power-head.R`: Figure 4 simulation (current code)
- `analysis/figure4/01-generate-power-sweep.R`: Figure 4 power sweep (sequential)
- `analysis/figure4/02-generate-power-sweep-original.R`: Publication reproduction
- `analysis/figure4/03-generate-power-by-commit.R`: Per-commit comparison
- `analysis/figure4/05-pd-diagnostics.R`: PD failure instrumentation
- `analysis/figure4/06-plot-power-heatmaps.R`: Heatmap plotting
- `analysis/figure5/01-generate-parameter-sensitivity.R`: Figure 5 simulation (revised)
- `analysis/figure5/02-generate-parameter-sensitivity-original.R`: Figure 5 (original)

### Key Results

- 1,000-replicate Figure 4 with zero PD failures
- Type I error: 2.5-5.8% (all nominal)
- N-of-1 carryover sensitivity confirmed at realistic
  half-lives (0.5-1.0 weeks)
- Design ranking preserved: N-of-1 > CO > OL+BDC > OL
- Runtime: 8 minutes (200 reps), 47 minutes (1,000 reps)
