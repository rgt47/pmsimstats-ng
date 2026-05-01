# pmsimstats Update Task Lists

Based on the audit and revision of the `orig` repository,
the following updates should be applied to the sibling
repositories. Each task references the specific finding
from the audit that motivates it.

---

## Repository: pmsimstats2025

Location: `~/Dropbox/prj/alz/01-pmsimstats/pmsimstats2025/`
File: `analysis/scripts/pm_functions.R`

### Priority 1: Statistical Correctness

#### Task 1.1: Remove carryover scale factor

**File:** `pm_functions.R`, `apply_carryover_to_component()`
(line ~260)

**Current code:**
```r
decay_factor <- (1/2)^(scale_factor * time_lag /
                       component_halflife)
```

**Change to:**
```r
decay_factor <- (1/2)^(time_lag / component_halflife)
```

Also remove `scale_factor` parameter from
`generate_data()` signature (line ~428) and
`build_sigma_matrix()` (lines ~704, 720, 729).

**Rationale:** The scale factor of 2 halves the effective
half-life without documentation or justification. The
parameter `carryover_t1half` should mean what it says.
See audit Section 4.1, Change 1.

---

#### Task 1.2: Fix timepoint-1 BM-BR correlation

**File:** `pm_functions.R`, `build_correlation_matrix()`
(lines ~392-414)

**Current code:**
```r
if (timepoint_idx > 1) {
  # ... set correlation
}
```

**Change to:**
```r
# No guard -- handle all timepoints
if (!is.null(bio_response_test) &&
    !bio_response_test[timepoint_idx]) {
  correlations["biomarker", name1] <- model_param$c.bm
  correlations[name1, "biomarker"] <- model_param$c.bm
} else if (trial_data$tsd[timepoint_idx] > 0 &&
           lambda_cor > 0) {
  decay <- exp(-lambda_cor * trial_data$tsd[timepoint_idx])
  correlations["biomarker", name1] <- model_param$c.bm * decay
  correlations[name1, "biomarker"] <- model_param$c.bm * decay
}
```

**Rationale:** The `p > 1` guard was a coding artifact to
avoid an index error in the proportional scaling formula.
Participants are on drug for 2.5-4 weeks at timepoint 1.
See audit Section 4.1, Change 2.

---

#### Task 1.3: Replace BM-BR correlation rule with
exponential decay

**File:** `pm_functions.R`, `build_correlation_matrix()`
(lines ~392-414) and `generate_data()` signature

**Current code (step function):**
```r
if (!bio_response_test[timepoint_idx]) {
  correlations["biomarker", name1] <- model_param$c.bm
} else {
  correlations["biomarker", name1] <- 0
}
```

**Change to:**
```r
if (trial_data$on_drug[timepoint_idx]) {
  correlations["biomarker", name1] <- model_param$c.bm
  correlations[name1, "biomarker"] <- model_param$c.bm
} else if (trial_data$tsd[timepoint_idx] > 0 &&
           lambda_cor > 0) {
  decay <- exp(-lambda_cor * trial_data$tsd[timepoint_idx])
  correlations["biomarker", name1] <- model_param$c.bm * decay
  correlations[name1, "biomarker"] <- model_param$c.bm * decay
}
```

Add `lambda_cor` parameter to `generate_data()` with
default `NA`, auto-computed as `log(2) / carryover_t1half`.

**Rationale:** Derived from pharmacokinetic first principles.
When each participant's drug effect decays by the same
fraction, the biomarker-specific spread decays
proportionally, preserving the coefficient of variation.
See correlation_decay_derivation.pdf.

---

#### Task 1.4: Switch autocorrelation from compound symmetry
to AR(1)

**File:** `pm_functions.R`, `build_correlation_matrix()`
(lines ~318-345)

**Current code:**
```r
correlations[name1[idx], name2[idx]] <- autocorrelation
```

**Change to:**
```r
time_gap <- abs(weeks[point_indices$p2[idx]] -
                weeks[point_indices$p1[idx]])
correlations[name1[idx], name2[idx]] <- autocorrelation^time_gap
```

Also apply AR(1) decay to cross-factor different-time
correlations (lines ~364-388):
```r
correlations[name1[idx], name2[idx]] <-
  model_param$c.cfct * autocorrelation^time_gap
```

**Rationale:** AR(1) expands the PD-feasible c.bm range
from 0.25 to 0.49+ and is more clinically realistic
(nearby timepoints more correlated than distant ones).
See pd_failure_white_paper.pdf.

---

#### Task 1.5: Switch analysis model from lmer to nlme with
corCAR1

**File:** `pm_functions.R`, `lme_analysis()` (lines ~1170-1298)

**Current code:**
```r
model <- lmer(formula, data = data_for_model)
```

**Change to:**
```r
model <- tryCatch(
  nlme::lme(formula, random = ~1|participant_id,
    correlation = nlme::corCAR1(
      form = ~t|participant_id),
    data = data_for_model,
    control = nlme::lmeControl(
      opt = "optim", maxIter = 200,
      msMaxIter = 200)),
  error = function(e) {
    nlme::lme(formula, random = ~1|participant_id,
      data = data_for_model,
      control = nlme::lmeControl(opt = "optim"))
  }
)
```

Update coefficient extraction from
`summary(model)$coefficients` to `summary(model)$tTable`,
with column names `Value`, `Std.Error`, `p-value`.

Add graceful handling for rank-deficient models (check
if coefficient exists before extracting).

Remove NA rows before fitting (`nlme` uses `na.fail`).

**Rationale:** The AR(1) DGP produces time-dependent
residual correlations that `lmer` (random intercept only)
cannot absorb, inflating Type I error to 13-17% for CO
and OL designs. `corCAR1` corrects this to nominal 5%.
See analysis_model_ar1.pdf.

---

#### Task 1.6: Also apply to Hendrickson-aligned functions

The Hendrickson-aligned functions in `pm_functions.R`
(lines ~1305-1620) implement a separate code path:

- `modgompertz_orig()` -- correct, no change needed
- `build_path_sigma()` -- needs AR(1) and lambda_cor
- `build_bm_br_correlations()` -- needs exponential decay
  instead of proportional scaling
- `compute_tsd_orig()`, `compute_tod_orig()`,
  `compute_tpb_orig()` -- correct, no change needed

Apply the same AR(1) and lambda_cor changes to
`build_path_sigma()` and replace the `build_bm_br_correlations()`
proportional scaling with exponential decay.

---

### Priority 2: Parameter Constraints

#### Task 2.1: Update default correlation parameters

**File:** `pm_functions.R`, `hendrickson_corr_params` and
simulation scripts

**Change:**
```r
# Was:
c.tv = 0.8, c.pb = 0.8, c.br = 0.8

# Now (AR(1) base rate):
c.tv = 0.7, c.pb = 0.7, c.br = 0.7
```

---

#### Task 2.2: Constrain c.bm range

**File:** Simulation scripts

**Change:** Limit `c.bm` to `{0, 0.25, 0.5}` instead of
`{0.2, 0.4}` or `{0, 0.3, 0.6}`. The maximum PD-feasible
value under AR(1) rho=0.7 is 0.49 for the OL design and
0.45 for OL+BDC.

---

#### Task 2.3: Update carryover half-life range

**File:** Simulation scripts and CLAUDE.md

**Change:** Replace `c(0, 0.1, 0.2)` weeks with
`c(0, 0.5, 1.0)` weeks. The publication values (0.1-0.2
weeks) are too short to produce genuine carryover effects.
Update CLAUDE.md which currently mandates the legacy values.

---

### Priority 3: Code Quality

#### Task 3.1: Remove deprecated function

Delete `calculate_carryover_adjusted_correlations()` (lines
~156-205) which is already marked deprecated.

---

#### Task 3.2: Prune unused dependencies

Remove from DESCRIPTION Imports: `DBI`, `RMySQL`,
`RPostgres`, `RSQLite`, `odbc`, `palmerpenguins`, `naniar`,
`visdat`, `janitor`, `skimr`, `covr`, `doParallel`,
`foreach`, `conflicted`, `sessioninfo`, `digest`,
`bookdown`.

---

#### Task 3.3: Archive exploratory scripts

Move the 30+ exploratory scripts
(`check_model_specification*.R`, `debug_*.R`,
`diagnose_*.R`, `test_*.R`) to an `analysis/archive/`
directory. Keep only the canonical simulation workflow.

---

#### Task 3.4: Resolve dual Gompertz functions

The codebase has two Gompertz implementations:
- `mod_gompertz()` (modernized, does NOT pass through origin)
- `modgompertz_orig()` (Hendrickson original, passes through
  origin)

These produce different trajectories. Remove
`mod_gompertz()` and use `modgompertz_orig()` exclusively,
or document when each should be used.

---

### Priority 4: Add nlme to Dependencies

**File:** DESCRIPTION

Add `nlme` to Imports (already listed but confirm it is
used in the analysis path).

---

## Repository: pmsimstats-simple

Location: `~/Dropbox/prj/alz/01b-pmsimstats-simple/pmsimstats-simple/`
File: `analysis/scripts/simulation.R`

### Priority 1: Bug Fixes

#### Task S1.1: Fix global br_rate reference in analysis

**File:** `simulation.R`, `analyze_trial()` (line ~284)

**Current code:**
```r
adjusted_response = response - carryover * br_rate
```

**Problem:** `br_rate` is a global variable from the DGP,
not available to the analysis function. In a real analysis,
the true drug effect rate would be unknown.

**Options:**
1. Remove the carryover adjustment from the analysis
   entirely (simplest, honest about what the analysis
   can know)
2. Estimate the carryover effect from the data (more
   complex, but realistic)
3. Document as an oracle analysis (the simulation assumes
   known carryover parameters, providing an upper bound
   on power)

**Recommended:** Option 3 -- add a comment documenting that
this is an oracle analysis, and note in the output that
power estimates are upper bounds.

---

#### Task S1.2: Fix stateless carryover

**File:** `simulation.R`, `generate_participant()` (lines
~150-162)

**Current code:**
```r
weeks_since <- w - last_drug_week
carryover_factors[w] <- calc_carryover(weeks_since,
                                       carryover_halflife)
```

**Problem:** Carryover is computed statelessly from
`weeks_since_last_drug`. In the original and revised DGP,
carryover propagates sequentially: the residual at time t
depends on the carryover-adjusted value at t-1, not the
raw on-drug value.

**Change to:**
```r
if (on_drug[w]) {
  last_effect <- br_rate * (1 + biomarker_mod * biomarker)
  carryover_factors[w] <- 0
} else if (last_drug_week > 0) {
  weeks_since <- w - last_drug_week
  carryover_factors[w] <- calc_carryover(weeks_since,
                                         carryover_halflife)
}
```

This is already close to the current implementation.
The key difference is whether the carryover decay uses
`last_drug_week` (current, stateless) or propagates
through intermediate off-drug timepoints (sequential).
For the simple repo's weekly timepoints with no
re-initiation, the two approaches are equivalent. Document
this equivalence.

---

### Priority 2: Documentation

#### Task S2.1: Document simplifications relative to orig

**File:** `docs/simplification-plan.md`

Add a section listing the specific simplifications and
their consequences, now that we understand the orig code
in detail:

| Simplification | Consequence |
|---------------|-------------|
| Direct mean moderation vs MVN correlation | No PD problem; no correlation decay issue; simpler but different statistical mechanism |
| OLS on phase means vs LME | No residual autocorrelation issue; potentially different power characteristics |
| Fixed sd vs Gompertz-derived sd | Less realistic variance structure |
| No expectancy scaling on variance | Avoids OL+BDC PD failures entirely |

---

#### Task S2.2: Add Type I error validation

**File:** `simulation.R`

Add a validation block that runs the simulation at
`biomarker_mod = 0` and reports the rejection rate. This
should be approximately 5%. Currently no Type I error
check exists.

---

#### Task S2.3: Add tests

**File:** `tests/testthat/test-basic.R`

The current test is `expect_true(TRUE)`. Add:
- Test that `calc_carryover()` produces expected decay
- Test that `generate_response()` returns 0 when off drug
  with no carryover
- Test that `analyze_trial()` returns a p-value between
  0 and 1
- Test that Type I error at `biomarker_mod = 0` is
  within binomial CI of 0.05 (at low Nreps)

---

### Priority 3: Parameter Updates

#### Task S3.1: Update carryover half-life range

**File:** `simulation.R` (line ~37)

**Current:** `carryover_halflife_levels <- c(0, 0.5, 1, 2)`

These are reasonable and more clinically realistic than the
publication's 0.1-0.2. No change needed.

---

#### Task S3.2: Seed handling

**File:** `simulation.R` (line ~18)

**Current:** `set.seed(42)` at the top with sequential
iteration through conditions.

**Problem:** The random number sequence is not independent
across conditions, which could introduce subtle
correlations in power estimates.

**Change to:** Set a per-condition seed:
```r
set.seed(42 + condition_index)
```

Or use `withr::with_seed()` for each condition.

---

## Estimated Effort

| Repository | Priority 1 | Priority 2 | Priority 3 | Total |
|-----------|:----------:|:----------:|:----------:|:-----:|
| pmsimstats2025 | 6 tasks (1-2 days) | 3 tasks (half day) | 4 tasks (half day) | ~3 days |
| pmsimstats-simple | 2 tasks (2 hours) | 3 tasks (2 hours) | 2 tasks (1 hour) | ~5 hours |

---

*Generated 2026-03-22 from the pmsimstats orig audit.*
*See pmsimstats_revision_summary.pdf for full context.*
