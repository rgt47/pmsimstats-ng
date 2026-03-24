# Code Review: analysis/2025/pm_functions.R

## Dead Code

### 1. `calculate_carryover()` (lines 24-77)

Multi-model carryover function supporting exponential, linear, and
Weibull decay. **Never called** by any function in the file. The
actual carryover is computed inline in `apply_carryover_to_component()`
using the exponential formula `(1/2)^(tsd/halflife)`. The linear
and Weibull models were never integrated into the DGP.

**Action:** Remove entirely.

### 2. `calculate_component_halflives()` (lines 138-145)

Returns a list with three identical values (`br_halflife`,
`pb_halflife`, `tv_halflife`), all set to `base_halflife`. Only
`br_halflife` is ever used (the function guards on
`component_name != "br"` in `apply_carryover_to_component()`).

**Action:** Inline the single value. Replace all calls with
`model_param$carryover_t1half` directly.

### 3. `calculate_bio_response_with_interaction()` (lines 95-126)

Thin wrapper that:

1. Calls `mod_gompertz()` on `trial_data$tod`
2. Sets `bio_response_test <- base_bio_response == 0`
3. Calls `apply_carryover_to_component()`
4. Returns a 3-element list

This adds indirection without value. The orig code does this inline
in `buildSigma()` in 10 lines. The `raw_bio_response_means` return
field is never used by any caller.

**Action:** Inline into `generate_data()` and `build_sigma_matrix()`.

### 4. `validate_correlation_structure()` (lines 778-813)

Diagnostic function. Calls `build_sigma_matrix()` and reports
eigenvalues and condition number. Not called by the simulation
pipeline. Was useful during development but is now redundant with
`validate_parameter_grid()`.

**Action:** Remove. The orig's `validateParameterGrid()` (ported
as Step 7c) covers this use case.

### 5. `create_sigma_cache_key()` (lines 815-822)

Creates a cache key from `design_name`, `params$n_participants`,
`params$biomarker_correlation`, `params$carryover_t1half`. **Never
called.** The actual caching in `generate_simulated_results()` uses
`paste(vpg indices, sep='_')` directly.

**Action:** Remove entirely.

### 6. `validate_parameter_grid()` (lines 828-991)

170-line validation function with extensive console output, emoji,
condition number statistics, and recommendations. References
`param_grid$biomarker_correlation` and `param_grid$n_participants`
-- column names from the old 2025 parameter format that no longer
exists after the alignment. Would crash if called.

**Action:** Remove. The orig's `validateParameterGrid()` and the
ported `generate_simulated_results()` handle pre-flight validation.

### 7. `report_parameter_validation()` (lines 997-1044)

Companion reporting function for `validate_parameter_grid()`. Same
issue -- references old column names. Never called by the pipeline.

**Action:** Remove.

### 8. Deprecated docstring block (lines 147-158)

Roxygen comment block for a function that was already removed.
Describes the deprecated `calculate_carryover_adjusted_correlations()`
which no longer exists.

**Action:** Remove.

## Duplication

### 9. `generate_data()` and `build_sigma_matrix()` share ~100 lines

Both functions independently compute: labels, standard_deviations,
means (including the Gompertz + carryover loop), and the correlation
matrix. The logic is copy-pasted. In orig, `buildSigma()` is a
standalone function and `generateData()` calls it (or uses a
cached version).

**Action:** Make `generate_data()` always call `build_sigma_matrix()`
when no cached sigma is provided, eliminating the duplicated sigma
construction path (lines 433-582 of generate_data). This matches
the orig architecture where `generateData()` calls `buildSigma()`.

### 10. `factor_types` and `factor_abbreviations` are identical

After the alignment, `factor_types = c("tv", "pb", "br")` and
`factor_abbreviations = c("tv", "pb", "br")` are the same vector.
Both are passed around as separate arguments and iterated in
parallel. The `current_factor` / `factor_abbrev` distinction no
longer exists.

**Action:** Collapse to a single `factors <- c("tv", "pb", "br")`
throughout.

## Simplifications

### 11. `apply_carryover_to_component()` over-engineered

Takes `component_halflives` (a list) and `component_name` (a
string), then immediately returns unchanged if `component_name !=
"br"`. This function is called three times in a loop -- twice it
returns immediately. The function could simply accept the halflife
value directly and only be called for BR.

**Action:** Call only for BR component. Pass
`model_param$carryover_t1half` directly instead of the halflives
list. Remove the `component_name` guard.

### 12. `prepare_trial_data()` fallback logic

Lines 238-246 handle a `week` column fallback that is never
exercised. All trial designs in this codebase use `t_wk`. The
`on_drug` column it adds duplicates information already available
from `tod > 0`.

**Action:** Simplify to: add `t_wk_cumulative <- cumsum(t_wk)` and
`on_drug <- (tod > 0)`. Remove the `week` fallback.

### 13. `build_correlation_matrix()` signature has 11 parameters

Several are vestigial: `bio_response_test`, `bio_response_means`,
`means` are only used for the BM-BR correlation block, and even
there `bio_response_test` is only checked when `trial_data` is
NULL (which never happens in practice since `trial_data` is always
passed). The `factor_types` parameter is now identical to
`factor_abbreviations`.

**Action:** Simplify signature. Remove `bio_response_test`,
`bio_response_means`, `means`, `factor_types`. Use `trial_data`
(always present) for the on-drug check.

### 14. `process_participant_data()` pre-allocates then overwrites

Lines 1054-1058 create a zero-filled tibble and bind it, then
lines 1060-1079 overwrite every column. The pre-allocation has no
benefit since `mutate()` creates new columns regardless.

**Action:** Remove the pre-allocation block. Let the `mutate()`
calls create the columns directly.

### 15. PD tracking in `generate_data()` is single-call accounting

`sigma_count` and `non_positive_definite_count` are initialized to
0, incremented at most once (the function builds one sigma), and
then attached as attributes. The `non_positive_definite_rate` is
always 0 or 1. This tracking made sense in a loop but is trivial
for a single-sigma function.

**Action:** Replace with a single boolean attribute
`was_pd_corrected`.

## Estimated Impact

| Category | Items | Lines removed | Lines simplified |
|----------|------:|:------------:|:----------------:|
| Dead code | 7 | ~350 | -- |
| Duplication | 2 | ~100 | ~50 |
| Simplifications | 5 | ~30 | ~60 |
| **Total** | **14** | **~480** | **~110** |

Current file: 1650 lines. After cleanup: ~1060 lines.

---
*Generated 2026-03-24.*
