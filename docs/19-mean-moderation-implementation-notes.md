# Mean Moderation Implementation Notes

## Summary

This document records the implementation of Architecture A
(direct mean moderation) for the biomarker-treatment
interaction in `pmsimstats-ng`, including two failed
approaches and the final correct formulation derived from
the original 2025 codebase.

## Background

The original R package code (`R/generateData.R`) uses
Architecture B: the BM-BR interaction is encoded as a
differential correlation in the covariance matrix of the
joint MVN draw. Under carryover, this correlation erodes as
the drug effect decays, producing the 40-60% power drops
documented in `05-revision-summary.pdf`.

The 2025 tidyverse codebase (pre-alignment version in
`~/prj/alz/01-pmsimstats/pmsimstats2025`, initial commit
957c8bd) used Architecture A: additive mean moderation
applied post-MVN draw. That version showed full power
preservation under carryover.

## Implementation Attempts

### Attempt 1: Multiplicative scaling (incorrect)

```r
dat[, (br_col) := get(br_col) * (1 + beta_bm * bm_centered)]
```

**Problem:** When the base BR value is near zero (off-drug
timepoints with minimal carryover), multiplying by
`(1 + small)` produces negligible moderation signal. The
analysis model sees strong biomarker signal on-drug but
almost none off-drug, mimicking Architecture B's correlation
erosion. Power dropped ~60% under carryover -- no
improvement over MVN.

### Attempt 2: Additive with drug exposure decay (incorrect)

```r
drug_exposure <- (1/2)^(tsd / t1half)  # on-drug=1, decays off-drug
dat[, (br_col) := get(br_col) + beta_bm * bm_centered *
  drug_exposure * br_sd]
```

**Problem:** Drug exposure decays rapidly during off-drug
periods (e.g., 0.18 at first off-drug timepoint, 0.001 by
the fourth). The additive signal is proportional to this
decay, so the off-drug moderation is still negligible
relative to noise (br_sd = 5). The `bm:Dbc` interaction
term sees minimal contrast. Power remained low under
carryover.

### Attempt 3: Additive with time-on-drug (correct)

```r
dat[, (br_col) := get(br_col) + bm * tod * beta_bm]
```

**Source:** Original 2025 codebase, initial commit 957c8bd,
`analysis/scripts/pm_functions.R` lines ~650-670:

```r
treatment_status <- trial_design$tod[timepoint_idx]
participant_data[[br_col]] <- participant_data[[br_col]] +
  biomarker * treatment_status * interaction_strength
```

**Why this works:** The moderation term is proportional to
`tod` (cumulative time on drug), not to the current drug
exposure level. On-drug timepoints accumulate `tod`
monotonically (2.5, 5.0, 7.5, ...), producing a strong,
growing biomarker signal. Off-drug timepoints have `tod = 0`
(never treated in that path), so no moderation is applied
there. The signal lives entirely in the mean structure and
is orthogonal to the covariance -- carryover cannot erode
it.

**Result:** Power = 1.00 across all carryover half-lives
(0, 0.5, 1.0) at N=70, c.bm=0.45 for both N-of-1 and
OL+BDC designs.

## Architectural Difference

The key distinction between Architectures A and B is where
the interaction signal lives:

- **Architecture B (MVN):** Signal is in the second moment
  (covariance). The analysis model infers the interaction
  from conditional correlations between bm and BR. Carryover
  erodes the off-drug correlation, compressing the contrast.

- **Architecture A (mean moderation):** Signal is in the
  first moment (mean). The moderation is additive and
  proportional to cumulative drug exposure. The analysis
  model detects a direct shift in mean outcome conditional
  on biomarker value. This shift is unaffected by carryover
  because it depends on `tod`, not on current drug status.

## Implementation Details

The `dgp_architecture` parameter is threaded through:

- `buildSigma()`: skips BM-BR correlation block when
  `dgp_architecture = "mean_moderation"`
- `generateData()`: applies additive `bm * tod * c.bm`
  post-draw
- `generateSimulatedResults()`: passes parameter to both
  sigma cache and data generation

Default is `"mvn"` to preserve backward compatibility. All
67 existing tests pass with the default.

## Files Modified

- `R/generateData.R`: `buildSigma()` and `generateData()`
- `R/generateSimulatedResults.R`: parameter threading

## Branch

`feature/mean-moderation-dgp`
