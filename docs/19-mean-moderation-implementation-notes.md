# Mean Moderation Implementation Notes

> **Note**: This is a quick-reference summary. For the full
> mathematical derivation, the positive-definiteness ceiling
> discussion, the `validateParameterGrid()` known limitation,
> and the empirical evaluation of the Hendrickson (2020)
> publication parameter grid under Architecture A, see the
> comprehensive document `19-mean-moderation-implementation-notes.tex`
> (rendered as `19-mean-moderation-implementation-notes.pdf`).

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

### Attempt 3: Additive with time-on-drug (incorrect scaling)

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

**Problem:** With uncentered biomarker (mean ~124) and
`tod` scaling with weeks (up to 10-20), the moderation
shift was ~650 units against BR noise SD of ~8. Power was
1.00 everywhere -- an artifact of an absurdly overpowered
signal, not genuine resilience to carryover. The original
2025 codebase had this same scaling issue.

### Attempt 4: Centered, on-drug binary, calibrated (correct)

```r
bm_z <- (dat$bm - bm_mean) / bm_sd
for (tp in 1:nP) {
  if (d[tp]$onDrug) {
    dat[, (br_col) := get(br_col) + c.bm * bm_z * br_sd]
  }
}
```

**Why this works:** The formulation matches Architecture B's
conditional expectation from the MVN model:

$$E[BR | bm, \text{on drug}] = \mu_{BR} + c_{bm} \cdot
\frac{\sigma_{BR}}{\sigma_{bm}} \cdot (bm - \mu_{bm})$$

The standardized biomarker `bm_z` centers the signal
(mean moderation = 0 in the population). The `br_sd`
scaling produces comparable effect sizes between
architectures: a 1-SD biomarker shift produces
`c.bm * sigma_br` units of BR shift. The binary on-drug
indicator ensures moderation is present whenever drug is
active and absent otherwise.

**Calibration (c.bm = 0.45, extracted parameters):**

- Architecture A: 3.36 BR units per SD of biomarker
- Architecture B: 3.48 BR units per SD of biomarker
- Signal-to-noise ratio: ~0.65-0.9 (realistic)

**Power (N-of-1, N=70, c.bm=0.45, Nreps=50):**

| $t_{1/2}$ | Architecture A | Architecture B |
|:---:|:---:|:---:|
| 0.0 | 0.74 | 0.82 |
| 0.5 | 0.68 | 0.64 |
| 1.0 | 0.72 | 0.50 |

Architecture A power is stable across carryover levels
(variation is Monte Carlo noise at Nreps=50). Architecture B
shows the characteristic 40% relative decline.

## Architectural Difference

The key distinction between Architectures A and B is where
the interaction signal lives:

- **Architecture B (MVN):** Signal is in the second moment
  (covariance). The analysis model infers the interaction
  from conditional correlations between bm and BR. Carryover
  erodes the off-drug correlation, compressing the contrast.

- **Architecture A (mean moderation):** Signal is in the
  first moment (mean). The moderation is a fixed additive
  shift to BR during on-drug timepoints. Off-drug timepoints
  carry no moderation signal regardless of carryover state.
  The analysis model detects the interaction from the
  differential mean between high-biomarker and
  low-biomarker participants when on drug.

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
