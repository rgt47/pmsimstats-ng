# Original Implementation (Architecture B Only)

## Collection Label: `original`

The corrected Hendrickson et al. (2020) publication code with seven documented DGP improvements applied. This collection implements Architecture B (MVN differential correlation) exclusively and serves as a baseline reference.

## Architecture Support

- **Architecture B only (MVN Differential Correlation):** Biomarker-response interaction via treatment-state-dependent correlation in the covariance structure

No `dgp_architecture` parameter — all functions operate in MVN mode.

## Coding Style

- data.table for efficient data manipulation
- Base R and data.table operations
- camelCase function names
- Parameter lists (list-based function signatures)

## File Structure

- `R/generateData.R` — MVN covariance construction and data generation
- `R/generateSimulatedResults.R` — Batch simulation orchestration
- `R/buildtrialdesign.R` — Trial design specification
- `R/lme_analysis.R` — LME analysis model (nlme::lme with corCAR1)
- `R/utilities.R` — Helper functions (Gompertz, cumulative effects)
- `R/censordata.R` — Missing data mechanisms
- `R/carryover_analysis.R` — Extended analysis (half-life characterization, clinical thresholds)

## Key Parameters

| Parameter | Type | Meaning |
|---|---|---|
| `c.bm` | numeric | BM-BR correlation coefficient (0-1) under MVN architecture |
| `carryover_t1half` | numeric | Carryover half-life (weeks). 0 = no carryover |
| `lambda_cor` | numeric | Correlation decay rate (per week). Auto-derived as ln(2)/carryover_t1half if NA |
| `dgp_architecture` | character | Not present; architecture is fixed to "mvn" |

## Biological Assumptions

Architecture B assumes the biomarker's predictive value emerges from **shared variance** with a drug-responsive physiological component. As the drug washes out, both the biomarker and drug-responsive component decay in parallel, so their correlation weakens. This is appropriate when:
- The biomarker reflects current drug pathway engagement (e.g., active metabolite concentration)
- The biomarker is a dynamic biomarker (changes with treatment status)

## Living Copy Note

This collection is maintained as a **static copy of master's R/ directory** that is manually updated when master changes. Use this collection to reference the baseline (Architecture B-only) code.

## Seven DGP Corrections Applied

1. AR(1) autocorrelation (instead of compound symmetry)
2. Exponential BM-BR correlation decay during off-drug periods
3. nlme::lme with corCAR1 residual structure in analysis model
4. Removed undocumented carryover scale factor
5. Fixed timepoint-1 correlation assignment
6. Corrected Gompertz offset for origin-passing behavior
7. Updated default parameters to reflect realistic carryover half-lives

See `docs/03-audit-and-revision-report.md` for detailed rationale.

## Cross-Reference

- **tidyverse:** Tidyverse reimplementation, Architecture A + B (both architectures available)
- **original-extended:** This codebase extended with Architecture A support (dgp_architecture parameter)

## Usage

See `original-extended/README.md` for usage examples. Functions are identical between `original` and `original-extended` when called with `dgp_architecture = "mvn"`.
