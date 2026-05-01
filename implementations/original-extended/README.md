# Original Extended (Architecture A + B)

## Collection Label: `original-extended`

The corrected Hendrickson et al. (2020) code extended with support for Architecture A (direct mean moderation). This is the most capable implementation, supporting both DGP architectures via the `dgp_architecture` parameter.

## Architecture Support

- **Architecture A (Mean Moderation):** Direct mean scaling of biomarker response (`dgp_architecture = "mean_moderation"`)
- **Architecture B (MVN Differential Correlation):** Biomarker-response interaction via covariance structure (`dgp_architecture = "mvn"`, default)

## Coding Style

- data.table for efficient data manipulation
- Base R and data.table operations
- camelCase function names
- Parameter lists (list-based function signatures)
- Identical to `original` collection except for `dgp_architecture` parameter threading

## File Structure

Same as `original` with Architecture A support:
- `R/generateData.R` — Covariance construction + post-hoc mean moderation
- `R/generateSimulatedResults.R` — Batch simulation with `dgp_architecture` parameter
- `R/buildtrialdesign.R` — Trial design
- `R/lme_analysis.R` — LME analysis
- `R/utilities.R` — Helpers
- `R/censordata.R` — Censoring
- `R/carryover_analysis.R` — Extended analysis

## Key Parameters

| Parameter | Type | Meaning |
|---|---|---|
| `c.bm` | numeric | Under MVN: BM-BR correlation. Under mean moderation: regression coefficient scaling BR by centered biomarker |
| `carryover_t1half` | numeric | Carryover half-life (weeks). 0 = no carryover |
| `lambda_cor` | numeric | Correlation decay rate (per week). Auto-derived as ln(2)/carryover_t1half if NA |
| `dgp_architecture` | character | `"mvn"` (default) or `"mean_moderation"` |

## Biological Assumptions

### Architecture B (default): Shared Variance Mechanism
The biomarker reflects current drug pathway engagement. Its predictive value depends on active drug presence. Appropriate when:
- The biomarker is a dynamic biomarker (changes with treatment status)
- The biomarker indexes a drug-responsive physiological state
- Example: blood pressure as a PTSD subtype marker in the Hendrickson prazosin study

### Architecture A: Deterministic Moderation
The biomarker governs the magnitude of the drug's biological effect for each individual. The relationship is preserved under carryover. Appropriate when:
- The biomarker determines effective drug dose (pharmacokinetics)
- The biomarker is a genetic or structural baseline (unchanging)
- Example: CYP450 metabolizer phenotype determining drug activation

## Effect Sizes Under Carryover

| Architecture | Design | t_half=0 | t_half=0.5 | t_half=1.0 | Power Loss |
|---|---|---|---|---|---|
| **A (Mean Mod)** | N-of-1 | 0.74 | 0.72 | 0.68 | 8% |
| **B (MVN)** | N-of-1 | 0.82 | 0.64 | 0.50 | 39% |

Architecture B is more sensitive to carryover because the differential correlation signal erodes.

## Usage Example

```r
source("implementations/original-extended/R/generateData.R")
source("implementations/original-extended/R/generateSimulatedResults.R")
source("implementations/original-extended/R/lme_analysis.R")

# Build trial design
td <- buildtrialdesign(nP = 20, design = "Hybrid", carryover_t1half = 1.0)

# Parameters
modelparam <- list(N = 50, c.bm = 0.45, carryover_t1half = 1.0)
respparam <- data.table(
  cat = rep(c("tv", "pb", "br"), each = 1),
  max = c(50, 20, 30),
  disp = c(4, 3, 2),
  rate = c(0.5, 0.3, 0.4),
  sd = c(10, 8, 5)
)
blparam <- data.table(
  cat = c("bm", "BL"),
  m = c(120, 60),
  sd = c(15, 10)
)

# Generate data (Architecture A)
dat_a <- generateData(
  modelparam, respparam, blparam, td,
  empirical = FALSE, makePositiveDefinite = TRUE,
  dgp_architecture = "mean_moderation"
)

# Generate data (Architecture B)
dat_b <- generateData(
  modelparam, respparam, blparam, td,
  empirical = FALSE, makePositiveDefinite = TRUE,
  dgp_architecture = "mvn"
)
```

## Batch Simulation

```r
# Batch simulation
results <- generateSimulatedResults(
  designparam = list(carryover_t1half = c(0, 0.5, 1.0)),
  respparamsets = respparamsets,
  blparamsets = blparamsets,
  modelparamsets = modelparamsets,
  n_reps = 1000,
  n_cores = -1,
  dgp_architecture = "mean_moderation"  # or "mvn"
)
```

## Installed Package

This collection is also the canonical `R/` directory of the installed pmsimstats package. Both implementations (`original` and `original-extended`) are available in a single package installation via the `dgp_architecture` parameter.

## Cross-Reference

- **original:** Architecture B only (baseline reference)
- **tidyverse:** Tidyverse reimplementation, Architecture A + B (alternative coding style)

## All Seven DGP Corrections Applied

Same corrections as `original` (see that README for details).
