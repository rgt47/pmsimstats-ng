# pmsimstats Implementations

This directory contains three parallel, complete implementations of the pmsimstats-ng simulation framework. Each collection supports different DGP architectures and uses different coding styles.

## Quick Reference

| Collection | Architecture | Style | Use Case |
|---|---|---|---|
| **original/** | B (MVN) only | data.table | Historical reference, backward compatibility baseline |
| **original-extended/** | A + B (dual) | data.table | Production use, default choice |
| **tidyverse/** | A + B (dual) | tidyverse | Modern alternative, type-stable code |

## Architecture Definitions

**Architecture A (Mean Moderation):** Biomarker directly scales the magnitude of drug response via additive mean shift. Power is preserved under carryover (8-13% loss). Appropriate when the biomarker determines the effective dose or biological responsiveness of each individual.

**Architecture B (MVN Differential Correlation):** Biomarker-response interaction encoded in the covariance structure. Power drops substantially under carryover (40-60% loss). Appropriate when the biomarker indexes a drug-responsive physiological state or subtype.

## Which One Should I Use?

### For new work
Use **`original-extended/`** — it's the production-ready implementation with both architectures available via the `dgp_architecture` parameter.

### For comparing coding styles
Use **`tidyverse/`** to see a modern tidyverse alternative while still having dual architecture support.

### For historical reference or backward compatibility
Use **`original/`** to match the exact corrected Hendrickson publication code (Architecture B only).

## Common Workflows

### Single trial simulation (Architecture A)
```r
source("implementations/original-extended/R/generateData.R")
source("implementations/original-extended/R/lme_analysis.R")

dat <- generateData(
  modelparam, respparam, blparam, trial_design,
  empirical = FALSE, makePositiveDefinite = TRUE,
  dgp_architecture = "mean_moderation"
)
```

### Batch simulation (both architectures)
```r
results_a <- generateSimulatedResults(
  ..., dgp_architecture = "mean_moderation"
)

results_b <- generateSimulatedResults(
  ..., dgp_architecture = "mvn"
)
```

### Tidyverse alternative
```r
source("implementations/tidyverse/R/functions.R")

dat <- generate_data(
  model_param, resp_param, baseline_param, trial_design,
  empirical = FALSE, make_positive_definite = TRUE,
  dgp_architecture = "mean_moderation"
)
```

## Core Parameters (All Collections)

| Parameter | Type | Meaning |
|---|---|---|
| `c.bm` | numeric | Under MVN: BM-BR correlation (0-1). Under mean moderation: regression coefficient scaling BR by centered biomarker |
| `carryover_t1half` | numeric | Half-life of carryover effect (weeks). 0 = no carryover |
| `lambda_cor` | numeric | Correlation decay rate (per week). Auto-derived as ln(2)/carryover_t1half if NA |
| `dgp_architecture` | character | `"mvn"` (default, Architecture B) or `"mean_moderation"` (Architecture A). Not present in `original/` |

## Architecture Support Matrix

| Collection | dgp_architecture parameter | Default | Architecture A | Architecture B |
|---|---|---|---|---|
| original/ | Not available | MVN only | ✗ | ✓ |
| original-extended/ | Available | mvn | ✓ | ✓ |
| tidyverse/ | Available | mvn | ✓ | ✓ |

## File Organization

Each collection contains:
```
<collection>/
  R/                  # Core simulation functions
    generateData.R
    generateSimulatedResults.R
    buildtrialdesign.R
    lme_analysis.R
    utilities.R
    censordata.R
    carryover_analysis.R
  README.md           # Collection-specific documentation
```

## Validation & Testing

- `implementations/tidyverse/R/test-alignment.R` validates that tidyverse and original collections produce identical results under `dgp_architecture = "mvn"`
- `implementations/original-extended/` produces identical results to `original/` when `dgp_architecture = "mvn"`
- Cross-architecture validation: both implementations produce same power estimates under both architectures

## Differences from the Main R/ Package

The `R/` directory at the repository root is the installable pmsimstats package (DESCRIPTION, NAMESPACE, roxygen documentation). It tracks `original-extended/` (dual architecture support).

The `implementations/` collections are standalone reference implementations without package infrastructure. Use `source()` to load functions directly.

## See Also

- Each collection's **README.md** for detailed parameter documentation and usage examples
- **CLAUDE.md** for project-wide conventions and guidance
- **docs/02-dgp-mean-moderation-vs-mvn.md** for theoretical comparison of the two architectures
