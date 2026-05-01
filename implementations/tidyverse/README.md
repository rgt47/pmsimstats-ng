# Tidyverse Implementation (Architecture A + B)

## Collection Label: `tidyverse`

A complete tidyverse-style reimplementation of the pmsimstats simulation framework, supporting both DGP architectures via the `dgp_architecture` parameter.

## Architecture Support

- **Architecture A (Mean Moderation):** Direct mean scaling of biomarker response (default `dgp_architecture = "mean_moderation"`)
- **Architecture B (MVN Differential Correlation):** Biomarker-response interaction via covariance structure (`dgp_architecture = "mvn"`)

## Coding Style

- Tidyverse packages: `tibble`, `dplyr`, `purrr`
- Native R pipe: `|>`
- Snake_case naming conventions
- No data.table usage

## File Structure

- `R/functions.R` — Complete simulation functions (tidyverse implementations)
- `R/test-alignment.R` — Cross-validation against `original` collection

## Key Functions

| Function | Purpose |
|---|---|
| `generate_data()` | Data generation with dual architecture support |
| `build_sigma_matrix()` | Covariance matrix construction |
| `generate_simulated_results()` | Batch simulation orchestration |
| `lme_analysis()` | LME analysis model with nlme::lme |
| `build_trial_design()` | Trial design specification |

## Usage Example

```r
source("implementations/tidyverse/R/functions.R")

# Build trial design
td <- build_trial_design(n_reps = 500, design = "Hybrid", carryover_t1half = 1.0)

# Define parameters
model_param <- list(N = 50, c.bm = 0.45, carryover_t1half = 1.0)
resp_param <- tibble::tibble(
  cat = rep(c("tv", "pb", "br"), each = 1),
  max = c(50, 20, 30),
  disp = c(4, 3, 2),
  rate = c(0.5, 0.3, 0.4),
  sd = c(10, 8, 5)
)
baseline_param <- tibble::tibble(
  cat = c("bm", "BL"),
  m = c(120, 60),
  sd = c(15, 10)
)

# Generate data (Architecture A - mean moderation)
dat <- generate_data(
  model_param, resp_param, baseline_param, td,
  empirical = FALSE, make_positive_definite = TRUE,
  dgp_architecture = "mean_moderation"
)

# Or use Architecture B (MVN)
dat_mvn <- generate_data(
  model_param, resp_param, baseline_param, td,
  empirical = FALSE, make_positive_definite = TRUE,
  dgp_architecture = "mvn"
)
```

## Validation

This implementation is validated to produce identical results to `implementations/original/` under `dgp_architecture = "mvn"` via `test-alignment.R`.

## Cross-Reference

- **original:** Data.table implementation, Architecture B only (corrected Hendrickson publication code)
- **original-extended:** Data.table implementation, Architecture A + B (feature branch code)

## Parameters

See `implementations/original/README.md` for complete parameter documentation. The tidyverse collection accepts the same parameters as the original implementations, with the same semantics for `c.bm`, `carryover_t1half`, and `lambda_cor`.
