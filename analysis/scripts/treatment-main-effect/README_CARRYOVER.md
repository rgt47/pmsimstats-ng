# Carryover Sensitivity Analysis Vignette

This directory contains files related to the carryover sensitivity analysis in clinical trial simulations.

## Key Files

- `carryover_sensitivity_analysis.Rmd`: The main vignette demonstrating carryover model comparisons
- `carryover_only.Rmd`: A simpler vignette focusing only on carryover effects in the Hybrid N-of-1 design
- `vig_carry.R`: A monolithic script containing the code from the carryover_only vignette
- `vig5.R`: A monolithic script containing the code from the monte_carlo_design_comparison vignette

## Rendering the Vignettes

To render the vignettes, use:

```r
rmarkdown::render("vignettes/carryover_sensitivity_analysis.Rmd")
rmarkdown::render("vignettes/carryover_only.Rmd")
```

## Modifications Made

1. **Removed pmsimstats package references**: The vignettes now directly source the required R functions instead of loading the pmsimstats package:
   - `source("../R/carryover_models.R")`
   - `source("../R/carryover_sensitivity.R")`
   - `source("../R/math_utilities.R")`
   - `source("../R/utility_functions.R")`
   - `source("../R/logging_utilities.R")`

2. **Fixed plot implementations**: Modified parameter effects plot implementations to avoid errors with the `effect` variable in ggplot.

3. **Disabled caching**: Set `cache = FALSE` in the setup chunk to avoid dependency issues.

4. **Improved performance**: Reduced number of iterations and parameter combinations for faster rendering.

## Monolithic Script

The `vig_carry.R` script contains all the code from the carryover vignette in a single standalone file. 
This allows you to run the analysis without having to render the entire vignette.

To run it in RStudio:

1. Open `vig_carry.R` in RStudio
2. Use "Run All" or select portions to run

This approach is suitable for quick testing or when you want to modify parameters without re-rendering the full vignette.