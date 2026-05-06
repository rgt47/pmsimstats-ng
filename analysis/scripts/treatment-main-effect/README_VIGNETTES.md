# Notes on Running the Vignettes

This document provides guidance on how to run the package vignettes with the refactored code.

## Background

The package has been fully refactored to use consistent `snake_case` naming for all functions and variables. The vignettes in this directory have been updated to work with the refactored code, but may require some additional configuration.

## Running the Vignettes

1. **Direct Package Install:** The best way to ensure the vignettes run correctly is to install the package first:

```r
# From the package directory
devtools::install()
# Then run a specific vignette
rmarkdown::render("vignettes/design_comparison_extended.Rmd")
```

2. **Using `load_all()`:** If you prefer to use `load_all()` instead of installing the package, make sure to use `export_all = TRUE`:

```r
# In your vignette setup chunk
library(devtools)
# When running from vignettes/ directory
load_all(".", export_all = TRUE)  # The '.' refers to the package root

# When testing in a different directory, use the absolute path
# load_all("/path/to/pmsimstats", export_all = TRUE)
```

3. **Troubleshooting Function Access:** If you encounter "function not found" errors:

```r
# Add this to your vignette
# List all available functions in the global environment
cat(paste(ls(envir = .GlobalEnv)[grep("^[a-z]", ls(envir = .GlobalEnv))], collapse = "\n"))
```

## Common Issues

1. **"function not found" errors**: This typically happens when functions aren't properly exported to the global environment. Ensure you're using `load_all(., export_all = TRUE)`.

2. **Parameter naming issues**: All functions now use `snake_case` parameter names. Check function calls to ensure they match the refactored parameter names.

3. **Parallel Processing Function Access**: The most common issue is functions not being available to parallel worker processes. Use the newly added `setup_vignette_parallel()` function:

```r
# In your vignette setup
setup_vignette_parallel(
  load_full_package = TRUE,
  functions = c(
    "generate_data", 
    "generate_simulated_results", 
    "lme_analysis", 
    "mod_gompertz",
    "cumulative",
    "apply_carryover_effects",
    "calculate_carryover"
  ),
  verbose = TRUE
)
```

4. **Manual Global Function Export**: If `setup_vignette_parallel()` is not available, manually export functions:

```r
# Create globals list with all necessary functions
globals_list <- list(
  generate_data = generate_data,
  generate_simulated_results = generate_simulated_results,
  lme_analysis = lme_analysis,
  mod_gompertz = mod_gompertz,
  cumulative = cumulative,
  apply_carryover_effects = apply_carryover_effects,
  calculate_carryover = calculate_carryover
)

# Set up basic multisession plan
future::plan(future::multisession)

# Increase globals size limit
options(future.globals.maxSize = 1000 * 1024^2)  # 1GB

# Make functions available to global environment for workers to find
for (func_name in names(globals_list)) {
  if (!is.null(globals_list[[func_name]])) {
    assign(func_name, globals_list[[func_name]], envir = .GlobalEnv)
  }
}
```

5. **Dependencies**: The vignettes rely on several packages. Make sure they're all installed:

```r
pkgs <- c("tidyverse", "future", "furrr", "progressr", "knitr", 
          "kableExtra", "ggplot2", "scales", "rlang", "lme4", 
          "lmerTest", "MASS", "corpcor")
install.packages(pkgs)
```

## Using Parallel Processing in Vignettes

The package now includes optimized utilities for parallel processing:

1. **For Vignettes**: Use `setup_vignette_parallel()` to ensure all necessary functions are made available to worker processes.

2. **For Custom Scripts**: Use `setup_parallel_processing()` to configure optimal parallel processing.

3. **Manual Configuration**: When needed, you can still manually configure future/furrr with globals:

```r
# For explicit control - proper approach
# Set up basic multisession plan with desired workers
future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))

# Increase globals size limit
options(future.globals.maxSize = 1000 * 1024^2)  # 1GB

# Create list of functions to export
globals_list <- list(
  generate_data = generate_data,
  generate_simulated_results = generate_simulated_results,
  # other functions needed
)

# Make functions available to global environment for workers to find
for (func_name in names(globals_list)) {
  if (!is.null(globals_list[[func_name]])) {
    assign(func_name, globals_list[[func_name]], envir = .GlobalEnv)
  }
}
```

## Changes to Note

1. All function names are in `snake_case`: `generate_data`, `run_carryover_sensitivity`, etc.
2. All parameter names are in `snake_case`: `model_param` (not `modelParam`)
3. All variable names are in `snake_case`: `participant_id` (not `ptID`), `on_drug` (not `onDrug`)
4. New parallel utilities provide improved function accessibility in parallel processes

These changes make the code more consistent with R style guidelines and improve maintainability.

## Troubleshooting

If you still encounter issues:

1. Check that all required packages are installed
2. Check for the use of any deprecated camelCase function names or parameters 
3. Try running with `verbose=TRUE` to see detailed output about function loading
4. Use `load_all(".", export_all = TRUE)` to ensure all package functions are available
5. For parallel processing issues, make functions available to worker processes by assigning them to the global environment: `future::plan(future::multisession); for (name in names(globals)) { assign(name, globals[[name]], envir = .GlobalEnv) }`

## Optimizing Monte Carlo Vignette Rendering

The `monte_carlo_design_comparison.Rmd` vignette performs Monte Carlo simulations that can be time-consuming. Several optimizations have been implemented to speed up rendering:

1. **Reduced Parameter Grid**: Uses a minimal parameter grid for demonstration purposes.

2. **Caching**: All computation-heavy chunks have caching enabled with `cache=TRUE`.

3. **Pre-computed Results**: Can use pre-computed simulation results to skip time-consuming computations.
   - Run `/vignettes/temp/create_precomputed_script.sh` once to generate cached data
   - Set `use_cached_results <- TRUE` in the vignette setup chunk to use pre-computed data

4. **Optimized Visualizations**: Uses simplified data subsets and efficient plotting settings.

5. **Parallel Processing**: Automatically configures optimal parallel processing based on available cores.

If the vignette rendering is taking too long or hitting LaTeX memory limits, try these steps:
1. Use the pre-computed data by running the script in the temp/ directory
2. Reduce the parameter grid by modifying the parameter ranges
3. Reduce the number of Monte Carlo iterations by setting base_params$n_iterations to a smaller value