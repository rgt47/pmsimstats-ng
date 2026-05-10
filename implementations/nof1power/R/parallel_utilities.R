#' Create global function list for parallel processing
#'
#' This function creates a list of global functions that are needed in 
#' parallel processes. It ensures that all necessary functions are available
#' to workers in future/furrr parallel processing.
#'
#' @param functions Vector of function names to export as globals
#' @param packages Vector of package names to load in workers
#' @param add_defaults Whether to include a default set of commonly used functions
#' @return List of functions for use with future::plan globals parameter
#' @importFrom rlang env_get env_has
#' @export
prepare_parallel_globals <- function(
  functions = NULL,
  packages = NULL,
  add_defaults = TRUE
) {
  globals_list <- list()
  
  # Get current environment and package namespace
  pkg_env <- asNamespace("nof1power")
  
  # Add requested functions
  if (!is.null(functions)) {
    for (func_name in functions) {
      # First try to get from package namespace
      if (exists(func_name, envir = pkg_env, inherits = FALSE)) {
        globals_list[[func_name]] <- get(func_name, envir = pkg_env)
      # Then try from global environment
      } else if (exists(func_name, envir = .GlobalEnv, inherits = FALSE)) {
        globals_list[[func_name]] <- get(func_name, envir = .GlobalEnv)
      # If not found, try to get from imported packages
      } else {
        found <- FALSE
        for (package in packages) {
          if (requireNamespace(package, quietly = TRUE)) {
            pkg_namespace <- asNamespace(package)
            if (exists(func_name, envir = pkg_namespace, inherits = FALSE)) {
              globals_list[[func_name]] <- get(func_name, envir = pkg_namespace)
              found <- TRUE
              break
            }
          }
        }
        if (!found) {
          warning("Function '", func_name, "' not found in package namespace or global environment.")
          globals_list[[func_name]] <- NULL
        }
      }
    }
  }
  
  # Add default functions used in simulation
  if (add_defaults) {
    default_functions <- c(
      "generate_data",
      "generate_simulated_results",
      "lme_analysis",
      "mod_gompertz",
      "apply_carryover_effects",
      "calculate_carryover",
      "log_progress",
      "logging_progress_bar",
      "update_logging_progress_bar",
      "close_logging_progress_bar",
      "cumulative"
    )
    
    for (func_name in default_functions) {
      # Skip if already added
      if (func_name %in% names(globals_list)) {
        next
      }
      
      # Try to get from package namespace
      if (exists(func_name, envir = pkg_env, inherits = FALSE)) {
        globals_list[[func_name]] <- get(func_name, envir = pkg_env)
      } else if (exists(func_name, envir = .GlobalEnv, inherits = FALSE)) {
        globals_list[[func_name]] <- get(func_name, envir = .GlobalEnv)
      }
    }
  }
  
  return(globals_list)
}

#' Set up parallel processing for simulations
#'
#' This function configures parallel processing with optimal settings
#' for simulation runs, ensuring all necessary functions are available.
#'
#' @param workers Number of worker processes to use (default: auto-detect)
#' @param additional_functions Vector of additional function names to include
#' @param additional_packages Vector of additional package names to include
#' @param strategy Parallel strategy ("multisession" or "multicore")
#' @return Invisibly returns the future plan object
#' @export
setup_parallel_processing <- function(
  workers = NULL,
  additional_functions = NULL,
  additional_packages = NULL,
  strategy = "multisession"
) {
  # Auto-detect number of cores if not specified
  if (is.null(workers)) {
    workers <- max(1, parallel::detectCores() - 1)
  }
  
  # Combine default and additional functions
  functions <- c(
    "generate_data",
    "generate_simulated_results",
    "lme_analysis",
    "mod_gompertz",
    "apply_carryover_effects",
    "calculate_carryover",
    "cumulative",
    "log_progress",
    "logging_progress_bar",
    "update_logging_progress_bar",
    "close_logging_progress_bar",
    additional_functions
  )
  
  # Combine default and additional packages
  packages <- c(
    "nof1power",
    "dplyr",
    "tidyr",
    "purrr",
    additional_packages
  )
  
  # Prepare globals list
  globals_list <- prepare_parallel_globals(
    functions = functions,
    packages = packages,
    add_defaults = TRUE
  )
  
  # Set up parallel processing plan
  parallel_strategy <- match.arg(
    strategy, 
    c("multisession", "multicore"),
    several.ok = FALSE
  )
  
  # Choose appropriate plan based on platform
  if (parallel_strategy == "multicore" && .Platform$OS.type == "windows") {
    # Multicore not available on Windows, fallback to multisession
    parallel_strategy <- "multisession"
  }

  # Just set up a basic plan - we'll handle globals differently
  # Use the appropriate future strategy based on name
  if (parallel_strategy == "multisession") {
    plan_obj <- future::plan(future::multisession, workers = workers)
  } else if (parallel_strategy == "multicore") {
    plan_obj <- future::plan(future::multicore, workers = workers)
  } else {
    # Default to sequential if unknown
    plan_obj <- future::plan(future::sequential)
  }

  # Increase globals size limit
  options(future.globals.maxSize = 1000 * 1024^2)  # 1GB

  # Make globals directly available in global environment
  # This allows worker processes to find them automatically
  for (func_name in names(globals_list)) {
    if (!is.null(globals_list[[func_name]])) {
      assign(func_name, globals_list[[func_name]], envir = .GlobalEnv)
    }
  }
  
  # Return plan object invisibly
  invisible(plan_obj)
}

#' Configure vignette for parallel processing
#'
#' This function is specifically designed to be called at the top of vignettes
#' to ensure all needed functions are available to workers when using
#' future/furrr for parallel processing.
#'
#' @param load_full_package Whether to attempt a devtools::load_all() operation
#' @param functions Vector of function names that should be made available to workers
#' @param verbose Whether to print status messages
#' @return Invisibly returns the future plan object
#' @export
setup_vignette_parallel <- function(
  load_full_package = TRUE,
  functions = NULL,
  verbose = TRUE
) {
  # First check if the package is already loaded
  package_loaded <- "nof1power" %in% (.packages())
  
  # Try to load the package using devtools if requested
  if (load_full_package && !package_loaded) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      if (verbose) cat("Attempting to load package with devtools::load_all()...\n")
      tryCatch({
        devtools::load_all(".", export_all = TRUE)
        if (verbose) cat("Successfully loaded package with devtools::load_all()\n")
        package_loaded <- TRUE
      }, error = function(e) {
        if (verbose) cat("Failed to load with devtools::load_all(): ", conditionMessage(e), "\n")
      })
    }

    if (!package_loaded) {
      if (verbose) cat("Trying to load the installed package instead...\n")
      tryCatch({
        library(nof1power)
        if (verbose) cat("Successfully loaded nof1power library\n")
        package_loaded <- TRUE
      }, error = function(e) {
        if (verbose) cat("Failed to load nof1power: ", conditionMessage(e), "\n")
        if (verbose) cat("Will continue with limited functionality\n")
      })
    }
  }
  
  # Default functions to export to workers
  default_functions <- c(
    "generate_data",
    "generate_simulated_results",
    "lme_analysis",
    "mod_gompertz",
    "apply_carryover_effects",
    "calculate_carryover",
    "cumulative",
    "log_progress",
    "logging_progress_bar",
    "update_logging_progress_bar",
    "close_logging_progress_bar"
  )
  
  # Combine with user-specified functions
  all_functions <- c(default_functions, functions)
  
  # Create globals list
  globals_list <- list()
  
  # First try to get functions from package namespace
  if (package_loaded) {
    pkg_env <- asNamespace("nof1power")
    for (func_name in all_functions) {
      if (exists(func_name, envir = pkg_env, inherits = FALSE)) {
        globals_list[[func_name]] <- get(func_name, envir = pkg_env)
      }
    }
  }
  
  # Then get any remaining functions from global environment
  for (func_name in all_functions) {
    if (!(func_name %in% names(globals_list)) && 
        exists(func_name, envir = .GlobalEnv, inherits = FALSE)) {
      globals_list[[func_name]] <- get(func_name, envir = .GlobalEnv)
    }
  }
  
  # Verify which functions were found
  if (verbose) {
    found_funcs <- names(globals_list)
    missing_funcs <- setdiff(all_functions, found_funcs)
    
    cat("Found", length(found_funcs), "functions for parallel export.\n")
    if (length(missing_funcs) > 0) {
      cat("Warning: Could not find these functions:", 
          paste(missing_funcs, collapse=", "), "\n")
    }
  }
  
  # Set up parallel processing with these globals
  workers <- max(1, parallel::detectCores() - 1)
  
  # Just set up a basic plan - we'll handle globals differently
  plan_obj <- future::plan(future::multisession, workers = workers)

  # Increase globals size limit
  options(future.globals.maxSize = 1000 * 1024^2)  # 1GB

  # Make globals directly available in global environment
  # This allows worker processes to find them automatically
  for (func_name in names(globals_list)) {
    if (!is.null(globals_list[[func_name]])) {
      assign(func_name, globals_list[[func_name]], envir = .GlobalEnv)
    }
  }
  
  if (verbose) {
    cat("Parallel processing configured with", workers, "workers.\n")
  }
  
  # Return the plan object invisibly
  invisible(plan_obj)
}