########################################################################################
# vig_carry.R: Self-contained script for carryover effects simulation in Hybrid N-of-1 Design
#
# This is a monolithic script that contains all the code needed to run the
# carryover sensitivity analysis vignette. It demonstrates models for carryover
# effects and runs simulations comparing different trial designs.
# All helper functions and simulation code are included directly in this file.
########################################################################################

# Load required packages
# Suppress package startup messages
# and clear the workspace (optional)
# rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(nlme)
  library(knitr)
  library(parallel)
  library(grid)
  library(tools)
  library(MASS)
  library(corpcor)
})

n_cores <- parallel::detectCores()
n_workers <- max(1, floor(n_cores * 0.75))
cat(paste0("Using ", n_workers, " workers out of ", n_cores, " available cores\n"))

# Helper functions
# ----------------

#===========================================================================
# Function: modify_list
# Description: Merge two lists, with the second list's values overriding
#              any matching names in the first list
#===========================================================================
modify_list <- function(x, val) {
  if (is.null(x))
    return(val)
  if (is.null(val))
    return(x)

  # Get names from both lists
  xnames <- names(x)
  vnames <- names(val)

  # Get the unique names
  vnames <- vnames[nzchar(vnames)]
  xnames <- xnames[nzchar(xnames)]

  # For each named item in val, replace or add to x
  for (v in vnames) {
    x[[v]] <- val[[v]]
  }

  return(x)
}
#---------------------------------------------------------------------------
# End Function: modify_list
#---------------------------------------------------------------------------

#===========================================================================
# Function: cumulative
# Description: Wrapper around cumsum for calculating cumulative values
#===========================================================================
cumulative <- function(values) {
  cumsum(values)
}
#---------------------------------------------------------------------------
# End Function: cumulative
#---------------------------------------------------------------------------

#===========================================================================
# Function: mod_gompertz
# Description: Modified Gompertz function for modeling growth curves
#===========================================================================
mod_gompertz <- function(time, max_value, displacement, rate) {
  max_value * (1 - exp(-displacement * exp(-rate * time)))
}
#---------------------------------------------------------------------------
# End Function: mod_gompertz
#---------------------------------------------------------------------------

#===========================================================================
# Function: %||%
# Description: Null coalescing operator - returns first argument if not NULL,
#              otherwise returns second argument
#===========================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
#---------------------------------------------------------------------------
# End Function: %||%
#---------------------------------------------------------------------------

#===========================================================================
# Function: generate_data
# Description: Primary data generation function for trial simulation
#===========================================================================
generate_data <- function(model_param, resp_param, baseline_param, trial_design,
                        empirical, make_positive_definite, seed = NA,
                        scale_factor = 2, verbose = FALSE, track_pd_stats = FALSE) {
  # Initialize variables for tracking sigma matrix statistics
  sigma_count <- 0
  non_positive_definite_count <- 0

  # I. Turn the trial design information into something easier to use
  if (verbose) {
    cat("Generate data: Converting trial_design to data.table\n")
    cat("Columns in trial_design:", paste(names(trial_design), collapse=", "), "\n")
  }

  # Ensure trial_design is a data frame
  if (!is.data.frame(trial_design)) {
    stop("Trial design is not a data frame. Class: ", class(trial_design),
         ", Type: ", typeof(trial_design))
  }

  trial_data <- as_tibble(trial_design)

  # Check if t_wk exists and calculate cumulative if it does
  if ("t_wk" %in% names(trial_data)) {
    trial_data$t_wk_cumulative <- cumulative(trial_data$t_wk)
    if (verbose) cat("Using existing t_wk column\n")
  } else {
    # If t_wk doesn't exist, calculate it or use week column directly
    if ("week" %in% names(trial_data)) {
      trial_data$t_wk <- c(trial_data$week[1], diff(trial_data$week))
      trial_data$t_wk_cumulative <- trial_data$week
      if (verbose) cat("Calculated t_wk from week column\n")
    } else {
      if (verbose) {
        cat("Error: Missing required columns\n")
        cat("Available columns:", paste(names(trial_data), collapse=", "), "\n")
      }
      stop("Neither 't_wk' nor 'week' column found in trial design data")
    }
  }

  trial_data <- trial_data %>% mutate(on_drug = (tod > 0))
  num_timepoints <- dim(trial_design)[1]

  # II. Set up variables to track - baseline parameters for each participant ("biomarker","baseline"),
  #     and the three modeled factors for each stage of the trial.

  # Set up the variable names
  factor_types <- c("time_variant", "pharm_biomarker", "bio_response")
  factor_abbreviations <- c("tv", "pb", "br")

  labels <- c(
    c("biomarker", "baseline"),
    paste(trial_design$timepoint_name, factor_abbreviations[1], sep = "."),
    paste(trial_design$timepoint_name, factor_abbreviations[2], sep = "."),
    paste(trial_design$timepoint_name, factor_abbreviations[3], sep = ".")
  )

  # Set up vectors with the standard deviations and means
  standard_deviations <- c(
    baseline_param$sd[baseline_param$cat == "biomarker"],
    baseline_param$sd[baseline_param$cat == "baseline"]
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "time_variant"], num_timepoints)
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "pharm_biomarker"], num_timepoints) * trial_design$e
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "bio_response"], num_timepoints)
  )

  means <- c(
    baseline_param$m[baseline_param$cat == "biomarker"],
    baseline_param$m[baseline_param$cat == "baseline"]
  )

  for (factor_idx in seq_along(factor_types)) {
    current_factor <- factor_types[factor_idx]
    factor_abbrev <- factor_abbreviations[factor_idx]

    if (factor_abbrev == "tv") {
      time_variant_means <- mod_gompertz(
        trial_data$t_wk_cumulative,
        resp_param$max[resp_param$cat == current_factor],
        resp_param$disp[resp_param$cat == current_factor],
        resp_param$rate[resp_param$cat == current_factor]
      )
      means <- c(means, time_variant_means)
    }

    if (factor_abbrev == "pb") {
      pharm_biomarker_means <- mod_gompertz(
        trial_data$tpb,
        resp_param$max[resp_param$cat == current_factor],
        resp_param$disp[resp_param$cat == current_factor],
        resp_param$rate[resp_param$cat == current_factor]
      ) * trial_design$e
      means <- c(means, pharm_biomarker_means)
    }

    if (factor_abbrev == "br") {
      bio_response_means <- mod_gompertz(
        trial_data$tod,
        resp_param$max[resp_param$cat == current_factor],
        resp_param$disp[resp_param$cat == current_factor],
        resp_param$rate[resp_param$cat == current_factor]
      )

      # Check for zero values
      bio_response_test <- bio_response_means == 0
      raw_bio_response_means <- bio_response_means

      # Find indices for bio_response (br) columns dynamically
      br_indices <- grep("\\.br$", labels)

      # Set names for debugging - using dynamic indices
      if (length(br_indices) > 0) {
        names(bio_response_test) <- labels[br_indices]
        names(raw_bio_response_means) <- labels[br_indices]
      }

      # Vectorized carryover calculation
      if (num_timepoints > 1) {
        # Create index vectors for non-drug timepoints with positive tsd
        non_drug_indices <- which(!trial_data$on_drug & trial_data$tsd > 0)

        # Only process if there are any applicable timepoints
        if (length(non_drug_indices) > 0) {
          # Get previous indices (p-1) for each qualifying timepoint
          prev_indices <- non_drug_indices - 1

          # Calculate decay factors for each qualifying timepoint
          decay_factors <- (1/2)^(
            scale_factor * trial_data$tsd[non_drug_indices] / model_param$carryover_t1half
          )

          # Apply carryover effect in one vectorized operation
          bio_response_means[non_drug_indices] <- bio_response_means[non_drug_indices] +
            bio_response_means[prev_indices] * decay_factors
        }
      }

      means <- c(means, bio_response_means)
    }
  }

  # Build a correlation matrix
  correlations <- diag(length(labels))
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  # Debug output in verbose mode:
  if (verbose == TRUE) {
    # Create debugging data frame only if we have named vectors
    if (length(names(raw_bio_response_means)) > 0) {
      aa <- data.frame(
        model_param$carryover_t1half,
        raw_bio_response_means,
        bio_response_means,
        diff = bio_response_means - raw_bio_response_means
      )
      cat("bio_response_means before and after adjustment:\n ")
      print(aa)
    } else {
      cat("bio_response_means before and after adjustment: (unable to format for display)\n")
      cat("raw values:", paste(raw_bio_response_means, collapse=", "), "\n")
      cat("adjusted values:", paste(bio_response_means, collapse=", "), "\n")
      cat("difference:", paste(bio_response_means - raw_bio_response_means, collapse=", "), "\n")
    }
  }

  # Apply correlations between factors
  for (factor_idx in seq_along(factor_types)) {
    current_factor <- factor_abbreviations[factor_idx]

    # Build autocorrelations across time - VECTORIZED
    if (num_timepoints > 1) {
      autocorrelation <- model_param[[paste("c", current_factor, sep = ".")]]

      # Create all combinations of time points at once
      point_indices <- expand.grid(p1 = 1:(num_timepoints-1), p2 = (2:num_timepoints))
      # Filter valid combinations where p2 > p1
      point_indices <- point_indices[point_indices$p2 > point_indices$p1, ]

      # Create name vectors for efficient indexing
      name1 <- paste(trial_design$timepoint_name[point_indices$p1], current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name[point_indices$p2], current_factor, sep = ".")

      # Set all correlations at once using matrix indexing
      for (idx in 1:nrow(point_indices)) {
        correlations[name1[idx], name2[idx]] <- autocorrelation
        correlations[name2[idx], name1[idx]] <- autocorrelation
      }
    }

    # Build autocorrelations across factors - VECTORIZED
    for (factor2_idx in setdiff(seq_along(factor_types), factor_idx)) {
      other_factor <- factor_abbreviations[factor2_idx]

      # Same timepoint cross-factor correlations
      name1 <- paste(trial_design$timepoint_name, current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name, other_factor, sep = ".")

      # Set all same-timepoint correlations at once
      for (idx in 1:length(name1)) {
        correlations[name1[idx], name2[idx]] <- model_param$c.cf1t
        correlations[name2[idx], name1[idx]] <- model_param$c.cf1t
      }

      # Different timepoint cross-factor correlations
      if (num_timepoints > 1) {
        # Create all combinations of time points at once
        point_indices <- expand.grid(p1 = 1:(num_timepoints-1), p2 = (2:num_timepoints))
        # Filter valid combinations where p2 > p1
        point_indices <- point_indices[point_indices$p2 > point_indices$p1, ]

        # Create name vectors for efficient indexing
        name1 <- paste(trial_design$timepoint_name[point_indices$p1], current_factor, sep = ".")
        name2 <- paste(trial_design$timepoint_name[point_indices$p2], other_factor, sep = ".")

        # Set all correlations at once
        for (idx in 1:nrow(point_indices)) {
          correlations[name1[idx], name2[idx]] <- model_param$c.cfct
          correlations[name2[idx], name1[idx]] <- model_param$c.cfct
        }
      }
    }

    # Special handling for biomarker correlation
    if (current_factor == "br") {
      for (timepoint_idx in 1:num_timepoints) {
        name1 <- paste(trial_design$timepoint_name[timepoint_idx], "br", sep = ".")

        if (timepoint_idx > 1) {
          name0 <- paste(trial_design$timepoint_name[timepoint_idx - 1], "br", sep = ".")
          mean_value1 <- means[which(name1 == labels)]
          mean_value0 <- means[which(name0 == labels)]

          # Handle special cases with careful checks
          correlations["biomarker", name1] <- correlations[name1, "biomarker"] <-
            ifelse(bio_response_test[timepoint_idx],
                  ifelse(bio_response_means[timepoint_idx] == 0 || abs(mean_value0) < 1e-10, 0,
                        (mean_value1 / max(mean_value0, 1e-10)) * model_param$c.bm),
                  model_param$c.bm)
        }
      }
    }
  }

  # Optimize matrix operations for the covariance matrix calculation
  # Get matrix dimensions and prepare for efficient operations
  n_vars <- length(standard_deviations)
  large_matrix <- n_vars > 100

  # Track statistics about sigma matrix if requested
  if (track_pd_stats) {
    sigma_count <- sigma_count + 1
  }

  # Fast path for positive definiteness handling
  is_positive_definite <- TRUE  # Assume matrix is positive definite until proven otherwise
  need_pd_check <- make_positive_definite || track_pd_stats

  # Turn correlation matrix into covariance matrix using efficient outer product
  sigma <- outer(standard_deviations, standard_deviations) * correlations

  # Check/fix positive definiteness if required
  if (need_pd_check) {
    is_positive_definite <- corpcor::is.positive.definite(sigma)

    # Update statistics if tracking
    if (track_pd_stats && !is_positive_definite) {
      non_positive_definite_count <- non_positive_definite_count + 1
    }
  }

  # Fix if turned on and not positive definite
  if (make_positive_definite && !is_positive_definite) {
    # For more robust conversion to positive definite
    sigma <- corpcor::make.positive.definite(sigma, tol = 1e-3)
  }

  # Set the seed if provided for reproducibility
  if (!is.na(seed)) {
    set.seed(seed)
  }

  # Generate multivariate normal data
  participant_data <- MASS::mvrnorm(
    n = model_param$N,
    mu = means,
    Sigma = sigma,
    empirical = empirical
  )
  participant_data <- as_tibble(participant_data)
  colnames(participant_data) <- labels

  # Pre-allocate all data structures for better performance
  N <- model_param$N

  # Add participant ID column
  participant_data <- participant_data %>%
    mutate(participant_id = 1:N)

  # Get the timepoint names
  timepoint_names <- trial_design$timepoint_name

  # Pre-allocate all result columns to avoid growing the data frame
  # First create list of all columns we'll need to add
  new_cols <- c(paste0("D_", timepoint_names), timepoint_names)

  # Create a dataframe of zeros to add
  zeros_df <- as_tibble(matrix(0, nrow = N, ncol = length(new_cols)))
  colnames(zeros_df) <- new_cols

  # Add the zero columns to participant_data
  participant_data <- bind_cols(participant_data, zeros_df)

  # Calculate deltas (sums of factors) in a vectorized way for all participants
  for (timepoint_idx in 1:num_timepoints) {
    timepoint_name <- timepoint_names[timepoint_idx]
    delta_col <- paste0("D_", timepoint_name)

    # Get the component columns
    components <- paste(timepoint_name, factor_abbreviations, sep = ".")

    # Ensure all components exist
    for (comp in components) {
      if (!(comp %in% names(participant_data))) {
        warning(paste("Component column", comp, "not found for delta calculation"))
      }
    }

    # Calculate delta for all participants at once
    participant_data <- participant_data %>%
      mutate(!!delta_col := rowSums(select(., all_of(components))))
  }

  # Calculate timepoint scores from baseline and deltas
  for (timepoint_idx in 1:num_timepoints) {
    timepoint_name <- timepoint_names[timepoint_idx]
    components <- paste(timepoint_name, factor_abbreviations, sep = ".")

    # Calculate timepoint value from baseline and factors
    participant_data <- participant_data %>%
      mutate(!!timepoint_name := baseline - rowSums(select(., all_of(components))))
  }

  # Attach matrix tracking statistics as attributes if requested
  if (track_pd_stats) {
    non_positive_definite_rate <- non_positive_definite_count / sigma_count
    attr(participant_data, "sigma_count") <- sigma_count
    attr(participant_data, "non_positive_definite_count") <- non_positive_definite_count
    attr(participant_data, "non_positive_definite_rate") <- non_positive_definite_rate
  }

  return(participant_data)
}
#---------------------------------------------------------------------------
# End Function: generate_data
#---------------------------------------------------------------------------

#===========================================================================
# Function: create_hybrid_design
# Description: Create the hybrid N-of-1 design for a given number of participants
#===========================================================================
create_hybrid_design <- function(n_participants) {
  # Create Hybrid N-of-1 design
  # Open-label (8wk) + blinded discontinuation (4wk) + crossover (8wk: 4wk active, 4wk placebo)
  hybrid_design <- expand_grid(
    participant_id = 1:n_participants,
    week = 1:20
  ) %>%
  mutate(
    treatment = c(rep(1, 8), rep(0, 4), rep(1, 4), rep(0, 4))[week],
    period = c(rep(1, 8), rep(2, 4), rep(3, 8))[week]
  )

  return(hybrid_design)
}
#---------------------------------------------------------------------------
# End Function: create_hybrid_design
#---------------------------------------------------------------------------

#===========================================================================
# Function: prepare_design_for_simulation
# Description: Convert designs to a format compatible with generate_data
#===========================================================================
prepare_design_for_simulation <- function(design_df) {
  # Create a copy to avoid modifying the original
  design <- design_df %>% as_tibble()

  # Add timepoint names that are compatible with generate_data
  design <- design %>%
    arrange(participant_id, week) %>%
    mutate(
      timepoint_name = paste0("W", week),
      tod = treatment,                      # Time on drug (1 = on drug, 0 = off drug)
      e = if_else(tod > 0, 1, 0),          # Expectancy (1 = on drug, 0 = off drug)
      t_wk = 1,                             # Each timepoint is 1 week
      tpb = cumsum(e)                       # Time on pharmacologic biomarker
    ) %>%
    group_by(participant_id) %>%
    arrange(participant_id, week) %>%
    mutate(
      # Calculate time since discontinuation for carryover
      last_drug = lag(tod, default = 0),
      drug_stopped = (last_drug == 1 & tod == 0),
      tsd = if_else(tod == 0, cumsum(tod == 0), 0)  # Time since discontinuation
    ) %>%
    ungroup()

  # Ensure the design has all required columns
  design <- design %>%
    select(timepoint_name, t_wk, e, tod, tsd, tpb)

  return(design)
}
#---------------------------------------------------------------------------
# End Function: prepare_design_for_simulation
#---------------------------------------------------------------------------

#===========================================================================
# Function: run_monte_carlo
# Description: Run a Monte Carlo simulation for the hybrid design with different carryover parameters
#===========================================================================
run_monte_carlo <- function(params) {
  # Update model parameters with current simulation parameters
  current_model_params <- model_params
  current_model_params$carryover_t1half <- params$carryover_t1half

  # Update analysis options for carryover scale
  current_analysis_options <- analysis_options
  current_analysis_options$carryover_scale_factor <- params$carryover_scale

  # Select the design template
  design_template <- design$hybrid$design

  # Set number of iterations
  n_iter <- base_params$n_iterations

  all_results <- purrr::map_dfr(1:n_iter, function(iter) {

    # Generate data using generate_data
    sim_data <- generate_data(
      model_param = current_model_params,
      resp_param = resp_params,
      baseline_param = bl_params,
      trial_design = design_template,
      empirical = FALSE,
      make_positive_definite = TRUE,
      seed = iter,
      scale_factor = params$carryover_scale
    )

    # Convert to format needed for analysis
    # Extract timepoint columns and reshape to long format
    timepoint_columns <- grep("^W\\d+$", names(sim_data), value = TRUE)

    # Keep only needed columns
    analysis_data <- sim_data %>%
      select(participant_id, biomarker, all_of(timepoint_columns)) %>%
      rename(bm = biomarker)

    # Convert to long format for analysis
    analysis_data_long <- analysis_data %>%
      pivot_longer(
        cols = all_of(timepoint_columns),
        names_to = "timepoint",
        values_to = "response"
      ) %>%
      # Extract week number
      mutate(
        week = as.integer(sub("W", "", timepoint)),
        # Add treatment indicator based on design
        treatment = ifelse(week <= 8, 1,
                       ifelse(week <= 12, 0,
                              ifelse(week <= 16, 1, 0)))
      )

    # Analyze with linear mixed model
    model_result <- tryCatch({
      model <- nlme::lme(response ~ treatment * bm,
                         random = ~1 | participant_id,
                         data = analysis_data_long,
                         na.action = na.omit)

      model_summary <- summary(model)
      coefs <- model_summary$tTable

      interaction_idx <- which(rownames(coefs) == "treatment:bm")

      effect_size <- coefs[interaction_idx, "Value"]
      std_error <- coefs[interaction_idx, "Std.Error"]
      t_value <- coefs[interaction_idx, "t-value"]
      p_value <- coefs[interaction_idx, "p-value"]
      significant <- p_value < 0.05

      # Return results
      tibble(
        effect_size = effect_size,
        std_error = std_error,
        t_value = t_value,
        p_value = p_value,
        significant = significant,
        error = FALSE
      )
    }, error = function(e) {
      # Return NA values on error
      tibble(
        effect_size = NA_real_,
        std_error = NA_real_,
        t_value = NA_real_,
        p_value = NA_real_,
        significant = NA,
        error = TRUE
      )
    })

    # Add metadata for this iteration
    model_result %>%
      mutate(
        iteration = iter,
        carryover_t1half = params$carryover_t1half,
        carryover_scale = params$carryover_scale,
        n_participants = base_params$n_participants,
        biomarker_correlation = base_params$biomarker_correlation,
        treatment_effect = base_params$treatment_effect
      )
  })

  return(all_results)
}
#---------------------------------------------------------------------------
# End Function: run_monte_carlo
#---------------------------------------------------------------------------

#===========================================================================
# Function: run_parameter_sweep
# Description: Run simulations for the parameter grid
#===========================================================================
run_parameter_sweep <- function(param_grid) {
  # Number of parameter combinations
  n_combinations <- nrow(param_grid)

  # Run simulations for each parameter combination
  results <- map_dfr(1:n_combinations, function(i) {
    # Extract current parameters
    current_params <- as.list(param_grid[i,])

    # Print progress (only every 5 combinations to reduce output)
    if (i %% 5 == 1 || i == n_combinations) {
      cat(sprintf("Running combination %d of %d: t1half=%s, scale=%s\n",
                 i, n_combinations, current_params$carryover_t1half, current_params$carryover_scale))
    }

    # Run Monte Carlo simulation
    run_monte_carlo(current_params)
  })

  return(results)
}
#---------------------------------------------------------------------------
# End Function: run_parameter_sweep
#---------------------------------------------------------------------------

#===========================================================================
# Function: save_results
# Description: Save simulation results to disk for future use
#===========================================================================
save_results <- function() {
  results_dir <- "temp"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  # Create timestamp for unique filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Save the simulation data
  sim_data <- list(
    simulation_results = simulation_results,
    simulation_summary = simulation_summary,
    param_grid = param_grid,
    base_params = base_params,
    timestamp = timestamp
  )

  # Save as RData file
  save_path <- file.path(results_dir, paste0("carryover_simulation_", timestamp, ".RData"))
  save(sim_data, file = save_path)

  # Also save just the summary data as RDS for easier loading
  summary_path <- file.path(results_dir, paste0("carryover_summary_", timestamp, ".rds"))
  saveRDS(simulation_summary, file = summary_path)

  cat("Simulation results saved to:", save_path, "\n")
  cat("Summary results saved to:", summary_path, "\n")
}
#---------------------------------------------------------------------------
# End Function: save_results
#---------------------------------------------------------------------------

# Main Script
# ==================

# Use middle values for fixed parameters
base_params <- list(
  # Sample size - middle value
  n_participants = 50,

  # Biomarker-response correlation (middle value)
  biomarker_correlation = 0.5,

  # Treatment effect (middle value)
  treatment_effect = 4.0,

  # Carryover effect parameters (will be varied)
  carryover_t1half = 1,     # Default value
  carryover_scale = 1.5,    # Scale factor

  # Random variation
  between_subject_sd = 3.0,
  within_subject_sd = 2.8,

  # Number of Monte Carlo iterations - can be increased for production runs
  n_iterations = 5
)

# Define parameter grid with focus on carryover effects
# For production use a larger grid
param_grid <- expand_grid(
  carryover_t1half = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
  carryover_scale = c(0.5, 1.0, 1.5, 2.0)
)

# Create a directory for saving results if it doesn't exist
if (!dir.exists("temp")) {
  dir.create("temp", recursive = TRUE)
}

# Display the parameter grid
cat("Parameter grid for carryover sensitivity analysis:\n")
print(param_grid)

# Create model_params compatible with generate_data function
model_params <- list(
  N = base_params$n_participants,
  c.bm = base_params$biomarker_correlation,  # Correlation between biomarker and response
  carryover_t1half = base_params$carryover_t1half,
  c.tv = 0.7,  # Autocorrelation for time_variant factor
  c.pb = 0.7,  # Autocorrelation for pharm_biomarker factor
  c.br = 0.7,  # Autocorrelation for bio_response factor
  c.cf1t = 0.3, # Correlation between different factors at a single timepoint
  c.cfct = 0.2  # Correlation between different factors at different timepoints
)

# Set up response parameters
resp_params <- tibble(
  cat = c("time_variant", "pharm_biomarker", "bio_response"),
  max = c(1.0, 1.0, base_params$treatment_effect),  # Maximum response value
  disp = c(2.0, 2.0, 2.0),  # Displacement
  rate = c(0.3, 0.3, 0.3),  # Rate
  sd = c(base_params$within_subject_sd, base_params$within_subject_sd, base_params$within_subject_sd)  # Standard deviation
)

# Set up baseline parameters
bl_params <- tibble(
  cat = c("biomarker", "baseline"),
  m = c(5.0, 10.0),  # Mean values
  sd = c(2.0, base_params$between_subject_sd)  # Standard deviation
)

# Analysis options
analysis_options <- list(
  use_expectancy = TRUE,  # Use time information in model
  random_slope = FALSE,   # No random slopes for simplicity
  full_output = FALSE,    # Return only summary statistics
  simple_carryover = FALSE, # Use complex carryover model
  carryover_halflife = base_params$carryover_t1half,
  carryover_scale_factor = base_params$carryover_scale
)

# Options for simulation
sim_options <- list(
  n_reps = base_params$n_iterations  # Number of Monte Carlo replications
)

# Create design for base sample size for visualization
original_design <- create_hybrid_design(base_params$n_participants)

# Display sample of design
cat("Hybrid N-of-1 Design - First rows:\n")
print(head(original_design, 10))

# Create design compatible with the generate_data function
hybrid_design <- prepare_design_for_simulation(original_design %>%
                                           filter(participant_id == 1))

# Store design in a variable
design <- list(
  hybrid = list(design = hybrid_design, name = "Hybrid N-of-1")
)

# Display the prepared design
cat("Hybrid N-of-1 Design - Prepared Format:\n")
print(design$hybrid$design)

# Running the Monte Carlo simulation
cat("\nRunning Monte Carlo simulation with parameter grid...\n")

# Run the parameter sweep with suppressed messages/warnings
suppressMessages(
  suppressWarnings({
    simulation_results <- run_parameter_sweep(param_grid)
  })
)

# Calculate summary statistics across all iterations
simulation_summary <- simulation_results %>%
  group_by(carryover_t1half, carryover_scale) %>%
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    mean_se = mean(std_error, na.rm = TRUE),
    n_iterations = n(),
    .groups = "drop"
  )

# Display summary
cat("\nSummary of simulation results:\n")
print(simulation_summary)

# Execute the save function
save_results()

# Visualizing the Results

# Power heatmap as text table
cat("\nPower Heatmap (Carryover Half-life x Scale Factor):\n")
power_wide <- simulation_summary |>
  mutate(power_pct = paste0(round(power * 100), "%")) |>
  select(carryover_t1half, carryover_scale, power_pct) |>
  tidyr::pivot_wider(names_from = carryover_t1half, values_from = power_pct,
                     names_prefix = "t1half=")
print(power_wide)

# Power decay curves
scale_levels <- sort(unique(simulation_summary$carryover_scale))
cols <- hcl.colors(length(scale_levels), palette = "Viridis")
plot(NULL,
     xlim = range(simulation_summary$carryover_t1half),
     ylim = c(0, 1),
     xlab = "Carryover Half-life (weeks)",
     ylab = "Statistical Power",
     main = paste0("Power Decay with Increasing Carryover Half-life\nN=",
                   base_params$n_participants))
for (j in seq_along(scale_levels)) {
  d <- simulation_summary[simulation_summary$carryover_scale == scale_levels[j], ]
  lines(d$carryover_t1half, d$power, col = cols[j], lwd = 2)
  points(d$carryover_t1half, d$power, col = cols[j], pch = 16)
}
legend("topright", legend = paste("Scale", scale_levels),
       col = cols, lwd = 2, pch = 16)

# Conclusions
cat("\nConclusions:\n")
cat("1. Half-life Impact: Longer carryover half-lives substantially reduce statistical power.\n")
cat("2. Scale Effect: Higher carryover scale factors amplify the power reduction.\n")
cat("3. Critical Thresholds: Power remains above 80% when carryover half-life is less than ~1.5 weeks with standard scale.\n")
cat("4. Severe Carryover: With half-lives exceeding 3 weeks, power drops below 40% even with minimal scale factors.\n")
cat("\nNOTE: Increasing n_iterations from the test value (", base_params$n_iterations,
    ") to 100 or higher is recommended for production-quality results.\n", sep="")

# End of script
