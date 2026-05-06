#' vig5.R: Self-contained Monte Carlo clinical trial design comparison
#'
#' This is a monolithic script that contains all the code needed to run the
#' Monte Carlo simulation from the monte_carlo_design_comparison vignette.
#' All helper functions and simulation code are included directly in this file.

# Load required packages
# Suppress package startup messages
# and clear the workspace
rm(list = ls())

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
#===========================================================================
# End Function: modify_list
#===========================================================================

#===========================================================================
# Function: cumulative
# Description: Wrapper around cumsum for calculating cumulative values
#===========================================================================
cumulative <- function(values) {
  cumsum(values)
}
#===========================================================================
# End Function: cumulative
#===========================================================================

#===========================================================================
# Function: mod_gompertz
# Description: Modified Gompertz function for modeling growth curves
#===========================================================================
mod_gompertz <- function(time, max_value, displacement, rate) {
  max_value * (1 - exp(-displacement * exp(-rate * time)))
}
#===========================================================================
# End Function: mod_gompertz
#===========================================================================

#===========================================================================
# Function: %||%
# Description: Null coalescing operator - returns first argument if not NULL,
#              otherwise returns second argument
#===========================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a
#===========================================================================
# End Function: %||%
#===========================================================================

#===========================================================================
# Function: analyze_trial
# Description: Analyze simulated trial data using linear mixed effects model
#===========================================================================
analyze_trial <- function(data, alpha = 0.05) {
  tryCatch({
    model <- nlme::lme(response ~ treatment * biomarker,
                       random = ~1 | participant_id,
                       data = data,
                       na.action = na.omit)

    model_summary <- summary(model)
    coefs <- model_summary$tTable

    interaction_idx <- which(rownames(coefs) == "treatment:biomarker")

    if (length(interaction_idx) > 0) {
      effect_size <- coefs[interaction_idx, "Value"]
      std_error <- coefs[interaction_idx, "Std.Error"]
      t_value <- coefs[interaction_idx, "t-value"]
      p_value <- coefs[interaction_idx, "p-value"]
      significant <- p_value < alpha

      results <- tibble(
        effect_size = effect_size,
        std_error = std_error,
        t_value = t_value,
        p_value = p_value,
        significant = significant
      )
    } else {
      results <- tibble(
        effect_size = NA_real_,
        std_error = NA_real_,
        t_value = NA_real_,
        p_value = NA_real_,
        significant = NA
      )
    }

    results
  }, error = function(e) {
    message("Model fitting error: ", e$message)
    tibble(
      effect_size = NA_real_,
      std_error = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      significant = NA
    )
  })
}
#===========================================================================
# End Function: analyze_trial
#===========================================================================

#===========================================================================
# Function: calculate_carryover
# Description: Calculate carryover effect using different decay models
#===========================================================================
calculate_carryover <- function(time_since_discontinuation, previous_effect, model = "exponential", params) {
  switch(model,
    "exponential" = {
      # Exponential decay: previous_effect * (1/2)^(scale_factor * time/halflife)
      scale_factor <- params$scale_factor %||% 2
      halflife <- params$carryover_t1half %||% 1
      previous_effect * (1/2)^(scale_factor * time_since_discontinuation / halflife)
    },

    "linear" = {
      # Linear decay to zero over total_time weeks
      total_time <- params$total_time %||% (2 * (params$carryover_t1half %||% 1))
      decay_factor <- pmax(0, 1 - time_since_discontinuation / total_time)
      previous_effect * decay_factor
    },

    "weibull" = {
      # Weibull decay: exp(-(time/lambda)^k)
      k <- params$k %||% 1  # k=1 gives exponential, k>1 gives accelerating decay
      halflife <- params$carryover_t1half %||% 1
      lambda <- params$lambda %||% (halflife / (log(2)^(1/k)))
      previous_effect * exp(-(time_since_discontinuation / lambda)^k)
    },

    stop("Unknown carryover model: ", model)
  )
}
#===========================================================================
# End Function: calculate_carryover
#===========================================================================

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

  # Additional verbose output
  if (verbose == TRUE) {
    cat("carryover: ")
    print(model_param$carryover_t1half)
    cat("br correlations: \n")
    br_indices <- grep("\\.br$", labels)
    if (length(br_indices) > 0) {
      print(correlations[br_indices, 1])
    } else {
      cat("No bio_response (br) columns found in correlations matrix\n")
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
#===========================================================================
# End Function: generate_data
#===========================================================================

#===========================================================================
# Function: lme_analysis
# Description: Linear mixed effects analysis of trial data
#===========================================================================
lme_analysis <- function(trial_design_set, data, options = list()) {
  # Set default options if not provided
  default_options <- list(
    use_expectancy = TRUE,
    random_slope = FALSE,
    full_output = FALSE,
    simple_carryover = FALSE,
    carryover_halflife = 0,
    carryover_scale_factor = 1
  )

  # Merge provided options with defaults
  options <- modify_list(default_options, options)

  # Validate incompatible options
  if ((options$carryover_halflife > 0) && (options$simple_carryover == TRUE)) {
    stop("Cannot use both simple_carryover and carryover_halflife models simultaneously")
  }

  # Initialize
  num_groups <- length(trial_design_set)

  # Pre-process all trial designs at once to avoid repetitive operations
  trial_designs_with_baseline <- map(1:num_groups, function(group_idx) {
    # Get this group's trial design
    trial_design <- trial_design_set[[group_idx]]

    # Add baseline row and process
    bind_rows(
      tibble(
        timepoint_names = "baseline",
        t_wk = 0,
        e = 0,
        tod = 0,
        tsd = 0,
        tpb = 0
      ),
      as_tibble(trial_design)
    ) %>%
    mutate(
      t = cumulative(t_wk),
      group = group_idx  # Add group identifier for joining later
    )
  })

  # For large datasets, use a map-reduce approach instead of growing data frames
  # Determine if it's a large dataset
  large_dataset <- nrow(data) > 10000

  if (large_dataset) {
    # Process groups in parallel using map and then combine results only once

    # First get unique participant IDs per group to calculate offsets
    group_sizes <- vector("numeric", num_groups)

    # Add path column if it doesn't exist
    if (!"path" %in% names(data)) {
      # If path column doesn't exist, assume all data belongs to group 1
      data$path <- 1
      message("Note: 'path' column not found in data. Adding path=1 to all rows.")
    }

    for (group_idx in 1:num_groups) {
      group_data <- data %>% filter(path == group_idx)
      group_sizes[group_idx] <- n_distinct(group_data$participant_id)
    }

    # Calculate participant_id offsets for each group to ensure uniqueness
    id_offsets <- c(0, cumsum(group_sizes)[-num_groups])

    # Process each group and then combine
    processed_groups <- map(1:num_groups, function(group_idx) {
      # Get this group's trial design and data
      trial_design <- trial_designs_with_baseline[[group_idx]]
      group_data <- data %>%
        filter(path == group_idx) %>%
        as_tibble()

      # Extract data for this group
      timepoint_names <- c("baseline", unique(trial_design$timepoint_names))

      # Check which timepoint names actually exist in the data to avoid missing values
      valid_timepoint_names <- timepoint_names[timepoint_names %in% names(group_data)]

      # Process group data efficiently
      group_processed <- group_data %>%
        select(participant_id, biomarker, all_of(valid_timepoint_names)) %>%
        # Adjust participant IDs to ensure uniqueness across groups
        mutate(participant_id = participant_id + id_offsets[group_idx])

      # Convert to long format
      group_long <- group_processed %>%
        pivot_longer(
          cols = all_of(valid_timepoint_names),
          names_to = "timepoint_names",
          values_to = "symptoms",
          values_drop_na = FALSE
        )

      # Merge with trial design information
      group_long %>%
        left_join(trial_design, by = "timepoint_names")
    })

    # Combine all group data efficiently
    all_data <- bind_rows(processed_groups)

  } else {
    # Original approach for smaller datasets
    all_data <- tibble()
    last_participant_id <- 0

    # Process each group in the trial design
    for (group_idx in 1:num_groups) {
      # Get this group's trial design and data
      trial_design <- trial_designs_with_baseline[[group_idx]]
      group_data <- data %>%
        filter(path == group_idx) %>%
        as_tibble()

      # Extract data for this group
      timepoint_names <- c("baseline", unique(trial_design$timepoint_names))

      # Check which timepoint names actually exist in the data to avoid missing values
      valid_timepoint_names <- timepoint_names[timepoint_names %in% names(group_data)]

      # Debug message if there are missing timepoints
      if (length(valid_timepoint_names) < length(timepoint_names)) {
        missing_timepoints <- setdiff(timepoint_names, valid_timepoint_names)
        message("Warning: Missing timepoint columns: ", paste(missing_timepoints, collapse=", "))
      }

      group_processed <- group_data %>%
        select(participant_id, biomarker, all_of(valid_timepoint_names)) %>%
        # Adjust participant IDs to ensure uniqueness across groups
        mutate(participant_id = participant_id + last_participant_id)

      # Update last participant ID for next group
      last_participant_id <- max(group_processed$participant_id)

      # Convert to long format
      group_long <- group_processed %>%
        pivot_longer(
          cols = all_of(valid_timepoint_names),
          names_to = "timepoint_names",
          values_to = "symptoms",
          values_drop_na = FALSE
        )

      # Merge with trial design information
      group_merged <- group_long %>%
        left_join(trial_design, by = "timepoint_names")

      # Add to combined data
      all_data <- bind_rows(all_data, group_merged)
    }
  }

  # Calculate derived variables for model
  data_for_model <- all_data %>%
    mutate(
      drug_binary = ifelse(tod > 0, 1, 0)
    )

  # Apply carryover effects if specified - optimized calculation
  if (options$carryover_halflife > 0) {
    # Pre-calculate the half-life factor for efficiency
    half_life_factor <- options$carryover_scale_factor / options$carryover_halflife

    # For large datasets, process each participant group separately
    if (large_dataset) {
      # Get unique participant IDs
      unique_ids <- unique(data_for_model$participant_id)
      n_participants <- length(unique_ids)

      # Create a lookup table for faster filtering
      participant_lookup <- vector("list", max(unique_ids))
      for (id in unique_ids) {
        participant_lookup[[id]] <- which(data_for_model$participant_id == id)
      }

      # Create a copy of drug_binary for updating
      data_for_model$drug_binary_new <- data_for_model$drug_binary

      # Process each participant separately
      for (id in unique_ids) {
        # Get indices for this participant
        idx <- participant_lookup[[id]]
        if (length(idx) < 2) next  # Skip if only one row

        # Get data for this participant
        participant_rows <- data_for_model[idx, ]

        # Order by time
        participant_rows <- participant_rows[order(participant_rows$t), ]

        # Process each timepoint
        for (i in 2:nrow(participant_rows)) {
          # Get current row
          current_row <- participant_rows[i, ]

          # Skip if on drug, NA, or invalid tsd
          if (is.na(current_row$tod) ||
              current_row$tod > 0 ||
              is.na(current_row$tsd) ||
              current_row$tsd <= 0) next

          # Calculate carryover
          prev_drug_binary <- participant_rows[i-1, "drug_binary"]
          decay <- (1/2)^(half_life_factor * current_row$tsd)
          new_drug_binary <- prev_drug_binary * decay

          # Update in the main dataset
          data_for_model$drug_binary_new[idx[i]] <- new_drug_binary
        }
      }

      # Replace drug_binary with new values
      data_for_model$drug_binary <- data_for_model$drug_binary_new
      data_for_model$drug_binary_new <- NULL
    } else {
      # For smaller datasets, use the tidyverse approach
      data_for_model <- data_for_model %>%
        group_by(participant_id) %>%
        mutate(
          # Calculate exponential falloff for carryover with optimized calculation
          # Handle NA values safely
          drug_binary = case_when(
            is.na(tod) ~ NA_real_,
            tod > 0 ~ 1,
            is.na(tsd) ~ 0,
            tsd <= 0 ~ 0,
            TRUE ~ dplyr::lag(drug_binary, default = 0) * (1/2)^(half_life_factor * tsd)
          )
        ) %>%
        ungroup()
    }
  } else if (options$simple_carryover) {
    # Simple carryover is already efficient
    data_for_model <- data_for_model %>%
      mutate(tsd_effect = tsd)
  }

  # Calculate group average at each timepoint after baseline - optimized
  # Avoid recalculation by computing means once per timepoint
  if (large_dataset) {
    # Calculate means by timepoint first
    timepoint_means <- data_for_model %>%
      group_by(timepoint_names) %>%
      summarize(mean_symptoms = mean(symptoms, na.rm = TRUE), .groups = "drop")

    # Join back to main data
    data_for_model <- data_for_model %>%
      left_join(timepoint_means, by = "timepoint_names")
  } else {
    # For smaller datasets, use the standard approach
    data_for_model <- data_for_model %>%
      group_by(timepoint_names) %>%
      mutate(mean_symptoms = mean(symptoms, na.rm = TRUE)) %>%
      ungroup()
  }

  # Ensure consistent column order
  data_for_model <- data_for_model %>%
    select(participant_id, biomarker, timepoint_names, t, symptoms, drug_binary, everything())

  # Build formula based on options
  formula_str <- "symptoms ~ biomarker"

  # Add drug effect (drug_binary)
  formula_str <- paste(formula_str, "+ drug_binary")

  # Add interaction term (biomarker * drug effect)
  formula_str <- paste(formula_str, "+ biomarker:drug_binary")

  # Add time effect if using it
  if (options$use_expectancy) {
    formula_str <- paste(formula_str, "+ t")
  }

  # Add simple carryover effect if using it
  if (options$simple_carryover) {
    formula_str <- paste(formula_str, "+ tsd_effect")
  }

  # Build random effects formula
  if (options$random_slope) {
    random_formula <- ~1 + drug_binary | participant_id
  } else {
    random_formula <- ~1 | participant_id
  }

  # Convert fixed effects to formula
  fixed_formula <- as.formula(formula_str)

  # Fit the model and capture any warnings
  fit_messages <- character()
  withCallingHandlers(
    model <- nlme::lme(fixed_formula, random = random_formula,
                       data = data_for_model, na.action = na.omit),
    warning = function(w) {
      fit_messages <<- c(fit_messages, w$message)
    }
  )

  # Extract relevant summary information
  model_summary <- summary(model)

  # Find the relevant coefficient (biomarker:drug_binary interaction)
  beta_row <- which(rownames(model_summary$tTable) == "biomarker:drug_binary")

  # Create output
  output <- tibble(
    beta = model_summary$tTable[beta_row, "Value"],
    beta_se = model_summary$tTable[beta_row, "Std.Error"],
    p_value = model_summary$tTable[beta_row, "p-value"],
    # Singular fit check not applicable to nlme
    is_singular = FALSE,
    warning = paste(fit_messages, collapse = "; ")
  )

  # Return appropriate output based on options
  if (options$full_output) {
    return(list(
      formula = formula,
      model = model,
      data = data_for_model,
      summary = output
    ))
  } else {
    return(output)
  }
}
#===========================================================================
# End Function: lme_analysis
#===========================================================================

#===========================================================================
# Function: run_monte_carlo
# Description: Run a Monte Carlo simulation for a design and parameters
#===========================================================================
run_monte_carlo <- function(design_name, params) {
  # Update model parameters with current simulation parameters
  current_model_params <- model_params
  current_model_params$N <- params$n_participants
  current_model_params$c.bm <- params$biomarker_correlation
  current_model_params$carryover_t1half <- params$carryover_t1half

  # Select the appropriate design template
  if (design_name == "hybrid") {
    design_template <- designs$hybrid$design
  } else {
    design_template <- designs$crossover$design
  }

  # Set number of iterations
  n_iter <- params$n_iterations

  # Initialize results storage
  all_results <- list()

  # Run Monte Carlo iterations
  cat(paste0("Running ", n_iter, " iterations for ", design_name, " design...\n"))

  all_results <- purrr::map_dfr(1:n_iter, function(iter) {
    if (iter %% 10 == 0) {
      cat(paste0("  Iteration ", iter, " of ", n_iter, "\n"))
    }

    # Generate data using generate_data
    sim_data <- generate_data(
      model_param = current_model_params,
      resp_param = resp_param,
      baseline_param = baseline_param,
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
        treatment = if (design_name == "hybrid") {
          ifelse(week <= 8, 1,
                 ifelse(week <= 12, 0,
                        ifelse(week <= 16, 1, 0)))
        } else {
          # For crossover design
          ifelse(participant_id %% 2 == 1,
                 ifelse(week <= 10, 1, 0),
                 ifelse(week <= 10, 0, 1))
        }
      )

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
        design = design_name,
        n_participants = params$n_participants,
        biomarker_correlation = params$biomarker_correlation,
        carryover_t1half = params$carryover_t1half
      )
  })

  return(all_results)
}
#===========================================================================
# End Function: run_monte_carlo
#===========================================================================

#===========================================================================
# Function: run_parameter_sweep
# Description: Run simulations for multiple parameter combinations
#===========================================================================
run_parameter_sweep <- function(param_grid) {
  # Number of parameter combinations
  n_combinations <- nrow(param_grid)

  # Run simulations for each parameter combination
  results <- map_dfr(1:n_combinations, function(i) {
    # Extract current parameters
    current_params <- as.list(param_grid[i,])
    current_params$n_iterations <- base_params$n_iterations
    current_params$treatment_effect <- base_params$treatment_effect
    current_params$carryover_scale <- base_params$carryover_scale
    current_params$between_subject_sd <- base_params$between_subject_sd
    current_params$within_subject_sd <- base_params$within_subject_sd

    # Run Monte Carlo for each design
    cat("Running parameter combination", i, "of", n_combinations, ":",
        "n =", current_params$n_participants,
        "corr =", current_params$biomarker_correlation,
        "carryover =", current_params$carryover_t1half, "\n")

    # Run Hybrid N-of-1 design
    cat("  - Running Hybrid N-of-1 design...\n")
    hybrid_results <- run_monte_carlo("hybrid", current_params)

    # Run Traditional Crossover design
    cat("  - Running Traditional Crossover design...\n")
    crossover_results <- run_monte_carlo("crossover", current_params)

    # Combine results
    bind_rows(hybrid_results, crossover_results)
  })

  return(results)
}
#===========================================================================
# End Function: run_parameter_sweep
#===========================================================================

#===========================================================================
# Function: create_designs
# Description: Create trial designs for a given number of participants
#===========================================================================
create_designs <- function(n_participants) {
  # Create Hybrid N-of-1 design (Design 4)
  # Open-label (8wk) + blinded discontinuation (4wk) + crossover (8wk: 4wk active, 4wk placebo)
  hybrid_design <- expand_grid(
    participant_id = 1:n_participants,
    week = 1:20
  ) %>%
  mutate(
    treatment = c(rep(1, 8), rep(0, 4), rep(1, 4), rep(0, 4))[week],
    period = c(rep(1, 8), rep(2, 4), rep(3, 8))[week]
  )

  # Create Traditional Crossover design (Design 3)
  # Half start with treatment then placebo, half with placebo then treatment
  crossover_design <- expand_grid(
    participant_id = 1:n_participants,
    week = 1:20
  ) %>%
  mutate(
    # Determine sequence based on participant ID (odd=AB, even=BA)
    sequence = if_else(participant_id %% 2 == 0, "AB", "BA"),
    # Set treatment based on sequence
    treatment = case_when(
      sequence == "AB" & week <= 10 ~ 1,
      sequence == "AB" & week > 10 ~ 0,
      sequence == "BA" & week <= 10 ~ 0,
      sequence == "BA" & week > 10 ~ 1,
      TRUE ~ NA_real_
    ),
    period = if_else(week <= 10, 1, 2)
  )

  # Return both designs in a list
  return(list(
    hybrid = hybrid_design,
    crossover = crossover_design
  ))
}
#===========================================================================
# End Function: create_designs
#===========================================================================

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
#===========================================================================
# End Function: prepare_design_for_simulation
#===========================================================================

# ==================
# Main script begins
# ==================

# Base set of parameters - these will be varied in the Monte Carlo simulation
# Parameters informed by Hendrickson et al. (2020) for realistic clinical scenarios
base_params <- list(
  # Sample size - modest for pragmatic trial designs
  n_participants = 30,

  # Biomarker-response correlation (default - will be varied)
  # Hendrickson et al. studied correlations in this range
  biomarker_correlation = 0.4,  # Moderate correlation as default (reduced from 0.6 for more realistic power)

  # Treatment effect (baseline)
  treatment_effect = 3.5,  # Reduced from 5 for more realistic power

  # Carryover effect parameters
  carryover_t1half = 1,     # Half-life of carryover effect in weeks
  carryover_scale = 1.5,    # Stronger scale factor (increases effect)

  # Random variation - increased to make signal detection more challenging
  between_subject_sd = 3.0,  # Increased between-subject variability
  within_subject_sd = 2.8,   # Increased within-subject variability

  # Number of Monte Carlo iterations
  n_iterations = 100  # Default number for full runs
)

# Define parameter ranges for sensitivity analysis
# Reduced from 144 combinations to 36 combinations for faster execution
param_grid <- expand_grid(
  n_participants = c(20, 70),  # Reduced from 3 to 2 sample sizes (smallest and largest)
  biomarker_correlation = c(0.2, 0.5, 0.8),  # Reduced from 4 to 3 correlation values
  carryover_t1half = c(0, 1.5, 2.5),  # Reduced from 4 to 3 carryover values
  treatment_effect = c(3.0, 5.0)  # Reduced from 3 to 2 treatment effect values
)

param_grid <- expand_grid(
  n_participants = c(70),  # Reduced from 3 to 2 sample sizes (smallest and largest)
  biomarker_correlation = c(0.2,  0.8),  # Reduced from 4 to 3 correlation values
  carryover_t1half = c(0, 1.5),  # Reduced from 4 to 3 carryover values
  treatment_effect = c(3.0)  # Reduced from 3 to 2 treatment effect values
)
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
resp_param <- tibble(
  cat = c("time_variant", "pharm_biomarker", "bio_response"),
  max = c(1.0, 1.0, base_params$treatment_effect),  # Maximum response value
  disp = c(2.0, 2.0, 2.0),  # Displacement
  rate = c(0.3, 0.3, 0.3),  # Rate
  sd = c(base_params$within_subject_sd, base_params$within_subject_sd, base_params$within_subject_sd)  # Standard deviation
)

# Set up baseline parameters
baseline_param <- tibble(
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

# Create designs for base sample size for visualization
original_designs <- create_designs(base_params$n_participants)

# Create designs compatible with the generate_data function
hybrid_design <- prepare_design_for_simulation(original_designs$hybrid %>%
                                            filter(participant_id == 1))
crossover_design <- prepare_design_for_simulation(original_designs$crossover %>%
                                               filter(participant_id == 1))

# Store designs in a list
designs <- list(
  hybrid = list(design = hybrid_design, name = "Hybrid N-of-1"),
  crossover = list(design = crossover_design, name = "Traditional Crossover")
)

# ==========================================
# Main simulation - use a reduced sample for quicker execution
# ==========================================

# Define a reduced demo grid for faster computation
demo_grid <- expand_grid(
  n_participants = c(20, 50, 70),  # Three sample sizes
  biomarker_correlation = c(0.2, 0.4, 0.6, 0.8),  # Four correlation values
  carryover_t1half = c(0, 0.5, 1.5, 2.5),  # Now includes zero carryover
  treatment_effect = c(3.0, 4.0, 5.0)  # Three treatment effect values
)
demo_grid <- param_grid
# Use a reduced number of iterations for faster execution
demo_params <- list(
  n_iterations = 10,  # Increased from 20
  treatment_effect = base_params$treatment_effect,
  carryover_scale = base_params$carryover_scale,
  between_subject_sd = base_params$between_subject_sd,
  within_subject_sd = base_params$within_subject_sd
)

# Use the full parameter grid instead of a subset
demo_subset <- demo_grid

# Initialize results storage
demo_results <- tibble()
iteration_tracker <- tibble()

# For each parameter combination in our subset
for (i in 1:nrow(demo_subset)) {
  # Extract current parameters
  current_params <- as.list(demo_subset[i,])
  current_params$n_iterations <- demo_params$n_iterations
  current_params$treatment_effect <- demo_params$treatment_effect
  current_params$carryover_scale <- demo_params$carryover_scale

  # Print current parameters
  cat("Running parameter combination", i, "of", nrow(demo_subset), ":",
      "n =", current_params$n_participants,
      "corr =", current_params$biomarker_correlation,
      "carryover =", current_params$carryover_t1half, "\n")

  # Run for each design using our Monte Carlo function
  for (design_name in c("hybrid", "crossover")) {
    # Run the Monte Carlo simulation
    cat("  - Running", design_name, "design...\n")
    design_results <- run_monte_carlo(design_name, current_params)

    # Add to overall results
    demo_results <- bind_rows(demo_results, design_results)

    # Calculate cumulative power at different iteration counts
    for (iter_cutoff in c(5, 10, 15)) {
      if (iter_cutoff <= demo_params$n_iterations) {
        # Calculate power using results up to this iteration
        cumulative_power <- design_results %>%
          filter(iteration <= iter_cutoff) %>%
          summarize(power = mean(significant, na.rm = TRUE)) %>%
          pull(power)

        # Add to iteration tracker
        iteration_tracker <- bind_rows(
          iteration_tracker,
          tibble(
            design = design_name,
            n_participants = current_params$n_participants,
            biomarker_correlation = current_params$biomarker_correlation,
            carryover_t1half = current_params$carryover_t1half,
            iterations = iter_cutoff,
            power = cumulative_power
          )
        )
      }
    }
  }
}

# Calculate summary statistics across all iterations
simulation_summary <- demo_results %>%
  group_by(design, n_participants, biomarker_correlation, carryover_t1half) %>%
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    mean_se = mean(std_error, na.rm = TRUE),
    n_iterations = n(),
    .groups = "drop"
  )

# Display summary statistics
cat("Summary of power values:", "\n")
print(summary(simulation_summary$power))
cat("Range of power values:", range(simulation_summary$power, na.rm = TRUE), "\n")
cat("Count of 'significant = TRUE' results:", sum(demo_results$significant, na.rm = TRUE), "\n")
cat("Total number of simulation runs:", nrow(demo_results), "\n")

# Display iteration stats
iteration_stats <- iteration_tracker %>%
  arrange(design, n_participants, biomarker_correlation, carryover_t1half, iterations)

# Save the simulation results for visualization and further analysis
sim_results <- list(
  demo_results = demo_results,
  iteration_tracker = iteration_tracker,
  simulation_summary = simulation_summary
)

# ==========================================
# Visualizations
# ==========================================

# Power comparison by design and parameters
unique_designs <- unique(simulation_summary$design)
design_cols <- hcl.colors(length(unique_designs), palette = "Set 2")

plot(NULL, xlim = range(simulation_summary$biomarker_correlation),
     ylim = c(0, 1),
     xlab = "Biomarker-Treatment Correlation",
     ylab = "Power",
     main = "Statistical Power by Design, Sample Size, and Carryover Effect")
legend_labels <- character()
legend_cols <- character()
legend_lty <- integer()
idx <- 0
for (j in seq_along(unique_designs)) {
  d_sub <- simulation_summary[simulation_summary$design == unique_designs[j], ]
  combos <- unique(d_sub[, c("n_participants", "carryover_t1half")])
  for (k in seq_len(nrow(combos))) {
    idx <- idx + 1
    dd <- d_sub[d_sub$n_participants == combos$n_participants[k] &
                d_sub$carryover_t1half == combos$carryover_t1half[k], ]
    lines(dd$biomarker_correlation, dd$power,
          col = design_cols[j], lwd = 2, lty = k)
    points(dd$biomarker_correlation, dd$power,
           col = design_cols[j], pch = 16)
    legend_labels <- c(legend_labels,
                       paste0(unique_designs[j], " n=", combos$n_participants[k],
                              " t1/2=", combos$carryover_t1half[k]))
    legend_cols <- c(legend_cols, design_cols[j])
    legend_lty <- c(legend_lty, k)
  }
}
if (length(legend_labels) <= 10) {
  legend("bottomright", legend = legend_labels, col = legend_cols,
         lwd = 2, lty = legend_lty, cex = 0.6)
}

# Distribution of effect estimates (histogram per design)
par(mfrow = c(1, length(unique_designs)))
for (des in unique_designs) {
  d <- demo_results$effect_size[demo_results$design == des]
  d <- d[!is.na(d)]
  if (length(d) > 0) {
    hist(d, breaks = 30, col = adjustcolor("steelblue", 0.5),
         main = paste("Effect Size:", des),
         xlab = "Estimated Interaction Effect")
  }
}
par(mfrow = c(1, 1))

# Carryover decay visualization
weeks <- 0:10
treatment_effect <- base_params$treatment_effect

carryover_scenarios <- expand_grid(
  model = c("exponential", "linear"),
  carryover_t1half = c(0, 0.5, 1.5, 2.5)
)

decay_data <- carryover_scenarios |>
  group_by(model, carryover_t1half) |>
  do({
    model_type <- .$model[1]
    halflife <- .$carryover_t1half[1]

    params <- list(
      carryover_t1half = halflife,
      scale_factor = 1,
      total_time = halflife * 4,
      k = 1.5
    )

    tibble(
      week = weeks,
      carryover_effect = sapply(weeks, function(w) {
        if (w == 0) return(treatment_effect)
        calculate_carryover(w, treatment_effect, model_type, params)
      })
    )
  }) |>
  ungroup()

unique_combos <- unique(decay_data[, c("model", "carryover_t1half")])
combo_cols <- hcl.colors(nrow(unique_combos), palette = "Set 2")

plot(NULL, xlim = range(weeks), ylim = c(0, treatment_effect * 1.1),
     xlab = "Weeks Since Discontinuation",
     ylab = "Carryover Effect Magnitude",
     main = "Carryover Effect Decay Patterns")
for (j in seq_len(nrow(unique_combos))) {
  dd <- decay_data[decay_data$model == unique_combos$model[j] &
                   decay_data$carryover_t1half == unique_combos$carryover_t1half[j], ]
  lty_val <- if (unique_combos$model[j] == "exponential") 1 else 2
  lines(dd$week, dd$carryover_effect, col = combo_cols[j], lwd = 2, lty = lty_val)
  points(dd$week, dd$carryover_effect, col = combo_cols[j], pch = 16, cex = 0.7)
}
legend("topright",
       legend = paste0(unique_combos$model, " t1/2=", unique_combos$carryover_t1half),
       col = combo_cols, lwd = 2,
       lty = ifelse(unique_combos$model == "exponential", 1, 2),
       cex = 0.6)

# Display summary statistics
cat("\nFinal Monte Carlo Simulation Results:\n")
print(simulation_summary)
# dim
dim(simulation_summary)
# ==========================================
# End of script
# ==========================================
