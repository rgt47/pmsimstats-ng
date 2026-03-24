#' pm_functions.R: Common functions for precision medicine trial simulation
#'
#' This script contains shared functions extracted from vig5.R for use in
#' multiple vignettes and analysis scripts.

#===========================================================================
# Function: mod_gompertz
# Origin-passing Gompertz (matches orig modgompertz exactly)
#===========================================================================
mod_gompertz <- function(t, maxr, disp, rate) {
  y <- maxr * exp(-disp * exp(-rate * t))
  vert_offset <- maxr * exp(-disp * exp(-rate * 0))
  y <- y - vert_offset
  y * (maxr / (maxr - vert_offset))
}



#===========================================================================
# Function: calculate_carryover
# Description: Calculate carryover effect using different decay
# models
#===========================================================================
calculate_carryover <- function(time_since_discontinuation,
                                previous_effect,
                                model = "exponential", params) {
  switch(model,
    "exponential" = {
      # Exponential decay: previous_effect *
      # (1/2)^(time/halflife)
      halflife <- if(is.null(params$carryover_t1half)) {
        1
      } else {
        params$carryover_t1half
      }
      previous_effect * (1/2)^(time_since_discontinuation /
                               halflife)
    },

    "linear" = {
      # Linear decay to zero over total_time weeks
      halflife_for_total <- if(is.null(params$carryover_t1half)) {
        1
      } else {
        params$carryover_t1half
      }
      total_time <- if(is.null(params$total_time)) {
        (2 * halflife_for_total)
      } else {
        params$total_time
      }
      decay_factor <- pmax(0, 1 - time_since_discontinuation /
                                   total_time)
      previous_effect * decay_factor
    },

    "weibull" = {
      # Weibull decay: exp(-(time/lambda)^k)
      # k=1 gives exponential, k>1 gives accelerating decay
      k <- if(is.null(params$k)) 1 else params$k
      halflife <- if(is.null(params$carryover_t1half)) {
        1
      } else {
        params$carryover_t1half
      }
      lambda <- if(is.null(params$lambda)) {
        (halflife / (log(2)^(1/k)))
      } else {
        params$lambda
      }
      previous_effect * exp(-(time_since_discontinuation /
                              lambda)^k)
    },

    stop("Unknown carryover model: ", model)
  )
}
#===========================================================================
# End Function: calculate_carryover
#===========================================================================

#===========================================================================
# Enhanced Carryover Functions
#===========================================================================

#' Calculate bio-response with biomarker×treatment interaction
#'
#' @param trial_data Trial design data with tod column
#' @param model_param Model parameters including c.bm interaction
#'   strength
#' @param resp_param Response parameters for bio_response category
#' @param component_halflives Component-specific half-lives for
#'   carryover
#' @return List with bio_response_means and bio_response_test
calculate_bio_response_with_interaction <- function(
    trial_data, model_param, resp_param,
    component_halflives) {
  # Base response follows mod_gompertz for treatment periods
  base_bio_response <- mod_gompertz(
    trial_data$tod,
    resp_param$max[resp_param$cat == "br"],
    resp_param$disp[resp_param$cat == "br"],
    resp_param$rate[resp_param$cat == "br"]
  )

  # NOTE: Biomarker×treatment interaction emerges from correlation structure
  # in the MVN draw, following Hendrickson et al. (2020) approach.
  # No population mean shift needed - interaction is created by differential
  # correlation between biomarker and bio_response by treatment status.
  # See: WHY_POST_MVN_ADJUSTMENT_ANALYSIS.md
  bio_response_means <- base_bio_response

  # Check for zero values (for correlation matrix logic)
  bio_response_test <- base_bio_response == 0

  # Apply enhanced carryover to bio-response component
  bio_response_means <- apply_carryover_to_component(
    bio_response_means, trial_data, component_halflives, "br"
  )

  return(list(
    bio_response_means = bio_response_means,
    bio_response_test = bio_response_test,
    raw_bio_response_means = base_bio_response
  ))
}

#' Calculate component-specific carryover half-lives
#'
#' @param base_halflife Base carryover half-life from model
#'   parameters
#' @return List of component-specific half-lives
#'
#' @details Following Hendrickson et al. (2020), all components use
#' the same half-life. Only BR component actually receives carryover
#' adjustment (tv and pb are returned unchanged by
#' apply_carryover_to_component).
calculate_component_halflives <- function(base_halflife) {
  list(
    # Uniform half-life for all components (Hendrickson approach)
    br_halflife = base_halflife,
    pb_halflife = base_halflife,
    tv_halflife = base_halflife
  )
}

#' DEPRECATED: Calculate carryover-adjusted correlation parameters
#'
#' This function has been DEPRECATED and is no longer used.
#'
#' REASON FOR DEPRECATION:
#' Based on theoretical analysis (see carryover_correlation_theory.tex):
#' - Carryover affects MEANS (systematic component), not CORRELATIONS
#' - Dynamic adjustment violates positive definite constraints
#' CURRENT APPROACH:
#' Use fixed correlation values (Hendrickson et al. 2020) that do not
#' vary with carryover parameters. Carryover is modeled in the mean
#' structure only.

#' Carryover calculation for bio-response component (Hendrickson method)
#'
#' @param component_means Current component means
#' @param trial_data Trial design data
#' @param component_halflives Component-specific half-lives
#' @param component_name Component abbreviation ("tv", "pb", "br")
#' @return Adjusted means with carryover effects
#'
#' @details Following Hendrickson et al. (2020), carryover is applied
#' ONLY to the bio-response (br) component when participants are off
#' drug. The formula is: mu[t] = base[t] + mu[t-1] * (1/2)^(tsd/t1half)
#' where tsd is time since discontinuation.
apply_carryover_to_component <- function(
    component_means, trial_data, component_halflives,
    component_name) {

  num_timepoints <- length(component_means)

  # ONLY apply carryover to BR component (Hendrickson approach)
  if (component_name != "br") {
    return(component_means)  # No carryover for tv or pb
  }

  if (num_timepoints > 1) {
    # Get component-specific half-life
    halflife_key <- paste0(component_name, "_halflife")
    component_halflife <- component_halflives[[halflife_key]]

    if (is.null(component_halflife) || component_halflife == 0) {
      return(component_means)  # No carryover if halflife is 0
    }

    # Bio-response: carryover when off drug and time since
    # discontinuation > 0
    carryover_indices <- which(!trial_data$on_drug &
                                trial_data$tsd > 0)

    if (length(carryover_indices) > 0) {
      for (idx in carryover_indices) {
        prev_idx <- idx - 1

        # Safety check: ensure prev_idx is valid
        if (prev_idx >= 1 &&
            prev_idx <= length(component_means) &&
            idx >= 1 && idx <= length(component_means)) {

          # Use time since discontinuation for bio-response
          time_lag <- trial_data$tsd[idx]

          decay_factor <- (1/2)^(time_lag / component_halflife)

          # Apply carryover effect
          component_means[idx] <- component_means[idx] +
            component_means[prev_idx] * decay_factor
        }
      }
    }
  }

  return(component_means)
}

#===========================================================================
# Helper Functions for Data Generation
#===========================================================================

prepare_trial_data <- function(trial_design) {
  # Ensure trial_design is a data frame
  if (!is.data.frame(trial_design)) {
    stop("Trial design is not a data frame. Class: ", class(trial_design),
         ", Type: ", typeof(trial_design))
  }

  trial_data <- as_tibble(trial_design)

  # Check if t_wk exists and calculate cumulative if it does
  if ("t_wk" %in% names(trial_data)) {
    trial_data$t_wk_cumulative <- cumsum(trial_data$t_wk)
  } else {
    # If t_wk doesn't exist, calculate it or use week column directly
    if ("week" %in% names(trial_data)) {
      trial_data$t_wk <- c(trial_data$week[1], diff(trial_data$week))
      trial_data$t_wk_cumulative <- trial_data$week
    } else {
      stop("Neither 't_wk' nor 'week' column found in trial design data")
    }
  }

  trial_data <- trial_data %>% mutate(on_drug = (tod > 0))
  return(trial_data)
}

build_correlation_matrix <- function(
    labels, trial_design, model_param, num_timepoints,
    factor_types, factor_abbreviations,
    bio_response_test = NULL, bio_response_means = NULL,
    means = NULL, trial_data = NULL, lambda_cor = 0) {
  # Build a correlation matrix
  correlations <- diag(length(labels))
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  # Extract cumulative week values for AR(1) time gaps
  if (!is.null(trial_data) && "t_wk_cumulative" %in% names(trial_data)) {
    weeks <- trial_data$t_wk_cumulative
  } else if ("week" %in% names(trial_design)) {
    weeks <- trial_design$week
  } else if ("t_wk" %in% names(trial_design)) {
    weeks <- cumsum(trial_design$t_wk)
  } else {
    weeks <- seq_len(num_timepoints)
  }

  # Apply correlations between factors
  for (factor_idx in seq_along(factor_types)) {
    current_factor <- factor_abbreviations[factor_idx]

    # Build autocorrelations across time using AR(1): rho^|t_i - t_j|
    if (num_timepoints > 1) {
      autocorrelation <- model_param[[paste("c", current_factor, sep = ".")]]

      # Create all combinations of time points at once
      point_indices <- expand.grid(
        p1 = 1:(num_timepoints-1),
        p2 = (2:num_timepoints)
      )
      # Filter valid combinations where p2 > p1
      point_indices <- point_indices[
        point_indices$p2 > point_indices$p1,
      ]

      # Create name vectors for efficient indexing
      name1 <- paste(trial_design$timepoint_name[
                       point_indices$p1],
                     current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name[
                       point_indices$p2],
                     current_factor, sep = ".")

      # AR(1) decay: rho^|t_i - t_j|
      for (idx in 1:nrow(point_indices)) {
        time_gap <- abs(weeks[point_indices$p2[idx]] -
                        weeks[point_indices$p1[idx]])
        correlations[name1[idx], name2[idx]] <- autocorrelation^time_gap
        correlations[name2[idx], name1[idx]] <- autocorrelation^time_gap
      }
    }

    # Build autocorrelations across factors - VECTORIZED
    for (factor2_idx in setdiff(seq_along(factor_types), factor_idx)) {
      other_factor <- factor_abbreviations[factor2_idx]

      # Same timepoint cross-factor correlations
      name1 <- paste(trial_design$timepoint_name,
                     current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name,
                     other_factor, sep = ".")

      # Set all same-timepoint correlations at once
      for (idx in 1:length(name1)) {
        correlations[name1[idx], name2[idx]] <- model_param$c.cf1t
        correlations[name2[idx], name1[idx]] <- model_param$c.cf1t
      }

      # Different timepoint cross-factor correlations
      if (num_timepoints > 1) {
        # Create all combinations of time points at once
        point_indices <- expand.grid(
          p1 = 1:(num_timepoints-1),
          p2 = (2:num_timepoints)
        )
        # Filter valid combinations where p2 > p1
        point_indices <- point_indices[
          point_indices$p2 > point_indices$p1,
        ]

        # Create name vectors for efficient indexing
        name1 <- paste(trial_design$timepoint_name[
                         point_indices$p1],
                       current_factor, sep = ".")
        name2 <- paste(trial_design$timepoint_name[
                         point_indices$p2],
                       other_factor, sep = ".")

        # AR(1) decay for cross-factor different-time correlations
        for (idx in 1:nrow(point_indices)) {
          time_gap <- abs(weeks[point_indices$p2[idx]] -
                          weeks[point_indices$p1[idx]])
          correlations[name1[idx], name2[idx]] <-
            model_param$c.cfct * autocorrelation^time_gap
          correlations[name2[idx], name1[idx]] <-
            model_param$c.cfct * autocorrelation^time_gap
        }
      }
    }

    # Special handling for biomarker correlation
    if (current_factor == "br") {
      for (timepoint_idx in 1:num_timepoints) {
        name1 <- paste(trial_design$timepoint_name[
                         timepoint_idx], "br", sep = ".")

        on_drug_now <- if (!is.null(trial_data)) {
          trial_data$on_drug[timepoint_idx]
        } else {
          !is.null(bio_response_test) &&
            !bio_response_test[timepoint_idx]
        }

        if (on_drug_now) {
          correlations["bm", name1] <- model_param$c.bm
          correlations[name1, "bm"] <- model_param$c.bm
        } else {
          tsd_now <- if (!is.null(trial_data)) {
            trial_data$tsd[timepoint_idx]
          } else {
            0
          }
          if (tsd_now > 0 && lambda_cor > 0) {
            decay <- exp(-lambda_cor * tsd_now)
            correlations["bm", name1] <-
              model_param$c.bm * decay
            correlations[name1, "bm"] <-
              model_param$c.bm * decay
          }
        }
      }
    }
  }

  return(correlations)
}

#===========================================================================
# Function: generate_data
# Description: Primary data generation function for trial
# simulation
#===========================================================================
generate_data <- function(
    model_param, resp_param, baseline_param, trial_design,
    empirical, make_positive_definite, seed = NA,
    lambda_cor = NA, verbose = FALSE, track_pd_stats = TRUE,
    cached_sigma = NULL) {
  # Auto-compute lambda_cor from carryover half-life when not provided
  if (is.na(lambda_cor)) {
    if (!is.null(model_param$carryover_t1half) &&
        model_param$carryover_t1half > 0) {
      lambda_cor <- log(2) / model_param$carryover_t1half
    } else {
      lambda_cor <- 0
    }
  }

  # Initialize variables for tracking sigma matrix statistics
  sigma_count <- 0
  non_positive_definite_count <- 0

  # Use cached sigma if provided, otherwise build from scratch
  if (!is.null(cached_sigma)) {
    # Use pre-computed sigma matrix
    sigma <- cached_sigma$sigma
    labels <- cached_sigma$labels
    standard_deviations <- cached_sigma$standard_deviations
    correlations <- cached_sigma$correlations

    # Extract needed variables for data generation
    trial_data <- prepare_trial_data(trial_design)
    num_timepoints <- dim(trial_design)[1]
    factor_types <- c("tv", "pb", "br")
    factor_abbreviations <- c("tv", "pb", "br")

    # Skip sigma construction, means will be calculated below

  } else {
    # Original sigma construction path
    # I. Turn the trial design information into something easier to use
    trial_data <- prepare_trial_data(trial_design)
    num_timepoints <- dim(trial_design)[1]

    # II. Set up variables to track - baseline parameters for
    # each participant ("bm","BL"), and the three
    # modeled factors for each stage of the trial.

    # Set up the variable names
    factor_types <- c("tv", "pb", "br")
    factor_abbreviations <- c("tv", "pb", "br")

    labels <- c(
      c("bm", "BL"),
      paste(trial_design$timepoint_name, factor_abbreviations[1], sep = "."),
      paste(trial_design$timepoint_name, factor_abbreviations[2], sep = "."),
      paste(trial_design$timepoint_name, factor_abbreviations[3], sep = ".")
    )

    # Set up vectors with the standard deviations and means
    standard_deviations <- c(
      baseline_param$sd[baseline_param$cat == "bm"],
      baseline_param$sd[baseline_param$cat == "BL"]
    )
    standard_deviations <- c(
      standard_deviations,
      rep(resp_param$sd[resp_param$cat == "tv"], num_timepoints)
    )
    standard_deviations <- c(
      standard_deviations,
      rep(resp_param$sd[resp_param$cat == "pb"], num_timepoints) * trial_design$e
    )
    standard_deviations <- c(
      standard_deviations,
      rep(resp_param$sd[resp_param$cat == "br"], num_timepoints)
    )

    # Note: bio_response variables will be calculated below, so use NULL for now
    # The correlation matrix will be rebuilt after means calculation
    correlations <- NULL

    # Sigma matrix will be built after means calculation

  } # End of else block for sigma construction

  # Calculate component-specific half-lives for enhanced carryover
  component_halflives <- calculate_component_halflives(model_param$carryover_t1half)

  # Calculate means with enhanced carryover (works for both cached and non-cached)
  means <- c(
    baseline_param$m[baseline_param$cat == "bm"],
    baseline_param$m[baseline_param$cat == "BL"]
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

      # Apply enhanced carryover to time-variant component
      time_variant_means <- apply_carryover_to_component(
        time_variant_means, trial_data, component_halflives, "tv"
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

      # Apply enhanced carryover to pharmacological-biomarker component
      pharm_biomarker_means <- apply_carryover_to_component(
        pharm_biomarker_means, trial_data, component_halflives, "pb"
      )

      means <- c(means, pharm_biomarker_means)
    }

    if (factor_abbrev == "br") {
      # Calculate bio-response with biomarker×treatment interaction using helper function
      br_result <- calculate_bio_response_with_interaction(
        trial_data, model_param, resp_param, component_halflives
      )
      
      bio_response_means <- br_result$bio_response_means
      bio_response_test <- br_result$bio_response_test

      # Find indices for bio_response (br) columns dynamically
      br_indices <- grep("\\.br$", labels)

      # Set names for debugging - using dynamic indices
      if (length(br_indices) > 0) {
        names(bio_response_test) <- labels[br_indices]
      }

      means <- c(means, bio_response_means)
    }
  }

  # Build correlation matrix for non-cached path (if needed)
  if (is.null(cached_sigma)) {
    correlations <- build_correlation_matrix(
      labels, trial_design, model_param, num_timepoints,
      factor_types, factor_abbreviations,
      bio_response_test, bio_response_means, means,
      trial_data = trial_data, lambda_cor = lambda_cor
    )

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
  }

  # Pre-compute Cholesky factor (matches orig's draw method)
  chol_sigma <- tryCatch(chol(sigma), error = function(e) {
    chol(corpcor::make.positive.definite(sigma, tol = 1e-3))
  })

  # Set the seed if provided for reproducibility
  if (!is.na(seed)) {
    set.seed(seed)
  }

  # Draw participants using Cholesky factor (matches orig exactly)
  n <- model_param$N
  p <- length(means)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  participant_data <- Z %*% chol_sigma +
    matrix(means, nrow = n, ncol = p, byrow = TRUE)
  participant_data <- as_tibble(participant_data)
  colnames(participant_data) <- labels

  # Process participant data using helper function
  participant_data <- process_participant_data(
    participant_data, trial_design$timepoint_name,
    factor_abbreviations, model_param$N
  )

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
# Helper function for building and caching sigma matrices
#===========================================================================
build_sigma_matrix <- function(model_param, resp_param, baseline_param, trial_design,
                              factor_types, factor_abbreviations, verbose = FALSE,
                              lambda_cor = NA) {
  # This function builds just the sigma matrix without generating data

  trial_data <- prepare_trial_data(trial_design)
  num_timepoints <- dim(trial_design)[1]

  # Auto-compute lambda_cor from carryover half-life
  if (is.na(lambda_cor)) {
    if (!is.null(model_param$carryover_t1half) &&
        model_param$carryover_t1half > 0) {
      lambda_cor <- log(2) / model_param$carryover_t1half
    } else {
      lambda_cor <- 0
    }
  }

  # Set up variable names (same as generate_data)
  labels <- c(
    c("bm", "BL"),
    paste(trial_design$timepoint_name, factor_abbreviations[1], sep = "."),
    paste(trial_design$timepoint_name, factor_abbreviations[2], sep = "."),
    paste(trial_design$timepoint_name, factor_abbreviations[3], sep = ".")
  )

  # Set up standard deviations (same as generate_data)
  standard_deviations <- c(
    baseline_param$sd[baseline_param$cat == "bm"],
    baseline_param$sd[baseline_param$cat == "BL"]
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "tv"], num_timepoints)
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "pb"], num_timepoints) * trial_design$e
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "br"], num_timepoints)
  )

  # Calculate component-specific half-lives for enhanced carryover (same as generate_data)
  component_halflives <- calculate_component_halflives(model_param$carryover_t1half)

  # Calculate means for biomarker correlation logic with enhanced carryover
  means <- c(
    baseline_param$m[baseline_param$cat == "bm"],
    baseline_param$m[baseline_param$cat == "BL"]
  )

  # Calculate means for each factor with enhanced carryover
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

      time_variant_means <- apply_carryover_to_component(
        time_variant_means, trial_data, component_halflives, "tv"
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

      pharm_biomarker_means <- apply_carryover_to_component(
        pharm_biomarker_means, trial_data, component_halflives, "pb"
      )

      means <- c(means, pharm_biomarker_means)
    }

    if (factor_abbrev == "br") {
      br_result <- calculate_bio_response_with_interaction(
        trial_data, model_param, resp_param, component_halflives
      )
      
      bio_response_means <- br_result$bio_response_means
      bio_response_test <- br_result$bio_response_test

      means <- c(means, bio_response_means)
    }
  }

  # Build correlation matrix with full biomarker logic
  correlations <- build_correlation_matrix(labels, trial_design, model_param, num_timepoints,
                                          factor_types, factor_abbreviations,
                                          bio_response_test, bio_response_means, means,
                                          trial_data = trial_data, lambda_cor = lambda_cor)

  # Convert to covariance matrix
  sigma <- outer(standard_deviations, standard_deviations) * correlations

  # Check positive definiteness - REJECT non-PD matrices
  is_pd <- corpcor::is.positive.definite(sigma)

  if (verbose) {
    cat("Sigma matrix positive definite:", is_pd, "\n")
  }

  if (!is_pd) {
    if (verbose) {
      eigenvalues <- eigen(sigma, only.values = TRUE)$values
      neg_evals <- eigenvalues[eigenvalues < 0]
      warning(sprintf(
        'Non-PD covariance matrix corrected (min eigenvalue: %.4f, %d negative)',
        min(eigenvalues), length(neg_evals)))
    }
    sigma <- corpcor::make.positive.definite(sigma, tol = 1e-3)
  }

  chol_sigma <- tryCatch(chol(sigma), error = function(e) {
    chol(corpcor::make.positive.definite(sigma, tol = 1e-3))
  })

  list(
    sigma = sigma,
    labels = labels,
    standard_deviations = standard_deviations,
    correlations = correlations,
    means = means,
    nP = num_timepoints,
    cl = factor_abbreviations,
    chol_sigma = chol_sigma
  )
}

#===========================================================================
# Function: validate_correlation_structure
# Description: Validate that correlation parameters produce a positive
# definite covariance matrix with good conditioning
#
# Returns: TRUE if valid, FALSE if non-positive definite
# Side effects: Prints diagnostic information about eigenvalues and
#               condition number
#===========================================================================
validate_correlation_structure <- function(model_params,
                                          resp_param,
                                          baseline_param,
                                          trial_design) {
  # Build sigma matrix
  sigma_result <- build_sigma_matrix(
    model_params, resp_param, baseline_param,
    trial_design,
    factor_types = c("tv", "pb",
                     "br"),
    factor_abbreviations = c("tv", "pb", "br"),
    verbose = TRUE
  )

  if (is.null(sigma_result)) {
    cat("FAILED: Non-positive definite\n")
    return(FALSE)
  }

  # Check condition number
  sigma <- sigma_result$sigma
  eigenvalues <- eigen(sigma, only.values = TRUE)$values
  condition_number <- max(eigenvalues) / min(eigenvalues)

  cat("Eigenvalue range: [", min(eigenvalues), ", ",
      max(eigenvalues), "]\n")
  cat("Condition number: ", condition_number, "\n")

  # Well-conditioned if condition number < 100
  if (condition_number > 100) {
    warning("Matrix is ill-conditioned (condition number = ",
            condition_number, ")")
  }

  return(TRUE)
}

create_sigma_cache_key <- function(design_name, params) {
  # Create unique key for each parameter/design combination
  paste(design_name,
        params$n_participants,
        params$biomarker_correlation,
        params$carryover_t1half,
        sep = "_")
}

#===========================================================================
# Function: validate_parameter_grid
# Description: Pre-simulation validation of all parameter combinations
#===========================================================================
validate_parameter_grid <- function(param_grid,
                                    trial_design,
                                    model_params,
                                    resp_param,
                                    baseline_param,
                                    verbose = TRUE) {
  #' Validate all parameter combinations before simulation starts
  #'
  #' This function tests each parameter combination for positive definiteness
  #' and returns a detailed report of problematic combinations.
  #'
  #' @param param_grid Tibble of parameter combinations to test
  #' @param trial_design Trial design tibble with timepoint information
  #' @param model_params Base model parameters (fixed correlations)
  #' @param resp_param Response parameters (variances, effects)
  #' @param baseline_param Baseline parameters (biomarker, baseline)
  #' @param verbose If TRUE, print detailed validation report
  #'
  #' @return List with:
  #'   - valid_combinations: Tibble of combinations that pass validation
  #'   - invalid_combinations: Tibble of combinations that fail
  #'   - summary: Validation summary statistics
  #'   - report: Character vector of validation messages

  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("PRE-SIMULATION PARAMETER VALIDATION\n")
  cat(strrep("=", 80), "\n\n")

  # Initialize tracking lists
  valid_rows <- c()
  invalid_rows <- c()
  invalid_reasons <- character()
  condition_numbers <- numeric()

  # Test each parameter combination
  for (i in 1:nrow(param_grid)) {
    current_params <- param_grid[i, ]

    # Update model_params with current biomarker correlation
    test_params <- model_params
    test_params$c.bm <- current_params$biomarker_correlation

    # Try to build sigma matrix
    sigma_result <- tryCatch({
      build_sigma_matrix(
        test_params, resp_param, baseline_param, trial_design,
        factor_types = c("tv", "pb", "br"),
        factor_abbreviations = c("tv", "pb", "br"),
        verbose = FALSE
      )
    }, error = function(e) {
      return(NULL)
    })

    # Check result
    if (is.null(sigma_result)) {
      invalid_rows <- c(invalid_rows, i)
      invalid_reasons <- c(
        invalid_reasons,
        sprintf("Row %d: Non-positive definite matrix (c.bm=%.2f)",
                i, current_params$biomarker_correlation)
      )
    } else {
      # Compute condition number
      sigma <- sigma_result$sigma
      eigenvalues <- eigen(sigma, only.values = TRUE)$values
      kappa <- max(eigenvalues) / max(min(eigenvalues), 1e-10)
      condition_numbers[i] <- kappa

      # Flag if ill-conditioned (kappa > 100)
      if (kappa > 100) {
        cat(sprintf(
          "⚠ WARNING [Row %d]: Ill-conditioned matrix (κ = %.1f, c.bm = %.2f)\n",
          i, kappa, current_params$biomarker_correlation
        ))
        valid_rows <- c(valid_rows, i)
      } else {
        valid_rows <- c(valid_rows, i)
      }
    }
  }

  # Compile results
  valid_grid <- param_grid[valid_rows, ]
  invalid_grid <- param_grid[invalid_rows, ]

  # Print summary
  cat("\n", strrep("=", 80), "\n")
  cat("VALIDATION SUMMARY\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("Total combinations tested: %d\n", nrow(param_grid)))
  cat(sprintf("✓ Valid combinations:       %d (%.1f%%)\n",
              nrow(valid_grid), 100*nrow(valid_grid)/nrow(param_grid)))
  cat(sprintf("✗ Invalid combinations:     %d (%.1f%%)\n\n",
              nrow(invalid_grid), 100*nrow(invalid_grid)/nrow(param_grid)))

  # Print invalid combinations with reasons
  if (nrow(invalid_grid) > 0) {
    cat("PROBLEMATIC COMBINATIONS:\n")
    cat(strrep("-", 80), "\n")
    for (j in 1:length(invalid_reasons)) {
      cat(sprintf("  %s\n", invalid_reasons[j]))
    }
    cat("\n")

    # Print the invalid combinations table
    cat("Details of invalid combinations:\n")
    print(invalid_grid)
    cat("\n")
  }

  # Condition number statistics
  if (length(condition_numbers) > 0) {
    cond_valid <- condition_numbers[valid_rows]
    cat(sprintf(
      "Condition number (κ) statistics for valid combinations:\n  Mean:   %.1f\n  Median: %.1f\n  Min:    %.1f (best conditioned)\n  Max:    %.1f (worst conditioned)\n\n",
      mean(cond_valid, na.rm = TRUE),
      median(cond_valid, na.rm = TRUE),
      min(cond_valid, na.rm = TRUE),
      max(cond_valid, na.rm = TRUE)
    ))
  }

  # Recommendations
  if (nrow(invalid_grid) > 0) {
    cat("RECOMMENDATIONS:\n")
    cat(strrep("-", 80), "\n")

    # Check if it's a biomarker correlation issue
    if (all(grepl("c.bm", invalid_reasons))) {
      cat("  • All failures are due to biomarker correlation (c.bm) values being too high\n")
      cat("  • Consider reducing the maximum biomarker_correlation in param_grid\n")
      cat("  • Current constraint: c.bm ≤ 0.6\n")
      cat("  • You may need to reduce to: c.bm ≤ 0.4 or lower\n\n")
    }

    cat("  • Alternatively, adjust correlation structure parameters:\n")
    cat("    - Reduce c.cf1t (same-time cross-correlation, currently 0.2)\n")
    cat("    - Reduce c.cfct (different-time cross-correlation, currently 0.1)\n")
    cat("    - Increase c.autocorr (would require modifying Hendrickson parameters)\n\n")
  }

  # Final status
  cat(strrep("=", 80), "\n")
  if (nrow(invalid_grid) == 0) {
    cat("✓ ALL PARAMETER COMBINATIONS ARE VALID\n")
    cat("  Ready to proceed with simulation!\n")
  } else {
    cat("⚠ SOME PARAMETER COMBINATIONS ARE INVALID\n")
    cat(sprintf("  %d combinations excluded from simulation\n", nrow(invalid_grid)))
  }
  cat(strrep("=", 80), "\n\n")

  # Return results
  return(list(
    valid_combinations = valid_grid,
    invalid_combinations = invalid_grid,
    n_valid = nrow(valid_grid),
    n_invalid = nrow(invalid_grid),
    condition_numbers = condition_numbers[valid_rows],
    invalid_reasons = invalid_reasons
  ))
}

#===========================================================================
# Function: report_parameter_validation
# Description: Print detailed validation report with recommendations
#===========================================================================
report_parameter_validation <- function(validation_result, param_grid) {
  #' Generate a detailed validation report
  #'
  #' @param validation_result Output from validate_parameter_grid()
  #' @param param_grid Original parameter grid tested

  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("DETAILED PARAMETER VALIDATION REPORT\n")
  cat(strrep("=", 80), "\n\n")

  cat("GRID COMPOSITION:\n")
  cat(sprintf("  • n_participants:        %s\n",
              paste(unique(param_grid$n_participants), collapse = ", ")))
  cat(sprintf("  • biomarker_correlation: %s\n",
              paste(unique(param_grid$biomarker_correlation), collapse = ", ")))
  cat(sprintf("  • carryover_t1half:      %s\n",
              paste(unique(param_grid$carryover_t1half), collapse = ", ")))

  cat("\nVALIDATION RESULTS:\n")
  cat(sprintf("  ✓ Valid:   %3d combinations (%5.1f%%)\n",
              validation_result$n_valid,
              100 * validation_result$n_valid / nrow(param_grid)))
  cat(sprintf("  ✗ Invalid: %3d combinations (%5.1f%%)\n",
              validation_result$n_invalid,
              100 * validation_result$n_invalid / nrow(param_grid)))

  if (validation_result$n_invalid > 0) {
    cat("\nINVALID COMBINATIONS:\n")
    for (reason in validation_result$invalid_reasons) {
      cat(sprintf("  • %s\n", reason))
    }
  }

  if (length(validation_result$condition_numbers) > 0) {
    cat("\nNUMERICAL STABILITY:\n")
    cond_nums <- validation_result$condition_numbers
    cat(sprintf("  • κ (condition number) range: [%.1f, %.1f]\n",
                min(cond_nums), max(cond_nums)))
    cat(sprintf("  • Mean κ: %.1f (geometric: %.1f)\n",
                mean(cond_nums),
                exp(mean(log(cond_nums)))))
    cat(sprintf("  • Status: %s\n",
                if(max(cond_nums) < 100) "✓ Well-conditioned" else "⚠ Some ill-conditioning detected"))
  }

  cat("\n")
}

#===========================================================================
# Helper function for processing participant data
#===========================================================================
process_participant_data <- function(participant_data, timepoint_names, factor_abbreviations, N) {
  # Add participant ID column
  participant_data <- participant_data %>%
    mutate(ptID = 1:N)

  # Pre-allocate all result columns to avoid growing the data frame
  new_cols <- c(paste0("D_", timepoint_names), timepoint_names)
  zeros_df <- as_tibble(matrix(0, nrow = N, ncol = length(new_cols)))
  colnames(zeros_df) <- new_cols
  participant_data <- bind_cols(participant_data, zeros_df)

  # Calculate deltas (sums of factors) for all timepoints
  for (timepoint_idx in 1:length(timepoint_names)) {
    timepoint_name <- timepoint_names[timepoint_idx]
    delta_col <- paste0("D_", timepoint_name)
    components <- paste(timepoint_name, factor_abbreviations, sep = ".")

    # Calculate delta for all participants at once
    participant_data <- participant_data %>%
      mutate(!!delta_col := rowSums(select(., all_of(components))))
  }

  # Calculate timepoint scores from baseline and factors
  for (timepoint_idx in 1:length(timepoint_names)) {
    timepoint_name <- timepoint_names[timepoint_idx]
    components <- paste(timepoint_name, factor_abbreviations, sep = ".")

    # Calculate timepoint value from baseline and factors
    participant_data <- participant_data %>%
      mutate(!!timepoint_name := BL - rowSums(select(., all_of(components))))
  }

  return(participant_data)
}

#===========================================================================
# Function: lme_analysis
# Matches orig lme_analysis conditional formula logic exactly
#===========================================================================
lme_analysis <- function(trial_design_set, dat, op) {

  if (missing(op)) {
    op <- list()
    op$useDE <- TRUE
    op$t_random_slope <- FALSE
    op$full_model_out <- FALSE
    op$carryover_t1half <- 0
    op$simplecarryover <- FALSE
    op$carryover_scalefactor <- 1
  } else if (!('simplecarryover' %in% names(op))) {
    op$simplecarryover <- FALSE
  }
  if (!('carryover_t1half' %in% names(op))) {
    op$carryover_t1half <- 0
    op$carryover_scalefactor <- 1
  }
  if ((op$carryover_t1half > 0) & (op$simplecarryover == TRUE)) {
    stop('Cannot use both simplecarryover and carryover_t1half')
  }

  n_groups <- length(trial_design_set)
  datout <- vector('list', n_groups)
  last_ptID <- 0

  for (g in seq_len(n_groups)) {
    td <- trial_design_set[[g]]

    dat_single <- dat |>
      filter(path == g) |>
      as_tibble()

    td_names <- if ('timeptnames' %in% names(td)) {
      td$timeptnames
    } else if ('timeptname' %in% names(td)) {
      td$timeptname
    } else {
      td$timepoint_name
    }

    td_with_bl <- bind_rows(
      tibble(
        timeptnames = 'BL', t_wk = 0, e = 0,
        tod = 0, tsd = 0, tpb = 0
      ),
      as_tibble(td) |>
        rename(any_of(c(timeptnames = 'timeptname',
                        timeptnames = 'timepoint_name')))
    ) |>
      mutate(t = cumsum(t_wk))

    all_names <- c('BL', td_names)
    valid_names <- all_names[all_names %in% names(dat_single)]

    data_wide <- dat_single |>
      select(ptID, bm, all_of(valid_names)) |>
      mutate(ptID = ptID + last_ptID)
    last_ptID <- max(data_wide$ptID)

    data_long <- data_wide |>
      pivot_longer(
        cols = all_of(valid_names),
        names_to = 'timeptnames',
        values_to = 'Sx',
        values_drop_na = FALSE
      ) |>
      left_join(
        td_with_bl |> select(timeptnames, t, De = e, tod, tsd),
        by = 'timeptnames'
      ) |>
      mutate(Db = (tod > 0))

    data_long <- data_long |>
      mutate(
        Dbc = case_when(
          Db ~ 1,
          TRUE ~ (1/2)^(op$carryover_scalefactor * tsd /
                        op$carryover_t1half)
        )
      )

    datout[[g]] <- data_long
  }

  datamerged <- bind_rows(datout)

  # Test 1: can we include expectancy?
  td_last <- trial_design_set[[n_groups]]
  e_vals <- if ('e' %in% names(td_last)) td_last$e else rep(0.5, nrow(td_last))
  var_in_exp <- length(unique(e_vals))

  # Test 2: does Db vary within participants?
  datamerged <- datamerged |>
    group_by(ptID) |>
    mutate(
      mean_Db = mean(Db[t > 0], na.rm = TRUE),
      Db_var = (mean_Db != 0) & (mean_Db != 1)
    ) |>
    ungroup()

  var_in_Db <- sum(datamerged$Db_var[datamerged$t > 0], na.rm = TRUE) > 0

  if (!var_in_Db) {
    ever_on_drug <- datamerged |>
      filter(t > 0, mean_Db == 1) |>
      pull(ptID) |>
      unique()
    datamerged <- datamerged |>
      filter(ptID %in% ever_on_drug)
  }

  # Test 3: is tsd ever non-zero?
  var_in_tsd <- any(datamerged$tsd != 0, na.rm = TRUE)

  # Build formula (matches orig exactly)
  model_base <- 'Sx ~ bm + t'
  if (var_in_Db) {
    model_base <- paste0(model_base, ' + Dbc + bm:Dbc')
  } else {
    model_base <- paste0(model_base, ' + bm:t')
  }
  if ((var_in_exp > 1) & (op$useDE == TRUE)) {
    model_base <- paste0(model_base, ' + De')
  }
  if (op$simplecarryover & var_in_tsd) {
    model_base <- paste0(model_base, ' + tsd')
  }
  form <- as.formula(model_base)

  # Remove NA rows for model variables
  model_vars <- c('Sx', 'bm', 't', 'ptID')
  if (var_in_Db) model_vars <- c(model_vars, 'Dbc')
  if (op$simplecarryover & var_in_tsd) model_vars <- c(model_vars, 'tsd')
  if ((var_in_exp > 1) & (op$useDE == TRUE)) model_vars <- c(model_vars, 'De')
  for (mv in model_vars) {
    datamerged <- datamerged |> filter(!is.na(.data[[mv]]))
  }

  # Fit nlme::lme with corCAR1 (matches orig)
  fit <- tryCatch(
    nlme::lme(form, random = ~1 | ptID,
      correlation = nlme::corCAR1(form = ~t | ptID),
      data = datamerged,
      control = nlme::lmeControl(
        opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) {
      tryCatch(
        nlme::lme(form, random = ~1 | ptID,
          data = datamerged,
          control = nlme::lmeControl(opt = 'optim')),
        error = function(e2) NULL
      )
    }
  )

  hold_warning <- as.character(NA)
  issingular <- FALSE

  if (is.null(fit)) {
    out <- tibble(
      beta = NA_real_, betaSE = NA_real_,
      p = NA_real_, issingular = NA,
      warning = 'model failed to converge'
    )
    if (op$full_model_out) {
      out <- list(form = form, fit = NULL,
                  datamerged = datamerged, stdout = out)
    }
    return(out)
  }

  # Extract coefficients (matches orig tTable extraction)
  cc <- summary(fit)$tTable
  coefnames <- rownames(cc)

  if (var_in_Db) {
    target <- intersect(c('bm:Dbc', 'Dbc:bm'), coefnames)
    if (length(target) == 0) {
      beta <- NA_real_; betaSE <- NA_real_; p <- NA_real_
    } else {
      target <- target[1]
      p <- cc[target, 'p-value']
      beta <- cc[target, 'Value']
      betaSE <- cc[target, 'Std.Error']
    }
  } else {
    target <- intersect(c('bm:t', 't:bm'), coefnames)
    if (length(target) == 0) {
      beta <- NA_real_; betaSE <- NA_real_; p <- NA_real_
    } else {
      target <- target[1]
      p <- cc[target, 'p-value']
      beta <- cc[target, 'Value']
      betaSE <- cc[target, 'Std.Error']
    }
  }

  out <- tibble(
    beta = beta, betaSE = betaSE, p = p,
    issingular = issingular, warning = hold_warning
  )

  if (op$full_model_out) {
    out <- list(form = form, fit = fit,
                datamerged = datamerged, stdout = out)
  }

  out
}


#===========================================================================
#===========================================================================
# Step 7a: build_trial_design (tidyverse port of buildtrialdesign)
#===========================================================================
build_trial_design <- function(name_longform = 'Trial Design 1',
                               name_shortform = name_longform,
                               timepoints,
                               timeptnames = paste0('V', seq_along(timepoints)),
                               expectancies,
                               ondrug) {

  n_paths <- length(ondrug)
  t_wk <- c(timepoints[1],
             diff(timepoints))

  trialpaths <- map(ondrug, function(od) {
    od_intervals <- t_wk * od
    tod <- od_intervals
    tsd <- t_wk * (od != 1)
    tpb <- t_wk * (expectancies > 0)

    for (i in 2:length(timepoints)) {
      if (od_intervals[i] > 0) {
        tod[i] <- tod[i] + tod[i - 1]
      } else {
        tsd[i] <- tsd[i] + tsd[i - 1]
      }
      if (tpb[i] > 0) {
        tpb[i] <- tpb[i] + tpb[i - 1]
      }
    }

    ever_on_drug <- (cumsum(od) > 0)
    tsd <- ever_on_drug * tsd

    tibble(
      timeptnames = timeptnames,
      t_wk = t_wk,
      e = expectancies,
      tod = tod,
      tsd = tsd,
      tpb = tpb
    )
  })

  metadata <- list(
    name_longform = name_longform,
    name_shortform = name_shortform,
    timepoints = timepoints,
    timeptnames = timeptnames,
    expectancies = expectancies,
    ondrug = ondrug
  )

  list(metadata = metadata, trialpaths = trialpaths)
}


#===========================================================================
# Step 7b: censor_data (tidyverse port of censordata)
#===========================================================================
censor_data <- function(dat, trialdesign, censorparam) {
  td <- as_tibble(trialdesign)
  td$t_wk_cumulative <- cumsum(td$t_wk)
  n_tp <- nrow(td)

  tp_names <- if ('timeptnames' %in% names(td)) {
    td$timeptnames
  } else if ('timeptname' %in% names(td)) {
    td$timeptname
  } else {
    td$timepoint_name
  }

  delta_cols <- paste0('D_', tp_names)
  cdt <- dat |> select(all_of(delta_cols))

  frac_NA <- censorparam$beta0
  frac_NA_biased <- censorparam$beta1
  fna1 <- frac_NA * (1 - frac_NA_biased)
  fna2 <- frac_NA * frac_NA_biased

  cdt_ps1 <- matrix(1, nrow = nrow(cdt), ncol = ncol(cdt))
  cdt_ps2 <- sapply(cdt, function(x) (x + 100)^censorparam$eb1)
  cdt_p1 <- t(t(cdt_ps1) * td$t_wk)
  cdt_p2 <- t(t(cdt_ps2) * td$t_wk)

  nr <- nrow(cdt_p1)
  nc <- ncol(cdt_p1)
  cdt_r1 <- runif(nr * nc, min = 0,
                   max = 2 * mean(cdt_p1) * (0.5 / fna1))
  cdt_r2 <- runif(nr * nc, min = 0,
                   max = 2 * mean(cdt_p2) * (0.5 / fna2))

  do1 <- cdt_p1 > matrix(cdt_r1, nrow = nr, ncol = nc)
  do2 <- cdt_p2 > matrix(cdt_r2, nrow = nr, ncol = nc)
  do_mat <- do1 | do2

  for (itp in seq_along(tp_names)) {
    for (tp_idx in 1:itp) {
      masking_col <- tp_idx
      masked_tp <- tp_names[itp]
      mask <- do_mat[, masking_col]
      dat[[masked_tp]][mask] <- NA
    }
  }

  dat
}


#===========================================================================
# Step 7c: generate_simulated_results
# (tidyverse port of generateSimulatedResults)
#===========================================================================
generate_simulated_results <- function(
    trialdesigns, respparamsets, blparamsets,
    censorparams, modelparams, simparam,
    analysisparams, rawdataout = FALSE,
    lambda_cor = NA, n_cores = 1) {

  if (missing(analysisparams)) {
    analysisparams <- list(useDE = TRUE, t_random_slope = FALSE,
                           full_model_out = FALSE)
  }

  if (n_cores < 1) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  initial_dir <- getwd()

  vpg_master <- expand.grid(
    trialdesign = seq_along(trialdesigns),
    respparamset = seq_along(respparamsets),
    blparamset = seq_along(blparamsets),
    modelparamset = seq_len(nrow(modelparams))
  )
  n_combos <- nrow(vpg_master)

  cat('Caching sigma matrices...\n')
  sigma_cache <- list()
  for (i_r in seq_len(n_combos)) {
    td <- trialdesigns[[vpg_master[i_r, 'trialdesign']]][[2]]
    rp <- respparamsets[[vpg_master[i_r, 'respparamset']]]$param
    bpp <- blparamsets[[vpg_master[i_r, 'blparamset']]]$param
    mpp <- modelparams[vpg_master[i_r, 'modelparamset'], ]
    n_paths <- length(td)

    for (i_p in seq_len(n_paths)) {
      cache_key <- paste(vpg_master[i_r, 'trialdesign'],
                         vpg_master[i_r, 'respparamset'],
                         vpg_master[i_r, 'blparamset'],
                         vpg_master[i_r, 'modelparamset'],
                         i_p, sep = '_')
      if (is.null(sigma_cache[[cache_key]])) {
        td_path <- td[[i_p]]
        if (!'timepoint_name' %in% names(td_path)) {
          if ('timeptnames' %in% names(td_path)) {
            td_path$timepoint_name <- td_path$timeptnames
          }
        }
        sigma_cache[[cache_key]] <- build_sigma_matrix(
          mpp, rp, bpp, td_path,
          factor_types = c('tv', 'pb', 'br'),
          factor_abbreviations = c('tv', 'pb', 'br'),
          lambda_cor = lambda_cor
        )
      }
    }
  }
  cat(sprintf('Cached %d unique sigma matrices\n\n',
              length(sigma_cache)))

  # Progressive save setup
  if (!simparam$progressiveSave) {
    n_large_loops <- 1
    ll_starts <- 1
    ll_stops <- n_combos
  } else {
    n_large_loops <- ceiling(n_combos / simparam$nRep2save)
    if (n_large_loops > 1) {
      ll_starts <- c(1, 1 + simparam$nRep2save *
                       (1:(n_large_loops - 1)))
      ll_stops <- c(ll_starts[2:n_large_loops] - 1, n_combos)
    } else {
      ll_starts <- 1
      ll_stops <- n_combos
    }
  }

  for (i_ll in simparam$saveunit2start:n_large_loops) {
    vpg <- vpg_master[ll_starts[i_ll]:ll_stops[i_ll], ,
                      drop = FALSE]
    i_paramset <- ll_starts[i_ll]
    n_r <- nrow(vpg)

    run_one <- function(i_r) {
      td <- trialdesigns[[vpg[i_r, 'trialdesign']]][[2]]
      rp <- respparamsets[[vpg[i_r, 'respparamset']]]$param
      bpp <- blparamsets[[vpg[i_r, 'blparamset']]]$param
      mpp <- modelparams[vpg[i_r, 'modelparamset'], ]

      n_paths <- length(td)
      Ns <- mpp$N %/% n_paths
      Ns <- Ns + c(rep(1, mpp$N %% n_paths),
                    rep(0, n_paths - mpp$N %% n_paths))
      NNs <- Ns * simparam$Nreps

      mpp_copy <- mpp
      mpp_copy$N <- NNs[[1]]
      ck <- paste(vpg[i_r, 'trialdesign'],
                  vpg[i_r, 'respparamset'],
                  vpg[i_r, 'blparamset'],
                  vpg[i_r, 'modelparamset'],
                  1, sep = '_')

      td_path <- td[[1]]
      if (!'timepoint_name' %in% names(td_path)) {
        if ('timeptnames' %in% names(td_path)) {
          td_path$timepoint_name <- td_path$timeptnames
        }
      }

      dat <- generate_data(mpp_copy, rp, bpp, td_path,
                           empirical = FALSE,
                           make_positive_definite = TRUE,
                           cached_sigma = sigma_cache[[ck]])
      dat$path <- 1
      dat$replicate <- rep(1:simparam$Nreps, Ns[1])

      if (n_paths > 1) {
        for (i_p in 2:n_paths) {
          mpp_copy$N <- NNs[[i_p]]
          ck <- paste(vpg[i_r, 'trialdesign'],
                      vpg[i_r, 'respparamset'],
                      vpg[i_r, 'blparamset'],
                      vpg[i_r, 'modelparamset'],
                      i_p, sep = '_')
          td_p <- td[[i_p]]
          if (!'timepoint_name' %in% names(td_p)) {
            if ('timeptnames' %in% names(td_p)) {
              td_p$timepoint_name <- td_p$timeptnames
            }
          }
          dat2 <- generate_data(mpp_copy, rp, bpp, td_p,
                                empirical = FALSE,
                                make_positive_definite = TRUE,
                                cached_sigma = sigma_cache[[ck]])
          dat2$path <- i_p
          dat2$replicate <- rep(1:simparam$Nreps, Ns[i_p])
          dat <- bind_rows(dat, dat2)
        }
      }
      mpp_copy$N <- sum(Ns)

      results_list <- list()
      ridx <- 0

      for (i_ap in seq_len(nrow(analysisparams))) {
        ap_row <- as.list(analysisparams[i_ap, ])
        et_out <- lme_analysis(td, dat, ap_row)
        for (i_s in seq_len(simparam$Nreps)) {
          dat_rep <- dat |> filter(replicate == i_s)
          a_out <- lme_analysis(td, dat_rep, ap_row)
          ridx <- ridx + 1
          results_list[[ridx]] <- bind_cols(
            as_tibble(vpg[i_r, , drop = FALSE]),
            as_tibble(as.list(mpp_copy)),
            tibble(
              censorparamset = 0,
              use_DE = ap_row$useDE,
              t_random_slope = ap_row$t_random_slope,
              irep = (i_r - 1) * simparam$Nreps + i_s,
              frac_NA = 0,
              ETbeta = et_out$beta,
              ETbetaSE = et_out$betaSE,
              beta = a_out$beta,
              betaSE = a_out$betaSE,
              p = a_out$p,
              issingular = a_out$issingular,
              warning = a_out$warning
            )
          )
        }

        nocensoring <- length(censorparams) < 2 &&
          is.na(censorparams)
        if (!nocensoring) {
          for (i_c in seq_len(nrow(censorparams))) {
            datc <- censor_data(dat, td[[1]], censorparams[i_c, ])
            for (i_s in seq_len(simparam$Nreps)) {
              datc_rep <- datc |> filter(replicate == i_s)
              frac_na <- sum(is.na(datc_rep)) /
                (mpp_copy$N * nrow(td[[1]]))
              a_out <- lme_analysis(td, datc_rep, ap_row)
              ridx <- ridx + 1
              results_list[[ridx]] <- bind_cols(
                as_tibble(vpg[i_r, , drop = FALSE]),
                as_tibble(as.list(mpp_copy)),
                tibble(
                  censorparamset = i_c,
                  use_DE = ap_row$useDE,
                  t_random_slope = ap_row$t_random_slope,
                  irep = (i_r - 1) * simparam$Nreps + i_s,
                  frac_NA = frac_na,
                  ETbeta = et_out$beta,
                  ETbetaSE = et_out$betaSE,
                  beta = a_out$beta,
                  betaSE = a_out$betaSE,
                  p = a_out$p,
                  issingular = a_out$issingular,
                  warning = a_out$warning
                )
              )
            }
          }
        }
      }
      bind_rows(results_list)
    }

    out_parts <- map(seq_len(n_r), function(i_r) {
      cat(sprintf('Parameter set %d of %d (%d of %d total)\n',
                  i_r, n_r, i_paramset + i_r - 1, n_combos))
      run_one(i_r)
    })
    out <- bind_rows(out_parts)

    outpt <- list(
      results = out,
      parameterselections = list(
        trialdesigns = trialdesigns,
        respparamsets = respparamsets,
        blparamsets = blparamsets,
        censorparams = censorparams,
        modelparams = modelparams,
        analysisparams = analysisparams,
        simparam = simparam
      )
    )

    if (simparam$progressiveSave) {
      setwd(simparam$savedir)
      saveRDS(outpt, paste(simparam$basesavename, i_ll,
                           sep = '_save'))
    }
  }

  setwd(initial_dir)
  if (!simparam$progressiveSave) outpt
}
