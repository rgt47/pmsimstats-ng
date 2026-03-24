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
# HENDRICKSON-ALIGNED FUNCTIONS (Phase 1)
# These functions replicate the exact DGP from Hendrickson et al. (2020)
# to enable faithful replication before extending with carryover modeling.
#===========================================================================

#' Modified Gompertz function (Hendrickson original)
#'
#' Offset-rescaled version that guarantees f(0) = 0 and asymptotes at
#' maxr. The pmsimstats2025 mod_gompertz() does NOT pass through the
#' origin; this version does, matching the original package exactly.
#'
#' @param t Time vector
#' @param maxr Maximum response value
#' @param disp Displacement parameter
#' @param rate Rate parameter
#' @return Numeric vector of response values, with f(0) = 0
modgompertz_orig <- function(t, maxr, disp, rate) {
  y <- maxr * exp(-disp * exp(-rate * t))
  vert_offset <- maxr * exp(-disp * exp(-rate * 0))
  y <- y - vert_offset
  y <- y * (maxr / (maxr - vert_offset))
  y
}

#' Hendrickson response and baseline parameters
#'
#' Extracted from the original pmsimstats package data objects
#' (extracted_rp.rda, extracted_bp.rda). These are the CAPS5/blood
#' pressure values used in the published simulations.
hendrickson_resp_params <- list(
  tv = list(max = 6.50647, disp = 5, rate = 0.35, sd = 10),
  pb = list(max = 6.50647, disp = 5, rate = 0.35, sd = 10),
  br = list(max = 10.98604, disp = 5, rate = 0.42, sd = 8)
)

hendrickson_bl_params <- list(
  BL = list(m = 83.069, sd = 18.483),
  bm = list(m = 124.328, sd = 15.362)
)

#' Hendrickson correlation defaults
hendrickson_corr_params <- list(
  c.tv = 0.7,
  c.pb = 0.7,
  c.br = 0.7,
  c.cf1t = 0.2,
  c.cfct = 0.1
)

#' Compute time-since-discontinuation (orig interval-based logic)
#'
#' In the original, tsd accumulates interval widths (t_wk) while
#' off drug, resets to 0 when on drug, and is zeroed for periods
#' before a participant was ever on drug (via everondrug mask).
#'
#' @param t_wk Numeric vector of interval widths (not cumulative)
#' @param ondrug Binary vector: 1 = on drug, 0 = off drug
#' @return Numeric vector of time-since-discontinuation values
compute_tsd_orig <- function(t_wk, ondrug) {
  nP <- length(ondrug)
  od_intervals <- t_wk * ondrug
  tsd <- t_wk * (ondrug != 1)
  for (iTP in 2:nP) {
    if (od_intervals[iTP] == 0) {
      tsd[iTP] <- tsd[iTP] + tsd[iTP - 1]
    }
  }
  everondrug <- (cumsum(ondrug) > 0)
  tsd * everondrug
}

#' Compute time-on-drug (orig interval-based logic)
#'
#' Accumulates interval widths while on drug, resets when off drug.
#' This differs from pmsimstats2025's cumsum(treatment) approach
#' which counts measurement occasions rather than elapsed time.
#'
#' @param t_wk Numeric vector of interval widths (not cumulative)
#' @param ondrug Binary vector: 1 = on drug, 0 = off drug
#' @return Numeric vector of time-on-drug values
compute_tod_orig <- function(t_wk, ondrug) {
  nP <- length(ondrug)
  od_intervals <- t_wk * ondrug
  tod <- od_intervals
  for (iTP in 2:nP) {
    if (od_intervals[iTP] > 0) {
      tod[iTP] <- tod[iTP] + tod[iTP - 1]
    }
  }
  tod
}

#' Compute time-in-placebo-belief (orig interval-based logic)
#'
#' Accumulates interval widths while expectancy > 0.
#'
#' @param t_wk Numeric vector of interval widths
#' @param expectancies Numeric vector of expectancy values (0-1)
#' @return Numeric vector of time-in-placebo-belief values
compute_tpb_orig <- function(t_wk, expectancies) {
  nP <- length(expectancies)
  tpb <- t_wk * (expectancies > 0)
  for (iTP in 2:nP) {
    if (tpb[iTP] > 0) {
      tpb[iTP] <- tpb[iTP] + tpb[iTP - 1]
    }
  }
  tpb
}

#' Biomarker-BR correlation with exponential decay
#'
#' On-drug timepoints receive the full c.bm correlation. Off-drug
#' timepoints with tsd > 0 receive exponentially decayed correlation:
#' c.bm * exp(-lambda_cor * tsd). Timepoints before any drug exposure
#' receive zero correlation.
#'
#' @param ondrug Binary vector: 1 = on drug, 0 = off drug
#' @param tsd Numeric vector of time-since-discontinuation
#' @param c.bm Base biomarker-BR correlation
#' @param lambda_cor Correlation decay rate
#' @param nP Number of timepoints
#' @return Numeric vector of per-timepoint BM-BR correlations
build_bm_br_correlations <- function(ondrug, tsd, c.bm,
                                     lambda_cor, nP) {
  corrs <- numeric(nP)
  for (p in 1:nP) {
    if (ondrug[p] == 1) {
      corrs[p] <- c.bm
    } else if (tsd[p] > 0 && lambda_cor > 0) {
      corrs[p] <- c.bm * exp(-lambda_cor * tsd[p])
    }
  }
  corrs
}

#' Build per-path sigma matrix (Hendrickson-aligned)
#'
#' Constructs the full covariance matrix for a single path through
#' the trial design, following the original's exact logic:
#' 1. Compute tod, tsd, tpb from ondrug vector
#' 2. Compute means via modgompertz_orig
#' 3. Apply carryover to BR means
#' 4. Build SDs with expectancy-scaled PB
#' 5. Build correlation matrix (compound symmetry within-factor,
#'    Ron Thomas BM-BR scaling)
#' 6. Sigma = outer(sds, sds) * correlations
#' 7. PD enforcement via corpcor
#'
#' @param timepoints Cumulative timepoint vector (e.g., c(4,8,9,...))
#' @param ondrug Binary vector for this path
#' @param expectancies Expectancy vector for this path
#' @param c.bm Biomarker-BR correlation
#' @param carryover_t1half Carryover half-life in weeks
#' @param lambda_cor Correlation decay rate (default NA, auto-computed)
#' @param resp_params Response parameters (default hendrickson_resp_params)
#' @param bl_params Baseline parameters (default hendrickson_bl_params)
#' @param corr_params Correlation parameters (default hendrickson_corr_params)
#' @param make_pd Force positive definiteness (default TRUE)
#' @return List with sigma, means, labels, sds, correlations, indices,
#'   and diagnostic info (brmeans, brtest, bm_br_corrs)
build_path_sigma <- function(
    timepoints,
    ondrug,
    expectancies,
    c.bm,
    carryover_t1half,
    lambda_cor = NA,
    resp_params = hendrickson_resp_params,
    bl_params = hendrickson_bl_params,
    corr_params = hendrickson_corr_params,
    make_pd = TRUE) {

  nP <- length(ondrug)

  # Auto-compute lambda_cor from carryover half-life
  if (is.na(lambda_cor)) {
    if (carryover_t1half > 0) {
      lambda_cor <- log(2) / carryover_t1half
    } else {
      lambda_cor <- 0
    }
  }

  t_wk <- c(timepoints[1],
            diff(timepoints))

  tod <- compute_tod_orig(t_wk, ondrug)
  tsd <- compute_tsd_orig(t_wk, ondrug)
  tpb <- compute_tpb_orig(t_wk, expectancies)
  t_cumulative <- cumsum(t_wk)

  timeptnames <- paste0('T', seq_len(nP))
  cl <- c('tv', 'pb', 'br')

  labels <- c(
    'bm', 'BL',
    paste(timeptnames, cl[1], sep = '.'),
    paste(timeptnames, cl[2], sep = '.'),
    paste(timeptnames, cl[3], sep = '.')
  )

  sds <- c(bl_params$bm$sd, bl_params$BL$sd)
  sds <- c(sds, rep(resp_params$tv$sd, nP))
  sds <- c(sds, rep(resp_params$pb$sd, nP) * expectancies)
  sds <- c(sds, rep(resp_params$br$sd, nP))

  tv_means <- modgompertz_orig(t_cumulative,
    resp_params$tv$max, resp_params$tv$disp, resp_params$tv$rate)
  pb_means <- modgompertz_orig(tpb,
    resp_params$pb$max, resp_params$pb$disp, resp_params$pb$rate) *
    expectancies

  brmeans <- modgompertz_orig(tod,
    resp_params$br$max, resp_params$br$disp, resp_params$br$rate)
  brtest <- brmeans == 0
  rawbrmeans <- brmeans

  if (nP > 1 && carryover_t1half > 0) {
    for (p in 2:nP) {
      if (ondrug[p] == 0 && tsd[p] > 0) {
        brmeans[p] <- brmeans[p] +
          brmeans[p - 1] * (1/2)^(tsd[p] / carryover_t1half)
      }
    }
  }

  means <- c(bl_params$bm$m, bl_params$BL$m,
             tv_means, pb_means, brmeans)

  correlations <- diag(length(labels))
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  for (cc in cl) {
    ac <- corr_params[[paste0('c.', cc)]]
    if (nP > 1) {
      for (p in 1:(nP - 1)) {
        for (p2 in (p + 1):nP) {
          n1 <- paste(timeptnames[p], cc, sep = '.')
          n2 <- paste(timeptnames[p2], cc, sep = '.')
          time_gap <- abs(t_cumulative[p2] - t_cumulative[p])
          correlations[n1, n2] <- ac^time_gap
          correlations[n2, n1] <- ac^time_gap
        }
      }
    }

    for (c2 in setdiff(cl, cc)) {
      for (p in 1:nP) {
        n1 <- paste(timeptnames[p], cc, sep = '.')
        n2 <- paste(timeptnames[p], c2, sep = '.')
        correlations[n1, n2] <- corr_params$c.cf1t
        correlations[n2, n1] <- corr_params$c.cf1t
      }
      if (nP > 1) {
        for (p in 1:(nP - 1)) {
          for (p2 in (p + 1):nP) {
            n1 <- paste(timeptnames[p], cc, sep = '.')
            n2 <- paste(timeptnames[p2], c2, sep = '.')
            time_gap <- abs(t_cumulative[p2] - t_cumulative[p])
            correlations[n1, n2] <- corr_params$c.cfct * ac^time_gap
            correlations[n2, n1] <- corr_params$c.cfct * ac^time_gap
          }
        }
      }
    }
  }

  bm_br_corrs <- build_bm_br_correlations(
    ondrug, tsd, c.bm = c.bm,
    lambda_cor = lambda_cor, nP = nP
  )
  for (p in 1:nP) {
    n1 <- paste(timeptnames[p], 'br', sep = '.')
    correlations['bm', n1] <- bm_br_corrs[p]
    correlations[n1, 'bm'] <- bm_br_corrs[p]
  }

  sigma <- outer(sds, sds) * correlations

  was_pd <- corpcor::is.positive.definite(sigma)
  if (make_pd && !was_pd) {
    sigma <- corpcor::make.positive.definite(sigma, tol = 1e-3)
  }

  bm_idx <- 1
  bl_idx <- 2
  tv_idx <- 3:(2 + nP)
  pb_idx <- (3 + nP):(2 + 2 * nP)
  br_idx <- (3 + 2 * nP):(2 + 3 * nP)

  list(
    sigma = sigma,
    means = means,
    labels = labels,
    sds = sds,
    correlations = correlations,
    indices = list(
      bm = bm_idx, bl = bl_idx,
      tv = tv_idx, pb = pb_idx, br = br_idx,
      nP = nP
    ),
    path_info = list(
      ondrug = ondrug,
      tod = tod,
      tsd = tsd,
      tpb = tpb,
      t_wk = t_wk,
      t_cumulative = t_cumulative,
      timeptnames = timeptnames,
      expectancies = expectancies
    ),
    diagnostics = list(
      brmeans = brmeans,
      rawbrmeans = rawbrmeans,
      brtest = brtest,
      bm_br_corrs = bm_br_corrs,
      was_pd = was_pd,
      eigenvalues = eigen(sigma, only.values = TRUE)$values
    )
  )
}
