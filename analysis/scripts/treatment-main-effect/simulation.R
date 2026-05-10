# N-of-1 vs RCT Simulation: Prazosin for PTSD Nightmares
# Comparing detection of treatment effects between study designs
# Enhanced version with Morris et al. (2019) compliant reporting
#
# Features added for Morris compliance:
# - Monte Carlo standard errors (MCSE) for all performance measures
# - Type I error evaluation under null hypothesis
# - nsim justification calculations
# - Explicit carryover modeling option (vs exclusion strategy)
# - Parameter validation functions

library(dplyr)
library(tidyr)
library(parallel)

set.seed(42)

# Detect available cores for parallel processing (use n-1 to leave one free)
n_cores <- max(1, detectCores() - 1)

#==========================================================================
# QUICK TEST MODE
#==========================================================================
# Set quick_test <- TRUE to run with reduced simulation counts for testing
# Set quick_test <- FALSE for full production runs.
#
# Default is FALSE so that running the script as-is reproduces the
# Morris-justified n_sim values quoted in the report. Setting to TRUE
# yields a small smoke run for development iteration.

quick_test <- FALSE

if (quick_test) {

  cat("*** QUICK TEST MODE ENABLED ***\n")
  cat("Using reduced simulation counts for faster execution.\n\n")
}

#==========================================================================
# MORRIS ET AL. (2019) COMPLIANT FUNCTIONS
#==========================================================================

#' Calculate Monte Carlo Standard Error for proportion (power/Type I error)
#'
#' @param p Estimated proportion (power or Type I error rate)
#' @param n_sim Number of simulations
#' @return Monte Carlo SE
calculate_mcse_proportion <- function(p, n_sim) {
  sqrt(p * (1 - p) / n_sim)
}

#' Calculate Monte Carlo Standard Error for mean (bias)
#'
#' @param estimates Vector of effect estimates
#' @return Monte Carlo SE for bias
calculate_mcse_bias <- function(estimates) {
  sqrt(var(estimates) / length(estimates))
}

#' Calculate Monte Carlo Standard Error for empirical SE
#'
#' @param estimates Vector of effect estimates
#' @return Monte Carlo SE for empirical SE
calculate_mcse_emp_se <- function(estimates) {
  n <- length(estimates)
  emp_se <- sd(estimates)
  emp_se / sqrt(2 * (n - 1))
}

#' Calculate required nsim for target MCSE (Morris Equation 1)
#'
#' @param expected_power Expected power (proportion)
#' @param target_mcse Target Monte Carlo SE
#' @return Required number of simulations
calculate_required_nsim <- function(expected_power, target_mcse = 0.01) {
  ceiling(expected_power * (1 - expected_power) / target_mcse^2)
}

#' Comprehensive simulation performance summary (Morris Table 6)
#'
#' @param estimates Vector of effect estimates
#' @param true_value True parameter value
#' @param p_values Vector of p-values
#' @param alpha Significance level
#' @param model_se Vector of model-based SEs (optional)
#' @return List with all performance measures and MCSEs
summarize_simulation_performance <- function(estimates, true_value, p_values,
                                             alpha = 0.05, model_se = NULL) {
  n_sim <- length(estimates)

  # Bias
  bias <- mean(estimates) - true_value
  bias_mcse <- calculate_mcse_bias(estimates)

empirical_se <- sd(estimates)
  emp_se_mcse <- calculate_mcse_emp_se(estimates)

  # MSE
  mse <- mean((estimates - true_value)^2)
  mse_mcse <- sqrt(var((estimates - true_value)^2) / n_sim)

  # Power (rejection rate)
  power <- mean(p_values < alpha)
  power_mcse <- calculate_mcse_proportion(power, n_sim)

  # 95% CI for power
  power_ci_lower <- power - 1.96 * power_mcse
  power_ci_upper <- power + 1.96 * power_mcse

  result <- list(
    n_sim = n_sim,
    true_value = true_value,
    mean_estimate = mean(estimates),
    bias = bias,
    bias_mcse = bias_mcse,
    empirical_se = empirical_se,
    emp_se_mcse = emp_se_mcse,
    mse = mse,
    mse_mcse = mse_mcse,
    power = power,
    power_mcse = power_mcse,
    power_ci_lower = power_ci_lower,
    power_ci_upper = power_ci_upper
  )

  # Coverage (if model SE provided)
  if (!is.null(model_se) && length(model_se) == n_sim) {
    ci_lower <- estimates - 1.96 * model_se
    ci_upper <- estimates + 1.96 * model_se
    coverage <- mean(ci_lower <= true_value & ci_upper >= true_value)
    coverage_mcse <- calculate_mcse_proportion(coverage, n_sim)
    result$coverage <- coverage
    result$coverage_mcse <- coverage_mcse
  }

  result
}

#' Validate parameter combinations for positive definiteness
#'
#' @param sigma_between Between-individual variance
#' @param sigma_within Within-individual variance
#' @param n_periods Number of periods
#' @return List with validation results
validate_parameters <- function(sigma_between, sigma_within, n_periods) {
  # Check basic parameter validity
  valid <- TRUE
  messages <- character(0)

  if (sigma_between <= 0) {
    valid <- FALSE
    messages <- c(messages, "sigma_between must be positive")
  }

  if (sigma_within <= 0) {
    valid <- FALSE
    messages <- c(messages, "sigma_within must be positive")
  }

  if (n_periods < 4) {
    valid <- FALSE
    messages <- c(messages, "n_periods should be at least 4 for N-of-1 design")
  }

  # Calculate variance ratio (ICC)
  if (valid) {
    total_var <- sigma_between^2 + sigma_within^2
    icc <- sigma_between^2 / total_var

    if (icc > 0.9) {
      messages <- c(messages,
        sprintf("Warning: Very high ICC (%.2f) may limit N-of-1 advantage",
                icc))
    }
  }

  list(
    valid = valid,
    messages = messages,
    icc = if (valid) icc else NA,
    total_variance = if (valid) total_var else NA
  )
}

#==========================================================================
# SIMULATION PARAMETERS
#==========================================================================

# Simulation counts (adjusted for quick_test mode)
if (quick_test) {
  n_simulations <- 50
  n_power_sims <- 20
  n_null_sims <- 200
  n_strategy_sims <- 50
  sensitivity_n_sims <- 50
} else {
  n_simulations <- 1000
  n_power_sims <- 70
  n_null_sims <- 500
  n_strategy_sims <- 200
  sensitivity_n_sims <- 200
}

true_effect <- -2  # Prazosin reduces nightmares by 2 points on average
baseline_nightmares <- 8  # Average nightmares per week at baseline
individual_variance <- 1.5  # Between-individual variance in baseline
measurement_error <- 1.0  # Within-individual measurement error

# Carryover parameters (pharmacokinetic modeling)
period_length_days <- 7  # Length of each period in days
carryover_halflife_days <- 3  # Half-life of carryover effect in days

#' Calculate carryover effect using multiple decay models
#'
#' Implements three decay models for sensitivity analysis:
#' - Exponential: Based on pharmacokinetic half-life (default)
#' - Linear: Simple linear decay to zero
#' - Weibull: Flexible decay with shape parameter
#'
#' @param time_since_treatment Time since treatment discontinuation (days)
#' @param treatment_effect Maximum treatment effect at discontinuation
#' @param model Decay model: "exponential", "linear", or "weibull"
#' @param halflife_days Half-life for exponential decay (days)
#' @param total_decay_time Total time for linear decay to reach zero (days)
#' @param weibull_scale Scale parameter (lambda) for Weibull decay
#' @param weibull_shape Shape parameter (k) for Weibull decay (k=1 is exponential)
#' @return Residual carryover effect
#'
#' @details
#' Decay model equations:
#' \itemize{
#'   \item Exponential: effect * exp(-ln(2) * t / halflife)
#'   \item Linear: effect * max(0, 1 - t / total_time)
#'   \item Weibull: effect * exp(-(t / lambda)^k)
#' }
#'
#' @examples
#' # Exponential decay with 3-day half-life
#' calculate_carryover_effect(7, -2, "exponential", halflife_days = 3)
#'
#' # Linear decay over 14 days
#' calculate_carryover_effect(7, -2, "linear", total_decay_time = 14)
#'
#' # Weibull decay with accelerating pattern (k > 1)
#' calculate_carryover_effect(7, -2, "weibull", weibull_scale = 5, weibull_shape = 2)
calculate_carryover_effect <- function(time_since_treatment,
                                       treatment_effect,
                                       model = "exponential",
                                       halflife_days = 3,
                                       total_decay_time = 14,
                                       weibull_scale = 5,
                                       weibull_shape = 1) {

  # Input validation

  if (!is.finite(treatment_effect) || time_since_treatment < 0) {
    return(0)
  }

  decay_fraction <- switch(model,

    "exponential" = {
      if (halflife_days <= 0) return(0)
      exp(-log(2) * time_since_treatment / halflife_days)
    },

    "linear" = {
      if (total_decay_time <= 0) return(0)
      max(0, 1 - time_since_treatment / total_decay_time)
    },

    "weibull" = {
      if (weibull_scale <= 0 || weibull_shape <= 0) return(0)
      exp(-(time_since_treatment / weibull_scale)^weibull_shape)
    },

    # Default to exponential if unknown model
    {
      warning(sprintf("Unknown decay model '%s', using exponential", model))
      if (halflife_days <= 0) return(0)
      exp(-log(2) * time_since_treatment / halflife_days)
    }
  )

  treatment_effect * decay_fraction
}

#' Legacy wrapper for backward compatibility
#'
#' @param halflife_days Half-life in days
#' @param period_days Period length in days
#' @param treatment_effect Maximum treatment effect
#' @return Residual effect after period (using exponential decay)
calculate_carryover_effect_legacy <- function(halflife_days, period_days,
                                              treatment_effect) {
  calculate_carryover_effect(
    time_since_treatment = period_days,
    treatment_effect = treatment_effect,
    model = "exponential",
    halflife_days = halflife_days
  )
}

#' Compare carryover decay across models
#'
#' @param time_points Vector of time points to evaluate
#' @param treatment_effect Treatment effect magnitude
#' @param halflife_days Half-life for exponential model
#' @param total_decay_time Total decay time for linear model
#' @param weibull_scale Scale for Weibull model
#' @param weibull_shape Shape for Weibull model
#' @return Data frame with decay values for each model
compare_decay_models <- function(time_points = seq(0, 14, by = 1),
                                 treatment_effect = -2,
                                 halflife_days = 3,
                                 total_decay_time = 14,
                                 weibull_scale = 5,
                                 weibull_shape = 1.5) {

  data.frame(
    time = time_points,
    exponential = sapply(time_points, function(t) {
      calculate_carryover_effect(t, treatment_effect, "exponential",
                                 halflife_days = halflife_days)
    }),
    linear = sapply(time_points, function(t) {
      calculate_carryover_effect(t, treatment_effect, "linear",
                                 total_decay_time = total_decay_time)
    }),
    weibull = sapply(time_points, function(t) {
      calculate_carryover_effect(t, treatment_effect, "weibull",
                                 weibull_scale = weibull_scale,
                                 weibull_shape = weibull_shape)
    })
  )
}

#' Create standard simulation result structure  
#'
#' @param ... Named arguments for result fields
#' @return Standardized list structure
create_sim_result <- function(...) {
  args <- list(...)
  required_fields <- c("p_value", "effect_estimate", "significant")
  
  # Ensure required fields exist
  for (field in required_fields) {
    if (!field %in% names(args)) {
      args[[field]] <- switch(field,
                             p_value = 1,
                             effect_estimate = 0, 
                             significant = FALSE)
    }
  }
  
  # Ensure significant field is logical
  args$significant <- as.logical(args$significant)
  args
}

# Calculate current carryover effect (using default exponential decay)
carryover_effect <- calculate_carryover_effect(
  time_since_treatment = period_length_days,
  treatment_effect = true_effect,
  model = "exponential",
  halflife_days = carryover_halflife_days
)

# Store default decay model for sensitivity analysis
default_decay_model <- "exponential"

# RCT parameters
# Primary comparison uses per-participant parity at N = 35; the
# N = 70 cell is a sensitivity row showing what the RCT achieves
# at twice the participant count (1:2 patient ratio), matching the
# Hendrickson 2020 reference grid.
rct_sample_size <- 35  # Per-participant parity vs. N-of-1
rct_sample_size_sensitivity <- 70  # 1:2 ratio sensitivity row

# N-of-1 parameters
n_of_1_periods <- 8  # 4 treatment, 4 placebo periods per individual
n_of_1_individuals <- 35  # Per-participant parity with parallel RCT

# Statistical constants
ALPHA_LEVEL <- 0.05  # Significance level
DF_APPROX <- 6  # Conservative estimate for degrees of freedom in N-of-1 trials
Z_CRITICAL <- 1.96  # Critical value for 95% confidence intervals

# Plot colors
PLOT_COLORS <- list(
  rct = "#2E86AB",
  n_of_1 = "#A23B72"
)

cat("=== SIMULATION SETUP ===\n")
cat("True treatment effect:", true_effect, "nightmares/week\n")
cat("Period length:", period_length_days, "days\n")
cat("Carryover half-life:", carryover_halflife_days, "days\n")
cat("Calculated carryover effect:", round(carryover_effect, 3), 
    "nightmares/week\n")
cat("RCT sample size:", rct_sample_size, "total participants\n")
cat("N-of-1 design:", n_of_1_individuals, "individuals ×", n_of_1_periods, 
    "periods each\n\n")

#==========================================================================
# FUNCTION: SIMULATE RCT
#==========================================================================

#' Simulate a randomized controlled trial
#'
#' @return Standardized simulation result with p_value, effect_estimate, significant
#' Draw a shared pool of individual patient latent baselines.
#'
#' Morris, White & Crowther (2019) §5.4: paired comparisons between
#' competing methods should operate on the SAME simulated population
#' within each replicate. This helper produces that shared population.
#' @param n Integer. Size of the patient pool.
#' @return Numeric vector of length `n` of per-patient baseline
#'   nightmare frequencies.
simulate_latent_population <- function(n) {
  stats::rnorm(n, baseline_nightmares, individual_variance)
}

simulate_rct <- function(latent_baselines = NULL) {
  # Morris §5.4: accept an optional shared latent pool so that paired
  # comparisons with N-of-1 operate on the SAME patients within a
  # replicate. When `latent_baselines` is supplied, the first
  # `rct_sample_size` values are used.
  individual_baselines <- if (is.null(latent_baselines)) {
    rnorm(rct_sample_size, baseline_nightmares,
          individual_variance)
  } else {
    stopifnot(length(latent_baselines) >= rct_sample_size)
    latent_baselines[seq_len(rct_sample_size)]
  }

  # Randomly assign to treatment (1) or placebo (0)
  treatment <- rbinom(rct_sample_size, 1, 0.5)

  # Generate outcomes
  outcomes <- individual_baselines +
              treatment * true_effect +
              rnorm(rct_sample_size, 0, measurement_error)
  
  # Perform t-test
  t_test <- t.test(outcomes[treatment == 1], outcomes[treatment == 0])
  
  create_sim_result(
    p_value = t_test$p.value,
    effect_estimate = mean(outcomes[treatment == 1]) - 
                     mean(outcomes[treatment == 0]),
    significant = t_test$p.value < ALPHA_LEVEL
  )
}

#==========================================================================
# FUNCTION: SIMULATE SINGLE N-OF-1 TRIAL WITH CARRY-OVER EFFECTS
#==========================================================================

#' Simulate a single N-of-1 trial with carry-over effects
#'
#' @return Standardized simulation result with additional periods_used field
simulate_n_of_1_individual <- function(latent_baseline = NULL) {
  # Morris §5.4 paired design: when `latent_baseline` is supplied, it
  # is the patient's true baseline drawn from the shared pool instead
  # of a fresh draw.
  individual_baseline <- if (is.null(latent_baseline)) {
    stats::rnorm(1, baseline_nightmares, individual_variance)
  } else {
    latent_baseline
  }
  
  # Randomize period order (ABAB vs BABA design)
  # A = placebo, B = treatment
  periods <- sample(c(rep("A", n_of_1_periods/2), 
                     rep("B", n_of_1_periods/2)))
  
  # Generate outcomes for each period with carry-over effects
  outcomes <- numeric(n_of_1_periods)
  for (i in 1:n_of_1_periods) {
    # Direct treatment effect
    treatment_effect <- ifelse(periods[i] == "B", true_effect, 0)
    
    # Carry-over effect: residual benefit from previous treatment period
    carryover <- 0
    if (i > 1 && periods[i-1] == "B" && periods[i] == "A") {
      # Previous period was treatment, current is placebo
      carryover <- carryover_effect
    }
    
    outcomes[i] <- individual_baseline + treatment_effect + carryover + 
                   rnorm(1, 0, measurement_error)
  }
  
  # Analyze with crossover design accounting for carry-over
  treatment_outcomes <- outcomes[periods == "B"]
  placebo_outcomes <- outcomes[periods == "A"]
  
  # Identify placebo periods that follow treatment (carry-over periods)
  carryover_periods <- numeric(0)
  for (i in 2:n_of_1_periods) {
    if (periods[i-1] == "B" && periods[i] == "A") {
      carryover_periods <- c(carryover_periods, i)
    }
  }
  
  # Separate placebo outcomes into pure placebo vs carry-over affected
  if (length(carryover_periods) > 0) {
    pure_placebo_outcomes <- placebo_outcomes[
      !((1:length(placebo_outcomes)) %in% 
        match(carryover_periods, which(periods == "A")))]
  } else {
    pure_placebo_outcomes <- placebo_outcomes
  }
  
  if (length(treatment_outcomes) > 0 & length(pure_placebo_outcomes) > 0) {
    # Compare treatment to pure placebo periods only
    mean_diff <- mean(treatment_outcomes) - mean(pure_placebo_outcomes)
    se_diff <- sqrt(var(treatment_outcomes)/length(treatment_outcomes) + 
                   var(pure_placebo_outcomes)/length(pure_placebo_outcomes))
    t_stat <- mean_diff / se_diff
    df <- length(treatment_outcomes) + length(pure_placebo_outcomes) - 2
    p_val <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)
    
    create_sim_result(
      p_value = p_val,
      effect_estimate = mean_diff,
      significant = p_val < ALPHA_LEVEL,
      periods_used = list(
        treatment = periods == "B", 
        pure_placebo = setdiff(which(periods == "A"), carryover_periods))
    )
  } else {
    create_sim_result(
      p_value = 1, 
      effect_estimate = 0, 
      significant = FALSE,
      periods_used = list(treatment = c(), pure_placebo = c())
    )
  }
}

#==========================================================================
# FUNCTION: SIMULATE MULTIPLE N-OF-1 TRIALS
#==========================================================================

#' Simulate multiple N-of-1 trials and perform meta-analysis
#'
#' Conducts multiple N-of-1 trials and combines results using 
#' DerSimonian-Laird random-effects meta-analysis, providing population-level
#' treatment effect estimates with heterogeneity assessment.
#'
#' @return A list containing:
#' \describe{
#'   \item{p_value}{Combined p-value from random-effects meta-analysis}
#'   \item{effect_estimate}{Pooled treatment effect estimate}
#'   \item{se}{Standard error of pooled estimate}
#'   \item{ci_lower}{Lower bound of 95% confidence interval}
#'   \item{ci_upper}{Upper bound of 95% confidence interval}
#'   \item{tau_squared}{Between-study variance (heterogeneity)}
#'   \item{i_squared}{I² statistic (percentage of heterogeneity)}
#'   \item{significant}{Logical indicating statistical significance}
#'   \item{individual_p_values}{Vector of individual trial p-values}
#'   \item{individual_effects}{Vector of individual effect estimates}
#'   \item{method}{Meta-analysis method used}
#' }
#'
#' @details
#' Meta-analysis process:
#' \itemize{
#'   \item Estimates standard errors from individual trial statistics
#'   \item Weights studies by inverse variance (fixed-effects)
#'   \item Calculates between-study heterogeneity (tau-squared, I²)
#'   \item Applies random-effects weights incorporating heterogeneity
#'   \item Provides pooled effect estimate with confidence intervals
#'   \item Handles missing/invalid results with appropriate defaults
#' }
#' 
#' Uses global parameter: n_of_1_individuals. Calls simulate_n_of_1_individual()
#' for each individual trial.
#'
#' @references
#' DerSimonian, R., & Laird, N. (1986). Meta-analysis in clinical trials. 
#' Controlled Clinical Trials, 7(3), 177-188.
#' 
#' Zucker, D. R., Ruthazer, R., & Schmid, C. H. (2010). Individual (N-of-1) 
#' trials can be combined to give population comparative treatment effect 
#' estimates: methodologic considerations. Journal of Clinical Epidemiology, 
#' 63(12), 1312-1323.
#'
#' @examples
#' \dontrun{
#' # Requires global parameters and simulate_n_of_1_individual() function
#' meta_result <- simulate_n_of_1_study()
#' print(paste("Pooled effect:", meta_result$effect_estimate))
#' print(paste("Heterogeneity I²:", meta_result$i_squared, "%"))
#' }
simulate_n_of_1_study <- function(latent_baselines = NULL) {
  # Morris §5.4 paired design: when `latent_baselines` is supplied,
  # the first `n_of_1_individuals` entries are used as each N-of-1
  # patient's true baseline, overlapping with the patients used by
  # `simulate_rct()` on the same replicate.
  if (is.null(latent_baselines)) {
    individual_results <- replicate(
      n_of_1_individuals,
      simulate_n_of_1_individual(),
      simplify = FALSE
    )
  } else {
    stopifnot(length(latent_baselines) >= n_of_1_individuals)
    individual_results <- lapply(
      seq_len(n_of_1_individuals),
      function(k) simulate_n_of_1_individual(latent_baselines[k])
    )
  }
  
  # Extract effect estimates and standard errors
  effect_estimates <- sapply(individual_results, 
                             function(x) x$effect_estimate)
  p_values <- sapply(individual_results, function(x) x$p_value)
  
  # Calculate standard errors from t-tests (approximate)
  std_errors <- numeric(n_of_1_individuals)
  for (i in 1:n_of_1_individuals) {
    if (!is.na(p_values[i]) && !is.na(effect_estimates[i]) &&
        p_values[i] < 0.999 && p_values[i] > 0.001) {
      # Approximate SE from p-value and effect estimate
      t_stat <- qt(p_values[i]/2, DF_APPROX, lower.tail = FALSE)
      if (t_stat > 0 && !is.na(t_stat) && effect_estimates[i] != 0) {
        std_errors[i] <- abs(effect_estimates[i]) / t_stat
      } else {
        std_errors[i] <- 1  # Default for failed calculations
      }
    } else {
      std_errors[i] <- 1  # Default for extreme p-values or NA effects
    }
  }
  
  # Handle missing or invalid effect estimates
  valid_indices <- !is.na(effect_estimates) & !is.na(std_errors) & 
                   std_errors > 0
  if (sum(valid_indices) == 0) {
    # Return default values if no valid estimates
    return(list(
      p_value = 1,
      effect_estimate = 0,
      se = 1,
      ci_lower = -1.96,
      ci_upper = 1.96,
      tau_squared = 0,
      i_squared = 0,
      significant = FALSE,
      individual_p_values = p_values,
      individual_effects = effect_estimates,
      method = "DerSimonian-Laird random-effects"
    ))
  }
  
  # Use only valid estimates
  valid_effects <- effect_estimates[valid_indices]
  valid_se <- std_errors[valid_indices]
  k_valid <- length(valid_effects)
  
  # Random-effects meta-analysis (DerSimonian & Laird method)
  # Weight by inverse variance
  weights <- 1 / (valid_se^2)
  weighted_effect <- sum(weights * valid_effects) / sum(weights)
  
  # Calculate heterogeneity
  Q <- sum(weights * (valid_effects - weighted_effect)^2)
  df_Q <- k_valid - 1
  
  # Ensure Q and df_Q are positive for heterogeneity calculation
  if (df_Q > 0 && Q > 0) {
    tau_squared <- max(0, (Q - df_Q) / 
                       (sum(weights) - sum(weights^2)/sum(weights)))
    i_squared <- max(0, (Q - df_Q) / Q * 100)
  } else {
    tau_squared <- 0
    i_squared <- 0
  }
  
  # Random-effects weights
  re_weights <- 1 / (valid_se^2 + tau_squared)
  re_effect <- sum(re_weights * valid_effects) / sum(re_weights)
  re_se <- sqrt(1 / sum(re_weights))
  
  # Test statistic and p-value
  z_stat <- re_effect / re_se
  combined_p_value <- 2 * (1 - pnorm(abs(z_stat)))
  
  # 95% confidence interval
  ci_lower <- re_effect - Z_CRITICAL * re_se
  ci_upper <- re_effect + Z_CRITICAL * re_se
  
  return(list(
    p_value = combined_p_value,
    effect_estimate = re_effect,
    se = re_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    tau_squared = tau_squared,
    i_squared = i_squared,
    significant = combined_p_value < 0.05,
    individual_p_values = p_values,
    individual_effects = effect_estimates,
    method = "DerSimonian-Laird random-effects"
  ))
}

#==========================================================================
# MAIN SIMULATION EXECUTION
#==========================================================================

cat("Running", n_simulations, "simulations...\n")

#' Run one paired replicate: draw a shared patient pool and pass it to
#' both designs so the RCT and the N-of-1 study operate on overlapping
#' people (Morris, White & Crowther 2019 §5.4, Question B pairing).
run_paired_replicate <- function() {
  pool_size <- max(rct_sample_size, n_of_1_individuals)
  latents <- simulate_latent_population(pool_size)
  list(
    rct  = simulate_rct(latent_baselines = latents),
    nof1 = simulate_n_of_1_study(latent_baselines = latents)
  )
}

paired_results <- replicate(n_simulations, run_paired_replicate(),
                            simplify = FALSE)
rct_results    <- lapply(paired_results, `[[`, "rct")
n_of_1_results <- lapply(paired_results, `[[`, "nof1")

# Paired-difference MCSE for rejection rates (McNemar-style): the
# variance of the paired proportion difference accounts for the
# within-replicate correlation introduced by the shared latent pool.
rct_reject  <- vapply(rct_results,    function(x) x$significant, logical(1))
nof1_reject <- vapply(n_of_1_results, function(x) x$significant, logical(1))
delta_reject     <- rct_reject - nof1_reject
paired_power_diff      <- mean(delta_reject)
paired_power_diff_mcse <- sqrt(stats::var(delta_reject) / n_simulations)
cat(sprintf(
  "Paired power difference (RCT - N-of-1): %.3f (MCSE = %.4f)\n",
  paired_power_diff, paired_power_diff_mcse
))

# Extract results
rct_power <- mean(sapply(rct_results, function(x) x$significant))
rct_effects <- sapply(rct_results, function(x) x$effect_estimate)

n_of_1_power <- mean(sapply(n_of_1_results, function(x) x$significant))
n_of_1_effects <- sapply(n_of_1_results, function(x) x$effect_estimate)

#==========================================================================
# RESULTS SUMMARY
#==========================================================================

cat("\n=== SIMULATION RESULTS ===\n")
cat("True effect size:", true_effect, "nightmares per week\n")
cat("Note: Negative values indicate benefit (nightmare reduction)\n\n")

cat("Randomized Controlled Trial (n =", rct_sample_size, "):\n")
cat("  Statistical Power:", round(rct_power * 100, 1), "%\n")
cat("  Mean Effect Estimate:", round(mean(rct_effects), 2), "\n")
cat("  Effect Estimate SD:", round(sd(rct_effects), 2), "\n")
cat("  Theoretical SE:", 
    round(sqrt(2 * (individual_variance^2 + measurement_error^2) / 
              (rct_sample_size/2)), 2), "\n\n")

cat("N-of-1 Study (", n_of_1_individuals, "individuals,", n_of_1_periods, 
    "periods each):\n")
cat("  Statistical Power:", round(n_of_1_power * 100, 1), "%\n")
cat("  Mean Effect Estimate:", round(mean(n_of_1_effects), 2), "\n")
cat("  Effect Estimate SD:", round(sd(n_of_1_effects), 2), "\n")
cat("  Meta-analysis method: DerSimonian-Laird random-effects\n")

# Calculate means safely
i_squared_values <- sapply(n_of_1_results, function(x) {
  if (is.null(x$i_squared) || is.na(x$i_squared)) return(0)
  return(as.numeric(x$i_squared))
})
tau_squared_values <- sapply(n_of_1_results, function(x) {
  if (is.null(x$tau_squared) || is.na(x$tau_squared)) return(0)
  return(as.numeric(x$tau_squared))
})

cat("  Mean I² (heterogeneity):", 
    round(mean(i_squared_values, na.rm = TRUE), 1), "%\n")
cat("  Mean τ² (between-study variance):", 
    round(mean(tau_squared_values, na.rm = TRUE), 3), "\n\n")

#==========================================================================
# VISUALIZATION: COMPARISON PLOTS
#==========================================================================

# Create comparison plots
results_df <- data.frame(
  Study_Type = rep(c("RCT", "N-of-1"), each = n_simulations),
  Effect_Estimate = c(rct_effects, n_of_1_effects),
  Significant = c(sapply(rct_results, function(x) x$significant),
                  sapply(n_of_1_results, function(x) x$significant))
)

# Plot 1: Distribution of effect estimates
par(mfrow = c(2, 1))
for (st in c("RCT", "N-of-1")) {
  d <- results_df$Effect_Estimate[results_df$Study_Type == st]
  hist(d, breaks = 30, col = adjustcolor("steelblue", 0.5),
       main = paste("Distribution of Effect Estimates:", st),
       xlab = "Effect Estimate (reduction in nightmares/week)",
       ylab = "Frequency")
  abline(v = true_effect, lty = 2, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))

# Plot 2: Power comparison across effect sizes and carryover half-lives
effect_sizes <- seq(0, -3, by = -0.5)
halflife_days_range <- c(0.5, 1, 2, 3, 5, 7)
# n_power_sims defined in SIMULATION PARAMETERS section

cat("Calculating power across effect sizes and carryover half-lives...\n")

# Pre-calculate carryover lookup table for efficiency
carryover_lookup <- expand.grid(
  Effect_Size = effect_sizes,
  Halflife_Days = halflife_days_range
)
carryover_lookup$Calculated_Carryover <- mapply(
  function(hl, eff) {
    calculate_carryover_effect(
      time_since_treatment = period_length_days,
      treatment_effect = eff,
      model = "exponential",
      halflife_days = hl
    )
  },
  hl = carryover_lookup$Halflife_Days,
  eff = carryover_lookup$Effect_Size
)

# Pre-allocate results data frame. Power_MCSE is computed via the
# existing `calculate_mcse_proportion()` helper per Morris, White,
# and Crowther (2019) Table 6 and is now reported alongside Power
# for every cell of the scenario grid (closing audit gap 3).
n_conditions <- nrow(carryover_lookup)
power_results <- data.frame(
  Effect_Size = rep(carryover_lookup$Effect_Size, each = 2),
  Halflife_Days = rep(carryover_lookup$Halflife_Days, each = 2),
  Calculated_Carryover = rep(carryover_lookup$Calculated_Carryover, each = 2),
  Study_Type = rep(c("RCT", "N-of-1"), n_conditions),
  Power = numeric(n_conditions * 2),
  Power_MCSE = numeric(n_conditions * 2),
  N_Sims = rep(0L, n_conditions * 2)
)

# Helper function for global variable management. Returns power and
# its Monte Carlo SE for both designs, plus the paired-difference
# MCSE. Morris, White & Crowther (2019) §5.4 paired design: the RCT
# and N-of-1 studies within each replicate share a single latent
# patient pool via `run_paired_replicate()`.
run_power_simulation <- function(effect_val, carryover_val, n_sims) {
  # Save current global values
  saved_effect <- true_effect
  saved_carryover <- carryover_effect

  # Set new values
  true_effect <<- effect_val
  carryover_effect <<- carryover_val

  # Run paired replicates: within each replicate, RCT and N-of-1
  # operate on the SAME simulated patient pool.
  paired <- replicate(n_sims, run_paired_replicate(),
                      simplify = FALSE)
  rct_sims <- lapply(paired, `[[`, "rct")
  n_of_1_sims <- lapply(paired, `[[`, "nof1")

  # Restore original values
  true_effect <<- saved_effect
  carryover_effect <<- saved_carryover

  # Return power values + MCSE per Morris Table 6
  rct_power <- mean(sapply(rct_sims, function(x) x$significant))
  nof1_power <- mean(sapply(n_of_1_sims, function(x) x$significant))
  c(
    RCT = rct_power,
    RCT_MCSE = calculate_mcse_proportion(rct_power, n_sims),
    N_of_1 = nof1_power,
    N_of_1_MCSE = calculate_mcse_proportion(nof1_power, n_sims),
    N_Sims = n_sims
  )
}

# Run simulations
for (i in 1:n_conditions) {
  power_vals <- run_power_simulation(
    carryover_lookup$Effect_Size[i],
    carryover_lookup$Calculated_Carryover[i],
    n_power_sims
  )

  # Fill in results
  idx_rct <- (i - 1) * 2 + 1
  idx_n_of_1 <- (i - 1) * 2 + 2
  power_results$Power[idx_rct] <- power_vals["RCT"]
  power_results$Power_MCSE[idx_rct] <- power_vals["RCT_MCSE"]
  power_results$N_Sims[idx_rct] <- as.integer(power_vals["N_Sims"])
  power_results$Power[idx_n_of_1] <- power_vals["N_of_1"]
  power_results$Power_MCSE[idx_n_of_1] <- power_vals["N_of_1_MCSE"]
  power_results$N_Sims[idx_n_of_1] <- as.integer(power_vals["N_Sims"])

  if (i %% 5 == 0) {
    cat(sprintf("  Progress: %d/%d conditions completed\n", i, n_conditions))
  }
}

# No need to restore values - helper function handles this automatically

# Create half-life labels for plotting
power_results$Halflife_Label <- paste0("Half-life: ", 
                                       power_results$Halflife_Days, " days")

# Plot 2A: Power curves by half-life (faceted)
unique_hl <- sort(unique(power_results$Halflife_Days))
n_panels <- length(unique_hl)
par(mfrow = c(2, ceiling(n_panels / 2)))
for (hl_val in unique_hl) {
  d <- power_results[power_results$Halflife_Days == hl_val, ]
  plot(NULL, xlim = c(0, 3), ylim = c(0, 100),
       xlab = "Effect Size (Absolute)", ylab = "Power (%)",
       main = paste0("Half-life: ", hl_val, " days"))
  abline(h = 5, lty = 3, col = "gray50")
  abline(h = 80, lty = 2, col = "gray50")
  for (st in c("RCT", "N-of-1")) {
    dd <- d[d$Study_Type == st, ]
    col_val <- if (st == "RCT") PLOT_COLORS$rct else PLOT_COLORS$n_of_1
    lty_val <- if (st == "RCT") 1 else 2
    lines(abs(dd$Effect_Size), dd$Power * 100, col = col_val, lwd = 2, lty = lty_val)
    points(abs(dd$Effect_Size), dd$Power * 100, col = col_val, pch = 16, cex = 0.8)
  }
  legend("topleft", legend = c("RCT", "N-of-1"),
         col = c(PLOT_COLORS$rct, PLOT_COLORS$n_of_1),
         lwd = 2, lty = c(1, 2), cex = 0.7)
}
par(mfrow = c(1, 1))

# Plot 2B: Power difference heatmap (N-of-1 minus RCT)
power_wide <- reshape(power_results[, c("Effect_Size", "Halflife_Days", 
                                        "Study_Type", "Power")],
                     idvar = c("Effect_Size", "Halflife_Days"),
                     timevar = "Study_Type", direction = "wide")
names(power_wide) <- c("Effect_Size", "Halflife_Days", 
                      "RCT_Power", "N_of_1_Power")
power_wide$Power_Difference <- (power_wide$N_of_1_Power - 
                               power_wide$RCT_Power) * 100

# Plot 2B: Power difference table (text-based, replacing heatmap)
cat("\nPower Advantage: N-of-1 vs RCT (percentage points)\n")
cat("Effect Size | Halflife | Difference\n")
cat("------------|----------|----------\n")
for (i in seq_len(nrow(power_wide))) {
  cat(sprintf("    %4.1f    |  %4.1f d  |  %+5.1f%%\n",
              abs(power_wide$Effect_Size[i]),
              power_wide$Halflife_Days[i],
              power_wide$Power_Difference[i]))
}

# Print summary statistics
cat("\nCarryover Effects by Half-life (for effect size -2.0):\n")
cat("Half-life (days) | Carryover Effect | % of Original Effect\n")
cat("-----------------|------------------|---------------------\n")
reference_effect <- -2.0
for (halflife_val in halflife_days_range) {
  carryover_val <- calculate_carryover_effect(
    time_since_treatment = period_length_days,
    treatment_effect = reference_effect,
    model = "exponential",
    halflife_days = halflife_val
  )
  percent_original <- abs(carryover_val / reference_effect) * 100
  cat(sprintf("      %4.1f      |      %6.3f      |       %5.1f%%\n", 
              halflife_val, carryover_val, percent_original))
}

# Summary of optimal conditions
cat("\nOptimal Design Conditions:\n")
max_n_of_1_advantage <- power_wide[which.max(power_wide$Power_Difference), ]
max_rct_advantage <- power_wide[which.min(power_wide$Power_Difference), ]

cat(sprintf("N-of-1 most advantageous: Effect size %4.1f, Half-life %3.1f days (N-of-1 +%4.1f%%)\n",
            abs(max_n_of_1_advantage$Effect_Size), 
            max_n_of_1_advantage$Halflife_Days,
            max_n_of_1_advantage$Power_Difference))

cat(sprintf("RCT most advantageous: Effect size %4.1f, Half-life %3.1f days (RCT +%4.1f%%)\n",
            abs(max_rct_advantage$Effect_Size), 
            max_rct_advantage$Halflife_Days,
            abs(max_rct_advantage$Power_Difference)))

#==========================================================================
# SAMPLE N-OF-1 TRIAL VISUALIZATION WITH CARRY-OVER
#==========================================================================

cat("\n=== SAMPLE N-of-1 TRIAL (with carry-over effects) ===\n")
sample_individual <- rnorm(1, baseline_nightmares, individual_variance)
sample_periods <- sample(c(rep("Placebo", n_of_1_periods/2), 
                          rep("Prazosin", n_of_1_periods/2)))
sample_outcomes <- numeric(n_of_1_periods)

for (i in 1:n_of_1_periods) {
  # Direct treatment effect
  treatment_effect <- ifelse(sample_periods[i] == "Prazosin", 
                             true_effect, 0)
  
  # Carry-over effect
  carryover <- 0
  if (i > 1 && sample_periods[i-1] == "Prazosin" && 
      sample_periods[i] == "Placebo") {
    carryover <- carryover_effect
  }
  
  sample_outcomes[i] <- sample_individual + treatment_effect + carryover + 
                       rnorm(1, 0, measurement_error)
}

# Identify carry-over affected periods
carryover_affected <- logical(n_of_1_periods)
for (i in 2:n_of_1_periods) {
  if (sample_periods[i-1] == "Prazosin" && sample_periods[i] == "Placebo") {
    carryover_affected[i] <- TRUE
  }
}

# Create enhanced visualization
sample_df <- data.frame(
  Period = 1:n_of_1_periods,
  Treatment = sample_periods,
  Nightmares = sample_outcomes,
  Carryover = carryover_affected,
  Period_Type = ifelse(carryover_affected, "Placebo (carry-over)", 
                      ifelse(sample_periods == "Prazosin", "Prazosin", 
                            "Placebo (pure)"))
)

plot(sample_df$Period, sample_df$Nightmares, type = "l", col = "gray50",
     xlab = "Time Period", ylab = "Nightmares per Week",
     main = "Sample N-of-1 Trial with Carry-over Effects\nOrange = carry-over affected")
pt_cols <- ifelse(sample_df$Period_Type == "Prazosin", "blue",
            ifelse(sample_df$Period_Type == "Placebo (carry-over)", "orange", "red"))
points(sample_df$Period, sample_df$Nightmares, pch = 16, cex = 1.5, col = pt_cols)
legend("topright",
       legend = c("Placebo (pure)", "Prazosin", "Placebo (carry-over)"),
       col = c("red", "blue", "orange"), pch = 16)

#==========================================================================
# SAMPLE N-OF-1 ANALYSIS OUTPUT
#==========================================================================

# Print sample data with carry-over information
print(sample_df)

cat("\nSample N-of-1 Analysis (accounting for carry-over):\n")
prazosin_periods <- sample_outcomes[sample_periods == "Prazosin"]
pure_placebo_periods <- sample_outcomes[sample_periods == "Placebo" & 
                                       !carryover_affected]
carryover_periods <- sample_outcomes[carryover_affected]

cat("  Prazosin periods mean:", round(mean(prazosin_periods), 2), "\n")
cat("  Pure placebo periods mean:", 
    round(mean(pure_placebo_periods), 2), "\n") 
if(length(carryover_periods) > 0) {
  cat("  Carry-over affected periods mean:", 
      round(mean(carryover_periods), 2), "\n")
}
cat("  Estimated treatment effect (vs pure placebo):", 
    round(mean(prazosin_periods) - mean(pure_placebo_periods), 2), "\n")
cat("  Carry-over effect estimate:",
    ifelse(length(carryover_periods) > 0,
           round(mean(carryover_periods) - mean(pure_placebo_periods), 2),
           "None detected"), "\n")

#==========================================================================
# MORRIS ET AL. (2019) COMPLIANT REPORTING
#==========================================================================

cat("\n=== MORRIS ET AL. (2019) COMPLIANT SUMMARY ===\n")

# Calculate comprehensive performance measures
rct_perf <- summarize_simulation_performance(
  estimates = rct_effects,
  true_value = true_effect,
  p_values = sapply(rct_results, function(x) x$p_value)
)

n_of_1_perf <- summarize_simulation_performance(
  estimates = n_of_1_effects,
  true_value = true_effect,
  p_values = sapply(n_of_1_results, function(x) x$p_value)
)

cat("\nRCT Performance Measures (with Monte Carlo SEs):\n")
cat(sprintf("  Bias: %.4f (MCSE = %.4f)\n",
            rct_perf$bias, rct_perf$bias_mcse))
cat(sprintf("  Empirical SE: %.4f (MCSE = %.4f)\n",
            rct_perf$empirical_se, rct_perf$emp_se_mcse))
cat(sprintf("  Power: %.1f%% (MCSE = %.1f%%, 95%% CI: %.1f%%-%.1f%%)\n",
            rct_perf$power * 100, rct_perf$power_mcse * 100,
            rct_perf$power_ci_lower * 100, rct_perf$power_ci_upper * 100))

cat("\nN-of-1 Performance Measures (with Monte Carlo SEs):\n")
cat(sprintf("  Bias: %.4f (MCSE = %.4f)\n",
            n_of_1_perf$bias, n_of_1_perf$bias_mcse))
cat(sprintf("  Empirical SE: %.4f (MCSE = %.4f)\n",
            n_of_1_perf$empirical_se, n_of_1_perf$emp_se_mcse))
cat(sprintf("  Power: %.1f%% (MCSE = %.1f%%, 95%% CI: %.1f%%-%.1f%%)\n",
            n_of_1_perf$power * 100, n_of_1_perf$power_mcse * 100,
            n_of_1_perf$power_ci_lower * 100, n_of_1_perf$power_ci_upper * 100))

# nsim justification
cat("\nSimulation Size Justification:\n")
expected_powers <- c(0.65, 0.70, 0.75, 0.80)
target_mcses <- c(0.02, 0.015, 0.01, 0.005)
cat("  Required nsim for target MCSE (assuming power ~ 0.75):\n")
for (mcse in target_mcses) {
  req_nsim <- calculate_required_nsim(0.75, mcse)
  cat(sprintf("    MCSE = %.1f%% requires nsim >= %d\n", mcse * 100, req_nsim))
}
cat(sprintf("  Current nsim = %d provides MCSE = %.2f%% for power = %.1f%%\n",
            n_simulations,
            calculate_mcse_proportion(n_of_1_perf$power, n_simulations) * 100,
            n_of_1_perf$power * 100))

#==========================================================================
# TYPE I ERROR EVALUATION (NULL HYPOTHESIS)
#==========================================================================

cat("\n=== TYPE I ERROR EVALUATION ===\n")
cat("Running simulations under null hypothesis (true effect = 0)...\n")

# Save original values
saved_true_effect <- true_effect
saved_carryover_effect <- carryover_effect

# Set null condition
true_effect <- 0
carryover_effect <- 0
# n_null_sims defined in SIMULATION PARAMETERS section

# Run null simulations under Morris §5.4 paired design: both designs
# share the same latent patient pool within each replicate.
set.seed(123)  # Different seed for Type I error evaluation
null_paired <- replicate(n_null_sims, run_paired_replicate(),
                         simplify = FALSE)
rct_null_results    <- lapply(null_paired, `[[`, "rct")
n_of_1_null_results <- lapply(null_paired, `[[`, "nof1")

# Calculate Type I error rates
rct_type1 <- mean(sapply(rct_null_results, function(x) x$significant))
n_of_1_type1 <- mean(sapply(n_of_1_null_results, function(x) x$significant))

rct_type1_mcse <- calculate_mcse_proportion(rct_type1, n_null_sims)
n_of_1_type1_mcse <- calculate_mcse_proportion(n_of_1_type1, n_null_sims)

# Paired-difference MCSE for Type I error (shared-latent variance).
rct_null_reject  <- vapply(rct_null_results,
                           function(x) x$significant, logical(1))
nof1_null_reject <- vapply(n_of_1_null_results,
                           function(x) x$significant, logical(1))
delta_null_reject    <- rct_null_reject - nof1_null_reject
paired_type1_diff    <- mean(delta_null_reject)
paired_type1_diff_mcse <- sqrt(
  stats::var(delta_null_reject) / n_null_sims
)
cat(sprintf(
  "Paired Type I difference (RCT - N-of-1): %.3f (MCSE = %.4f)\n",
  paired_type1_diff, paired_type1_diff_mcse
))

cat(sprintf("\nRCT Type I Error: %.1f%% (MCSE = %.1f%%, 95%% CI: %.1f%%-%.1f%%)\n",
            rct_type1 * 100, rct_type1_mcse * 100,
            (rct_type1 - 1.96 * rct_type1_mcse) * 100,
            (rct_type1 + 1.96 * rct_type1_mcse) * 100))
cat(sprintf("N-of-1 Type I Error: %.1f%% (MCSE = %.1f%%, 95%% CI: %.1f%%-%.1f%%)\n",
            n_of_1_type1 * 100, n_of_1_type1_mcse * 100,
            (n_of_1_type1 - 1.96 * n_of_1_type1_mcse) * 100,
            (n_of_1_type1 + 1.96 * n_of_1_type1_mcse) * 100))

# Check if Type I error is within acceptable range
nominal_alpha <- 0.05
acceptable_range <- c(0.025, 0.075)  # Reasonable range for nsim=500

cat("\nType I Error Control Assessment:\n")
if (rct_type1 >= acceptable_range[1] && rct_type1 <= acceptable_range[2]) {
  cat("  RCT: ACCEPTABLE (within 2.5%-7.5% range)\n")
} else {
  cat(sprintf("  RCT: OUTSIDE ACCEPTABLE RANGE (%.1f%%)\n", rct_type1 * 100))
}

if (n_of_1_type1 >= acceptable_range[1] && n_of_1_type1 <= acceptable_range[2]) {
  cat("  N-of-1: ACCEPTABLE (within 2.5%-7.5% range)\n")
} else {
  cat(sprintf("  N-of-1: OUTSIDE ACCEPTABLE RANGE (%.1f%%)\n",
              n_of_1_type1 * 100))
}

# Restore original values
true_effect <- saved_true_effect
carryover_effect <- saved_carryover_effect

#==========================================================================
# EXPLICIT CARRYOVER MODELING (ALTERNATIVE TO EXCLUSION)
#==========================================================================

#' Simulate N-of-1 with explicit carryover modeling in analysis
#'
#' Instead of excluding carryover-affected periods, this function
#' includes a carryover term in the analysis model, following
#' Granholm et al. (2021) recommendations.
#'
#' @return Standardized simulation result
simulate_n_of_1_explicit_carryover <- function() {
  # Individual's true baseline
  individual_baseline <- rnorm(1, baseline_nightmares, individual_variance)

  # Randomize period order
  periods <- sample(c(rep("A", n_of_1_periods/2),
                     rep("B", n_of_1_periods/2)))

  # Generate outcomes with carryover
  outcomes <- numeric(n_of_1_periods)
  carryover_term <- numeric(n_of_1_periods)

  for (i in 1:n_of_1_periods) {
    treatment_effect_i <- ifelse(periods[i] == "B", true_effect, 0)

    # Calculate carryover from previous periods
    if (i > 1 && periods[i-1] == "B" && periods[i] == "A") {
      carryover_term[i] <- carryover_effect
    }

    outcomes[i] <- individual_baseline + treatment_effect_i +
                   carryover_term[i] + rnorm(1, 0, measurement_error)
  }

  # Analysis with EXPLICIT carryover term (not exclusion)
  treatment_indicator <- as.numeric(periods == "B")

  # Simple regression with carryover term
  if (sum(carryover_term != 0) > 0) {
    # Fit model: outcome ~ treatment + carryover
    fit <- lm(outcomes ~ treatment_indicator + carryover_term)
    coef_summary <- summary(fit)$coefficients

    if (nrow(coef_summary) >= 2) {
      treatment_est <- coef_summary["treatment_indicator", "Estimate"]
      treatment_se <- coef_summary["treatment_indicator", "Std. Error"]
      treatment_t <- coef_summary["treatment_indicator", "t value"]
      treatment_p <- coef_summary["treatment_indicator", "Pr(>|t|)"]

      return(create_sim_result(
        p_value = treatment_p,
        effect_estimate = treatment_est,
        significant = treatment_p < ALPHA_LEVEL,
        method = "explicit_carryover"
      ))
    }
  }

  # Fallback: simple comparison without carryover term
  treatment_outcomes <- outcomes[periods == "B"]
  placebo_outcomes <- outcomes[periods == "A"]

  mean_diff <- mean(treatment_outcomes) - mean(placebo_outcomes)
  se_diff <- sqrt(var(treatment_outcomes)/length(treatment_outcomes) +
                 var(placebo_outcomes)/length(placebo_outcomes))
  t_stat <- mean_diff / se_diff
  df <- length(treatment_outcomes) + length(placebo_outcomes) - 2
  p_val <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)

  create_sim_result(
    p_value = p_val,
    effect_estimate = mean_diff,
    significant = p_val < ALPHA_LEVEL,
    method = "simple_comparison"
  )
}

# Compare exclusion vs explicit modeling strategies
cat("\n=== CARRYOVER MODELING STRATEGY COMPARISON ===\n")
cat("Comparing period exclusion vs explicit carryover modeling...\n")
# n_strategy_sims defined in SIMULATION PARAMETERS section

set.seed(42)

# Exclusion strategy (current approach)
exclusion_results <- replicate(n_strategy_sims, simulate_n_of_1_study(),
                               simplify = FALSE)
exclusion_power <- mean(sapply(exclusion_results, function(x) x$significant))
exclusion_effects <- sapply(exclusion_results, function(x) x$effect_estimate)

# Explicit modeling strategy
explicit_results <- replicate(n_strategy_sims,
                              simulate_n_of_1_explicit_carryover(),
                              simplify = FALSE)
explicit_power <- mean(sapply(explicit_results, function(x) x$significant))
explicit_effects <- sapply(explicit_results, function(x) x$effect_estimate)

cat(sprintf("\nExclusion Strategy:\n"))
cat(sprintf("  Power: %.1f%% (MCSE = %.1f%%)\n",
            exclusion_power * 100,
            calculate_mcse_proportion(exclusion_power, n_strategy_sims) * 100))
cat(sprintf("  Bias: %.3f\n", mean(exclusion_effects) - true_effect))

cat(sprintf("\nExplicit Carryover Modeling:\n"))
cat(sprintf("  Power: %.1f%% (MCSE = %.1f%%)\n",
            explicit_power * 100,
            calculate_mcse_proportion(explicit_power, n_strategy_sims) * 100))
cat(sprintf("  Bias: %.3f\n", mean(explicit_effects) - true_effect))

power_diff <- (explicit_power - exclusion_power) * 100
cat(sprintf("\nPower difference (explicit - exclusion): %+.1f percentage points\n",
            power_diff))
if (power_diff > 0) {
  cat("Recommendation: Explicit carryover modeling may improve power\n")
} else {
  cat("Note: Exclusion strategy performs comparably in this scenario\n")
}

#==========================================================================
# CARRYOVER DECAY MODEL SENSITIVITY ANALYSIS
#==========================================================================

cat("\n=== CARRYOVER DECAY MODEL SENSITIVITY ANALYSIS ===\n")
cat("Comparing exponential, linear, and Weibull decay models...\n\n")

# Parameters for sensitivity analysis
# sensitivity_n_sims defined in SIMULATION PARAMETERS section
sensitivity_effect_size <- -2.0

# Decay model parameters calibrated to approximately match at day 7
# (end of washout period)
decay_params <- list(
  exponential = list(halflife_days = 3),
  linear = list(total_decay_time = 14),
  weibull = list(weibull_scale = 5, weibull_shape = 1.5)
)

# 1. Visualize decay curves
cat("1. Decay Curve Comparison\n")
decay_comparison <- compare_decay_models(
  time_points = seq(0, 14, by = 0.5),
  treatment_effect = sensitivity_effect_size,
  halflife_days = decay_params$exponential$halflife_days,
  total_decay_time = decay_params$linear$total_decay_time,
  weibull_scale = decay_params$weibull$weibull_scale,
  weibull_shape = decay_params$weibull$weibull_shape
)

# Print key time points
cat("\n   Residual carryover effect by model (treatment effect = -2.0):\n")
cat("   Time (days) | Exponential |   Linear   |  Weibull   |\n")
cat("   ------------|-------------|------------|------------|\n")
key_times <- c(0, 3, 7, 10, 14)
for (t in key_times) {
  idx <- which.min(abs(decay_comparison$time - t))
  cat(sprintf("       %2d      |    %6.3f   |   %6.3f   |   %6.3f   |\n",
              t,
              decay_comparison$exponential[idx],
              decay_comparison$linear[idx],
              decay_comparison$weibull[idx]))
}

# Plot decay curves
decay_long <- tidyr::pivot_longer(
  decay_comparison,
  cols = c("exponential", "linear", "weibull"),
  names_to = "model",
  values_to = "effect"
)
decay_long$model <- factor(decay_long$model,
                           levels = c("exponential", "linear", "weibull"),
                           labels = c("Exponential", "Linear", "Weibull"))

decay_models <- levels(decay_long$model)
decay_cols <- c("#2E86AB", "#A23B72", "#E8871E")
plot(NULL, xlim = range(decay_long$time), ylim = range(decay_long$effect),
     xlab = "Days Since Treatment Discontinuation",
     ylab = "Residual Carryover Effect",
     main = sprintf("Carryover Decay Models Comparison\nTreatment effect = %.1f, Period length = %d days",
                    sensitivity_effect_size, period_length_days))
abline(h = 0, lty = 2, col = "gray50")
abline(v = period_length_days, lty = 3, col = "gray30")
for (j in seq_along(decay_models)) {
  dd <- decay_long[decay_long$model == decay_models[j], ]
  lines(dd$time, dd$effect, col = decay_cols[j], lwd = 2)
}
legend("bottomright", legend = decay_models, col = decay_cols, lwd = 2)

# 2. Run simulations under each decay model
cat("\n2. Power Sensitivity to Decay Model\n")

#' Simulate N-of-1 study with specified decay model
#'
#' @param decay_model One of "exponential", "linear", "weibull"
#' @param params List of model-specific parameters
#' @return Simulation result
simulate_n_of_1_with_decay_model <- function(decay_model, params) {
  # Calculate carryover effect for this model
  local_carryover <- calculate_carryover_effect(
    time_since_treatment = period_length_days,
    treatment_effect = true_effect,
    model = decay_model,
    halflife_days = params$halflife_days %||% 3,
    total_decay_time = params$total_decay_time %||% 14,
    weibull_scale = params$weibull_scale %||% 5,
    weibull_shape = params$weibull_shape %||% 1
  )

  # Individual's true baseline
  individual_baseline <- rnorm(1, baseline_nightmares, individual_variance)

  # Randomize period order
  periods <- sample(c(rep("A", n_of_1_periods/2),
                     rep("B", n_of_1_periods/2)))

  # Generate outcomes
  outcomes <- numeric(n_of_1_periods)
  for (i in 1:n_of_1_periods) {
    treatment_effect_i <- ifelse(periods[i] == "B", true_effect, 0)
    carryover <- 0
    if (i > 1 && periods[i-1] == "B" && periods[i] == "A") {
      carryover <- local_carryover
    }
    outcomes[i] <- individual_baseline + treatment_effect_i + carryover +
                   rnorm(1, 0, measurement_error)
  }

  # Analysis (excluding carryover periods)
  treatment_outcomes <- outcomes[periods == "B"]
  placebo_outcomes <- outcomes[periods == "A"]

  carryover_periods <- numeric(0)
  for (i in 2:n_of_1_periods) {
    if (periods[i-1] == "B" && periods[i] == "A") {
      carryover_periods <- c(carryover_periods, i)
    }
  }

  if (length(carryover_periods) > 0) {
    pure_placebo_outcomes <- placebo_outcomes[
      !((1:length(placebo_outcomes)) %in%
        match(carryover_periods, which(periods == "A")))]
  } else {
    pure_placebo_outcomes <- placebo_outcomes
  }

  if (length(treatment_outcomes) > 0 & length(pure_placebo_outcomes) > 0) {
    mean_diff <- mean(treatment_outcomes) - mean(pure_placebo_outcomes)
    se_diff <- sqrt(var(treatment_outcomes)/length(treatment_outcomes) +
                   var(pure_placebo_outcomes)/length(pure_placebo_outcomes))
    t_stat <- mean_diff / se_diff
    df <- length(treatment_outcomes) + length(pure_placebo_outcomes) - 2
    p_val <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)

    create_sim_result(
      p_value = p_val,
      effect_estimate = mean_diff,
      significant = p_val < ALPHA_LEVEL,
      carryover_value = local_carryover,
      decay_model = decay_model
    )
  } else {
    create_sim_result(
      p_value = 1,
      effect_estimate = 0,
      significant = FALSE,
      carryover_value = local_carryover,
      decay_model = decay_model
    )
  }
}

# Null coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x

# Run simulations for each decay model (parallelized)
set.seed(42)
decay_model_names <- c("exponential", "linear", "weibull")
decay_results <- list()

cat(sprintf("   Using %d cores for parallel processing\n", n_cores))

#' Run single N-of-1 study simulation for sensitivity analysis
#'
#' @param sim_id Simulation replicate ID (for RNG)
#' @param model_name Decay model name
#' @param params Model parameters
#' @param n_individuals Number of individuals per study
#' @param alpha Significance level
#' @return List with significance, effect estimate, and p-value
run_single_sensitivity_sim <- function(sim_id, model_name, params,
                                       n_individuals, alpha) {
  # Run individual N-of-1 simulations
  individual_results <- replicate(n_individuals,
                                 simulate_n_of_1_with_decay_model(model_name,
                                                                 params),
                                 simplify = FALSE)

  # Extract effects and p-values
  effects <- sapply(individual_results, function(x) x$effect_estimate)
  p_vals <- sapply(individual_results, function(x) x$p_value)

  # Simplified meta-analysis
  valid_idx <- !is.na(effects) & !is.na(p_vals)
  if (sum(valid_idx) < 2) {
    return(list(significant = FALSE, effect_estimate = 0, p_value = 1))
  }

  pooled_effect <- mean(effects[valid_idx])
  pooled_se <- sd(effects[valid_idx]) / sqrt(sum(valid_idx))
  z_stat <- pooled_effect / pooled_se
  combined_p <- 2 * (1 - pnorm(abs(z_stat)))

  list(
    significant = combined_p < alpha,
    effect_estimate = pooled_effect,
    p_value = combined_p
  )
}

for (model_name in decay_model_names) {
  cat(sprintf("   Running %s model simulations (parallel)...\n", model_name))
  params <- decay_params[[model_name]]

  # Set unique seeds for each simulation to ensure reproducibility
  sim_seeds <- 42 + (1:sensitivity_n_sims) +
               which(decay_model_names == model_name) * 1000

  # Parallel simulation using mclapply (Unix) or parLapply (Windows)
  if (.Platform$OS.type == "unix") {
    model_results <- mclapply(
      1:sensitivity_n_sims,
      function(i) {
        set.seed(sim_seeds[i])
        run_single_sensitivity_sim(i, model_name, params,
                                   n_of_1_individuals, ALPHA_LEVEL)
      },
      mc.cores = n_cores,
      mc.set.seed = FALSE
    )
  } else {
    # Windows: use parLapply with explicit cluster
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("simulate_n_of_1_with_decay_model",
                       "calculate_carryover_effect", "create_sim_result",
                       "period_length_days", "true_effect",
                       "baseline_nightmares", "individual_variance",
                       "n_of_1_periods", "measurement_error",
                       "ALPHA_LEVEL", "%||%"))
    model_results <- parLapply(
      cl,
      1:sensitivity_n_sims,
      function(i) {
        set.seed(sim_seeds[i])
        run_single_sensitivity_sim(i, model_name, params,
                                   n_of_1_individuals, ALPHA_LEVEL)
      }
    )
    stopCluster(cl)
  }

  decay_results[[model_name]] <- model_results
}

# 3. Summarize results
cat("\n3. Summary Statistics by Decay Model\n\n")

sensitivity_summary <- data.frame(
  Model = character(),
  Carryover_at_Washout = numeric(),
  Power_Pct = numeric(),
  Power_MCSE = numeric(),
  Mean_Effect = numeric(),
  Bias = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in decay_model_names) {
  results <- decay_results[[model_name]]
  power_val <- mean(sapply(results, function(x) x$significant))
  effects <- sapply(results, function(x) x$effect_estimate)

  # Calculate carryover at end of washout
  params <- decay_params[[model_name]]
  carryover_at_washout <- calculate_carryover_effect(
    time_since_treatment = period_length_days,
    treatment_effect = sensitivity_effect_size,
    model = model_name,
    halflife_days = params$halflife_days %||% 3,
    total_decay_time = params$total_decay_time %||% 14,
    weibull_scale = params$weibull_scale %||% 5,
    weibull_shape = params$weibull_shape %||% 1
  )

  sensitivity_summary <- rbind(sensitivity_summary, data.frame(
    Model = tools::toTitleCase(model_name),
    Carryover_at_Washout = carryover_at_washout,
    Power_Pct = power_val * 100,
    Power_MCSE = calculate_mcse_proportion(power_val, sensitivity_n_sims) * 100,
    Mean_Effect = mean(effects),
    Bias = mean(effects) - sensitivity_effect_size
  ))
}

cat("   Model       | Carryover | Power (%) | MCSE (%) | Mean Est |  Bias  |\n")
cat("   ------------|-----------|-----------|----------|----------|--------|\n")
for (i in 1:nrow(sensitivity_summary)) {
  cat(sprintf("   %-11s |   %6.3f  |   %5.1f   |   %4.1f   |  %6.3f  | %+6.3f |\n",
              sensitivity_summary$Model[i],
              sensitivity_summary$Carryover_at_Washout[i],
              sensitivity_summary$Power_Pct[i],
              sensitivity_summary$Power_MCSE[i],
              sensitivity_summary$Mean_Effect[i],
              sensitivity_summary$Bias[i]))
}

# 4. Statistical comparison
cat("\n4. Model Comparison Analysis\n")

power_range <- max(sensitivity_summary$Power_Pct) -
               min(sensitivity_summary$Power_Pct)
cat(sprintf("   Power range across models: %.1f percentage points\n", power_range))

if (power_range < 5) {
  cat("   Interpretation: Results are ROBUST to decay model choice.\n")
  cat("   Power estimates vary by < 5 percentage points across models.\n")
} else if (power_range < 10) {
  cat("   Interpretation: MODERATE sensitivity to decay model choice.\n")
  cat("   Consider reporting results under multiple decay assumptions.\n")
} else {
  cat("   Interpretation: HIGH sensitivity to decay model choice.\n")
  cat("   Decay model selection substantially affects conclusions.\n")
}

# 5. Create visualization of power by model
power_df <- data.frame(
  Model = factor(sensitivity_summary$Model,
                levels = c("Exponential", "Linear", "Weibull")),
  Power = sensitivity_summary$Power_Pct,
  MCSE = sensitivity_summary$Power_MCSE
)

bar_cols <- c("#2E86AB", "#A23B72", "#E8871E")
bp <- barplot(power_df$Power, names.arg = as.character(power_df$Model),
              col = bar_cols, ylim = c(0, 100),
              main = sprintf("Power Sensitivity to Carryover Decay Model\nN-of-1: %d individuals x %d periods, %d sims",
                             n_of_1_individuals, n_of_1_periods, sensitivity_n_sims),
              xlab = "Decay Model", ylab = "Statistical Power (%)")
arrows(bp, power_df$Power - 1.96 * power_df$MCSE,
       bp, power_df$Power + 1.96 * power_df$MCSE,
       angle = 90, code = 3, length = 0.1)
abline(h = 80, lty = 2, col = "gray50")

# Store sensitivity analysis results for manuscript
decay_sensitivity_results <- list(
  summary = sensitivity_summary,
  decay_curves = decay_comparison,
  params = decay_params,
  n_sims = sensitivity_n_sims,
  power_range = power_range,
  robust = power_range < 5
)

#==========================================================================
# CARRYOVER MAGNITUDE SENSITIVITY ANALYSIS
#==========================================================================

cat("\n=== CARRYOVER MAGNITUDE SENSITIVITY ANALYSIS ===\n")
cat("Assessing N-of-1 power robustness across carryover levels...\n\n")

# Test carryover half-lives from minimal to severe
carryover_halflife_levels <- c(0.5, 1, 2, 3, 5, 7, 10, 14)
carryover_n_sims <- sensitivity_n_sims

# Use smaller effect size to avoid ceiling effects
carryover_test_effect <- -0.75  # Smaller effect to see power differences
cat(sprintf("Testing with effect size = %.2f (smaller to avoid ceiling)\n\n",
            carryover_test_effect))

# Calculate residual carryover at each level
cat("Carryover levels being tested:\n")
cat("Half-life | Residual Effect | % of Treatment |\n")
cat("----------|-----------------|----------------|\n")
for (hl in carryover_halflife_levels) {
  residual <- calculate_carryover_effect(
    time_since_treatment = period_length_days,
    treatment_effect = carryover_test_effect,
    model = "exponential",
    halflife_days = hl
  )
  pct <- abs(residual / carryover_test_effect) * 100
  cat(sprintf("  %4.1f d  |      %6.3f     |     %5.1f%%     |\n",
              hl, residual, pct))
}

# Function to run N-of-1 simulation with specific carryover
simulate_n_of_1_with_carryover <- function(carryover_hl, test_effect) {
  local_carryover <- calculate_carryover_effect(
    time_since_treatment = period_length_days,
    treatment_effect = test_effect,
    model = "exponential",
    halflife_days = carryover_hl
  )

  # Simulate multiple individuals
  individual_results <- lapply(1:n_of_1_individuals, function(ind) {
    individual_baseline <- rnorm(1, baseline_nightmares, individual_variance)
    periods <- sample(c(rep("A", n_of_1_periods/2),
                       rep("B", n_of_1_periods/2)))

    outcomes <- numeric(n_of_1_periods)
    for (i in 1:n_of_1_periods) {
      treatment_effect_i <- ifelse(periods[i] == "B", test_effect, 0)
      carryover <- 0
      if (i > 1 && periods[i-1] == "B" && periods[i] == "A") {
        carryover <- local_carryover
      }
      outcomes[i] <- individual_baseline + treatment_effect_i + carryover +
                     rnorm(1, 0, measurement_error)
    }

    # Analysis excluding carryover periods
    treatment_outcomes <- outcomes[periods == "B"]
    placebo_indices <- which(periods == "A")
    carryover_indices <- numeric(0)
    for (i in 2:n_of_1_periods) {
      if (periods[i-1] == "B" && periods[i] == "A") {
        carryover_indices <- c(carryover_indices, i)
      }
    }
    pure_placebo_indices <- setdiff(placebo_indices, carryover_indices)
    pure_placebo_outcomes <- outcomes[pure_placebo_indices]

    if (length(treatment_outcomes) > 1 && length(pure_placebo_outcomes) > 1) {
      mean_diff <- mean(treatment_outcomes) - mean(pure_placebo_outcomes)
      se_diff <- sqrt(var(treatment_outcomes)/length(treatment_outcomes) +
                     var(pure_placebo_outcomes)/length(pure_placebo_outcomes))
      list(effect = mean_diff, se = se_diff)
    } else {
      list(effect = NA, se = NA)
    }
  })

  # Meta-analysis
  effects <- sapply(individual_results, function(x) x$effect)
  ses <- sapply(individual_results, function(x) x$se)
  valid <- !is.na(effects) & !is.na(ses) & ses > 0

  if (sum(valid) < 2) {
    return(list(significant = FALSE, effect = 0, p_value = 1))
  }

  effects <- effects[valid]
  ses <- ses[valid]

  # DerSimonian-Laird
  weights <- 1 / ses^2
  pooled_effect <- sum(weights * effects) / sum(weights)
  Q <- sum(weights * (effects - pooled_effect)^2)
  k <- length(effects)
  tau2 <- max(0, (Q - (k - 1)) / (sum(weights) - sum(weights^2)/sum(weights)))
  re_weights <- 1 / (ses^2 + tau2)
  re_effect <- sum(re_weights * effects) / sum(re_weights)
  re_se <- sqrt(1 / sum(re_weights))
  z_stat <- re_effect / re_se
  p_val <- 2 * (1 - pnorm(abs(z_stat)))

  list(significant = p_val < ALPHA_LEVEL, effect = re_effect, p_value = p_val)
}

# Run simulations across carryover levels (parallelized)
cat("\nRunning simulations across carryover levels...\n")
cat(sprintf("Using %d cores for parallel processing\n", n_cores))

carryover_power_results <- data.frame(
  Halflife_Days = numeric(),
  Residual_Carryover = numeric(),
  Pct_of_Treatment = numeric(),
  N_of_1_Power = numeric(),
  Power_MCSE = numeric(),
  Mean_Effect = numeric(),
  Bias = numeric(),
  stringsAsFactors = FALSE
)

for (hl in carryover_halflife_levels) {
  cat(sprintf("   Half-life = %.1f days...\n", hl))

  # Set seeds for reproducibility
  sim_seeds <- 42 + (1:carryover_n_sims) + hl * 1000

  if (.Platform$OS.type == "unix") {
    results <- mclapply(
      1:carryover_n_sims,
      function(i) {
        set.seed(sim_seeds[i])
        simulate_n_of_1_with_carryover(hl, carryover_test_effect)
      },
      mc.cores = n_cores,
      mc.set.seed = FALSE
    )
  } else {
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("simulate_n_of_1_with_carryover",
                       "calculate_carryover_effect", "period_length_days",
                       "carryover_test_effect", "n_of_1_individuals",
                       "n_of_1_periods", "baseline_nightmares",
                       "individual_variance", "measurement_error",
                       "ALPHA_LEVEL"))
    results <- parLapply(cl, 1:carryover_n_sims, function(i) {
      set.seed(sim_seeds[i])
      simulate_n_of_1_with_carryover(hl, carryover_test_effect)
    })
    stopCluster(cl)
  }

  power_val <- mean(sapply(results, function(x) x$significant))
  effects <- sapply(results, function(x) x$effect)
  residual <- calculate_carryover_effect(period_length_days, carryover_test_effect,
                                         "exponential", halflife_days = hl)

  carryover_power_results <- rbind(carryover_power_results, data.frame(
    Halflife_Days = hl,
    Residual_Carryover = residual,
    Pct_of_Treatment = abs(residual / carryover_test_effect) * 100,
    N_of_1_Power = power_val * 100,
    Power_MCSE = calculate_mcse_proportion(power_val, carryover_n_sims) * 100,
    Mean_Effect = mean(effects),
    Bias = mean(effects) - carryover_test_effect
  ))
}

# Add RCT power for comparison (unaffected by carryover)
# Need to temporarily set true_effect to carryover_test_effect
cat("   Running RCT baseline for comparison...\n")
saved_effect_for_rct <- true_effect
true_effect <- carryover_test_effect
set.seed(42)
rct_baseline_results <- replicate(carryover_n_sims, simulate_rct(),
                                  simplify = FALSE)
rct_baseline_power <- mean(sapply(rct_baseline_results, function(x) x$significant))
rct_baseline_mcse <- calculate_mcse_proportion(rct_baseline_power, carryover_n_sims)
true_effect <- saved_effect_for_rct  # Restore

# Display results
cat("\n=== N-of-1 Power Across Carryover Levels ===\n\n")
cat("Half-life | Residual | % Trt | N-of-1 Power | MCSE  |  Bias  | vs RCT |\n")
cat("----------|----------|-------|--------------|-------|--------|--------|\n")
for (i in 1:nrow(carryover_power_results)) {
  r <- carryover_power_results[i, ]
  power_diff <- r$N_of_1_Power - (rct_baseline_power * 100)
  cat(sprintf("  %4.1f d  |  %6.3f  | %5.1f | %10.1f%%  | %4.1f%% | %+6.3f | %+5.1f%% |\n",
              r$Halflife_Days, r$Residual_Carryover, r$Pct_of_Treatment,
              r$N_of_1_Power, r$Power_MCSE, r$Bias, power_diff))
}
cat(sprintf("\nRCT baseline power: %.1f%% (MCSE = %.1f%%)\n",
            rct_baseline_power * 100, rct_baseline_mcse * 100))

# Plot power vs carryover
plot(carryover_power_results$Halflife_Days, carryover_power_results$N_of_1_Power,
     type = "b", col = "#A23B72", lwd = 2, pch = 16,
     xlim = range(carryover_halflife_levels), ylim = c(0, 100),
     xlab = "Carryover Half-life (days)", ylab = "Statistical Power (%)",
     main = sprintf("N-of-1 Power Robustness to Carryover Effect\nN = %d, %d periods, true effect = %.1f",
                    n_of_1_individuals, n_of_1_periods, true_effect))
arrows(carryover_power_results$Halflife_Days,
       carryover_power_results$N_of_1_Power - 1.96 * carryover_power_results$Power_MCSE,
       carryover_power_results$Halflife_Days,
       carryover_power_results$N_of_1_Power + 1.96 * carryover_power_results$Power_MCSE,
       angle = 90, code = 3, length = 0.05, col = "#A23B72")
abline(h = rct_baseline_power * 100, lty = 2, col = "#2E86AB", lwd = 2)
abline(h = 80, lty = 3, col = "gray50")
legend("topright", legend = c("N-of-1", "RCT Power"),
       col = c("#A23B72", "#2E86AB"), lwd = 2, lty = c(1, 2))

# Summary interpretation
power_range_carryover <- max(carryover_power_results$N_of_1_Power) -
                         min(carryover_power_results$N_of_1_Power)
min_power_hl <- carryover_power_results$Halflife_Days[
  which.min(carryover_power_results$N_of_1_Power)]

cat(sprintf("\nSummary:\n"))
cat(sprintf("  Power range across carryover levels: %.1f percentage points\n",
            power_range_carryover))
cat(sprintf("  Minimum N-of-1 power: %.1f%% at half-life = %.1f days\n",
            min(carryover_power_results$N_of_1_Power), min_power_hl))

if (min(carryover_power_results$N_of_1_Power) > rct_baseline_power * 100) {
  cat("  Conclusion: N-of-1 maintains power ADVANTAGE over RCT across all carryover levels\n")
} else {
  crossover_hl <- carryover_power_results$Halflife_Days[
    which(carryover_power_results$N_of_1_Power < rct_baseline_power * 100)[1]]
  cat(sprintf("  Conclusion: N-of-1 loses advantage to RCT when half-life >= %.1f days\n",
              crossover_hl))
}

# Store results
carryover_sensitivity_results <- list(
  power_results = carryover_power_results,
  rct_power = rct_baseline_power,
  rct_mcse = rct_baseline_mcse,
  power_range = power_range_carryover,
  n_sims = carryover_n_sims
)

cat("\n=== SIMULATION COMPLETE ===\n")
