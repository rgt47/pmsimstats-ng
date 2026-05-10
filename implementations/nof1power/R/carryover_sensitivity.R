# carryover_sensitivity_refactored.R - Sensitivity analysis for carryover effects

#' Run carryover sensitivity analysis
#'
#' This function tests how sensitive trial design comparisons are to
#' different carryover model assumptions. It runs multiple simulations
#' using different carryover models and parameter settings, then compares
#' the results across models and scenarios.
#'
#' @param trial_designs List of trial designs to compare
#' @param base_model_params Base model parameters (without carryover model specification)
#' @param resp_param_sets Response parameter sets
#' @param bl_param_sets Baseline parameter sets
#' @param sim_options Simulation options (reduced n_reps recommended for sensitivity testing)
#' @param analysis_options Analysis options
#' @param carryover_models Vector of carryover models to test
#' @param carryover_scenarios List of parameter scenarios for each model
#' @param log_file Optional file path to write progress logs
#' @return List containing sensitivity analysis results
#' @export
run_carryover_sensitivity <- function(
  trial_designs,
  base_model_params,
  resp_param_sets,
  bl_param_sets,
  sim_options = list(n_reps = 25),  # Reduced for sensitivity testing
  analysis_options = list(use_expectancy = TRUE, random_slope = FALSE),
  carryover_models = c("exponential", "linear", "weibull", "power"),
  carryover_scenarios = NULL,
  log_file = NULL
) {

  # Define default scenarios if not provided
  if (is.null(carryover_scenarios)) {
    carryover_scenarios <- list(
      exponential = list(
        list(carryover_t1half = 0.5, scalefactor = 2),
        list(carryover_t1half = 1.0, scalefactor = 2),
        list(carryover_t1half = 2.0, scalefactor = 2)
      ),
      linear = list(
        list(carryover_t1half = 1.0, total_time = 2),
        list(carryover_t1half = 1.0, total_time = 4),
        list(carryover_t1half = 1.0, total_time = 6)
      ),
      weibull = list(
        list(carryover_t1half = 1.0, k = 0.8),  # Slower initial decay
        list(carryover_t1half = 1.0, k = 1.0),  # Exponential (equivalent)
        list(carryover_t1half = 1.0, k = 1.5)   # Faster initial decay
      ),
      power = list(
        list(carryover_t1half = 1.0, power = 1.0),  # Linear-like
        list(carryover_t1half = 1.0, power = 1.5),  # Moderate curvature
        list(carryover_t1half = 1.0, power = 2.0)   # Quadratic
      )
    )
  }

  # Create expanded parameter grid
  sensitivity_results <- list()

  # Progress tracking
  total_scenarios <- sum(sapply(carryover_scenarios, length))
  current_scenario <- 0

  log_progress("Running carryover sensitivity analysis...\n", log_file = log_file)
  log_progress("Testing ", length(carryover_models), " models with ", total_scenarios, " total scenarios\n\n", log_file = log_file)

  # Create a progress bar that can log to file
  pb <- logging_progress_bar(min = 0, max = total_scenarios, style = 3, log_file = log_file)

  for (model in carryover_models) {
    model_results <- list()
    scenarios <- carryover_scenarios[[model]]

    log_progress("Testing ", model, " model with ", length(scenarios), " scenarios...\n", log_file = log_file)

    for (scenario_idx in seq_along(scenarios)) {
      current_scenario <- current_scenario + 1
      scenario_params <- scenarios[[scenario_idx]]

      log_progress("  Scenario ", scenario_idx, " of ", length(scenarios),
          " (", current_scenario, " of ", total_scenarios, ")...", log_file = log_file)

      # Create model parameters for this scenario. The new
      # build_sigma_matrix / apply_carryover_to_component pipeline
      # reads the decay form from `carryover_form`, not from
      # `carryover_model`; translate on assignment.
      model_params <- base_model_params
      model_params$carryover_form <- model

      # Add scenario-specific parameters, translating shape/weibull
      # keys to the canonical names where needed.
      for (param_name in names(scenario_params)) {
        canon <- switch(param_name,
          "k" = "weibull_shape",
          "scalefactor" = "scale_factor",
          param_name
        )
        model_params[[canon]] <- scenario_params[[param_name]]
      }

      # Run simulation
      start_time <- Sys.time()

      # Add log file to sim_options if provided
      if (!is.null(log_file)) {
        sim_options$log_file <- log_file
      }

      # Set up parallel processing with correct function globals if possible
      if (requireNamespace("future", quietly = TRUE) &&
          requireNamespace("furrr", quietly = TRUE)) {

        if (exists("setup_parallel_processing", mode = "function")) {
          # Use the optimized parallel setup function
          setup_parallel_processing(
            additional_functions = c(
              "generate_data",
              "generate_simulated_results",
              "lme_analysis",
              "calculate_carryover",
              "apply_carryover_effects",
              "mod_gompertz",
              "cumulative",
              "log_progress"
            )
          )
        } else {
          # Manual setup as fallback
          future_globals <- list(
            generate_data = if(exists("generate_data")) generate_data else NULL,
            generate_simulated_results = if(exists("generate_simulated_results")) generate_simulated_results else NULL,
            lme_analysis = if(exists("lme_analysis")) lme_analysis else NULL,
            calculate_carryover = if(exists("calculate_carryover")) calculate_carryover else NULL,
            apply_carryover_effects = if(exists("apply_carryover_effects")) apply_carryover_effects else NULL,
            mod_gompertz = if(exists("mod_gompertz")) mod_gompertz else NULL,
            cumulative = if(exists("cumulative")) cumulative else NULL
          )

          # Clean up NULL entries
          future_globals <- future_globals[!sapply(future_globals, is.null)]

          # Set up basic multisession plan
          n_workers <- min(parallel::detectCores() - 1, 8)
          future::plan(future::multisession, workers = n_workers)

          # Increase globals size limit
          options(future.globals.maxSize = 1000 * 1024^2)  # 1GB

          # Make globals directly available in global environment
          # This allows worker processes to find them automatically
          for (func_name in names(future_globals)) {
            if (!is.null(future_globals[[func_name]])) {
              assign(func_name, future_globals[[func_name]], envir = .GlobalEnv)
            }
          }
        }
      }

      # Adapt legacy sim_options / analysis_options to the
      # pmsimstats-ng tidyverse generate_simulated_results signature.
      simparam_new <- list(
        Nreps = sim_options$n_reps %||% 25,
        progressiveSave = isTRUE(sim_options$progressive_save),
        basesavename = sim_options$base_save_name %||% "carryover",
        nRep2save = sim_options$n_rep_to_save %||% .Machine$integer.max,
        saveunit2start = sim_options$save_unit_to_start %||% 1,
        savedir = sim_options$save_dir %||% "."
      )
      analysisparams_new <- if (is.data.frame(analysis_options)) {
        analysis_options
      } else {
        do.call(data.frame, c(analysis_options,
                              stringsAsFactors = FALSE))
      }

      dgp_arch <- sim_options$dgp_architecture %||% "mvn"

      results <- generate_simulated_results(
        trialdesigns = trial_designs,
        respparamsets = resp_param_sets,
        blparamsets = bl_param_sets,
        censorparams = NA,
        modelparams = model_params,
        simparam = simparam_new,
        analysisparams = analysisparams_new,
        rawdataout = FALSE,
        dgp_architecture = dgp_arch
      )

      end_time <- Sys.time()
      log_progress(" completed in ", round(difftime(end_time, start_time, units = "secs"), 1), " seconds\n", log_file = log_file)

      # Update progress bar
      update_logging_progress_bar(pb, current_scenario)

      # Store results with metadata
      model_results[[paste0("scenario_", scenario_idx)]] <- list(
        results = results$results,
        model = model,
        parameters = scenario_params,
        scenario_id = scenario_idx
      )
    }

    sensitivity_results[[model]] <- model_results
  }

  # Close progress bar
  close_logging_progress_bar(pb)

  log_progress("\nSensitivity analysis complete!\n", log_file = log_file)

  # Add comparison summary
  sensitivity_results$comparison_summary <- summarize_carryover_sensitivity(sensitivity_results, trial_designs)

  return(sensitivity_results)
}

#' Summarize carryover sensitivity results
#'
#' Creates a summary of sensitivity analysis results across different
#' carryover models and scenarios, formatted as a data frame for easy
#' visualization and analysis.
#'
#' @param sensitivity_results Results from run_carryover_sensitivity
#' @param trial_designs Original trial designs for naming
#' @return Data frame with summary statistics
#' @importFrom dplyr group_by summarise bind_rows
#' @importFrom tidyr %>%
#' @export
summarize_carryover_sensitivity <- function(sensitivity_results, trial_designs) {

  # Extract design names
  design_names <- sapply(trial_designs, function(x) x$metadata$name_shortform)

  summary_data <- data.frame(stringsAsFactors = FALSE)

  for (model_name in names(sensitivity_results)) {
    if (model_name == "comparison_summary") next

    model_data <- sensitivity_results[[model_name]]

    for (scenario_name in names(model_data)) {
      scenario_data <- model_data[[scenario_name]]$results
      scenario_params <- model_data[[scenario_name]]$parameters

      # Calculate summary statistics for each design
      design_summary <- scenario_data %>%
        group_by(trialdesign) %>%
        summarise(
          model = model_name,
          scenario = scenario_name,
          design_name = design_names[trialdesign[1]],
          power = mean(p < 0.05, na.rm = TRUE),
          mean_beta = mean(beta, na.rm = TRUE),
          sd_beta = sd(beta, na.rm = TRUE),
          mean_beta_se = mean(betaSE, na.rm = TRUE),
          .groups = "drop"
        )

      # Add scenario parameters as columns
      for (param_name in names(scenario_params)) {
        design_summary[[paste0("param_", param_name)]] <- scenario_params[[param_name]]
      }

      summary_data <- bind_rows(summary_data, design_summary)
    }
  }

  return(summary_data)
}

#' Plot carryover sensitivity results
#'
#' Creates visualizations to compare how different carryover models
#' and parameters affect key outcome metrics across trial designs.
#'
#' @param sensitivity_summary Output from summarize_carryover_sensitivity
#' @param metric Metric to plot ("power", "mean_beta", "sd_beta")
#' @return invisible(NULL)
#' @importFrom rlang sym
#' @export
plot_carryover_sensitivity <- function(sensitivity_summary, metric = "power") {

  models <- unique(sensitivity_summary$model)
  n_models <- length(models)
  design_names <- unique(sensitivity_summary$design_name)
  cols <- seq_along(design_names)

  old_par <- par(mfrow = c(1, n_models), mar = c(6, 4, 3, 1),
                 oma = c(0, 0, 3, 0))
  on.exit(par(old_par), add = TRUE)

  metric_label <- gsub("_", " ", metric)
  metric_label <- paste0(toupper(substring(metric_label, 1, 1)),
                         substring(metric_label, 2))

  for (m in models) {
    sub <- sensitivity_summary[sensitivity_summary$model == m, ]
    scenarios <- unique(sub$scenario)
    scenario_idx <- seq_along(scenarios)

    y_range <- range(sub[[metric]], na.rm = TRUE)

    plot(range(scenario_idx), y_range, type = "n",
         xlab = "", ylab = metric_label,
         main = paste0(toupper(substring(m, 1, 1)), substring(m, 2)),
         xaxt = "n")
    axis(1, at = scenario_idx, labels = scenarios, las = 2, cex.axis = 0.8)

    for (di in seq_along(design_names)) {
      dn <- design_names[di]
      d_sub <- sub[sub$design_name == dn, ]
      x_pos <- match(d_sub$scenario, scenarios)
      lines(x_pos, d_sub[[metric]], col = cols[di], lwd = 1.5)
      points(x_pos, d_sub[[metric]], col = cols[di], pch = 16, cex = 1.2)
    }

    if (metric == "power") {
      axis_labels <- pretty(y_range)
      axis(2, at = axis_labels, labels = paste0(axis_labels * 100, "%"))
    }
  }

  title(paste("Carryover Model Sensitivity:", metric_label), outer = TRUE)

  # Add legend to last panel
  legend("bottomright", legend = design_names, col = cols,
         lty = 1, pch = 16, cex = 0.7, title = "Trial Design")

  invisible(NULL)
}

#' Compare carryover models against observed data
#'
#' Fits different carryover models to observed data with known
#' carryover effects, and evaluates which model best explains
#' the observed patterns.
#'
#' @param observed_data Observed data with carryover effects
#' @param time_column Name of column containing time values
#' @param response_column Name of column containing response values
#' @param carryover_models Vector of carryover models to test
#' @param log_file Optional file path to write progress logs
#' @return List with model comparisons and fit statistics
#' @importFrom stats nls AIC BIC
#' @importFrom purrr map_dfr
#' @export
compare_carryover_models <- function(
  observed_data,
  time_column = "time",
  response_column = "response",
  carryover_models = c("exponential", "linear", "weibull", "power"),
  log_file = NULL
) {
  # Extract data
  time_values <- observed_data[[time_column]]
  response_values <- observed_data[[response_column]]

  # Make sure we have valid input
  if (length(time_values) < 3 || length(response_values) < 3) {
    stop("Insufficient data points for model fitting (need at least 3)")
  }

  # Fit each carryover model
  model_fits <- list()
  model_stats <- data.frame(stringsAsFactors = FALSE)

  # Initial parameter estimates
  y0 <- max(response_values, na.rm = TRUE)  # Initial effect
  t_half_estimate <- max(time_values, na.rm = TRUE) / 2

  log_progress("Fitting carryover models to observed data...\n", log_file = log_file)

  for (model_name in carryover_models) {
    log_progress("  Fitting ", model_name, " model... ", log_file = log_file)

    tryCatch({
      # Define model formula based on carryover type
      model_formula <- switch(model_name,
        "exponential" = response ~ y0 * (1/2)^(time/t_half),
        "linear" = response ~ y0 * pmax(0, 1 - time/t_max),
        "weibull" = response ~ y0 * exp(-(time/lambda)^k),
        "power" = response ~ y0 * pmax(0, (1 - time/t_max)^p),
        stop("Unknown model type: ", model_name)
      )

      # Set up starting parameters
      start_params <- switch(model_name,
        "exponential" = list(y0 = y0, t_half = t_half_estimate),
        "linear" = list(y0 = y0, t_max = max(time_values) * 1.5),
        "weibull" = list(y0 = y0, lambda = t_half_estimate, k = 1),
        "power" = list(y0 = y0, t_max = max(time_values) * 1.5, p = 1.5)
      )

      # Fit model
      fit <- nls(
        model_formula,
        data = data.frame(time = time_values, response = response_values),
        start = start_params,
        control = nls.control(maxiter = 100)
      )

      # Store the fit
      model_fits[[model_name]] <- fit

      # Calculate fit statistics
      model_stats <- bind_rows(model_stats, data.frame(
        model = model_name,
        aic = AIC(fit),
        bic = BIC(fit),
        rss = sum(residuals(fit)^2),
        parameters = length(coef(fit)),
        converged = TRUE,
        stringsAsFactors = FALSE
      ))

      log_progress("success (AIC = ", round(AIC(fit), 2), ")\n", log_file = log_file)

    }, error = function(e) {
      # Record failed fit
      log_progress("failed: ", as.character(e), "\n", log_file = log_file)

      model_stats <<- bind_rows(model_stats, data.frame(
        model = model_name,
        aic = NA,
        bic = NA,
        rss = NA,
        parameters = NA,
        converged = FALSE,
        error = as.character(e),
        stringsAsFactors = FALSE
      ))
    })
  }

  # Add relative model performance measure (smaller is better)
  if (nrow(model_stats) > 0 && any(!is.na(model_stats$aic))) {
    min_aic <- min(model_stats$aic, na.rm = TRUE)
    model_stats <- model_stats %>%
      mutate(delta_aic = aic - min_aic)

    # Find best model
    best_model <- model_stats %>%
      filter(!is.na(aic)) %>%
      arrange(aic) %>%
      slice(1) %>%
      pull(model)

    log_progress("\nBest fitting model: ", best_model, " (lowest AIC)\n", log_file = log_file)
  } else {
    log_progress("\nNo models converged successfully\n", log_file = log_file)
    best_model <- NULL
  }

  # Return model comparison results
  return(list(
    model_stats = model_stats,
    best_model = best_model,
    model_fits = model_fits
  ))
}

#' Create comparison plot of carryover models
#'
#' Visualizes different carryover models on the same plot to compare
#' their decay profiles over time.
#'
#' @param time_range Vector of time values to plot
#' @param initial_effect Initial effect value at time 0
#' @param model_params List of model parameters for each model type
#' @param models Vector of models to include
#' @return invisible(NULL)
#' @export
plot_carryover_models <- function(
  time_range = seq(0, 10, by = 0.1),
  initial_effect = 10,
  model_params = list(
    exponential = list(carryover_t1half = 2, scalefactor = 2),
    linear = list(total_time = 8),
    weibull = list(carryover_t1half = 2, k = 1.5),
    power = list(max_time = 8, power = 1.5)
  ),
  models = c("exponential", "linear", "weibull", "power")
) {
  # Generate data for plotting
  plot_data <- list()

  for (model in models) {
    if (!(model %in% names(model_params))) next

    model_values <- sapply(time_range, function(t) {
      if (t == 0) return(initial_effect)
      calculate_carryover(t, initial_effect, model, model_params[[model]])
    })

    plot_data[[model]] <- data.frame(
      time = time_range,
      effect = model_values,
      stringsAsFactors = FALSE
    )
  }

  # Create visualization
  cols <- seq_along(models)

  y_range <- range(unlist(lapply(plot_data, function(d) d$effect)),
                   na.rm = TRUE)

  old_par <- par(mar = c(5, 4, 4, 2))
  on.exit(par(old_par), add = TRUE)

  plot(range(time_range), y_range, type = "n",
       xlab = "Time Since Discontinuation",
       ylab = "Remaining Effect",
       main = "Comparison of Carryover Effect Models")

  for (i in seq_along(names(plot_data))) {
    m <- names(plot_data)[i]
    lines(plot_data[[m]]$time, plot_data[[m]]$effect,
          col = cols[i], lwd = 2)
  }

  legend("topright", legend = names(plot_data), col = cols,
         lty = 1, lwd = 2, title = "Model Type")

  invisible(NULL)
}

# Example usage script
if (FALSE) {
  # Example of how to use the sensitivity analysis

  # Set up base parameters (without carryover model specification)
  base_model_params <- data.frame(
    N = 30,
    c.bm = 0.7,
    carryover_t1half = 1,  # This will be overridden by scenarios
    c.tv = 0.9,
    c.pb = 0.8,
    c.br = 0.85,
    c.cf1t = 0.3,
    c.cfct = 0.1,
    stringsAsFactors = FALSE
  )

  # Define custom carryover scenarios if desired
  custom_scenarios <- list(
    exponential = list(
      list(carryover_t1half = 0.25, scale_factor = 2),
      list(carryover_t1half = 0.5, scale_factor = 2),
      list(carryover_t1half = 1.0, scale_factor = 2),
      list(carryover_t1half = 2.0, scale_factor = 2)
    ),
    linear = list(
      list(carryover_t1half = 1.0, total_time = 1),
      list(carryover_t1half = 1.0, total_time = 2),
      list(carryover_t1half = 1.0, total_time = 4),
      list(carryover_t1half = 1.0, total_time = 8)
    ),
    weibull = list(
      list(carryover_t1half = 1.0, k = 0.5),  # Very slow initial decay
      list(carryover_t1half = 1.0, k = 1.0),  # Exponential
      list(carryover_t1half = 1.0, k = 1.5),  # Moderate acceleration
      list(carryover_t1half = 1.0, k = 2.0)   # Strong acceleration
    ),
    power = list(
      list(carryover_t1half = 1.0, power = 0.5),  # Slower initial drop
      list(carryover_t1half = 1.0, power = 1.0),  # Linear
      list(carryover_t1half = 1.0, power = 1.5),  # Moderate curve
      list(carryover_t1half = 1.0, power = 2.0)   # Quadratic decay
    )
  )

  # Set up a log file for progress reporting
  log_file <- "carryover_sensitivity_analysis.log"

  # Run sensitivity analysis
  sensitivity_results <- run_carryover_sensitivity(
    trial_designs = designs,
    base_model_params = base_model_params,
    resp_param_sets = resp_param_sets,
    bl_param_sets = bl_param_sets,
    sim_options = list(n_reps = 50, log_file = log_file),
    carryover_scenarios = custom_scenarios,
    log_file = log_file
  )

  # Summarize results
  summary <- sensitivity_results$comparison_summary

  # Create plots
  power_plot <- plot_carryover_sensitivity(summary, "power")
  beta_plot <- plot_carryover_sensitivity(summary, "mean_beta")
  variance_plot <- plot_carryover_sensitivity(summary, "sd_beta")

  # Visualize model shapes
  model_comparison <- plot_carryover_models(
    time_range = seq(0, 8, by = 0.1),
    initial_effect = 10,
    model_params = list(
      exponential = list(carryover_t1half = 2, scale_factor = 2),
      linear = list(total_time = 8),
      weibull = list(carryover_t1half = 2, k = 1.5),
      power = list(max_time = 8, power = 1.5)
    )
  )

  # Print summary table
  print(summary)
}

#' Backward compatibility wrapper for run_carryover_sensitivity
#'
#' @param trial_designs List of trial designs to compare
#' @param base_model_params Base model parameters
#' @param resp_param_sets Response parameter sets
#' @param bl_param_sets Baseline parameter sets
#' @param sim_options Simulation options
#' @param analysis_options Analysis options
#' @param carryover_models Vector of carryover models to test
#' @param carryover_scenarios List of parameter scenarios for each model
#' @param log_file Optional file path to write progress logs
#' @return List containing sensitivity analysis results
#' @export
run_carryover_sensitivity_original <- function(
  trial_designs_old,
  base_model_params_old,
  resp_param_sets_old,
  bl_param_sets_old,
  sim_options_old = list(n_reps = 25),
  analysis_options_old = list(use_expectancy = TRUE, random_slope = FALSE),
  carryover_models_old = c("exponential", "linear", "weibull", "power"),
  carryover_scenarios_old = NULL,
  log_file = NULL
) {
  # Call the new function with renamed parameters
  run_carryover_sensitivity(
    trial_designs = trial_designs_old,
    base_model_params = base_model_params_old,
    resp_param_sets = resp_param_sets_old,
    bl_param_sets = bl_param_sets_old,
    sim_options = sim_options_old,
    analysis_options = analysis_options_old,
    carryover_models = carryover_models_old,
    carryover_scenarios = carryover_scenarios_old,
    log_file = log_file
  )
}
