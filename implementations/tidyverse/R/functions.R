#' functions.R: Tidyverse implementation of the pmsimstats simulation framework
#'
#' Complete reimplementation supporting both DGP architectures (mvn and
#' mean_moderation) via the dgp_architecture parameter.

#===========================================================================
# Function: mod_gompertz
# Origin-passing Gompertz (matches orig modgompertz exactly)
#===========================================================================
mod_gompertz <- function(t, maxr, disp, rate) {
  # Guard: maxr=0 implies y identically zero; otherwise the final
  # rescale divides by (maxr - vert_offset) = 0 and returns NaN.
  if (isTRUE(all.equal(maxr, 0))) return(rep(0, length(t)))
  y <- maxr * exp(-disp * exp(-rate * t))
  vert_offset <- maxr * exp(-disp * exp(-rate * 0))
  y <- y - vert_offset
  y * (maxr / (maxr - vert_offset))
}

#' Residual drug-effect decay multiplier
#'
#' Returns the fraction of residual drug effect present at time-since-
#' discontinuation `tsd`, for a carryover with given `halflife` and
#' functional form.
#'
#' @param tsd Numeric, time since drug discontinuation (same time
#'   units as `halflife`). May be a vector.
#' @param halflife Numeric, carryover half-life. The value of `tsd`
#'   at which the multiplier equals 0.5 under the exponential and
#'   Weibull forms, and at which the multiplier equals 0.5 under the
#'   linear form (which reaches zero at `2 * halflife`).
#' @param form Character, one of `"exponential"` (default),
#'   `"linear"`, or `"weibull"`.
#' @param shape Numeric, Weibull shape parameter `k`. Ignored unless
#'   `form = "weibull"`. `k = 1` recovers the exponential form;
#'   `k < 1` gives a heavier tail (slow elimination); `k > 1` gives
#'   an accelerated washout.
#' @return Numeric vector of the same length as `tsd`, with values
#'   in `[0, 1]`.
#'
#' @details
#' - Exponential: \eqn{\phi(t) = \exp(-\lambda t)} with
#'   \eqn{\lambda = \ln 2 / t_{1/2}}.
#' - Linear: \eqn{\phi(t) = \max(0, 1 - t / (2 t_{1/2}))}.
#' - Weibull: \eqn{\phi(t) = \exp(-(\lambda_w t)^k)} with
#'   \eqn{\lambda_w = (\ln 2)^{1/k} / t_{1/2}}.
#'
#' All three forms satisfy \eqn{\phi(0) = 1} and
#' \eqn{\phi(t_{1/2}) = 0.5}.
carryover_decay <- function(tsd, halflife,
                            form = c("exponential", "linear", "weibull"),
                            shape = 1) {
  form <- match.arg(form)
  if (is.null(halflife) || halflife == 0) return(rep(0, length(tsd)))
  switch(form,
    exponential = exp(-log(2) / halflife * tsd),
    linear      = pmax(0, 1 - tsd / (2 * halflife)),
    weibull     = {
      if (!is.numeric(shape) || shape <= 0) {
        stop("Weibull `shape` must be a positive number.")
      }
      lambda_w <- (log(2))^(1 / shape) / halflife
      exp(-(lambda_w * tsd)^shape)
    }
  )
}

resolve_lambda_cor <- function(lambda_cor, model_param) {
  if (!is.na(lambda_cor)) return(lambda_cor)
  if (!is.null(model_param$carryover_t1half) &&
      model_param$carryover_t1half > 0) {
    log(2) / model_param$carryover_t1half
  } else {
    0
  }
}

ensure_timepoint_name <- function(td) {
  if (!'timepoint_name' %in% names(td)) {
    nm <- intersect(c('timeptnames', 'timeptname'), names(td))
    if (length(nm)) td$timepoint_name <- td[[nm[1]]]
  }
  td
}

ensure_timeptnames <- function(td) {
  nms <- names(td)
  if ('timeptnames' %in% nms) return(td)
  src <- intersect(c('timeptname', 'timepoint_name'), nms)
  if (length(src)) {
    td$timeptnames <- td[[src[1]]]
    td[[src[1]]] <- NULL
  }
  td
}

get_timepoint_names <- function(td) {
  nms <- names(td)
  if ('timeptnames' %in% nms) return(td$timeptnames)
  if ('timeptname' %in% nms) return(td$timeptname)
  td$timepoint_name
}

compute_means <- function(baseline_param, resp_param, trial_data,
                          trial_design, factors, halflife,
                          carryover_form, weibull_shape) {
  gompertz_by_cat <- function(cat, t_vec) {
    mod_gompertz(
      t_vec,
      resp_param$max[resp_param$cat == cat],
      resp_param$disp[resp_param$cat == cat],
      resp_param$rate[resp_param$cat == cat]
    )
  }

  tv_means <- gompertz_by_cat("tv", trial_data$t_wk_cumulative)
  pb_means <- gompertz_by_cat("pb", trial_data$tpb) * trial_design$e
  br_means <- gompertz_by_cat("br", trial_data$tod)
  br_means <- apply_carryover_to_component(
    br_means, trial_data, halflife,
    form = carryover_form, shape = weibull_shape
  )

  c(
    baseline_param$m[baseline_param$cat == "bm"],
    baseline_param$m[baseline_param$cat == "BL"],
    tv_means, pb_means, br_means
  )
}

#' Carryover calculation for bio-response component (Hendrickson method)
#'
#' @param component_means Current component means
#' @param trial_data Trial design data
#' @param halflife Carryover half-life value
#' @param form Decay form; see [carryover_decay()].
#' @param shape Weibull shape; see [carryover_decay()].
#' @return Adjusted means with carryover effects
#'
#' @details Following Hendrickson et al. (2020), carryover is applied
#' ONLY to the bio-response (br) component when participants are off
#' drug. The recurrence is
#' `mu[t] = base[t] + mu[t-1] * phi(tsd; halflife, form, shape)`
#' where `phi` is the decay multiplier from [carryover_decay()]. The
#' default `form = "exponential"` reproduces the Hendrickson
#' `(1/2)^(tsd/t1half)` behaviour exactly.
apply_carryover_to_component <- function(
    component_means, trial_data, halflife,
    form = "exponential", shape = 1) {

  num_timepoints <- length(component_means)

  if (is.null(halflife) || halflife == 0 || num_timepoints <= 1) {
    return(component_means)
  }

  carryover_indices <- which(!trial_data$on_drug & trial_data$tsd > 0)

  if (length(carryover_indices) > 0) {
    for (idx in carryover_indices) {
      prev_idx <- idx - 1
      time_lag <- trial_data$tsd[idx]
      decay_factor <- carryover_decay(time_lag, halflife,
                                      form = form, shape = shape)
      component_means[idx] <- component_means[idx] +
        component_means[prev_idx] * decay_factor
    }
  }

  component_means
}

#===========================================================================
# Helper Functions for Data Generation
#===========================================================================

prepare_trial_data <- function(trial_design) {
  if (!is.data.frame(trial_design)) {
    stop("Trial design is not a data frame. Class: ", class(trial_design),
         ", Type: ", typeof(trial_design))
  }

  trial_data <- as_tibble(trial_design)

  if (!("t_wk" %in% names(trial_data))) {
    stop("'t_wk' column not found in trial design data")
  }

  trial_data$t_wk_cumulative <- cumsum(trial_data$t_wk)
  trial_data$on_drug <- (trial_data$tod > 0)
  trial_data
}

build_correlation_matrix <- function(
    labels, trial_design, model_param, num_timepoints,
    factors, trial_data, lambda_cor = 0,
    carryover_form = "exponential", weibull_shape = 1,
    dgp_architecture = "mvn") {
  correlations <- diag(length(labels))
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  weeks <- trial_data$t_wk_cumulative

  if (num_timepoints > 1) {
    upper_pairs <- expand.grid(
      p1 = 1:(num_timepoints - 1),
      p2 = 2:num_timepoints
    )
    upper_pairs <- upper_pairs[upper_pairs$p2 > upper_pairs$p1, ]
  }

  for (factor_idx in seq_along(factors)) {
    current_factor <- factors[factor_idx]

    if (num_timepoints > 1) {
      autocorrelation <- model_param[[paste("c", current_factor, sep = ".")]]
      time_gaps <- abs(weeks[upper_pairs$p2] - weeks[upper_pairs$p1])
      ar_vals <- autocorrelation^time_gaps

      n1 <- paste(trial_design$timepoint_name[upper_pairs$p1],
                  current_factor, sep = ".")
      n2 <- paste(trial_design$timepoint_name[upper_pairs$p2],
                  current_factor, sep = ".")
      idx1 <- cbind(match(n1, labels), match(n2, labels))
      idx2 <- cbind(idx1[, 2], idx1[, 1])
      correlations[idx1] <- ar_vals
      correlations[idx2] <- ar_vals
    }

    for (factor2_idx in setdiff(seq_along(factors), factor_idx)) {
      other_factor <- factors[factor2_idx]

      n1_same <- paste(trial_design$timepoint_name,
                       current_factor, sep = ".")
      n2_same <- paste(trial_design$timepoint_name,
                       other_factor, sep = ".")
      idx1_same <- cbind(match(n1_same, labels), match(n2_same, labels))
      idx2_same <- cbind(idx1_same[, 2], idx1_same[, 1])
      correlations[idx1_same] <- model_param$c.cf1t
      correlations[idx2_same] <- model_param$c.cf1t

      if (num_timepoints > 1) {
        n1_lag <- paste(trial_design$timepoint_name[upper_pairs$p1],
                        current_factor, sep = ".")
        n2_lag <- paste(trial_design$timepoint_name[upper_pairs$p2],
                        other_factor, sep = ".")
        lag_vals <- model_param$c.cfct * ar_vals
        idx1_lag <- cbind(match(n1_lag, labels), match(n2_lag, labels))
        idx2_lag <- cbind(idx1_lag[, 2], idx1_lag[, 1])
        correlations[idx1_lag] <- lag_vals
        correlations[idx2_lag] <- lag_vals
      }
    }

    # BM-BR differential correlation (Architecture B / MVN and Architecture C)
    # Under "mvn" the weight is c.bm; under "combined" it is c.bm_b.
    if (current_factor == "br" && dgp_architecture %in% c("mvn", "combined")) {
      c_bm_cov <- if (dgp_architecture == "combined") model_param$c.bm_b
                  else model_param$c.bm
      for (timepoint_idx in seq_len(num_timepoints)) {
        name1 <- paste(trial_design$timepoint_name[
                         timepoint_idx], "br", sep = ".")

        on_drug_now <- trial_data$on_drug[timepoint_idx]

        if (on_drug_now) {
          correlations["bm", name1] <- c_bm_cov
          correlations[name1, "bm"] <- c_bm_cov
        } else {
          tsd_now <- trial_data$tsd[timepoint_idx]
          if (tsd_now > 0 && lambda_cor > 0) {
            halflife_for_decay <- log(2) / lambda_cor
            decay <- carryover_decay(
              tsd_now, halflife_for_decay,
              form = carryover_form, shape = weibull_shape
            )
            correlations["bm", name1] <- c_bm_cov * decay
            correlations[name1, "bm"] <- c_bm_cov * decay
          }
        }
      }
    }
  }

  correlations
}

#===========================================================================
# Helper function for building and caching sigma matrices
#===========================================================================
build_sigma_matrix <- function(model_param, resp_param, baseline_param,
                               trial_design,
                               verbose = FALSE,
                               lambda_cor = NA,
                               dgp_architecture = "mvn",
                               ...) {

  dgp_architecture <- match.arg(dgp_architecture, c("mvn", "mean_moderation", "combined"))

  factors <- c("tv", "pb", "br")

  trial_data <- prepare_trial_data(trial_design)
  num_timepoints <- dim(trial_design)[1]

  lambda_cor <- resolve_lambda_cor(lambda_cor, model_param)

  labels <- c(
    c("bm", "BL"),
    paste(trial_design$timepoint_name, factors[1], sep = "."),
    paste(trial_design$timepoint_name, factors[2], sep = "."),
    paste(trial_design$timepoint_name, factors[3], sep = ".")
  )

  standard_deviations <- c(
    baseline_param$sd[baseline_param$cat == "bm"],
    baseline_param$sd[baseline_param$cat == "BL"],
    rep(resp_param$sd[resp_param$cat == "tv"], num_timepoints),
    rep(resp_param$sd[resp_param$cat == "pb"], num_timepoints) * trial_design$e,
    rep(resp_param$sd[resp_param$cat == "br"], num_timepoints)
  )

  halflife <- model_param$carryover_t1half
  carryover_form <- model_param$carryover_form %||% "exponential"
  weibull_shape <- model_param$weibull_shape %||% 1

  means <- compute_means(
    baseline_param, resp_param, trial_data,
    trial_design, factors, halflife,
    carryover_form, weibull_shape
  )

  correlations <- build_correlation_matrix(
    labels, trial_design, model_param, num_timepoints,
    factors, trial_data = trial_data, lambda_cor = lambda_cor,
    carryover_form = carryover_form, weibull_shape = weibull_shape,
    dgp_architecture = dgp_architecture
  )

  sigma <- outer(standard_deviations, standard_deviations) * correlations

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
    cl = factors,
    chol_sigma = chol_sigma,
    was_pd_corrected = !is_pd
  )
}

#===========================================================================
# Function: generate_data
# Primary data generation function for trial simulation
#===========================================================================
generate_data <- function(
    model_param, resp_param, baseline_param, trial_design,
    empirical, make_positive_definite, seed = NA,
    lambda_cor = NA, verbose = FALSE, track_pd_stats = TRUE,
    cached_sigma = NULL, dgp_architecture = "mvn") {

  dgp_architecture <- match.arg(dgp_architecture, c("mvn", "mean_moderation", "combined"))

  lambda_cor <- resolve_lambda_cor(lambda_cor, model_param)

  factors <- c("tv", "pb", "br")

  trial_data <- prepare_trial_data(trial_design)

  if (!is.null(cached_sigma)) {
    labels <- cached_sigma$labels
    means <- cached_sigma$means
    chol_sigma <- cached_sigma$chol_sigma
    was_pd_corrected <- FALSE
  } else {
    sigma_result <- build_sigma_matrix(
      model_param, resp_param, baseline_param, trial_design,
      lambda_cor = lambda_cor, verbose = verbose,
      dgp_architecture = dgp_architecture
    )
    labels <- sigma_result$labels
    means <- sigma_result$means
    chol_sigma <- sigma_result$chol_sigma
    was_pd_corrected <- sigma_result$was_pd_corrected
  }

  if (!is.na(seed)) {
    set.seed(seed)
  }

  n <- model_param$N
  p <- length(means)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  participant_data <- Z %*% chol_sigma +
    matrix(means, nrow = n, ncol = p, byrow = TRUE)
  participant_data <- as_tibble(participant_data)
  colnames(participant_data) <- labels

  # Architecture A / C (mean-moderation channel): additive biomarker
  # moderation of BR. Under "mean_moderation" weight = c.bm;
  # under "combined" weight = c.bm_a.
  if (dgp_architecture %in% c("mean_moderation", "combined")) {
    beta_bm <- if (dgp_architecture == "combined") model_param$c.bm_a
               else model_param$c.bm
    bm_mean <- baseline_param |> filter(cat == "bm") |> pull(m)
    bm_sd <- baseline_param |> filter(cat == "bm") |> pull(sd)
    br_sd <- resp_param |> filter(cat == "br") |> pull(sd)
    bm_z <- (participant_data$bm - bm_mean) / bm_sd

    for (tp in seq_len(nrow(trial_design))) {
      if (trial_data$on_drug[tp]) {
        br_col <- paste(trial_design$timepoint_name[tp], "br", sep = ".")
        participant_data[[br_col]] <- participant_data[[br_col]] +
          beta_bm * bm_z * br_sd
      }
    }
  }

  participant_data <- process_participant_data(
    participant_data, trial_design$timepoint_name,
    factors, model_param$N
  )

  if (track_pd_stats) {
    was_corrected <- if (!is.null(cached_sigma)) FALSE else was_pd_corrected
    attr(participant_data, "sigma_count") <- 1L
    attr(participant_data, "non_positive_definite_count") <-
      as.integer(was_corrected)
    attr(participant_data, "non_positive_definite_rate") <-
      as.numeric(was_corrected)
  }

  participant_data
}

#===========================================================================
# Helper function for processing participant data
#===========================================================================
process_participant_data <- function(participant_data, timepoint_names,
                                     factors, N) {
  participant_data <- participant_data |>
    mutate(ptID = 1:N)

  for (timepoint_idx in seq_along(timepoint_names)) {
    tp <- timepoint_names[timepoint_idx]
    components <- paste(tp, factors, sep = ".")
    delta_col <- paste0("D_", tp)
    participant_data <- participant_data |>
      mutate(!!delta_col := rowSums(pick(all_of(components))),
             !!tp := BL - .data[[delta_col]])
  }

  participant_data
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
      dplyr::filter(path == g) |>
      as_tibble()

    td_names <- get_timepoint_names(td)

    td_with_bl <- bind_rows(
      tibble(
        timeptnames = 'BL', t_wk = 0, e = 0,
        tod = 0, tsd = 0, tpb = 0
      ),
      as_tibble(ensure_timeptnames(td))
    ) |>
      mutate(t = cumsum(t_wk))

    all_names <- c('BL', td_names)
    valid_names <- all_names[all_names %in% names(dat_single)]

    data_wide <- dat_single |>
      dplyr::select(ptID, bm, all_of(valid_names)) |>
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
        td_with_bl |> dplyr::select(timeptnames, t, De = e, tod, tsd),
        by = 'timeptnames'
      ) |>
      mutate(
        Db = (tod > 0),
        # Guard against 0/0 = NaN when carryover_t1half == 0 or
        # tsd == 0 for never-on-drug subjects; both collapse to
        # Dbc=0, recovering the binary on/off indicator and
        # preserving the never-on-drug arm for main-effect tests.
        Dbc = case_when(
          Db ~ 1,
          op$carryover_t1half == 0 ~ 0,
          tsd <= 0 ~ 0,
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
      dplyr::filter(t > 0, mean_Db == 1) |>
      pull(ptID) |>
      unique()
    datamerged <- datamerged |>
      dplyr::filter(ptID %in% ever_on_drug)
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
  datamerged <- datamerged |>
    dplyr::filter(if_all(all_of(model_vars), ~ !is.na(.x)))

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

  target_names <- if (var_in_Db) c('bm:Dbc', 'Dbc:bm') else c('bm:t', 't:bm')
  target <- intersect(target_names, coefnames)
  if (length(target) == 0) {
    beta <- NA_real_; betaSE <- NA_real_; p <- NA_real_
  } else {
    target <- target[1]
    p <- cc[target, 'p-value']
    beta <- cc[target, 'Value']
    betaSE <- cc[target, 'Std.Error']
  }

  out <- tibble(
    beta = beta, betaSE = betaSE, p = p,
    issingular = FALSE, warning = NA_character_
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

  t_wk <- c(timepoints[1],
             diff(timepoints))

  trialpaths <- map(ondrug, function(od) {
    od_intervals <- t_wk * od
    tod <- od_intervals
    tsd <- t_wk * (od != 1)
    tpb <- t_wk * (expectancies > 0)

    if (length(timepoints) > 1) {
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
  tp_names <- get_timepoint_names(td)

  delta_cols <- paste0('D_', tp_names)
  cdt <- dat |> dplyr::select(all_of(delta_cols))

  frac_NA <- censorparam$beta0
  frac_NA_biased <- censorparam$beta1
  fna1 <- frac_NA * (1 - frac_NA_biased)
  fna2 <- frac_NA * frac_NA_biased

  cdt_ps2 <- sapply(cdt, function(x) (x + 100)^censorparam$eb1)
  cdt_p1 <- matrix(td$t_wk, nrow = nrow(cdt), ncol = ncol(cdt),
                    byrow = TRUE)
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
    lambda_cor = NA, n_cores = 1, dgp_architecture = "mvn") {

  dgp_architecture <- match.arg(dgp_architecture, c("mvn", "mean_moderation", "combined"))

  if (missing(analysisparams)) {
    analysisparams <- list(useDE = TRUE, t_random_slope = FALSE,
                           full_model_out = FALSE)
  }

  if (n_cores < 1) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

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
        td_path <- ensure_timepoint_name(td[[i_p]])
        sigma_cache[[cache_key]] <- build_sigma_matrix(
          mpp, rp, bpp, td_path,
          lambda_cor = lambda_cor,
          dgp_architecture = dgp_architecture
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

      td_path <- ensure_timepoint_name(td[[1]])

      dat <- generate_data(mpp_copy, rp, bpp, td_path,
                           empirical = FALSE,
                           make_positive_definite = TRUE,
                           cached_sigma = sigma_cache[[ck]],
                           dgp_architecture = dgp_architecture)
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
          td_p <- ensure_timepoint_name(td[[i_p]])
          dat2 <- generate_data(mpp_copy, rp, bpp, td_p,
                                empirical = FALSE,
                                make_positive_definite = TRUE,
                                cached_sigma = sigma_cache[[ck]],
                                dgp_architecture = dgp_architecture)
          dat2$path <- i_p
          dat2$replicate <- rep(1:simparam$Nreps, Ns[i_p])
          dat <- bind_rows(dat, dat2)
        }
      }
      mpp_copy$N <- sum(Ns)

      make_result_row <- function(ap_row, et_out, a_out, i_s,
                                    censorparamset, frac_na) {
        bind_cols(
          as_tibble(vpg[i_r, , drop = FALSE]),
          as_tibble(as.list(mpp_copy)),
          tibble(
            censorparamset = censorparamset,
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

      results_list <- list()
      ridx <- 0
      nocensoring <- length(censorparams) < 2 && is.na(censorparams)

      for (i_ap in seq_len(nrow(analysisparams))) {
        ap_row <- as.list(analysisparams[i_ap, ])
        et_out <- lme_analysis(td, dat, ap_row)
        for (i_s in seq_len(simparam$Nreps)) {
          dat_rep <- dat |> dplyr::filter(replicate == i_s)
          a_out <- lme_analysis(td, dat_rep, ap_row)
          ridx <- ridx + 1
          results_list[[ridx]] <- make_result_row(
            ap_row, et_out, a_out, i_s, 0, 0
          )
        }

        if (!nocensoring) {
          for (i_c in seq_len(nrow(censorparams))) {
            datc <- censor_data(dat, td[[1]], censorparams[i_c, ])
            for (i_s in seq_len(simparam$Nreps)) {
              datc_rep <- datc |> dplyr::filter(replicate == i_s)
              frac_na <- sum(is.na(datc_rep)) /
                (mpp_copy$N * nrow(td[[1]]))
              a_out <- lme_analysis(td, datc_rep, ap_row)
              ridx <- ridx + 1
              results_list[[ridx]] <- make_result_row(
                ap_row, et_out, a_out, i_s, i_c, frac_na
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
      saveRDS(outpt, file.path(
        simparam$savedir,
        paste(simparam$basesavename, i_ll, sep = '_save')
      ))
    }
  }

  if (!simparam$progressiveSave) outpt
}
