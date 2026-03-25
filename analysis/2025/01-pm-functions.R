#' 01-pm-functions.R: Common functions for precision medicine trial simulation
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

#' Carryover calculation for bio-response component (Hendrickson method)
#'
#' @param component_means Current component means
#' @param trial_data Trial design data
#' @param halflife Carryover half-life value
#' @return Adjusted means with carryover effects
#'
#' @details Following Hendrickson et al. (2020), carryover is applied
#' ONLY to the bio-response (br) component when participants are off
#' drug. The formula is: mu[t] = base[t] + mu[t-1] * (1/2)^(tsd/t1half)
#' where tsd is time since discontinuation.
apply_carryover_to_component <- function(
    component_means, trial_data, halflife) {

  num_timepoints <- length(component_means)

  if (is.null(halflife) || halflife == 0 || num_timepoints <= 1) {
    return(component_means)
  }

  carryover_indices <- which(!trial_data$on_drug & trial_data$tsd > 0)

  if (length(carryover_indices) > 0) {
    for (idx in carryover_indices) {
      prev_idx <- idx - 1
      if (prev_idx >= 1 &&
          prev_idx <= length(component_means) &&
          idx <= length(component_means)) {
        time_lag <- trial_data$tsd[idx]
        decay_factor <- (1/2)^(time_lag / halflife)
        component_means[idx] <- component_means[idx] +
          component_means[prev_idx] * decay_factor
      }
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
    factors, trial_data, lambda_cor = 0) {
  correlations <- diag(length(labels))
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  weeks <- trial_data$t_wk_cumulative

  for (factor_idx in seq_along(factors)) {
    current_factor <- factors[factor_idx]

    if (num_timepoints > 1) {
      autocorrelation <- model_param[[paste("c", current_factor, sep = ".")]]

      point_indices <- expand.grid(
        p1 = 1:(num_timepoints-1),
        p2 = (2:num_timepoints)
      )
      point_indices <- point_indices[
        point_indices$p2 > point_indices$p1,
      ]

      name1 <- paste(trial_design$timepoint_name[
                       point_indices$p1],
                     current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name[
                       point_indices$p2],
                     current_factor, sep = ".")

      for (idx in 1:nrow(point_indices)) {
        time_gap <- abs(weeks[point_indices$p2[idx]] -
                        weeks[point_indices$p1[idx]])
        correlations[name1[idx], name2[idx]] <- autocorrelation^time_gap
        correlations[name2[idx], name1[idx]] <- autocorrelation^time_gap
      }
    }

    for (factor2_idx in setdiff(seq_along(factors), factor_idx)) {
      other_factor <- factors[factor2_idx]

      name1 <- paste(trial_design$timepoint_name,
                     current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name,
                     other_factor, sep = ".")

      for (idx in 1:length(name1)) {
        correlations[name1[idx], name2[idx]] <- model_param$c.cf1t
        correlations[name2[idx], name1[idx]] <- model_param$c.cf1t
      }

      if (num_timepoints > 1) {
        point_indices <- expand.grid(
          p1 = 1:(num_timepoints-1),
          p2 = (2:num_timepoints)
        )
        point_indices <- point_indices[
          point_indices$p2 > point_indices$p1,
        ]

        name1 <- paste(trial_design$timepoint_name[
                         point_indices$p1],
                       current_factor, sep = ".")
        name2 <- paste(trial_design$timepoint_name[
                         point_indices$p2],
                       other_factor, sep = ".")

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

    if (current_factor == "br") {
      for (timepoint_idx in 1:num_timepoints) {
        name1 <- paste(trial_design$timepoint_name[
                         timepoint_idx], "br", sep = ".")

        on_drug_now <- trial_data$on_drug[timepoint_idx]

        if (on_drug_now) {
          correlations["bm", name1] <- model_param$c.bm
          correlations[name1, "bm"] <- model_param$c.bm
        } else {
          tsd_now <- trial_data$tsd[timepoint_idx]
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

  correlations
}

#===========================================================================
# Helper function for building and caching sigma matrices
#===========================================================================
build_sigma_matrix <- function(model_param, resp_param, baseline_param,
                               trial_design,
                               factor_types = NULL,
                               factor_abbreviations = NULL,
                               verbose = FALSE,
                               lambda_cor = NA) {

  factors <- c("tv", "pb", "br")

  trial_data <- prepare_trial_data(trial_design)
  num_timepoints <- dim(trial_design)[1]

  if (is.na(lambda_cor)) {
    if (!is.null(model_param$carryover_t1half) &&
        model_param$carryover_t1half > 0) {
      lambda_cor <- log(2) / model_param$carryover_t1half
    } else {
      lambda_cor <- 0
    }
  }

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

  means <- c(
    baseline_param$m[baseline_param$cat == "bm"],
    baseline_param$m[baseline_param$cat == "BL"]
  )

  for (f in factors) {
    if (f == "tv") {
      tv_means <- mod_gompertz(
        trial_data$t_wk_cumulative,
        resp_param$max[resp_param$cat == "tv"],
        resp_param$disp[resp_param$cat == "tv"],
        resp_param$rate[resp_param$cat == "tv"]
      )
      means <- c(means, tv_means)
    }

    if (f == "pb") {
      pb_means <- mod_gompertz(
        trial_data$tpb,
        resp_param$max[resp_param$cat == "pb"],
        resp_param$disp[resp_param$cat == "pb"],
        resp_param$rate[resp_param$cat == "pb"]
      ) * trial_design$e
      means <- c(means, pb_means)
    }

    if (f == "br") {
      br_means <- mod_gompertz(
        trial_data$tod,
        resp_param$max[resp_param$cat == "br"],
        resp_param$disp[resp_param$cat == "br"],
        resp_param$rate[resp_param$cat == "br"]
      )
      br_means <- apply_carryover_to_component(
        br_means, trial_data, halflife
      )
      means <- c(means, br_means)
    }
  }

  correlations <- build_correlation_matrix(
    labels, trial_design, model_param, num_timepoints,
    factors, trial_data = trial_data, lambda_cor = lambda_cor
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
    cached_sigma = NULL) {

  if (is.na(lambda_cor)) {
    if (!is.null(model_param$carryover_t1half) &&
        model_param$carryover_t1half > 0) {
      lambda_cor <- log(2) / model_param$carryover_t1half
    } else {
      lambda_cor <- 0
    }
  }

  factors <- c("tv", "pb", "br")

  if (!is.null(cached_sigma)) {
    sigma <- cached_sigma$sigma
    labels <- cached_sigma$labels
    standard_deviations <- cached_sigma$standard_deviations
    correlations <- cached_sigma$correlations

    trial_data <- prepare_trial_data(trial_design)
    num_timepoints <- dim(trial_design)[1]

  } else {
    sigma_result <- build_sigma_matrix(
      model_param, resp_param, baseline_param, trial_design,
      lambda_cor = lambda_cor, verbose = verbose
    )
    sigma <- sigma_result$sigma
    labels <- sigma_result$labels
    standard_deviations <- sigma_result$standard_deviations
    correlations <- sigma_result$correlations

    trial_data <- prepare_trial_data(trial_design)
    num_timepoints <- dim(trial_design)[1]
  }

  halflife <- model_param$carryover_t1half

  means <- c(
    baseline_param$m[baseline_param$cat == "bm"],
    baseline_param$m[baseline_param$cat == "BL"]
  )

  for (f in factors) {
    if (f == "tv") {
      tv_means <- mod_gompertz(
        trial_data$t_wk_cumulative,
        resp_param$max[resp_param$cat == "tv"],
        resp_param$disp[resp_param$cat == "tv"],
        resp_param$rate[resp_param$cat == "tv"]
      )
      means <- c(means, tv_means)
    }

    if (f == "pb") {
      pb_means <- mod_gompertz(
        trial_data$tpb,
        resp_param$max[resp_param$cat == "pb"],
        resp_param$disp[resp_param$cat == "pb"],
        resp_param$rate[resp_param$cat == "pb"]
      ) * trial_design$e
      means <- c(means, pb_means)
    }

    if (f == "br") {
      br_means <- mod_gompertz(
        trial_data$tod,
        resp_param$max[resp_param$cat == "br"],
        resp_param$disp[resp_param$cat == "br"],
        resp_param$rate[resp_param$cat == "br"]
      )
      br_means <- apply_carryover_to_component(
        br_means, trial_data, halflife
      )
      means <- c(means, br_means)
    }
  }

  if (is.null(cached_sigma)) {
    was_pd_corrected <- sigma_result$was_pd_corrected

    if (make_positive_definite && was_pd_corrected) {
      sigma <- corpcor::make.positive.definite(sigma, tol = 1e-3)
    }
  }

  chol_sigma <- tryCatch(chol(sigma), error = function(e) {
    chol(corpcor::make.positive.definite(sigma, tol = 1e-3))
  })

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

  for (timepoint_idx in 1:length(timepoint_names)) {
    timepoint_name <- timepoint_names[timepoint_idx]
    delta_col <- paste0("D_", timepoint_name)
    components <- paste(timepoint_name, factors, sep = ".")

    participant_data <- participant_data |>
      mutate(!!delta_col := rowSums(pick(all_of(components))))
  }

  for (timepoint_idx in 1:length(timepoint_names)) {
    timepoint_name <- timepoint_names[timepoint_idx]
    components <- paste(timepoint_name, factors, sep = ".")

    participant_data <- participant_data |>
      mutate(!!timepoint_name := BL - rowSums(pick(all_of(components))))
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
  for (mv in model_vars) {
    datamerged <- datamerged |> dplyr::filter(!is.na(.data[[mv]]))
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
  cdt <- dat |> dplyr::select(all_of(delta_cols))

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
          dat_rep <- dat |> dplyr::filter(replicate == i_s)
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
              datc_rep <- datc |> dplyr::filter(replicate == i_s)
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
