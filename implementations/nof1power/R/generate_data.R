#' Build the correlation matrix for multivariate normal sampling
#'
#' Constructs a square correlation matrix with AR(1) within-factor
#' autocorrelation indexed by `t_wk_cumulative`, cross-factor
#' correlations at matched and lagged timepoints, and (under
#' Architecture B, `dgp_architecture = "mvn"`) a biomarker-bio-response
#' differential correlation that decays during off-drug periods. The
#' implementation is a direct port of the pmsimstats-ng tidyverse
#' reference.
#'
#' @param labels Character vector of variable labels.
#' @param trial_design Trial-design tibble (single path).
#' @param model_param Model-parameter tibble with correlation
#'   coefficients `c.tv`, `c.pb`, `c.br`, `c.cf1t`, `c.cfct`, `c.bm`.
#' @param num_timepoints Number of timepoints.
#' @param factors Character vector of factor abbreviations (default
#'   `c("tv", "pb", "br")`).
#' @param trial_data Augmented trial data from [prepare_trial_data()].
#' @param lambda_cor Correlation decay rate for off-drug BM-BR
#'   coupling.
#' @param carryover_form Decay form for the BM-BR correlation under
#'   Architecture B; see [carryover_decay()].
#' @param weibull_shape Weibull shape parameter; see
#'   [carryover_decay()].
#' @param dgp_architecture Either `"mvn"` or `"mean_moderation"`.
#' @return A square correlation matrix with dimensions matching
#'   `labels`.
#' @keywords internal
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
      autocorrelation <- model_param[[
        paste("c", current_factor, sep = ".")
      ]]
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
      idx1_same <- cbind(match(n1_same, labels),
                         match(n2_same, labels))
      idx2_same <- cbind(idx1_same[, 2], idx1_same[, 1])
      correlations[idx1_same] <- model_param$c.cf1t
      correlations[idx2_same] <- model_param$c.cf1t

      if (num_timepoints > 1) {
        n1_lag <- paste(trial_design$timepoint_name[upper_pairs$p1],
                        current_factor, sep = ".")
        n2_lag <- paste(trial_design$timepoint_name[upper_pairs$p2],
                        other_factor, sep = ".")
        lag_vals <- model_param$c.cfct * ar_vals
        idx1_lag <- cbind(match(n1_lag, labels),
                          match(n2_lag, labels))
        idx2_lag <- cbind(idx1_lag[, 2], idx1_lag[, 1])
        correlations[idx1_lag] <- lag_vals
        correlations[idx2_lag] <- lag_vals
      }
    }

    if (current_factor == "br" && dgp_architecture == "mvn") {
      for (timepoint_idx in seq_len(num_timepoints)) {
        name1 <- paste(trial_design$timepoint_name[timepoint_idx],
                       "br", sep = ".")
        on_drug_now <- trial_data$on_drug[timepoint_idx]

        if (on_drug_now) {
          correlations["bm", name1] <- model_param$c.bm
          correlations[name1, "bm"] <- model_param$c.bm
        } else {
          tsd_now <- trial_data$tsd[timepoint_idx]
          if (tsd_now > 0 && lambda_cor > 0) {
            halflife_for_decay <- log(2) / lambda_cor
            decay <- carryover_decay(
              tsd_now, halflife_for_decay,
              form = carryover_form, shape = weibull_shape
            )
            correlations["bm", name1] <- model_param$c.bm * decay
            correlations[name1, "bm"] <- model_param$c.bm * decay
          }
        }
      }
    }
  }

  correlations
}

#' Compute the mean vector for the MVN draw
#'
#' @param baseline_param Baseline-parameter tibble with `cat`, `m`,
#'   `sd`. `cat` values expected: `"bm"`, `"BL"`.
#' @param resp_param Response-parameter tibble with `cat`, `max`,
#'   `disp`, `rate`, `sd`. `cat` values expected: `"tv"`, `"pb"`,
#'   `"br"`.
#' @param trial_data Augmented trial data from [prepare_trial_data()].
#' @param trial_design Single-path trial-design tibble.
#' @param factors Character vector of factor abbreviations.
#' @param halflife Carryover half-life.
#' @param carryover_form Decay form for carryover.
#' @param weibull_shape Weibull shape parameter.
#' @return Numeric vector of length 2 + 3 * num_timepoints.
#' @keywords internal
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

#' Build Sigma (covariance) matrix for MVN sampling
#'
#' Composes [build_correlation_matrix()] with the standard-deviation
#' vector and performs positive-definiteness checking, correcting
#' non-PD matrices via `corpcor::make.positive.definite(tol = 1e-3)`.
#' Returns the Cholesky factor for efficient sampling.
#'
#' @param model_param Model-parameter tibble.
#' @param resp_param Response-parameter tibble.
#' @param baseline_param Baseline-parameter tibble.
#' @param trial_design Single-path trial-design tibble.
#' @param verbose Logical; emit diagnostic messages.
#' @param lambda_cor Correlation decay rate for off-drug BM-BR
#'   coupling. `NA` auto-derives from `carryover_t1half`.
#' @param dgp_architecture Either `"mvn"` or `"mean_moderation"`.
#' @param ... Reserved for future use.
#' @return A list with elements `sigma`, `labels`,
#'   `standard_deviations`, `correlations`, `means`, `nP`, `cl`,
#'   `chol_sigma`, and `was_pd_corrected`.
#' @export
build_sigma_matrix <- function(model_param, resp_param, baseline_param,
                               trial_design,
                               verbose = FALSE,
                               lambda_cor = NA,
                               dgp_architecture = "mvn",
                               ...) {

  dgp_architecture <- match.arg(dgp_architecture,
                                c("mvn", "mean_moderation"))

  if (!is.null(model_param$carryover_t1half)) {
    stopifnot(
      "`carryover_t1half` must be non-negative" =
        is.na(model_param$carryover_t1half) ||
          model_param$carryover_t1half >= 0
    )
  }

  factors <- c("tv", "pb", "br")

  trial_design <- ensure_timepoint_name(trial_design)
  trial_data <- prepare_trial_data(trial_design)
  num_timepoints <- nrow(trial_design)

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
    rep(resp_param$sd[resp_param$cat == "pb"], num_timepoints) *
      trial_design$e,
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

  sigma <- outer(standard_deviations, standard_deviations) *
    correlations

  is_pd <- corpcor::is.positive.definite(sigma)

  if (verbose) {
    cat("Sigma matrix positive definite:", is_pd, "\n")
  }

  if (!is_pd) {
    if (verbose) {
      eigenvalues <- eigen(sigma, only.values = TRUE)$values
      neg_evals <- eigenvalues[eigenvalues < 0]
      warning(sprintf(
        paste0("Non-PD covariance matrix corrected (min eigenvalue: ",
               "%.4f, %d negative)"),
        min(eigenvalues), length(neg_evals)
      ))
    }
    sigma <- corpcor::make.positive.definite(sigma, tol = 1e-3)
  }

  chol_sigma <- tryCatch(
    chol(sigma),
    error = function(e) {
      chol(corpcor::make.positive.definite(sigma, tol = 1e-3))
    }
  )

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

#' Append delta and timepoint-symptom columns to participant data
#'
#' Adds `ptID`, per-timepoint delta columns `D_<tp>`, and the
#' symptom-scale columns `<tp>` computed as `BL - D_<tp>`.
#'
#' @param participant_data Participant data tibble from the MVN draw.
#' @param timepoint_names Character vector of timepoint names.
#' @param factors Factor abbreviations (`c("tv","pb","br")`).
#' @param N Number of participants.
#' @return A tibble augmented with `ptID`, `D_<tp>`, and `<tp>` cols.
#' @keywords internal
process_participant_data <- function(participant_data, timepoint_names,
                                     factors, N) {
  participant_data <- participant_data |>
    dplyr::mutate(ptID = seq_len(N))

  for (timepoint_idx in seq_along(timepoint_names)) {
    tp <- timepoint_names[timepoint_idx]
    components <- paste(tp, factors, sep = ".")
    delta_col <- paste0("D_", tp)
    participant_data <- participant_data |>
      dplyr::mutate(
        !!delta_col := rowSums(dplyr::pick(dplyr::all_of(components))),
        !!tp := .data$BL - .data[[delta_col]]
      )
  }

  participant_data
}

#' Generate simulated N-of-1 trial data
#'
#' Core data-generating function for the simulation pipeline. Samples
#' participant data from a multivariate normal with covariance
#' constructed by [build_sigma_matrix()]. Under Architecture A
#' (`dgp_architecture = "mean_moderation"`), applies an additive
#' biomarker-scaled shift to the bio-response columns at on-drug
#' timepoints. Under Architecture B (`"mvn"`, default), the
#' biomarker-response interaction is encoded via differential
#' correlation in the covariance matrix.
#'
#' @param model_param Model-parameter tibble (single row) with
#'   elements including `N`, `c.bm`, `c.tv`, `c.pb`, `c.br`,
#'   `c.cf1t`, `c.cfct`, `carryover_t1half`.
#' @param resp_param Response-parameter tibble with columns `cat`
#'   (in `{"tv","pb","br"}`), `max`, `disp`, `rate`, `sd`.
#' @param baseline_param Baseline-parameter tibble with columns `cat`
#'   (in `{"bm","BL"}`), `m`, `sd`.
#' @param trial_design Single-path trial-design tibble with columns
#'   `timepoint_name` (or `timeptnames`), `t_wk`, `e`, `tod`, `tsd`,
#'   `tpb`.
#' @param empirical Logical; `TRUE` forces the sample correlation to
#'   exactly match the target. Typically `FALSE`.
#' @param make_positive_definite Logical; force PD correction when
#'   covariance is non-PD. Usually `TRUE`.
#' @param seed Randomization seed or `NA`.
#' @param lambda_cor BM-BR correlation decay rate. `NA` auto-derives.
#' @param verbose Logical; emit diagnostic messages.
#' @param track_pd_stats Logical; attach PD-correction attributes to
#'   the returned tibble (`sigma_count`,
#'   `non_positive_definite_count`, `non_positive_definite_rate`).
#' @param cached_sigma Pre-built Sigma list from
#'   [build_sigma_matrix()] to skip construction.
#' @param dgp_architecture Either `"mvn"` (Architecture B, default)
#'   or `"mean_moderation"` (Architecture A).
#' @return A tibble with `bm`, `BL`, one column per factor-timepoint
#'   pair (`<tp>.tv`, `<tp>.pb`, `<tp>.br`), `ptID`, delta columns
#'   `D_<tp>`, and symptom columns `<tp>`.
#' @export
generate_data <- function(
    model_param, resp_param, baseline_param, trial_design,
    empirical, make_positive_definite, seed = NA,
    lambda_cor = NA, verbose = FALSE, track_pd_stats = TRUE,
    cached_sigma = NULL, dgp_architecture = "mvn") {

  dgp_architecture <- match.arg(dgp_architecture,
                                c("mvn", "mean_moderation"))

  if (!is.null(model_param$carryover_t1half)) {
    stopifnot(
      "`carryover_t1half` must be non-negative" =
        is.na(model_param$carryover_t1half) ||
          model_param$carryover_t1half >= 0
    )
  }

  lambda_cor <- resolve_lambda_cor(lambda_cor, model_param)

  factors <- c("tv", "pb", "br")

  trial_design <- ensure_timepoint_name(trial_design)
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

  if (!is.na(seed)) set.seed(seed)

  n <- model_param$N
  p <- length(means)

  if (isTRUE(empirical)) {
    sigma_mat <- if (!is.null(cached_sigma)) cached_sigma$sigma else
      sigma_result$sigma
    participant_matrix <- MASS::mvrnorm(
      n = n, mu = means, Sigma = sigma_mat, empirical = TRUE
    )
  } else {
    Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    participant_matrix <- Z %*% chol_sigma +
      matrix(means, nrow = n, ncol = p, byrow = TRUE)
  }

  participant_data <- tibble::as_tibble(participant_matrix,
                                       .name_repair = "minimal")
  colnames(participant_data) <- labels

  if (dgp_architecture == "mean_moderation") {
    beta_bm <- model_param$c.bm
    bm_mean <- baseline_param$m[baseline_param$cat == "bm"]
    bm_sd <- baseline_param$sd[baseline_param$cat == "bm"]
    br_sd <- resp_param$sd[resp_param$cat == "br"]
    bm_z <- (participant_data$bm - bm_mean) / bm_sd

    for (tp in seq_len(nrow(trial_design))) {
      if (trial_data$on_drug[tp]) {
        br_col <- paste(trial_design$timepoint_name[tp],
                        "br", sep = ".")
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
    was_corrected <- if (!is.null(cached_sigma)) FALSE else
      was_pd_corrected
    attr(participant_data, "sigma_count") <- 1L
    attr(participant_data, "non_positive_definite_count") <-
      as.integer(was_corrected)
    attr(participant_data, "non_positive_definite_rate") <-
      as.numeric(was_corrected)
  }

  participant_data
}
