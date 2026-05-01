#' Characterize carryover from observed trial data
#'
#' Fits the analysis model across a range of candidate carryover
#' half-lives and compares model fit to identify the best-fitting
#' half-life. Also extracts the main drug effect and the
#' biomarker-drug interaction at each half-life.
#'
#' @param trialdesign_set Trial design paths (from buildtrialdesign)
#' @param dat Data as produced by generateData or actual trial data
#' @param half_lives Numeric vector of candidate half-lives to test
#'   (in weeks). Default: c(0.1, 0.25, 0.5, 1, 2, 4, 8)
#' @param op Analysis options (same as lme_analysis). If missing,
#'   defaults are used.
#' @return A list with:
#'   \itemize{
#'     \item{\code{comparison}}{  Data.table with one row per
#'       half-life: t_half, AIC, BIC, logLik, drug_beta, drug_se,
#'       drug_p, interaction_beta, interaction_se, interaction_p,
#'       ar1_phi}
#'     \item{\code{best_t_half}}{  The half-life with the lowest AIC}
#'     \item{\code{best_model}}{  The full model fit at the best
#'       half-life}
#'     \item{\code{best_data}}{  The analysis data at the best
#'       half-life}
#'   }
#' @export
characterize_carryover <- function(trialdesign_set, dat,
                                    half_lives = c(0.1, 0.25, 0.5,
                                                   1, 2, 4, 8),
                                    op) {
  if (missing(op)) {
    op <- list(useDE = FALSE, t_random_slope = FALSE,
               full_model_out = TRUE, simplecarryover = FALSE,
               carryover_t1half = 0, carryover_scalefactor = 1)
  } else {
    op$full_model_out <- TRUE
  }

  results <- vector("list", length(half_lives))

  for (i in seq_along(half_lives)) {
    hl <- half_lives[i]
    op_i <- op
    op_i$carryover_t1half <- hl
    op_i$carryover_scalefactor <- 1

    fit_result <- tryCatch(
      lme_analysis(trialdesign_set, dat, op_i),
      error = function(e) NULL
    )

    if (is.null(fit_result) || is.null(fit_result$fit)) {
      results[[i]] <- data.table(
        t_half = hl,
        aic = NA_real_, bic = NA_real_, loglik = NA_real_,
        drug_beta = NA_real_, drug_se = NA_real_,
        drug_p = NA_real_,
        interaction_beta = NA_real_, interaction_se = NA_real_,
        interaction_p = NA_real_,
        ar1_phi = NA_real_,
        converged = FALSE
      )
      next
    }

    fit <- fit_result$fit
    s <- summary(fit)$tTable
    coefnames <- rownames(s)

    # Extract main drug effect (Dbc)
    drug_target <- intersect(c("Dbc", "drug_binary"), coefnames)
    if (length(drug_target) > 0) {
      drug_target <- drug_target[1]
      drug_beta <- s[drug_target, "Value"]
      drug_se <- s[drug_target, "Std.Error"]
      drug_p <- s[drug_target, "p-value"]
    } else {
      drug_beta <- drug_se <- drug_p <- NA_real_
    }

    # Extract interaction (bm:Dbc)
    int_target <- intersect(c("bm:Dbc", "Dbc:bm",
                              "biomarker:drug_binary",
                              "drug_binary:biomarker"),
                            coefnames)
    if (length(int_target) > 0) {
      int_target <- int_target[1]
      int_beta <- s[int_target, "Value"]
      int_se <- s[int_target, "Std.Error"]
      int_p <- s[int_target, "p-value"]
    } else {
      int_beta <- int_se <- int_p <- NA_real_
    }

    # Extract AR(1) parameter if available
    ar1_phi <- tryCatch({
      coef(fit$modelStruct$corStruct, unconstrained = FALSE)
    }, error = function(e) NA_real_)

    results[[i]] <- data.table(
      t_half = hl,
      aic = AIC(fit),
      bic = BIC(fit),
      loglik = as.numeric(logLik(fit)),
      drug_beta = drug_beta,
      drug_se = drug_se,
      drug_p = drug_p,
      interaction_beta = int_beta,
      interaction_se = int_se,
      interaction_p = int_p,
      ar1_phi = as.numeric(ar1_phi),
      converged = TRUE
    )
  }

  comparison <- rbindlist(results)
  best_idx <- which.min(comparison$aic)
  best_hl <- comparison$t_half[best_idx]

  # Refit best model to return full output
  op_best <- op
  op_best$carryover_t1half <- best_hl
  op_best$carryover_scalefactor <- 1
  best_fit <- lme_analysis(trialdesign_set, dat, op_best)

  list(
    comparison = comparison,
    best_t_half = best_hl,
    best_model = best_fit$fit,
    best_data = best_fit$datamerged
  )
}


#' Analyze trial for main drug effect and biomarker interaction
#'
#' Extended version of lme_analysis that extracts both the main
#' drug effect (does the drug work?) and the biomarker-drug
#' interaction (does the biomarker predict who responds?).
#' Optionally computes predicted response curves and clinical
#' decision thresholds.
#'
#' @param trialdesign_set Trial design paths
#' @param dat Data as produced by generateData or actual trial data
#' @param op Analysis options (same as lme_analysis)
#' @param threshold Clinical significance threshold for the drug
#'   effect (e.g., 5 CAPS points). Used to compute the biomarker
#'   value at which predicted benefit exceeds this threshold.
#'   Default: NULL (no threshold computation).
#' @return A list with:
#'   \itemize{
#'     \item{\code{drug_effect}}{  Data.table with beta, SE, p-value
#'       for the main drug effect (Dbc or bm:t)}
#'     \item{\code{interaction}}{  Data.table with beta, SE, p-value
#'       for the biomarker-drug interaction}
#'     \item{\code{coefficients}}{  Full coefficient table from the
#'       model}
#'     \item{\code{variance_components}}{  Random intercept variance,
#'       residual variance, AR(1) parameter}
#'     \item{\code{predicted_effect}}{  Function that takes a
#'       biomarker value and returns the predicted drug effect}
#'     \item{\code{threshold_bm}}{  Biomarker value at which the
#'       predicted drug effect equals the clinical threshold
#'       (if threshold is provided)}
#'     \item{\code{model}}{  The full nlme model object}
#'     \item{\code{data}}{  The analysis-ready data}
#'   }
#' @export
analyze_trial_extended <- function(trialdesign_set, dat, op,
                                    threshold = NULL) {
  if (missing(op)) {
    op <- list(useDE = FALSE, t_random_slope = FALSE,
               full_model_out = TRUE, simplecarryover = FALSE,
               carryover_t1half = 0, carryover_scalefactor = 1)
  } else {
    op$full_model_out <- TRUE
  }

  fit_result <- lme_analysis(trialdesign_set, dat, op)

  if (is.null(fit_result$fit)) {
    return(list(
      drug_effect = data.table(beta = NA, se = NA, p = NA),
      interaction = data.table(beta = NA, se = NA, p = NA),
      coefficients = NULL,
      variance_components = NULL,
      predicted_effect = NULL,
      threshold_bm = NA,
      model = NULL,
      data = NULL
    ))
  }

  fit <- fit_result$fit
  s <- summary(fit)$tTable
  coefnames <- rownames(s)

  # Main drug effect
  drug_target <- intersect(c("Dbc", "drug_binary"), coefnames)
  if (length(drug_target) > 0) {
    dt <- drug_target[1]
    drug_effect <- data.table(
      beta = s[dt, "Value"],
      se = s[dt, "Std.Error"],
      p = s[dt, "p-value"],
      ci_lower = s[dt, "Value"] - 1.96 * s[dt, "Std.Error"],
      ci_upper = s[dt, "Value"] + 1.96 * s[dt, "Std.Error"]
    )
  } else {
    drug_effect <- data.table(
      beta = NA, se = NA, p = NA,
      ci_lower = NA, ci_upper = NA)
  }

  # Biomarker-drug interaction
  int_target <- intersect(c("bm:Dbc", "Dbc:bm", "bm:t", "t:bm",
                            "biomarker:drug_binary"),
                          coefnames)
  if (length(int_target) > 0) {
    it <- int_target[1]
    interaction <- data.table(
      beta = s[it, "Value"],
      se = s[it, "Std.Error"],
      p = s[it, "p-value"],
      ci_lower = s[it, "Value"] - 1.96 * s[it, "Std.Error"],
      ci_upper = s[it, "Value"] + 1.96 * s[it, "Std.Error"]
    )
  } else {
    interaction <- data.table(
      beta = NA, se = NA, p = NA,
      ci_lower = NA, ci_upper = NA)
  }

  # Variance components
  vc <- tryCatch({
    v <- VarCorr(fit)
    ar1_phi <- tryCatch(
      as.numeric(coef(fit$modelStruct$corStruct,
                      unconstrained = FALSE)),
      error = function(e) NA_real_)
    data.table(
      random_intercept_sd = as.numeric(v[1, "StdDev"]),
      residual_sd = as.numeric(v[2, "StdDev"]),
      ar1_phi = ar1_phi
    )
  }, error = function(e) {
    data.table(random_intercept_sd = NA, residual_sd = NA,
               ar1_phi = NA)
  })

  # Biomarker mean from data
  bm_mean <- mean(fit_result$datamerged$bm, na.rm = TRUE)

  # Predicted drug effect as function of biomarker
  # effect(bm) = beta_Dbc + beta_interaction * (bm - bm_mean)
  if (!is.na(drug_effect$beta) && !is.na(interaction$beta)) {
    predicted_effect <- function(bm) {
      drug_effect$beta + interaction$beta * (bm - bm_mean)
    }
  } else {
    predicted_effect <- NULL
  }

  # Clinical threshold: at what biomarker value does
  # predicted benefit exceed the threshold?
  threshold_bm <- NA_real_
  if (!is.null(threshold) && !is.null(predicted_effect) &&
      !is.na(interaction$beta) && interaction$beta != 0) {
    threshold_bm <- bm_mean +
      (threshold - drug_effect$beta) / interaction$beta
  }

  list(
    drug_effect = drug_effect,
    interaction = interaction,
    coefficients = as.data.table(s, keep.rownames = "term"),
    variance_components = vc,
    predicted_effect = predicted_effect,
    biomarker_mean = bm_mean,
    threshold_bm = threshold_bm,
    threshold_value = threshold,
    model = fit,
    data = fit_result$datamerged
  )
}


#' Print summary of carryover characterization
#'
#' @param result Output from characterize_carryover()
#' @export
print_carryover_summary <- function(result) {
  cat("=== Carryover Characterization ===\n\n")

  comp <- result$comparison[converged == TRUE]
  if (nrow(comp) == 0) {
    cat("No models converged.\n")
    return(invisible(NULL))
  }

  cat(sprintf("%-8s  %8s  %8s  %10s  %8s  %8s\n",
    "t_half", "AIC", "BIC", "Drug p",
    "Int p", "AR1 phi"))
  cat(strrep("-", 60), "\n")

  for (i in 1:nrow(comp)) {
    r <- comp[i]
    best <- ifelse(r$t_half == result$best_t_half, " <--", "")
    cat(sprintf("%6.2f wk  %8.1f  %8.1f  %10.4f  %8.4f  %8.3f%s\n",
      r$t_half, r$aic, r$bic,
      r$drug_p, r$interaction_p, r$ar1_phi, best))
  }

  cat(sprintf("\nBest half-life (by AIC): %.2f weeks\n",
    result$best_t_half))
  cat(sprintf("Drug effect at best: beta = %.3f, p = %.4f\n",
    comp[t_half == result$best_t_half]$drug_beta,
    comp[t_half == result$best_t_half]$drug_p))
  cat(sprintf("Interaction at best: beta = %.4f, p = %.4f\n",
    comp[t_half == result$best_t_half]$interaction_beta,
    comp[t_half == result$best_t_half]$interaction_p))
}


#' Print summary of extended trial analysis
#'
#' @param result Output from analyze_trial_extended()
#' @export
print_trial_summary <- function(result) {
  cat("=== Trial Analysis Summary ===\n\n")

  cat("Drug Effect (does the drug work?):\n")
  cat(sprintf("  beta = %.3f (95%% CI: %.3f to %.3f)\n",
    result$drug_effect$beta,
    result$drug_effect$ci_lower,
    result$drug_effect$ci_upper))
  cat(sprintf("  p = %.4f %s\n\n",
    result$drug_effect$p,
    ifelse(result$drug_effect$p < 0.05,
           "(significant)", "(not significant)")))

  cat("Biomarker-Drug Interaction (does biomarker predict response?):\n")
  cat(sprintf("  beta = %.4f (95%% CI: %.4f to %.4f)\n",
    result$interaction$beta,
    result$interaction$ci_lower,
    result$interaction$ci_upper))
  cat(sprintf("  p = %.4f %s\n\n",
    result$interaction$p,
    ifelse(result$interaction$p < 0.05,
           "(significant)", "(not significant)")))

  cat("Variance Components:\n")
  cat(sprintf("  Random intercept SD: %.2f\n",
    result$variance_components$random_intercept_sd))
  cat(sprintf("  Residual SD: %.2f\n",
    result$variance_components$residual_sd))
  cat(sprintf("  AR(1) phi: %.3f\n\n",
    result$variance_components$ar1_phi))

  if (!is.null(result$predicted_effect)) {
    bm <- result$biomarker_mean
    cat("Predicted Drug Effect by Biomarker:\n")
    cat(sprintf("  At mean BM (%.1f): %.2f CAPS points\n",
      bm, result$predicted_effect(bm)))
    cat(sprintf("  At BM +1 SD:       %.2f CAPS points\n",
      result$predicted_effect(bm + 15.4)))
    cat(sprintf("  At BM -1 SD:       %.2f CAPS points\n",
      result$predicted_effect(bm - 15.4)))
  }

  if (!is.na(result$threshold_bm) && !is.null(result$threshold_value)) {
    cat(sprintf("\nClinical Threshold (%.0f CAPS points):\n",
      result$threshold_value))
    cat(sprintf("  Biomarker value needed: %.1f mmHg\n",
      result$threshold_bm))
  }
}
