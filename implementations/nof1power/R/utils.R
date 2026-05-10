#' Create AR(1) correlation matrix
#'
#' Convenience helper for a discrete-time AR(1) correlation matrix
#' used in diagnostic work. The simulation pipeline builds AR(1)
#' structure directly inside [build_correlation_matrix()] against the
#' trial-design time vector; this helper is retained for exploratory
#' analysis.
#'
#' @param n_times Number of timepoints.
#' @param rho Lag-one autocorrelation.
#' @return A square correlation matrix with AR(1) structure.
#' @export
create_ar1_corr <- function(n_times, rho) {
  matrix(rho^abs(outer(1:n_times, 1:n_times, "-")), n_times, n_times)
}

#' Extract working correlation matrix from a `gls`/`lme` fit
#'
#' @param glsob A `gls` or `lme` object with a structured correlation
#'   component.
#' @param cov Logical; if `TRUE`, return covariance (default).
#' @param ... Reserved for future use.
#' @return A correlation matrix.
#' @export
corandcov <- function(glsob, cov = TRUE, ...) {
  cc <- nlme::corMatrix(glsob$modelStruct$corStruct)[[5]]
  varstruct <- glsob$modelStruct$varStruct
  varests <- stats::coef(varstruct, uncons = FALSE, allCoef = TRUE)
  rownames(cc) <- names(varests)
  colnames(cc) <- names(varests)
  cc
}

#' Resolve lambda_cor from carryover half-life
#'
#' Auto-derives the BM-BR correlation decay rate `lambda_cor` from
#' the carryover half-life `carryover_t1half` when `lambda_cor` is
#' `NA`. Returns `0` when the half-life is missing or non-positive.
#'
#' @param lambda_cor Supplied decay rate, or `NA` to auto-derive.
#' @param model_param Model-parameter tibble or list with optional
#'   `carryover_t1half`.
#' @return Numeric decay rate (per week).
#' @keywords internal
resolve_lambda_cor <- function(lambda_cor, model_param) {
  if (!is.na(lambda_cor)) return(lambda_cor)
  if (!is.null(model_param$carryover_t1half) &&
      model_param$carryover_t1half > 0) {
    log(2) / model_param$carryover_t1half
  } else {
    0
  }
}

#' Ensure a trial-design tibble has a `timepoint_name` column
#'
#' Reconciles tidyverse (`timepoint_name`), pmsimstats short-form
#' (`timeptnames`), and hybrid (`timeptname`) column conventions
#' by adding `timepoint_name` if missing.
#'
#' @param td Trial-design data frame or tibble.
#' @return The same object with a `timepoint_name` column.
#' @keywords internal
ensure_timepoint_name <- function(td) {
  if (!"timepoint_name" %in% names(td)) {
    nm <- intersect(c("timeptnames", "timeptname"), names(td))
    if (length(nm)) td$timepoint_name <- td[[nm[1]]]
  }
  td
}

#' Ensure a trial-design tibble has a `timeptnames` column
#'
#' Converts `timepoint_name` or `timeptname` to `timeptnames` if
#' present; the pmsimstats short-form column name used by
#' [lme_analysis()].
#'
#' @param td Trial-design data frame or tibble.
#' @return The same object with a `timeptnames` column.
#' @keywords internal
ensure_timeptnames <- function(td) {
  nms <- names(td)
  if ("timeptnames" %in% nms) return(td)
  src <- intersect(c("timeptname", "timepoint_name"), nms)
  if (length(src)) {
    td$timeptnames <- td[[src[1]]]
    td[[src[1]]] <- NULL
  }
  td
}

#' Return the timepoint-name column of a trial-design tibble
#'
#' @param td Trial-design data frame or tibble.
#' @return Character vector of timepoint names.
#' @keywords internal
get_timepoint_names <- function(td) {
  nms <- names(td)
  if ("timeptnames" %in% nms) return(td$timeptnames)
  if ("timeptname" %in% nms) return(td$timeptname)
  td$timepoint_name
}

#' Augment a trial-design path with cumulative-time and on-drug columns
#'
#' @param trial_design A single-path trial design data frame (one row
#'   per timepoint).
#' @return A tibble with the same rows and added `t_wk_cumulative`
#'   and `on_drug` columns.
#' @keywords internal
prepare_trial_data <- function(trial_design) {
  if (!is.data.frame(trial_design)) {
    stop("Trial design is not a data frame. Class: ",
         class(trial_design), ", Type: ", typeof(trial_design))
  }

  trial_data <- tibble::as_tibble(trial_design)

  if (!("t_wk" %in% names(trial_data))) {
    stop("'t_wk' column not found in trial design data")
  }

  trial_data$t_wk_cumulative <- cumsum(trial_data$t_wk)
  trial_data$on_drug <- (trial_data$tod > 0)
  trial_data
}
