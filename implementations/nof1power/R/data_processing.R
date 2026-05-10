#' Apply dropout-based censoring to simulated trial data
#'
#' Simulates participant dropout as a mixture of time-weighted
#' uniform censoring and symptom-trajectory-biased censoring. Once a
#' timepoint is censored for a participant, all subsequent timepoints
#' are also censored (carry-forward masking). Direct port of the
#' pmsimstats-ng tidyverse `censor_data`.
#'
#' @param dat Data tibble produced by [generate_data()] with symptom
#'   columns `<tp>` and delta columns `D_<tp>`.
#' @param trialdesign Single-path trial-design tibble with
#'   `timeptnames` (or `timepoint_name`), `t_wk`.
#' @param censorparam One-row data frame or list with numeric entries
#'   `beta0` (total censoring rate), `beta1` (fraction of censoring
#'   that is trajectory-biased), and `eb1` (exponent for the biased
#'   component).
#' @return The input `dat` with symptom cells set to `NA` at censored
#'   timepoints.
#' @examples
#' # See vignettes for integrated usage.
#' @export
censor_data <- function(dat, trialdesign, censorparam) {
  td <- tibble::as_tibble(trialdesign)
  tp_names <- get_timepoint_names(td)

  delta_cols <- paste0("D_", tp_names)
  cdt <- dat |> dplyr::select(dplyr::all_of(delta_cols))

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
  cdt_r1 <- stats::runif(
    nr * nc, min = 0,
    max = 2 * mean(cdt_p1) * (0.5 / fna1)
  )
  cdt_r2 <- stats::runif(
    nr * nc, min = 0,
    max = 2 * mean(cdt_p2) * (0.5 / fna2)
  )

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
