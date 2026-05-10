#' Construct a trial-design object as a list of per-path tibbles
#'
#' Returns a list `list(metadata, trialpaths)` where `trialpaths` is
#' a list of per-path tibbles with one row per timepoint and columns
#' `timeptnames`, `t_wk`, `e` (expectancy), `tod` (cumulative time on
#' drug), `tsd` (cumulative time since discontinuation), and `tpb`
#' (cumulative time with positive expectancy). Direct port of the
#' pmsimstats-ng tidyverse `build_trial_design`.
#'
#' @param name_longform Human-readable design name.
#' @param name_shortform Short tag for plots and filenames.
#' @param timepoints Numeric vector of visit weeks (cumulative).
#' @param timeptnames Character vector of visit labels. Default
#'   `paste0("V", seq_along(timepoints))`.
#' @param expectancies Numeric vector of per-visit expectancy values
#'   (same length as `timepoints`).
#' @param ondrug List of binary vectors (same length as `timepoints`)
#'   indicating on-drug status per path. One element per trial path.
#' @return A list with elements `metadata` (design name, timepoints,
#'   expectancies, ondrug) and `trialpaths` (a list of per-path
#'   tibbles).
#' @examples
#' td <- build_trial_design(
#'   name_longform = "Hybrid",
#'   name_shortform = "H",
#'   timepoints = c(2, 4, 6, 8),
#'   expectancies = c(0.5, 0.5, 0, 0),
#'   ondrug = list(c(1, 1, 0, 0), c(0, 0, 1, 1))
#' )
#' @export
build_trial_design <- function(
    name_longform = "Trial Design 1",
    name_shortform = name_longform,
    timepoints,
    timeptnames = paste0("V", seq_along(timepoints)),
    expectancies, ondrug) {

  t_wk <- c(timepoints[1], diff(timepoints))

  trialpaths <- purrr::map(ondrug, function(od) {
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

    tibble::tibble(
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

#' Process a single trial-design path
#'
#' Convenience helper that wraps [prepare_trial_data()] to attach
#' `t_wk_cumulative` and `on_drug` columns. Accepts either a flat
#' single-path tibble or a `list(metadata, trialpaths)` object (in
#' which case the first path is processed).
#'
#' @param trial_design Either a single-path tibble or a list returned
#'   by [build_trial_design()].
#' @return The augmented single-path tibble.
#' @export
process_trial_design <- function(trial_design) {
  if (is.list(trial_design) && "trialpaths" %in% names(trial_design)) {
    path <- trial_design$trialpaths[[1]]
  } else {
    path <- trial_design
  }
  prepare_trial_data(ensure_timepoint_name(path))
}
