#' Cumulative sum of a vector
#'
#' Thin wrapper around [cumsum()] used throughout the simulation
#' pipeline to emphasize that time grids and drug-on-time variables
#' are accumulated rather than differenced.
#'
#' @param values A numeric vector.
#' @return A numeric vector of the same length containing the
#'   cumulative sum.
#' @examples
#' cumulative(c(1, 2, 3, 4))
#' @export
cumulative <- function(values) {
  cumsum(values)
}

#' Modified Gompertz response curve
#'
#' Origin-passing Gompertz curve used to describe each latent factor
#' (time-variant, pharm-biomarker, bio-response) as a function of
#' cumulative time in the relevant phase. The implementation matches
#' the corrected Hendrickson formulation adopted in the pmsimstats-ng
#' tidyverse reference.
#'
#' @param t Numeric time vector (typically `t_wk_cumulative`, `tpb`,
#'   or `tod`).
#' @param maxr Asymptotic maximum of the response curve.
#' @param disp Displacement parameter controlling curve shape.
#' @param rate Rate parameter controlling curve steepness.
#' @return Numeric vector of the same length as `t`. `y(0) = 0` and
#'   `y(t) -> maxr` as `t -> Inf`.
#'
#' @details The curve is constructed from the base Gompertz form
#' `maxr * exp(-disp * exp(-rate * t))` with a vertical offset
#' subtraction so the curve passes through the origin, then rescaled
#' to preserve the asymptote `maxr`.
#' @examples
#' mod_gompertz(seq(0, 10, by = 1), maxr = 10, disp = 1, rate = 0.5)
#' @export
mod_gompertz <- function(t, maxr, disp, rate) {
  if (isTRUE(all.equal(maxr, 0))) {
    return(rep(0, length(t)))
  }
  y <- maxr * exp(-disp * exp(-rate * t))
  vert_offset <- maxr * exp(-disp * exp(-rate * 0))
  y <- y - vert_offset
  y * (maxr / (maxr - vert_offset))
}

#' Backward-compatibility wrapper for the modified Gompertz curve
#'
#' @param t Numeric time vector.
#' @param max Asymptotic maximum (maps to `maxr`).
#' @param disp Displacement parameter.
#' @param rate Rate parameter.
#' @return Numeric vector; see [mod_gompertz()].
#' @export
modgompertz <- function(t, max, disp, rate) {
  mod_gompertz(t = t, maxr = max, disp = disp, rate = rate)
}
