#' Residual drug-effect decay multiplier
#'
#' Returns the fraction of residual drug effect present at
#' time-since-discontinuation `tsd` for a carryover with given
#' `halflife` and functional form. The multiplier equals 1 at
#' `tsd = 0` and 0.5 at `tsd = halflife` for all supported forms.
#'
#' @param tsd Time since drug discontinuation (same time units as
#'   `halflife`). May be a vector.
#' @param halflife Carryover half-life. The value of `tsd` at which
#'   the multiplier equals 0.5.
#' @param form One of `"exponential"` (default), `"linear"`,
#'   `"weibull"`, or `"power"`.
#' @param shape Weibull shape parameter `k`. Ignored unless
#'   `form = "weibull"`. `k = 1` recovers the exponential form;
#'   `k < 1` gives a heavier tail; `k > 1` gives accelerated washout.
#' @param power Power exponent. Ignored unless `form = "power"`.
#' @return Numeric vector of the same length as `tsd`, with values
#'   in `[0, 1]`.
#'
#' @details
#' - Exponential: `phi(t) = exp(-lambda t)` with
#'   `lambda = ln(2) / halflife`.
#' - Linear: `phi(t) = max(0, 1 - t / (2 * halflife))`.
#' - Weibull: `phi(t) = exp(-(lambda_w t)^k)` with
#'   `lambda_w = (ln 2)^{1/k} / halflife`.
#' - Power: `phi(t) = max(0, (1 - t / max_time)^power)` with
#'   `max_time = 3 * halflife`; supplied as a nof1_power extension.
#'
#' @examples
#' carryover_decay(c(0, 0.5, 1, 2), halflife = 1)
#' carryover_decay(c(0, 0.5, 1, 2), halflife = 1, form = "weibull",
#'                 shape = 2)
#' @export
carryover_decay <- function(tsd, halflife,
                            form = c("exponential", "linear",
                                     "weibull", "power"),
                            shape = 1, power = 1.5) {
  form <- match.arg(form)
  if (is.null(halflife) || halflife == 0) {
    return(rep(0, length(tsd)))
  }
  switch(
    form,
    exponential = exp(-log(2) / halflife * tsd),
    linear      = pmax(0, 1 - tsd / (2 * halflife)),
    weibull     = {
      if (!is.numeric(shape) || shape <= 0) {
        stop("Weibull `shape` must be a positive number.")
      }
      lambda_w <- (log(2))^(1 / shape) / halflife
      exp(-(lambda_w * tsd)^shape)
    },
    power = {
      if (!is.numeric(power) || power <= 0) {
        stop("Power exponent must be a positive number.")
      }
      max_time <- 3 * halflife
      pmax(0, (1 - tsd / max_time)^power)
    }
  )
}

#' Apply carryover decay to a bio-response mean trajectory
#'
#' Hendrickson-style carryover: at each off-drug timepoint, add the
#' previous timepoint's bio-response mean scaled by the decay
#' multiplier from [carryover_decay()]. Carryover is applied only to
#' the bio-response (br) component and only at timepoints with
#' `tsd > 0` and `on_drug == FALSE`.
#'
#' @param component_means Numeric vector of component means per
#'   timepoint (bio-response component).
#' @param trial_data Trial design tibble augmented by
#'   [prepare_trial_data()] with `on_drug` and `tsd` columns.
#' @param halflife Carryover half-life (weeks). `NULL` or `0` leaves
#'   `component_means` unchanged.
#' @param form Decay form passed to [carryover_decay()].
#' @param shape Weibull shape parameter; see [carryover_decay()].
#' @param power Power exponent; see [carryover_decay()].
#' @return Numeric vector of adjusted component means.
#' @export
apply_carryover_to_component <- function(
    component_means, trial_data, halflife,
    form = "exponential", shape = 1, power = 1.5) {

  num_timepoints <- length(component_means)

  if (is.null(halflife) || halflife == 0 || num_timepoints <= 1) {
    return(component_means)
  }

  carryover_indices <- which(!trial_data$on_drug & trial_data$tsd > 0)

  if (length(carryover_indices) > 0) {
    for (idx in carryover_indices) {
      prev_idx <- idx - 1
      if (prev_idx < 1) next
      time_lag <- trial_data$tsd[idx]
      decay_factor <- carryover_decay(
        time_lag, halflife, form = form,
        shape = shape, power = power
      )
      component_means[idx] <- component_means[idx] +
        component_means[prev_idx] * decay_factor
    }
  }

  component_means
}

#' Standalone carryover calculation from prior effect and decay form
#'
#' Convenience wrapper around [carryover_decay()] that scales a prior
#' effect by the decay multiplier. Retained for sensitivity-analysis
#' code that operates on scalar prior effects rather than full mean
#' trajectories.
#'
#' @param time_since_discontinuation Numeric, time since discontinuation.
#' @param previous_effect Numeric, prior effect to be scaled.
#' @param model Decay form; see [carryover_decay()].
#' @param params Named list of model parameters with elements
#'   `carryover_t1half`, optional `scale_factor`, and optional
#'   `k` (Weibull) or `power`.
#' @return Numeric vector of decayed effects.
#' @examples
#' calculate_carryover(2, 10, "exponential",
#'                     list(carryover_t1half = 1, scale_factor = 2))
#' calculate_carryover(2, 10, "weibull",
#'                     list(carryover_t1half = 1, k = 1.5))
#' @export
calculate_carryover <- function(time_since_discontinuation,
                                previous_effect,
                                model = "exponential", params) {
  halflife <- params$carryover_t1half %||% 1
  scale <- params$scale_factor %||% 1
  effective_tsd <- scale * time_since_discontinuation
  shape <- params$k %||% 1
  power <- params$power %||% 1.5

  decay <- carryover_decay(
    effective_tsd, halflife, form = model,
    shape = shape, power = power
  )
  previous_effect * decay
}

#' Null-coalescing operator
#'
#' Returns `a` unless `a` is `NULL`, in which case it returns `b`.
#' Used internally by the simulation pipeline for default parameters.
#'
#' @param a First value.
#' @param b Default value if `a` is `NULL`.
#' @return `a` if not `NULL`, otherwise `b`.
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
