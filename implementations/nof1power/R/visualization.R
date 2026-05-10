#' Plot a power heatmap across a parameter grid
#'
#' Produces a `ggplot2` heatmap from the `results` tibble returned
#' by [generate_simulated_results()]. Rows of `results` are
#' per-replicate analysis outputs with columns `trialdesign`,
#' `respparamset`, `blparamset`, `modelparamset`, parameters from
#' `modelparams` (typically `N`, `c.bm`, `carryover_t1half`), and
#' replicate-level quantities `beta`, `betaSE`, `p`, `irep`.
#'
#' The caller specifies four varying parameters (`param2vary`),
#' parameters to hold fixed (`param2hold`), and which aggregate
#' quantity to plot (`param2plot`). The first two `param2vary`
#' entries form the facet grid; the last two form the tile x and y
#' axes. Aggregation is computed with `dplyr::group_by` over the
#' four varying parameters.
#'
#' @param data Either the full list returned by
#'   [generate_simulated_results()] (with `$results` and
#'   `$parameterselections`) or just the `results` tibble.
#' @param param2plot Aggregate quantity; one of `"power"`,
#'   `"mean_beta"`, `"mean_betaSE"`, `"sd_beta"`, or
#'   `"mean_frac_NA"`.
#' @param param2vary Character vector of length 4 naming the
#'   varying parameters (columns of the results tibble).
#' @param param2hold Named list of parameters to hold fixed.
#' @param alpha Significance threshold for the `"power"` metric.
#' @return A `ggplot` object.
#' @export
plot_modeling_results <- function(data, param2plot = "power",
                                  param2vary, param2hold = list(),
                                  alpha = 0.05) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_modeling_results()")
  }

  if (is.list(data) && "results" %in% names(data)) {
    d <- data$results
  } else {
    d <- data
  }

  if (length(param2vary) != 4) {
    stop("param2vary must be a character vector of length 4.")
  }

  missing_vary <- setdiff(param2vary, names(d))
  if (length(missing_vary)) {
    stop("Missing columns in results: ",
         paste(missing_vary, collapse = ", "))
  }

  for (nm in names(param2hold)) {
    if (nm %in% names(d)) {
      d <- d[d[[nm]] == param2hold[[nm]], , drop = FALSE]
    }
  }

  agg <- d |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(param2vary))
    ) |>
    dplyr::summarise(
      power = mean(.data$p < alpha, na.rm = TRUE),
      mean_beta = mean(.data$beta, na.rm = TRUE),
      sd_beta = stats::sd(.data$beta, na.rm = TRUE),
      mean_betaSE = mean(.data$betaSE, na.rm = TRUE),
      mean_frac_NA = mean(.data$frac_NA, na.rm = TRUE),
      .groups = "drop"
    )

  if (!(param2plot %in% names(agg))) {
    stop("Unknown param2plot: ", param2plot)
  }

  x_var <- param2vary[3]
  y_var <- param2vary[4]
  row_facet <- param2vary[2]
  col_facet <- param2vary[1]

  ggplot2::ggplot(
    agg,
    ggplot2::aes(
      x = factor(.data[[x_var]]),
      y = factor(.data[[y_var]]),
      fill = .data[[param2plot]]
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data[[param2plot]], 2)),
      size = 3
    ) +
    ggplot2::scale_fill_distiller(
      type = "seq", palette = "RdYlBu", name = param2plot
    ) +
    ggplot2::facet_grid(
      stats::reformulate(col_facet, row_facet)
    ) +
    ggplot2::labs(x = x_var, y = y_var) +
    ggplot2::theme_minimal()
}

#' Plot mean latent-factor trajectories for a trial design
#'
#' Takes the `trialpaths` list of a [build_trial_design()] result
#' and a parameter set for [build_sigma_matrix()], and plots the
#' mean trajectories of the three latent factors (`tv`, `pb`,
#' `br`) across timepoints for each path. The drug-on intervals are
#' overlaid at the bottom of each facet.
#'
#' @param trial_design The `list(metadata, trialpaths)` object
#'   returned by [build_trial_design()].
#' @param model_param Model-parameter tibble (single row).
#' @param resp_param Response-parameter tibble.
#' @param baseline_param Baseline-parameter tibble.
#' @param dgp_architecture Either `"mvn"` or `"mean_moderation"`.
#' @return A `ggplot` object.
#' @export
plot_factor_trajectories <- function(trial_design, model_param,
                                     resp_param, baseline_param,
                                     dgp_architecture = "mvn") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_factor_trajectories()")
  }

  if (!is.list(trial_design) ||
      !all(c("metadata", "trialpaths") %in% names(trial_design))) {
    stop("trial_design must be the output of build_trial_design().")
  }

  paths <- trial_design$trialpaths

  per_path <- purrr::map2(
    paths, seq_along(paths),
    function(path, path_idx) {
      path <- ensure_timepoint_name(path)
      sig <- build_sigma_matrix(
        model_param, resp_param, baseline_param, path,
        dgp_architecture = dgp_architecture
      )

      means <- sig$means
      labels <- sig$labels

      n_skip <- 2  # bm, BL come first
      tp_names <- path$timepoint_name
      nTP <- length(tp_names)

      tv_means <- means[(n_skip + 1):(n_skip + nTP)]
      pb_means <- means[(n_skip + nTP + 1):(n_skip + 2 * nTP)]
      br_means <- means[(n_skip + 2 * nTP + 1):(n_skip + 3 * nTP)]

      tibble::tibble(
        path = path_idx,
        timepoint = rep(tp_names, 3),
        week = rep(cumsum(path$t_wk), 3),
        factor = rep(c("tv", "pb", "br"), each = nTP),
        value = c(tv_means, pb_means, br_means)
      )
    }
  )

  traj <- dplyr::bind_rows(per_path)

  ondrug_segments <- purrr::map2(
    paths, seq_along(paths),
    function(path, path_idx) {
      starts <- c(0, cumsum(path$t_wk)[-nrow(path)])
      stops <- cumsum(path$t_wk)
      on_flag <- path$tod > 0
      idx <- which(on_flag)
      if (!length(idx)) return(NULL)
      tibble::tibble(
        path = path_idx,
        start = starts[idx],
        stop = stops[idx]
      )
    }
  ) |>
    dplyr::bind_rows()

  p <- ggplot2::ggplot(traj) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data$week, y = .data$value,
                   color = .data$factor)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$week, y = .data$value,
                   color = .data$factor),
      size = 1.2
    ) +
    ggplot2::facet_wrap(~ path, labeller = ggplot2::label_both) +
    ggplot2::scale_color_manual(
      name = "Factor",
      values = c(br = "red", pb = "blue", tv = "brown"),
      labels = c(br = "BR", pb = "ER", tv = "TR"),
      breaks = c("br", "pb", "tv")
    ) +
    ggplot2::labs(
      x = "Weeks",
      y = "Mean factor value",
      title = paste0(
        "Factor trajectories (",
        trial_design$metadata$name_shortform, ", ",
        dgp_architecture, ")"
      )
    ) +
    ggplot2::theme_minimal()

  if (nrow(ondrug_segments) > 0) {
    y_min <- min(traj$value) - 0.05 * diff(range(traj$value))
    p <- p + ggplot2::geom_segment(
      data = ondrug_segments,
      ggplot2::aes(x = .data$start, xend = .data$stop,
                   y = y_min, yend = y_min),
      linewidth = 1.2, color = "black"
    )
  }

  p
}

#' Plot a two-factor summary across a simulation parameter grid
#'
#' Lightweight tile plot that displays the mean of `param_to_plot`
#' across two varying parameters. A thin wrapper around
#' `ggplot2::geom_tile` that works directly on the results tibble.
#'
#' @param data Either the full list returned by
#'   [generate_simulated_results()] or just the `results` tibble.
#' @param param_to_plot One of `"power"`, `"mean_beta"`,
#'   `"mean_betaSE"`, `"sd_beta"`, or `"mean_frac_NA"`.
#' @param params_to_vary Character vector of length 2.
#' @param params_to_hold Named list of parameters to hold fixed.
#' @param alpha Significance threshold for the `"power"` metric.
#' @return A `ggplot` object.
#' @export
plot_parameter_space <- function(data, param_to_plot = "power",
                                 params_to_vary,
                                 params_to_hold = list(),
                                 alpha = 0.05) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_parameter_space()")
  }

  if (is.list(data) && "results" %in% names(data)) {
    d <- data$results
  } else {
    d <- data
  }

  if (length(params_to_vary) != 2) {
    stop("params_to_vary must be a character vector of length 2.")
  }

  for (nm in names(params_to_hold)) {
    if (nm %in% names(d)) {
      d <- d[d[[nm]] == params_to_hold[[nm]], , drop = FALSE]
    }
  }

  agg <- d |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(params_to_vary))
    ) |>
    dplyr::summarise(
      power = mean(.data$p < alpha, na.rm = TRUE),
      mean_beta = mean(.data$beta, na.rm = TRUE),
      sd_beta = stats::sd(.data$beta, na.rm = TRUE),
      mean_betaSE = mean(.data$betaSE, na.rm = TRUE),
      mean_frac_NA = mean(.data$frac_NA, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot(
    agg,
    ggplot2::aes(
      x = factor(.data[[params_to_vary[1]]]),
      y = factor(.data[[params_to_vary[2]]]),
      fill = .data[[param_to_plot]]
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = round(.data[[param_to_plot]], 2)),
      size = 3
    ) +
    ggplot2::scale_fill_distiller(
      type = "seq", palette = "RdYlBu", name = param_to_plot
    ) +
    ggplot2::labs(x = params_to_vary[1], y = params_to_vary[2]) +
    ggplot2::theme_minimal()
}
