#' Linear mixed-effects analysis of a single simulated trial
#'
#' Fits the canonical pmsimstats analysis model via `nlme::lme` with
#' a `corCAR1(form = ~t | ptID)` residual correlation and a random
#' intercept per participant. The fixed-effects formula is assembled
#' adaptively based on three properties of the input data:
#'
#' - `var_in_Db`: whether the binary drug indicator varies within
#'   participant. If yes, the formula includes `Dbc + bm:Dbc`; if
#'   not, it falls back to `bm:t`.
#' - `var_in_exp`: whether the expectancy vector varies. If yes and
#'   `op$useDE = TRUE`, adds `+ De`.
#' - `var_in_tsd`: whether `tsd` ever varies. If yes and
#'   `op$simplecarryover = TRUE`, adds `+ tsd`.
#'
#' Direct port of the pmsimstats-ng tidyverse `lme_analysis`.
#'
#' @param trial_design_set A list of per-path tibbles with
#'   timepoint/time/expectancy/drug columns (the `trialpaths`
#'   element of a [build_trial_design()] result).
#' @param dat Simulated data tibble produced by [generate_data()].
#' @param op Options list with components:
#'   \itemize{
#'     \item `useDE` (default `TRUE`): include `+ De` when
#'       expectancy varies.
#'     \item `t_random_slope` (default `FALSE`): reserved for
#'       future extension.
#'     \item `full_model_out` (default `FALSE`): return the full
#'       list including the `nlme::lme` fit.
#'     \item `simplecarryover` (default `FALSE`): add `+ tsd` as a
#'       nuisance term when carryover is being modelled linearly.
#'     \item `carryover_t1half` (default `0`): carryover half-life
#'       used to compute the continuous drug indicator `Dbc`.
#'     \item `carryover_scalefactor` (default `1`): scale factor in
#'       the `Dbc` decay.
#'   }
#' @return A tibble with columns `beta`, `betaSE`, `p`, `issingular`,
#'   `warning`. If `op$full_model_out = TRUE`, returns a list with
#'   `form`, `fit`, `datamerged`, and the `stdout` tibble.
#' @export
lme_analysis <- function(trial_design_set, dat, op) {

  if (missing(op)) {
    op <- list(
      useDE = TRUE, t_random_slope = FALSE,
      full_model_out = FALSE, carryover_t1half = 0,
      simplecarryover = FALSE, carryover_scalefactor = 1,
      target_coef = "interaction"
    )
  } else if (!("simplecarryover" %in% names(op))) {
    op$simplecarryover <- FALSE
  }
  if (!("carryover_t1half" %in% names(op))) {
    op$carryover_t1half <- 0
    op$carryover_scalefactor <- 1
  }
  if (!("useDE" %in% names(op))) op$useDE <- TRUE
  if (!("full_model_out" %in% names(op))) op$full_model_out <- FALSE
  if (!("target_coef" %in% names(op))) op$target_coef <- "interaction"
  op$target_coef <- match.arg(
    op$target_coef, c("interaction", "main")
  )

  if ((op$carryover_t1half > 0) && (op$simplecarryover == TRUE)) {
    stop("Cannot use both simplecarryover and carryover_t1half")
  }

  stopifnot(
    "`carryover_t1half` must be non-negative" =
      op$carryover_t1half >= 0
  )

  n_groups <- length(trial_design_set)
  datout <- vector("list", n_groups)
  last_ptID <- 0

  for (g in seq_len(n_groups)) {
    td <- trial_design_set[[g]]
    td_names <- get_timepoint_names(td)

    dat_single <- dat |>
      dplyr::filter(.data$path == g) |>
      tibble::as_tibble()

    td_with_bl <- dplyr::bind_rows(
      tibble::tibble(
        timeptnames = "BL", t_wk = 0, e = 0,
        tod = 0, tsd = 0, tpb = 0
      ),
      tibble::as_tibble(ensure_timeptnames(td))
    ) |>
      dplyr::mutate(t = cumsum(.data$t_wk))

    all_names <- c("BL", td_names)
    valid_names <- all_names[all_names %in% names(dat_single)]

    data_wide <- dat_single |>
      dplyr::select(.data$ptID, .data$bm,
                    dplyr::all_of(valid_names)) |>
      dplyr::mutate(ptID = .data$ptID + last_ptID)
    last_ptID <- max(data_wide$ptID)

    data_long <- data_wide |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(valid_names),
        names_to = "timeptnames",
        values_to = "Sx",
        values_drop_na = FALSE
      ) |>
      dplyr::left_join(
        td_with_bl |>
          dplyr::select(.data$timeptnames, .data$t,
                        De = .data$e, .data$tod, .data$tsd),
        by = "timeptnames"
      ) |>
      dplyr::mutate(
        Db = (.data$tod > 0),
        Dbc = dplyr::case_when(
          .data$Db ~ 1,
          op$carryover_t1half == 0 ~ 0,
          .data$tsd <= 0 ~ 0,
          TRUE ~ (1 / 2)^(
            op$carryover_scalefactor * .data$tsd /
              op$carryover_t1half
          )
        )
      )

    datout[[g]] <- data_long
  }

  datamerged <- dplyr::bind_rows(datout)

  td_last <- trial_design_set[[n_groups]]
  e_vals <- if ("e" %in% names(td_last)) td_last$e else
    rep(0.5, nrow(td_last))
  var_in_exp <- length(unique(e_vals))

  datamerged <- datamerged |>
    dplyr::group_by(.data$ptID) |>
    dplyr::mutate(
      mean_Db = mean(.data$Db[.data$t > 0], na.rm = TRUE),
      Db_var = (.data$mean_Db != 0) & (.data$mean_Db != 1)
    ) |>
    dplyr::ungroup()

  var_in_Db <- sum(
    datamerged$Db_var[datamerged$t > 0],
    na.rm = TRUE
  ) > 0

  if (!var_in_Db && op$target_coef == "interaction") {
    ever_on_drug <- datamerged |>
      dplyr::filter(.data$t > 0, .data$mean_Db == 1) |>
      dplyr::pull(.data$ptID) |>
      unique()
    datamerged <- datamerged |>
      dplyr::filter(.data$ptID %in% ever_on_drug)
  }

  var_in_tsd <- any(datamerged$tsd != 0, na.rm = TRUE)

  # Drug indicator (Dbc) is included as a fixed effect whenever
  # it varies at all: within subject for N-of-1 / crossover, or
  # between subject for parallel RCT. The biomarker-by-treatment
  # interaction term bm:Dbc is added only when the target is the
  # interaction and Dbc varies within subject.
  any_Db_variation <- var_in_Db ||
    (length(unique(datamerged$Db)) > 1)

  model_base <- "Sx ~ bm + t"
  if (any_Db_variation) {
    model_base <- paste0(model_base, " + Dbc")
    if (op$target_coef == "interaction" && var_in_Db) {
      model_base <- paste0(model_base, " + bm:Dbc")
    }
  } else if (op$target_coef == "interaction") {
    model_base <- paste0(model_base, " + bm:t")
  }
  if ((var_in_exp > 1) && (op$useDE == TRUE)) {
    model_base <- paste0(model_base, " + De")
  }
  if (op$simplecarryover && var_in_tsd) {
    model_base <- paste0(model_base, " + tsd")
  }
  form <- stats::as.formula(model_base)

  model_vars <- c("Sx", "bm", "t", "ptID")
  if (any_Db_variation) model_vars <- c(model_vars, "Dbc")
  if (op$simplecarryover && var_in_tsd) {
    model_vars <- c(model_vars, "tsd")
  }
  if ((var_in_exp > 1) && (op$useDE == TRUE)) {
    model_vars <- c(model_vars, "De")
  }
  datamerged <- datamerged |>
    dplyr::filter(dplyr::if_all(
      dplyr::all_of(model_vars), ~ !is.na(.x)
    ))

  fit <- tryCatch(
    nlme::lme(
      form, random = ~ 1 | ptID,
      correlation = nlme::corCAR1(form = ~ t | ptID),
      data = datamerged,
      control = nlme::lmeControl(
        opt = "optim", maxIter = 200, msMaxIter = 200
      )
    ),
    error = function(e) {
      tryCatch(
        nlme::lme(
          form, random = ~ 1 | ptID,
          data = datamerged,
          control = nlme::lmeControl(opt = "optim")
        ),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(fit)) {
    out <- tibble::tibble(
      beta = NA_real_, betaSE = NA_real_,
      p = NA_real_, issingular = NA,
      warning = "model failed to converge"
    )
    if (op$full_model_out) {
      out <- list(form = form, fit = NULL,
                  datamerged = datamerged, stdout = out)
    }
    return(out)
  }

  cc <- summary(fit)$tTable
  coefnames <- rownames(cc)

  target_names <- if (op$target_coef == "main") {
    c("Dbc")
  } else if (var_in_Db) {
    c("bm:Dbc", "Dbc:bm")
  } else {
    c("bm:t", "t:bm")
  }
  target <- intersect(target_names, coefnames)
  if (length(target) == 0) {
    beta <- NA_real_; betaSE <- NA_real_; p <- NA_real_
  } else {
    target <- target[1]
    p <- cc[target, "p-value"]
    beta <- cc[target, "Value"]
    betaSE <- cc[target, "Std.Error"]
  }

  vc <- nlme::VarCorr(fit)
  re_var <- suppressWarnings(as.numeric(
    vc[rownames(vc) == "(Intercept)", "Variance"]
  ))
  issingular <- !is.na(re_var) && re_var < 1e-4

  out <- tibble::tibble(
    beta = beta, betaSE = betaSE, p = p,
    issingular = issingular, warning = NA_character_
  )

  if (op$full_model_out) {
    out <- list(form = form, fit = fit,
                datamerged = datamerged, stdout = out)
  }

  out
}
