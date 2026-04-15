## analysis/02-carryover-sensitivity/simulation-core.R
##
## Shared simulation helpers for the carryover-sensitivity factorial
## study. Sourced by 01-run-factorial.R.
##
## Builds on implementations/tidyverse/R/functions.R:
##   - generate_data(): single-path MVN data generation with
##     carryover_form / weibull_shape plumbed through
##   - build_trial_design(): single-path trial design object
##
## Design presets follow the canonical Hendrickson (2020) forms
## reproduced in analysis/figure5/01-generate-parameter-sensitivity.R.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(nlme)
})

## -----------------------------------------------------------------
## Trial-design presets
## -----------------------------------------------------------------

design_preset <- function(name) {
  name <- match.arg(name, c('OL', 'CO', 'Hybrid', 'OLBDC'))
  switch(name,
    OL = list(
      timepoints    = cumsum(rep(2.5, 8)),
      timeptnames   = paste0('OL', 1:8),
      expectancies  = rep(1, 8),
      ondrug        = list(pathA = rep(1, 8))
    ),
    CO = list(
      timepoints    = cumsum(rep(2.5, 8)),
      timeptnames   = c(paste0('COa', 1:4), paste0('COb', 1:4)),
      expectancies  = rep(0.5, 8),
      ondrug        = list(
        pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
        pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
      )
    ),
    OLBDC = list(
      timepoints    = c(4, 8, 12, 16, 17, 18, 19, 20),
      timeptnames   = c('OL1', 'OL2', 'OL3', 'OL4',
                        'BD1', 'BD2', 'BD3', 'BD4'),
      expectancies  = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
      ondrug        = list(
        pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
        pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
      )
    ),
    Hybrid = list(
      timepoints    = c(4, 8, 9, 10, 11, 12, 16, 20),
      timeptnames   = c('OL1', 'OL2', 'BD1', 'BD2',
                        'BD3', 'BD4', 'COd', 'COp'),
      expectancies  = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
      ondrug        = list(
        pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
        pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
        pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
        pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
      )
    )
  )
}

build_design_set <- function(name) {
  spec <- design_preset(name)
  map(spec$ondrug, function(od) {
    build_trial_design(
      timepoints   = spec$timepoints,
      timeptnames  = spec$timeptnames,
      expectancies = spec$expectancies,
      ondrug       = list(od)
    )$trialpaths[[1]] |>
      dplyr::rename(any_of(c(timepoint_name = 'timeptnames')))
  })
}

## -----------------------------------------------------------------
## Default response / baseline parameters
## -----------------------------------------------------------------

default_resp_param <- function() {
  ## br sd = 8 matches the docs/19 calibration target
  ## (c_bm * sigma_BR = 0.45 * 8 = 3.6 BR units per SD of biomarker).
  ## tv / pb values follow the Hendrickson (2020) defaults; tv sd is
  ## bumped to 6 so the TV trajectory is not dominated by BR noise.
  tibble::tibble(
    cat  = c('tv', 'pb', 'br'),
    max  = c(10.98604, 6.50647, 10.98604),
    disp = c(5, 5, 5),
    rate = c(0.42, 0.35, 0.42),
    sd   = c(6, 2, 8)
  )
}

default_baseline_param <- function() {
  tibble::tibble(
    cat = c('bm', 'BL'),
    m   = c(0, 70),
    sd  = c(1, 10)
  )
}

## -----------------------------------------------------------------
## Generate data across paths for a multi-path design
## -----------------------------------------------------------------

generate_data_multi_path <- function(model_param, resp_param,
                                     baseline_param, design_set) {
  n_paths <- length(design_set)
  N_per_path <- rep(floor(model_param$N / n_paths), n_paths)
  N_per_path[seq_len(model_param$N %% n_paths)] <-
    N_per_path[seq_len(model_param$N %% n_paths)] + 1

  map2_dfr(design_set, seq_along(design_set), function(td_path, i) {
    mp_i <- model_param
    mp_i$N <- N_per_path[i]
    if (mp_i$N == 0) return(tibble())
    dat <- generate_data(
      mp_i, resp_param, baseline_param, td_path,
      empirical = FALSE, make_positive_definite = TRUE,
      dgp_architecture = mp_i$dgp_architecture %||% 'mvn'
    )
    dat$path <- i
    dat
  })
}

## `%||%` operator for default values
`%||%` <- function(x, y) if (is.null(x)) y else x

## -----------------------------------------------------------------
## Build long-form data with Db / Dbc / L columns
## -----------------------------------------------------------------

prepare_long_data <- function(dat, design_set, carryover_t1half,
                              carryover_form = 'exponential',
                              weibull_shape = 1) {

  design_set <- map(seq_along(design_set), function(i) {
    td <- design_set[[i]]
    td$path <- i
    td$t_cum <- cumsum(td$t_wk)
    td
  })

  td_all <- bind_rows(design_set) |>
    group_by(path) |>
    mutate(
      D_prev = lag(tod > 0, default = FALSE),
      Db_lg  = tod > 0,
      L_lg   = D_prev & !Db_lg
    ) |>
    ungroup() |>
    mutate(
      Db  = as.numeric(Db_lg),
      L   = as.numeric(L_lg),
      Dbc = dplyr::case_when(
        Db_lg ~ 1,
        TRUE  ~ carryover_decay(tsd, carryover_t1half,
                                form = carryover_form,
                                shape = weibull_shape)
      )
    ) |>
    dplyr::select(-Db_lg, -L_lg, -D_prev)

  tp_names <- unique(td_all$timepoint_name)

  dat_long <- dat |>
    dplyr::select(ptID, path, bm, all_of(tp_names)) |>
    pivot_longer(all_of(tp_names),
                 names_to = 'timepoint_name',
                 values_to = 'Sx') |>
    inner_join(
      td_all |> dplyr::select(path, timepoint_name,
                              t = t_cum, Db, Dbc, L, tsd),
      by = c('path', 'timepoint_name')
    ) |>
    mutate(ptID = paste(path, ptID, sep = '_')) |>
    filter(!is.na(Sx))

  dat_long
}

## -----------------------------------------------------------------
## Fit three analysis-model carryover specifications
## -----------------------------------------------------------------

fit_spec <- function(dat_long, spec) {
  spec <- match.arg(spec, c('A1', 'A2', 'A3'))
  ## A3 follows the Jones & Kenward (2014) crossover pattern: the
  ## lagged-treatment indicator L is a nuisance covariate; the
  ## biomarker interaction of interest remains bm:Db.
  form <- switch(spec,
    A1 = as.formula('Sx ~ bm + t + Db  + bm:Db'),
    A2 = as.formula('Sx ~ bm + t + Dbc + bm:Dbc'),
    A3 = as.formula('Sx ~ bm + t + Db  + bm:Db + L')
  )

  fit <- tryCatch(
    nlme::lme(form, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dat_long,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) {
      tryCatch(nlme::lme(form, random = ~1 | ptID, data = dat_long,
                         control = nlme::lmeControl(opt = 'optim')),
               error = function(e2) NULL)
    }
  )

  if (is.null(fit)) {
    return(tibble(spec = spec, estimate = NA_real_,
                  p_value = NA_real_, converged = FALSE))
  }

  cc <- summary(fit)$tTable
  target <- switch(spec,
    A1 = intersect(c('bm:Db',  'Db:bm'),  rownames(cc)),
    A2 = intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc)),
    A3 = intersect(c('bm:Db',  'Db:bm'),  rownames(cc))
  )

  if (length(target) == 0) {
    return(tibble(spec = spec, estimate = NA_real_,
                  p_value = NA_real_, converged = FALSE))
  }

  tibble(spec = spec,
         estimate = cc[target[1], 'Value'],
         p_value  = cc[target[1], 'p-value'],
         converged = TRUE)
}

fit_three_specs <- function(dat_long) {
  has_off_drug <- any(dat_long$Db == 0)
  has_lagged   <- any(dat_long$L  == 1)
  na_row <- function(spec, reason) {
    tibble(spec = spec, estimate = NA_real_,
           p_value = NA_real_, converged = FALSE,
           reason = reason)
  }
  bind_rows(
    if (has_off_drug) fit_spec(dat_long, 'A1') |>
      mutate(reason = NA_character_)
    else na_row('A1', 'no off-drug observations'),
    fit_spec(dat_long, 'A2') |> mutate(reason = NA_character_),
    if (has_off_drug && has_lagged) fit_spec(dat_long, 'A3') |>
      mutate(reason = NA_character_)
    else na_row('A3', 'no lagged-on timepoints')
  )
}

## -----------------------------------------------------------------
## Single-cell simulator
## -----------------------------------------------------------------

simulate_cell <- function(cell, n_reps) {
  design_set <- build_design_set(cell$design)
  resp_param <- default_resp_param()
  baseline_param <- default_baseline_param()

  model_param <- list(
    N = cell$N,
    c.bm = cell$c_bm,
    carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.2, c.cfct = 0.1
  )

  map_dfr(seq_len(n_reps), function(rep_i) {
    dat <- tryCatch(
      generate_data_multi_path(model_param, resp_param,
                               baseline_param, design_set),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) {
      return(tibble(rep = rep_i, spec = c('A1', 'A2', 'A3'),
                    estimate = NA_real_, p_value = NA_real_,
                    converged = FALSE))
    }

    dat_long <- prepare_long_data(
      dat, design_set,
      carryover_t1half = cell$t1half,
      carryover_form = cell$carryover_form,
      weibull_shape = cell$weibull_shape
    )

    fit_three_specs(dat_long) |>
      mutate(rep = rep_i, .before = 1)
  })
}
