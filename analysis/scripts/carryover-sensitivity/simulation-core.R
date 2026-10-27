## analysis/scripts/carryover-sensitivity/simulation-core.R
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
## stop_on_nonpd defaults to TRUE here and in every simulate_cell*()
## below (changed 2026-10-26; the shared build_sigma_matrix() in
## implementations/tidyverse keeps FALSE for other callers). Under
## the FALSE default, jobs 28 and 30 silently repaired non-PD
## matrices at prognostic-confound strength 0.6 and reported the
## attenuated interaction as a power loss (docs/60, 2026-10-26
## update; 30-check-prognostic-confound-pd.R). A non-PD cell now
## surfaces as a failed replicate instead. Pass stop_on_nonpd =
## FALSE explicitly to reproduce runs made before this date.

generate_data_multi_path <- function(model_param, resp_param,
                                     baseline_param, design_set,
                                     stop_on_nonpd = TRUE) {
  n_paths <- length(design_set)
  N_per_path <- rep(floor(model_param$N / n_paths), n_paths)
  ## Deterministic path allocation (RH's original approach, retained
  ## in the refined process; path allocation is not randomized).
  ## When N is not evenly divisible by n_paths, the remainder goes to
  ## the first `remainder` paths.
  remainder <- model_param$N %% n_paths
  if (remainder > 0) {
    N_per_path[seq_len(remainder)] <- N_per_path[seq_len(remainder)] + 1
  }

  map2_dfr(design_set, seq_along(design_set), function(td_path, i) {
    mp_i <- model_param
    mp_i$N <- N_per_path[i]
    if (mp_i$N == 0) return(tibble())
    dat <- generate_data(
      mp_i, resp_param, baseline_param, td_path,
      empirical = FALSE, make_positive_definite = TRUE,
      dgp_architecture = mp_i$dgp_architecture %||% 'mvn',
      stop_on_nonpd = stop_on_nonpd
    )
    dat$path <- i
    dat
  })
}

## default operator: returns y when x is NULL or NA (scalar)
default_val <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (is.na(x)) return(y)
  x
}

## `%||%` operator for NULL-only defaults (base-R convention)
`%||%` <- function(x, y) if (is.null(x)) y else x

## -----------------------------------------------------------------
## Build long-form data with Db / Dbc / L columns
## -----------------------------------------------------------------

## Exposure-weighted drug indicator: 1 on drug, the decay weight off
## drug after exposure, 0 off drug before any exposure. The design sets
## tsd = 0 before first exposure, and carryover_decay(0, h) = 1 for
## h > 0, so without the tsd > 0 guard a placebo-first CO participant
## would be coded fully exposed at every visit (docs/45, Section 8).
exposure_weight <- function(Db, tsd, halflife, form = 'exponential',
                            shape = 1) {
  w <- carryover_decay(tsd, halflife, form = form, shape = shape)
  ifelse(Db == 1, 1, ifelse(tsd > 0, w, 0))
}

prepare_long_data <- function(dat, design_set, carryover_t1half,
                              carryover_form = 'exponential',
                              weibull_shape = 1,
                              analysis_t1half = NULL,
                              analysis_form = NULL,
                              analysis_shape = NULL) {
  ## Dbc in the analysis model uses analysis_* parameters if
  ## supplied; otherwise falls back to the DGP parameters
  ## (analyst perfectly specifies truth).
  if (is.null(analysis_t1half)) analysis_t1half <- carryover_t1half
  if (is.null(analysis_form))   analysis_form   <- carryover_form
  if (is.null(analysis_shape))  analysis_shape  <- weibull_shape


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
      Dbc = exposure_weight(Db, tsd, analysis_t1half,
                            form = analysis_form,
                            shape = analysis_shape)
    ) |>
    dplyr::select(-Db_lg, -L_lg, -D_prev)

  tp_names <- unique(td_all$timepoint_name)

  dat_long <- dat |>
    dplyr::select(ptID, path, bm, all_of(tp_names)) |>
    pivot_longer(all_of(tp_names),
                 names_to = 'timepoint_name',
                 values_to = 'Sx') |>
    inner_join(
      ## De: the design-expectancy covariate (F-family models,
      ## docs/57/docs/62; De = e from build_trial_design()).
      td_all |> dplyr::select(path, timepoint_name,
                              t = t_cum, Db, Dbc, L, tsd, De = e),
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
  spec <- match.arg(spec, c('E1', 'E2', 'E3', 'E4'))
  ## E3 follows the Jones & Kenward (2014) crossover pattern: the
  ## lagged-treatment indicator L is a nuisance covariate; the
  ## biomarker interaction of interest remains bm:Db.
  ##
  ## E4 is the alternative offered by the reference implementation
  ## (implementations/original/R/lme_analysis.R, simplecarryover =
  ## TRUE): time since discontinuation enters linearly as a nuisance
  ## main effect. Like E3 it consumes no half-life, but unlike E3 it
  ## is indexed by elapsed time rather than by position in the visit
  ## sequence, so it can track a graded decline across every off-drug
  ## occasion. Its interaction target is bm:Db, the same estimand as
  ## E1 and E3.
  form <- switch(spec,
    E1 = as.formula('Sx ~ bm + t + Db  + bm:Db'),
    E2 = as.formula('Sx ~ bm + t + Dbc + bm:Dbc'),
    E3 = as.formula('Sx ~ bm + t + Db  + bm:Db + L'),
    E4 = as.formula('Sx ~ bm + t + Db  + bm:Db + tsd')
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
                  std_error = NA_real_,
                  p_value = NA_real_, converged = FALSE))
  }

  cc <- summary(fit)$tTable
  target <- switch(spec,
    E1 = intersect(c('bm:Db',  'Db:bm'),  rownames(cc)),
    E2 = intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc)),
    E3 = intersect(c('bm:Db',  'Db:bm'),  rownames(cc)),
    E4 = intersect(c('bm:Db',  'Db:bm'),  rownames(cc))
  )

  if (length(target) == 0) {
    return(tibble(spec = spec, estimate = NA_real_,
                  std_error = NA_real_,
                  p_value = NA_real_, converged = FALSE))
  }

  tibble(spec = spec,
         estimate  = cc[target[1], 'Value'],
         std_error = cc[target[1], 'Std.Error'],
         p_value   = cc[target[1], 'p-value'],
         converged = TRUE)
}

## -----------------------------------------------------------------
## Robust-inference fitters (used when robust = TRUE)
##
## Both apply a working-covariance-agnostic standard error to the A2
## (exposure-weighted Dbc) interaction. fit_cr2 keeps the conditional
## lme fit and replaces its standard error with a CR2 bias-reduced
## cluster-robust sandwich (clubSandwich). fit_gee_md fits the
## marginal GEE counterpart, with an exchangeable working
## correlation, and uses the Mancl-DeRouen bias-corrected sandwich
## (geesmv). An exchangeable working correlation is used because it
## is stable (one parameter) and handles the unbalanced clusters
## dropout produces; an AR(1) form fails under dropout and an
## unstructured form (28 parameters on eight occasions) is unstable
## at these cluster counts. See analysis/report/
## 10-interaction-test-calibration/ for the rationale.
## -----------------------------------------------------------------

fit_cr2 <- function(dat_long) {
  na_out <- tibble(spec = 'lme+CR2', estimate = NA_real_,
                   std_error = NA_real_, p_value = NA_real_,
                   converged = FALSE)
  fit <- tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dat_long,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_out)
  cc <- summary(fit)$tTable
  tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc))
  if (length(tgt) == 0) return(na_out)
  tgt <- tgt[1]
  ct <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2',
                            cluster = dat_long$ptID,
                            test = 'Satterthwaite'),
    error = function(e) NULL)
  if (is.null(ct)) return(na_out)
  cf  <- if (!is.null(ct$Coef)) as.character(ct$Coef)
         else rownames(ct)
  idx <- which(cf == tgt)
  if (length(idx) == 0) return(na_out)
  row  <- ct[idx[1], , drop = FALSE]
  grab <- function(patterns) {
    for (p in patterns) {
      hit <- grep(p, names(row), value = TRUE, ignore.case = TRUE)
      if (length(hit)) return(as.numeric(row[[hit[1]]]))
    }
    NA_real_
  }
  tibble(spec = 'lme+CR2',
         estimate  = cc[tgt, 'Value'],
         std_error = grab('^SE$'),
         p_value   = grab(c('^p_Satt', '^p_val', '^p\\.',
                            '^p$', '^p')),
         converged = TRUE)
}

fit_gee_md <- function(dat_long) {
  na_out <- tibble(spec = 'GEE+MD', estimate = NA_real_,
                   std_error = NA_real_, p_value = NA_real_,
                   converged = FALSE)
  ## geeglm/geesmv need contiguous, time-ordered clusters, a numeric
  ## cluster id, and a plain data frame.
  dl <- dat_long[order(dat_long$ptID, dat_long$t), ]
  dl$id <- as.integer(factor(dl$ptID))
  dl <- as.data.frame(dl)
  fit <- tryCatch(
    geepack::geeglm(Sx ~ bm + t + Dbc + bm:Dbc, family = gaussian,
                    data = dl, id = id, corstr = 'exchangeable'),
    error = function(e) NULL)
  if (is.null(fit)) return(na_out)
  sm  <- summary(fit)$coefficients
  tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(sm))
  if (length(tgt) == 0) return(na_out)
  tgt <- tgt[1]
  est <- sm[tgt, 'Estimate']
  md <- tryCatch({
    utils::capture.output(
      m <- geesmv::GEE.var.md(Sx ~ bm + t + Dbc + bm:Dbc,
                              id = 'id', family = gaussian,
                              data = dl, corstr = 'exchangeable'))
    m
  }, error = function(e) NULL)
  if (is.null(md)) return(na_out)
  idx <- match(tgt, names(coef(fit)))
  if (is.na(idx)) return(na_out)
  md_se <- sqrt(md$cov.beta[idx])
  n_id  <- length(unique(dl$id))
  df    <- n_id - length(coef(fit))
  tibble(spec = 'GEE+MD', estimate = est,
         std_error = md_se,
         p_value   = 2 * pt(-abs(est / md_se), df),
         converged = TRUE)
}

fit_all_specs <- function(dat_long, robust = FALSE) {
  has_off_drug <- any(dat_long$Db == 0)
  has_lagged   <- any(dat_long$L  == 1)
  has_washout  <- any(dat_long$tsd > 0)
  na_row <- function(spec, reason) {
    tibble(spec = spec, estimate = NA_real_,
           std_error = NA_real_,
           p_value = NA_real_, converged = FALSE,
           reason = reason)
  }
  out <- bind_rows(
    if (has_off_drug) fit_spec(dat_long, 'E1') |>
      mutate(reason = NA_character_)
    else na_row('E1', 'no off-drug observations'),
    fit_spec(dat_long, 'E2') |> mutate(reason = NA_character_),
    if (has_off_drug && has_lagged) fit_spec(dat_long, 'E3') |>
      mutate(reason = NA_character_)
    else na_row('E3', 'no lagged-on timepoints'),
    if (has_off_drug && has_washout) fit_spec(dat_long, 'E4') |>
      mutate(reason = NA_character_)
    else na_row('E4', 'no post-discontinuation timepoints')
  )
  if (robust) {
    out <- bind_rows(
      out,
      fit_cr2(dat_long)    |> mutate(reason = NA_character_),
      fit_gee_md(dat_long) |> mutate(reason = NA_character_)
    )
  }
  out
}

## -----------------------------------------------------------------
## Single-cell simulator
## -----------------------------------------------------------------

simulate_cell <- function(cell, n_reps, robust = FALSE,
                          resp_param = NULL,
                          baseline_param = NULL,
                          stop_on_nonpd = TRUE) {
  design_set <- build_design_set(cell$design)
  if (is.null(resp_param))     resp_param     <- default_resp_param()
  if (is.null(baseline_param)) baseline_param <- default_baseline_param()

  ## rho: common within-factor AR(1) autocorrelation. Falls back
  ## to 0.7 if not provided (Tier 1 default).
  rho <- default_val(cell[['rho']], 0.7)

  ## analyst-assumed carryover parameters for the A2 analysis
  ## predictor. Default to DGP truth (perfectly specified).
  analysis_t1half <- default_val(cell[['analysis_t1half']], cell$t1half)
  analysis_form   <- default_val(cell[['analysis_form']],   cell$carryover_form)
  analysis_shape  <- default_val(cell[['analysis_shape']],  cell$weibull_shape)

  ## optional dropout block (S3)
  dropout_rate <- default_val(cell[['dropout_rate']], 0)
  dropout_mech <- default_val(cell[['dropout_mech']], 'MCAR')

  ## For combined architecture, c_bm_a and c_bm_b are supplied
  ## separately; c_bm is unused. For mvn/mean_moderation the single
  ## c_bm value is used as before.
  is_combined <- isTRUE(cell$dgp_arch == 'combined')
  model_param <- list(
    N = cell$N,
    carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.tv = rho, c.pb = rho, c.br = rho,
    ## c.cf1t (same timepoint, lag 0) and c.cfct (other timepoints,
    ## scaled by rho^lag) together implement the refined process's
    ## single separable cross-factor correlation c_1 rho^|w_t-w_s|
    ## (docs/series/ch10-refined-process.md; paper 01 Table
    ## tab:spec, c_1 = 0.2) only when both equal c_1: at lag 0,
    ## rho^0 = 1, so c.cf1t alone must carry c_1; at lag > 0 the
    ## formula is c_1 * rho^lag, carried by c.cfct. Previously
    ## c.cfct = 0.1 here, which is not c_1 = 0.2 and so did not
    ## match the documented spec at any lag beyond zero.
    c.cf1t = 0.2, c.cfct = 0.2
  )
  if (is_combined) {
    model_param$c.bm_a <- cell$c_bm_a
    model_param$c.bm_b <- cell$c_bm_b
  } else {
    model_param$c.bm <- cell$c_bm
  }

  ## optional prognostic-baseline coupling (docs/57 Section 7.3,
  ## stress-test cells 40-41; cell field matches that convention,
  ## `cl$c_bl`). Absent or zero leaves bm and BL uncorrelated,
  ## reproducing all existing cells exactly.
  model_param$c.bm.bl <- default_val(cell[['c_bl']], 0)

  ## optional prognostic-TV coupling (prog-TV-const cell, docs/57
  ## Section 7.3; cell field `cl$c_tv`). Absent or zero leaves bm and
  ## TV uncorrelated, reproducing all existing cells exactly.
  model_param$c.bm.tv <- default_val(cell[['c_tv']], 0)

  ## optional prognostic-TV mean shift (prog-TV-grow cell, docs/57
  ## Section 7.3; cell field `cl$b_tv`). Absent or zero leaves TV
  ## unshifted, reproducing all existing cells exactly.
  model_param$b.bm.tv <- default_val(cell[['b_tv']], 0)

  ## optional prognostic-PB coupling (prog-PB cell, docs/57 Section
  ## 7.3; cell field `cl$c_pb`, matching
  ## 11-strawman-stress-test.R's convention). Not to be confused
  ## with model_param$c.pb (the AR(1) PB autocorrelation, from
  ## cell$rho above) or with the trial-design field of the same
  ## name used elsewhere. Absent or zero leaves bm and PB
  ## uncorrelated, reproducing all existing cells exactly.
  model_param$c.bm.pb <- default_val(cell[['c_pb']], 0)

  map_dfr(seq_len(n_reps), function(rep_i) {
    ## NOTE: with stop_on_nonpd = TRUE, generate_data_multi_path()
    ## raises an error on a non-PD matrix rather than repairing it,
    ## but this outer tryCatch still catches it per replicate (so a
    ## long run doesn't crash on one bad draw) and records it as a
    ## non-converged/failed replicate below, the same as any other
    ## generation failure. This differs from the package's stop()
    ## semantics (which propagate to the caller) but is the
    ## appropriate adaptation for a Monte Carlo loop: a parameter
    ## combination that is persistently non-PD will show up as an
    ## elevated non-convergence rate for its cell, visible in the
    ## usual summary, rather than as a silently-repaired matrix.
    dat <- tryCatch(
      generate_data_multi_path(model_param, resp_param,
                               baseline_param, design_set,
                               stop_on_nonpd = stop_on_nonpd),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) {
      specs <- if (robust)
        c('E1', 'E2', 'E3', 'E4', 'lme+CR2', 'GEE+MD')
      else c('E1', 'E2', 'E3', 'E4')
      return(tibble(rep = rep_i, spec = specs,
                    estimate = NA_real_,
                    std_error = NA_real_,
                    p_value = NA_real_,
                    converged = FALSE))
    }

    if (dropout_rate > 0) {
      dat <- apply_dropout(dat, design_set,
                           rate = dropout_rate,
                           mechanism = dropout_mech)
    }

    dat_long <- prepare_long_data(
      dat, design_set,
      carryover_t1half = cell$t1half,
      carryover_form = cell$carryover_form,
      weibull_shape = cell$weibull_shape,
      analysis_t1half = analysis_t1half,
      analysis_form = analysis_form,
      analysis_shape = analysis_shape
    )

    fit_all_specs(dat_long, robust = robust) |>
      mutate(rep = rep_i, .before = 1)
  })
}

## -----------------------------------------------------------------
## Dropout (S3)
##
## apply_dropout imposes missingness on timepoint outcome columns
## under two mechanisms:
##   - MCAR: uniform random probability `rate` per observation
##   - MAR:  probability scales with baseline severity (BL column),
##           with participants in the upper half of BL about twice
##           as likely to drop out as those in the lower half; the
##           marginal rate still averages to `rate`.
## Dropout is monotone: once a participant misses a timepoint,
## subsequent timepoints are also missing, matching clinical
## trial dropout conventions.
## -----------------------------------------------------------------

apply_dropout <- function(dat, design_set, rate, mechanism) {
  if (rate == 0) return(dat)
  mechanism <- match.arg(mechanism, c('MCAR', 'MAR'))

  tp_names <- design_set |>
    map(~ .x$timepoint_name) |>
    unlist() |>
    unique()

  n <- nrow(dat)
  n_tp <- length(tp_names)

  drop_prob <- if (mechanism == 'MCAR') {
    matrix(rate, nrow = n, ncol = n_tp)
  } else {
    bl_rank <- rank(dat$BL) / n
    base <- rate * 2 * bl_rank
    base <- pmin(base, 1)
    matrix(base, nrow = n, ncol = n_tp)
  }

  draws <- matrix(runif(n * n_tp), nrow = n, ncol = n_tp)
  drop_mask <- draws < drop_prob

  for (i in seq_len(n)) {
    first_drop <- which(drop_mask[i, ])[1]
    if (!is.na(first_drop)) {
      dat[i, tp_names[first_drop:n_tp]] <- NA
    }
  }

  dat
}

## -----------------------------------------------------------------
## Biomarker-interacted carryover terms and their replacements (S7)
##
## E1 is refit here (not reused from fit_all_specs) so this block is
## self-contained and does not disturb the Tier 1 pipeline. E5/E6
## (a carryover term crossed with bm) were tested in an earlier
## version of this block and found to underperform E1 at every cell
## examined (Blocks S7-S9, structurally over-parameterized relative
## to a DGP with a single true decay parameter; Section 4.5). They
## were replaced by two specifications answering different
## questions about the same open problem (the analyst's unknown
## true decay shape/half-life), rather than another interacted-term
## variant expected to fail for the same reason:
##
##   E7 = AIC-selected half-life. Refits E2's formula
##        (Sx ~ bm + t + Dbc + bm:Dbc) across a half-life grid and
##        keeps the model with the lowest AIC, mirroring
##        R/carryover_analysis.R::characterize_carryover()'s
##        candidate-grid-plus-AIC-selection logic (same
##        carryover_decay() helper used everywhere else in this
##        manuscript), reimplemented against this driver's
##        precomputed dat_long rather than calling
##        characterize_carryover()/lme_analysis() directly, whose
##        wide-data/trialdesign_set interface is incompatible with
##        this pipeline. Asks whether ESTIMATING the half-life
##        recovers what assuming it correctly would have given you
##        (Section 3.6, S2's question from the other direction).
##
##   E9 = paired-difference regression. For each patient with both
##        on- and off-drug observations, compute their own mean
##        on-drug minus mean off-drug Sx (no repeated-measures model,
##        no correlation structure at all), then regress that
##        patient-level difference on bm across patients
##        (diff ~ bm, plain OLS). A biomarker-interaction extension
##        of the paired t-test method Senn, Julious & Araujo (2014)
##        found outperforms carryover-adjusted mixed models for
##        ordinary N-of-1 treatment-effect estimation; included as a
##        genuinely different inferential paradigm (two-stage,
##        correlation-agnostic) rather than another elaboration
##        within the mixed-model family E1-E8 all belonged to.
##
## E7 targets bm:Dbc (single-coefficient test, like E2). E9 targets
## the diff~bm slope, a different estimand again (not bm:Db or
## bm:Dbc; estimate_bmDb is left NA for both). Neither adds a second
## interaction coefficient the way E5/E6 did.
## -----------------------------------------------------------------

## CR2 cluster-robust extraction shared by the E1cr2/E3cr2/E7cr2
## (G6/G7/G9) branches below. Given an already-fitted lme object and
## the target coefficient name(s) (first match wins), returns the
## CR2-corrected SE and Satterthwaite p-value, point estimate
## unchanged from the model fit (same convention as the existing
## fit_cr2() used for G8/E2's CR2 variant).
cr2_extract <- function(fit, dat_long, target_candidates) {
  cc <- summary(fit)$tTable
  tgt <- intersect(target_candidates, rownames(cc))
  if (length(tgt) == 0) return(NULL)
  tgt <- tgt[1]
  ct <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2',
                            cluster = dat_long$ptID,
                            test = 'Satterthwaite'),
    error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  cf  <- if (!is.null(ct$Coef)) as.character(ct$Coef) else rownames(ct)
  idx <- which(cf == tgt)
  if (length(idx) == 0) return(NULL)
  row  <- ct[idx[1], , drop = FALSE]
  grab <- function(patterns) {
    for (p in patterns) {
      hit <- grep(p, names(row), value = TRUE, ignore.case = TRUE)
      if (length(hit)) return(as.numeric(row[[hit[1]]]))
    }
    NA_real_
  }
  list(target = tgt, estimate = cc[tgt, 'Value'],
      std_error = grab('^SE$'),
      p_value = grab(c('^p_Satt', '^p_val', '^p\\.', '^p$', '^p')))
}

## AIC selection of the analysis half-life, shared by E7 (G4) and
## E7cr2 (G9). The candidates differ in the Dbc column, i.e. in their
## fixed-effect design, and REML log-likelihoods are not comparable
## across different fixed-effect designs, so the candidates are fitted
## by ML for the comparison and the selected model is refitted by REML
## for inference. Returns NULL if no candidate or the refit fails.
select_half_life_aic <- function(dat_long,
                                 half_lives = c(0.25, 0.5, 1.0, 2.0),
                                 random = ~1 | ptID) {
  fit_dbc <- function(dl, method) tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc, random = random,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl, method = method,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  cands <- purrr::map(half_lives, function(hl) {
    dl <- dat_long
    dl$Dbc <- exposure_weight(dl$Db, dl$tsd, hl)
    f <- fit_dbc(dl, 'ML')
    list(hl = hl, dl = dl, aic = if (is.null(f)) NA_real_ else AIC(f))
  })
  aics <- purrr::map_dbl(cands, ~ .x$aic)
  if (all(is.na(aics))) return(NULL)
  best <- cands[[which.min(aics)]]
  fit <- fit_dbc(best$dl, 'REML')
  if (is.null(fit)) return(NULL)
  list(hl = best$hl, fit = fit, dl = best$dl)
}

## G9 with a richer random-effects structure (29-run-g9-random-
## effects.R). Same AIC half-life selection and CR2 extraction as
## E7cr2; only `random` differs. Job 28 showed every specification,
## G9 included, losing power under a prog-PB or prog-TV-const
## confound with no size inflation, i.e. the confound adds patient-
## level variance correlated with bm rather than bias. The F-family's
## fixed expectancy term is identical across patients and cannot
## absorb that, so these variants model it as variance components:
## a random slope on t (patient-specific natural-history trend,
## the TV mechanism), a random coefficient on De (patient-specific
## placebo responsiveness, the PB mechanism), or both.
G9_RANDOM_VARIANTS <- list(
  G9     = ~1 | ptID,
  G9rs   = ~1 + t | ptID,
  G9re   = ~1 + De | ptID,
  G9rse  = ~1 + t + De | ptID
)

fit_g9_random <- function(dat_long, variant) {
  random <- G9_RANDOM_VARIANTS[[variant]]
  na_out <- tibble(spec = variant, estimate = NA_real_,
                   std_error = NA_real_, p_value = NA_real_,
                   best_t_half = NA_real_, converged = FALSE)
  best <- select_half_life_aic(dat_long, random = random)
  if (is.null(best)) return(na_out)
  r <- cr2_extract(best$fit, best$dl, c('bm:Dbc', 'Dbc:bm'))
  if (is.null(r)) return(na_out)
  tibble(spec = variant, estimate = r$estimate,
         std_error = r$std_error, p_value = r$p_value,
         best_t_half = best$hl, converged = TRUE)
}

fit_spec_s7 <- function(dat_long, spec) {
  spec <- match.arg(spec, c('E1', 'E2', 'E3', 'E7', 'E9',
                            'E1cr2', 'E3cr2', 'E2cr2', 'E7cr2'))

  na_out <- tibble(spec = spec, estimate = NA_real_,
                   std_error = NA_real_, p_value = NA_real_,
                   estimate_bmDb = NA_real_, p_value_bmDb = NA_real_,
                   best_t_half = NA_real_, converged = FALSE)

  ## G6/G7/G8: E1/E3/E2 refit with a CR2 cluster-robust SE instead of
  ## the model-based SE (bm:Db for E1/E3, bm:Dbc for E2).
  if (spec %in% c('E1cr2', 'E3cr2', 'E2cr2')) {
    form <- switch(spec,
      E1cr2 = as.formula('Sx ~ bm + t + Db + bm:Db'),
      E3cr2 = as.formula('Sx ~ bm + t + Db + bm:Db + L'),
      E2cr2 = as.formula('Sx ~ bm + t + Dbc + bm:Dbc'))
    fit <- tryCatch(
      nlme::lme(form, random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = dat_long,
                control = nlme::lmeControl(
                  opt = 'optim', maxIter = 200, msMaxIter = 200)),
      error = function(e) NULL)
    if (is.null(fit)) return(na_out)
    if (spec == 'E2cr2') {
      r <- cr2_extract(fit, dat_long, c('bm:Dbc', 'Dbc:bm'))
      if (is.null(r)) return(na_out)
      return(tibble(spec = spec, estimate = r$estimate,
                    std_error = r$std_error, p_value = r$p_value,
                    estimate_bmDb = NA_real_, p_value_bmDb = NA_real_,
                    best_t_half = NA_real_, converged = TRUE))
    }
    r <- cr2_extract(fit, dat_long, c('bm:Db', 'Db:bm'))
    if (is.null(r)) return(na_out)
    return(tibble(spec = spec, estimate = r$estimate,
                  std_error = r$std_error, p_value = r$p_value,
                  estimate_bmDb = r$estimate, p_value_bmDb = r$p_value,
                  best_t_half = NA_real_, converged = TRUE))
  }

  ## G9: E7 (AIC-selected half-life) refit with a CR2 cluster-robust
  ## SE on bm:Dbc at the selected half-life, instead of the
  ## model-based SE. Tests whether CR2 offsets any part of the
  ## post-selection inference risk AIC selection carries (untested
  ## elsewhere in this manuscript; CR2 was not designed for that
  ## specific problem and should not be assumed to fix it).
  if (spec == 'E7cr2') {
    best <- select_half_life_aic(dat_long)
    if (is.null(best)) return(na_out)
    r <- cr2_extract(best$fit, best$dl, c('bm:Dbc', 'Dbc:bm'))
    if (is.null(r)) return(na_out)
    return(tibble(spec = spec, estimate = r$estimate,
                  std_error = r$std_error, p_value = r$p_value,
                  estimate_bmDb = NA_real_, p_value_bmDb = NA_real_,
                  best_t_half = best$hl, converged = TRUE))
  }

  if (spec == 'E7') {
    best <- select_half_life_aic(dat_long)
    if (is.null(best)) return(na_out)
    cc <- summary(best$fit)$tTable
    tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc))
    if (length(tgt) == 0) return(na_out)
    tgt <- tgt[1]
    return(tibble(spec = spec,
                  estimate  = cc[tgt, 'Value'],
                  std_error = cc[tgt, 'Std.Error'],
                  p_value   = cc[tgt, 'p-value'],
                  estimate_bmDb = NA_real_,
                  p_value_bmDb  = NA_real_,
                  best_t_half = best$hl,
                  converged = TRUE))
  }

  ## E9: two-stage paired-difference regression. No repeated-measures
  ## model and no correlation structure at all; patients with only
  ## on-drug or only off-drug observations (e.g. CO's
  ## never-discontinuing path) cannot contribute a within-patient
  ## difference and are dropped, mirroring the "only analyze folks
  ## who are ever on drug" convention already used in
  ## R/lme_analysis.R for the analogous varInDb == FALSE case.
  if (spec == 'E9') {
    pd <- dat_long |>
      dplyr::group_by(ptID) |>
      dplyr::summarise(
        bm = dplyr::first(bm),
        has_on  = any(Db == 1),
        has_off = any(Db == 0),
        diff = mean(Sx[Db == 1], na.rm = TRUE) -
               mean(Sx[Db == 0], na.rm = TRUE),
        .groups = 'drop'
      ) |>
      dplyr::filter(has_on, has_off)
    if (nrow(pd) < 3) return(na_out)
    fit9 <- tryCatch(stats::lm(diff ~ bm, data = pd), error = function(e) NULL)
    if (is.null(fit9)) return(na_out)
    cc9 <- tryCatch(summary(fit9)$coefficients, error = function(e) NULL)
    if (is.null(cc9) || !('bm' %in% rownames(cc9))) return(na_out)
    return(tibble(spec = spec,
                  estimate  = cc9['bm', 'Estimate'],
                  std_error = cc9['bm', 'Std. Error'],
                  p_value   = cc9['bm', 'Pr(>|t|)'],
                  estimate_bmDb = NA_real_,
                  p_value_bmDb  = NA_real_,
                  best_t_half = NA_real_,
                  converged = TRUE))
  }

  form <- switch(spec,
    E1 = as.formula('Sx ~ bm + t + Db  + bm:Db'),
    E2 = as.formula('Sx ~ bm + t + Dbc + bm:Dbc'),
    E3 = as.formula('Sx ~ bm + t + Db  + bm:Db + L')
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
  if (is.null(fit)) return(na_out)

  cc <- summary(fit)$tTable

  ## E2 targets bm:Dbc, a different regressor/estimand from the
  ## other three, which all target bm:Db; estimate_bmDb is left NA
  ## for E2 rather than silently comparing across the estimand
  ## divide (Section 2.6 of the manuscript).
  if (spec == 'E2') {
    tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc))
    if (length(tgt) == 0) return(na_out)
    tgt <- tgt[1]
    return(tibble(spec = spec,
                  estimate  = cc[tgt, 'Value'],
                  std_error = cc[tgt, 'Std.Error'],
                  p_value   = cc[tgt, 'p-value'],
                  estimate_bmDb = NA_real_,
                  p_value_bmDb  = NA_real_,
                  best_t_half = NA_real_,
                  converged = TRUE))
  }

  bmdb <- intersect(c('bm:Db', 'Db:bm'), rownames(cc))
  if (length(bmdb) == 0) return(na_out)
  bmdb <- bmdb[1]

  tibble(spec = spec,
         estimate  = cc[bmdb, 'Value'],
         std_error = cc[bmdb, 'Std.Error'],
         p_value   = cc[bmdb, 'p-value'],
         estimate_bmDb = cc[bmdb, 'Value'],
         p_value_bmDb  = cc[bmdb, 'p-value'],
         best_t_half = NA_real_,
         converged = TRUE)
}

## E0: Hendrickson's actual original analysis (github.com/
## rchendrickson/pmsimstats, commit 58b32a9, R/lme_analysis.R;
## *verified* by fetching that file directly via `gh api`, 2026-10-23
## -- NOT the vendored implementations/original/, which is a later,
## already-modified copy that has since picked up nlme::lme +
## corCAR1; the actual 58b32a9 source never had that). Fits with
## lme4::lmer (via lmerTest for Satterthwaite p-values), random
## intercept only, NO residual correlation structure at all --
## corCAR1/AR(1), used by every other spec in this file (G1-G9 and
## F1-F1A0B alike), is a project addition the original paper never
## had. Binary Db (not the exposure-decayed Dbc any other spec
## uses). Model switches on whether Db varies within a patient after
## t > 0: bm:Db if it does (and patients who are never on drug are
## dropped, mirroring orig's own `iEverOnDrug` filter exactly), else
## bm:t; De is added whenever the design's expectancy varies, both
## mirroring orig's op$useDE = TRUE default. This translation uses
## this driver's long-format data (prepare_long_data()), which has
## no explicit baseline (BL, t = 0) outcome row the way orig's own
## data munging does -- a limitation shared with every other spec
## compared against it here, not unique to E0.
fit_spec_e0 <- function(dat_long) {
  na_out <- tibble(spec = 'E0', estimate = NA_real_, std_error = NA_real_,
                   p_value = NA_real_, estimate_bmDb = NA_real_,
                   p_value_bmDb = NA_real_, best_t_half = NA_real_,
                   converged = FALSE)
  dd <- dat_long
  db_tab <- dd |> dplyr::filter(t > 0) |> dplyr::group_by(ptID) |>
    dplyr::summarise(mean_db = mean(Db), .groups = 'drop')
  var_in_db <- any(db_tab$mean_db != 0 & db_tab$mean_db != 1)
  if (!var_in_db) {
    ever_on <- db_tab$ptID[db_tab$mean_db == 1]
    dd <- dplyr::filter(dd, ptID %in% ever_on)
  }
  var_in_exp <- length(unique(dd$De[dd$t > 0])) > 1
  form <- if (var_in_db) {
    if (var_in_exp) Sx ~ bm + De + Db + t + bm:Db else Sx ~ bm + Db + t + bm:Db
  } else {
    if (var_in_exp) Sx ~ bm + De + t + bm:t else Sx ~ bm + t + bm:t
  }
  fit <- tryCatch(
    lmerTest::lmer(update(form, . ~ . + (1 | ptID)), data = dd),
    error = function(e) NULL)
  if (is.null(fit)) return(na_out)
  cc <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(cc)) return(na_out)
  tgt <- if (var_in_db) intersect(c('bm:Db', 'Db:bm'), rownames(cc))
         else intersect(c('bm:t', 't:bm'), rownames(cc))
  if (length(tgt) == 0) return(na_out)
  tgt <- tgt[1]
  tibble(spec = 'E0', estimate = cc[tgt, 'Estimate'],
        std_error = cc[tgt, 'Std. Error'], p_value = cc[tgt, 'Pr(>|t|)'],
        estimate_bmDb = NA_real_, p_value_bmDb = NA_real_,
        best_t_half = NA_real_, converged = TRUE)
}

simulate_cell_s7 <- function(cell, n_reps,
                             specs = c('E1', 'E7', 'E9'),
                             stop_on_nonpd = TRUE) {
  design_set <- build_design_set(cell$design)
  resp_param     <- default_resp_param()
  baseline_param <- default_baseline_param()
  rho <- default_val(cell[['rho']], 0.7)

  ## analyst-assumed carryover parameters for E2/E2cr2's fixed Dbc
  ## predictor. Default to DGP truth (perfectly specified), same
  ## convention as simulate_cell(). E7/E7cr2 build their own Dbc via
  ## internal AIC grid search and ignore these, so a mismatch here
  ## affects only the fixed-assumption specifications, as intended.
  analysis_t1half <- default_val(cell[['analysis_t1half']], cell$t1half)
  analysis_form   <- default_val(cell[['analysis_form']],   cell$carryover_form)
  analysis_shape  <- default_val(cell[['analysis_shape']],  cell$weibull_shape)

  ## optional dropout block, same convention as simulate_cell()
  dropout_rate <- default_val(cell[['dropout_rate']], 0)
  dropout_mech <- default_val(cell[['dropout_mech']], 'MCAR')

  model_param <- list(
    N = cell$N,
    carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.bm = cell$c_bm,
    c.tv = rho, c.pb = rho, c.br = rho,
    ## See simulate_cell()'s comment above: both constants must
    ## equal c_1 = 0.2 to match the refined process's documented
    ## single separable cross-factor correlation c_1 rho^|w_t-w_s|.
    c.cf1t = 0.2, c.cfct = 0.2
  )

  ## optional prognostic-baseline coupling, same convention as
  ## simulate_cell().
  model_param$c.bm.bl <- default_val(cell[['c_bl']], 0)

  ## optional prognostic-TV coupling, same convention as
  ## simulate_cell().
  model_param$c.bm.tv <- default_val(cell[['c_tv']], 0)

  ## optional prognostic-TV mean shift, same convention as
  ## simulate_cell().
  model_param$b.bm.tv <- default_val(cell[['b_tv']], 0)

  ## optional prognostic-PB coupling, same convention as
  ## simulate_cell().
  model_param$c.bm.pb <- default_val(cell[['c_pb']], 0)

  furrr::future_map_dfr(seq_len(n_reps), function(rep_i) {
    ## See simulate_cell()'s note on stop_on_nonpd: caught per
    ## replicate here too, recorded as converged = FALSE.
    dat <- tryCatch(
      generate_data_multi_path(model_param, resp_param,
                               baseline_param, design_set,
                               stop_on_nonpd = stop_on_nonpd),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) {
      return(tibble(rep = rep_i, spec = specs,
                    estimate = NA_real_, std_error = NA_real_,
                    p_value = NA_real_, estimate_bmDb = NA_real_,
                    p_value_bmDb = NA_real_, best_t_half = NA_real_,
                    converged = FALSE))
    }

    if (dropout_rate > 0) {
      dat <- apply_dropout(dat, design_set,
                           rate = dropout_rate,
                           mechanism = dropout_mech)
    }

    dat_long <- prepare_long_data(
      dat, design_set,
      carryover_t1half = cell$t1half,
      carryover_form = cell$carryover_form,
      weibull_shape = cell$weibull_shape,
      analysis_t1half = analysis_t1half,
      analysis_form = analysis_form,
      analysis_shape = analysis_shape
    )

    map_dfr(specs, ~ fit_spec_s7(dat_long, .x)) |>
      mutate(rep = rep_i, .before = 1)
  }, .options = furrr::furrr_options(seed = TRUE))
}

## -----------------------------------------------------------------
## F-family analysis models (docs/57, docs/62): the richer model
## Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De. bm is already centered
## at 0 by the DGP's baseline_param (bm_m = 0), so bm serves as bmc
## directly. De requires prepare_long_data() above to carry the
## expectancy column, which it now does. F1 fixes t1/2^a at the
## project's policy value (0.2 week, CLAUDE.md); F1A selects it by
## AIC with the interaction already in the model (the liberal
## pattern); F1A0 selects it by AIC on the null (no-interaction)
## model and refits with the interaction at the selected value;
## F1A0B calibrates F1A0's test via a warp-speed parametric
## bootstrap, one resample per Monte Carlo replicate, pooled across
## the cell (docs/62, "Bootstrap calibration").
##
## This driver's 'CO' design preset (expectancies = rep(0.5, 8)) has
## no variation in De, so De and bm:De are exactly collinear with the
## intercept and bm under CO: FM1_F/FM0_F are not identifiable there
## and every fit returns NA (verified). Use 'OLBDC' instead for the
## design-robustness arm; its expectancy varies (1.0 open-label vs.
## 0.5 blinded discontinuation).
## -----------------------------------------------------------------

F_GRID <- c(0, 0.2, 0.5, 1, 2)
FM1_F <- Sx ~ bm + t + Dbc + bm:Dbc + De + bm:De
FM0_F <- Sx ~ bm + t + Dbc + De + bm:De

fit_lme_f <- function(form, dat_long, method = 'REML') {
  tryCatch(
    nlme::lme(form, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dat_long, method = method,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
}

ll_or_na_f <- function(f) if (is.null(f)) NA_real_ else as.numeric(logLik(f))

with_dbc_f <- function(dat_long, h, analysis_form = 'exponential',
                       analysis_shape = 1) {
  dat_long$Dbc <- exposure_weight(dat_long$Db, dat_long$tsd, h,
                                  form = analysis_form, shape = analysis_shape)
  dat_long
}

## ML log-likelihoods with (ll1) and without (ll0) the interaction at
## every grid half-life, plus the REML + CR2 test of bm:Dbc at FM1_F.
## Refits FM1_F by REML at every candidate (not just the one
## eventually selected) so F1/F1A/F1A0 share one pass over the grid;
## costs 2 ML fits + 1 REML fit per grid point, matching docs/62's
## reported 5-6x-a-single-fit cost for the data-driven methods.
grid_fits_f <- function(dat_long, grid = F_GRID,
                        analysis_form = 'exponential', analysis_shape = 1) {
  map_dfr(grid, function(h) {
    dd <- with_dbc_f(dat_long, h, analysis_form, analysis_shape)
    fit1_reml <- fit_lme_f(FM1_F, dd, 'REML')
    r <- if (is.null(fit1_reml)) NULL else cr2_extract(fit1_reml, dd, c('bm:Dbc', 'Dbc:bm'))
    tibble(h = h,
          ll1 = ll_or_na_f(fit_lme_f(FM1_F, dd, 'ML')),
          ll0 = ll_or_na_f(fit_lme_f(FM0_F, dd, 'ML')),
          estimate = if (is.null(r)) NA_real_ else r$estimate,
          std_error = if (is.null(r)) NA_real_ else r$std_error,
          p_value = if (is.null(r)) NA_real_ else r$p_value)
  })
}

pick_f <- function(g, k, name) {
  if (length(k) == 0 || is.na(k)) {
    return(tibble(method = name, h = NA_real_, estimate = NA_real_,
                  std_error = NA_real_, p_value = NA_real_))
  }
  tibble(method = name, h = g$h[k], estimate = g$estimate[k],
        std_error = g$std_error[k], p_value = g$p_value[k])
}
wmax_f <- function(x) if (all(is.na(x))) NA_integer_ else which.max(x)

grid_methods_f <- function(g, h_fixed = 0.2) {
  bind_rows(
    pick_f(g, which(g$h == h_fixed), 'F1'),
    pick_f(g, wmax_f(g$ll1), 'F1A'),
    pick_f(g, wmax_f(g$ll0), 'F1A0')
  )
}

## Null-model parametric resample at h0 (F1A0's own ML choice): FM0_F
## refit by REML, simulated with the fitted random-intercept +
## corCAR1 covariance structure, on the same participants, visits,
## biomarker values and paths. Mirrors
## analysis/scripts/data-chosen-half-life/57-common.R's sim_null().
sim_null_f <- function(fit, dat_long, h0, analysis_form = 'exponential',
                       analysis_shape = 1) {
  dd <- with_dbc_f(dat_long, h0, analysis_form, analysis_shape)
  mu <- drop(model.matrix(FM0_F[-2], dd) %*% nlme::fixef(fit))
  tau <- as.numeric(nlme::VarCorr(fit)['(Intercept)', 'StdDev'])
  phi <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  sig <- fit$sigma
  y <- mu
  for (idx in split(seq_len(nrow(dd)), dd$ptID)) {
    tt <- dd$t[idx]
    S <- tau^2 + sig^2 * phi^abs(outer(tt, tt, '-'))
    y[idx] <- mu[idx] + drop(crossprod(chol(S), rnorm(length(idx))))
  }
  ds <- dd
  ds$Sx <- y
  ds
}

## F1A0's |t| statistic on one resample, rerunning the entire
## selection (not just the test).
boot_stat_f1a0 <- function(ds, grid = F_GRID, analysis_form = 'exponential',
                           analysis_shape = 1) {
  g <- grid_fits_f(ds, grid, analysis_form, analysis_shape)
  m <- grid_methods_f(g)
  row <- dplyr::filter(m, method == 'F1A0')
  if (nrow(row) == 0 || is.na(row$estimate[1]) || is.na(row$std_error[1])) {
    return(NA_real_)
  }
  abs(row$estimate[1] / row$std_error[1])
}

## One trial, all four F-family methods. Returns a 3-row table (F1,
## F1A, F1A0, each with t = estimate/std_error) plus one warp-speed
## bootstrap resample's F1A0 |t| statistic (boot_t_f1a0, identical
## across the three rows) for later per-cell pooling into F1A0B's
## calibrated size/power (docs/62; ws_size-style comparison of the
## observed |t| against the pooled resampled |t| of the same cell).
## `grid` overrides F_GRID (docs/62 criterion 3's "what if the grid
## is denser" stress test, 2026-10-23: a denser grid gives AIC more
## candidate half-lives to overfit to under the null, which is the
## scenario F1A0B's bootstrap correction is meant to guard against
## but that the 5-point project grid never actually exercises).
fit_f_family <- function(dat_long, n_boot = 1L, analysis_form = 'exponential',
                         analysis_shape = 1, grid = F_GRID) {
  g <- grid_fits_f(dat_long, grid = grid, analysis_form = analysis_form,
                   analysis_shape = analysis_shape)
  m <- grid_methods_f(g)
  m <- mutate(m, converged = !is.na(estimate) & !is.na(std_error),
             t = estimate / std_error)
  h0_rows <- dplyr::filter(m, method == 'F1A0')
  h0 <- if (nrow(h0_rows) == 1) h0_rows$h[1] else NA_real_
  boot_t <- NA_real_
  if (n_boot > 0 && length(h0) == 1 && is.finite(h0)) {
    nf <- fit_lme_f(FM0_F, with_dbc_f(dat_long, h0, analysis_form, analysis_shape), 'REML')
    if (!is.null(nf)) {
      boot_vals <- vapply(seq_len(n_boot), function(b)
        boot_stat_f1a0(sim_null_f(nf, dat_long, h0, analysis_form, analysis_shape),
                       grid = grid, analysis_form = analysis_form,
                       analysis_shape = analysis_shape),
        numeric(1))
      boot_t <- mean(boot_vals, na.rm = TRUE)
    }
  }
  mutate(m, boot_t_f1a0 = boot_t)
}

## Warp-speed calibrated rejection rate: fraction of observed |t|
## values (a) exceeding the quantile of the pooled resampled |t|
## values (s) implied by alpha, p* = (1 + #{s >= a}) / (n + 1) < alpha.
ws_size <- function(a, s, alpha = 0.05) {
  s <- sort(s[is.finite(s)]); a <- a[is.finite(a)]; n <- length(s)
  if (n == 0 || length(a) == 0) return(NA_real_)
  mean((1 + n - findInterval(a, s, left.open = TRUE)) / (n + 1) < alpha)
}

## F-family counterpart to simulate_cell(): same cell convention
## (design, N, t1half, carryover_form, weibull_shape, dgp_arch, c_bm,
## plus the optional prognostic fields c_bl/c_tv/b_tv), fitting the
## four F-family methods instead of E1-E4. n_boot is the number of
## warp-speed resamples per replicate (1 is the standard warp-speed
## shortcut; more reduces resample noise at proportional extra cost).
simulate_cell_f <- function(cell, n_reps, n_boot = 1L,
                            resp_param = NULL, baseline_param = NULL,
                            stop_on_nonpd = TRUE) {
  design_set <- build_design_set(cell$design)
  if (is.null(resp_param)) resp_param <- default_resp_param()
  if (is.null(baseline_param)) baseline_param <- default_baseline_param()

  rho <- default_val(cell[['rho']], 0.7)
  analysis_t1half <- default_val(cell[['analysis_t1half']], cell$t1half)
  analysis_form   <- default_val(cell[['analysis_form']], cell$carryover_form)
  analysis_shape  <- default_val(cell[['analysis_shape']], cell$weibull_shape)
  f_grid          <- if (is.null(cell[['f_grid']])) F_GRID else cell[['f_grid']]
  is_combined <- isTRUE(cell$dgp_arch == 'combined')

  model_param <- list(
    N = cell$N, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.tv = rho, c.pb = rho, c.br = rho,
    c.cf1t = 0.2, c.cfct = 0.2
  )
  if (is_combined) {
    model_param$c.bm_a <- cell$c_bm_a
    model_param$c.bm_b <- cell$c_bm_b
  } else {
    model_param$c.bm <- cell$c_bm
  }
  model_param$c.bm.bl <- default_val(cell[['c_bl']], 0)
  model_param$c.bm.tv <- default_val(cell[['c_tv']], 0)
  model_param$b.bm.tv <- default_val(cell[['b_tv']], 0)
  model_param$c.bm.pb <- default_val(cell[['c_pb']], 0)

  map_dfr(seq_len(n_reps), function(rep_i) {
    dat <- tryCatch(
      generate_data_multi_path(model_param, resp_param, baseline_param,
                               design_set, stop_on_nonpd = stop_on_nonpd),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) {
      return(tibble(rep = rep_i,
                    method = c('F1', 'F1A', 'F1A0'),
                    h = NA_real_, estimate = NA_real_, std_error = NA_real_,
                    p_value = NA_real_, converged = FALSE, t = NA_real_,
                    boot_t_f1a0 = NA_real_))
    }
    dat_long <- prepare_long_data(
      dat, design_set,
      carryover_t1half = cell$t1half,
      carryover_form = cell$carryover_form,
      weibull_shape = cell$weibull_shape,
      analysis_t1half = analysis_t1half,
      analysis_form = analysis_form,
      analysis_shape = analysis_shape
    )
    fit_f_family(dat_long, n_boot = n_boot,
                analysis_form = analysis_form,
                analysis_shape = analysis_shape,
                grid = f_grid) |>
      as_tibble() |>
      mutate(rep = rep_i, .before = 1)
  })
}

## -----------------------------------------------------------------
## G-family (G1-G9, via fit_spec_s7's E1-E9/+cr2 specs) vs F-family
## head-to-head: one DGP draw per replicate, both families fit to
## the *same* simulated trial (not independently resampled), so the
## comparison is not diluted by extra Monte Carlo noise between
## families. docs/62 Part 3 names this comparison as missing
## entirely ("no single document currently compares the two
## families directly") and names criterion 7 (prognostic biomarker)
## as untested for every G-family specification. Added 2026-10-23.
## Also fits E0, Hendrickson's actual original analysis (see
## fit_spec_e0() above) alongside the nine G-family specs, labeled
## family = 'G' for convenience though it predates that catalog and
## is not really part of it.
## Same cell convention as simulate_cell_f()/simulate_cell_s7();
## does not support the dual-channel ('combined') architecture or
## dropout, neither used by any cell this function has been run on.
simulate_cell_combined <- function(cell, n_reps, n_boot = 1L,
                                   g_specs = c('E1', 'E2', 'E3', 'E7', 'E9',
                                              'E1cr2', 'E3cr2', 'E2cr2', 'E7cr2'),
                                   stop_on_nonpd = TRUE) {
  design_set <- build_design_set(cell$design)
  resp_param     <- default_resp_param()
  baseline_param <- default_baseline_param()
  rho <- default_val(cell[['rho']], 0.7)
  analysis_t1half <- default_val(cell[['analysis_t1half']], cell$t1half)
  analysis_form   <- default_val(cell[['analysis_form']],   cell$carryover_form)
  analysis_shape  <- default_val(cell[['analysis_shape']],  cell$weibull_shape)
  f_grid <- if (is.null(cell[['f_grid']])) F_GRID else cell[['f_grid']]

  model_param <- list(
    N = cell$N, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.bm = cell$c_bm,
    c.tv = rho, c.pb = rho, c.br = rho,
    c.cf1t = 0.2, c.cfct = 0.2
  )
  model_param$c.bm.bl <- default_val(cell[['c_bl']], 0)
  model_param$c.bm.tv <- default_val(cell[['c_tv']], 0)
  model_param$b.bm.tv <- default_val(cell[['b_tv']], 0)
  model_param$c.bm.pb <- default_val(cell[['c_pb']], 0)
  ## optional dropout block, same convention as simulate_cell()
  ## (added 2026-10-26 for 32-run-g9-gap-fill.R)
  dropout_rate <- default_val(cell[['dropout_rate']], 0)
  dropout_mech <- default_val(cell[['dropout_mech']], 'MCAR')

  map_dfr(seq_len(n_reps), function(rep_i) {
    dat <- tryCatch(
      generate_data_multi_path(model_param, resp_param, baseline_param,
                               design_set, stop_on_nonpd = stop_on_nonpd),
      error = function(e) NULL
    )
    ## A failed generation (a non-PD cell under stop_on_nonpd = TRUE)
    ## is recorded as a non-converged row per method so the cell's
    ## nonconv rate shows it, rather than dropping the replicate.
    if (is.null(dat) || nrow(dat) == 0) {
      return(tibble(rep = rep_i,
                    family = rep(c('G', 'F'), c(length(g_specs) + 1L, 3L)),
                    method = c(g_specs, 'E0', 'F1', 'F1A', 'F1A0'),
                    h = NA_real_, estimate = NA_real_, std_error = NA_real_,
                    p_value = NA_real_, converged = FALSE, t = NA_real_,
                    boot_t_f1a0 = NA_real_))
    }
    if (dropout_rate > 0) {
      dat <- apply_dropout(dat, design_set,
                           rate = dropout_rate,
                           mechanism = dropout_mech)
    }

    dat_long <- prepare_long_data(
      dat, design_set,
      carryover_t1half = cell$t1half,
      carryover_form = cell$carryover_form,
      weibull_shape = cell$weibull_shape,
      analysis_t1half = analysis_t1half,
      analysis_form = analysis_form,
      analysis_shape = analysis_shape
    )

    g_res <- bind_rows(map_dfr(g_specs, ~ fit_spec_s7(dat_long, .x)),
                      fit_spec_e0(dat_long)) |>
      transmute(family = 'G', method = spec, h = best_t_half,
               estimate, std_error, p_value, converged,
               t = estimate / std_error, boot_t_f1a0 = NA_real_)

    f_res <- fit_f_family(dat_long, n_boot = n_boot,
                          analysis_form = analysis_form,
                          analysis_shape = analysis_shape,
                          grid = f_grid) |>
      as_tibble() |>
      transmute(family = 'F', method, h, estimate, std_error,
               p_value, converged, t, boot_t_f1a0)

    bind_rows(g_res, f_res) |> mutate(rep = rep_i, .before = 1)
  })
}

## G9 and its random-effects variants (G9_RANDOM_VARIANTS) on the
## same simulated trial. Data generation is the same sequence of
## calls as simulate_cell_combined(), so under the same seed and
## platform the trial is identical to job 28's for a matching cell.
simulate_cell_g9_random <- function(cell, n_reps,
                                    variants = names(G9_RANDOM_VARIANTS),
                                    stop_on_nonpd = TRUE) {
  design_set <- build_design_set(cell$design)
  resp_param     <- default_resp_param()
  baseline_param <- default_baseline_param()
  rho <- default_val(cell[['rho']], 0.7)

  model_param <- list(
    N = cell$N, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.bm = cell$c_bm,
    c.tv = rho, c.pb = rho, c.br = rho,
    c.cf1t = 0.2, c.cfct = 0.2
  )
  model_param$c.bm.bl <- default_val(cell[['c_bl']], 0)
  model_param$c.bm.tv <- default_val(cell[['c_tv']], 0)
  model_param$b.bm.tv <- default_val(cell[['b_tv']], 0)
  model_param$c.bm.pb <- default_val(cell[['c_pb']], 0)

  map_dfr(seq_len(n_reps), function(rep_i) {
    dat <- tryCatch(
      generate_data_multi_path(model_param, resp_param, baseline_param,
                               design_set, stop_on_nonpd = stop_on_nonpd),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) {
      return(tibble(rep = rep_i, method = variants, h = NA_real_,
                    estimate = NA_real_, std_error = NA_real_,
                    p_value = NA_real_, converged = FALSE, t = NA_real_))
    }

    dat_long <- prepare_long_data(
      dat, design_set,
      carryover_t1half = cell$t1half,
      carryover_form = cell$carryover_form,
      weibull_shape = cell$weibull_shape
    )

    map_dfr(variants, ~ fit_g9_random(dat_long, .x)) |>
      transmute(method = spec, h = best_t_half, estimate, std_error,
                p_value, converged, t = estimate / std_error) |>
      mutate(rep = rep_i, .before = 1)
  })
}
