#' 18-null-shift-check.R
#'
#' Check of the mean of the main-effect estimate under the null.
#'
#' The 2026-10-29 rerun (zeus job 51) showed that the alternating
#' design and the two-period crossover have a mean estimate below zero
#' under the null (S1, S12). The analysis model of a blinded
#' within-subject design is `Sx ~ bm + t + Dbc`. The baseline row has
#' expectancy 0 and every follow-up visit 0.5, but no expectancy term
#' enters, because the engine counts expectancy values over the
#' follow-up visits only. This script tests whether the shift comes
#' from that mean model. One large null dataset per design is fitted
#' four ways: the engine formula; the engine formula plus a follow-up
#' indicator (1 at t > 0); the engine formula without the baseline
#' row; and visit as a factor in place of linear time.
#'
#' One dataset per design, 3,000 participants per path (2,000 for the
#' Hybrid), baseline configuration, analysis half-life 0.2 week. The
#' coefficient of one large fit approximates the large-sample value.
#' Takes about two minutes on one core. No simulation sweep is run.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 18L
SWEEP_NAME <- 'null-shift-check'
set.seed(42 + SWEEP_ID)

ensure_tp <- get('ensure_timepoint_name', envir = asNamespace('nof1power'))

fit_ways <- function(td, label, n_per_path) {
  paths <- td$trialpaths
  rp <- baseline_resp_param(br_max = 0)
  bp <- baseline_bl_param()
  mp <- baseline_model_param(N = n_per_path)
  dat <- dplyr::bind_rows(lapply(seq_along(paths), function(i) {
    d <- generate_data(mp, rp, bp, ensure_tp(paths[[i]]),
                       empirical = FALSE, make_positive_definite = TRUE)
    d$path <- i
    d
  }))
  ap <- as.list(default_analysisparams())
  ap$full_model_out <- TRUE
  out <- lme_analysis(paths, dat, ap)
  dm <- out$datamerged
  dm$fu <- as.numeric(dm$t > 0)
  ctl <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)
  base_rhs <- paste(attr(stats::terms(out$form), 'term.labels'),
                    collapse = ' + ')
  fit <- function(rhs, data) {
    tryCatch(
      nlme::lme(stats::as.formula(paste('Sx ~', rhs)),
                random = ~ 1 | ptID,
                correlation = nlme::corCAR1(form = ~ t | ptID),
                data = data, control = ctl),
      error = function(e) NULL)
  }
  # Expectancy is a function of visit, so De is dropped when visit
  # enters as a factor.
  no_t <- gsub('(^| )t( \\+|$)', '\\1factor(t)\\2', base_rhs)
  no_t <- gsub(' \\+ De', '', no_t)
  fits <- list(
    engine = NULL,
    follow_up = fit(paste(base_rhs, '+ fu'), dm),
    no_baseline_row = fit(base_rhs, dm[dm$t > 0, ]),
    visit_factor = fit(no_t, dm)
  )
  get_d <- function(f) {
    if (is.null(f)) return(c(NA_real_, NA_real_))
    summary(f)$tTable['Dbc', c('Value', 'Std.Error')]
  }
  rows <- lapply(names(fits), function(nm) {
    v <- if (nm == 'engine') {
      c(out$stdout$beta, out$stdout$betaSE)
    } else {
      get_d(fits[[nm]])
    }
    tibble::tibble(design_name = label, fit = nm, beta = v[[1]],
                   se = v[[2]], n_participants = n_per_path * length(paths),
                   engine_formula = deparse(out$form))
  })
  dplyr::bind_rows(rows)
}

res <- list(results = dplyr::bind_rows(
  fit_ways(design_crossover(), 'XO', 3000),
  fit_ways(design_nof1_alternating(), 'Alt', 3000),
  fit_ways(design_hybrid(), 'Hybrid', 2000)
))
print(res$results)
save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = 1L,
           notes = paste('One large null dataset per design; Dbc',
                         'coefficient under four mean models.'))
