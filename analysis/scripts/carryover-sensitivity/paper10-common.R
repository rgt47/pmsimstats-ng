## analysis/scripts/carryover-sensitivity/paper10-common.R
##
## Shared set-up for the paper 10 diagnostic drivers
## (diagnose-null-type1*.R, diagnose-alt-cell-cr2.R; manuscript
## analysis/report/10-interaction-test-calibration/). Sourced after
## implementations/tidyverse/R/functions.R and simulation-core.R.
##
## DGP (corrected 2026-10-29, referee report 2026-10-28, Major 2):
##   - c.cfct = 0.2 = c.cf1t, the single cross-component constant c_1
##     of the refined process. The drivers ran with c.cfct = 0.1 until
##     2026-10-29, the mismatch documented in docs/61.
##   - carryover on the drug-response mean decays from the last on-drug
##     mean (the tidyverse apply_carryover_to_component() since
##     2026-10-08); the 2026-05-20 runs used the original compounding
##     rule.
##   - stop_on_nonpd = TRUE: a non-positive-definite covariance stops
##     generation and the replicate is recorded as failed.
## Outputs go to output/paper10-refined/. The pre-correction outputs
## in output/ (diag-*.rds, 2026-05-20, relabeled 2026-07-31) are kept
## unchanged for comparison.

P10_CFCT <- 0.2

p10_out_dir <- function(repo_root) {
  d <- file.path(repo_root,
    'analysis/scripts/carryover-sensitivity/output/paper10-refined')
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

p10_model_param <- function(N, c_bm, t1half = 1.0, rho = 0.7,
                            carryover_form = 'exponential',
                            weibull_shape = 1, dgp_arch = 'mvn') {
  list(
    N = N, c.bm = c_bm,
    carryover_t1half = t1half,
    carryover_form   = carryover_form,
    weibull_shape    = weibull_shape,
    dgp_architecture = dgp_arch,
    c.tv = rho, c.pb = rho, c.br = rho,
    c.cf1t = 0.2, c.cfct = P10_CFCT)
}

p10_dgp_meta <- function(repo_root) {
  commit <- tryCatch(
    system2('git', c('-C', repo_root, 'rev-parse', '--short', 'HEAD'),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_)
  dirty <- tryCatch(
    length(system2('git', c('-C', repo_root, 'status', '--porcelain',
                            'implementations/tidyverse/R/functions.R',
                            'analysis/scripts/carryover-sensitivity'),
                   stdout = TRUE, stderr = FALSE)) > 0,
    error = function(e) NA)
  list(c.cf1t = 0.2, c.cfct = P10_CFCT,
       carryover_rule = 'decays from the last on-drug mean',
       stop_on_nonpd = TRUE, commit = commit,
       working_tree_dirty = dirty, run_date = as.character(Sys.time()),
       r_version = R.version.string)
}

## Generate one dataset in long form, or NULL when generation fails
## (a non-PD covariance under stop_on_nonpd = TRUE).
p10_gen_long <- function(model_param, design_set, t1half = 1.0,
                         carryover_form = 'exponential',
                         weibull_shape = 1) {
  dat <- tryCatch(
    generate_data_multi_path(model_param, default_resp_param(),
                             default_baseline_param(), design_set,
                             stop_on_nonpd = TRUE),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  prepare_long_data(dat, design_set,
                    carryover_t1half = t1half,
                    carryover_form   = carryover_form,
                    weibull_shape    = weibull_shape)
}

## True mean of Sx at each path and visit under the DGP, for the
## correct-mean (oracle) analysis: Sx minus its true path mean.
p10_true_means <- function(model_param, design_set) {
  rp <- default_resp_param()
  bp <- default_baseline_param()
  purrr::map_dfr(seq_along(design_set), function(i) {
    s <- build_sigma_matrix(model_param, rp, bp, design_set[[i]],
                            dgp_architecture = 'mvn',
                            stop_on_nonpd = TRUE)
    tpn <- design_set[[i]]$timepoint_name
    m <- s$means
    lab <- s$labels
    mu <- vapply(tpn, function(tp) m[lab == 'BL'] -
      sum(m[lab %in% paste(tp, c('tv', 'pb', 'br'), sep = '.')]),
      numeric(1))
    tibble::tibble(path = i, timepoint_name = tpn, mu_true = unname(mu))
  })
}

p10_spec_formula <- function(spec, mean_model = 'fitted') {
  dterm <- if (spec == 'E2') 'Dbc' else 'Db'
  rhs <- switch(mean_model,
    fitted  = paste('bm + t +', dterm, '+ bm:', dterm),
    expect  = paste('bm + t +', dterm, '+ bm:', dterm, '+ De + bm:De'),
    satur   = paste('0 + pv + bm + bm:', dterm))
  if (spec == 'E3' && mean_model != 'satur') rhs <- paste(rhs, '+ L')
  as.formula(paste('Sx ~', gsub(': ', ':', rhs)))
}

p10_target <- function(spec) {
  if (spec == 'E2') c('bm:Dbc', 'Dbc:bm') else c('bm:Db', 'Db:bm')
}

p10_ctrl <- function() nlme::lmeControl(opt = 'optim', maxIter = 200,
                                        msMaxIter = 200)

## Fit an lme; if the optim optimizer fails, refit the SAME model with
## nlminb. The model is never simplified (no fallback that drops the
## corCAR1 term). Returns NULL if both fail.
p10_lme <- function(form, data, random = ~1 | ptID, correlation = NULL,
                    method = 'REML') {
  fit <- tryCatch(
    nlme::lme(form, random = random, correlation = correlation,
              data = data, method = method, control = p10_ctrl()),
    error = function(e) NULL)
  if (!is.null(fit)) return(fit)
  tryCatch(
    nlme::lme(form, random = random, correlation = correlation,
              data = data, method = method,
              control = nlme::lmeControl(opt = 'nlminb', maxIter = 200,
                                         msMaxIter = 200)),
    error = function(e) NULL)
}

## Binomial summaries used by every driver.
p10_mcse <- function(p, n) sqrt(p * (1 - p) / n)
