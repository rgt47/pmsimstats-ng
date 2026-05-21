suppressPackageStartupMessages({
  library(dplyr); library(geepack); library(geesmv)
})
repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

ds <- build_design_set('Hybrid')
rp <- default_resp_param(); bp <- default_baseline_param()
mp <- function(cbm) list(N = 70, c.bm = cbm, carryover_t1half = 1.0,
  carryover_form = 'exponential', weibull_shape = 1,
  dgp_architecture = 'mvn', c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.2, c.cfct = 0.1)

gee_exch <- function(dl) {
  dl <- dl[order(dl$ptID, dl$t), ]
  dl$id <- as.integer(factor(dl$ptID))
  dl <- as.data.frame(dl)
  fit <- tryCatch(geepack::geeglm(Sx ~ bm + t + Dbc + bm:Dbc,
    family = gaussian, data = dl, id = id,
    corstr = 'exchangeable'), error = function(e) NULL)
  if (is.null(fit)) return(c(est = NA, se = NA, p = NA))
  sm <- summary(fit)$coefficients
  tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(sm))[1]
  est <- sm[tgt, 'Estimate']
  md <- tryCatch({
    utils::capture.output(m <- geesmv::GEE.var.md(
      Sx ~ bm + t + Dbc + bm:Dbc, id = 'id', family = gaussian,
      data = dl, corstr = 'exchangeable'))
    m
  }, error = function(e) NULL)
  if (is.null(md)) return(c(est = est, se = NA, p = NA))
  idx <- match(tgt, names(coef(fit)))
  se  <- sqrt(md$cov.beta[idx])
  df  <- length(unique(dl$id)) - length(coef(fit))
  c(est = est, se = se, p = 2 * pt(-abs(est / se), df))
}

run <- function(cbm, dropout, tag, N = 120) {
  set.seed(20260520)
  res <- t(vapply(seq_len(N), function(i) {
    dat <- generate_data_multi_path(mp(cbm), rp, bp, ds)
    if (dropout > 0)
      dat <- apply_dropout(dat, ds, rate = dropout,
                           mechanism = 'MCAR')
    gee_exch(prepare_long_data(dat, ds, carryover_t1half = 1.0,
      carryover_form = 'exponential', weibull_shape = 1))
  }, c(est = 0, se = 0, p = 0)))
  cat(sprintf('=== %s (%d reps) ===\n', tag, N))
  cat(sprintf('  SE: median=%.3f mean=%.3f max=%.3f | emp_se=%.3f',
      median(res[, 'se'], na.rm = TRUE),
      mean(res[, 'se'], na.rm = TRUE),
      max(res[, 'se'], na.rm = TRUE),
      sd(res[, 'est'], na.rm = TRUE)))
  cat(sprintf(' | rej=%.3f NA=%.3f\n',
      mean(res[, 'p'] < 0.05, na.rm = TRUE),
      mean(is.na(res[, 'p']))))
}

run(0.45, 0,   'exchangeable, ref cell c_bm=0.45')
run(0.00, 0,   'exchangeable, null c_bm=0')
run(0.45, 0.3, 'exchangeable, dropout 0.3 MCAR')
