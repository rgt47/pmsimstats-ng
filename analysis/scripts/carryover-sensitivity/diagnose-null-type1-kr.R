## analysis/scripts/carryover-sensitivity/diagnose-null-type1-kr.R
##
## Paper 10 (analysis/report/10-interaction-test-calibration/),
## Kenward-Roger check. Kenward-Roger is available in R for marginal
## models through the mmrm package, which has no random intercept.
## The fits here are therefore of the ladder rung M2, not of the
## production model M3:
##   sp_exp : continuous-time exponential (CAR(1)) covariance in
##            calendar weeks, the marginal model of rung M2
##   us     : unstructured 8 x 8 covariance (36 parameters)
## Each is fitted with Kenward-Roger and with Satterthwaite degrees of
## freedom (the Satterthwaite fit keeps the unadjusted model-based
## covariance and is the model-based baseline). On the same datasets
## the M2 rung is also fitted with nlme::gls, to compare the two
## engines on nominally the same model.
##
## Corrections of 2026-10-29 (referee report 2026-10-28): refined DGP
## (paper10-common.R), the fits are labeled M2, and a gls M2 fit with
## its REML log-likelihood is recorded beside the mmrm sp_exp fit.
##
## The seed matches diagnose-null-type1.R, so the datasets are
## identical to the other null-cell drivers.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1-kr.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(nlme)
  library(mmrm)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/paper10-common.R'))

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
N_REPS <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 2000L
SEED <- 20260520L

cell <- list(
  c_bm = 0, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', N = 70, rho = 0.7,
  dgp_arch = 'mvn')

design_set  <- build_design_set(cell$design)
model_param <- p10_model_param(cell$N, cell$c_bm, cell$t1half, cell$rho)

rhs_of <- function(spec) switch(spec,
  E1 = 'bm + t + Db  + bm:Db',
  E2 = 'bm + t + Dbc + bm:Dbc',
  E3 = 'bm + t + Db  + bm:Db + L')

cov_term <- c(sp_exp = 'sp_exp(t | ptID)', us = 'us(visit | ptID)')

fit_one <- function(dl, spec, cov_label, method) {
  f <- tryCatch(
    mmrm(as.formula(paste('Sx ~', rhs_of(spec), '+',
                          cov_term[[cov_label]])),
         data = dl, method = method, reml = TRUE),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- tryCatch(summary(f)$coefficients, error = function(e) NULL)
  if (is.null(s)) return(NULL)
  tgt <- intersect(p10_target(spec), rownames(s))
  if (length(tgt) == 0) return(NULL)
  tgt <- tgt[1]
  tibble(spec = spec, cov = cov_label, method = method,
         estimate = s[tgt, 'Estimate'], se = s[tgt, 'Std. Error'],
         df = s[tgt, 'df'], p = s[tgt, 'Pr(>|t|)'],
         loglik = as.numeric(logLik(f)))
}

fit_gls <- function(dl, spec) {
  f <- tryCatch(
    nlme::gls(as.formula(paste('Sx ~', rhs_of(spec))),
              correlation = nlme::corCAR1(form = ~t | ptID), data = dl,
              control = nlme::glsControl(maxIter = 200, msMaxIter = 200,
                                         tolerance = 1e-8,
                                         msTol = 1e-10)),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- summary(f)$tTable
  tgt <- intersect(p10_target(spec), rownames(s))[1]
  tibble(spec = spec, cov = 'gls_car1', method = 'model-based',
         estimate = s[tgt, 'Value'], se = s[tgt, 'Std.Error'],
         df = NA_real_, p = s[tgt, 'p-value'],
         loglik = as.numeric(logLik(f)))
}

grid <- expand_grid(
  spec   = c('E1', 'E2', 'E3'),
  cov    = c('sp_exp', 'us'),
  method = c('Kenward-Roger', 'Satterthwaite'))

one_rep <- function(i) {
  dl <- p10_gen_long(model_param, design_set, cell$t1half)
  if (is.null(dl)) return(NULL)
  dl$visit <- factor(dl$t)
  dl$ptID <- factor(dl$ptID)
  bind_rows(
    pmap_dfr(grid, function(spec, cov, method)
      fit_one(dl, spec, cov, method)),
    map_dfr(c('E1', 'E2', 'E3'), function(sp) fit_gls(dl, sp))) |>
    mutate(rep = i, .before = 1)
}

cat(sprintf('KR check (M2 rung): %s, N=%d, c_bm=0, c.cfct=%.1f, %d reps\n',
            cell$design, cell$N, P10_CFCT, N_REPS))

plan(multicore, workers = max(1, parallelly::availableCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

summ <- res |>
  filter(!is.na(p)) |>
  group_by(cov, method, spec) |>
  summarise(
    n        = n(),
    emp_se   = sd(estimate),
    mean_se  = mean(se),
    se_ratio = mean(se) / sd(estimate),
    mean_df  = mean(df),
    type1    = mean(p < 0.05),
    .groups  = 'drop') |>
  mutate(mcse = p10_mcse(type1, n))

print(as.data.frame(summ), digits = 4)

out_dir <- p10_out_dir(repo_root)
saveRDS(list(results = res, summary = summ, cell = cell,
             n_reps = N_REPS, seed = SEED, dgp = p10_dgp_meta(repo_root),
             elapsed_secs = elapsed),
        file.path(out_dir, 'diag-null-type1-kr.rds'))
cat(sprintf('Wrote %s\n', file.path(out_dir, 'diag-null-type1-kr.rds')))
