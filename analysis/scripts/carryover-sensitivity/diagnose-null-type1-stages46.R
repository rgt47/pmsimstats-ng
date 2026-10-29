## analysis/scripts/carryover-sensitivity/diagnose-null-type1-stages46.R
##
## Paper 10 (analysis/report/10-interaction-test-calibration/),
## stages 4 and 6 of type1-examination-plan.md.
##
## Stage 4  Two-factor ladder at the primary null cell. The same
##          datasets are refitted under a mean model crossed with a
##          working covariance, for specifications E1 and E2:
##            mean   fitted : bm + t + D + bm:D (the E specification)
##                   expect : fitted plus De + bm:De (the F1 mean
##                            structure at the E drug term)
##                   visit  : factor(t) + bm + D + bm:D (a free mean
##                            per visit, common to all paths)
##                   satur  : a free mean per path-by-visit cell
##                            (0 + pv) + bm + bm:D
##            cov    M0 lm (independence)
##                   M1 lme, random intercept
##                   M2 gls, corCAR1(~t | ptID)
##                   M3 lme, random intercept + corCAR1 (production)
##          M0 and M2 are fitted with the fitted mean only; M1 and M3
##          with all four means.
##
## Stage 6  Sample-size gradient. M3 with the fitted mean at
##          N in {35, 70, 140, 280}.
##
## Corrections of 2026-10-29 (referee report 2026-10-28): refined DGP
## (paper10-common.R), mean-model rungs added, M2 gls fit records its
## REML log-likelihood and optimizer convergence (the 2026-05-20 gls
## M2 rung disagreed with mmrm sp_exp by about 2 percent in mean SE).
##
## Usage:
##   Rscript diagnose-null-type1-stages46.R --stage 4 [--reps 5000]
##   Rscript diagnose-null-type1-stages46.R --stage 6 [--reps 1500]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(nlme)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/paper10-common.R'))

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) && i < length(args)) args[i + 1] else default
}
STAGE <- as.integer(arg_val('--stage', '4'))
SEED  <- 20260520L

cell <- list(
  c_bm = 0, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', rho = 0.7, dgp_arch = 'mvn')

design_set <- build_design_set(cell$design)

gen_long <- function(N) {
  dl <- p10_gen_long(p10_model_param(N, cell$c_bm, cell$t1half,
                                     cell$rho),
                     design_set, cell$t1half)
  if (is.null(dl)) return(NULL)
  dl$pv <- factor(paste(sub('_.*', '', dl$ptID), dl$timepoint_name,
                        sep = ':'))
  dl$visit <- factor(dl$t)
  dl
}

formula_for <- function(spec, mean_model) {
  if (mean_model == 'visit') {
    dterm <- if (spec == 'E2') 'Dbc' else 'Db'
    return(as.formula(sprintf('Sx ~ 0 + visit + bm + %s + bm:%s',
                              dterm, dterm)))
  }
  p10_spec_formula(spec, mean_model)
}

na_row <- function(model, spec, mean_model) tibble(
  model = model, spec = spec, mean_model = mean_model,
  estimate = NA_real_, std_error = NA_real_, p_value = NA_real_,
  loglik = NA_real_, converged = FALSE)

pick <- function(tab, spec) {
  hit <- intersect(p10_target(spec), rownames(tab))
  if (length(hit) == 0) NULL else tab[hit[1], ]
}

fit_M0 <- function(dl, spec, mean_model) {
  fit <- tryCatch(lm(formula_for(spec, mean_model), data = dl),
                  error = function(e) NULL)
  if (is.null(fit)) return(na_row('M0', spec, mean_model))
  r <- pick(summary(fit)$coefficients, spec)
  if (is.null(r)) return(na_row('M0', spec, mean_model))
  tibble(model = 'M0', spec = spec, mean_model = mean_model,
         estimate = r['Estimate'], std_error = r['Std. Error'],
         p_value = r['Pr(>|t|)'], loglik = NA_real_, converged = TRUE)
}

fit_lme <- function(dl, spec, mean_model, model, correlation) {
  fit <- p10_lme(formula_for(spec, mean_model), dl,
                 correlation = correlation)
  if (is.null(fit)) return(na_row(model, spec, mean_model))
  r <- pick(summary(fit)$tTable, spec)
  if (is.null(r)) return(na_row(model, spec, mean_model))
  tibble(model = model, spec = spec, mean_model = mean_model,
         estimate = r['Value'], std_error = r['Std.Error'],
         p_value = r['p-value'], loglik = as.numeric(logLik(fit)),
         converged = TRUE)
}

fit_M1 <- function(dl, spec, mean_model)
  fit_lme(dl, spec, mean_model, 'M1', NULL)
fit_M3 <- function(dl, spec, mean_model)
  fit_lme(dl, spec, mean_model, 'M3', nlme::corCAR1(form = ~t | ptID))

## gls with the default nlminb optimizer and a tight tolerance; the
## 2026-05-20 version used opt = 'optim' with default tolerances.
fit_M2 <- function(dl, spec, mean_model) {
  fit <- tryCatch(
    nlme::gls(formula_for(spec, mean_model),
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl,
              control = nlme::glsControl(maxIter = 200,
                                         msMaxIter = 200,
                                         tolerance = 1e-8,
                                         msTol = 1e-10)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row('M2', spec, mean_model))
  r <- pick(summary(fit)$tTable, spec)
  if (is.null(r)) return(na_row('M2', spec, mean_model))
  tibble(model = 'M2', spec = spec, mean_model = mean_model,
         estimate = r['Value'], std_error = r['Std.Error'],
         p_value = r['p-value'], loglik = as.numeric(logLik(fit)),
         converged = TRUE)
}

summarise_ladder <- function(res, ...) {
  res |>
    filter(converged) |>
    group_by(...) |>
    summarise(
      n_conv   = n(),
      mean_est = mean(estimate),
      emp_se   = sd(estimate),
      mod_se   = mean(std_error),
      se_ratio = mean(std_error) / sd(estimate),
      type1_05 = mean(p_value < 0.05),
      .groups = 'drop') |>
    mutate(mc_se = p10_mcse(type1_05, n_conv),
           mcse_kappa = se_ratio / sqrt(2 * (n_conv - 1)))
}

out_dir <- p10_out_dir(repo_root)

if (STAGE == 4) {
  N_REPS <- as.integer(arg_val('--reps', '5000'))
  cat(sprintf('Stage 4 two-factor ladder: N=70, c_bm=0, %d reps\n',
              N_REPS))
  one_rep <- function(i) {
    dl <- gen_long(70)
    if (is.null(dl)) return(tibble(rep = i, gen_ok = FALSE))
    map_dfr(c('E1', 'E2'), function(sp) bind_rows(
      fit_M0(dl, sp, 'fitted'),
      fit_M2(dl, sp, 'fitted'),
      map_dfr(c('fitted', 'expect', 'visit', 'satur'), function(mm)
        bind_rows(fit_M1(dl, sp, mm), fit_M3(dl, sp, mm))))) |>
      mutate(rep = i, gen_ok = TRUE, .before = 1)
  }
  plan(multicore, workers = max(1, parallelly::availableCores() - 1))
  t0 <- Sys.time()
  res <- future_map_dfr(seq_len(N_REPS), one_rep,
                        .options = furrr_options(seed = SEED))
  elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
  cat(sprintf('Completed in %.0f s.\n\n', elapsed))
  summ <- summarise_ladder(res, spec, mean_model, model)
  conv <- res |> filter(gen_ok) |>
    group_by(spec, mean_model, model) |>
    summarise(n = n(), conv = mean(converged), .groups = 'drop')
  print(as.data.frame(summ), digits = 4)
  saveRDS(list(results = res, summary = summ, conv = conv,
               n_reps = N_REPS, seed = SEED, cell = cell,
               dgp = p10_dgp_meta(repo_root), elapsed_secs = elapsed),
          file.path(out_dir, 'diag-null-type1-stage4.rds'))
  cat(sprintf('Wrote %s\n',
      file.path(out_dir, 'diag-null-type1-stage4.rds')))
}

if (STAGE == 6) {
  N_REPS <- as.integer(arg_val('--reps', '1500'))
  N_GRID <- c(35L, 70L, 140L, 280L)
  cat(sprintf('Stage 6 sample-size gradient: M3, N in {%s}, %d reps\n',
              paste(N_GRID, collapse = ','), N_REPS))
  one_rep <- function(i) {
    map_dfr(N_GRID, function(NN) {
      dl <- gen_long(NN)
      if (is.null(dl)) return(tibble(N = NN, gen_ok = FALSE))
      map_dfr(c('E1', 'E2', 'E3'),
              function(sp) fit_M3(dl, sp, 'fitted') |>
                mutate(N = NN, gen_ok = TRUE))
    }) |> mutate(rep = i)
  }
  plan(multicore, workers = max(1, parallelly::availableCores() - 1))
  t0 <- Sys.time()
  res <- future_map_dfr(seq_len(N_REPS), one_rep,
                        .options = furrr_options(seed = SEED + 1L))
  elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
  cat(sprintf('Completed in %.0f s.\n\n', elapsed))
  summ <- summarise_ladder(res, N, spec) |> arrange(spec, N)
  print(as.data.frame(summ), digits = 4)
  saveRDS(list(results = res, summary = summ, n_reps = N_REPS,
               seed = SEED + 1L, cell = cell,
               dgp = p10_dgp_meta(repo_root), elapsed_secs = elapsed),
          file.path(out_dir, 'diag-null-type1-stage6.rds'))
  cat(sprintf('Wrote %s\n',
      file.path(out_dir, 'diag-null-type1-stage6.rds')))
}
