## analysis/scripts/carryover-sensitivity/diagnose-null-type1-kr.R
##
## Empirical Kenward-Roger check at the primary null cell. The
## manuscript (Sections 4.3 and 5.2) argues that a Kenward-Roger
## finite-sample correction is the wrong instrument for the observed
## conservatism, because that conservatism is asymptotic rather than a
## finite-sample artefact. This script supplies the empirical leg of
## that argument.
##
## Tooling note. Kenward-Roger is not implemented for the production
## nlme::lme model (random intercept + corCAR1) in any standard R
## package: pbkrtest applies only to lme4 objects, which cannot carry a
## corCAR1 residual. We therefore apply KR through the mmrm package to
## the faithful marginal representation of the implicated autoregressive
## residual layer: a continuous-time spatial-exponential (sp_exp)
## covariance is the marginal equivalent of a corCAR1 residual
## (manuscript ladder rung M2). An unstructured (us) covariance is
## fitted alongside as the correctly specified comparator. For each we
## report both Kenward-Roger and Satterthwaite degrees of freedom; the
## Satterthwaite fit uses the unadjusted asymptotic covariance and so
## serves as the model-based baseline, while the Kenward-Roger fit adds
## the finite-sample variance inflation.
##
## Primary null cell and RNG seed match diagnose-null-type1-cr2.R, so
## the simulated datasets are identical to the CR2 verification.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1-kr.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(mmrm)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
N_REPS <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 3000L
SEED <- 20260520L

cell <- list(
  c_bm = 0, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', N = 70, rho = 0.7,
  dgp_arch = 'mvn')

resp_param     <- default_resp_param()
baseline_param <- default_baseline_param()
design_set     <- build_design_set(cell$design)

model_param <- list(
  N = cell$N, c.bm = cell$c_bm,
  carryover_t1half = cell$t1half,
  carryover_form   = cell$carryover_form,
  weibull_shape    = cell$weibull_shape,
  dgp_architecture = cell$dgp_arch,
  c.tv = cell$rho, c.pb = cell$rho, c.br = cell$rho,
  c.cf1t = 0.2, c.cfct = 0.1)

gen_long <- function() {
  dat <- tryCatch(
    generate_data_multi_path(model_param, resp_param,
                             baseline_param, design_set),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dl <- prepare_long_data(
    dat, design_set,
    carryover_t1half = cell$t1half,
    carryover_form   = cell$carryover_form,
    weibull_shape    = cell$weibull_shape)
  dl$visit <- factor(dl$t)
  dl
}

spec_formula <- function(spec, cov_term) {
  rhs <- switch(spec,
    A1 = 'bm + t + Db  + bm:Db',
    A2 = 'bm + t + Dbc + bm:Dbc',
    A3 = 'bm + t + Db  + bm:Db + L')
  as.formula(paste('Sx ~', rhs, '+', cov_term))
}

target_names <- function(spec)
  if (spec == 'A2') c('bm:Dbc', 'Dbc:bm') else c('bm:Db', 'Db:bm')

cov_term <- c(sp_exp = 'sp_exp(t | ptID)', us = 'us(visit | ptID)')

fit_one <- function(dl, spec, cov_label, method) {
  f <- tryCatch(
    mmrm(spec_formula(spec, cov_term[[cov_label]]), data = dl,
         method = method, reml = TRUE),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- tryCatch(summary(f)$coefficients, error = function(e) NULL)
  if (is.null(s)) return(NULL)
  tgt <- intersect(target_names(spec), rownames(s))
  if (length(tgt) == 0) return(NULL)
  tgt <- tgt[1]
  tibble(
    spec = spec, cov = cov_label, method = method,
    estimate = s[tgt, 'Estimate'], se = s[tgt, 'Std. Error'],
    df = s[tgt, 'df'], p = s[tgt, 'Pr(>|t|)'])
}

grid <- expand.grid(
  spec   = c('A1', 'A2', 'A3'),
  cov    = c('sp_exp', 'us'),
  method = c('Kenward-Roger', 'Satterthwaite'),
  stringsAsFactors = FALSE)

one_rep <- function(i) {
  dl <- gen_long()
  if (is.null(dl)) return(NULL)
  pmap_dfr(grid, function(spec, cov, method)
    fit_one(dl, spec, cov, method)) |>
    mutate(rep = i, .before = 1)
}

cat(sprintf('KR check: %s, N=%d, c_bm=0, exp DGP, %d reps\n',
            cell$design, cell$N, N_REPS))

plan(multicore, workers = max(1, parallel::detectCores() - 1))
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
    .groups  = 'drop')

mcse <- sqrt(0.05 * 0.95 / max(summ$n))

cat('=== Kenward-Roger vs Satterthwaite at the primary null cell ===\n')
cat('cov sp_exp : continuous AR(1) residual (corCAR1 layer, rung M2)\n')
cat('cov us     : unstructured (correctly specified comparator)\n')
cat('se_ratio   : mean model SE / empirical SD of the estimate\n')
cat('type1      : realised type-I error at alpha = 0.05\n')
cat(sprintf('MC SE on a type-I estimate near 0.05: ~%.4f\n\n', mcse))
print(as.data.frame(summ), digits = 4)

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
saveRDS(list(results = res, summary = summ, cell = cell,
             n_reps = N_REPS, seed = SEED, elapsed_secs = elapsed),
        file.path(out_dir, 'diag-null-type1-kr.rds'))

cat(sprintf('\nWrote %s\n',
    file.path(out_dir, 'diag-null-type1-kr.rds')))
