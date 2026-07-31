## analysis/scripts/carryover-sensitivity/diagnose-null-type1.R
##
## Stages 1 to 3 of type1-examination-plan.md. Examines why the
## empirical type-I error of the biomarker-by-treatment interaction
## test runs below nominal 0.05 in the null (c_bm = 0) cells.
##
## Stage 1  Instrument the per-replicate output: an augmented fitter
##          records, beyond estimate/std_error/p_value, the lme
##          degrees of freedom, the corCAR1 phi, and the random-
##          intercept and residual standard deviations.
## Stage 2  Decompose the conservatism into four channels: point-
##          estimate bias, SE calibration, statistic dispersion, and
##          p-value uniformity.
## Stage 3  Separate the degrees-of-freedom channel from the SE
##          channel by recomputing p-values against a normal
##          reference.
##
## Primary null cell: c_bm = 0, exponential DGP, Hybrid design,
## N = 70, t1half = 1.0, rho = 0.7. Architecture is irrelevant at
## the null (see check-null-architecture-equivalence.R).
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(ggplot2)
  library(nlme)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
N_REPS <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 5000L
SEED <- 20260520L

## ------------------------------------------------------------------
## Primary null cell
## ------------------------------------------------------------------
cell <- list(
  c_bm = 0, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', N = 70, rho = 0.7,
  dgp_arch = 'mvn'
)

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
  c.cf1t = 0.2, c.cfct = 0.1
)

## ------------------------------------------------------------------
## Stage 1 ancillary: positive-definiteness of the null covariance.
## At c_bm = 0 the sigma is fixed per path, so a single build per
## path is definitive.
## ------------------------------------------------------------------
pd_corrected <- map_lgl(design_set, function(td) {
  build_sigma_matrix(model_param, resp_param, baseline_param, td,
                     dgp_architecture = cell$dgp_arch)$was_pd_corrected
})

## ------------------------------------------------------------------
## Stage 1: augmented fitter
## ------------------------------------------------------------------
na_spec_row <- function(spec) {
  tibble(spec = spec, estimate = NA_real_, std_error = NA_real_,
         df = NA_real_, t_value = NA_real_, p_value = NA_real_,
         phi = NA_real_, sd_ranef = NA_real_, sd_resid = NA_real_,
         has_corcar1 = NA, converged = FALSE)
}

fit_spec_diag <- function(dat_long, spec) {
  spec <- match.arg(spec, c('E1', 'E2', 'E3'))
  form <- switch(spec,
    E1 = as.formula('Sx ~ bm + t + Db  + bm:Db'),
    E2 = as.formula('Sx ~ bm + t + Dbc + bm:Dbc'),
    E3 = as.formula('Sx ~ bm + t + Db  + bm:Db + L'))

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
    })

  if (is.null(fit)) return(na_spec_row(spec))

  cc <- summary(fit)$tTable
  target <- switch(spec,
    E1 = intersect(c('bm:Db',  'Db:bm'),  rownames(cc)),
    E2 = intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc)),
    E3 = intersect(c('bm:Db',  'Db:bm'),  rownames(cc)))
  if (length(target) == 0) return(na_spec_row(spec))

  has_cor <- !is.null(fit$modelStruct$corStruct)
  phi <- if (has_cor)
    as.numeric(coef(fit$modelStruct$corStruct,
                    unconstrained = FALSE))[1] else NA_real_
  sd_ranef <- suppressWarnings(
    as.numeric(nlme::VarCorr(fit)['(Intercept)', 'StdDev']))

  tibble(
    spec        = spec,
    estimate    = cc[target[1], 'Value'],
    std_error   = cc[target[1], 'Std.Error'],
    df          = cc[target[1], 'DF'],
    t_value     = cc[target[1], 't-value'],
    p_value     = cc[target[1], 'p-value'],
    phi         = phi,
    sd_ranef    = sd_ranef,
    sd_resid    = fit$sigma,
    has_corcar1 = has_cor,
    converged   = TRUE)
}

## ------------------------------------------------------------------
## One replicate: generate null data, fit all three specs
## ------------------------------------------------------------------
one_rep <- function(i) {
  dat <- tryCatch(
    generate_data_multi_path(model_param, resp_param,
                             baseline_param, design_set),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) {
    return(bind_rows(na_spec_row('E1'), na_spec_row('E2'),
                     na_spec_row('E3')) |> mutate(rep = i, .before = 1))
  }
  dat_long <- prepare_long_data(
    dat, design_set,
    carryover_t1half = cell$t1half,
    carryover_form   = cell$carryover_form,
    weibull_shape    = cell$weibull_shape)
  bind_rows(
    fit_spec_diag(dat_long, 'E1'),
    fit_spec_diag(dat_long, 'E2'),
    fit_spec_diag(dat_long, 'E3')
  ) |> mutate(rep = i, .before = 1)
}

cat(sprintf('Primary null cell: %s, N=%d, c_bm=0, exp DGP, t1half=%.1f\n',
            cell$design, cell$N, cell$t1half))
cat(sprintf('Replicates: %d\n', N_REPS))
cat(sprintf('PD correction triggered on %d of %d %s paths.\n',
            sum(pd_corrected), length(pd_corrected), cell$design))

plan(multicore, workers = max(1, parallel::detectCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

## ------------------------------------------------------------------
## Stage 3 column: p-value under a standard-normal reference
## ------------------------------------------------------------------
res <- res |> mutate(p_value_z = 2 * pnorm(-abs(t_value)))

## ------------------------------------------------------------------
## Stage 2: four-channel decomposition
## ------------------------------------------------------------------
decomp <- res |>
  filter(converged) |>
  group_by(spec) |>
  summarise(
    n_conv    = n(),
    mean_est  = mean(estimate),
    bias_z    = mean(estimate) / (sd(estimate) / sqrt(n())),
    emp_se    = sd(estimate),
    mod_se    = mean(std_error),
    se_ratio  = mean(std_error) / sd(estimate),
    sd_tstat  = sd(t_value),
    type1_01  = mean(p_value < 0.01),
    type1_05  = mean(p_value < 0.05),
    type1_10  = mean(p_value < 0.10),
    ks_p      = suppressWarnings(
      ks.test(p_value, 'punif')$p.value),
    .groups = 'drop')

## ------------------------------------------------------------------
## Stage 3: degrees-of-freedom versus standard-error channel
## ------------------------------------------------------------------
stage3 <- res |>
  filter(converged) |>
  group_by(spec) |>
  summarise(
    mean_df    = mean(df),
    min_df     = min(df),
    max_df     = max(df),
    type1_05_t = mean(p_value   < 0.05),
    type1_05_z = mean(p_value_z < 0.05),
    .groups = 'drop')

conv_rate <- res |>
  group_by(spec) |>
  summarise(conv = mean(converged),
            corcar1 = mean(has_corcar1, na.rm = TRUE),
            .groups = 'drop')

mcse_05 <- sqrt(0.05 * 0.95 / N_REPS)

cat('=== Convergence ===\n')
print(as.data.frame(conv_rate), digits = 4)
cat(sprintf('\nMC SE on a 0.05 type-I estimate at n=%d: %.4f\n\n',
            N_REPS, mcse_05))

cat('=== Stage 2: four-channel decomposition ===\n')
cat('mean_est : mean interaction estimate (should be ~0)\n')
cat('bias_z   : mean_est / SE(mean_est); |bias_z| > ~2 flags bias\n')
cat('emp_se   : SD(estimate) across reps (the true SE)\n')
cat('mod_se   : mean model-reported std_error\n')
cat('se_ratio : mod_se / emp_se; > 1.05 flags SE over-estimation\n')
cat('sd_tstat : SD(t-statistic) under H0; < 1 flags conservatism\n')
cat('ks_p     : KS test of p-values vs Uniform(0,1)\n\n')
print(as.data.frame(decomp), digits = 4)

cat('\n=== Stage 3: DF channel vs SE channel ===\n')
cat('type1_05_t : type-I under the lme t reference (production)\n')
cat('type1_05_z : type-I under a normal reference (same stat)\n')
cat('If type1_05_z rises to ~0.05, the DF rule is the cause.\n')
cat('If mean_df is large, t and z coincide and DF is excluded.\n\n')
print(as.data.frame(stage3), digits = 4)

## ------------------------------------------------------------------
## p-value ECDF figure
## ------------------------------------------------------------------
out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')

p_ecdf <- res |>
  filter(converged) |>
  ggplot(aes(p_value, colour = spec)) +
  stat_ecdf(geom = 'step', linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dashed', colour = 'grey50') +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(
    values = c(E1 = '#1f78b4', E2 = '#33a02c', E3 = '#e31a1c')) +
  labs(x = 'Null p-value', y = 'Empirical CDF', colour = 'Spec',
       title = 'Null p-value ECDF versus Uniform(0,1)',
       subtitle = sprintf(
         'Hybrid, N=70, c_bm=0, exp DGP, t1half=1.0, %d reps',
         N_REPS)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, 'diag-null-type1-ecdf.pdf'),
       p_ecdf, width = 5.5, height = 5.2)

saveRDS(
  list(results = res, decomp = decomp, stage3 = stage3,
       conv = conv_rate, cell = cell, n_reps = N_REPS,
       seed = SEED, pd_corrected = pd_corrected,
       elapsed_secs = elapsed),
  file.path(out_dir, 'diag-null-type1.rds'))

cat(sprintf('\nWrote %s\n', file.path(out_dir, 'diag-null-type1-ecdf.pdf')))
cat(sprintf('Wrote %s\n', file.path(out_dir, 'diag-null-type1.rds')))
