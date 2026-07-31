## analysis/scripts/carryover-sensitivity/diagnose-null-type1-stages46.R
##
## Stages 4 and 6 of type1-examination-plan.md. Follows the
## Stage 1 to 3 finding (diagnose-null-type1.R) that the interaction
## test is conservative because the model-reported standard error
## over-states the true sampling SD by 6 to 10 percent.
##
## Stage 4  Model ladder. Refit the primary null cell under four
##          correlation specifications and locate which structure
##          inflates the interaction SE:
##            M0  lm                          (no within-subject cor)
##            M1  lme, random ~1|ptID          (random intercept only)
##            M2  gls, corCAR1(~t|ptID)        (AR(1) residual only)
##            M3  lme, ~1|ptID + corCAR1       (production model)
##
## Stage 6  Sample-size gradient. Run the M3 null decomposition at
##          N in {35,70,140,280}. If type-I rises toward 0.05 as N
##          grows the conservatism is benign finite-sample
##          behaviour; if it stays flat it is an asymptotic
##          mis-specification.
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
  library(ggplot2)
  library(nlme)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) && i < length(args)) args[i + 1] else default
}
STAGE  <- as.integer(arg_val('--stage', '4'))
SEED   <- 20260520L

cell <- list(
  c_bm = 0, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', rho = 0.7, dgp_arch = 'mvn'
)

resp_param     <- default_resp_param()
baseline_param <- default_baseline_param()
design_set     <- build_design_set(cell$design)

base_model_param <- function(N) {
  list(
    N = N, c.bm = cell$c_bm,
    carryover_t1half = cell$t1half,
    carryover_form   = cell$carryover_form,
    weibull_shape    = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.tv = cell$rho, c.pb = cell$rho, c.br = cell$rho,
    c.cf1t = 0.2, c.cfct = 0.1)
}

gen_long <- function(N) {
  dat <- tryCatch(
    generate_data_multi_path(base_model_param(N), resp_param,
                             baseline_param, design_set),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  prepare_long_data(
    dat, design_set,
    carryover_t1half = cell$t1half,
    carryover_form   = cell$carryover_form,
    weibull_shape    = cell$weibull_shape)
}

spec_formula <- function(spec) switch(spec,
  E1 = Sx ~ bm + t + Db  + bm:Db,
  E2 = Sx ~ bm + t + Dbc + bm:Dbc,
  E3 = Sx ~ bm + t + Db  + bm:Db + L)

target_names <- function(spec)
  if (spec == 'E2') c('bm:Dbc', 'Dbc:bm') else c('bm:Db', 'Db:bm')

na_row <- function(model, spec) tibble(
  model = model, spec = spec, estimate = NA_real_,
  std_error = NA_real_, p_value = NA_real_, converged = FALSE)

pick <- function(tab, spec) {
  hit <- intersect(target_names(spec), rownames(tab))
  if (length(hit) == 0) NULL else tab[hit[1], ]
}

## ------------------------------------------------------------------
## The four correlation specifications
## ------------------------------------------------------------------
fit_M0 <- function(dl, spec) {
  fit <- tryCatch(lm(spec_formula(spec), data = dl),
                  error = function(e) NULL)
  if (is.null(fit)) return(na_row('M0', spec))
  r <- pick(summary(fit)$coefficients, spec)
  if (is.null(r)) return(na_row('M0', spec))
  tibble(model = 'M0', spec = spec, estimate = r['Estimate'],
         std_error = r['Std. Error'], p_value = r['Pr(>|t|)'],
         converged = TRUE)
}

fit_lme <- function(dl, spec, model, correlation) {
  fit <- tryCatch(
    nlme::lme(spec_formula(spec), random = ~1 | ptID,
              correlation = correlation, data = dl,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row(model, spec))
  r <- pick(summary(fit)$tTable, spec)
  if (is.null(r)) return(na_row(model, spec))
  tibble(model = model, spec = spec, estimate = r['Value'],
         std_error = r['Std.Error'], p_value = r['p-value'],
         converged = TRUE)
}

fit_M1 <- function(dl, spec) fit_lme(dl, spec, 'M1', NULL)
fit_M3 <- function(dl, spec)
  fit_lme(dl, spec, 'M3', nlme::corCAR1(form = ~t | ptID))

fit_M2 <- function(dl, spec) {
  fit <- tryCatch(
    nlme::gls(spec_formula(spec),
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl,
              control = nlme::glsControl(opt = 'optim')),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row('M2', spec))
  r <- pick(summary(fit)$tTable, spec)
  if (is.null(r)) return(na_row('M2', spec))
  tibble(model = 'M2', spec = spec, estimate = r['Value'],
         std_error = r['Std.Error'], p_value = r['p-value'],
         converged = TRUE)
}

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')

theme_paper <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'grey92', colour = NA),
        legend.position = 'top')

## ==================================================================
## Stage 4
## ==================================================================
if (STAGE == 4) {
  N_REPS <- as.integer(arg_val('--reps', '5000'))
  cat(sprintf('Stage 4 model ladder: %s, N=70, c_bm=0, %d reps\n',
              cell$design, N_REPS))

  one_rep <- function(i) {
    dl <- gen_long(70)
    if (is.null(dl)) return(NULL)
    map_dfr(c('E1', 'E2'), function(sp)
      bind_rows(fit_M0(dl, sp), fit_M1(dl, sp),
                fit_M2(dl, sp), fit_M3(dl, sp))) |>
      mutate(rep = i)
  }

  plan(multicore, workers = max(1, parallel::detectCores() - 1))
  t0 <- Sys.time()
  res <- future_map_dfr(seq_len(N_REPS), one_rep,
                        .options = furrr_options(seed = SEED))
  elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
  cat(sprintf('Completed in %.0f s.\n\n', elapsed))

  summ <- res |>
    filter(converged) |>
    group_by(spec, model) |>
    summarise(
      n_conv   = n(),
      mean_est = mean(estimate),
      emp_se   = sd(estimate),
      mod_se   = mean(std_error),
      se_ratio = mean(std_error) / sd(estimate),
      type1_05 = mean(p_value < 0.05),
      mc_se    = sqrt(mean(p_value < 0.05) *
                      (1 - mean(p_value < 0.05)) / n()),
      .groups = 'drop') |>
    mutate(model = factor(model, levels = c('M0','M1','M2','M3'))) |>
    arrange(spec, model)

  cat('=== Stage 4: model ladder ===\n')
  cat('M0 lm | M1 lme RI | M2 gls corCAR1 | M3 lme RI+corCAR1\n')
  cat('se_ratio > 1 = SE over-estimated (conservative);',
      '< 1 = anti-conservative\n\n')
  print(as.data.frame(summ), digits = 4)

  conv <- res |> group_by(spec, model) |>
    summarise(conv = mean(converged), .groups = 'drop')
  cat('\nConvergence:\n')
  print(as.data.frame(conv), digits = 4)

  ref <- tibble(metric = c('type1_05', 'se_ratio'),
                yint = c(0.05, 1.0))
  s4_long <- summ |>
    select(spec, model, type1_05, se_ratio) |>
    pivot_longer(c(type1_05, se_ratio),
                 names_to = 'metric', values_to = 'value')

  p4 <- ggplot(s4_long, aes(model, value, colour = spec,
                            group = spec)) +
    geom_hline(data = ref, aes(yintercept = yint),
               linetype = 'dashed', colour = 'grey50') +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = 'free_y') +
    scale_colour_manual(values = c(E1 = '#1f78b4', E2 = '#33a02c')) +
    labs(x = 'Correlation specification', y = NULL, colour = 'Spec',
         title = 'Stage 4: interaction-SE calibration by model',
         subtitle = sprintf(
           'Null cell: Hybrid, N=70, c_bm=0, exp DGP, %d reps',
           N_REPS)) +
    theme_paper

  ggsave(file.path(out_dir, 'diag-null-type1-ladder.pdf'),
         p4, width = 7, height = 3.6)

  saveRDS(list(results = res, summary = summ, conv = conv,
               n_reps = N_REPS, seed = SEED, elapsed_secs = elapsed),
          file.path(out_dir, 'diag-null-type1-stage4.rds'))
  cat(sprintf('\nWrote %s\n',
      file.path(out_dir, 'diag-null-type1-ladder.pdf')))
  cat(sprintf('Wrote %s\n',
      file.path(out_dir, 'diag-null-type1-stage4.rds')))
}

## ==================================================================
## Stage 6
## ==================================================================
if (STAGE == 6) {
  N_REPS <- as.integer(arg_val('--reps', '1500'))
  N_GRID <- c(35L, 70L, 140L, 280L)
  cat(sprintf('Stage 6 sample-size gradient: M3, N in {%s}, %d reps\n',
              paste(N_GRID, collapse = ','), N_REPS))

  one_rep <- function(i) {
    map_dfr(N_GRID, function(NN) {
      dl <- gen_long(NN)
      if (is.null(dl)) return(NULL)
      map_dfr(c('E1', 'E2', 'E3'),
              function(sp) fit_M3(dl, sp) |> mutate(N = NN))
    }) |> mutate(rep = i)
  }

  plan(multicore, workers = max(1, parallel::detectCores() - 1))
  t0 <- Sys.time()
  res <- future_map_dfr(seq_len(N_REPS), one_rep,
                        .options = furrr_options(seed = SEED + 1L))
  elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
  cat(sprintf('Completed in %.0f s.\n\n', elapsed))

  summ <- res |>
    filter(converged) |>
    group_by(N, spec) |>
    summarise(
      n_conv   = n(),
      emp_se   = sd(estimate),
      mod_se   = mean(std_error),
      se_ratio = mean(std_error) / sd(estimate),
      type1_05 = mean(p_value < 0.05),
      mc_se    = sqrt(mean(p_value < 0.05) *
                      (1 - mean(p_value < 0.05)) / n()),
      .groups = 'drop') |>
    arrange(spec, N)

  cat('=== Stage 6: sample-size gradient (model M3) ===\n')
  cat('type-I -> 0.05 as N grows = benign finite-sample;',
      'flat = asymptotic\n\n')
  print(as.data.frame(summ), digits = 4)

  ref <- tibble(metric = c('type1_05', 'se_ratio'),
                yint = c(0.05, 1.0))
  s6_long <- summ |>
    select(N, spec, type1_05, se_ratio) |>
    pivot_longer(c(type1_05, se_ratio),
                 names_to = 'metric', values_to = 'value')

  p6 <- ggplot(s6_long, aes(N, value, colour = spec, group = spec)) +
    geom_hline(data = ref, aes(yintercept = yint),
               linetype = 'dashed', colour = 'grey50') +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = 'free_y') +
    scale_x_continuous(breaks = N_GRID) +
    scale_colour_manual(
      values = c(E1 = '#1f78b4', E2 = '#33a02c', E3 = '#e31a1c')) +
    labs(x = 'Sample size N', y = NULL, colour = 'Spec',
         title = 'Stage 6: type-I and SE calibration versus N',
         subtitle = sprintf(
           'Production model M3, Hybrid, c_bm=0, exp DGP, %d reps',
           N_REPS)) +
    theme_paper

  ggsave(file.path(out_dir, 'diag-null-type1-ngradient.pdf'),
         p6, width = 7, height = 3.6)

  saveRDS(list(results = res, summary = summ,
               n_reps = N_REPS, seed = SEED + 1L,
               elapsed_secs = elapsed),
          file.path(out_dir, 'diag-null-type1-stage6.rds'))
  cat(sprintf('\nWrote %s\n',
      file.path(out_dir, 'diag-null-type1-ngradient.pdf')))
  cat(sprintf('Wrote %s\n',
      file.path(out_dir, 'diag-null-type1-stage6.rds')))
}
