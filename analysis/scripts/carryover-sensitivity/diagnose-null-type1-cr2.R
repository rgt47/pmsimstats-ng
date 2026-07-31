## analysis/scripts/carryover-sensitivity/diagnose-null-type1-cr2.R
##
## Step 1 of the remediation sequence in type1-examination-plan.md
## (Option A). Stages 1 to 6 established that the production model
## M3 (lme with a random intercept and a corCAR1 residual)
## over-estimates the interaction standard error by 6 to 10 percent,
## and that the bias is asymptotic.
##
## This script verifies the proposed fix: replace the model-based
## standard error of the interaction with a CR2 bias-reduced
## cluster-robust standard error (clustered on ptID), computed with
## clubSandwich::coef_test(). A sandwich estimator is consistent for
## the true sampling variance irrespective of working-covariance
## mis-specification, so the prediction is se_ratio -> 1 and
## type-I -> 0.05.
##
## Primary null cell: c_bm = 0, exponential DGP, Hybrid design,
## N = 70, t1half = 1.0, rho = 0.7. The RNG seed matches
## diagnose-null-type1.R so the datasets are identical.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1-cr2.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(ggplot2)
  library(nlme)
  library(clubSandwich)
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

na_row <- function(spec) tibble(
  spec = spec, estimate = NA_real_, mod_se = NA_real_,
  mod_p = NA_real_, cr2_se = NA_real_, cr2_p = NA_real_,
  cr2_df = NA_real_, converged = FALSE)

## Pull the CR2 row for the target coefficient, robust to
## clubSandwich column-naming across versions.
cr2_extract <- function(fit, tgt, cluster_vec) {
  ct <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2', cluster = cluster_vec,
                            test = 'Satterthwaite'),
    error = function(e) NULL)
  if (is.null(ct)) return(list(se = NA_real_, p = NA_real_,
                               df = NA_real_))
  cf <- if (!is.null(ct$Coef)) as.character(ct$Coef) else rownames(ct)
  idx <- which(cf == tgt)
  if (length(idx) == 0) return(list(se = NA_real_, p = NA_real_,
                                    df = NA_real_))
  row <- ct[idx[1], , drop = FALSE]
  grab <- function(patterns) {
    for (p in patterns) {
      hit <- grep(p, names(row), value = TRUE, ignore.case = TRUE)
      if (length(hit)) return(as.numeric(row[[hit[1]]]))
    }
    NA_real_
  }
  list(se = grab('^SE$'),
       p  = grab(c('^p_Satt', '^p_val', '^p\\.', '^p$', '^p')),
       df = grab(c('df_Satt', '^df')))
}

fit_M3_cr2 <- function(dl, spec) {
  fit <- tryCatch(
    nlme::lme(spec_formula(spec), random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row(spec))

  cc <- summary(fit)$tTable
  tgt <- intersect(target_names(spec), rownames(cc))
  if (length(tgt) == 0) return(na_row(spec))
  tgt <- tgt[1]

  cr2 <- cr2_extract(fit, tgt, dl$ptID)

  tibble(
    spec      = spec,
    estimate  = cc[tgt, 'Value'],
    mod_se    = cc[tgt, 'Std.Error'],
    mod_p     = cc[tgt, 'p-value'],
    cr2_se    = cr2$se,
    cr2_p     = cr2$p,
    cr2_df    = cr2$df,
    converged = TRUE)
}

one_rep <- function(i) {
  dl <- gen_long()
  if (is.null(dl)) return(NULL)
  map_dfr(c('E1', 'E2', 'E3'),
          function(sp) fit_M3_cr2(dl, sp)) |>
    mutate(rep = i, .before = 1)
}

cat(sprintf('CR2 verification: %s, N=%d, c_bm=0, exp DGP, %d reps\n',
            cell$design, cell$N, N_REPS))

plan(multicore, workers = max(1, parallel::detectCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

summ <- res |>
  filter(converged, !is.na(cr2_p)) |>
  group_by(spec) |>
  summarise(
    n            = n(),
    emp_se       = sd(estimate),
    mod_se       = mean(mod_se),
    cr2_se       = mean(cr2_se),
    mod_se_ratio = mean(mod_se) / sd(estimate),
    cr2_se_ratio = mean(cr2_se) / sd(estimate),
    type1_mod    = mean(mod_p < 0.05),
    type1_cr2    = mean(cr2_p < 0.05),
    cr2_df       = mean(cr2_df),
    .groups = 'drop')

mcse <- sqrt(0.05 * 0.95 / unique(summ$n)[1])

cat('=== CR2 verification: model-based vs cluster-robust SE ===\n')
cat('mod_se_ratio : production SE / empirical SE (expect ~1.07-1.10)\n')
cat('cr2_se_ratio : CR2 robust SE / empirical SE (expect ~1.00)\n')
cat('type1_mod    : type-I, production Wald test (expect ~0.03)\n')
cat('type1_cr2    : type-I, CR2 Satterthwaite test (expect ~0.05)\n')
cat(sprintf('MC SE on a type-I estimate: ~%.4f\n\n', mcse))
print(as.data.frame(summ), digits = 4)

verdict <- all(abs(summ$cr2_se_ratio - 1) < 0.03) &&
           all(abs(summ$type1_cr2 - 0.05) < 3 * mcse)
cat(sprintf('\nVERDICT: CR2 calibration %s\n',
    if (verdict) 'CONFIRMED (se_ratio ~1, type-I ~0.05).'
    else 'NOT confirmed - inspect table above.'))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')

ecdf_dat <- res |>
  filter(converged, !is.na(cr2_p)) |>
  select(spec, mod_p, cr2_p) |>
  pivot_longer(c(mod_p, cr2_p), names_to = 'source',
               values_to = 'p') |>
  mutate(source = recode(source,
    mod_p = 'Model-based Wald', cr2_p = 'CR2 cluster-robust'))

p_ecdf <- ggplot(ecdf_dat, aes(p, colour = source)) +
  stat_ecdf(geom = 'step', linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 'dashed', colour = 'grey50') +
  facet_wrap(~ spec) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c(
    'Model-based Wald' = '#e31a1c',
    'CR2 cluster-robust' = '#1f78b4')) +
  labs(x = 'Null p-value', y = 'Empirical CDF', colour = NULL,
       title = 'CR2 cluster-robust vs model-based null p-values',
       subtitle = sprintf(
         'Model M3, Hybrid, N=70, c_bm=0, exp DGP, %d reps',
         N_REPS)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'top')

ggsave(file.path(out_dir, 'diag-null-type1-cr2-ecdf.pdf'),
       p_ecdf, width = 7.5, height = 3.4)

saveRDS(list(results = res, summary = summ, cell = cell,
             n_reps = N_REPS, seed = SEED, elapsed_secs = elapsed),
        file.path(out_dir, 'diag-null-type1-cr2.rds'))

cat(sprintf('\nWrote %s\n',
    file.path(out_dir, 'diag-null-type1-cr2-ecdf.pdf')))
cat(sprintf('Wrote %s\n',
    file.path(out_dir, 'diag-null-type1-cr2.rds')))
