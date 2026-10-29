## analysis/scripts/carryover-sensitivity/diagnose-null-type1-cr2.R
##
## Paper 10 (analysis/report/10-interaction-test-calibration/),
## Section "A cluster-robust standard error restores the size".
## At the primary null cell, each replicate is analyzed by
##   - M3 (lme, random intercept + corCAR1) for E1, E2, E3, tested
##     with the model-based SE and with the CR2 cluster-robust SE
##     (clustered on participant, Satterthwaite df; clubSandwich).
##     In paper 02's codes these are G1/G6 (E1), G3/G8 (E2).
##   - G9: E2 with the analysis half-life chosen by ML AIC over
##     {0.25, 0.5, 1, 2} weeks, REML refit, CR2 test (paper 02's
##     recommended analysis; simulation-core.R fit_spec_s7('E7cr2')).
##   - F1: bm + t + Dbc + bm:Dbc + De + bm:De at a fixed analysis
##     half-life of 0.2 week, CR2 test (simulation-core.R
##     fit_f1_cr2()).
## Failures of the CR2 computation are recorded, not dropped.
##
## Corrections of 2026-10-29 (referee report 2026-10-28): refined DGP
## (paper10-common.R), G9 and F1 arms added, CR2 failures recorded.
##
## The seed matches diagnose-null-type1.R, so the datasets are
## identical and every comparison can be paired.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1-cr2.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(nlme)
  library(clubSandwich)
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
N_REPS <- as.integer(arg_val('--reps', '3000'))
C_BM   <- as.numeric(arg_val('--cbm', '0'))
SEED   <- as.integer(arg_val('--seed', '20260520'))
OUT    <- arg_val('--out', 'diag-null-type1-cr2.rds')

cell <- list(
  c_bm = C_BM, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', N = 70, rho = 0.7,
  dgp_arch = 'mvn')

design_set  <- build_design_set(cell$design)
model_param <- p10_model_param(cell$N, cell$c_bm, cell$t1half, cell$rho)

na_row <- function(spec) tibble(
  spec = spec, estimate = NA_real_, mod_se = NA_real_,
  mod_p = NA_real_, cr2_se = NA_real_, cr2_p = NA_real_,
  cr2_df = NA_real_, best_t_half = NA_real_, fit_ok = FALSE,
  cr2_ok = FALSE)

cr2_full <- function(fit, dl, targets) {
  ct <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2', cluster = dl$ptID,
                            test = 'Satterthwaite'),
    error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  cf <- if (!is.null(ct$Coef)) as.character(ct$Coef) else rownames(ct)
  idx <- which(cf %in% targets)
  if (length(idx) == 0) return(NULL)
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
  fit <- p10_lme(p10_spec_formula(spec), dl,
                 correlation = nlme::corCAR1(form = ~t | ptID))
  if (is.null(fit)) return(na_row(spec))
  cc <- summary(fit)$tTable
  tgt <- intersect(p10_target(spec), rownames(cc))
  if (length(tgt) == 0) return(na_row(spec))
  tgt <- tgt[1]
  cr2 <- cr2_full(fit, dl, tgt)
  tibble(spec = spec, estimate = cc[tgt, 'Value'],
         mod_se = cc[tgt, 'Std.Error'], mod_p = cc[tgt, 'p-value'],
         cr2_se = if (is.null(cr2)) NA_real_ else cr2$se,
         cr2_p  = if (is.null(cr2)) NA_real_ else cr2$p,
         cr2_df = if (is.null(cr2)) NA_real_ else cr2$df,
         best_t_half = NA_real_, fit_ok = TRUE,
         cr2_ok = !is.null(cr2) && !is.na(cr2$p))
}

fit_g9 <- function(dl) {
  best <- select_half_life_aic(dl)
  if (is.null(best)) return(na_row('G9'))
  cr2 <- cr2_full(best$fit, best$dl, c('bm:Dbc', 'Dbc:bm'))
  cc <- summary(best$fit)$tTable
  tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc))[1]
  tibble(spec = 'G9', estimate = cc[tgt, 'Value'],
         mod_se = cc[tgt, 'Std.Error'], mod_p = cc[tgt, 'p-value'],
         cr2_se = if (is.null(cr2)) NA_real_ else cr2$se,
         cr2_p  = if (is.null(cr2)) NA_real_ else cr2$p,
         cr2_df = if (is.null(cr2)) NA_real_ else cr2$df,
         best_t_half = best$hl, fit_ok = TRUE,
         cr2_ok = !is.null(cr2) && !is.na(cr2$p))
}

fit_f1 <- function(dl) {
  dd <- with_dbc_f(dl, 0.2)
  fit <- p10_lme(FM1_F, dd, correlation = nlme::corCAR1(form = ~t | ptID))
  if (is.null(fit)) return(na_row('F1'))
  cc <- summary(fit)$tTable
  tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(cc))[1]
  cr2 <- cr2_full(fit, dd, tgt)
  tibble(spec = 'F1', estimate = cc[tgt, 'Value'],
         mod_se = cc[tgt, 'Std.Error'], mod_p = cc[tgt, 'p-value'],
         cr2_se = if (is.null(cr2)) NA_real_ else cr2$se,
         cr2_p  = if (is.null(cr2)) NA_real_ else cr2$p,
         cr2_df = if (is.null(cr2)) NA_real_ else cr2$df,
         best_t_half = 0.2, fit_ok = TRUE,
         cr2_ok = !is.null(cr2) && !is.na(cr2$p))
}

one_rep <- function(i) {
  dl <- p10_gen_long(model_param, design_set, cell$t1half)
  if (is.null(dl)) {
    return(map_dfr(c('E1', 'E2', 'E3', 'G9', 'F1'), na_row) |>
             mutate(rep = i, gen_ok = FALSE, .before = 1))
  }
  bind_rows(map_dfr(c('E1', 'E2', 'E3'), function(sp) fit_M3_cr2(dl, sp)),
            fit_g9(dl), fit_f1(dl)) |>
    mutate(rep = i, gen_ok = TRUE, .before = 1)
}

cat(sprintf('CR2 run: %s, N=%d, c_bm=%.2f, c.cfct=%.1f, %d reps\n',
            cell$design, cell$N, cell$c_bm, P10_CFCT, N_REPS))

plan(multicore, workers = max(1, parallelly::availableCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

summ <- res |>
  filter(fit_ok, cr2_ok) |>
  group_by(spec) |>
  summarise(
    n            = n(),
    emp_se       = sd(estimate),
    mod_se_ratio = mean(mod_se) / sd(estimate),
    cr2_se_ratio = mean(cr2_se) / sd(estimate),
    cv_cr2_se    = sd(cr2_se) / mean(cr2_se),
    rej_mod      = mean(mod_p < 0.05),
    rej_cr2      = mean(cr2_p < 0.05),
    cr2_df_q05   = quantile(cr2_df, 0.05),
    cr2_df_q95   = quantile(cr2_df, 0.95),
    cr2_df       = mean(cr2_df),
    .groups = 'drop') |>
  mutate(mcse_mod = p10_mcse(rej_mod, n), mcse_cr2 = p10_mcse(rej_cr2, n))

fails <- res |>
  group_by(spec) |>
  summarise(n = n(), gen_fail = sum(!gen_ok), fit_fail = sum(gen_ok & !fit_ok),
            cr2_fail = sum(fit_ok & !cr2_ok), .groups = 'drop')

print(as.data.frame(summ), digits = 4)
print(as.data.frame(fails))

out_dir <- p10_out_dir(repo_root)
saveRDS(list(results = res, summary = summ, failures = fails, cell = cell,
             n_reps = N_REPS, seed = SEED, dgp = p10_dgp_meta(repo_root),
             elapsed_secs = elapsed),
        file.path(out_dir, OUT))
cat(sprintf('Wrote %s\n', file.path(out_dir, OUT)))
