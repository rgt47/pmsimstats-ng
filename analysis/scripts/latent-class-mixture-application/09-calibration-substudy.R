## Cross-study effect-size calibration sub-study.
##
## Phase A final infrastructure deliverable for paper 03 path
## (b). Per the ADEMP pre-registration §"Effect-size
## calibration":
##
##   "A cross-study calibration sub-study expresses each DGP's
##   biomarker-treatment moderation parameter as a marginal
##   Cohen's f^2-style standardised effect size at a reference
##   time (week 8) and on-drug state. The calibration sub-study
##   runs once during Phase 1 and produces calibration-table.rds,
##   which subsequent cells consult to align their effect-size
##   axes."
##
## This script computes, for each of the five DGP families used
## across Studies 1-5, the marginal Cor(B, BR) and Cohen's f^2
## at the reference timepoint (week 8, on-drug) for a sweep of
## moderation-strength parameter values. The output table lets
## downstream Studies translate between native parameters
## (c.bm, gating_slope, alpha, etc.) and the standardised
## marginal-correlation effect-size axis.
##
## Scope: large-N (n = 500), single-replicate-per-cell (the
## marginal correlation is a population moment, not a sampling
## quantity, so single replicate suffices once N is large).
## Wall-clock target: ~3 minutes total for the full sweep.
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/09-calibration-substudy.R
##
## Output:
##   analysis/data/quick-sim/03-calibration-table.rds

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
})

source('analysis/scripts/latent-class-mixture-application/05-marginal-mixture-dgp.R')
source('analysis/scripts/latent-class-mixture-application/06-combined-moderation-dgp.R')
source('analysis/scripts/latent-class-mixture-application/07-true-heterogeneous-slopes-dgp.R')

OUT_DIR <- file.path('analysis', 'data', 'quick-sim')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PATH <- file.path(OUT_DIR, '03-calibration-table.rds')

set.seed(20260509)

## ---- Common reference cell ----------------------------------------
N_CALIB <- 500            ## n_subjects for calibration draws
REF_TIMEPOINT <- 'OL8'    ## reference (week 8) on-drug timepoint

td <- buildtrialdesign(
  name_longform  = 'open label',
  name_shortform = 'OL',
  timepoints     = cumulative(rep(2.5, 8)),
  timeptnames    = paste0('OL', 1:8),
  expectancies   = rep(1, 8),
  ondrug         = list(pathA = rep(1, 8))
)
trialdesign <- td$trialpaths[[1]]

base_modelparam <- data.table(
  N = N_CALIB, c.bm = 0.45,
  carryover_t1half = 1.0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.1, c.cfct = 0.05
)
respparam <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)
blparam <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)

## ---- Helper: compute marginal Cor(B, BR) and f^2 ----
##
## Cor(B, BR) at the reference timepoint is the canonical
## standardised effect-size measure for biomarker-treatment
## moderation in the pmsimstats framework. Cohen's f^2 is the
## variance-explained version: f^2 = R^2 / (1 - R^2) where
## R^2 is the coefficient of determination of BR on B.
##
## The reference timepoint is on-drug at week 8 ('OL8'). The
## column name in the long-format data is paste0(REF_TIMEPOINT,
## '.br').
calibration_metrics <- function(dat) {
  br_col <- paste0(REF_TIMEPOINT, '.br')
  if (!(br_col %in% names(dat))) {
    return(list(achieved_cor = NA_real_, achieved_f2 = NA_real_))
  }
  bm <- dat$bm
  br <- dat[[br_col]]
  if (length(bm) == 0 || var(bm) == 0 || var(br) == 0) {
    return(list(achieved_cor = NA_real_, achieved_f2 = NA_real_))
  }
  cor_BR <- cor(bm, br)
  R2 <- cor_BR^2
  f2 <- R2 / pmax(1 - R2, 1e-9)
  list(achieved_cor = cor_BR, achieved_f2 = f2)
}

## ---- Sweep 1: Architecture A (mean moderation) ----
cat('Calibrating Architecture A across c.bm values...\n')
arch_a_grid <- data.table(
  dgp_family = 'mean_moderation',
  c.bm = c(0.0, 0.10, 0.20, 0.30, 0.45, 0.60)
)
arch_a_rows <- vector('list', nrow(arch_a_grid))
for (i in seq_len(nrow(arch_a_grid))) {
  mp <- copy(base_modelparam)
  mp$c.bm <- arch_a_grid$c.bm[i]
  dat <- generateData(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, empirical = FALSE,
    makePositiveDefinite = TRUE, seed = 100 + i,
    dgp_architecture = 'mean_moderation'
  )
  m <- calibration_metrics(dat)
  arch_a_rows[[i]] <- data.table(
    dgp_family = 'mean_moderation',
    parameter = 'c.bm', parameter_value = mp$c.bm,
    achieved_cor = m$achieved_cor, achieved_f2 = m$achieved_f2
  )
}
arch_a_table <- rbindlist(arch_a_rows)

## ---- Sweep 2: Architecture B (MVN) ----
cat('Calibrating Architecture B across c.bm values...\n')
arch_b_grid <- data.table(
  c.bm = c(0.0, 0.10, 0.20, 0.30, 0.45)
)
arch_b_rows <- vector('list', nrow(arch_b_grid))
for (i in seq_len(nrow(arch_b_grid))) {
  mp <- copy(base_modelparam)
  mp$c.bm <- arch_b_grid$c.bm[i]
  dat <- generateData(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, empirical = FALSE,
    makePositiveDefinite = TRUE, seed = 200 + i,
    dgp_architecture = 'mvn'
  )
  m <- calibration_metrics(dat)
  arch_b_rows[[i]] <- data.table(
    dgp_family = 'mvn', parameter = 'c.bm',
    parameter_value = mp$c.bm,
    achieved_cor = m$achieved_cor, achieved_f2 = m$achieved_f2
  )
}
arch_b_table <- rbindlist(arch_b_rows)

## ---- Sweep 3: Marginal-mixture DGP across gating slope ----
cat('Calibrating marginal-mixture DGP across gating_slope...\n')
mm_grid <- data.table(
  gating_slope = c(0, 0.5, 1.0, 1.5, 2.0)
)
mm_rows <- vector('list', nrow(mm_grid))
for (i in seq_len(nrow(mm_grid))) {
  mp <- copy(base_modelparam)
  mp$gating_slope <- mm_grid$gating_slope[i]
  mp$beta_R_factor <- 1.0
  mp$beta_NR_factor <- 0.0
  dat <- generate_marginal_mixture_data(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, seed = 300 + i
  )
  m <- calibration_metrics(dat)
  mm_rows[[i]] <- data.table(
    dgp_family = 'marginal_mixture', parameter = 'gating_slope',
    parameter_value = mp$gating_slope,
    achieved_cor = m$achieved_cor, achieved_f2 = m$achieved_f2
  )
}
mm_table <- rbindlist(mm_rows)

## ---- Sweep 4: Combined-moderation DGP across (alpha, c.bm) ----
cat('Calibrating combined-moderation DGP across (alpha, c.bm)...\n')
combined_grid <- CJ(
  alpha = c(0.0, 0.25, 0.5, 0.75, 1.0),
  c.bm = c(0.30, 0.45)
)
combined_rows <- vector('list', nrow(combined_grid))
for (i in seq_len(nrow(combined_grid))) {
  mp <- copy(base_modelparam)
  mp$alpha <- combined_grid$alpha[i]
  mp$c.bm <- combined_grid$c.bm[i]
  ## Naive split for this calibration pass; the achieved_cor
  ## column tells downstream Studies how naive maps to actual.
  mp$c_bm_mean <- (1 - mp$alpha) * mp$c.bm
  mp$c_bm_cov <- mp$alpha * mp$c.bm
  dat <- generate_combined_moderation_data(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, seed = 400 + i
  )
  m <- calibration_metrics(dat)
  combined_rows[[i]] <- data.table(
    dgp_family = 'combined_moderation',
    parameter = 'alpha_c.bm',
    parameter_value = paste0('alpha=', mp$alpha,
                              ',c.bm=', mp$c.bm),
    alpha = mp$alpha, c_bm_target = mp$c.bm,
    c_bm_mean = mp$c_bm_mean, c_bm_cov = mp$c_bm_cov,
    achieved_cor = m$achieved_cor, achieved_f2 = m$achieved_f2
  )
}
combined_table <- rbindlist(combined_rows)

## ---- Sweep 5: True random-slopes DGP across mu_R/sigma_R ----
cat('Calibrating true-random-slopes DGP across (mu_R, sigma_R)...\n')
slopes_grid <- CJ(
  mu_R = c(0.5, 1.0, 1.5),
  sigma_R = c(0.1, 0.3, 0.5)
)
slopes_rows <- vector('list', nrow(slopes_grid))
for (i in seq_len(nrow(slopes_grid))) {
  mp <- copy(base_modelparam)
  mp$mu_R <- slopes_grid$mu_R[i]
  mp$sigma_R <- slopes_grid$sigma_R[i]
  mp$mu_NR <- 0.0
  mp$sigma_NR <- mp$sigma_R
  mp$gating_slope <- 1.0
  dat <- generate_true_random_slopes_data(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, seed = 500 + i
  )
  m <- calibration_metrics(dat)
  slopes_rows[[i]] <- data.table(
    dgp_family = 'true_random_slopes',
    parameter = 'mu_R_sigma_R',
    parameter_value = paste0('mu_R=', mp$mu_R,
                              ',sigma_R=', mp$sigma_R),
    mu_R = mp$mu_R, sigma_R = mp$sigma_R,
    achieved_cor = m$achieved_cor, achieved_f2 = m$achieved_f2
  )
}
slopes_table <- rbindlist(slopes_rows)

## ---- Combine and save ----
calib_table <- rbindlist(list(arch_a_table, arch_b_table,
                                mm_table, combined_table,
                                slopes_table),
                          fill = TRUE)

attr(calib_table, 'reference_timepoint') <- REF_TIMEPOINT
attr(calib_table, 'N_calibration') <- N_CALIB
attr(calib_table, 'date') <- format(Sys.time(),
                                       '%Y-%m-%d %H:%M %Z')

saveRDS(calib_table, OUT_PATH)

cat('\nCalibration table written to:', OUT_PATH, '\n')
cat('Rows:', nrow(calib_table), '\n')
cat('Reference timepoint:', REF_TIMEPOINT, '\n\n')
cat('Sample rows:\n')
print(head(calib_table, 10))
