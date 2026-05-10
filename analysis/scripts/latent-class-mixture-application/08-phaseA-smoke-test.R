## Phase A smoke test — verify the four DGP code paths and the
## lcmm_analysis wrapper run end-to-end at one cell each, and
## produce non-trivial output (non-zero variance in the symptom
## score, non-degenerate analysis output).
##
## Phase A success criterion (per simgo-paper03-phasing.md):
## "Every code path runs to completion at one smoke-test cell
## with valid output."
##
## Run from the project root:
##   Rscript analysis/scripts/latent-class-mixture-application/08-phaseA-smoke-test.R
##
## Output: a TRUE/FALSE pass for each DGP and a brief summary.
## Failures print full traceback for debugging.

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
})

source('analysis/scripts/latent-class-mixture-application/01-lcmm-wrapper.R')
source('analysis/scripts/latent-class-mixture-application/05-marginal-mixture-dgp.R')
source('analysis/scripts/latent-class-mixture-application/06-combined-moderation-dgp.R')
source('analysis/scripts/latent-class-mixture-application/07-true-heterogeneous-slopes-dgp.R')

set.seed(20260509)

## Build a single trial-design + parameter set for the smoke
## test, matching the conventions in 03-study5-pilot.R.
td <- buildtrialdesign(
  name_longform  = 'open label',
  name_shortform = 'OL',
  timepoints     = cumulative(rep(2.5, 8)),
  timeptnames    = paste0('OL', 1:8),
  expectancies   = rep(1, 8),
  ondrug         = list(pathA = rep(1, 8))
)
trialdesign <- td$trialpaths[[1]]

## Standard reference parameters (Hendrickson 2020 calibration).
modelparam <- data.table(
  N = 35, c.bm = 0.45,
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

## ---- Test 1: marginal-mixture DGP (Study 1) ----
cat('Test 1: marginal-mixture DGP ... ')
result_mm <- tryCatch({
  dat <- generate_marginal_mixture_data(
    modelparam = modelparam, respparam = respparam,
    blparam = blparam, trialdesign = trialdesign, seed = 1
  )
  ok <- nrow(dat) == 35 &&
        'latent_class' %in% names(dat) &&
        'gating_prob' %in% names(dat) &&
        var(dat[[trialdesign$timeptnames[1]]]) > 0
  list(ok = ok, n_class1 = sum(dat$latent_class),
       n_total = nrow(dat))
}, error = function(e) {
  list(ok = FALSE, err = conditionMessage(e))
})
cat(if (result_mm$ok) 'PASS' else 'FAIL', '\n')
if (result_mm$ok) {
  cat(sprintf('   class 1 prevalence: %d / %d\n',
              result_mm$n_class1, result_mm$n_total))
} else {
  cat('   error:', result_mm$err, '\n')
}

## ---- Test 2: marginal-mixture DGP with PB-class-correlation ----
cat('Test 2: marginal-mixture DGP, PB-correlated cell ... ')
mp_corr <- copy(modelparam)
mp_corr$pb_class_correlation <- 'correlated'
mp_corr$pb_class_target <- 0.3
result_mm_corr <- tryCatch({
  dat <- generate_marginal_mixture_data(
    modelparam = mp_corr, respparam = respparam,
    blparam = blparam, trialdesign = trialdesign, seed = 2
  )
  ok <- nrow(dat) == 35 &&
        'pb_class_shift' %in% names(dat) &&
        unique(dat$pb_class_shift)[1] != 0
  list(ok = ok,
       pb_shift = unique(dat$pb_class_shift)[1])
}, error = function(e) {
  list(ok = FALSE, err = conditionMessage(e))
})
cat(if (result_mm_corr$ok) 'PASS' else 'FAIL', '\n')
if (result_mm_corr$ok) {
  cat(sprintf('   pb_class_shift = %.3f\n', result_mm_corr$pb_shift))
} else {
  cat('   error:', result_mm_corr$err, '\n')
}

## ---- Test 3: combined-moderation DGP (Study 2), alpha = 0.5 ----
cat('Test 3: combined-moderation DGP, alpha = 0.5 ... ')
mp_combined <- copy(modelparam)
mp_combined$alpha <- 0.5
result_combined <- tryCatch({
  dat <- generate_combined_moderation_data(
    modelparam = mp_combined, respparam = respparam,
    blparam = blparam, trialdesign = trialdesign, seed = 3
  )
  ok <- nrow(dat) == 35 &&
        'alpha' %in% names(dat) &&
        unique(dat$alpha) == 0.5
  list(ok = ok)
}, error = function(e) {
  list(ok = FALSE, err = conditionMessage(e))
})
cat(if (result_combined$ok) 'PASS' else 'FAIL', '\n')
if (!result_combined$ok) {
  cat('   error:', result_combined$err, '\n')
}

## ---- Test 4: true heterogeneous-random-slopes DGP (Study 3) ----
cat('Test 4: true heterogeneous-random-slopes DGP ... ')
mp_slopes <- copy(modelparam)
mp_slopes$mu_R <- 1.0
mp_slopes$mu_NR <- 0.0
mp_slopes$sigma_R <- 0.3
mp_slopes$sigma_NR <- 0.3
result_slopes <- tryCatch({
  dat <- generate_true_random_slopes_data(
    modelparam = mp_slopes, respparam = respparam,
    blparam = blparam, trialdesign = trialdesign, seed = 4
  )
  ok <- nrow(dat) == 35 &&
        'beta_slope' %in% names(dat) &&
        var(dat$beta_slope) > 0
  list(ok = ok,
       n_class1 = sum(dat$latent_class),
       slope_var = var(dat$beta_slope))
}, error = function(e) {
  list(ok = FALSE, err = conditionMessage(e))
})
cat(if (result_slopes$ok) 'PASS' else 'FAIL', '\n')
if (result_slopes$ok) {
  cat(sprintf('   beta_slope variance = %.4f\n',
              result_slopes$slope_var))
} else {
  cat('   error:', result_slopes$err, '\n')
}

## ---- Test 5: lcmm_analysis on a marginal-mixture draw ----
cat('Test 5: lcmm_analysis on marginal-mixture data ... ')
result_lcmm <- tryCatch({
  dat <- generate_marginal_mixture_data(
    modelparam = modelparam, respparam = respparam,
    blparam = blparam, trialdesign = trialdesign, seed = 5
  )
  dat[, path := 1]
  fit <- lcmm_analysis(td$trialpaths, dat,
                       op = list(carryover_t1half = 1.0))
  ok <- !is.null(fit) && nrow(fit) == 1 &&
        is.finite(fit$beta) && is.finite(fit$delta_bic) &&
        fit$conv == 1
  list(ok = ok,
       beta = if (ok) fit$beta else NA,
       p = if (ok) fit$p else NA,
       delta_bic = if (ok) fit$delta_bic else NA,
       conv = if (ok) fit$conv else NA)
}, error = function(e) {
  list(ok = FALSE, err = conditionMessage(e))
})
cat(if (result_lcmm$ok) 'PASS' else 'FAIL', '\n')
if (result_lcmm$ok) {
  cat(sprintf(
    '   beta=%.3f p=%.3f delta_bic=%.2f conv=%d\n',
    result_lcmm$beta, result_lcmm$p,
    result_lcmm$delta_bic, result_lcmm$conv))
} else {
  cat('   error:', result_lcmm$err, '\n')
}

## ---- Phase A summary ----
all_ok <- all(c(result_mm$ok, result_mm_corr$ok,
                result_combined$ok, result_slopes$ok,
                result_lcmm$ok))
cat('\n')
cat('Phase A smoke-test summary: ',
    if (all_ok) 'ALL PASS' else 'FAILURES PRESENT', '\n')
cat('Date:', format(Sys.time(), '%Y-%m-%d %H:%M %Z'), '\n')
