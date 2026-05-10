## Study 5 pilot: small-N identifiability and type I error
## under a single-component gaussian DGP (c.bm = 0).
##
## This is a SCALED-DOWN pilot run (n_reps = 100, single N value)
## to verify that:
##   (1) the lcmm_analysis wrapper achieves convergence at >= 90%
##       of replicates,
##   (2) per-replicate timing matches the smoke-test estimate
##       (~1.5-2 s),
##   (3) the gating-coefficient null distribution under c.bm = 0
##       is centred and has approximately nominal type I error.
##
## Production Study 5 (n_reps = 5000, three N values, three K
## values for class-count selection, mildly non-gaussian DGP cell,
## true-two-class DGP cell) is committed only after this pilot
## confirms feasibility.
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/03-study5-pilot.R
##
## Pilot output is written to:
##   analysis/scripts/latent-class-mixture-application/output/
##     03-study5-pilot.rds

library(devtools)
load_all('.', quiet = TRUE)
library(data.table)

source('analysis/scripts/latent-class-mixture-application/01-lcmm-wrapper.R')

OUT_DIR <- file.path('analysis', 'scripts',
                     'latent-class-mixture-application', 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## Pilot configuration -----------------------------------------------
N_REPS <- as.integer(Sys.getenv('NREPS', '100'))
N_PT   <- as.integer(Sys.getenv('N',     '70'))
SEED_BASE <- 20260507L

## Trial design (open-label, eight weekly timepoints).
td <- buildtrialdesign(
  name_longform  = 'open label',
  name_shortform = 'OL',
  timepoints     = cumulative(rep(2.5, 8)),
  timeptnames    = paste0('OL', 1:8),
  expectancies   = rep(1, 8),
  ondrug         = list(pathA = rep(1, 8))
)

mp_null <- data.table(
  N = N_PT, c.bm = 0,                ## NULL: no biomarker moderation
  carryover_t1half = 1.0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.1, c.cfct = 0.05
)
rp <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)
bp <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)

op <- list(
  useDE = TRUE, t_random_slope = FALSE, full_model_out = FALSE,
  carryover_t1half = 1.0, simplecarryover = FALSE,
  carryover_scalefactor = 1
)

## Worker for one replicate. Returns a one-row data.table with
## both lme and lcmm results plus per-replicate timing.
run_one <- function(rep_idx) {
  set.seed(SEED_BASE + rep_idx)
  dat <- generateData(
    modelparam = mp_null, respparam = rp, blparam = bp,
    trialdesign = td$trialpaths[[1]],
    empirical = FALSE, makePositiveDefinite = TRUE,
    dgp_architecture = 'mvn'
  )
  dat[, path := 1]

  t_lme_start <- Sys.time()
  out_lme <- tryCatch(
    lme_analysis(td$trialpaths, dat, op),
    error = function(e) data.table(beta = NA_real_, betaSE = NA_real_,
                                   p = NA_real_, issingular = NA,
                                   warning = paste0('err: ',
                                                    conditionMessage(e)))
  )
  t_lme <- as.numeric(difftime(Sys.time(), t_lme_start, units = 'secs'))

  t_lcmm_start <- Sys.time()
  out_lcmm <- tryCatch(
    lcmm_analysis(td$trialpaths, dat, op),
    error = function(e) data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      beta_within = NA_real_, betaSE_within = NA_real_,
      p_within = NA_real_,
      delta_bic = NA_real_, lrt_stat = NA_real_,
      lrt_df = NA_integer_,
      bic_ng1 = NA_real_, bic_ng2 = NA_real_,
      issingular = NA, warning = paste0('err: ',
                                        conditionMessage(e)),
      entropy = NA_real_, bic = NA_real_, conv = 0L,
      conv_ng1 = 0L
    )
  )
  t_lcmm <- as.numeric(difftime(Sys.time(), t_lcmm_start,
                                 units = 'secs'))

  data.table(
    rep_idx = rep_idx,
    lme_beta   = out_lme$beta,
    lme_betaSE = out_lme$betaSE,
    lme_p      = out_lme$p,
    lme_t      = t_lme,
    lcmm_gating_beta   = out_lcmm$beta,
    lcmm_gating_betaSE = out_lcmm$betaSE,
    lcmm_gating_p      = out_lcmm$p,
    lcmm_within_beta   = out_lcmm$beta_within,
    lcmm_within_betaSE = out_lcmm$betaSE_within,
    lcmm_within_p      = out_lcmm$p_within,
    lcmm_delta_bic = out_lcmm$delta_bic,
    lcmm_lrt_stat  = out_lcmm$lrt_stat,
    lcmm_lrt_df    = out_lcmm$lrt_df,
    lcmm_bic_ng1   = out_lcmm$bic_ng1,
    lcmm_bic_ng2   = out_lcmm$bic_ng2,
    lcmm_entropy = out_lcmm$entropy,
    lcmm_conv    = out_lcmm$conv,
    lcmm_conv_ng1 = out_lcmm$conv_ng1,
    lcmm_t       = t_lcmm
  )
}

## Pilot run --------------------------------------------------------
cat(sprintf('Study 5 pilot: c.bm = 0, N = %d, n_reps = %d\n',
            N_PT, N_REPS))
cat(sprintf('Started at %s\n', format(Sys.time())))

t0 <- Sys.time()
results_list <- vector('list', N_REPS)
for (r in seq_len(N_REPS)) {
  results_list[[r]] <- run_one(r)
  if (r %% 25 == 0 || r == N_REPS) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
    cat(sprintf('  rep %d/%d, cumulative elapsed %.0f s\n',
                r, N_REPS, elapsed))
  }
}
results <- rbindlist(results_list)

## Diagnostics -------------------------------------------------------
cat('\n=== Pilot diagnostics ===\n')
cat(sprintf('lcmm convergence rate: %.1f%%\n',
            100 * mean(results$lcmm_conv == 1, na.rm = TRUE)))
cat(sprintf('lme convergence rate:  %.1f%%\n',
            100 * mean(!is.na(results$lme_p))))
cat(sprintf('Mean lme    fit time: %.2f s\n',
            mean(results$lme_t)))
cat(sprintf('Mean lcmm   fit time: %.2f s\n',
            mean(results$lcmm_t)))
cat(sprintf('Median lcmm fit time: %.2f s\n',
            median(results$lcmm_t)))

cat('\nNull-DGP rejection rate at alpha = 0.05',
    sprintf('(MCSE for p=0.05: %.3f):\n',
            sqrt(0.05 * 0.95 / N_REPS)))
cat(sprintf('  lme bm:Dbc Wald t-test:               %.3f\n',
            mean(results$lme_p < 0.05, na.rm = TRUE)))
cat(sprintf('  lcmm gating-on-bm Wald z-test:        %.3f\n',
            mean(results$lcmm_gating_p < 0.05, na.rm = TRUE)))
cat(sprintf('  lcmm within-class bm:Dbc Wald z-test: %.3f\n',
            mean(results$lcmm_within_p < 0.05, na.rm = TRUE)))

## Class-enumeration model comparisons.
cat('\nClass-enumeration model comparison (under null, ng=1 is correct):\n')
cat(sprintf('  P(BIC selects ng=2 over ng=1):  %.3f\n',
            mean(results$lcmm_delta_bic > 0, na.rm = TRUE)))
cat(sprintf('  P(BIC strongly prefers ng=2, dBIC>6): %.3f\n',
            mean(results$lcmm_delta_bic > 6, na.rm = TRUE)))
lrt_p <- pchisq(results$lcmm_lrt_stat,
                df = results$lcmm_lrt_df, lower.tail = FALSE)
cat(sprintf('  P(LRT chi-sq(df=%d) p < 0.05):  %.3f\n',
            median(results$lcmm_lrt_df, na.rm = TRUE),
            mean(lrt_p < 0.05, na.rm = TRUE)))

## Conditional procedures (BIC-gated).
sel_ng2 <- results$lcmm_delta_bic > 0
cat('\nConditional class-aware moderation procedures:\n')
cat(sprintf('  P(BIC selects ng=2 AND gating Wald rejects): %.3f\n',
            mean(sel_ng2 & results$lcmm_gating_p < 0.05,
                 na.rm = TRUE)))
cat(sprintf('  P(BIC selects ng=2 AND within-class rejects): %.3f\n',
            mean(sel_ng2 & results$lcmm_within_p < 0.05,
                 na.rm = TRUE)))

out_path <- file.path(OUT_DIR, '03-study5-pilot.rds')
saveRDS(list(
  config = list(N_REPS = N_REPS, N_PT = N_PT, SEED_BASE = SEED_BASE,
                started = format(t0)),
  results = results
), out_path)
cat(sprintf('\nPilot output written to %s\n', out_path))
