## Pilot 4: power profile of the linear-mixed baseline and the
## BIC-conditional class-aware procedure on the heterogeneous-
## random-slopes DGP.
##
## Sweeps the gating slope (the strength with which the biomarker
## predicts latent-class membership) at fixed responder / non-
## responder multipliers (1.0 / 0.0). At each gating slope we
## report the power of:
##   - lme bm:Dbc Wald t-test (the linear-mixed baseline),
##   - lcmm gating-on-bm Wald z-test (unconditional, for record),
##   - lcmm within-class bm:Dbc Wald z-test (unconditional),
##   - BIC selection rate (P(ng=2 preferred over ng=1)),
##   - BIC-conditional gating procedure (P(BIC selects ng=2 AND
##     the gating Wald rejects)),
##   - BIC-conditional within-class procedure.
##
## The Type I cell from Pilot 3 (Study 5 null DGP) is included
## here as a calibration check; its rejection rates appear in
## the gating_slope = 0 row for direct comparison.
##
## Run from the package root:
##   NREPS=100 Rscript analysis/scripts/latent-class-mixture-application/04-pilot-power-random-slopes.R
##
## Output:
##   analysis/scripts/latent-class-mixture-application/output/
##     04-pilot-power-random-slopes.rds

library(devtools)
load_all('.', quiet = TRUE)
library(data.table)

source('analysis/scripts/latent-class-mixture-application/01-lcmm-wrapper.R')
source('analysis/scripts/latent-class-mixture-application/02-heterogeneous-slopes-dgp.R')

OUT_DIR <- file.path('analysis', 'scripts',
                     'latent-class-mixture-application', 'output')
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_REPS <- as.integer(Sys.getenv('NREPS', '100'))
N_PT   <- as.integer(Sys.getenv('N',     '70'))
SEED_BASE <- 20260507L

td <- buildtrialdesign(
  name_longform  = 'open label',
  name_shortform = 'OL',
  timepoints     = cumulative(rep(2.5, 8)),
  timeptnames    = paste0('OL', 1:8),
  expectancies   = rep(1, 8),
  ondrug         = list(pathA = rep(1, 8))
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

## Sweep grid: gating-slope axis at fixed beta_R=1, beta_NR=0.
GATING_SLOPES <- c(0, 0.5, 1.0, 1.5, 2.0)

run_one_cell <- function(gating_slope, n_reps) {
  out <- vector('list', n_reps)
  for (r in seq_len(n_reps)) {
    set.seed(SEED_BASE + r + round(gating_slope * 1000))
    if (gating_slope == 0) {
      ## Calibration cell: matches the Pilot 3 null DGP (single
      ## component, c.bm = 0). Use generateData directly.
      mp_null <- data.table(
        N = N_PT, c.bm = 0, carryover_t1half = 1.0,
        c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
        c.cf1t = 0.1, c.cfct = 0.05
      )
      dat <- generateData(
        modelparam = mp_null, respparam = rp, blparam = bp,
        trialdesign = td$trialpaths[[1]],
        empirical = FALSE, makePositiveDefinite = TRUE,
        dgp_architecture = 'mvn'
      )
    } else {
      mp_rs <- data.table(
        N = N_PT, c.bm = 0, carryover_t1half = 1.0,
        c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
        c.cf1t = 0.1, c.cfct = 0.05,
        beta_R_factor = 1.0, beta_NR_factor = 0.0,
        gating_intercept = 0.0, gating_slope = gating_slope
      )
      dat <- generate_random_slopes_data(
        modelparam = mp_rs, respparam = rp, blparam = bp,
        trialdesign = td$trialpaths[[1]],
        empirical = FALSE, makePositiveDefinite = TRUE
      )
    }
    dat[, path := 1]

    out_lme <- tryCatch(
      lme_analysis(td$trialpaths, dat, op),
      error = function(e) data.table(beta = NA_real_,
                                     betaSE = NA_real_,
                                     p = NA_real_,
                                     issingular = NA,
                                     warning = paste0('err: ',
                                                      conditionMessage(e)))
    )
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

    out[[r]] <- data.table(
      gating_slope = gating_slope, rep_idx = r,
      lme_p   = out_lme$p,
      lme_beta = out_lme$beta,
      lcmm_gating_p   = out_lcmm$p,
      lcmm_gating_beta = out_lcmm$beta,
      lcmm_within_p   = out_lcmm$p_within,
      lcmm_within_beta = out_lcmm$beta_within,
      delta_bic       = out_lcmm$delta_bic,
      lcmm_entropy    = out_lcmm$entropy,
      lcmm_conv       = out_lcmm$conv
    )
  }
  rbindlist(out)
}

cat(sprintf('Pilot 4: power on heterogeneous-random-slopes DGP\n'))
cat(sprintf('  N = %d, n_reps per cell = %d, %d cells\n',
            N_PT, N_REPS, length(GATING_SLOPES)))
cat(sprintf('  Started %s\n', format(Sys.time())))

t0 <- Sys.time()
results_list <- vector('list', length(GATING_SLOPES))
for (i in seq_along(GATING_SLOPES)) {
  gs <- GATING_SLOPES[i]
  cat(sprintf('  cell %d/%d: gating_slope = %.2f ...',
              i, length(GATING_SLOPES), gs))
  t_cell <- Sys.time()
  results_list[[i]] <- run_one_cell(gs, N_REPS)
  cat(sprintf(' done in %.0f s\n',
              as.numeric(difftime(Sys.time(), t_cell,
                                   units = 'secs'))))
}
results <- rbindlist(results_list)

## Summarise -------------------------------------------------------
alpha <- 0.05
mcse <- function(p, n) sqrt(p * (1 - p) / n)
summary_tbl <- results[, .(
  n_reps = .N,
  conv_rate = mean(lcmm_conv == 1, na.rm = TRUE),
  power_lme = mean(lme_p < alpha, na.rm = TRUE),
  power_lcmm_gating = mean(lcmm_gating_p < alpha, na.rm = TRUE),
  power_lcmm_within = mean(lcmm_within_p < alpha, na.rm = TRUE),
  bic_select_ng2 = mean(delta_bic > 0, na.rm = TRUE),
  bic_strong_ng2 = mean(delta_bic > 6, na.rm = TRUE),
  power_bic_cond_gating = mean(delta_bic > 0 &
                                 lcmm_gating_p < alpha,
                                 na.rm = TRUE),
  power_bic_cond_within = mean(delta_bic > 0 &
                                 lcmm_within_p < alpha,
                                 na.rm = TRUE),
  mean_entropy = mean(lcmm_entropy, na.rm = TRUE)
), by = gating_slope][order(gating_slope)]

cat('\n=== Pilot 4 summary ===\n')
print(summary_tbl, digits = 3)
cat(sprintf('\nMCSE for proportion p = 0.05: %.3f\n',
            mcse(0.05, N_REPS)))
cat(sprintf('MCSE for proportion p = 0.50: %.3f\n',
            mcse(0.50, N_REPS)))

out_path <- file.path(OUT_DIR, '04-pilot-power-random-slopes.rds')
saveRDS(list(
  config = list(N_REPS = N_REPS, N_PT = N_PT,
                SEED_BASE = SEED_BASE,
                GATING_SLOPES = GATING_SLOPES,
                started = format(t0)),
  results = results,
  summary = summary_tbl
), out_path)
cat(sprintf('\nOutput written to %s\n', out_path))
