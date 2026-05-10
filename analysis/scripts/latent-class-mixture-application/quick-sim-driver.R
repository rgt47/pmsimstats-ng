## Quick-sim prototype driver for paper 03-latent-class-mixture-
## application. Three gating-slope cells x lme/lcmm analyses.
## Per-replicate cost is dominated by lcmm; budget is 4 minutes
## wall-clock with 30 reps per cell, halving to 15 if first 15
## reps in cell 1 take >50 s.
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/quick-sim-driver.R

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
})

source('analysis/scripts/latent-class-mixture-application/01-lcmm-wrapper.R')
source('analysis/scripts/latent-class-mixture-application/02-heterogeneous-slopes-dgp.R')

DATA_DIR <- file.path('analysis', 'data', 'quick-sim')
FIG_DIR  <- file.path('analysis', 'figures', 'quick-sim')
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

SEED_BASE <- 20260507L
GATING_SLOPES <- c(0, 1.0, 1.5)
N_REPS_INITIAL <- 30L
N_REPS_FALLBACK <- 15L
PROBE_BUDGET_SEC <- 50

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

make_modelparam <- function(gating_slope) {
  if (gating_slope == 0) {
    data.table(
      N = 70, c.bm = 0,
      carryover_t1half = 1.0,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.1, c.cfct = 0.05
    )
  } else {
    data.table(
      N = 70, c.bm = 0,
      carryover_t1half = 1.0,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.1, c.cfct = 0.05,
      beta_R_factor = 1.0, beta_NR_factor = 0.0,
      gating_intercept = 0.0, gating_slope = gating_slope
    )
  }
}

simulate_dat <- function(gating_slope, rep_idx) {
  set.seed(SEED_BASE + rep_idx + 1000L * which(GATING_SLOPES == gating_slope))
  mp <- make_modelparam(gating_slope)
  if (gating_slope == 0) {
    dat <- generateData(
      modelparam = mp, respparam = rp, blparam = bp,
      trialdesign = td$trialpaths[[1]],
      empirical = FALSE, makePositiveDefinite = TRUE,
      dgp_architecture = 'mvn'
    )
  } else {
    dat <- generate_random_slopes_data(
      modelparam = mp, respparam = rp, blparam = bp,
      trialdesign = td$trialpaths[[1]],
      empirical = FALSE, makePositiveDefinite = TRUE
    )
  }
  dat[, path := 1]
  dat
}

run_one <- function(gating_slope, rep_idx) {
  dat <- simulate_dat(gating_slope, rep_idx)

  t_lme0 <- Sys.time()
  out_lme <- tryCatch(
    lme_analysis(td$trialpaths, dat, op),
    error = function(e) data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      issingular = NA, warning = paste0('err: ', conditionMessage(e))
    )
  )
  t_lme <- as.numeric(difftime(Sys.time(), t_lme0, units = 'secs'))

  t_lcmm0 <- Sys.time()
  out_lcmm <- tryCatch(
    lcmm_analysis(td$trialpaths, dat, op),
    error = function(e) data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      beta_within = NA_real_, betaSE_within = NA_real_,
      p_within = NA_real_, delta_bic = NA_real_,
      lrt_stat = NA_real_, lrt_df = NA_integer_,
      bic_ng1 = NA_real_, bic_ng2 = NA_real_,
      issingular = NA, warning = paste0('err: ', conditionMessage(e)),
      entropy = NA_real_, bic = NA_real_, conv = 0L, conv_ng1 = 0L
    )
  )
  t_lcmm <- as.numeric(difftime(Sys.time(), t_lcmm0, units = 'secs'))

  list(
    lme = data.table(
      gating_slope = gating_slope, analysis = 'lme', rep_idx = rep_idx,
      beta = out_lme$beta, betaSE = out_lme$betaSE, p = out_lme$p,
      gating_p = NA_real_, within_p = NA_real_, delta_bic = NA_real_,
      conv = as.integer(!is.na(out_lme$p)), elapsed_s = t_lme
    ),
    lcmm = data.table(
      gating_slope = gating_slope, analysis = 'lcmm', rep_idx = rep_idx,
      beta = out_lcmm$beta, betaSE = out_lcmm$betaSE, p = out_lcmm$p,
      gating_p = out_lcmm$p, within_p = out_lcmm$p_within,
      delta_bic = out_lcmm$delta_bic,
      conv = as.integer(out_lcmm$conv == 1), elapsed_s = t_lcmm
    )
  )
}

t_global <- Sys.time()
cat(sprintf('Quick-sim driver started at %s\n', format(t_global)))
cat(sprintf('Cells: gating_slope = %s; initial reps per cell = %d\n',
            paste(GATING_SLOPES, collapse = ', '), N_REPS_INITIAL))

## Probe phase: run first 15 reps in cell 1 (gating=0).
n_per_cell <- N_REPS_INITIAL
probe_n <- 15L
cat(sprintf('Probe: first %d reps in cell gating=%g...\n',
            probe_n, GATING_SLOPES[1]))

probe_results <- vector('list', probe_n)
t_probe0 <- Sys.time()
for (r in seq_len(probe_n)) {
  probe_results[[r]] <- run_one(GATING_SLOPES[1], r)
}
t_probe <- as.numeric(difftime(Sys.time(), t_probe0, units = 'secs'))
cat(sprintf('  probe elapsed: %.1f s for %d reps\n', t_probe, probe_n))

if (t_probe > PROBE_BUDGET_SEC) {
  cat(sprintf('  probe exceeded %d s budget; halving to %d reps/cell\n',
              PROBE_BUDGET_SEC, N_REPS_FALLBACK))
  n_per_cell <- N_REPS_FALLBACK
}

## Continue cell 1 to n_per_cell, then run remaining cells.
all_results <- list()
all_results <- c(all_results, probe_results)

if (n_per_cell > probe_n) {
  for (r in (probe_n + 1L):n_per_cell) {
    all_results[[length(all_results) + 1L]] <- run_one(GATING_SLOPES[1], r)
  }
}

for (g in GATING_SLOPES[-1]) {
  cat(sprintf('Cell gating_slope = %g: running %d reps...\n', g, n_per_cell))
  t_cell0 <- Sys.time()
  for (r in seq_len(n_per_cell)) {
    all_results[[length(all_results) + 1L]] <- run_one(g, r)
  }
  t_cell <- as.numeric(difftime(Sys.time(), t_cell0, units = 'secs'))
  cat(sprintf('  cell elapsed: %.1f s\n', t_cell))
}

t_total <- as.numeric(difftime(Sys.time(), t_global, units = 'secs'))
cat(sprintf('\nTotal wall time: %.1f s (%.2f min)\n', t_total, t_total / 60))

## Stitch per-replicate data.
per_rep <- rbindlist(c(
  lapply(all_results, `[[`, 'lme'),
  lapply(all_results, `[[`, 'lcmm')
), use.names = TRUE, fill = TRUE)
setorder(per_rep, gating_slope, analysis, rep_idx)

saveRDS(per_rep, file.path(DATA_DIR, '03-latent-replicates.rds'))
cat(sprintf('Per-replicate output -> %s (rows=%d)\n',
            file.path(DATA_DIR, '03-latent-replicates.rds'),
            nrow(per_rep)))

## Morris summary table -----------------------------------------
##
## Three test channels:
##   lme_bmDbc                       (analysis = 'lme', use p)
##   lcmm_gating_unconditional       (analysis = 'lcmm', use gating_p)
##   lcmm_gating_BIC_conditional     ('lcmm', delta_bic > 0 AND gating_p<.05)

binom_mcse <- function(rate, n) sqrt(rate * (1 - rate) / n)

summarise_cell <- function(df_lme, df_lcmm, gating_slope) {
  dgp_label <- if (gating_slope == 0) 'null' else 'random-slopes'

  ## lme bm:Dbc
  n_lme <- sum(!is.na(df_lme$p))
  conv_lme <- mean(df_lme$conv == 1, na.rm = TRUE)
  rej_lme <- mean(df_lme$p < 0.05, na.rm = TRUE)

  ## lcmm gating unconditional
  n_lcmm <- sum(!is.na(df_lcmm$gating_p))
  conv_lcmm <- mean(df_lcmm$conv == 1, na.rm = TRUE)
  rej_uncond <- mean(df_lcmm$gating_p < 0.05, na.rm = TRUE)

  ## lcmm BIC-conditional gating: rejection counted only when
  ## delta_bic > 0 AND gating_p < 0.05, with denominator = total
  ## convergent lcmm reps (Morris ADEMP convention: denominator
  ## fixed across tests within cell so rejection rates compose).
  rej_cond <- mean(
    (df_lcmm$delta_bic > 0) & (df_lcmm$gating_p < 0.05),
    na.rm = TRUE
  )

  rbind(
    data.table(
      dgp = dgp_label, gating_slope = gating_slope,
      analysis_test = 'lme_bmDbc',
      n_reps = n_lme, rejection_rate = rej_lme,
      mcse_rate = binom_mcse(rej_lme, n_lme),
      conv_rate = conv_lme
    ),
    data.table(
      dgp = dgp_label, gating_slope = gating_slope,
      analysis_test = 'lcmm_gating_unconditional',
      n_reps = n_lcmm, rejection_rate = rej_uncond,
      mcse_rate = binom_mcse(rej_uncond, n_lcmm),
      conv_rate = conv_lcmm
    ),
    data.table(
      dgp = dgp_label, gating_slope = gating_slope,
      analysis_test = 'lcmm_gating_BIC_conditional',
      n_reps = n_lcmm, rejection_rate = rej_cond,
      mcse_rate = binom_mcse(rej_cond, n_lcmm),
      conv_rate = conv_lcmm
    )
  )
}

summary_tab <- rbindlist(lapply(GATING_SLOPES, function(g) {
  df_lme  <- per_rep[gating_slope == g & analysis == 'lme']
  df_lcmm <- per_rep[gating_slope == g & analysis == 'lcmm']
  summarise_cell(df_lme, df_lcmm, g)
}))

summary_path <- file.path(DATA_DIR, '03-latent-summary.txt')
fwrite(summary_tab, summary_path, sep = '\t')
cat(sprintf('Morris summary -> %s\n', summary_path))
print(summary_tab)

## Figure --------------------------------------------------------
plot_dat <- copy(summary_tab)
plot_dat[, ci_lo := pmax(0, rejection_rate - 1.96 * mcse_rate)]
plot_dat[, ci_hi := pmin(1, rejection_rate + 1.96 * mcse_rate)]
plot_dat[, analysis_test := factor(
  analysis_test,
  levels = c('lme_bmDbc',
             'lcmm_gating_unconditional',
             'lcmm_gating_BIC_conditional'),
  labels = c('lme bm:Dbc',
             'lcmm gating Wald (unconditional)',
             'lcmm gating, BIC-conditional')
)]

p_fig <- ggplot(plot_dat,
                aes(x = gating_slope, y = rejection_rate)) +
  geom_hline(yintercept = 0.05, linetype = 'dashed',
             colour = 'grey50') +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.08) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~ analysis_test, ncol = 3) +
  scale_x_continuous(breaks = GATING_SLOPES) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = 'Gating slope (heterogeneous-slopes DGP intensity)',
    y = 'Rejection rate (alpha = 0.05)',
    title = sprintf(
      'Quick-sim: class-aware vs lme moderation tests (n_reps = %d/cell)',
      n_per_cell
    ),
    caption = sprintf(
      'Reference line at 0.05; error bars are 95%% binomial CIs. Total wall: %.1f min.',
      t_total / 60
    )
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 11),
        plot.caption = element_text(size = 8))

fig_path <- file.path(FIG_DIR, '03-latent-power.pdf')
ggsave(fig_path, p_fig, width = 9, height = 3.5)
cat(sprintf('Figure -> %s\n', fig_path))

cat('\nDone.\n')
