## quick-sim-prototype-1k.R
##
## Higher-rep prototype simulation for paper 02-carryover-sensitivity.
## Identical cell structure to quick-sim-prototype.R; only n_reps_target
## and the wall-clock budget are bumped (200 -> 1000, 4 min -> 14 min).

suppressPackageStartupMessages({
  library(devtools)
  library(data.table)
  library(ggplot2)
})

t_start <- Sys.time()
budget_secs <- 14 * 60

load_all('.', quiet = TRUE)

set.seed(20260507)

td <- buildtrialdesign(
  name_longform = 'OL+BDC',
  name_shortform = 'OLBDC',
  timepoints   = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames  = c('OL1', 'OL2', 'OL3', 'OL4',
                   'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

N_total <- 35
N_per_path <- c(ceiling(N_total / 2), floor(N_total / 2))

resp_param <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)

baseline_param <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)

op_specs <- list(
  A1_binary = list(useDE = TRUE, t_random_slope = FALSE,
                   full_model_out = FALSE, carryover_t1half = 0,
                   simplecarryover = FALSE, carryover_scalefactor = 1),
  A2_matched = list(useDE = TRUE, t_random_slope = FALSE,
                    full_model_out = FALSE, carryover_t1half = 1.0,
                    simplecarryover = FALSE, carryover_scalefactor = 1),
  A3_lagged = list(useDE = TRUE, t_random_slope = FALSE,
                   full_model_out = FALSE, carryover_t1half = 0,
                   simplecarryover = TRUE, carryover_scalefactor = 1)
)

c_bm_levels <- c(0, 0.45)
n_reps_target <- 1000L

simulate_one_dataset <- function(c_bm, seed) {
  set.seed(seed)
  mp1 <- data.table(N = N_per_path[1], c.bm = c_bm,
                    carryover_t1half = 1.0,
                    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
                    c.cf1t = 0.1, c.cfct = 0.05)
  mp2 <- data.table(N = N_per_path[2], c.bm = c_bm,
                    carryover_t1half = 1.0,
                    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
                    c.cf1t = 0.1, c.cfct = 0.05)
  d1 <- generateData(mp1, resp_param, baseline_param,
                     td$trialpaths[[1]],
                     empirical = FALSE,
                     makePositiveDefinite = TRUE)
  d1[, path := 1]
  d2 <- generateData(mp2, resp_param, baseline_param,
                     td$trialpaths[[2]],
                     empirical = FALSE,
                     makePositiveDefinite = TRUE)
  d2[, path := 2]
  d2[, ptID := ptID + max(d1$ptID)]
  rbind(d1, d2)
}

results <- list()
row_i <- 0L
aborted <- FALSE
reps_completed <- list()

cat('Starting simulation: target = ', n_reps_target,
    ' reps x ', length(c_bm_levels), ' c.bm levels x ',
    length(op_specs), ' specs.\n', sep = '')
cat('Wall-clock budget: ', budget_secs, ' s.\n', sep = '')

for (cb in c_bm_levels) {
  reps_completed[[as.character(cb)]] <- 0L
  for (rep_idx in seq_len(n_reps_target)) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
    if (elapsed > budget_secs) {
      cat('Wall-clock budget exceeded at rep ', rep_idx,
          ' c.bm = ', cb, '; aborting cleanly.\n', sep = '')
      aborted <- TRUE
      break
    }
    seed_i <- 1000000L * which(c_bm_levels == cb) + rep_idx
    dat <- tryCatch(
      simulate_one_dataset(cb, seed_i),
      error = function(e) NULL
    )
    if (is.null(dat)) next
    for (spec_name in names(op_specs)) {
      fit <- tryCatch(
        lme_analysis(td$trialpaths, dat, op_specs[[spec_name]]),
        error = function(e) NULL
      )
      row_i <- row_i + 1L
      if (is.null(fit)) {
        results[[row_i]] <- data.table(
          analysis_spec = spec_name, c.bm = cb,
          rep_idx = rep_idx,
          beta = NA_real_, betaSE = NA_real_, p = NA_real_,
          converged = FALSE
        )
      } else {
        results[[row_i]] <- data.table(
          analysis_spec = spec_name, c.bm = cb,
          rep_idx = rep_idx,
          beta = fit$beta, betaSE = fit$betaSE, p = fit$p,
          converged = !is.na(fit$beta)
        )
      }
    }
    reps_completed[[as.character(cb)]] <- rep_idx
    if (rep_idx %% 100L == 0L) {
      cat('  c.bm = ', cb, ': ', rep_idx, ' / ', n_reps_target,
          ' reps; elapsed = ', round(elapsed, 1), ' s.\n', sep = '')
    }
  }
  if (aborted) break
}

replicates <- rbindlist(results)

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat('Simulation finished. Elapsed: ', round(elapsed_total, 1),
    ' s. Aborted: ', aborted, '.\n', sep = '')
for (cb in names(reps_completed)) {
  cat('  c.bm = ', cb, ': ', reps_completed[[cb]], ' reps completed.\n',
      sep = '')
}

dir.create('analysis/data/quick-sim', showWarnings = FALSE,
           recursive = TRUE)
dir.create('analysis/figures/quick-sim', showWarnings = FALSE,
           recursive = TRUE)

saveRDS(replicates,
        'analysis/data/quick-sim/02-carryover-replicates.rds')

summary_dt <- replicates[, {
  conv <- !is.na(beta)
  n_total <- .N
  n_conv  <- sum(conv)
  pwr     <- mean(p[conv] < 0.05, na.rm = TRUE)
  mcse_pwr <- sqrt(pwr * (1 - pwr) / max(n_conv, 1))
  m_b   <- mean(beta[conv], na.rm = TRUE)
  s_b   <- sd(beta[conv], na.rm = TRUE)
  mcse_b <- s_b / sqrt(max(n_conv, 1))
  list(n_reps = n_total,
       power = pwr, mcse_power = mcse_pwr,
       mean_beta = m_b, sd_beta = s_b,
       mcse_mean_beta = mcse_b,
       conv_rate = n_conv / n_total)
}, by = .(analysis_spec, c.bm)]

setorder(summary_dt, c.bm, analysis_spec)

write.table(summary_dt,
            file = 'analysis/data/quick-sim/02-carryover-summary.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

cat('\nMorris summary table:\n')
print(summary_dt)

plot_dt <- copy(summary_dt)
plot_dt[, c_bm_lbl := paste0('c.bm = ', c.bm)]
plot_dt[, lo := pmax(0, power - 1.96 * mcse_power)]
plot_dt[, hi := pmin(1, power + 1.96 * mcse_power)]

p <- ggplot(plot_dt,
            aes(x = analysis_spec, y = power,
                ymin = lo, ymax = hi)) +
  geom_col(fill = 'grey75', colour = 'grey30') +
  geom_errorbar(width = 0.2) +
  geom_text(aes(label = sprintf('%.3f', power)),
            vjust = -0.6, size = 3.2) +
  facet_wrap(~ c_bm_lbl) +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(x = 'Analysis-side carryover specification',
       y = 'Power (or type I error at c.bm = 0)',
       title = 'Power by analysis-side carryover specification',
       subtitle = paste0(
         'OL+BDC, N = 35, DGP t1/2 = 1.0 wk, ',
         'reps = ', max(summary_dt$n_reps),
         ' (target ', n_reps_target, ')')) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave('analysis/figures/quick-sim/02-carryover-power.pdf',
       p, width = 7, height = 4)

cat('\nFigure written to analysis/figures/quick-sim/',
    '02-carryover-power.pdf\n', sep = '')

invisible(NULL)
