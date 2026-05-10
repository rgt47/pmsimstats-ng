## quick-sim-prototype-30min.R
##
## 30-minute prototype simulation for paper 02-carryover-sensitivity.
## Expands the 1k driver from 6 cells (3 spec x 2 c.bm) to 24 cells
## (3 spec x 2 c.bm x 2 N x 2 DGP carryover half-life). Wall-clock
## budget 1700 s with graceful abort.

suppressPackageStartupMessages({
  library(devtools)
  library(data.table)
  library(ggplot2)
})

t_start <- Sys.time()
budget_secs <- 1700

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

## Analysis-side specifications. A2's analysis-side carryover
## half-life intentionally tracks the DGP value within each cell so
## that 'matched' remains correctly matched as t1half_dgp varies.
make_op_specs <- function(t1half_dgp) {
  list(
    A1_binary = list(useDE = TRUE, t_random_slope = FALSE,
                     full_model_out = FALSE, carryover_t1half = 0,
                     simplecarryover = FALSE,
                     carryover_scalefactor = 1),
    A2_matched = list(useDE = TRUE, t_random_slope = FALSE,
                      full_model_out = FALSE,
                      carryover_t1half = t1half_dgp,
                      simplecarryover = FALSE,
                      carryover_scalefactor = 1),
    A3_lagged = list(useDE = TRUE, t_random_slope = FALSE,
                     full_model_out = FALSE, carryover_t1half = 0,
                     simplecarryover = TRUE,
                     carryover_scalefactor = 1)
  )
}

c_bm_levels    <- c(0, 0.45)
N_levels       <- c(35L, 70L)
t1half_levels  <- c(0.5, 1.0)
n_reps_target  <- 800L

dgp_grid <- CJ(c_bm = c_bm_levels,
               N = N_levels,
               t1half_dgp = t1half_levels,
               sorted = FALSE)
n_dgp_cells <- nrow(dgp_grid)

cat('Starting simulation: target = ', n_reps_target,
    ' reps x ', n_dgp_cells,
    ' DGP cells x 3 specs = ',
    n_reps_target * n_dgp_cells * 3L, ' fits.\n', sep = '')
cat('Wall-clock budget: ', budget_secs, ' s.\n', sep = '')

simulate_one_dataset <- function(c_bm, N_total, t1half_dgp, seed) {
  set.seed(seed)
  N_per_path <- c(ceiling(N_total / 2), floor(N_total / 2))
  mp1 <- data.table(N = N_per_path[1], c.bm = c_bm,
                    carryover_t1half = t1half_dgp,
                    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
                    c.cf1t = 0.1, c.cfct = 0.05)
  mp2 <- data.table(N = N_per_path[2], c.bm = c_bm,
                    carryover_t1half = t1half_dgp,
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
fits_done <- 0L
aborted <- FALSE
reps_completed <- vector('list', n_dgp_cells)
names(reps_completed) <- paste(dgp_grid$c_bm, dgp_grid$N,
                               dgp_grid$t1half_dgp, sep = '_')

for (cell_i in seq_len(n_dgp_cells)) {
  cb         <- dgp_grid$c_bm[cell_i]
  Nv         <- dgp_grid$N[cell_i]
  th         <- dgp_grid$t1half_dgp[cell_i]
  cell_key   <- names(reps_completed)[cell_i]
  op_specs   <- make_op_specs(th)
  reps_completed[[cell_key]] <- 0L

  cat(sprintf('Cell %d/%d: c.bm = %.2f, N = %d, t1half_dgp = %.2f\n',
              cell_i, n_dgp_cells, cb, Nv, th))

  for (rep_idx in seq_len(n_reps_target)) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
    if (elapsed > budget_secs) {
      cat('Wall-clock budget exceeded at cell ', cell_i,
          ' rep ', rep_idx, '; aborting cleanly.\n', sep = '')
      aborted <- TRUE
      break
    }
    seed_i <- 1000000L * cell_i + rep_idx
    dat <- tryCatch(
      simulate_one_dataset(cb, Nv, th, seed_i),
      error = function(e) NULL
    )
    if (is.null(dat)) next
    for (spec_name in names(op_specs)) {
      fit <- tryCatch(
        lme_analysis(td$trialpaths, dat, op_specs[[spec_name]]),
        error = function(e) NULL
      )
      row_i <- row_i + 1L
      fits_done <- fits_done + 1L
      if (is.null(fit)) {
        results[[row_i]] <- data.table(
          analysis_spec = spec_name, c.bm = cb,
          N = Nv, t1half_dgp = th,
          rep_idx = rep_idx,
          beta = NA_real_, betaSE = NA_real_, p = NA_real_,
          converged = FALSE
        )
      } else {
        results[[row_i]] <- data.table(
          analysis_spec = spec_name, c.bm = cb,
          N = Nv, t1half_dgp = th,
          rep_idx = rep_idx,
          beta = fit$beta, betaSE = fit$betaSE, p = fit$p,
          converged = !is.na(fit$beta)
        )
      }
      if (fits_done %% 500L == 0L) {
        cat('  fits = ', fits_done, '; cell ', cell_i,
            ' rep ', rep_idx, '; elapsed = ', round(elapsed, 1),
            ' s.\n', sep = '')
      }
    }
    reps_completed[[cell_key]] <- rep_idx
  }
  if (aborted) break
}

replicates <- rbindlist(results)

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat('Simulation finished. Elapsed: ', round(elapsed_total, 1),
    ' s. Aborted: ', aborted, '. Total fits: ', fits_done, '.\n',
    sep = '')
for (k in names(reps_completed)) {
  cat('  cell ', k, ': ', reps_completed[[k]], ' reps completed.\n',
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
}, by = .(analysis_spec, c.bm, N, t1half_dgp)]

setorder(summary_dt, c.bm, N, t1half_dgp, analysis_spec)

write.table(summary_dt,
            file = 'analysis/data/quick-sim/02-carryover-summary.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

cat('\nMorris summary table:\n')
print(summary_dt)

plot_dt <- copy(summary_dt[c.bm == 0.45])
plot_dt[, lo := pmax(0, power - 1.96 * mcse_power)]
plot_dt[, hi := pmin(1, power + 1.96 * mcse_power)]
plot_dt[, N_lbl := paste0('N = ', N)]
plot_dt[, th_lbl := paste0('DGP t1/2 = ', t1half_dgp, ' wk')]

p <- ggplot(plot_dt,
            aes(x = analysis_spec, y = power,
                ymin = lo, ymax = hi)) +
  geom_col(fill = 'grey75', colour = 'grey30') +
  geom_errorbar(width = 0.2) +
  geom_text(aes(label = sprintf('%.3f', power)),
            vjust = -0.6, size = 3.0) +
  facet_grid(N_lbl ~ th_lbl) +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(x = 'Analysis-side carryover specification',
       y = 'Power',
       title = paste0(
         'Power by analysis-side carryover specification ',
         '(c.bm = 0.45)'),
       subtitle = paste0(
         'OL+BDC; faceted by N x DGP carryover half-life; ',
         'reps per cell = ', max(summary_dt$n_reps),
         ' (target ', n_reps_target, ')')) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave('analysis/figures/quick-sim/02-carryover-power.pdf',
       p, width = 8, height = 6)

cat('\nFigure written to analysis/figures/quick-sim/',
    '02-carryover-power.pdf\n', sep = '')

invisible(NULL)
