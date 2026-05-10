# 05-design-prototype.R
#
# Expanded prototype simulation for paper
# 05-nof1-design-sensitivity.
#
# Sweep S8 expanded to a 5 x 3 factorial grid:
#   c.br levels x N levels
#     c.br in {0, 0.15, 0.30, 0.45, 0.60}
#     N    in {25, 35, 50}
# Target 500 reps per cell (MCSE ~ 0.022 at power = 0.5).
# Estimand is the Dbc treatment main-effect coefficient.
#
# Trial design: hybrid OL+BDC (Hendrickson 2020 design 4).
#
# Time budget: 1700 s (30 min target with 100 s headroom).
# Parallelisation: parallel::mclapply with detectCores() - 1
# workers. Progress message every 200 completed fits.
#
# Outputs (overwritten):
#   analysis/data/quick-sim/05-design-replicates.rds
#   analysis/data/quick-sim/05-design-summary.txt
#   analysis/figures/quick-sim/05-design-power.pdf

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

t_start <- Sys.time()
budget_secs <- 1700

set.seed(20260507)

n_cores <- max(1L, parallel::detectCores() - 1L)
cat('Cores:', n_cores, '\n')

# ---------------------------------------------------------------
# Trial design: hybrid OL+BDC (Hendrickson 2020 design 4)
# ---------------------------------------------------------------

td <- buildtrialdesign(
  name_longform = 'OL + BDC + crossover (hybrid)',
  name_shortform = 'Hybrid',
  timepoints = seq(1, 8, 1),
  timeptnames = paste0('W', 1:8),
  expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 1),
    pathB = c(1, 1, 1, 1, 1, 1, 0, 0)
  )
)

# ---------------------------------------------------------------
# Parameter sets (matched to test_lme_analysis.R conventions)
# ---------------------------------------------------------------

bp <- data.table(
  cat = c('bm', 'BL'),
  m = c(0, 70),
  sd = c(1, 10)
)

rp <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, -2.0),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)

# ---------------------------------------------------------------
# Sweep grid: c.br x N
# ---------------------------------------------------------------

c_br_levels <- c(0, 0.15, 0.30, 0.45, 0.60)
N_levels    <- c(25L, 35L, 50L)
n_reps      <- 500L

cell_grid <- CJ(c_br = c_br_levels, N = N_levels, sorted = FALSE)
n_cells <- nrow(cell_grid)
cat('Grid:', n_cells, 'cells x', n_reps, 'reps =',
    n_cells * n_reps, 'fits\n')

# ---------------------------------------------------------------
# Single-rep worker
# ---------------------------------------------------------------

run_one_rep <- function(c_br_val, N_val, seed) {
  set.seed(seed)

  mp <- data.table(
    N = as.integer(N_val),
    c.bm = 0,
    carryover_t1half = 1.0,
    c.tv = 0.7, c.pb = 0.7, c.br = c_br_val,
    c.cf1t = 0.1, c.cfct = 0.05
  )

  out_na <- function(reason) {
    data.table(
      beta_Dbc = NA_real_, betaSE_Dbc = NA_real_,
      p_Dbc = NA_real_, converged = FALSE,
      reason = reason
    )
  }

  tryCatch(
    {
      dat1 <- generateData(
        modelparam = mp, respparam = rp, blparam = bp,
        trialdesign = td$trialpaths[[1]],
        empirical = FALSE, makePositiveDefinite = TRUE
      )
      dat1[, path := 1L]
      dat2 <- generateData(
        modelparam = mp, respparam = rp, blparam = bp,
        trialdesign = td$trialpaths[[2]],
        empirical = FALSE, makePositiveDefinite = TRUE
      )
      dat2[, path := 2L]
      dat <- rbind(dat1, dat2, fill = TRUE)

      op <- list(
        useDE = TRUE,
        t_random_slope = FALSE,
        full_model_out = TRUE,
        carryover_t1half = 1.0,
        simplecarryover = FALSE,
        carryover_scalefactor = 1
      )

      fit_obj <- lme_analysis(td$trialpaths, dat, op)
      if (is.null(fit_obj$fit)) {
        out_na('fit failed')
      } else {
        tt <- summary(fit_obj$fit)$tTable
        if ('Dbc' %in% rownames(tt)) {
          data.table(
            beta_Dbc   = tt['Dbc', 'Value'],
            betaSE_Dbc = tt['Dbc', 'Std.Error'],
            p_Dbc      = tt['Dbc', 'p-value'],
            converged  = TRUE,
            reason     = NA_character_
          )
        } else {
          out_na('no Dbc row')
        }
      }
    },
    error = function(e) out_na(conditionMessage(e))
  )
}

# ---------------------------------------------------------------
# Main loop: per-cell mclapply, abort gracefully on time budget
# ---------------------------------------------------------------

total_fits <- 0L
all_rows <- list()
aborted <- FALSE
abort_cell <- NA_integer_

for (i in seq_len(n_cells)) {
  c_br_val <- cell_grid$c_br[i]
  N_val    <- cell_grid$N[i]

  remaining <- budget_secs -
    as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
  if (remaining <= 0) {
    aborted <- TRUE
    abort_cell <- i
    cat(sprintf('\nBudget exhausted before cell %d.\n', i))
    break
  }

  cat(sprintf('\n[%4.0fs left] Cell %d/%d: c.br=%.2f N=%d\n',
              remaining, i, n_cells, c_br_val, N_val))

  cell_seeds <- 100000L * i + seq_len(n_reps)

  # Run reps in chunks so we can check the time budget mid-cell
  chunk_size <- min(n_reps, n_cores * 25L)
  cell_rows <- vector('list', n_reps)
  k_done <- 0L

  while (k_done < n_reps) {
    elapsed_now <- as.numeric(difftime(Sys.time(), t_start,
                                       units = 'secs'))
    if (elapsed_now > budget_secs) {
      cat('  budget exceeded mid-cell at',
          k_done, '/', n_reps, 'reps\n')
      aborted <- TRUE
      abort_cell <- i
      break
    }

    chunk_end <- min(k_done + chunk_size, n_reps)
    idx <- (k_done + 1L):chunk_end
    seeds_chunk <- cell_seeds[idx]

    chunk_res <- mclapply(
      seeds_chunk,
      function(s) run_one_rep(c_br_val, N_val, s),
      mc.cores = n_cores,
      mc.preschedule = TRUE
    )

    cell_rows[idx] <- chunk_res
    k_done <- chunk_end
    total_fits <- total_fits + length(idx)

    if (total_fits %% 200L < length(idx)) {
      el <- as.numeric(difftime(Sys.time(), t_start,
                                units = 'secs'))
      cat(sprintf('  ...%d fits done (cell %d, %d/%d), %.0fs\n',
                  total_fits, i, k_done, n_reps, el))
    }
  }

  cell_rows <- cell_rows[!vapply(cell_rows, is.null, logical(1))]
  if (length(cell_rows) == 0) next

  cell_dt <- rbindlist(cell_rows, fill = TRUE)
  cell_dt[, `:=`(c_br = c_br_val,
                 N = N_val,
                 rep_idx = seq_len(.N))]
  all_rows[[i]] <- cell_dt

  if (aborted) break
}

reps_dt <- rbindlist(all_rows, fill = TRUE)
setcolorder(reps_dt,
            c('c_br', 'N', 'rep_idx',
              'beta_Dbc', 'betaSE_Dbc', 'p_Dbc',
              'converged', 'reason'))

t_total <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat(sprintf('\nTotal wall-clock: %.1f s (budget %d s)\n',
            t_total, budget_secs))
cat(sprintf('Total reps stored: %d (aborted=%s)\n',
            nrow(reps_dt), aborted))

# ---------------------------------------------------------------
# Persist replicates
# ---------------------------------------------------------------

dir.create('analysis/data/quick-sim', showWarnings = FALSE,
           recursive = TRUE)
dir.create('analysis/figures/quick-sim', showWarnings = FALSE,
           recursive = TRUE)

reps_out <- reps_dt[, .(c_br, N, rep_idx,
                        beta_Dbc, betaSE_Dbc, p_Dbc,
                        converged)]
saveRDS(reps_out,
        'analysis/data/quick-sim/05-design-replicates.rds')

# ---------------------------------------------------------------
# Morris-style summary table (ADEMP-aligned)
# ---------------------------------------------------------------

alpha <- 0.05
summary_dt <- reps_dt[, {
  conv <- converged & !is.na(p_Dbc)
  n_used <- sum(conv)
  rej <- if (n_used > 0) sum(p_Dbc[conv] < alpha) else 0L
  rate <- if (n_used > 0) rej / n_used else NA_real_
  mcse_rate <- if (n_used > 0)
    sqrt(rate * (1 - rate) / n_used) else NA_real_
  mb <- if (n_used > 0) mean(beta_Dbc[conv]) else NA_real_
  sb <- if (n_used > 1) sd(beta_Dbc[conv]) else NA_real_
  mcse_mb <- if (n_used > 1) sb / sqrt(n_used) else NA_real_
  list(
    n_reps         = .N,
    n_converged    = n_used,
    conv_rate      = n_used / .N,
    power          = rate,
    mcse_power     = mcse_rate,
    mean_beta      = mb,
    sd_beta        = sb,
    mcse_mean_beta = mcse_mb
  )
}, by = .(c_br, N)]

setorder(summary_dt, N, c_br)

fwrite(summary_dt,
       'analysis/data/quick-sim/05-design-summary.txt',
       sep = '\t')

cat('\n--- Summary ---\n')
print(summary_dt)

# ---------------------------------------------------------------
# Power figure: faceted by N
# ---------------------------------------------------------------

plot_dt <- copy(summary_dt)
plot_dt[, ci_lo := pmax(0, power - 1.96 * mcse_power)]
plot_dt[, ci_hi := pmin(1, power + 1.96 * mcse_power)]
plot_dt[, N_lab := factor(paste0('N = ', N),
                          levels = paste0('N = ', sort(N_levels)))]

p_power <- ggplot(plot_dt,
                  aes(x = c_br, y = power)) +
  geom_hline(yintercept = 0.80, linetype = 'dotted',
             colour = 'grey50') +
  geom_hline(yintercept = 0.05, linetype = 'dashed',
             colour = 'grey60') +
  geom_line() +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.025) +
  facet_wrap(~ N_lab) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = c_br_levels) +
  labs(x = 'AR(1) within-factor BR correlation (c.br)',
       y = 'Power: P(p_Dbc < 0.05)',
       title = paste0('Treatment main-effect power vs AR(1) ',
                      'BR correlation, by N'),
       subtitle = sprintf(
         'Hybrid OL+BDC, br_max=-2, target %d reps/cell',
         n_reps),
       caption = paste0('Error bars: 1.96 x MCSE (binomial). ',
                        'Dotted = 0.80 power; dashed = nominal ',
                        'alpha. Estimand: Dbc main-effect.')) +
  theme_bw(base_size = 11)

ggsave('analysis/figures/quick-sim/05-design-power.pdf',
       p_power, width = 9, height = 4.5)

cat('\nWrote:\n')
cat('  analysis/data/quick-sim/05-design-replicates.rds\n')
cat('  analysis/data/quick-sim/05-design-summary.txt\n')
cat('  analysis/figures/quick-sim/05-design-power.pdf\n')

if (aborted) {
  cat(sprintf('\nNOTE: aborted at cell %d due to time budget. ',
              abort_cell))
  cat('Some cells have fewer than target reps.\n')
}
