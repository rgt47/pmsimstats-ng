## quick-sim-rewrite-summary.R
##
## Rewrite analysis/data/quick-sim/02-carryover-summary.txt with a
## proper column-header row (the original driver omitted col.names
## when appending to the metadata header). Re-derives the Morris
## summary from the per-replicate file. Run from the package root.

suppressPackageStartupMessages({
  library(data.table)
})

binom_mcse <- function(rate, n) {
  if (is.na(rate) || is.na(n) || n == 0) return(NA_real_)
  sqrt(rate * (1 - rate) / n)
}

rep <- readRDS('analysis/data/quick-sim/02-carryover-replicates.rds')

summary_dt <- rep[, {
  conv <- !is.na(beta)
  n_total <- .N
  n_conv  <- sum(conv)
  pwr     <- if (n_conv > 0) mean(p[conv] < 0.05, na.rm = TRUE)
             else NA_real_
  mcse_pwr <- binom_mcse(pwr, n_conv)
  m_b     <- if (n_conv > 0) mean(beta[conv], na.rm = TRUE)
             else NA_real_
  s_b     <- if (n_conv > 1) sd(beta[conv], na.rm = TRUE)
             else NA_real_
  mcse_b  <- if (n_conv > 0) s_b / sqrt(n_conv) else NA_real_
  list(n_reps = n_total,
       power = pwr, mcse_power = mcse_pwr,
       mean_beta = m_b, sd_beta = s_b,
       mcse_mean_beta = mcse_b,
       conv_rate = if (n_total > 0) n_conv / n_total
                   else NA_real_)
}, by = .(analysis_spec, c.bm, N, t1half_dgp)]

setorder(summary_dt, c.bm, N, t1half_dgp, analysis_spec)

summary_path <- 'analysis/data/quick-sim/02-carryover-summary.txt'
header_lines <- c(
  '# Morris ADEMP summary: 02-carryover-sensitivity',
  '# run_started: 2026-05-08 05:51:56',
  '# total_wall_sec: 435.9',
  '# mode_used: mclapply mc.cores=8',
  '# datasets_target: 4800',
  '# datasets_done: 4800',
  '# fits_done: 14400',
  '# time_was_binding: FALSE',
  '#'
)
writeLines(header_lines, summary_path)
fwrite(summary_dt, summary_path, sep = '\t', append = TRUE,
       col.names = TRUE)

cat('Rewrote ', summary_path, ' with column header row.\n', sep = '')
print(summary_dt)
