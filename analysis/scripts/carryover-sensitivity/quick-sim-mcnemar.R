## quick-sim-mcnemar.R
##
## Pairwise McNemar tests on the 02-carryover-sensitivity per-rep
## replicate file, restricted to c.bm = 0.45. Reports A2 vs A1 and
## A2 vs A3 paired comparisons by (N, t1half_dgp) cell, plus the
## A3 gap-recovery fraction (A3-A1)/(A2-A1).

suppressPackageStartupMessages({
  library(data.table)
})

rep <- readRDS('analysis/data/quick-sim/02-carryover-replicates.rds')
rep_alt <- rep[c.bm == 0.45]
rep_alt[, sig := p < 0.05]

pair_mcnemar <- function(d, spec_a, spec_b) {
  da <- d[analysis_spec == spec_a, .(rep_idx, sig_a = sig)]
  db <- d[analysis_spec == spec_b, .(rep_idx, sig_b = sig)]
  m  <- merge(da, db, by = 'rep_idx')
  b01 <- sum(!m$sig_a &  m$sig_b, na.rm = TRUE)
  b10 <- sum( m$sig_a & !m$sig_b, na.rm = TRUE)
  n   <- nrow(m)
  diff_pp <- mean(m$sig_b, na.rm = TRUE) -
             mean(m$sig_a, na.rm = TRUE)
  if ((b01 + b10) > 0) {
    chi <- (abs(b01 - b10) - 1)^2 / (b01 + b10)
    pval <- pchisq(chi, df = 1, lower.tail = FALSE)
  } else {
    pval <- 1
  }
  list(n = n, b_a_only = b10, b_b_only = b01,
       diff_pp = diff_pp, mcnemar_p = pval)
}

cells <- unique(rep_alt[, .(N, t1half_dgp)])
setorder(cells, N, t1half_dgp)

cat('\n=== A2 vs A1 (paired, by cell, c.bm = 0.45) ===\n')
res21 <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  N0 <- cells$N[i]; th0 <- cells$t1half_dgp[i]
  d  <- rep_alt[N == N0 & t1half_dgp == th0]
  out <- pair_mcnemar(d, 'A1_binary', 'A2_matched')
  data.table(N = N0, t1half_dgp = th0,
             A1_only = out$b_a_only, A2_only = out$b_b_only,
             diff_A2_minus_A1 = out$diff_pp,
             mcnemar_p = out$mcnemar_p)
}))
print(res21)

cat('\n=== A2 vs A3 (paired, by cell, c.bm = 0.45) ===\n')
res23 <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  N0 <- cells$N[i]; th0 <- cells$t1half_dgp[i]
  d  <- rep_alt[N == N0 & t1half_dgp == th0]
  out <- pair_mcnemar(d, 'A3_lagged', 'A2_matched')
  data.table(N = N0, t1half_dgp = th0,
             A3_only = out$b_a_only, A2_only = out$b_b_only,
             diff_A2_minus_A3 = out$diff_pp,
             mcnemar_p = out$mcnemar_p)
}))
print(res23)

cat('\n=== Gap recovery: (A3-A1)/(A2-A1) by cell ===\n')
gap_dt <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  N0 <- cells$N[i]; th0 <- cells$t1half_dgp[i]
  d  <- rep_alt[N == N0 & t1half_dgp == th0]
  p1 <- mean(d[analysis_spec == 'A1_binary',  sig], na.rm = TRUE)
  p2 <- mean(d[analysis_spec == 'A2_matched', sig], na.rm = TRUE)
  p3 <- mean(d[analysis_spec == 'A3_lagged',  sig], na.rm = TRUE)
  data.table(N = N0, t1half_dgp = th0,
             pwr_A1 = p1, pwr_A2 = p2, pwr_A3 = p3,
             gap_A2_A1 = p2 - p1,
             gap_A3_A1 = p3 - p1,
             recovery = (p3 - p1) / (p2 - p1))
}))
print(gap_dt)
