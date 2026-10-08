## Numeric checks for the propositions of docs/41, Section 4, using the
## implied outcome covariance of the published construct (Hybrid, c_bm = 0).
setwd('/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng')
suppressMessages(library(data.table))
dgp <- new.env()
keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
for (e in parse('analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R')) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
h <- dgp$designs$Hybrid
out <- rbindlist(lapply(seq_along(h$paths), function(k) {
  pi <- dgp$path_info(h$t_wk, h$e, h$paths[[k]])
  b <- dgp$build(pi, 'A', 0, 0)
  n <- b$n
  L <- matrix(0, n + 1, 2 + 3 * n); L[, 2] <- 1
  for (t in 1:n) L[t + 1, c(2 + t, 2 + n + t, 2 + 2 * n + t)] <- -1
  V <- L %*% b$Sigma %*% t(L)
  s00 <- V[1, 1]; sp0 <- V[-1, 1]; Spp <- V[-1, -1]
  on <- pi$on; n_on <- sum(on); n_off <- sum(!on)
  a <- ifelse(on, 1 / n_on, -1 / n_off)                       # E9 weights
  ab <- c(-1 / (n_off + 1), ifelse(on, 1 / n_on, -1 / (n_off + 1)))  # E9b
  Spp0 <- Spp - tcrossprod(sp0) / s00                          # conditional on BL
  data.table(path = k, n_on = n_on, n_off = n_off,
             gamma_range = paste(range(round(sp0 / s00, 6)), collapse = '..'),
             a_sp0 = sum(a * sp0),
             var_E9 = drop(t(a) %*% Spp %*% a),
             var_E9_given_BL = drop(t(a) %*% Spp0 %*% a),
             var_E9b = drop(t(ab) %*% V %*% ab),
             shared_S_cov = mean(Spp[upper.tri(Spp)]) - s00)
}))
print(out, digits = 4)
## Classic pre-post efficiencies for a between-participant contrast with
## one post visit, correlation r with baseline, equal variances.
r <- cov2cor(matrix(c(341.6, 341.6, 341.6, 599), 2))[1, 2]
cat(sprintf('\nBaseline-blinded visit correlation r = %.3f\n', r))
cat(sprintf('Relative variance: post only 1, change 2(1-r) = %.3f (equal-variance formula), ANCOVA 1-r^2 = %.3f\n',
            2 * (1 - r), 1 - r^2))
## exact change-score variance with unequal variances
cat(sprintf('Exact Var(post) = 599, Var(change) = %.1f, Var(post | BL) = %.1f\n',
            599 + 341.6 - 2 * 341.6, 599 - 341.6^2 / 341.6))
