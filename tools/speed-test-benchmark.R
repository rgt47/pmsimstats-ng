#!/usr/bin/env Rscript
# Cross-host compute speed benchmark. Pure base R plus the `parallel`
# package (bundled with every R install, not a CRAN dependency), so
# the same script runs identically via plain `Rscript` on julia,
# zeus, and the Expanse login node. Defaults to single-threaded, for
# a per-core comparison; pass a second argument above 1 to compare
# horizontal scale instead (relevant on Expanse, where per-core speed
# lags but core count does not).
#
# Usage: Rscript speed-test-benchmark.R <n_reps> [n_cores]

library(parallel)

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 20000L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1L
n_cores <- max(1L, min(n_cores, detectCores()))

n_units <- 40L
n_time <- 30L

run_one_rep <- function(seed) {
  set.seed(seed)
  trt <- rep(c(0, 1), length.out = n_time)
  subject <- rep(seq_len(n_units), each = n_time)
  time <- rep(seq_len(n_time), times = n_units)
  trt_full <- rep(trt, times = n_units)
  baseline <- rnorm(n_units, mean = 50, sd = 10)[subject]
  noise <- rnorm(n_units * n_time, sd = 5)
  y <- baseline + 2.5 * trt_full + 0.1 * time + noise
  fit <- lm(y ~ trt_full + time + factor(subject))
  unname(coef(fit)["trt_full"])
}

start <- Sys.time()
if (n_cores > 1) {
  estimates <- unlist(mclapply(seq_len(n_reps), run_one_rep,
                                mc.cores = n_cores))
} else {
  estimates <- vapply(seq_len(n_reps), run_one_rep, numeric(1))
}
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

cat(sprintf("host: %s\n", Sys.info()[["nodename"]]))
cat(sprintf("r_version: %s\n", R.version.string))
cat(sprintf("n_reps: %d\n", n_reps))
cat(sprintf("n_cores: %d\n", n_cores))
cat(sprintf("elapsed_sec: %.2f\n", elapsed))
cat(sprintf("reps_per_sec: %.2f\n", n_reps / elapsed))
cat(sprintf("mean_estimate: %.4f\n", mean(estimates)))
