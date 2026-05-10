library(tinytest)

#---------------------------------------------------------------------
# Type I error check under the null (c.bm = 0)
#
# Runs a modest batch of replicates through generate_data and
# lme_analysis on a hybrid design at c.bm = 0, under both
# architectures, and asserts that the empirical rejection rate at
# alpha = 0.05 is in [0.02, 0.10]. The upper bound is generous
# (the tolerance at 100 replicates is wide); assertion tightens
# naturally at higher n_reps.
#
# This test is moderately expensive (seconds) and is skipped when
# the environment variable NOF1POWER_SKIP_SLOW is set.
#---------------------------------------------------------------------

if (nzchar(Sys.getenv("NOF1POWER_SKIP_SLOW"))) {
  exit_file("slow tests skipped (NOF1POWER_SKIP_SLOW set)")
}

set.seed(2026)

td <- nof1power::build_trial_design(
  name_longform = "Hybrid",
  name_shortform = "H",
  timepoints = c(2, 4, 6, 8),
  expectancies = c(0.5, 0.5, 0, 0),
  ondrug = list(c(1, 1, 0, 0), c(0, 0, 1, 1))
)

model_param_null <- data.frame(
  N = 40, c.bm = 0, carryover_t1half = 0,
  c.tv = 0.6, c.pb = 0.6, c.br = 0.6,
  c.cf1t = 0.3, c.cfct = 0.1,
  stringsAsFactors = FALSE
)

resp_param <- data.frame(
  cat = c("tv", "pb", "br"),
  max = c(10, 5, 8), disp = c(1, 1, 1),
  rate = c(0.5, 0.5, 0.5), sd = c(2, 1, 3),
  stringsAsFactors = FALSE
)

baseline_param <- data.frame(
  cat = c("bm", "BL"),
  m = c(120, 50), sd = c(15, 10),
  stringsAsFactors = FALSE
)

N_REPS <- 100
ARCHES <- c("mvn", "mean_moderation")
target_lo <- 0.02
target_hi <- 0.10

for (arch in ARCHES) {
  p_values <- numeric(N_REPS)
  for (i in seq_len(N_REPS)) {
    # Split participants evenly across the two paths.
    paths <- td$trialpaths
    n_paths <- length(paths)
    N_total <- model_param_null$N
    per_path <- N_total %/% n_paths
    leftover <- N_total %% n_paths
    ns <- per_path + c(rep(1, leftover),
                       rep(0, n_paths - leftover))

    dat_paths <- purrr::map2(paths, ns, function(p, n) {
      mp <- model_param_null
      mp$N <- n
      d <- nof1power::generate_data(
        mp, resp_param, baseline_param, p,
        empirical = FALSE, make_positive_definite = TRUE,
        dgp_architecture = arch
      )
      d
    })
    for (j in seq_along(dat_paths)) {
      dat_paths[[j]]$path <- j
    }
    dat <- do.call(rbind, dat_paths)

    out <- nof1power::lme_analysis(paths, dat, op = list(useDE = TRUE))
    p_values[i] <- out$p
  }
  rejection_rate <- mean(p_values < 0.05, na.rm = TRUE)
  expect_true(
    rejection_rate >= target_lo && rejection_rate <= target_hi,
    info = sprintf(
      "Architecture %s: empirical rejection rate %.3f not in [%.2f, %.2f]",
      arch, rejection_rate, target_lo, target_hi
    )
  )
}
