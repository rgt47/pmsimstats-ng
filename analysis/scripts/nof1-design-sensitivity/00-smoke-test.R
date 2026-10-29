#' 00-smoke-test.R
#'
#' Single-cell end-to-end check of the engine and the shared helpers.
#' Runs before any of the sweep scripts. A pass here is the minimum
#' green light for S1-S12. It is a structural check, not a test of
#' the statistical properties (see 14-rct-null-check.R for those).
#'
#' On success: prints timing, dimensions of the results tibble, and
#' a one-line summary of the beta/p from a handful of replicates.
#' On failure: stops with an informative message.

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

set.seed(42)

message("Smoke test: building fixtures...")

td <- design_nof1()
mp <- baseline_model_param(N = 20)
rp <- baseline_resp_param()
bp <- baseline_bl_param()

message("Smoke test: generate_data on single path (mvn)...")
t0 <- Sys.time()
dat_single <- generate_data(
  mp, rp, bp,
  td$trialpaths[[1]],
  empirical = FALSE, make_positive_definite = TRUE,
  seed = 42, dgp_architecture = "mvn"
)
stopifnot(is.data.frame(dat_single))
stopifnot("ptID" %in% names(dat_single))
stopifnot("bm" %in% names(dat_single))
stopifnot("OL1.br" %in% names(dat_single))

message(sprintf(
  "  dim(dat_single) = %d x %d; elapsed %.2fs",
  nrow(dat_single), ncol(dat_single),
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))

message("Smoke test: lme_analysis on four-path data...")
dat <- dplyr::bind_rows(lapply(seq_along(td$trialpaths), function(i) {
  generate_data(mp, rp, bp, td$trialpaths[[i]],
                empirical = FALSE, make_positive_definite = TRUE,
                seed = 42 + i, dgp_architecture = "mvn") |>
    dplyr::mutate(path = i)
}))

out <- lme_analysis(td$trialpaths, dat,
                    op = list(useDE = TRUE, target_coef = "main"))
stopifnot(is.data.frame(out))
stopifnot(all(c("beta", "betaSE", "p") %in% names(out)))
stopifnot(!is.na(out$beta))

message(sprintf(
  "  beta = %.3f; betaSE = %.3f; p = %.4g",
  out$beta, out$betaSE, out$p
))

message("Smoke test: generate_simulated_results (Nreps = 3)...")

t0 <- Sys.time()
res <- generate_simulated_results(
  trialdesigns = list(td),
  respparamsets = wrap_param_set(rp, "resp_default"),
  blparamsets = wrap_param_set(bp, "bl_default"),
  censorparams = NA,
  modelparams = mp,
  simparam = default_simparam(n_reps = 3),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = 1
)
stopifnot("results" %in% names(res))
stopifnot(nrow(res$results) >= 3)

message(sprintf(
  "  generate_simulated_results OK: %d rows, elapsed %.2fs",
  nrow(res$results),
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))

stopifnot(all(c("fit_formula", "delta_on", "min_eig", "pd_repaired",
                "analysis_t1half") %in% names(res$results)))
stopifnot(!any(res$results$pd_repaired))

message("Smoke test PASSED.")
