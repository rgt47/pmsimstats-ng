#' 00-smoke-test.R
#'
#' Single-cell end-to-end check of the pipeline swapped in from
#' pmsimstats-ng tidyverse (see review-pmsimstats-sync-2026-04-17-post.md
#' item R3). Runs before any of the sweep scripts. A pass here is
#' the minimum green light for S1-S12.
#'
#' On success: prints timing, dimensions of the results tibble, and
#' a one-line summary of the beta/p from a handful of replicates.
#' On failure: stops with an informative message.

source(here::here(
  "analysis", "scripts", "sensitivity", "00-common.R"
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
stopifnot("V1.br" %in% names(dat_single))

message(sprintf(
  "  dim(dat_single) = %d x %d; elapsed %.2fs",
  nrow(dat_single), ncol(dat_single),
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
))

message("Smoke test: lme_analysis on two-path data...")
dat <- dplyr::bind_rows(
  dat_single |> dplyr::mutate(path = 1),
  generate_data(mp, rp, bp, td$trialpaths[[2]],
                empirical = FALSE, make_positive_definite = TRUE,
                seed = 43, dgp_architecture = "mvn") |>
    dplyr::mutate(path = 2)
)

out <- lme_analysis(td$trialpaths, dat,
                    op = list(useDE = TRUE))
stopifnot(is.data.frame(out))
stopifnot(all(c("beta", "betaSE", "p") %in% names(out)))

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

message("Smoke test PASSED.")
