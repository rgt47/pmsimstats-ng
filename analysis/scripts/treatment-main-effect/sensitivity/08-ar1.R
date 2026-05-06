#' S8. AR(1) serial-correlation sweep (with/without corCAR1).
#'
#' Precedent: Wang-Schork 2019.
#' Varies the within-factor AR(1) correlations rho = c.tv = c.pb =
#' c.br. Crossed with the analysis-model switch: fit with vs
#' without corCAR1(form = ~t | ptID). Quantifies the cost of
#' ignoring serial correlation.
#'
#' Note: in the current lme_analysis implementation, corCAR1 is
#' enabled by default. Fitting without corCAR1 requires a small
#' analysis-level hook that is not yet exposed. This script runs
#' the DGP side of the sweep and leaves analysis toggling for
#' a follow-up once the `op$use_corCAR1` option is added.

source(here::here(
  "analysis", "scripts", "sensitivity", "00-common.R"
))

set.seed(42 + 8)

SWEEP_ID <- 8L
SWEEP_NAME <- "ar1"

rho_grid <- c(0, 0.3, 0.5, 0.7, 0.9)

designs <- list(design_nof1(), design_rct())
rp <- baseline_resp_param()
bp <- baseline_bl_param()

model_params <- do.call(rbind, lapply(rho_grid, function(rho) {
  mp <- baseline_model_param()
  mp$c.tv <- rho
  mp$c.pb <- rho
  mp$c.br <- rho
  mp
}))

message("S8: running ", nrow(model_params) * length(designs),
        " AR(1) cells...")

res <- generate_simulated_results(
  trialdesigns = designs,
  respparamsets = wrap_param_set(rp, "resp_default"),
  blparamsets = wrap_param_set(bp, "bl_default"),
  censorparams = NA,
  modelparams = model_params,
  simparam = default_simparam(n_reps = 500),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = max(1, parallel::detectCores() - 1)
)

res$results$rho <- rho_grid[res$results$modelparamset]
res$results$design_name <- c("Nof1", "RCT")[res$results$trialdesign]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
