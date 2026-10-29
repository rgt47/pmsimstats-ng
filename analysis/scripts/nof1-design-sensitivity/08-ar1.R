#' S8. AR(1) serial-correlation sweep.
#'
#' Varies rho = c.tv = c.pb = c.br. The analysis always fits
#' corCAR1 (no with-versus-without comparison; see the Limitations).
#' Cells whose covariance matrix is not positive definite are
#' dropped, never repaired, and listed in the `infeasible` element
#' of the artifact with the smallest eigenvalue. The drug-response
#' ceiling is set per design to the common `delta_on`.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 8L
SWEEP_NAME <- 'ar1'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

rho_grid <- c(0, 0.3, 0.5, 0.7, 0.9)

designs <- list(design_hybrid(), design_nof1_alternating(), design_rct())
design_names <- c('Hybrid', 'Alt', 'RCT')

model_params <- do.call(rbind, lapply(rho_grid, function(rho) {
  baseline_model_param(N = BASE_N, c.tv = rho, c.pb = rho, c.br = rho)
}))

message('S8: ', nrow(model_params) * length(designs), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = model_params,
  n_reps = N_REPS
)

res$results$rho <- res$results$c.br

checks <- function(res) {
  run_rho <- unique(res$results$c.br)
  drop_rho <- rho_grid[res$infeasible$modelparamset]
  assert_stored_grid(data.frame(rho = c(run_rho, drop_rho)), 'rho',
                     rho_grid)
  assert_cell_count(res, nrow(model_params) * length(designs))
  stopifnot(!any(res$results$pd_repaired))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
