#' S10. Biomarker-response correlation crossed with the DGP
#' architecture, for the main-effect target.
#'
#' The sweep stores the main-effect coefficient (`Dbc`), not the
#' biomarker-by-drug interaction. Four values of c_bm, two carryover
#' settings and the two encodings of the interaction (the correlation
#' encoding `mvn` and the mean-shift encoding `mean_moderation`).
#' Cells whose covariance matrix is not positive definite are dropped,
#' never repaired, and listed in `infeasible` with the smallest
#' eigenvalue. The drug-response ceiling is set per design to the
#' common `delta_on`.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 10L
SWEEP_NAME <- 'biomarker'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

c_bm_grid <- c(0, 0.2, 0.4, 0.6)
halflife_days_grid <- c(0, 3)
arch_grid <- c('mvn', 'mean_moderation')

designs <- list(design_rct(), design_nof1_alternating(), design_hybrid())
design_names <- c('RCT', 'Alt', 'Hybrid')

model_params <- do.call(rbind, lapply(c_bm_grid, function(cb) {
  do.call(rbind, lapply(halflife_days_grid, function(hd) {
    baseline_model_param(
      N = BASE_N, c.bm = cb, carryover_t1half = hd / 7
    )
  }))
}))

parts <- lapply(arch_grid, function(arch) {
  message('S10: dgp_architecture = ', arch)
  sweep_call(
    designs, design_names,
    respparamsets = list(resp_set()),
    blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
    modelparams = model_params,
    n_reps = N_REPS,
    dgp_architecture = arch
  )
})
res <- combine_sweeps(parts)

res$results$c_bm <- res$results$c.bm
res$results$halflife_days <- round(res$results$carryover_t1half * 7, 6)

checks <- function(res) {
  run_cbm <- unique(res$results$c.bm)
  drop_cbm <- model_params$c.bm[res$infeasible$modelparamset]
  assert_stored_grid(data.frame(c_bm = c(run_cbm, drop_cbm)), 'c_bm',
                     c_bm_grid)
  assert_stored_grid(res$results, 'carryover_t1half',
                     halflife_days_grid / 7)
  assert_stored_grid(res$results, 'dgp_architecture', arch_grid)
  stopifnot(!any(res$results$pd_repaired))
  n_run <- nrow(dplyr::distinct(
    res$results, trialdesign, modelparamset, dgp_architecture))
  stopifnot(n_run + nrow(res$infeasible) ==
              nrow(model_params) * length(designs) * length(arch_grid))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
