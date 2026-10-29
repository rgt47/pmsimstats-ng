#' S7. Total participant count.
#'
#' Varies N over six values. The drug-response ceiling is set per
#' design to the common `delta_on`.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 7L
SWEEP_NAME <- 'patients'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

N_grid <- c(16, 24, 35, 60, 80, 120)

designs <- list(design_hybrid(), design_nof1_alternating(), design_rct())
design_names <- c('Hybrid', 'Alt', 'RCT')

model_params <- do.call(rbind, lapply(N_grid, function(n) {
  baseline_model_param(N = n)
}))

message('S7: ', nrow(model_params) * length(designs), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = model_params,
  n_reps = N_REPS
)

res$results$n_total <- res$results$N

checks <- function(res) {
  assert_stored_grid(res$results, 'N', N_grid)
  assert_cell_count(res, nrow(model_params) * length(designs))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
