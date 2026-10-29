#' S5. Within-participant variability.
#'
#' Multiplies the SDs of all three components (natural history,
#' placebo response, drug response) by a common factor, so the
#' within-participant SD of the outcome scales with the factor. The
#' 2026-05 runs scaled the drug-response SD only. The drug-response
#' ceiling is set per design to the common `delta_on`.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 5L
SWEEP_NAME <- 'sigma2'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

scale_grid <- c(0.5, 0.75, 1, 1.25, 1.5)

designs <- list(design_hybrid(), design_nof1_alternating(), design_rct())
design_names <- c('Hybrid', 'Alt', 'RCT')

rp_list <- lapply(scale_grid, function(s) {
  rp <- baseline_resp_param()
  rp$sd <- rp$sd * s
  resp_set(rp, name = sprintf('sd_scale=%g', s))
})

message('S5: ', length(scale_grid) * length(designs), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = rp_list,
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = baseline_model_param(N = BASE_N),
  n_reps = N_REPS
)

res$results$sd_scale <- res$results$br_sd /
  baseline_resp_param()$sd[3]

checks <- function(res) {
  assert_stored_grid(res$results, 'sd_scale', scale_grid)
  assert_cell_count(res, length(scale_grid) * length(designs))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
