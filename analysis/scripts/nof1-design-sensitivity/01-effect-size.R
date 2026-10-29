#' S1. Effect-size sweep on the common estimand.
#'
#' Varies the mean drug response on drug (`delta_on`, outcome units)
#' over a seven-point grid. The drug-response ceiling is set per
#' design so that every design has the same `delta_on`; the grid is
#' therefore a grid of one estimand, not of design-specific ceilings.
#' The null cell (delta_on = 0) gives Type I error, bias and
#' calibration for every design. RCT-OL is the parallel group with
#' open-label expectancy (the assumption of the 2026-05 runs).

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 1L
SWEEP_NAME <- 'effect-size'
N_REPS <- sens_reps(1000L)
set.seed(42 + SWEEP_ID)

delta_grid <- c(0, -0.25, -0.5, -0.75, -1, -1.5, -2)
designs <- list(
  design_hybrid(), design_compact2(), design_nof1_alternating(),
  design_crossover(), design_rct(), design_rct_openlabel()
)
design_names <- c('Hybrid', 'Compact2', 'Alt', 'XO', 'RCT', 'RCT-OL')

message('S1: ', length(delta_grid) * length(designs), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = lapply(delta_grid, function(d) resp_set(
    target_delta = d)),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = baseline_model_param(N = BASE_N),
  n_reps = N_REPS
)

checks <- function(res) {
  assert_stored_grid(res$results, 'delta_on', delta_grid)
  assert_stored_grid(res$results, 'N', BASE_N)
  assert_cell_count(res, length(delta_grid) * length(designs))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
