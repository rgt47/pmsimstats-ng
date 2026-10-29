#' S6. Number of periods at fixed trial length.
#'
#' The alternating N-of-1 design has 16 weekly visits and 2, 4, 8 or
#' 16 equal periods (1, 2, 4 or 8 on-off cycles), so trial length and
#' number of assessments are fixed and only the number of switches
#' varies. The parallel group has the same 16 weekly visits. The
#' drug-response ceiling is set per design to the common `delta_on`,
#' so the longer exposure of long periods is not an effect-size
#' difference. The 2026-05 runs varied the number of weekly visits,
#' which confounded cycles with trial length and effect.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 6L
SWEEP_NAME <- 'cycles'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

period_grid <- c(2L, 4L, 8L, 16L)

designs <- c(
  lapply(period_grid, design_nof1_periods),
  list(design_rct_weekly())
)
design_names <- c(sprintf('Alt-P%d', period_grid), 'RCT')

message('S6: ', length(designs), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = baseline_model_param(N = BASE_N),
  n_reps = N_REPS
)

res$results$n_periods <- ifelse(
  res$results$design_name == 'RCT', NA_integer_,
  res$results$n_switch + 1L
)

checks <- function(res) {
  assert_stored_grid(res$results, 'n_switch', c(0, period_grid - 1))
  assert_stored_grid(res$results, 'n_visits', 16)
  assert_cell_count(res, length(designs))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
