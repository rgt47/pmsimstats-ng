#' S9. Period length crossed with carryover half-life.
#'
#' Five period lengths (alternating N-of-1, eight visits, one per
#' period) crossed with three half-lives: 15 cells, each run at its
#' own half-life. The 2026-05 version built 15 designs but passed
#' only the first model-parameter row, so every cell ran at the
#' first half-life. The drug-response ceiling is set per design to
#' the common `delta_on`.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 9L
SWEEP_NAME <- 'period'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

period_days_grid <- c(3, 5, 7, 10, 14)
halflife_days_grid <- c(0.5, 3, 7)

designs <- lapply(period_days_grid, function(p) {
  weeks <- seq(p / 7, 8 * p / 7, by = p / 7)
  design_nof1_alternating(weeks = weeks)
})
design_names <- sprintf('Alt-%dd', period_days_grid)

model_params <- do.call(rbind, lapply(halflife_days_grid / 7, function(h) {
  baseline_model_param(N = BASE_N, carryover_t1half = h)
}))

message('S9: ', length(designs) * nrow(model_params), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = model_params,
  n_reps = N_REPS
)

res$results$period_days <-
  period_days_grid[match(res$results$design_name, design_names)]
res$results$halflife_days <- round(res$results$carryover_t1half * 7, 6)
res$results$ratio <- res$results$halflife_days / res$results$period_days

checks <- function(res) {
  assert_stored_grid(res$results, 'carryover_t1half',
                     halflife_days_grid / 7)
  assert_stored_grid(res$results, 'period_days', period_days_grid)
  assert_stored_grid(res$results, 'N', BASE_N)
  assert_cell_count(res, length(designs) * nrow(model_params))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
