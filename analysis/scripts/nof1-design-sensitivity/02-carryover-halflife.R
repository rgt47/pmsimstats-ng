#' S2. Carryover half-life sweep, crossed with the analysis half-life.
#'
#' Varies the data-generating half-life over nine values from 0 to 14
#' days. Each design is analyzed twice: with the binary on/off
#' indicator (analysis half-life 0) and with the decaying indicator
#' at the policy value of 0.2 week. The drug-response ceiling is set
#' per design to the common `delta_on`.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 2L
SWEEP_NAME <- 'carryover-halflife'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

halflife_days <- c(0, 0.5, 1, 2, 3, 5, 7, 10, 14)
analysis_t1half <- c(0, 0.2)

designs <- list(
  design_rct(), design_crossover(),
  design_nof1_alternating(), design_hybrid()
)
design_names <- c('RCT', 'XO', 'Alt', 'Hybrid')

model_params <- do.call(rbind, lapply(
  halflife_days / 7,
  function(h) baseline_model_param(N = BASE_N, carryover_t1half = h)
))

message('S2: ', nrow(model_params) * length(designs), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = model_params,
  n_reps = N_REPS,
  analysisparams = default_analysisparams(
    analysis_t1half = analysis_t1half)
)

res$results$halflife_days <- round(res$results$carryover_t1half * 7, 6)

checks <- function(res) {
  assert_stored_grid(res$results, 'carryover_t1half', halflife_days / 7)
  assert_stored_grid(res$results, 'analysis_t1half', analysis_t1half)
  assert_cell_count(res, nrow(model_params) * length(designs))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
