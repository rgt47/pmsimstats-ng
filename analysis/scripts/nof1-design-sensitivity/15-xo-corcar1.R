#' S12b. Two-period crossover and alternating design under the null,
#' with and without the corCAR1 residual correlation.
#'
#' The referee report (comment 5.23) asked whether the Type I
#' inflation of the two-period crossover in S12 comes from fitting
#' corCAR1 to two follow-up occasions per participant. The
#' crossover and the alternating design are refitted with the random
#' intercept only. Null effect, baseline configuration, 2,000
#' replicates, analysis half-life 0.2 week.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 15L
SWEEP_NAME <- 'xo-corcar1'
N_REPS <- sens_reps(2000L)
set.seed(42 + SWEEP_ID)

designs <- list(design_crossover(), design_nof1_alternating())
design_names <- c('XO', 'Alt')

ap <- rbind(
  transform(default_analysisparams(), use_corCAR1 = TRUE),
  transform(default_analysisparams(), use_corCAR1 = FALSE)
)

res <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set(target_delta = 0)),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = baseline_model_param(N = BASE_N),
  n_reps = N_REPS,
  analysisparams = ap
)
res$results$corCAR1 <- c(TRUE, FALSE)[res$results$analysis_row]

checks <- function(res) {
  assert_stored_grid(res$results, 'delta_on', 0)
  stopifnot(setequal(unique(res$results$corCAR1), c(TRUE, FALSE)))
  assert_cell_count(res, length(designs))
}
save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
