#' S3. Decay form of the data-generating carryover, crossed with the
#' analysis half-life.
#'
#' The data-generating carryover takes one of five forms. The analysis
#' always builds the decaying indicator from an exponential, at one of
#' three half-lives: 0 (binary indicator), 0.2 week (policy value) and
#' 3/7 week (the data-generating half-life). Only the exponential DGP
#' form is matched in form at 3/7 week.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 3L
SWEEP_NAME <- 'decay-form'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

form_grid <- list(
  list(form = 'exponential', shape = 1),
  list(form = 'linear', shape = 1),
  list(form = 'weibull', shape = 0.8),
  list(form = 'weibull', shape = 1.5),
  list(form = 'power', shape = 1.5)
)
form_labels <- vapply(form_grid, function(fg) {
  if (fg$form == 'weibull') sprintf('weibull_k=%g', fg$shape) else fg$form
}, character(1))
analysis_t1half <- c(0, 0.2, 3 / 7)

designs <- list(design_nof1_alternating(), design_hybrid(), design_rct())
design_names <- c('Alt', 'Hybrid', 'RCT')

model_params <- do.call(rbind, lapply(form_grid, function(fg) {
  baseline_model_param(
    N = BASE_N, carryover_t1half = 3 / 7,
    carryover_form = fg$form, weibull_shape = fg$shape
  )
}))

message('S3: ', nrow(model_params) * length(designs), ' cells')

res <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = model_params,
  n_reps = N_REPS,
  analysisparams = default_analysisparams(
    analysis_t1half = analysis_t1half)
)

res$results$decay_form <- ifelse(
  res$results$carryover_form == 'weibull',
  sprintf('weibull_k=%g', res$results$weibull_shape),
  res$results$carryover_form
)

checks <- function(res) {
  assert_stored_grid(res$results, 'decay_form', form_labels)
  assert_stored_grid(res$results, 'analysis_t1half', analysis_t1half)
  assert_cell_count(res, nrow(model_params) * length(designs))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
