#' S11. Misspecification factorial: data-generating carryover and
#' AR(1), crossed with the analysis half-life.
#'
#' DGP carryover (off, 3 days) by DGP AR(1) (off, rho = 0.7), at the
#' null and at the common baseline effect, for three designs. Each
#' dataset is analyzed with the analysis half-life 0 (binary
#' indicator), 0.2 week (policy value) and 3/7 week (the
#' data-generating half-life of the carryover cells).

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 11L
SWEEP_NAME <- 'misspec'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

grid <- expand.grid(
  carryover = c('off', 'on'),
  ar1 = c('off', 'high'),
  stringsAsFactors = FALSE
)
analysis_t1half <- c(0, 0.2, 3 / 7)

designs <- list(design_hybrid(), design_nof1_alternating(), design_rct())
design_names <- c('Hybrid', 'Alt', 'RCT')

model_params <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  t_half <- if (grid$carryover[i] == 'on') 3 / 7 else 0
  rho <- if (grid$ar1[i] == 'high') 0.7 else 0
  baseline_model_param(
    N = BASE_N, carryover_t1half = t_half,
    c.tv = rho, c.pb = rho, c.br = rho
  )
}))

message('S11: ', nrow(grid) * length(designs) * 2, ' cells')

parts <- lapply(c(null = 0, alt = BASE_DELTA), function(delta) {
  sweep_call(
    designs, design_names,
    respparamsets = list(resp_set(target_delta = delta)),
    blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
    modelparams = model_params,
    n_reps = N_REPS,
    analysisparams = default_analysisparams(
      analysis_t1half = analysis_t1half)
  )
})
res <- combine_sweeps(parts, tag = 'delta_label')

res$results$carryover <- grid$carryover[res$results$modelparamset]
res$results$ar1 <- grid$ar1[res$results$modelparamset]

checks <- function(res) {
  assert_stored_grid(res$results, 'carryover_t1half', c(0, 3 / 7))
  assert_stored_grid(res$results, 'c.br', c(0, 0.7))
  assert_stored_grid(res$results, 'analysis_t1half', analysis_t1half)
  assert_stored_grid(res$results, 'delta_on', c(0, BASE_DELTA))
  stopifnot(!any(res$results$pd_repaired))
  n_run <- nrow(dplyr::distinct(
    res$results, delta_label, trialdesign, modelparamset))
  stopifnot(n_run + nrow(res$infeasible) ==
              nrow(grid) * length(designs) * 2)
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
