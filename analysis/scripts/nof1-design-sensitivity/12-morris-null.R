#' S12. Null-effect calibration with 2,000 replicates per cell.
#'
#' The drug-response ceiling is 0 (the biomarker-response correlation
#' c_bm keeps its baseline value of 0.3). Four configurations of the
#' baseline SD, the drug-response SD and the carryover half-life, four
#' designs, and two analysis half-lives (0 and 0.2 week).

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 12L
SWEEP_NAME <- 'morris-null'
N_REPS <- sens_reps(2000L)
set.seed(42 + SWEEP_ID)

configs <- list(
  list(bl_sd = 1.5, br_sd = 1.0, halflife_days = 0),
  list(bl_sd = 1.5, br_sd = 1.0, halflife_days = 3),
  list(bl_sd = 0.75, br_sd = 1.5, halflife_days = 3),
  list(bl_sd = 2.0, br_sd = 0.5, halflife_days = 3)
)
analysis_t1half <- c(0, 0.2)

designs <- list(
  design_rct(), design_crossover(),
  design_nof1_alternating(), design_hybrid()
)
design_names <- c('RCT', 'XO', 'Alt', 'Hybrid')

parts <- lapply(seq_along(configs), function(i) {
  cfg <- configs[[i]]
  message(sprintf(
    'S12 config %d/%d: bl_sd=%g br_sd=%g halflife_days=%g',
    i, length(configs), cfg$bl_sd, cfg$br_sd, cfg$halflife_days
  ))
  bp <- baseline_bl_param()
  bp$sd[bp$cat == 'BL'] <- cfg$bl_sd
  rp <- baseline_resp_param()
  rp$sd[rp$cat == 'br'] <- cfg$br_sd
  sweep_call(
    designs, design_names,
    respparamsets = list(resp_set(rp, target_delta = 0)),
    blparamsets = wrap_param_set(bp, 'bl_config'),
    modelparams = baseline_model_param(
      N = BASE_N, carryover_t1half = cfg$halflife_days / 7),
    n_reps = N_REPS,
    analysisparams = default_analysisparams(
      analysis_t1half = analysis_t1half)
  )
})
res <- combine_sweeps(parts, tag = 'config')
res$results$config <- as.integer(res$results$config)

checks <- function(res) {
  assert_stored_grid(res$results, 'bl_sd',
                     vapply(configs, `[[`, 0, 'bl_sd'))
  assert_stored_grid(res$results, 'br_sd',
                     vapply(configs, `[[`, 0, 'br_sd'))
  assert_stored_grid(res$results, 'carryover_t1half',
                     vapply(configs, `[[`, 0, 'halflife_days') / 7)
  assert_stored_grid(res$results, 'delta_on', 0)
  stopifnot(nrow(dplyr::distinct(res$results, config, trialdesign)) ==
              length(configs) * length(designs))
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
