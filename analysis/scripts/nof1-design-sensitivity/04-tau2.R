#' S4. Treatment-effect heterogeneity and baseline SD.
#'
#' Two one-axis blocks, each with the drug-response ceiling set per
#' design to the common `delta_on`.
#'
#'   tau    SD of participant-level heterogeneity of the drug effect
#'          (`tau` in the model parameters; each participant's drug
#'          response is shifted by v_i times the unit response curve,
#'          v_i ~ N(0, tau^2)). This is the patient-by-treatment
#'          variance component of Senn (2019). The 2026-05 runs had
#'          no such component.
#'   bl_sd  SD of the baseline symptom score, a participant-level
#'          random intercept that enters every visit equally. It is
#'          not treatment heterogeneity.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

SWEEP_ID <- 4L
SWEEP_NAME <- 'tau2'
N_REPS <- sens_reps(500L)
set.seed(42 + SWEEP_ID)

tau_grid <- c(0, 0.25, 0.5, 0.75, 1)
bl_sd_grid <- c(0.5, 1, 1.5, 2, 3)

designs <- list(design_hybrid(), design_nof1_alternating(), design_rct())
design_names <- c('Hybrid', 'Alt', 'RCT')

message('S4: ', (length(tau_grid) + length(bl_sd_grid)) *
          length(designs), ' cells')

res_tau <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = do.call(rbind, lapply(tau_grid, function(t) {
    baseline_model_param(N = BASE_N, tau = t)
  })),
  n_reps = N_REPS
)

bp_list <- lapply(bl_sd_grid, function(s) {
  bp <- baseline_bl_param()
  bp$sd[bp$cat == 'BL'] <- s
  list(name = sprintf('bl_sd=%g', s), param = bp)
})
res_bl <- sweep_call(
  designs, design_names,
  respparamsets = list(resp_set()),
  blparamsets = bp_list,
  modelparams = baseline_model_param(N = BASE_N),
  n_reps = N_REPS
)

res <- combine_sweeps(list(tau = res_tau, bl_sd = res_bl), tag = 'axis')

checks <- function(res) {
  tau_rows <- res$results$axis == 'tau'
  assert_stored_grid(res$results, 'tau', tau_grid, among = tau_rows)
  assert_stored_grid(res$results, 'bl_sd', bl_sd_grid,
                     among = !tau_rows)
  stopifnot(
    nrow(dplyr::distinct(res$results[tau_rows, ], trialdesign,
                         modelparamset)) ==
      length(tau_grid) * length(designs),
    nrow(dplyr::distinct(res$results[!tau_rows, ], trialdesign,
                         blparamset)) ==
      length(bl_sd_grid) * length(designs)
  )
}

save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)
