#' 14. Check of the parallel-group analyses under the null and an
#' alternative.
#'
#' The referee report (comment 4.5) found that the mixed-model
#' analysis of the parallel group (a baseline row at t = 0 with the
#' drug indicator zero, one shared linear time term) has a nonzero
#' mean under the null and a miscalibrated SE. This script compares
#' that analysis ('legacy') with the baseline-adjusted ANCOVA arm
#' contrast ('ancova') at three AR(1) values, for the blinded and the
#' open-label parallel group. The zero-mean criterion is
#' |mean estimate| < 2 MCSE under the null. The calibration criterion is
#' kappa (model SE / empirical SD) within 5% of 1.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))
source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-summarize.R'
))

SWEEP_ID <- 14L
SWEEP_NAME <- 'rct-null-check'
N_REPS <- sens_reps(1000L)
set.seed(42 + SWEEP_ID)

rho_grid <- c(0, 0.3, 0.7)
delta_grid <- c(0, BASE_DELTA)

designs <- list(design_rct(), design_rct_openlabel())
design_names <- c('RCT', 'RCT-OL')

model_params <- do.call(rbind, lapply(rho_grid, function(rho) {
  baseline_model_param(N = BASE_N, c.tv = rho, c.pb = rho, c.br = rho)
}))

ap <- rbind(
  transform(default_analysisparams(analysis_t1half = 0),
            parallel_analysis = 'legacy'),
  default_analysisparams(analysis_t1half = 0)
)

res <- sweep_call(
  designs, design_names,
  respparamsets = lapply(delta_grid, function(d) resp_set(
    target_delta = d)),
  blparamsets = wrap_param_set(baseline_bl_param(), 'bl_default'),
  modelparams = model_params,
  n_reps = N_REPS,
  analysisparams = ap
)
res$results$analysis <- c('legacy', 'ancova')[res$results$analysis_row]
res$results$rho <- res$results$c.br

checks <- function(res) {
  assert_stored_grid(res$results, 'c.br', rho_grid)
  assert_stored_grid(res$results, 'delta_on', delta_grid)
  assert_cell_count(res, length(designs) * length(delta_grid) *
                      length(rho_grid))
}
save_sweep(res, SWEEP_ID, SWEEP_NAME, n_reps = N_REPS,
           checks = checks)

tab <- summarize_cells(
  res$results, c('design_name', 'analysis', 'delta_on', 'rho')
)
tab$z_estimand <- tab$bias_estimand / tab$mcse_bias
print(as.data.frame(tab[, c(
  'design_name', 'analysis', 'delta_on', 'rho', 'n_fit', 'mean_beta',
  'bias_estimand', 'mcse_bias', 'z_estimand', 'bias', 'emp_se',
  'model_se', 'kappa', 'power', 'coverage')]), digits = 3)
