#' S9. Period-length sweep (crossed with carryover half-life).
#'
#' Precedent: Wang-Schork 2019.
#' Varies the N-of-1 period length (days between measurements) and
#' the carryover half-life jointly; the key quantity is the ratio
#' t_half / period.

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

set.seed(42 + 9)

SWEEP_ID <- 9L
SWEEP_NAME <- "period"

period_days_grid <- c(3, 5, 7, 10, 14)
halflife_days_grid <- c(0.5, 3, 7)

grid <- expand.grid(
  period_days = period_days_grid,
  halflife_days = halflife_days_grid
)

td_list <- lapply(seq_len(nrow(grid)), function(i) {
  p <- grid$period_days[i]
  weeks <- seq(p / 7, 8 * p / 7, by = p / 7)
  design_nof1_alternating(weeks = weeks)
})

rp <- baseline_resp_param()
bp <- baseline_bl_param()

model_params <- data.frame(
  N = 40, c.bm = 0.3,
  carryover_t1half = grid$halflife_days / 7,
  c.tv = 0.3, c.pb = 0.3, c.br = 0.3,
  c.cf1t = 0.1, c.cfct = 0.05,
  carryover_form = "exponential",
  weibull_shape = 1,
  stringsAsFactors = FALSE
)

message("S9: running ", nrow(grid), " period-by-halflife cells...")

res <- generate_simulated_results(
  trialdesigns = td_list,
  respparamsets = wrap_param_set(rp, "resp_default"),
  blparamsets = wrap_param_set(bp, "bl_default"),
  censorparams = NA,
  modelparams = model_params[1, , drop = FALSE],
  simparam = default_simparam(n_reps = 300),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = max(1, parallel::detectCores() - 1)
)

res$results$period_days <- grid$period_days[res$results$trialdesign]
res$results$halflife_days <-
  grid$halflife_days[res$results$trialdesign]
res$results$ratio <-
  res$results$halflife_days / res$results$period_days

save_sweep(res, SWEEP_ID, SWEEP_NAME)
