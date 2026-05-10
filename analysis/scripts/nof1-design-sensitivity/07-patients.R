#' S7. Patient-count (n) sweep.
#'
#' Precedent: Senn 2018; Blackston 2019.
#' Varies N-of-1 patient count and RCT participant count
#' concurrently so that the total number of observations per cell
#' is matched between designs.

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

set.seed(42 + 7)

SWEEP_ID <- 7L
SWEEP_NAME <- "patients"

N_grid <- c(16, 24, 35, 60, 80, 120)

td_nof1 <- design_nof1()
td_rct <- design_rct()
rp <- baseline_resp_param()
bp <- baseline_bl_param()

model_params <- do.call(rbind, lapply(N_grid, function(N) {
  baseline_model_param(N = N)
}))

message("S7: running ", nrow(model_params) * 2, " cells...")

res <- generate_simulated_results(
  trialdesigns = list(td_nof1, td_rct),
  respparamsets = wrap_param_set(rp, "resp_default"),
  blparamsets = wrap_param_set(bp, "bl_default"),
  censorparams = NA,
  modelparams = model_params,
  simparam = default_simparam(n_reps = 500),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = max(1, parallel::detectCores() - 1)
)

res$results$n_total <- res$results$N
res$results$design_name <- c("Nof1", "RCT")[res$results$trialdesign]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
