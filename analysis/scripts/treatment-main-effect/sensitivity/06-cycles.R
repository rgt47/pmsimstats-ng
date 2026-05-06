#' S6. Cycles-per-patient (k) sweep for N-of-1.
#'
#' Precedent: Senn 2018.
#' Varies the number of measurement timepoints per patient k by
#' reconstructing the N-of-1 design for each k.
#' Senn's variance decomposition predicts concave gains in k.

source(here::here(
  "analysis", "scripts", "sensitivity", "00-common.R"
))

set.seed(42 + 6)

SWEEP_ID <- 6L
SWEEP_NAME <- "cycles"

k_grid <- c(2, 4, 6, 8, 12, 16)

build_designs_by_k <- function(k_values) {
  lapply(k_values, function(k) {
    weeks <- seq_len(k)
    design_nof1_alternating(weeks = weeks)
  })
}

build_rcts_by_k <- function(k_values) {
  lapply(k_values, function(k) {
    weeks <- seq_len(k)
    design_rct(weeks = weeks)
  })
}

designs <- c(
  build_designs_by_k(k_grid),
  build_rcts_by_k(k_grid)
)
design_labels <- c(
  sprintf("Nof1_k=%d", k_grid),
  sprintf("RCT_k=%d", k_grid)
)

rp <- baseline_resp_param()
bp <- baseline_bl_param()
mp <- baseline_model_param()

message("S6: running ", length(designs), " design-by-k cells...")

res <- generate_simulated_results(
  trialdesigns = designs,
  respparamsets = wrap_param_set(rp, "resp_default"),
  blparamsets = wrap_param_set(bp, "bl_default"),
  censorparams = NA,
  modelparams = mp,
  simparam = default_simparam(n_reps = 500),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = max(1, parallel::detectCores() - 1)
)

res$results$design_k <- design_labels[res$results$trialdesign]
res$results$design_name <- ifelse(
  res$results$trialdesign <= length(k_grid), "Nof1-alt", "RCT"
)
res$results$k <- k_grid[
  ((res$results$trialdesign - 1) %% length(k_grid)) + 1
]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
