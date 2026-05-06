#' S2. Carryover half-life sweep.
#'
#' Precedent: Blackston 2019; Hendrickson 2020.
#' Varies carryover half-life across nine values including
#' beyond-period half-lives (10 and 14 days; 8-week periods are 1
#' week each). Four designs are compared.

source(here::here(
  "analysis", "scripts", "sensitivity", "00-common.R"
))

set.seed(42 + 2)

SWEEP_ID <- 2L
SWEEP_NAME <- "carryover-halflife"

halflife_grid_days <- c(0, 0.5, 1, 2, 3, 5, 7, 10, 14)
halflife_grid_weeks <- halflife_grid_days / 7

designs <- list(
  design_rct(),
  design_crossover(),
  design_nof1_alternating(),
  design_hybrid()
)

rp <- baseline_resp_param()
bp <- baseline_bl_param()

model_params <- do.call(rbind, lapply(
  halflife_grid_weeks,
  function(h) baseline_model_param(carryover_t1half = h)
))

message("S2: running ", nrow(model_params) * length(designs),
        " cells across 4 designs...")

res <- generate_simulated_results(
  trialdesigns = designs,
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

res$results$halflife_days <-
  halflife_grid_days[res$results$modelparamset]
res$results$design_name <-
  c("RCT", "XO", "Nof1-alt", "Hybrid")[res$results$trialdesign]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
