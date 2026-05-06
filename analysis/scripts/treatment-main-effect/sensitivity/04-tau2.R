#' S4. Between-patient heterogeneity (tau^2) sweep.
#'
#' Precedent: Senn 2018; Zucker 2010.
#' Varies the biomarker (bm) baseline SD, which in pmsimstats
#' short-form encodes the between-patient heterogeneity component
#' that drives parallel-RCT inefficiency in the Senn variance
#' decomposition.

source(here::here(
  "analysis", "scripts", "sensitivity", "00-common.R"
))

set.seed(42 + 4)

SWEEP_ID <- 4L
SWEEP_NAME <- "tau2"

tau_grid <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)

designs <- list(design_nof1(), design_rct())
rp <- baseline_resp_param()
bp_base <- baseline_bl_param()

bp_list <- lapply(tau_grid, function(t) {
  bp <- bp_base
  bp$sd[bp$cat == "BL"] <- t
  list(name = sprintf("tau=%s", t), param = bp)
})

mp <- baseline_model_param()

message("S4: running ", length(tau_grid) * length(designs),
        " cells across tau^2 levels...")

res <- generate_simulated_results(
  trialdesigns = designs,
  respparamsets = wrap_param_set(rp, "resp_default"),
  blparamsets = bp_list,
  censorparams = NA,
  modelparams = mp,
  simparam = default_simparam(n_reps = 500),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = max(1, parallel::detectCores() - 1)
)

res$results$tau <- tau_grid[res$results$blparamset]
res$results$design_name <- c("Nof1", "RCT")[res$results$trialdesign]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
