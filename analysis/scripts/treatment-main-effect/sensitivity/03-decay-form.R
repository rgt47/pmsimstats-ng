#' S3. Carryover decay-form sweep with analysis misspecification.
#'
#' Precedent: Gaertner 2023; Section 9 gap of the lit review.
#' DGP carryover uses one of four forms (exponential, linear,
#' Weibull with k=0.8, Weibull with k=1.5, power), while the
#' analysis model always assumes exponential decay for Dbc.
#' This quantifies the bias and power cost of form mis-specification.

source(here::here(
  "analysis", "scripts", "sensitivity", "00-common.R"
))

set.seed(42 + 3)

SWEEP_ID <- 3L
SWEEP_NAME <- "decay-form"

form_grid <- list(
  list(form = "exponential", shape = 1),
  list(form = "linear",      shape = 1),
  list(form = "weibull",     shape = 0.8),
  list(form = "weibull",     shape = 1.5),
  list(form = "power",       shape = 1.5)
)

designs <- list(design_nof1_alternating(), design_rct())
rp <- baseline_resp_param()
bp <- baseline_bl_param()

model_params <- do.call(rbind, lapply(form_grid, function(fg) {
  mp <- baseline_model_param(carryover_t1half = 3 / 7)
  mp$carryover_form <- fg$form
  mp$weibull_shape <- fg$shape
  mp
}))

message("S3: running ", nrow(model_params) * length(designs),
        " cells (DGP decay forms)...")

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

form_labels <- vapply(form_grid, function(fg) {
  if (fg$form == "weibull") {
    sprintf("weibull_k=%g", fg$shape)
  } else {
    fg$form
  }
}, character(1))

res$results$decay_form <- form_labels[res$results$modelparamset]
res$results$design_name <- c("Nof1-alt", "RCT")[res$results$trialdesign]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
