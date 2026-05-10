#' S11. Carryover and AR(1) misspecification factorial.
#'
#' Precedent: Gaertner 2023; Blackston 2019 Section 7.3.
#' Factorial across DGP carryover on/off and DGP AR(1) high/off,
#' with analysis model always fitting corCAR1 and never-carryover
#' Dbc. Quantifies the Type I error inflation reported in
#' Blackston 2019 when moderate-to-strong carryover is ignored.

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

set.seed(42 + 11)

SWEEP_ID <- 11L
SWEEP_NAME <- "misspec"

grid <- expand.grid(
  carryover = c("off", "on"),
  ar1 = c("off", "high"),
  stringsAsFactors = FALSE
)

designs <- list(design_nof1(), design_rct())
rp <- baseline_resp_param()
bp <- baseline_bl_param()

build_model_params <- function() {
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    t_half <- if (grid$carryover[i] == "on") 3 / 7 else 0
    rho <- if (grid$ar1[i] == "high") 0.7 else 0
    mp <- baseline_model_param(carryover_t1half = t_half)
    mp$c.tv <- rho
    mp$c.pb <- rho
    mp$c.br <- rho
    mp
  }))
}

message("S11: running ", nrow(grid) * length(designs),
        " misspecification cells (null + non-null)...")

results_all <- list()

for (delta_label in c("null", "alt")) {
  rp_label <- rp
  if (delta_label == "null") {
    rp_label$max[rp_label$cat == "br"] <- 0
  }
  model_params <- build_model_params()

  res <- generate_simulated_results(
    trialdesigns = designs,
    respparamsets = wrap_param_set(
      rp_label, sprintf("resp_%s", delta_label)
    ),
    blparamsets = wrap_param_set(bp, "bl_default"),
    censorparams = NA,
    modelparams = model_params,
    simparam = default_simparam(n_reps = 500),
    analysisparams = default_analysisparams(),
    rawdataout = FALSE,
    dgp_architecture = "mvn",
    n_cores = max(1, parallel::detectCores() - 1)
  )

  res$results$delta_label <- delta_label
  res$results$carryover <- grid$carryover[res$results$modelparamset]
  res$results$ar1 <- grid$ar1[res$results$modelparamset]
  res$results$design_name <- c("Nof1", "RCT")[res$results$trialdesign]

  results_all[[delta_label]] <- res$results
}

combined <- do.call(rbind, results_all)

save_sweep(
  list(results = combined,
       parameterselections = list(grid = grid)),
  SWEEP_ID, SWEEP_NAME
)
