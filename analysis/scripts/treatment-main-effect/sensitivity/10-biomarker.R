#' S10. Biomarker interaction and DGP-architecture sweep.
#'
#' Precedent: Hendrickson 2020; Section 9 gap of the lit review.
#' Varies c.bm (biomarker-treatment interaction strength) across
#' four levels and dgp_architecture across the two supported modes.
#' Three designs are compared: parallel RCT, aggregated N-of-1,
#' hybrid OL+BDC+crossover. Crossed with two carryover settings
#' (0 and 3 days). Highest-cost sweep; run last.

source(here::here(
  "analysis", "scripts", "sensitivity", "00-common.R"
))

set.seed(42 + 10)

SWEEP_ID <- 10L
SWEEP_NAME <- "biomarker"

c_bm_grid <- c(0, 0.2, 0.4, 0.6)
halflife_grid_days <- c(0, 3)
arch_grid <- c("mvn", "mean_moderation")

designs <- list(design_rct(),
                design_nof1_alternating(),
                design_hybrid())
rp <- baseline_resp_param()
bp <- baseline_bl_param()

results_all <- list()

for (arch in arch_grid) {
  message("S10: dgp_architecture = ", arch)

  model_params <- do.call(rbind, lapply(c_bm_grid, function(cb) {
    do.call(rbind, lapply(halflife_grid_days, function(hd) {
      baseline_model_param(
        c.bm = cb, carryover_t1half = hd / 7
      )
    }))
  }))

  res <- generate_simulated_results(
    trialdesigns = designs,
    respparamsets = wrap_param_set(rp, "resp_default"),
    blparamsets = wrap_param_set(bp, "bl_default"),
    censorparams = NA,
    modelparams = model_params,
    simparam = default_simparam(n_reps = 500),
    analysisparams = default_analysisparams(),
    rawdataout = FALSE,
    dgp_architecture = arch,
    n_cores = max(1, parallel::detectCores() - 1)
  )

  res$results$dgp_architecture <- arch
  res$results$c_bm <- rep(c_bm_grid,
                           each = length(halflife_grid_days))[
    res$results$modelparamset
  ]
  res$results$halflife_days <- rep(halflife_grid_days,
                                    length(c_bm_grid))[
    res$results$modelparamset
  ]
  res$results$design_name <- c("RCT", "Nof1-alt", "Hybrid")[
    res$results$trialdesign
  ]

  results_all[[arch]] <- res$results
}

combined <- do.call(rbind, results_all)

save_sweep(
  list(results = combined,
       parameterselections = list(
         c_bm_grid = c_bm_grid,
         halflife_grid_days = halflife_grid_days,
         arch_grid = arch_grid
       )),
  SWEEP_ID, SWEEP_NAME
)
