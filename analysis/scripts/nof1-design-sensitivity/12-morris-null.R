#' S12. Morris null-effect calibration.
#'
#' Precedent: Morris 2019; nof1_power report Section 4.3.
#' Higher-replicate Type I error check under delta = 0 across four
#' representative (tau, sigma, halflife) settings for all four
#' designs. Target MCSE ~0.5% at nominal 0.05.

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

set.seed(42 + 12)

SWEEP_ID <- 12L
SWEEP_NAME <- "morris-null"

configs <- list(
  list(tau = 1.5,  sigma = 1.0, halflife_days = 0),
  list(tau = 1.5,  sigma = 1.0, halflife_days = 3),
  list(tau = 0.75, sigma = 1.5, halflife_days = 3),
  list(tau = 2.0,  sigma = 0.5, halflife_days = 3)
)

designs <- list(
  design_rct(), design_crossover(),
  design_nof1_alternating(), design_hybrid()
)
design_names <- c("RCT", "XO", "Nof1-alt", "Hybrid")

rp_base <- baseline_resp_param()
rp_base$max[rp_base$cat == "br"] <- 0  # main-effect null

results_all <- list()

for (i in seq_along(configs)) {
  cfg <- configs[[i]]
  message(sprintf(
    "S12 config %d/%d: tau=%g sigma=%g halflife_days=%g",
    i, length(configs), cfg$tau, cfg$sigma, cfg$halflife_days
  ))

  bp <- baseline_bl_param()
  bp$sd[bp$cat == "BL"] <- cfg$tau

  rp <- rp_base
  rp$sd[rp$cat == "br"] <- cfg$sigma

  mp <- baseline_model_param(
    carryover_t1half = cfg$halflife_days / 7
  )

  res <- generate_simulated_results(
    trialdesigns = designs,
    respparamsets = wrap_param_set(rp, "resp_brmax0"),
    blparamsets = wrap_param_set(bp, "bl_config"),
    censorparams = NA,
    modelparams = mp,
    simparam = default_simparam(n_reps = 2000),
    analysisparams = default_analysisparams(),
    rawdataout = FALSE,
    dgp_architecture = "mvn",
    n_cores = max(1, parallel::detectCores() - 1)
  )

  res$results$config <- i
  res$results$tau <- cfg$tau
  res$results$sigma <- cfg$sigma
  res$results$halflife_days <- cfg$halflife_days
  res$results$design_name <- design_names[res$results$trialdesign]

  results_all[[i]] <- res$results
}

combined <- do.call(rbind, results_all)

save_sweep(
  list(results = combined,
       parameterselections = list(configs = configs)),
  SWEEP_ID, SWEEP_NAME
)
