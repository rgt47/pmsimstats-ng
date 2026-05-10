#' S1. Effect-size sweep.
#'
#' Precedent: Blackston et al. 2019 (DOI 10.3390/healthcare7040137);
#' replicated in the current prazosin-for-PTSD report at a single
#' effect size.
#'
#' Varies the true effect size Delta over a seven-point grid while
#' holding carryover, variance components, and design geometry fixed
#' at the baseline. Two designs are compared: aggregated N-of-1 and
#' parallel RCT at matched observation totals.

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

set.seed(42 + 1)

SWEEP_ID <- 1L
SWEEP_NAME <- "effect-size"

delta_grid <- c(-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0)

td_nof1 <- design_nof1()
td_rct <- design_rct()
rp <- baseline_resp_param()
bp <- baseline_bl_param()

model_params <- baseline_model_param(
  N = 40, c.bm = 0.3, carryover_t1half = 3 / 7
)

rp_by_delta <- lapply(delta_grid, function(d) {
  rp_d <- rp
  rp_d$max[rp_d$cat == "br"] <- abs(d)
  list(name = sprintf("delta=%s", d), param = rp_d)
})

message("S1: running ", length(rp_by_delta) * 2,
        " effect-size cells across 2 designs...")

res <- generate_simulated_results(
  trialdesigns = list(td_nof1, td_rct),
  respparamsets = rp_by_delta,
  blparamsets = wrap_param_set(bp, "bl_default"),
  censorparams = NA,
  modelparams = model_params,
  simparam = default_simparam(n_reps = 1000),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = max(1, parallel::detectCores() - 1)
)

res$results$delta <- delta_grid[res$results$respparamset]
res$results$design_name <- c("Nof1", "RCT")[res$results$trialdesign]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
