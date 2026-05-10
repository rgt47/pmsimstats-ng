#' S13. Cycle-period design sweep (interaction-test target).
#'
#' Tier 1 of the cycle-period sensitivity sweep specified in
#' docs/27-cycle-period-design-sweep-plan.md. Walks the
#' (k, T_on, T_off) cycle-structure grid, crossed with carryover
#' half-life, biomarker effect size, and DGP architecture.
#'
#' Held fixed: T_runin = 0, T_openlabel = 8, T_runout = 0,
#' density = 1 obs/wk, N = 70, analysis spec A2 (matched-decay
#' continuous Dbc). Target coefficient: bm:Dbc interaction.
#'
#' Total subject-time per participant is
#'   T_openlabel + k * (T_on + T_off)
#' Cells with total > 32 weeks are filtered out (feasibility cap
#' for chronic-condition N-of-1 trials).

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

suppressPackageStartupMessages({
  library(tidyr)
  library(purrr)
})

set.seed(42 + 13)

SWEEP_ID <- 13L
SWEEP_NAME <- "cycle-period"

#=====================================================================
# Mode
#=====================================================================

args <- commandArgs(trailingOnly = TRUE)
MODE <- if ("--smoke" %in% args) {
  "smoke"
} else if ("--dev" %in% args) {
  "dev"
} else {
  "production"
}

n_reps <- switch(MODE, smoke = 5, dev = 50, production = 500)

#=====================================================================
# Cycle-structure grid
#=====================================================================

if (MODE == "smoke") {
  k_levels    <- c(1, 2)
  T_on_levels <- c(4)
  T_off_levels <- c(4)
  t_half_levels <- c(0.5)
  c_bm_levels  <- c(0.45)
  arch_levels  <- c("mvn")
} else {
  k_levels      <- c(1, 2, 3, 4)
  T_on_levels   <- c(2, 4, 6)
  T_off_levels  <- c(2, 4, 6)
  t_half_levels <- c(0, 0.5, 1.0)
  c_bm_levels   <- c(0, 0.10, 0.20, 0.30, 0.45)
  arch_levels   <- c("mvn", "mean_moderation")
}

T_OPENLABEL <- 8
DURATION_CAP <- if (MODE == "smoke") 1000 else 32

design_grid <- expand.grid(
  k = k_levels,
  T_on = T_on_levels,
  T_off = T_off_levels,
  KEEP.OUT.ATTRS = FALSE
) |>
  dplyr::mutate(
    duration = T_OPENLABEL + k * (T_on + T_off)
  ) |>
  dplyr::filter(duration <= DURATION_CAP)

message(sprintf(
  "S13 mode=%s: %d (k, T_on, T_off) cells under %d-week cap, n_reps=%d",
  MODE, nrow(design_grid), DURATION_CAP, n_reps
))

#=====================================================================
# Build trial designs
#=====================================================================

designs <- lapply(seq_len(nrow(design_grid)), function(i) {
  row <- design_grid[i, ]
  build_cycle_design(
    T_runin = 0,
    T_openlabel = T_OPENLABEL,
    k = row$k,
    T_on = row$T_on,
    T_off = row$T_off,
    T_runout = 0,
    order = "on_first",
    density = 1
  )
})

#=====================================================================
# Model parameter grid
#=====================================================================

mp_grid <- expand.grid(
  t_half = t_half_levels,
  c_bm = c_bm_levels,
  KEEP.OUT.ATTRS = FALSE
)
mp <- do.call(rbind, lapply(seq_len(nrow(mp_grid)), function(i) {
  baseline_model_param(
    N = 70,
    c.bm = mp_grid$c_bm[i],
    carryover_t1half = mp_grid$t_half[i]
  )
}))
mp$mp_index <- seq_len(nrow(mp))

#=====================================================================
# Sweep loop over architectures
#=====================================================================

rp <- baseline_resp_param()
bp <- baseline_bl_param()

results_per_arch <- list()
for (arch in arch_levels) {
  message(sprintf("S13 arch=%s: running %d designs x %d mp rows = %d cells, %d reps each",
    arch, length(designs), nrow(mp),
    length(designs) * nrow(mp), n_reps))

  res_arch <- generate_simulated_results(
    trialdesigns   = designs,
    respparamsets  = wrap_param_set(rp, "resp_default"),
    blparamsets    = wrap_param_set(bp, "bl_default"),
    censorparams   = NA,
    modelparams    = mp,
    simparam       = default_simparam(n_reps = n_reps,
                                       basesavename = sprintf(
                                         "s13-%s", arch)),
    analysisparams = default_analysisparams(target_coef = "interaction"),
    rawdataout     = FALSE,
    dgp_architecture = arch,
    n_cores        = max(1, parallel::detectCores() - 1)
  )

  ## Annotate results with the cycle-design axes and mp axes.
  res_arch$results <- res_arch$results |>
    dplyr::mutate(
      arch = arch,
      k     = design_grid$k[trialdesign],
      T_on  = design_grid$T_on[trialdesign],
      T_off = design_grid$T_off[trialdesign],
      duration = design_grid$duration[trialdesign],
      t_half = mp_grid$t_half[modelparamset],
      c_bm   = mp_grid$c_bm[modelparamset]
    )

  results_per_arch[[arch]] <- res_arch$results
}

results_combined <- dplyr::bind_rows(results_per_arch)

#=====================================================================
# Save
#=====================================================================

out <- list(
  results = results_combined,
  design_grid = design_grid,
  mp_grid = mp_grid,
  meta = list(
    script = "13-cycle-period.R",
    mode = MODE,
    date = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    seed = 42 + SWEEP_ID,
    n_reps = n_reps,
    T_openlabel = T_OPENLABEL,
    duration_cap = DURATION_CAP,
    r_version = R.version.string,
    elapsed_secs = NA  # populated below
  )
)

tic_start <- Sys.time()
out$meta$elapsed_secs <- as.numeric(
  difftime(Sys.time(), tic_start, units = "secs"))

save_sweep(out, SWEEP_ID, SWEEP_NAME)
message(sprintf("S13 done: %d total replicate rows",
                nrow(results_combined)))
