#' S5. Within-patient variance (sigma^2) sweep.
#'
#' Precedent: Senn 2018; Wang-Schork 2019.
#' Varies the bio-response (br) SD which encodes the within-patient
#' measurement error component in the Senn variance decomposition.

source(here::here(
  "analysis", "scripts", "nof1-design-sensitivity", "00-common.R"
))

set.seed(42 + 5)

SWEEP_ID <- 5L
SWEEP_NAME <- "sigma2"

sigma_grid <- c(0.5, 0.75, 1.0, 1.25, 1.5)

designs <- list(design_nof1(), design_rct())
rp_base <- baseline_resp_param()
bp <- baseline_bl_param()

rp_list <- lapply(sigma_grid, function(s) {
  rp <- rp_base
  rp$sd[rp$cat == "br"] <- s
  list(name = sprintf("sigma=%s", s), param = rp)
})

mp <- baseline_model_param()

message("S5: running ", length(sigma_grid) * length(designs),
        " cells across sigma^2 levels...")

res <- generate_simulated_results(
  trialdesigns = designs,
  respparamsets = rp_list,
  blparamsets = wrap_param_set(bp, "bl_default"),
  censorparams = NA,
  modelparams = mp,
  simparam = default_simparam(n_reps = 500),
  analysisparams = default_analysisparams(),
  rawdataout = FALSE,
  dgp_architecture = "mvn",
  n_cores = max(1, parallel::detectCores() - 1)
)

res$results$sigma <- sigma_grid[res$results$respparamset]
res$results$design_name <- c("Nof1", "RCT")[res$results$trialdesign]

save_sweep(res, SWEEP_ID, SWEEP_NAME)
