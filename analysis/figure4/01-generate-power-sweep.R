# 01-generate-power-sweep.R
#
# Reproduces the results_core simulation data from
# Hendrickson et al. (2020), Figure 4.
#
# This script runs the full parameter sweep over 4 trial designs,
# 2 sample sizes, 3 biomarker correlations, 3 carryover half-lives,
# and 5 censoring conditions. Output is saved as an .rds file in
# analysis/output/ (not overwriting the bundled results_core.rda).
#
# Runtime: ~8-12 hours at Nreps=1000 on a modern machine.
#          ~30-60 minutes at Nreps=50 for a quick validation run.
#
# Usage:
#   Rscript analysis/figure4/01-generate-power-sweep.R
#   NREPS=50 Rscript analysis/figure4/01-generate-power-sweep.R

library(devtools)
load_all(".")

output_dir <- file.path("analysis", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv("NREPS", unset = "1000"))
cat("Running with Nreps =", Nreps, "\n")

# ---------------------------------------------------------------
# 1. Trial designs (identical to vignette 1)
# ---------------------------------------------------------------

tdOL <- buildtrialdesign(
  name_longform = "open label",
  name_shortform = "OL",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste("OL", 1:8, sep = ""),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)

tdOLBDC <- buildtrialdesign(
  name_longform = "open label+blinded discontinuation",
  name_shortform = "OL+BDC",
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptname = c("OL1", "OL2", "OL3", "OL4",
                 "BD1", "BD2", "BD3", "BD4"),
  expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

tdCO <- buildtrialdesign(
  name_longform = "traditional crossover",
  name_shortform = "CO",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = c(paste("COa", 1:4, sep = ""),
                 paste("COb", 1:4, sep = "")),
  expectancies = rep(.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)

tdNof1 <- buildtrialdesign(
  name_longform = "primary N-of-1 design",
  name_shortform = "N-of-1",
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptname = c("OL1", "OL2", "BD1", "BD2",
                 "BD3", "BD4", "COd", "COp"),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)

trialdesigns <- list(OL = tdOL, OLBDC = tdOLBDC,
                     CO = tdCO, Nof1 = tdNof1)

# ---------------------------------------------------------------
# 2. Baseline parameters (extracted from pilot data)
# ---------------------------------------------------------------

data(extracted_bp)
TRblparam <- list(
  name = "TR",
  verbaldesc = "Extracted from data with blank slate assumptions",
  param = extracted_bp
)
blparamsets <- list(TRblparam)

# ---------------------------------------------------------------
# 3. Response parameters (tabula rasa)
# ---------------------------------------------------------------

data(extracted_rp)
TRrespparam <- list(
  name = "TR",
  verbaldesc = "Extracted from data with blank slate assumptions",
  param = extracted_rp
)
respparamsets <- list(TRrespparam)

# ---------------------------------------------------------------
# 4. Censoring parameters
# ---------------------------------------------------------------

censorparams <- data.table(
  names = c("balanced", "more of flat",
            "more of biased", "high dropout"),
  beta0 = c(.05, .05, .05, .15),
  beta1 = c(.5, .2, .8, .5),
  eb1 = c(2, 2, 2, 2)
)

# ---------------------------------------------------------------
# 5. Model parameters (core grid for Figure 4)
# ---------------------------------------------------------------

coremodelparams <- expand.grid(
  N = c(35, 70),
  c.bm = c(0, .3, .6),
  carryover_t1half = c(0, .1, .2),
  c.tv = .8, c.pb = .8, c.br = .8,
  c.cf1t = .2, c.cfct = .1
)

cat("Model parameter grid: ", nrow(coremodelparams), "combinations\n")
cat("Trial designs: ", length(trialdesigns), "\n")
cat("Total cells: ",
    nrow(coremodelparams) * length(trialdesigns), "\n")
cat("Total simulation runs: ",
    nrow(coremodelparams) * length(trialdesigns) * Nreps, "\n\n")

# ---------------------------------------------------------------
# 6. Analysis parameters
# ---------------------------------------------------------------

analysisparams <- expand.grid(
  useDE = FALSE,
  t_random_slope = FALSE,
  full_model_out = FALSE
)

# ---------------------------------------------------------------
# 7. Run the simulation
# ---------------------------------------------------------------

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
save_basename <- paste0("results_core_", timestamp)

cat("Starting simulation at", format(Sys.time()), "\n")

use_progressive <- (Nreps > 50)

if (use_progressive) {
  cat("Using progressive save (basename:", save_basename, ")\n\n")

  generateSimulatedResults(
    trialdesigns = list(trialdesigns$OL, trialdesigns$OLBDC,
                        trialdesigns$CO, trialdesigns$Nof1),
    respparamsets = respparamsets,
    blparamsets = blparamsets,
    censorparams = censorparams,
    modelparams = coremodelparams,
    simparam = list(
      Nreps = Nreps,
      progressiveSave = TRUE,
      basesavename = save_basename,
      nRep2save = 5,
      saveunit2start = 1,
      savedir = output_dir
    ),
    analysisparams = analysisparams,
    rawdataout = FALSE
  )

  cat("Recombining progressive save chunks...\n")
  simresults <- reknitsimresults(save_basename, output_dir)

} else {
  cat("Using single-pass mode (no progressive save)\n\n")

  simresults <- generateSimulatedResults(
    trialdesigns = list(trialdesigns$OL, trialdesigns$OLBDC,
                        trialdesigns$CO, trialdesigns$Nof1),
    respparamsets = respparamsets,
    blparamsets = blparamsets,
    censorparams = censorparams,
    modelparams = coremodelparams,
    simparam = list(
      Nreps = Nreps,
      progressiveSave = FALSE,
      basesavename = save_basename,
      nRep2save = 5,
      saveunit2start = 1,
      savedir = output_dir
    ),
    analysisparams = analysisparams,
    rawdataout = FALSE
  )
}

# ---------------------------------------------------------------
# 8. Save final combined results
# ---------------------------------------------------------------

outfile <- file.path(output_dir,
                     paste0(save_basename, "_final.rds"))
saveRDS(simresults, outfile)

cat("\nSimulation complete at", format(Sys.time()), "\n")
cat("Results saved to:", outfile, "\n")
cat("Results dimensions:", dim(simresults$results), "\n")
