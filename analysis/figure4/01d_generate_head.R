# 01d_generate_head.R
#
# Runs the Figure 4 parameter sweep using the current HEAD code
# loaded as a proper package via devtools::load_all().
#
# Usage:
#   NREPS=200 Rscript analysis/scripts/01d_generate_head.R

library(devtools)
load_all(".")

output_dir <- file.path("analysis", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv("NREPS", unset = "200"))
cat("=== Commit: HEAD (package load) ===\n")
cat("Nreps:", Nreps, "\n\n")

# ---------------------------------------------------------------
# Trial designs
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
# Parameters
# ---------------------------------------------------------------

data(extracted_bp)
blparamsets <- list(list(name = "TR",
  verbaldesc = "Extracted from data", param = extracted_bp))

data(extracted_rp)
respparamsets <- list(list(name = "TR",
  verbaldesc = "Extracted from data", param = extracted_rp))

censorparams <- data.table(
  names = c("balanced", "more of flat",
            "more of biased", "high dropout"),
  beta0 = c(.05, .05, .05, .15),
  beta1 = c(.5, .2, .8, .5),
  eb1 = c(2, 2, 2, 2)
)

coremodelparams <- expand.grid(
  N = c(35, 70),
  c.bm = c(0, .25, .45),
  carryover_t1half = c(0, .5, 1.0),
  c.tv = .7, c.pb = .7, c.br = .7,
  c.cf1t = .2, c.cfct = .1
)

analysisparams <- expand.grid(
  useDE = FALSE,
  t_random_slope = FALSE,
  full_model_out = FALSE
)

# ---------------------------------------------------------------
# Run
# ---------------------------------------------------------------

n_cores <- as.integer(Sys.getenv("NCORES", unset = "-1"))
if(n_cores < 0) n_cores <- max(1, parallel::detectCores() - 1)
cat("Starting simulation at", format(Sys.time()), "\n")
cat("Cores:", n_cores, "\n\n")

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
    basesavename = "HEAD",
    nRep2save = 5,
    saveunit2start = 1,
    savedir = output_dir
  ),
  analysisparams = analysisparams,
  rawdataout = FALSE,
  n_cores = n_cores
)

# ---------------------------------------------------------------
# Save
# ---------------------------------------------------------------

outfile <- file.path(output_dir,
  paste0("results_HEAD_n", Nreps, ".rds"))
saveRDS(simresults, outfile)

cat("\nComplete at", format(Sys.time()), "\n")
cat("Saved:", outfile, "\n")
cat("Dimensions:", dim(simresults$results), "\n")
