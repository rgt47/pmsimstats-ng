# 03-generate-power-by-commit.R
#
# Parameterized simulation script that runs the Figure 4 parameter
# sweep using R source files extracted from a specific git commit.
#
# Environment variables:
#   COMMIT_DIR  - path to directory with extracted R files (required)
#   COMMIT_HASH - short hash for labeling output (required)
#   NREPS       - number of replicates per cell (default: 200)
#
# Usage:
#   COMMIT_DIR=analysis/scripts/commit_42ac030 \
#   COMMIT_HASH=42ac030 \
#   NREPS=200 \
#   Rscript analysis/figure4/03-generate-power-by-commit.R

library(data.table)
library(lme4)
library(lmerTest)
library(corpcor)
library(MASS)
library(tictoc)

commit_dir <- Sys.getenv("COMMIT_DIR", unset = "")
commit_hash <- Sys.getenv("COMMIT_HASH", unset = "")
Nreps <- as.integer(Sys.getenv("NREPS", unset = "200"))

if (!nzchar(commit_dir) || !nzchar(commit_hash)) {
  stop("COMMIT_DIR and COMMIT_HASH must be set")
}

cat("=== Commit:", commit_hash, "===\n")
cat("Source dir:", commit_dir, "\n")
cat("Nreps:", Nreps, "\n\n")

source(file.path(commit_dir, "utilities.R"))
source(file.path(commit_dir, "buildtrialdesign.R"))
source(file.path(commit_dir, "generateData.R"))
source(file.path(commit_dir, "censordata.R"))
source(file.path(commit_dir, "lme_analysis.R"))
source(file.path(commit_dir, "generateSimulatedResults.R"))

output_dir <- file.path("analysis", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

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

extracted_bp <- data.table(
  cat = c("BL", "bm"),
  m = c(83.06897, 124.32759),
  sd = c(18.48267, 15.36159)
)
blparamsets <- list(list(name = "TR",
  verbaldesc = "Extracted from data", param = extracted_bp))

extracted_rp <- data.table(
  cat = c("tv", "pb", "br"),
  max = c(6.50647, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.35, 0.35, 0.42),
  sd = c(10, 10, 8)
)
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
  c.bm = c(0, .3, .6),
  carryover_t1half = c(0, .1, .2),
  c.tv = .8, c.pb = .8, c.br = .8,
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

cat("Starting simulation at", format(Sys.time()), "\n\n")

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
    basesavename = commit_hash,
    nRep2save = 5,
    saveunit2start = 1,
    savedir = output_dir
  ),
  analysisparams = analysisparams,
  rawdataout = FALSE
)

# ---------------------------------------------------------------
# Save
# ---------------------------------------------------------------

outfile <- file.path(output_dir,
  paste0("results_", commit_hash, "_n", Nreps, ".rds"))
saveRDS(simresults, outfile)

cat("\nComplete at", format(Sys.time()), "\n")
cat("Saved:", outfile, "\n")
cat("Dimensions:", dim(simresults$results), "\n")
