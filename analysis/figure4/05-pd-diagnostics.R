# 05-pd-diagnostics.R
#
# Measures the rate at which the covariance matrix fails the
# positive definiteness check and requires correction via
# corpcor::make.positive.definite() during simulation.
#
# Strategy: monkey-patch corpcor::is.positive.definite to
# record every call and its result, then run the simulation
# and report statistics. No commit-specific files are modified.
#
# Environment variables:
#   COMMIT_DIR  - path to directory with extracted R files
#   COMMIT_HASH - short hash for labeling output
#   NREPS       - number of replicates (default: 50)
#
# Usage:
#   COMMIT_DIR=analysis/scripts/commit_42ac030 \
#   COMMIT_HASH=42ac030 NREPS=50 \
#   Rscript analysis/figure4/05-pd-diagnostics.R

library(data.table)
library(lme4)
library(lmerTest)
library(corpcor)
library(MASS)
library(tictoc)

commit_dir <- Sys.getenv("COMMIT_DIR", unset = "")
commit_hash <- Sys.getenv("COMMIT_HASH", unset = "")
Nreps <- as.integer(Sys.getenv("NREPS", unset = "50"))

if (!nzchar(commit_dir) || !nzchar(commit_hash)) {
  stop("COMMIT_DIR and COMMIT_HASH must be set")
}

cat("=== PD Diagnostics: Commit", commit_hash, "===\n")
cat("Source dir:", commit_dir, "\n")
cat("Nreps:", Nreps, "\n\n")

# ---------------------------------------------------------------
# Instrumentation: patch is.positive.definite and
# make.positive.definite to record calls
# ---------------------------------------------------------------

pd_log <- new.env(parent = emptyenv())
pd_log$checks <- 0L
pd_log$failures <- 0L
pd_log$fixes <- 0L
pd_log$eigenvalues_at_failure <- list()
pd_log$condition_numbers <- numeric()

original_is_pd <- corpcor::is.positive.definite
original_make_pd <- corpcor::make.positive.definite

instrumented_is_pd <- function(m, tol = 1e-8, method = c("eigen", "chol")) {
  pd_log$checks <- pd_log$checks + 1L
  result <- original_is_pd(m, tol = tol)

  evals <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  kappa <- max(evals) / max(min(evals), .Machine$double.eps)
  pd_log$condition_numbers <- c(pd_log$condition_numbers, kappa)

  if (!result) {
    pd_log$failures <- pd_log$failures + 1L
    pd_log$eigenvalues_at_failure[[
      length(pd_log$eigenvalues_at_failure) + 1]] <- evals
  }
  result
}

instrumented_make_pd <- function(m, tol = 1e-3) {
  pd_log$fixes <- pd_log$fixes + 1L
  original_make_pd(m, tol = tol)
}

# Replace in both the corpcor namespace and the global environment
# (sourced code resolves unqualified calls in the global env)
environment(instrumented_is_pd) <- asNamespace("corpcor")
environment(instrumented_make_pd) <- asNamespace("corpcor")
assignInNamespace("is.positive.definite", instrumented_is_pd,
                  ns = "corpcor")
assignInNamespace("make.positive.definite", instrumented_make_pd,
                  ns = "corpcor")
assign("is.positive.definite", instrumented_is_pd,
       envir = .GlobalEnv)
assign("make.positive.definite", instrumented_make_pd,
       envir = .GlobalEnv)

# ---------------------------------------------------------------
# Source commit-specific code
# ---------------------------------------------------------------

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
# Run simulation
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
    basesavename = paste0("pd_diag_", commit_hash),
    nRep2save = 5,
    saveunit2start = 1,
    savedir = output_dir
  ),
  analysisparams = analysisparams,
  rawdataout = FALSE
)

# ---------------------------------------------------------------
# Report PD diagnostics
# ---------------------------------------------------------------

cat("\n")
cat(strrep("=", 65), "\n")
cat("POSITIVE DEFINITENESS DIAGNOSTICS: Commit", commit_hash, "\n")
cat(strrep("=", 65), "\n\n")

cat("Total is.positive.definite() checks:", pd_log$checks, "\n")
cat("Failures (non-PD matrices):         ",
    pd_log$failures, "\n")
cat("make.positive.definite() calls:     ", pd_log$fixes, "\n")
cat("Failure rate:                        ",
    sprintf("%.4f (%d / %d)\n",
            pd_log$failures / max(pd_log$checks, 1),
            pd_log$failures, pd_log$checks))

if (length(pd_log$condition_numbers) > 0) {
  kappas <- pd_log$condition_numbers
  cat("\nCondition number statistics (all matrices):\n")
  cat("  Min:    ", sprintf("%.1f", min(kappas)), "\n")
  cat("  Median: ", sprintf("%.1f", median(kappas)), "\n")
  cat("  Mean:   ", sprintf("%.1f", mean(kappas)), "\n")
  cat("  Max:    ", sprintf("%.1f", max(kappas)), "\n")
  cat("  >100:   ", sum(kappas > 100), "of",
      length(kappas), "\n")
  cat("  >1000:  ", sum(kappas > 1000), "of",
      length(kappas), "\n")
}

if (length(pd_log$eigenvalues_at_failure) > 0) {
  cat("\nEigenvalue analysis of non-PD matrices:\n")
  for (i in seq_along(pd_log$eigenvalues_at_failure)) {
    evals <- pd_log$eigenvalues_at_failure[[i]]
    neg_evals <- evals[evals < 0]
    cat(sprintf(
      "  Failure %d: %d negative eigenvalues, min = %.6f, ",
      i, length(neg_evals), min(evals)))
    cat(sprintf("max = %.2f, kappa = %.1f\n",
      max(evals), max(evals) / max(min(evals),
                                    .Machine$double.eps)))
  }

  if (length(pd_log$eigenvalues_at_failure) > 5) {
    all_min_evals <- sapply(pd_log$eigenvalues_at_failure, min)
    cat("\n  Summary of minimum eigenvalues across failures:\n")
    cat("    Min of min:   ", sprintf("%.6f", min(all_min_evals)), "\n")
    cat("    Median of min:", sprintf("%.6f", median(all_min_evals)), "\n")
    cat("    Max of min:   ", sprintf("%.6f", max(all_min_evals)), "\n")
  }
}

# ---------------------------------------------------------------
# Save diagnostics
# ---------------------------------------------------------------

diag_result <- list(
  commit = commit_hash,
  nreps = Nreps,
  checks = pd_log$checks,
  failures = pd_log$failures,
  fixes = pd_log$fixes,
  failure_rate = pd_log$failures / max(pd_log$checks, 1),
  condition_numbers = pd_log$condition_numbers,
  eigenvalues_at_failure = pd_log$eigenvalues_at_failure
)

outfile <- file.path(output_dir,
  paste0("pd_diagnostics_", commit_hash, "_n", Nreps, ".rds"))
saveRDS(diag_result, outfile)

cat("\nDiagnostics saved to:", outfile, "\n")
cat("Simulation complete at", format(Sys.time()), "\n")

# Restore original functions
assignInNamespace("is.positive.definite", original_is_pd,
                  ns = "corpcor")
assignInNamespace("make.positive.definite", original_make_pd,
                  ns = "corpcor")
