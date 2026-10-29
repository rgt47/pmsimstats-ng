#' run-all.R
#'
#' Orchestrator for the sensitivity sweeps (protocol:
#' 00-ademp-protocol.md).
#'
#' Run from the repository root with the corrected engine installed
#' (devtools::install_local('implementations/nof1power')). Execution
#' order: smoke test, the parallel-group check (14), the feasibility
#' table (17, no simulation), then S1, S12,
#' S12b (15), S2, S6, S7, S4, S5, S8, S9, S3, S10, S11, and last the
#' table of cell summaries (16). The commands, the expected host time
#' and the check list are in analysis/report/05-nof1-design-
#' sensitivity/README.md.
#'
#' Each sweep is sourced in its own R environment (via local()) so
#' that state does not leak between sweeps. Elapsed time is printed
#' per sweep.

script_dir <- file.path("analysis", "scripts", "nof1-design-sensitivity")

order <- c(
  "00-smoke-test.R",
  "14-rct-null-check.R",
  "17-feasibility.R",
  "01-effect-size.R",
  "12-morris-null.R",
  "15-xo-corcar1.R",
  "02-carryover-halflife.R",
  "06-cycles.R",
  "07-patients.R",
  "04-tau2.R",
  "05-sigma2.R",
  "08-ar1.R",
  "09-period.R",
  "03-decay-form.R",
  "10-biomarker.R",
  "11-misspec.R",
  "16-summarize-all.R"
)

total_start <- Sys.time()

for (script in order) {
  path <- file.path(script_dir, script)
  if (!file.exists(path)) {
    warning("Missing script: ", path)
    next
  }
  cat(sprintf(
    "\n=== %s ===\nStart: %s\n",
    script, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ))
  t0 <- Sys.time()
  local({
    source(path, local = TRUE)
  })
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  cat(sprintf("Done: %.2f min\n", as.numeric(elapsed)))
}

total_elapsed <- difftime(Sys.time(), total_start, units = "mins")
cat(sprintf(
  "\n== All sensitivity sweeps complete: %.2f min total ==\n",
  as.numeric(total_elapsed)
))
