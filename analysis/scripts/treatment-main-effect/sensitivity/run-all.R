#' run-all.R
#'
#' Orchestrator for the sensitivity plan
#' (docs/sensitivity-sweep-plan-2026-04-17.md).
#'
#' Execution order matches plan Section 5: smoke test first, then
#' S1, S12, S2, S6, S7, S4, S5, S8, S9, S3, S10, S11.
#'
#' Each sweep is sourced in its own R environment (via local()) so
#' that state does not leak between sweeps. Elapsed time is printed
#' per sweep.

script_dir <- file.path("analysis", "scripts", "sensitivity")

order <- c(
  "00-smoke-test.R",
  "01-effect-size.R",
  "12-morris-null.R",
  "02-carryover-halflife.R",
  "06-cycles.R",
  "07-patients.R",
  "04-tau2.R",
  "05-sigma2.R",
  "08-ar1.R",
  "09-period.R",
  "03-decay-form.R",
  "10-biomarker.R",
  "11-misspec.R"
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
