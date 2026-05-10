#' run-primary-simulation.R
#'
#' Wrapper that sources `simulation.R` and persists the resulting
#' workspace as `analysis/data/derived_data/sim_workspace.RData`,
#' which `analysis/report/treatment-main-effect/report.Rmd` loads
#' at knit time. Run from the project root:
#'
#'   Rscript analysis/scripts/treatment-main-effect/run-primary-simulation.R
#'
#' Wall-clock: approximately 5-15 min on a single core depending on
#' system load. The script honours the `quick_test` flag inside
#' `simulation.R`; the production path keeps `quick_test <- FALSE`
#' (n_simulations = 1000 for the primary RCT-vs-N-of-1 comparison).

simulation_path <- file.path(
  "analysis", "scripts", "treatment-main-effect", "simulation.R"
)
if (!file.exists(simulation_path)) {
  stop(sprintf("simulation.R not found at %s; run from project root.",
               simulation_path))
}

t0 <- Sys.time()
source(simulation_path, local = FALSE)
elapsed <- difftime(Sys.time(), t0, units = "mins")
cat(sprintf("\nsimulation.R completed in %.2f min.\n",
            as.numeric(elapsed)))

out_dir <- file.path("analysis", "data", "derived_data")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_path <- file.path(out_dir, "sim_workspace.RData")

save.image(file = out_path)
cat(sprintf("Workspace saved to %s.\n", out_path))
