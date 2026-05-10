# Integration Test: Data Pipeline
# Tests that research data files exist and are accessible

library(tinytest)
library(here)

# Raw data files exist
raw_data_dir <- here("analysis", "data", "raw_data")
expect_true(dir.exists(raw_data_dir))

ct_data_path <- file.path(raw_data_dir, "CTdata.rda")
expect_true(file.exists(ct_data_path))

# Derived data files exist
derived_data_dir <- here("analysis", "data", "derived_data")
expect_true(dir.exists(derived_data_dir))

expected_files <- c(
  "extracted_bp.rda",
  "extracted_rp.rda",
  "results_core.rda",
  "results_maxes.rda",
  "results_rates.rda",
  "results_trajectories.rda"
)

for (file in expected_files) {
  file_path <- file.path(derived_data_dir, file)
  expect_true(file.exists(file_path), info = paste("Missing:", file))
}

# Simulation results cache exists
derived_data_dir <- here("analysis", "data", "derived_data")

cache_file <- file.path(derived_data_dir, "simulation_results_20250827.RData")
expect_true(file.exists(cache_file))

summary_file <- file.path(derived_data_dir, "simulation_summary_20250827.rds")
expect_true(file.exists(summary_file))

# CSV summaries are available
derived_data_dir <- here("analysis", "data", "derived_data")

expect_true(file.exists(file.path(derived_data_dir, "quick_sim_raw.csv")))
expect_true(file.exists(file.path(derived_data_dir, "quick_sim_summary.csv")))
expect_true(file.exists(file.path(derived_data_dir, "biomarker_pvalues_table.csv")))
