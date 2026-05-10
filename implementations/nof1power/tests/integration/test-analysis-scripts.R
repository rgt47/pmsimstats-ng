# Integration Test: Analysis Scripts
# Tests that key analysis scripts exist and are accessible

library(tinytest)
library(here)

# Main simulation script exists
script_path <- here("analysis", "scripts", "sim.R")
expect_true(file.exists(script_path))

# Monte Carlo design comparison vignette exists
script_path <- here("analysis", "scripts", "vig5.R")
expect_true(file.exists(script_path))

# Carryover analysis scripts exist
scripts <- c(
  here("analysis", "scripts", "vig_carry.R"),
  here("analysis", "scripts", "vig_carryover.R")
)

for (script in scripts) {
  expect_true(file.exists(script), info = paste("Missing:", script))
}

# Analysis outputs directories exist
figures_dir <- here("analysis", "figures")
expect_true(dir.exists(figures_dir))

test_plot_path <- file.path(figures_dir, "test_plot.png")
png(test_plot_path, width = 800, height = 600)
plot(1:10, 1:10, main = "Test Plot")
dev.off()

expect_true(file.exists(test_plot_path))

if (file.exists(test_plot_path)) {
  unlink(test_plot_path)
}

tables_dir <- here("analysis", "tables")
expect_true(dir.exists(tables_dir))
