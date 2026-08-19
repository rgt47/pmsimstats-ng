## analysis/scripts/carryover-sensitivity/15-run-sensitivity-s10.R
##
## Tier 2 sensitivity block S10: does a CR2 cluster-robust standard
## error change the ranking among the five base specifications
## (G1-G5), at the Hybrid, Architecture B reference cell (N = 70,
## c_bm = 0.45, t1/2 = 1.0, exponential DGP)?
##
## Compares eight specifications (a ninth, G8 = E2+CR2, already
## exists from Block S6's independent run and is merged in at
## report-writing time rather than rerun here; a tenth, G5+CR2, does
## not apply since E9's paired-difference regression has no
## repeated-measures cluster structure left for CR2 to correct):
##
##   G1 = E1      Unadjusted
##   G2 = E3      Lag-adjusted
##   G3 = E2      Exposure-weighted
##   G4 = E7      AIC-selected half-life
##   G5 = E9      Paired-difference
##   G6 = E1cr2   Unadjusted + CR2
##   G7 = E3cr2   Lag-adjusted + CR2
##   G9 = E7cr2   AIC-selected half-life + CR2
##
## See simulation-core.R, "Biomarker-interacted carryover terms and
## their replacements (S7)", for the E1cr2/E3cr2/E7cr2 fitters
## (cr2_extract() helper), and spec-labels.R for the G-code display
## mapping (g_code, g_labels).
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/15-run-sensitivity-s10.R [--dev] [--reps N]
##
## --dev    : 50 reps (development); default is 500 (production).
## --reps N : override the rep count for the current mode.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(purrr)
  library(furrr)
})

repo_root <- here::here()
source(file.path(repo_root,
  'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
dev_mode <- '--dev' %in% args
reps_idx <- which(args == '--reps')
n_reps_override <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else NA_integer_
n_reps <- if (!is.na(n_reps_override)) n_reps_override else
          if (dev_mode) 50 else 500

seed <- 20260415L
set.seed(seed)

ref <- tibble::tibble(
  carryover_form = 'exponential', weibull_shape = 1.0, dgp_arch = 'mvn',
  t1half = 1.0, design = 'Hybrid', N = 70, c_bm = 0.45, block = 'S10'
)

specs_g <- c('E1', 'E2', 'E3', 'E7', 'E9', 'E1cr2', 'E3cr2', 'E7cr2')

cat(sprintf('S10 (G1-G9 minus G8): 1 cell, n_reps = %d, specs = %s\n',
            n_reps, paste(specs_g, collapse = ', ')))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()
results <- simulate_cell_s7(ref, n_reps, specs = specs_g)
t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

summary_g <- results |>
  group_by(spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    mean_best_t_half = mean(best_t_half, na.rm = TRUE),
    .groups = 'drop'
  )
print(as.data.frame(summary_g))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '15-run-sensitivity-s10.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_g,
  reference = as.list(ref),
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

out_file <- if (dev_mode) '15-sensitivity-s10-g-dev.rds' else '15-sensitivity-s10-g.rds'
saveRDS(list(results = results, summary = summary_g, ref = ref, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows)',
        file.path(out_dir, out_file), nrow(results)))
