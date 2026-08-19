## analysis/scripts/carryover-sensitivity/16-run-sensitivity-s11.R
##
## Tier 2 sensitivity block S11: Type I error at the null (c_bm = 0)
## for the full G1-G5 set plus G9 (E7cr2), at the Hybrid, Architecture
## B reference cell (N = 70, exponential DGP, t1/2 = 1.0).
##
## This closes a gap flagged during S7/S8/S9/S10 development: every
## power comparison for G4 (E7, AIC-selected half-life) and G9
## (E7+CR2) was run at c_bm = 0.45, and neither specification's Type
## I error had been checked anywhere in this manuscript. G4 carries a
## specific, unaddressed risk that the other specifications do not:
## it selects the half-life by AIC and then tests the interaction
## coefficient from the selected model, a post-selection inference
## problem (the same data drive both the selection and the test) that
## can inflate Type I error beyond nominal. This block tests that
## directly rather than leaving it as a documented but unverified
## concern.
##
## Six specifications (G1, G2, G3, G4, G5, G9); G6/G7/G8 are not
## included here since their Type I behavior is expected, on
## structural grounds, to track G1/G2/G3's respectively (CR2 on an
## unselected model has no post-selection inference issue to begin
## with), and Block S6 already establishes CR2 restores nominal size
## for the model-based test's conservatism at this same cell for G3.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/16-run-sensitivity-s11.R [--dev] [--reps N]
##
## --dev    : 50 reps (development); default is 2000 (production, for
##            tighter MCSE on a rate expected near 0.05).
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
          if (dev_mode) 50 else 2000

seed <- 20260415L
set.seed(seed)

ref_null <- tibble::tibble(
  carryover_form = 'exponential', weibull_shape = 1.0, dgp_arch = 'mvn',
  t1half = 1.0, design = 'Hybrid', N = 70, c_bm = 0, block = 'S11'
)

specs_s11 <- c('E1', 'E3', 'E2', 'E7', 'E9', 'E7cr2')

cat(sprintf('S11: null cell (c_bm = 0), n_reps = %d, specs = %s\n',
            n_reps, paste(specs_s11, collapse = ', ')))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()
results <- simulate_cell_s7(ref_null, n_reps, specs = specs_s11)
t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

summary_s11 <- results |>
  group_by(spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    type1 = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_type1 = sqrt(type1 * (1 - type1) / n_converged),
    mean_best_t_half = mean(best_t_half, na.rm = TRUE),
    .groups = 'drop'
  )
print(as.data.frame(summary_s11))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '16-run-sensitivity-s11.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_s11,
  reference = as.list(ref_null),
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

out_file <- if (dev_mode) '16-sensitivity-s11-dev.rds' else '16-sensitivity-s11.rds'
saveRDS(list(results = results, summary = summary_s11,
             ref = ref_null, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows)',
        file.path(out_dir, out_file), nrow(results)))
