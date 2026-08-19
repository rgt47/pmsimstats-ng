## analysis/scripts/carryover-sensitivity/18-run-sensitivity-s6-g9.R
##
## Extends Block S6 (originally 06-cr2-recalibration.R, G1-G3 model-
## based vs CR2 only) to the full nine-specification G1-G9 set, via
## simulate_cell_s7()/fit_spec_s7() instead of that script's bespoke
## fit_spec_both(). Same five stress cells (Reference, Small N, High
## autocorrelation, Half-life mis-specification, 30% MCAR dropout),
## each run at both c_bm = 0.45 (power) and c_bm = 0 (Type I error),
## so both the S6 power/gap comparison and the null-cell Type I
## comparison come from one unified run with common random numbers
## across all nine specifications within each cell (unlike the
## original 06 script and Block S10's G8, which came from
## independently-seeded runs).
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/18-run-sensitivity-s6-g9.R [--dev] [--reps N]
##
## --dev    : 50 reps/cell (development); default is 500 (production).
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

specs_g9 <- c('E1', 'E3', 'E2', 'E7', 'E9',
             'E1cr2', 'E3cr2', 'E2cr2', 'E7cr2')

base_cell <- tibble::tibble(
  design = 'Hybrid', N = 70, c_bm = 0.45, t1half = 1.0,
  carryover_form = 'exponential', weibull_shape = 1, rho = 0.7,
  dgp_arch = 'mvn', analysis_t1half = NA_real_,
  dropout_rate = NA_real_, dropout_mech = NA_character_
)

mk_cell <- function(label, ...) {
  d <- base_cell
  overrides <- list(...)
  for (nm in names(overrides)) d[[nm]] <- overrides[[nm]]
  d$cell <- label
  d
}

stress_cells <- dplyr::bind_rows(
  mk_cell('Reference (N=70)'),
  mk_cell('Small N (N=35)', N = 35),
  mk_cell('High autocorr (rho=0.9)', rho = 0.9),
  mk_cell('Half-life mis-spec (assume 0.25)', analysis_t1half = 0.25),
  mk_cell('30% MCAR dropout', dropout_rate = 0.30, dropout_mech = 'MCAR')
)

## Power cells (c_bm = 0.45) and null cells (c_bm = 0) mirroring each
## stress cell, so the null rejection rate is the Type I error.
grid <- dplyr::bind_rows(
  stress_cells |> dplyr::mutate(arm = 'power'),
  stress_cells |> dplyr::mutate(c_bm = 0, arm = 'null')
)

cat(sprintf('S6 (G1-G9): %d cells (5 stress x power/null), n_reps = %d\n',
            nrow(grid), n_reps))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()

results <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
  cell <- grid[i, ]
  cell_result <- simulate_cell_s7(cell, n_reps, specs = specs_g9)
  dplyr::bind_cols(cell[rep(1, nrow(cell_result)), ], cell_result)
})

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '18-run-sensitivity-s6-g9.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_g9,
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

summary_g9 <- results |>
  group_by(cell, arm, spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    .groups = 'drop'
  )

print(as.data.frame(summary_g9) |> dplyr::filter(arm == 'power') |>
       dplyr::select(cell, spec, power) |> dplyr::arrange(cell, spec))

out_file <- if (dev_mode) '18-sensitivity-s6-g9-dev.rds' else '18-sensitivity-s6-g9.rds'
saveRDS(list(results = results, summary = summary_g9, grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
