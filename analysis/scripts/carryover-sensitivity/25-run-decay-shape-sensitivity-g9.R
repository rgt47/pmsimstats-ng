## analysis/scripts/carryover-sensitivity/25-run-decay-shape-sensitivity-g9.R
##
## Extends the decay-shape sensitivity grid behind Figure 4
## (fig-hendrickson-c, Tier 1 grid decay-shape axis) from the
## original three specifications (G1-G3, 23-run-decay-shape-
## sensitivity.R) to the full nine-specification G1-G9 set, via
## simulate_cell_s7() (same pipeline as 19-run-tier1-hendrickson-g9.R,
## which performed the analogous extension for Panels A/B).
##
## Grid: carryover_form in {exponential, weibull(k=0.25),
## weibull(k=0.5), weibull(k=2.0), weibull(k=4.0)} x design in
## {CO, Hybrid, OLBDC} x N in {35, 70}, at c_bm = 0.45, t1half = 1.0
## (matching Figure 4's existing reference cell and grid). Same seed
## and grid construction as 23-run-decay-shape-sensitivity.R so the
## G1-G3 rows reproduce that script's results under common random
## numbers; only the specs argument differs.
##
## Progressive save: each of the 30 cells is checkpointed to its own
## file under a checkpoints/ subdirectory as soon as it completes,
## keyed by (mode, n_reps) so a --dev and a production run never
## collide. A rerun with the same mode/n_reps skips any cell whose
## checkpoint already exists, so a run killed partway through (by an
## interrupt, a session timeout, or a machine sleep) resumes from
## where it left off rather than restarting the full grid. Delete
## the checkpoint directory (or pass --restart) to force a clean
## re-run.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/25-run-decay-shape-sensitivity-g9.R [--dev] [--reps N] [--restart]
##
## --dev     : 3 reps/cell (smoke test); default is 500 (production).
## --reps N  : override the rep count for the current mode.
## --restart : ignore and overwrite any existing checkpoints for this
##             mode/n_reps combination.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
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
restart <- '--restart' %in% args
reps_idx <- which(args == '--reps')
n_reps_override <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else NA_integer_
n_reps <- if (!is.na(n_reps_override)) n_reps_override else
          if (dev_mode) 3 else 500

seed <- 20260415L
set.seed(seed)

specs_g9 <- c('E1', 'E3', 'E2', 'E7', 'E9',
             'E1cr2', 'E3cr2', 'E2cr2', 'E7cr2')

forms <- tidyr::expand_grid(
  design = c('CO', 'Hybrid', 'OLBDC'),
  N      = c(35, 70)
) |>
  dplyr::mutate(c_bm = 0.45, t1half = 1.0, dgp_arch = 'mvn')

grid <- dplyr::bind_rows(
  forms |> dplyr::mutate(carryover_form = 'exponential', weibull_shape = 1.0),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 0.25),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 0.5),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 2.0),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 4.0)
)

cat(sprintf('Decay-shape sensitivity (k=0.25, 0.5, 2.0, 4.0; G1-G9): %d cells, n_reps = %d\n',
            nrow(grid), n_reps))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

out_dir <- file.path(repo_root, 'analysis/data')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mode_tag <- if (dev_mode) 'dev' else 'prod'
ckpt_dir <- file.path(out_dir,
  sprintf('checkpoints-decay-shape-g9-%s-n%d', mode_tag, n_reps))
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

if (restart) {
  unlink(list.files(ckpt_dir, full.names = TRUE))
  cat(sprintf('--restart: cleared existing checkpoints in %s\n', ckpt_dir))
}

t_start <- Sys.time()

for (i in seq_len(nrow(grid))) {
  ckpt_file <- file.path(ckpt_dir, sprintf('cell-%02d.rds', i))
  if (file.exists(ckpt_file)) {
    cat(sprintf('[%2d/%d] skipping (checkpoint exists): %s\n',
                i, nrow(grid), ckpt_file))
    next
  }
  cell <- grid[i, ]
  cell_t0 <- Sys.time()
  cell_result <- simulate_cell_s7(cell, n_reps, specs = specs_g9)
  cell_result <- dplyr::bind_cols(cell[rep(1, nrow(cell_result)), ], cell_result)
  saveRDS(cell_result, ckpt_file)
  cat(sprintf('[%2d/%d] wrote %s (%.1f sec)\n',
              i, nrow(grid), ckpt_file,
              as.numeric(Sys.time() - cell_t0, units = 'secs')))
}

results <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
  readRDS(file.path(ckpt_dir, sprintf('cell-%02d.rds', i)))
})

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

meta <- list(
  script = '25-run-decay-shape-sensitivity-g9.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_g9,
  weibull_shapes = c(0.25, 0.5, 2.0, 4.0),
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

summary_tab <- results |>
  group_by(design, N, carryover_form, weibull_shape, spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    .groups = 'drop'
  )

print(as.data.frame(summary_tab) |>
       dplyr::select(design, N, carryover_form, weibull_shape, spec, power) |>
       dplyr::arrange(design, N, carryover_form, weibull_shape, spec))

out_file <- if (dev_mode) '02-decay-shape-sensitivity-g9-dev.rds' else '02-decay-shape-sensitivity-g9.rds'
saveRDS(list(results = results, summary = summary_tab, grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
