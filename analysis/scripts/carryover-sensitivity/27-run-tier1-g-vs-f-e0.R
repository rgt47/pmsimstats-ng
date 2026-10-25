## analysis/scripts/carryover-sensitivity/27-run-tier1-g-vs-f-e0.R
##
## Extends 19-run-tier1-hendrickson-g9.R's 30-cell Tier 1 grid
## (Figures 1-2, Panels A/B: exponential DGP, covariance moderation
## only, design x {CO, Hybrid, OLBDC}, N x {35, 70}) from G1-G9 alone
## to the full comparison this project has lacked a shared-grid
## version of (docs/62 Part 3, criterion 3): G1-G9, F1/F1A/F1A0/
## F1A0B, and E0 (Hendrickson's actual original analysis, fetched
## directly from GitHub commit 58b32a9; see report.Rmd Block S12),
## all fit to the *same* simulated trial each replicate, under the
## corrected c.cf1t/c.cfct = 0.2 DGP constants (docs/61, fixed
## 2026-10-21; simulation-core.R's model-param constructors already
## carry the fix, so no DGP-side change was needed here -- this
## script only widens which specifications are fit per cell).
##
## This replaces 19-run-tier1-hendrickson-g9.R's output for this grid
## (02-grid-summary-hendrickson-g9.rds) rather than supplementing it;
## that file's G1-G9-only numbers were run under the pre-fix DGP
## constants and are now superseded for this grid. Figures 1-2's
## G1-G9 panels should be regenerated from this script's output once
## it completes (12-hendrickson-heatmap.R reads the old file; update
## it to read this one once this run is verified).
##
## Scope decision (2026-10-24): the manuscript's full headline grid
## (540 cells / the 216-cell covariance-moderation slice) combined
## with G1-G9 + F + E0 was estimated at ~44 hours on zeus's 8 cores
## (11.7 core-sec/replicate, calibrated from simrun job 21's 15-cell,
## 11-spec run) -- not practical as one job. This 30-cell grid
## matches the scope already established for the G1-G9 Tier 1
## figures and is estimated at ~6.1 hours at the same per-replicate
## cost, rather than inventing new scope.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/27-run-tier1-g-vs-f-e0.R \
##     [--reps R] [--n-boot B] [--cores K] [--out-dir PATH]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
REPS <- as.integer(argv('--reps', '500'))
N_BOOT <- as.integer(argv('--n-boot', '1'))
CORES <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
OUT_DIR <- argv('--out-dir',
  'analysis/data/carryover-sensitivity/tier1-g-vs-f-e0')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

## Same grid as 19-run-tier1-hendrickson-g9.R: Panel A (c_bm x 3,
## t1half = 1.0) plus Panel B (t1half x 2, c_bm = 0.45), design x
## {CO, Hybrid, OLBDC}, N x {35, 70}, exponential DGP, mvn
## architecture. 30 unique cells.
panel_a <- tidyr::expand_grid(
  design = c('CO', 'Hybrid', 'OLBDC'),
  N      = c(35, 70),
  c_bm   = c(0, 0.30, 0.45)
) |>
  dplyr::mutate(t1half = 1.0, panel = 'A')

panel_b <- tidyr::expand_grid(
  design = c('CO', 'Hybrid', 'OLBDC'),
  N      = c(35, 70),
  t1half = c(0, 0.5)
) |>
  dplyr::mutate(c_bm = 0.45, panel = 'B')

grid <- dplyr::bind_rows(panel_a, panel_b) |>
  dplyr::mutate(carryover_form = 'exponential', weibull_shape = 1.0,
               dgp_arch = 'mvn') |>
  dplyr::mutate(cell_id = sprintf('panel%s-design%s-N%d-cbm%.2f-t%.2f',
    panel, design, N, c_bm, t1half))

cat(sprintf('Tier 1, G1-G9 + F1/F1A/F1A0/F1A0B + E0: %d cells, n_reps = %d\n',
            nrow(grid), REPS))

jobs <- tidyr::expand_grid(cell_i = seq_len(nrow(grid)), rep = seq_len(REPS))

one_job <- function(i) {
  cell_i <- jobs$cell_i[i]; r <- jobs$rep[i]
  cell_id <- grid$cell_id[cell_i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(24102600L + r)
  res <- tryCatch(
    simulate_cell_combined(as.list(grid[cell_i, ]), n_reps = 1,
                           n_boot = N_BOOT),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    res <- dplyr::mutate(res, rep = r, cell_id = cell_id, .before = 1)
    saveRDS(res, f)
  }
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job,
                   mc.cores = CORES, mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

res <- bind_rows(lapply(list.files(ck_dir, full.names = TRUE), readRDS))
res <- dplyr::left_join(res, grid, by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
n_done <- length(unique(paste(res$cell_id, res$rep)))
cat(sprintf(
  '\n%d of %d (cell, rep) jobs present in %.1f min (this run; checkpointed jobs from earlier runs are not re-timed)\n',
  n_done, nrow(jobs), elapsed_min))

## ---------------- summary: nominal rejection rate per method, plus warp-speed F1A0B ----------------
summarize_cell <- function(d) {
  nominal <- d |>
    filter(!(family == 'F' & method == 'F1A0')) |>
    group_by(family, method) |>
    summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
              n = sum(!is.na(p_value)), .groups = 'drop')
  boot_pool <- d |> filter(family == 'F', method == 'F1A0') |>
    distinct(rep, boot_t_f1a0) |> pull(boot_t_f1a0)
  obs_f1a0 <- d |> filter(family == 'F', method == 'F1A0') |> pull(t)
  f1a0_nominal <- d |> filter(family == 'F', method == 'F1A0') |>
    summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
              n = sum(!is.na(p_value)))
  bind_rows(
    nominal,
    tibble(family = 'F', method = 'F1A0', reject = f1a0_nominal$reject,
          n = f1a0_nominal$n),
    tibble(family = 'F', method = 'F1A0B',
          reject = ws_size(abs(obs_f1a0), boot_pool),
          n = sum(!is.na(obs_f1a0)))
  )
}

summary_tab <- res |>
  group_by(cell_id, panel, design, N, c_bm, t1half) |>
  group_modify(~ summarize_cell(.x)) |>
  ungroup()

cat('\nRejection rate (Type I error at c_bm=0 cells, power elsewhere):\n')
print(as.data.frame(summary_tab), max = 2000)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))

meta <- list(
  script = '27-run-tier1-g-vs-f-e0.R',
  date = Sys.time(),
  n_reps = REPS,
  n_boot = N_BOOT,
  elapsed_min_this_run = elapsed_min,
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
