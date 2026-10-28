## analysis/scripts/carryover-sensitivity/36-run-referee-variants.R
##
## Simulation variants requested by the paper 02 referee report of
## 2026-10-27 (analysis/report/02-carryover-sensitivity/
## referee-report-2026-10-27.md), modeled on 32-run-g9-gap-fill.R.
## Four blocks, each submittable as its own simrun job:
##
##   N  null cells carrying the recommendation (M1, M4): Hybrid,
##      N in {35, 70}, true t1/2 in {0.5, 1}, c_bm = 0, with and
##      without a prognostic PB coupling of 0.3 (8 cells x REPS_NULL,
##      default 5,000)
##   A  alternative cells (M9): c_bm in {0.25, 0.45}, Hybrid,
##      N in {35, 70}, true t1/2 in {0.5, 2} (8 cells x REPS_ALT,
##      default 1,000)
##   T  prog-TV-trend pathway (M7): b_tv_trend = 0.3 with c_bm in
##      {0, 0.45}, Hybrid, N = 70, t1/2 = 0.5 (2 cells x REPS_OTHER,
##      default 1,000)
##   S  patient-specific slope (M12): sd_slope = 0.3 outcome units per
##      week with c_bm in {0, 0.45}, Hybrid, N = 70, t1/2 = 0.5
##      (2 cells x REPS_OTHER)
##
## Blocks N and A fit, to the same simulated trial each replicate,
## G8 (E2cr2), G9 (E7cr2), the null-model selection variant
## (E7null_cr2, M1), F1, and G8 at the fixed analysis half-lives 0.25
## and 1 week whatever the truth (E3fix_cr2, M9). Blocks T and S fit
## G8, G9 and F1. The half-life each method used is stored per
## replicate in `h` (selected for G9 and G9-null, fixed otherwise),
## so the selection distribution of M9 can be tabulated.
## simulate_cell_referee() in simulation-core.R does the fitting.
##
## Seeds (M4): every (cell, replicate) job seeds from a hash of its
## own label, label_seed(sprintf('%s-rep-%d', cell_id, r)), so cells
## do not share a random-number stream. Consequently no cell here is
## paired with any stored run or with any other cell; comparisons
## across cells are between independent estimates, comparisons across
## methods within a cell remain paired.
##
## Cost (timed locally 2026-10-27, one worker, Apple Silicon): 1.08 s
## per replicate for the five-specification fit at N = 70 (the
## heaviest cell, null-t1.0-N70-cpb0.3), 0.49 s at N = 35, and 0.56 s
## for the three-specification fit of blocks T and S. 52,000
## replicates in all, about 40,000 local worker-seconds (11 h on one
## worker). The same machine runs simulate_cell_combined() at 3.7 s
## per replicate against the 6 to 7 core-seconds measured on the
## simulation host (32-gap-fill-protocol.md), so the host is taken as
## 1.8 times slower per thread: about 20 host core-hours. On the
## host's 2 physical cores (4 threads) that is roughly 8 h with 4
## threads if hyperthreading yields about 2.5 cores' worth of
## throughput, 10 h if it yields none: block N about 6.5 h (8 cells
## x 5,000), A about 1.3 h, T and S about 15 min each. Estimates,
## not measurements: no job has run this driver on the host.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/36-run-referee-variants.R \
##     [--block N|A|T|S|all] [--reps-null R] [--reps-alt R] \
##     [--reps-other R] [--cores K] [--out-dir PATH] \
##     [--preflight] [--smoke] [--time]
##   --preflight checks every cell's covariance matrix and exits;
##   --smoke runs two replicates of every cell; --time runs 5
##   replicates of the heaviest cell on one worker and reports the
##   per-replicate wall time.
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
BLOCK      <- argv('--block', 'all')
REPS_NULL  <- as.integer(argv('--reps-null', '5000'))
REPS_ALT   <- as.integer(argv('--reps-alt', '1000'))
REPS_OTHER <- as.integer(argv('--reps-other', '1000'))
CORES      <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
PREFLIGHT  <- '--preflight' %in% args
SMOKE      <- '--smoke' %in% args
TIME       <- '--time' %in% args
OUT_DIR <- argv('--out-dir',
  if (SMOKE) 'analysis/data/carryover-sensitivity/referee-variants-smoke'
  else 'analysis/data/carryover-sensitivity/referee-variants')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

HL_FIXED <- c(0.25, 1)
SPECS_FULL  <- c('E2cr2', 'E7cr2', 'E7null_cr2', 'F1', 'E3fix_cr2')
SPECS_SHORT <- c('E2cr2', 'E7cr2', 'F1')

base_cell <- list(design = 'Hybrid', N = 70, t1half = 0.5,
                  carryover_form = 'exponential', weibull_shape = 1,
                  dgp_arch = 'mvn', c_bm = 0, hl_analysis = HL_FIXED)

cell_row <- function(block, id, reps, specs, ...) {
  tibble(block = block, cell_id = id, reps = reps, specs = list(specs),
         cell = list(modifyList(base_cell, list(...))))
}

grid_n <- pmap_dfr(
  expand_grid(t1half = c(0.5, 1), N = c(35, 70), c_pb = c(0, 0.3)),
  function(t1half, N, c_pb)
    cell_row('N', sprintf('null-t%.1f-N%d-cpb%.1f', t1half, N, c_pb),
             REPS_NULL, SPECS_FULL, t1half = t1half, N = N, c_bm = 0,
             c_pb = c_pb))

grid_a <- pmap_dfr(
  expand_grid(t1half = c(0.5, 2), N = c(35, 70), c_bm = c(0.25, 0.45)),
  function(t1half, N, c_bm)
    cell_row('A', sprintf('alt-t%.1f-N%d-cbm%.2f', t1half, N, c_bm),
             REPS_ALT, SPECS_FULL, t1half = t1half, N = N, c_bm = c_bm))

grid_t <- map_dfr(c(0, 0.45), function(c_bm)
  cell_row('T', sprintf('trend0.3-cbm%.2f-N70', c_bm), REPS_OTHER,
           SPECS_SHORT, c_bm = c_bm, N = 70, b_tv_trend = 0.3))

grid_s <- map_dfr(c(0, 0.45), function(c_bm)
  cell_row('S', sprintf('slope0.3-cbm%.2f-N70', c_bm), REPS_OTHER,
           SPECS_SHORT, c_bm = c_bm, N = 70, sd_slope = 0.3))

grid <- bind_rows(grid_n, grid_a, grid_t, grid_s)
if (BLOCK != 'all') grid <- grid |> filter(block == BLOCK)
stopifnot(nrow(grid) > 0, !anyDuplicated(grid$cell_id))
if (SMOKE) grid <- grid |> mutate(reps = 2L)

## ---------------- pre-flight: every cell must be PD on every path ----------------
## The trend shift and the patient slope are additive after the MVN
## draw and do not enter the covariance matrix, so the check is the
## same as in 32-run-g9-gap-fill.R.
preflight_cell <- function(cell) {
  design_set <- build_design_set(cell$design)
  rho <- default_val(cell[['rho']], 0.7)
  mp <- list(N = cell$N, carryover_t1half = cell$t1half,
             carryover_form = cell$carryover_form,
             weibull_shape = cell$weibull_shape,
             dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
             c.tv = rho, c.pb = rho, c.br = rho, c.cf1t = 0.2, c.cfct = 0.2,
             c.bm.bl = default_val(cell[['c_bl']], 0),
             c.bm.tv = default_val(cell[['c_tv']], 0),
             b.bm.tv = default_val(cell[['b_tv']], 0),
             c.bm.pb = default_val(cell[['c_pb']], 0))
  map_dfr(seq_along(design_set), function(i) {
    b <- build_sigma_matrix(mp, default_resp_param(), default_baseline_param(),
                            design_set[[i]], stop_on_nonpd = FALSE)
    raw <- outer(b$standard_deviations, b$standard_deviations) * b$correlations
    pos <- b$standard_deviations > 0
    ev <- eigen(raw[pos, pos], symmetric = TRUE, only.values = TRUE)$values
    tibble(path = i, min_eigenvalue = min(ev), feasible = min(ev) > 0)
  })
}
pre <- grid |>
  mutate(pf = map(cell, preflight_cell)) |>
  select(block, cell_id, pf) |>
  unnest(pf) |>
  group_by(block, cell_id) |>
  summarise(min_eigenvalue = min(min_eigenvalue),
            feasible = all(feasible), .groups = 'drop')
saveRDS(pre, file.path(OUT_DIR, 'preflight.rds'))
cat(sprintf('Pre-flight: %d cells, %d feasible, min eigenvalue %.4f\n',
            nrow(pre), sum(pre$feasible), min(pre$min_eigenvalue)))
if (!all(pre$feasible)) {
  print(as.data.frame(pre |> filter(!feasible)), row.names = FALSE)
  stop('Infeasible cells present; adjust the grid before running.')
}
if (PREFLIGHT) {
  print(as.data.frame(pre |> mutate(min_eigenvalue = round(min_eigenvalue, 4))),
        row.names = FALSE)
  quit(save = 'no')
}

## ---------------- --time: per-replicate cost of the heaviest cell ----------------
if (TIME) {
  heavy <- grid_n |> filter(cell_id == 'null-t1.0-N70-cpb0.3')
  cell <- heavy$cell[[1]]
  set.seed(label_seed(heavy$cell_id))
  t0 <- Sys.time()
  res <- simulate_cell_referee(cell, n_reps = 5, specs = SPECS_FULL)
  secs <- as.numeric(Sys.time() - t0, units = 'secs')
  cat(sprintf('%s: 5 replicates, %d method rows, %.2f s per replicate on one worker\n',
              heavy$cell_id, nrow(res), secs / 5))
  quit(save = 'no')
}

cat(sprintf('Referee variants, block %s: %d cells, %d (cell, rep) jobs\n',
            BLOCK, nrow(grid), sum(grid$reps)))

jobs <- grid |> select(cell_id, reps) |>
  mutate(rep = map(reps, seq_len)) |> select(-reps) |> unnest(rep)

one_job <- function(i) {
  cell_id <- jobs$cell_id[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  k <- match(cell_id, grid$cell_id)
  cell <- grid$cell[[k]]
  set.seed(label_seed(sprintf('%s-rep-%d', cell_id, r)))
  res <- tryCatch(
    simulate_cell_referee(cell, n_reps = 1, specs = grid$specs[[k]]),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    res <- mutate(res, rep = r, cell_id = cell_id, .before = 1)
    saveRDS(res, f)
  }
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job,
                   mc.cores = CORES, mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

files <- list.files(ck_dir, full.names = TRUE)
files <- files[sub('-rep-\\d+\\.rds$', '', basename(files)) %in% grid$cell_id]
res <- bind_rows(lapply(files, readRDS)) |>
  left_join(grid |> select(block, cell_id), by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
n_done <- length(unique(paste(res$cell_id, res$rep)))
cat(sprintf('\n%d of %d (cell, rep) jobs present in %.1f min (this run)\n',
            n_done, nrow(jobs), elapsed_min))

## ---------------- summary ----------------
## Per (cell, method, fixed half-life): rejection rate with its MCSE,
## non-convergence, kappa. E3fix_cr2 is split by the fixed half-life
## it used (`h_fixed`); for the selecting methods h_fixed is NA and
## the selection distribution is tabulated separately.
cell_meta <- grid |>
  mutate(N = map_dbl(cell, 'N'), c_bm = map_dbl(cell, 'c_bm'),
         t1half = map_dbl(cell, 't1half'),
         c_pb = map_dbl(cell, ~ default_val(.x[['c_pb']], 0)),
         b_tv_trend = map_dbl(cell, ~ default_val(.x[['b_tv_trend']], 0)),
         sd_slope = map_dbl(cell, ~ default_val(.x[['sd_slope']], 0))) |>
  select(-cell, -reps, -specs)

summary_tab <- res |>
  mutate(h_fixed = if_else(method == 'E3fix_cr2', h, NA_real_)) |>
  group_by(cell_id, method, h_fixed) |>
  summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
            n = sum(!is.na(p_value)),
            mcse = sqrt(reject * (1 - reject) / n),
            nonconv = mean(!converged),
            mean_est = mean(estimate, na.rm = TRUE),
            mean_se = mean(std_error, na.rm = TRUE),
            emp_se = sd(estimate, na.rm = TRUE),
            .groups = 'drop') |>
  mutate(kappa = mean_se / emp_se) |>
  left_join(cell_meta, by = 'cell_id') |>
  relocate(block, .before = cell_id)

selection_tab <- res |>
  filter(method %in% c('E7cr2', 'E7null_cr2'), converged) |>
  count(cell_id, method, h) |>
  group_by(cell_id, method) |>
  mutate(share = n / sum(n)) |>
  ungroup() |>
  left_join(cell_meta |> select(cell_id, t1half), by = 'cell_id')

## Paired contrasts against G8 on the same trials (referee M4(c)):
## per-replicate rejection differences and their paired SE.
paired_tab <- res |>
  filter(method == 'E2cr2') |>
  select(cell_id, rep, p_g8 = p_value) |>
  inner_join(res |> filter(method != 'E2cr2') |>
               mutate(h_fixed = if_else(method == 'E3fix_cr2', h, NA_real_)) |>
               select(cell_id, rep, method, h_fixed, p_value),
             by = c('cell_id', 'rep')) |>
  mutate(d = (p_value < 0.05) - (p_g8 < 0.05)) |>
  group_by(cell_id, method, h_fixed) |>
  summarise(diff_vs_g8 = mean(d, na.rm = TRUE),
            paired_se = sd(d, na.rm = TRUE) / sqrt(sum(!is.na(d))),
            n = sum(!is.na(d)), .groups = 'drop')

cat('\nRejection rate per (cell, method, fixed half-life); kappa = mean SE / empirical SE:\n')
print(as.data.frame(summary_tab |>
  select(block, cell_id, method, h_fixed, reject, mcse, n, nonconv, kappa) |>
  mutate(across(c(reject, mcse, nonconv, kappa), ~ round(.x, 3)))), max = 4000)
cat('\nSelected half-life shares (G9 and G9-null):\n')
print(as.data.frame(selection_tab |> mutate(share = round(share, 3))), max = 2000)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))
saveRDS(selection_tab, file.path(OUT_DIR, 'selection.rds'))
saveRDS(paired_tab, file.path(OUT_DIR, 'paired.rds'))

meta <- list(
  script = '36-run-referee-variants.R', block = BLOCK, date = Sys.time(),
  reps_null = REPS_NULL, reps_alt = REPS_ALT, reps_other = REPS_OTHER,
  hl_fixed = HL_FIXED, smoke = SMOKE, elapsed_min_this_run = elapsed_min,
  r_version = R.version.string,
  git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                     error = function(e) NA_character_)
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
cat(sprintf('\nWritten: %s\n', OUT_DIR))
