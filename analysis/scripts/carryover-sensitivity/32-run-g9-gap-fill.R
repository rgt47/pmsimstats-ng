## analysis/scripts/carryover-sensitivity/32-run-g9-gap-fill.R
##
## Gap-fill simulation for paper 02's G9 recommendation. Protocol:
## 32-gap-fill-protocol.md (same directory). Three blocks, each
## submittable as its own simrun job:
##
##   A  dropout (24 cells x REPS) and a corrected-code size
##      diagnostic (4 null cells x REPS_SIZE)
##   B  DGP scenarios: Weibull decay shape, true half-life outside
##      the candidate grid, mean moderation, autocorrelation
##      (38 cells x REPS)
##   C  feasible prognostic-confound ladder at N = 70 (14 cells x
##      REPS) and prog-PB at N in {140, 280} (4 cells x REPS_LARGE)
##
## Every cell fits G1-G9, E0 and the F-family to the same simulated
## trial each replicate via simulate_cell_combined() (now with
## dropout plumbing), under the corrected DGP constants and the
## stop-on-non-PD default. A pre-flight eigenvalue check over every
## cell runs before any simulation and aborts if a cell is
## infeasible, so the silent-repair artifact of docs/65 cannot
## recur here.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/32-run-g9-gap-fill.R \
##     [--block A|B|C|all] [--reps R] [--reps-size R] [--reps-large R] \
##     [--n-boot B] [--cores K] [--out-dir PATH] [--preflight] [--smoke]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
BLOCK      <- argv('--block', 'all')
REPS       <- as.integer(argv('--reps', '500'))
REPS_SIZE  <- as.integer(argv('--reps-size', '2000'))
REPS_LARGE <- as.integer(argv('--reps-large', '1000'))
N_BOOT     <- as.integer(argv('--n-boot', '1'))
CORES      <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
PREFLIGHT  <- '--preflight' %in% args
SMOKE      <- '--smoke' %in% args
OUT_DIR <- argv('--out-dir',
  if (SMOKE) 'analysis/data/carryover-sensitivity/g9-gap-fill-smoke'
  else 'analysis/data/carryover-sensitivity/g9-gap-fill')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

base_cell <- list(design = 'Hybrid', N = 70, t1half = 0.5,
                  carryover_form = 'exponential', weibull_shape = 1,
                  dgp_arch = 'mvn', c_bm = 0)

cell_row <- function(block, sub, id, reps, ...) {
  tibble(block = block, sub = sub, cell_id = id, reps = reps,
         cell = list(modifyList(base_cell, list(...))))
}

## ---------------- block A ----------------
grid_a <- bind_rows(
  pmap_dfr(expand_grid(rate = c(0.1, 0.2, 0.3), mech = c('MCAR', 'MAR'),
                       c_bm = c(0, 0.45), N = c(35, 70)),
           function(rate, mech, c_bm, N)
             cell_row('A', 'A1-dropout',
                      sprintf('drop-%s%.1f-cbm%.2f-N%d', mech, rate, c_bm, N),
                      REPS, dropout_rate = rate, dropout_mech = mech,
                      c_bm = c_bm, N = N)),
  pmap_dfr(expand_grid(t1half = c(0.5, 1), N = c(35, 70)),
           function(t1half, N)
             cell_row('A', 'A2-size',
                      sprintf('size-t%.1f-cbm0.00-N%d', t1half, N),
                      REPS_SIZE, t1half = t1half, c_bm = 0, N = N))
)

## ---------------- block B ----------------
grid_b <- bind_rows(
  pmap_dfr(expand_grid(k = c(0.5, 2), c_bm = c(0, 0.25, 0.45), N = c(35, 70)),
           function(k, c_bm, N)
             cell_row('B', 'B1-shape',
                      sprintf('shape-k%.1f-cbm%.2f-N%d', k, c_bm, N),
                      REPS, carryover_form = 'weibull', weibull_shape = k,
                      analysis_form = 'exponential', analysis_shape = 1,
                      c_bm = c_bm, N = N)),
  pmap_dfr(expand_grid(t1half = c(0.2, 2, 4), c_bm = c(0, 0.25), N = c(35, 70)),
           function(t1half, c_bm, N)
             cell_row('B', 'B2-halflife',
                      sprintf('hl-t%.1f-cbm%.2f-N%d', t1half, c_bm, N),
                      REPS, t1half = t1half, c_bm = c_bm, N = N)),
  pmap_dfr(expand_grid(c_bm = c(0, 0.25, 0.45), N = c(35, 70)),
           function(c_bm, N)
             cell_row('B', 'B3-meanmod',
                      sprintf('meanmod-cbm%.2f-N%d', c_bm, N),
                      REPS, dgp_arch = 'mean_moderation', c_bm = c_bm, N = N)),
  ## c_bm = 0.30, not 0.45: the coupling ceiling falls with rho and
  ## 0.45 is infeasible at rho = 0.9 (pre-flight, 2026-10-26).
  pmap_dfr(expand_grid(rho = c(0.5, 0.9), c_bm = c(0, 0.30), N = c(35, 70)),
           function(rho, c_bm, N)
             cell_row('B', 'B4-rho',
                      sprintf('rho%.1f-cbm%.2f-N%d', rho, c_bm, N),
                      REPS, rho = rho, c_bm = c_bm, N = N))
)

## ---------------- block C ----------------
ladder <- tribble(
  ~field, ~strength, ~c_bm_levels,
  'c_pb', 0.40, c(0, 0.30),
  'c_pb', 0.45, c(0, 0.30),
  'c_tv', 0.40, c(0, 0.30),
  'c_tv', 0.45, c(0, 0.30),
  'c_bl', 0.60, c(0, 0.30),
  'c_bl', 0.80, c(0, 0.30),
  'c_bl', 0.45, c(0, 0.45)
)
path_name <- c(c_pb = 'progpb', c_tv = 'progtvconst', c_bl = 'progbl')
grid_c <- bind_rows(
  pmap_dfr(ladder, function(field, strength, c_bm_levels)
    map_dfr(c_bm_levels, function(c_bm) {
      extra <- list(c_bm = c_bm, N = 70); extra[[field]] <- strength
      do.call(cell_row, c(list('C', 'C1-ladder',
        sprintf('%s-conf%.2f-cbm%.2f-N70', path_name[[field]], strength, c_bm),
        REPS), extra))
    })),
  pmap_dfr(expand_grid(N = c(140, 280), c_pb = c(0.3, 0.45)),
           function(N, c_pb)
             cell_row('C', 'C2-largeN',
                      sprintf('progpb-conf%.2f-cbm0.00-N%d', c_pb, N),
                      REPS_LARGE, c_pb = c_pb, c_bm = 0, N = N))
)

grid <- bind_rows(grid_a, grid_b, grid_c)
if (BLOCK != 'all') grid <- grid |> filter(block == BLOCK)
if (SMOKE) {
  grid <- grid |> group_by(block) |> slice(c(1, n())) |> ungroup() |>
    mutate(reps = 2L)
}
stopifnot(!anyDuplicated(grid$cell_id))

## ---------------- pre-flight: every cell must be PD on every path ----------------
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
  select(block, sub, cell_id, pf) |>
  unnest(pf) |>
  group_by(block, sub, cell_id) |>
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

cat(sprintf('G9 gap fill, block %s: %d cells, %d (cell, rep) jobs, n_boot = %d\n',
            BLOCK, nrow(grid), sum(grid$reps), N_BOOT))

jobs <- grid |> select(cell_id, reps) |>
  mutate(rep = map(reps, seq_len)) |> select(-reps) |> unnest(rep)

one_job <- function(i) {
  cell_id <- jobs$cell_id[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  cell <- grid$cell[[match(cell_id, grid$cell_id)]]
  set.seed(26102700L + r)
  res <- tryCatch(
    simulate_cell_combined(cell, n_reps = 1, n_boot = N_BOOT),
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
  left_join(grid |> select(block, sub, cell_id), by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
n_done <- length(unique(paste(res$cell_id, res$rep)))
cat(sprintf('\n%d of %d (cell, rep) jobs present in %.1f min (this run)\n',
            n_done, nrow(jobs), elapsed_min))

## ---------------- summary: rejection, kappa, nonconv; F1A0B pooled ----------------
summarize_cell <- function(d) {
  nominal <- d |>
    filter(!(family == 'F' & method == 'F1A0')) |>
    group_by(family, method) |>
    summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
              n = sum(!is.na(p_value)),
              nonconv = mean(!converged),
              mean_est = mean(estimate, na.rm = TRUE),
              mean_se = mean(std_error, na.rm = TRUE),
              emp_se = sd(estimate, na.rm = TRUE),
              .groups = 'drop') |>
    mutate(kappa = mean_se / emp_se)
  boot_pool <- d |> filter(family == 'F', method == 'F1A0') |>
    distinct(rep, boot_t_f1a0) |> pull(boot_t_f1a0)
  obs_f1a0 <- d |> filter(family == 'F', method == 'F1A0') |> pull(t)
  f1a0 <- d |> filter(family == 'F', method == 'F1A0') |>
    summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
              n = sum(!is.na(p_value)), nonconv = mean(!converged),
              mean_est = mean(estimate, na.rm = TRUE),
              mean_se = mean(std_error, na.rm = TRUE),
              emp_se = sd(estimate, na.rm = TRUE)) |>
    mutate(kappa = mean_se / emp_se)
  bind_rows(
    nominal,
    tibble(family = 'F', method = 'F1A0') |> bind_cols(f1a0),
    tibble(family = 'F', method = 'F1A0B',
           reject = ws_size(abs(obs_f1a0), boot_pool),
           n = sum(!is.na(obs_f1a0)))
  )
}

cell_meta <- grid |>
  mutate(N = map_dbl(cell, 'N'), c_bm = map_dbl(cell, 'c_bm'),
         t1half = map_dbl(cell, 't1half'),
         rho = map_dbl(cell, ~ default_val(.x[['rho']], 0.7)),
         dgp_arch = map_chr(cell, 'dgp_arch'),
         carryover_form = map_chr(cell, 'carryover_form'),
         weibull_shape = map_dbl(cell, 'weibull_shape'),
         dropout_rate = map_dbl(cell, ~ default_val(.x[['dropout_rate']], 0)),
         dropout_mech = map_chr(cell, ~ default_val(.x[['dropout_mech']], 'none')),
         c_pb = map_dbl(cell, ~ default_val(.x[['c_pb']], 0)),
         c_tv = map_dbl(cell, ~ default_val(.x[['c_tv']], 0)),
         c_bl = map_dbl(cell, ~ default_val(.x[['c_bl']], 0))) |>
  select(-cell, -reps)

summary_tab <- res |>
  group_by(cell_id) |>
  group_modify(~ summarize_cell(.x)) |>
  ungroup() |>
  left_join(cell_meta, by = 'cell_id') |>
  relocate(block, sub, .before = cell_id)

cat('\nRejection rate per (cell, method); kappa = mean SE / empirical SE:\n')
print(as.data.frame(summary_tab |>
  filter(method %in% c('E1', 'E2', 'E7', 'E2cr2', 'E7cr2', 'E0', 'F1', 'F1A0B')) |>
  select(block, cell_id, method, reject, n, nonconv, kappa) |>
  mutate(across(c(reject, nonconv, kappa), ~ round(.x, 3)))), max = 4000)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))

meta <- list(
  script = '32-run-g9-gap-fill.R', block = BLOCK, date = Sys.time(),
  n_reps = REPS, n_reps_size = REPS_SIZE, n_reps_large = REPS_LARGE,
  n_boot = N_BOOT, smoke = SMOKE, elapsed_min_this_run = elapsed_min,
  r_version = R.version.string,
  git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                     error = function(e) NA_character_)
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
cat(sprintf('\nWritten: %s\n', OUT_DIR))
