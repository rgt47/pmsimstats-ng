## analysis/scripts/carryover-sensitivity/41-run-g10-backfill.R
##
## G10 backfill for paper 02 (split plan
## analysis/report/02-carryover-sensitivity/revisions/
## split-plan-prognostic-2026-10-29.md, Section 2.8). G10 (E7null_cr2:
## half-life chosen by ML AIC on the null model Sx ~ bm + t + Dbc,
## interaction model refitted by REML, bm:Dbc tested with CR2) has so
## far been fitted only in the referee-variant blocks N and A. This
## driver fits it in every other family paper 02 keeps, on the same
## trials as the stored G8 and G9 fits:
##
##   source  family  stored job and driver                      seed base
##   t1      1       tier1-g-vs-f-e0 (job 27, script 27)        24102600 + r
##   gvf     2       g-vs-f-comparison (job 21, script 06),     23102600 + r
##                   the four Stage 4 cells only
##   gf      2, 3, 4 g9-gap-fill (script 32): B1-shape and       26102700 + r
##                   B2-halflife (job 32), B4-rho (job 32),
##                   A1-dropout patient-level (job 33), and the
##                   A2-size diagnostic (job 31)
##   rvS     slope   referee-variants block S (job 37, script   label_seed()
##                   36)
##
## Each (cell, replicate) is regenerated with its stored job's seed by
## simulate_cell_g10() (simulation-core.R), which repeats the stored
## drivers' data-generation calls, and fitted with G9 (E7cr2) and G10
## (E7null_cr2). G9 is the pairing guard: its estimate must reproduce
## the stored G9 estimate to GUARD_TOL in every replicate. A source whose
## guard fails cannot be paired with its stored G8 and G9; per the
## plan, G8, G9 and G10 must then be rerun together for it.
##
## Not included: B3-meanmod and job 31's per-occasion A1-dropout cells
## (neither is in a paper 02 family), and the prognostic cells, which
## move to paper 15.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/41-run-g10-backfill.R \
##     [--source t1|gvf|gf|rvS|all] [--sub SUB] [--max-reps R] \
##     [--cores K] [--out-dir PATH] [--smoke] [--time]
##   --sub restricts the gf source to one sub-block (A1-dropout,
##   A2-size, B1-shape, B2-halflife, B4-rho); --max-reps caps the
##   replicates per cell; --smoke runs two replicates of every cell
##   into a separate directory; --time fits a few replicates of one
##   cell per (source, N) on one worker and projects the host time.
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
SOURCE   <- argv('--source', 'all')
SUB      <- argv('--sub', 'all')
MAX_REPS <- as.integer(argv('--max-reps', NA))
CORES    <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
SMOKE    <- '--smoke' %in% args
TIME     <- '--time' %in% args
DATA_DIR <- 'analysis/data/carryover-sensitivity'
OUT_DIR <- argv('--out-dir',
  file.path(DATA_DIR, if (SMOKE) 'g10-backfill-smoke' else 'g10-backfill'))
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
## The smoke run (2026-10-31, local) reproduced every stored trial with
## refitted-G9 differences up to 4.5e-8: optimizer and BLAS rounding on
## identical data. A different trial differs by the order of a standard
## error (about 0.1), so 1e-6 separates the two by five orders.
GUARD_TOL <- 1e-6

## ---------------- cell grids, copied from the stored drivers ----------------
## Each grid reproduces its driver's cell list and cell ids exactly;
## a copying error shows up as a cell id with no stored match or as a
## guard failure.

## Script 27: Panel A and Panel B, 30 cells, 500 replicates.
t1_grid <- bind_rows(
  expand_grid(design = c('CO', 'Hybrid', 'OLBDC'), N = c(35, 70),
              c_bm = c(0, 0.30, 0.45)) |>
    mutate(t1half = 1.0, panel = 'A'),
  expand_grid(design = c('CO', 'Hybrid', 'OLBDC'), N = c(35, 70),
              t1half = c(0, 0.5)) |>
    mutate(c_bm = 0.45, panel = 'B')
) |>
  mutate(carryover_form = 'exponential', weibull_shape = 1.0,
         dgp_arch = 'mvn',
         cell_id = sprintf('panel%s-design%s-N%d-cbm%.2f-t%.2f',
                           panel, design, N, c_bm, t1half))
t1_cells <- tibble(
  source = 't1', family = '1', sub = NA_character_,
  cell_id = t1_grid$cell_id, reps = 500L,
  cell = map(seq_len(nrow(t1_grid)), ~ as.list(t1_grid[.x, ])))

## Script 06: the four Stage 4 cells paper 02 uses, 500 replicates.
gvf_base <- list(N = 70, t1half = 0.5, carryover_form = 'exponential',
                 weibull_shape = 1, dgp_arch = 'mvn', design = 'Hybrid')
gvf_def <- list(
  t1hifhi_null  = list(c_bm = 0, t1half = 2),
  t1hifhi_alt   = list(c_bm = 0.25, t1half = 2),
  shapemis_null = list(c_bm = 0, carryover_form = 'weibull',
                       weibull_shape = 2, analysis_form = 'exponential',
                       analysis_shape = 1),
  shapemis_alt  = list(c_bm = 0.25, carryover_form = 'weibull',
                       weibull_shape = 2, analysis_form = 'exponential',
                       analysis_shape = 1)
)
gvf_cells <- tibble(
  source = 'gvf', family = '2', sub = NA_character_,
  cell_id = names(gvf_def), reps = 500L,
  cell = map(gvf_def, ~ modifyList(gvf_base, .x)))

## Script 32: sub-blocks A1 (patient-level), A2, B1, B2 and B4.
gf_base <- list(design = 'Hybrid', N = 70, t1half = 0.5,
                carryover_form = 'exponential', weibull_shape = 1,
                dgp_arch = 'mvn', c_bm = 0)
gf_row <- function(sub, family, id, reps, ...) {
  tibble(source = 'gf', family = family, sub = sub, cell_id = id,
         reps = as.integer(reps),
         cell = list(modifyList(gf_base, list(...))))
}
gf_cells <- bind_rows(
  pmap_dfr(expand_grid(rate = c(0.1, 0.2, 0.3), mech = c('MCAR', 'MAR'),
                       c_bm = c(0, 0.45), N = c(35, 70)),
           function(rate, mech, c_bm, N)
             gf_row('A1-dropout', '4',
                    sprintf('drop-%s%.1f-cbm%.2f-N%d-pt', mech, rate, c_bm, N),
                    500, dropout_rate = rate, dropout_mech = mech,
                    dropout_per = 'patient', c_bm = c_bm, N = N)),
  pmap_dfr(expand_grid(t1half = c(0.5, 1), N = c(35, 70)),
           function(t1half, N)
             gf_row('A2-size', 'A2',
                    sprintf('size-t%.1f-cbm0.00-N%d', t1half, N),
                    2000, t1half = t1half, c_bm = 0, N = N)),
  pmap_dfr(expand_grid(k = c(0.5, 2), c_bm = c(0, 0.25, 0.45), N = c(35, 70)),
           function(k, c_bm, N)
             gf_row('B1-shape', '2',
                    sprintf('shape-k%.1f-cbm%.2f-N%d', k, c_bm, N),
                    500, carryover_form = 'weibull', weibull_shape = k,
                    analysis_form = 'exponential', analysis_shape = 1,
                    c_bm = c_bm, N = N)),
  pmap_dfr(expand_grid(t1half = c(0.2, 2, 4), c_bm = c(0, 0.25), N = c(35, 70)),
           function(t1half, c_bm, N)
             gf_row('B2-halflife', '2',
                    sprintf('hl-t%.1f-cbm%.2f-N%d', t1half, c_bm, N),
                    500, t1half = t1half, c_bm = c_bm, N = N)),
  pmap_dfr(expand_grid(rho = c(0.5, 0.9), c_bm = c(0, 0.30), N = c(35, 70)),
           function(rho, c_bm, N)
             gf_row('B4-rho', '3',
                    sprintf('rho%.1f-cbm%.2f-N%d', rho, c_bm, N),
                    500, rho = rho, c_bm = c_bm, N = N))
)

## Script 36, block S: patient-specific slope, 1,000 replicates.
rv_base <- list(design = 'Hybrid', N = 70, t1half = 0.5,
                carryover_form = 'exponential', weibull_shape = 1,
                dgp_arch = 'mvn', c_bm = 0, hl_analysis = c(0.25, 1))
rvs_cells <- map_dfr(c(0, 0.45), function(c_bm)
  tibble(source = 'rvS', family = 'slope', sub = NA_character_,
         cell_id = sprintf('slope0.3-cbm%.2f-N70', c_bm), reps = 1000L,
         cell = list(modifyList(rv_base, list(c_bm = c_bm, N = 70,
                                              sd_slope = 0.3)))))

grid <- bind_rows(t1_cells, gvf_cells, gf_cells, rvs_cells) |>
  mutate(N = map_dbl(cell, 'N'), c_bm = map_dbl(cell, 'c_bm'),
         design = map_chr(cell, 'design'),
         key = paste(source, cell_id, sep = '|'))
if (SOURCE != 'all') grid <- grid |> filter(source == SOURCE)
if (SUB != 'all') grid <- grid |> filter(sub == SUB)
stopifnot(nrow(grid) > 0, !anyDuplicated(grid$key))
if (!is.na(MAX_REPS)) grid <- grid |> mutate(reps = pmin(reps, MAX_REPS))
if (SMOKE) grid <- grid |> mutate(reps = 2L)

seed_for <- function(source, cell_id, r) {
  switch(source,
    t1  = 24102600L + r,
    gvf = 23102600L + r,
    gf  = 26102700L + r,
    rvS = label_seed(sprintf('%s-rep-%d', cell_id, r)))
}

## ---------------- stored G8 and G9 fits, for the guard and pairing ----------------
load_stored <- function() {
  t1 <- readRDS(file.path(DATA_DIR, 'tier1-g-vs-f-e0/all.rds')) |>
    filter(family == 'G') |> mutate(source = 't1')
  gvf <- readRDS(file.path(DATA_DIR, 'g-vs-f-comparison/all.rds')) |>
    filter(family == 'G') |> rename(cell_id = cell) |>
    mutate(source = 'gvf')
  gf <- map_dfr(c(31, 32, 33), function(j)
    readRDS(file.path(DATA_DIR, 'g9-gap-fill', paste0('job-', j),
                      'all.rds')) |>
      filter(family == 'G') |> mutate(source = 'gf', job = j))
  ## Job 31's A1 cells are the per-occasion run (no '-pt' suffix);
  ## keep its A2 cells only, and job 33's A1 and job 32's B cells.
  gf <- gf |> filter((job == 31 & sub == 'A2-size') | job %in% c(32, 33))
  rvs <- readRDS(file.path(DATA_DIR, 'referee-variants/block-S/all.rds')) |>
    mutate(source = 'rvS')
  bind_rows(t1, gvf, gf, rvs) |>
    filter(method %in% c('E2cr2', 'E7cr2')) |>
    transmute(source, cell_id, rep, method, estimate, std_error, p_value,
              converged)
}

## ---------------- --time: per-replicate cost and host projection ----------------
## The host rate is calibrated against job 37 (block S, script 36):
## 2,000 replicates of the E2cr2 + E7cr2 + F1 fit in 481 s on zeus's 8
## cores. The same fit and this driver's fit are timed here on one
## local worker, on the same block S cell, and their ratio scales the
## host rate.
if (TIME) {
  n_time <- 3L
  time_fit <- function(f) {
    t0 <- Sys.time(); f()
    as.numeric(Sys.time() - t0, units = 'secs') / n_time
  }
  s_cell <- rvs_cells$cell[[2]]
  set.seed(1L)
  ref_s <- time_fit(function()
    simulate_cell_referee(s_cell, n_reps = n_time,
                          specs = c('E2cr2', 'E7cr2', 'F1')))
  set.seed(1L)
  g10_s <- time_fit(function() simulate_cell_g10(s_cell, n_reps = n_time))
  host_ref <- 481 * 8 / 2000
  host_per_local <- host_ref / ref_s
  cat(sprintf(paste0('Calibration on %s: script-36 fit %.2f s, G9 + G10 ',
                     '%.2f s per replicate locally; host %.2f core-s per ',
                     'local s\n'),
              rvs_cells$cell_id[2], ref_s, g10_s, host_per_local))

  probe <- grid |> group_by(source, N, design) |> slice(1) |> ungroup()
  probe$secs <- map_dbl(seq_len(nrow(probe)), function(i) {
    set.seed(seed_for(probe$source[i], probe$cell_id[i], 1L))
    time_fit(function() simulate_cell_g10(probe$cell[[i]], n_reps = n_time))
  })
  load <- grid |> count(source, N, design, wt = reps, name = 'reps') |>
    left_join(probe |> select(source, N, design, secs),
              by = c('source', 'N', 'design')) |>
    mutate(local_core_h = reps * secs / 3600,
           host_core_h = local_core_h * host_per_local)
  print(as.data.frame(load |>
    mutate(across(c(secs, local_core_h, host_core_h), ~ round(.x, 2)))),
    row.names = FALSE)
  cat(sprintf(paste0('\n%d replicates. Projected: %.1f local core-hours; ',
                     '%.1f host core-hours, %.1f h wall on 8 cores\n'),
              sum(load$reps), sum(load$local_core_h),
              sum(load$host_core_h), sum(load$host_core_h) / 8))
  quit(save = 'no')
}

cat(sprintf('G10 backfill, source %s, sub %s: %d cells, %d (cell, rep) jobs\n',
            SOURCE, SUB, nrow(grid), sum(grid$reps)))

jobs <- grid |> select(source, cell_id, reps) |>
  mutate(rep = map(reps, seq_len)) |> select(-reps) |> unnest(rep)

one_job <- function(i) {
  src <- jobs$source[i]; cell_id <- jobs$cell_id[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s--%s-rep-%04d.rds', src, cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  cell <- grid$cell[[match(paste(src, cell_id, sep = '|'), grid$key)]]
  set.seed(seed_for(src, cell_id, r))
  res <- tryCatch(simulate_cell_g10(cell, n_reps = 1),
                  error = function(e) NULL)
  if (!is.null(res)) {
    res <- mutate(res, source = src, cell_id = cell_id, rep = r,
                  .before = 1)
    saveRDS(res, f)
  }
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job,
                   mc.cores = CORES, mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

files <- list.files(ck_dir, full.names = TRUE)
files <- files[sub('-rep-\\d+\\.rds$', '', basename(files)) %in%
                 paste(grid$source, grid$cell_id, sep = '--')]
res <- bind_rows(lapply(files, readRDS)) |>
  left_join(grid |> select(source, family, sub, cell_id, N, c_bm, design),
            by = c('source', 'cell_id'))
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
n_done <- nrow(distinct(res, source, cell_id, rep))
cat(sprintf('\n%d of %d (cell, rep) jobs present in %.1f min (this run)\n',
            n_done, nrow(jobs), elapsed_min))

## ---------------- pairing guard: refitted G9 against stored G9 ----------------
stored <- load_stored()
guard <- res |>
  filter(method == 'E7cr2') |>
  select(source, cell_id, rep, est_new = estimate) |>
  left_join(stored |> filter(method == 'E7cr2') |>
              select(source, cell_id, rep, est_stored = estimate),
            by = c('source', 'cell_id', 'rep')) |>
  mutate(
    abs_diff = abs(est_new - est_stored),
    match = (is.na(est_new) & is.na(est_stored)) |
      (!is.na(abs_diff) & abs_diff <= GUARD_TOL))
guard_tab <- guard |>
  left_join(grid |> select(source, cell_id, sub), by = c('source', 'cell_id')) |>
  group_by(source, sub) |>
  summarise(reps = n(), matched = sum(match),
            max_abs_diff = suppressWarnings(max(abs_diff, na.rm = TRUE)),
            .groups = 'drop') |>
  mutate(pass = matched == reps)
cat(sprintf('\nPairing guard (refitted G9 against stored G9, tolerance %g):\n',
            GUARD_TOL))
print(as.data.frame(guard_tab), row.names = FALSE)
saveRDS(guard, file.path(OUT_DIR, 'guard.rds'))

## ---------------- summary and paired contrasts against stored G8 and G9 ----------------
summary_tab <- res |>
  group_by(source, family, sub, cell_id, N, c_bm, design, method) |>
  summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
            n = sum(!is.na(p_value)),
            mcse = sqrt(reject * (1 - reject) / n),
            nonconv = mean(!converged),
            mean_se = mean(std_error, na.rm = TRUE),
            emp_se = sd(estimate, na.rm = TRUE),
            .groups = 'drop') |>
  mutate(kappa = mean_se / emp_se)

## Paired only where the guard passed for that replicate.
paired_tab <- res |>
  filter(method == 'E7null_cr2') |>
  semi_join(guard |> filter(match), by = c('source', 'cell_id', 'rep')) |>
  select(source, cell_id, rep, p_g10 = p_value) |>
  inner_join(stored |> select(source, cell_id, rep, method, p_value) |>
               pivot_wider(names_from = method, values_from = p_value) |>
               rename(p_g8 = E2cr2, p_g9 = E7cr2),
             by = c('source', 'cell_id', 'rep')) |>
  mutate(d_g8 = (p_g10 < 0.05) - (p_g8 < 0.05),
         d_g9 = (p_g10 < 0.05) - (p_g9 < 0.05)) |>
  group_by(source, cell_id) |>
  summarise(n = sum(!is.na(d_g8)),
            g10_minus_g8 = mean(d_g8, na.rm = TRUE),
            se_g8 = sd(d_g8, na.rm = TRUE) / sqrt(n),
            g10_minus_g9 = mean(d_g9, na.rm = TRUE),
            se_g9 = sd(d_g9, na.rm = TRUE) / sqrt(n),
            .groups = 'drop')

selection_tab <- res |>
  filter(converged) |>
  count(source, cell_id, method, h) |>
  group_by(source, cell_id, method) |>
  mutate(share = n / sum(n)) |>
  ungroup()

cat('\nG10 rejection per cell, with paired differences against stored G8 and G9:\n')
print(as.data.frame(summary_tab |>
  filter(method == 'E7null_cr2') |>
  left_join(paired_tab |> select(-n), by = c('source', 'cell_id')) |>
  select(source, sub, cell_id, reject, mcse, n, nonconv, kappa,
         g10_minus_g8, se_g8, g10_minus_g9, se_g9) |>
  mutate(across(where(is.double), ~ round(.x, 4)))), max = 4000)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))
saveRDS(paired_tab, file.path(OUT_DIR, 'paired.rds'))
saveRDS(selection_tab, file.path(OUT_DIR, 'selection.rds'))

meta <- list(
  script = '41-run-g10-backfill.R', source = SOURCE, sub = SUB,
  max_reps = MAX_REPS, smoke = SMOKE, date = Sys.time(),
  guard_tol = GUARD_TOL, guard_pass = all(guard_tab$pass),
  elapsed_min_this_run = elapsed_min, r_version = R.version.string,
  git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                     error = function(e) NA_character_)
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
cat(sprintf('\nWritten: %s\n', OUT_DIR))
