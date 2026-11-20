## analysis/scripts/carryover-sensitivity/55-run-a4-baseline.R
##
## A4 with and without the baseline row, on the same trials.
##
## The published analysis of Hendrickson et al. prepends a baseline
## (BL, t = 0) outcome row per patient; none of paper 02's ladder
## analyses has it. Script 54 showed that the row carries much of the
## published analysis's power (E0pnc against E1). This driver asks
## whether the same row helps the recommended analysis A4 (E7null_cr2)
## and keeps its size. The baseline row is the only change: same
## formula, half-life grid, ML selection on the null model, REML refit
## and CR2 test, and every patient kept (the published exclusion rule
## is not used).
##
## Fits per trial, plain and with the baseline row ('bl' suffix):
##   A2  E1cr2       Sx ~ bm + t + Db + bm:Db, CR2
##   A3  E2cr2       Sx ~ bm + t + Dbc + bm:Dbc at the true half-life, CR2
##   A4  E7null_cr2  half-life by ML-AIC on Sx ~ bm + t + Dbc, REML refit
##                   with bm:Dbc, CR2
## The baseline row is Sx = BL, t = 0, Db = Dbc = L = tsd = De = 0, so
## Dbc rebuilt at any candidate half-life is 0 there.
##
## Trials: the 30 cells of job 27 (tier1-g-vs-f-e0), regenerated with
## its seeds (24102600 + r) exactly as script 54 does. Guard: the plain
## refits must match the stored E1cr2, E2cr2 and E7null_cr2 estimates
## (and the stored A4 half-life); replicates that fail are excluded.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/55-run-a4-baseline.R \
##     [--check] [--smoke] [--time] [--max-reps R] [--cores K] \
##     [--params published|simcore2026]
## --check fits one trial per design and verifies the baseline-row
## construction and that row order does not change the CR2 result.
.args0 <- commandArgs(trailingOnly = TRUE)
PARAMS <- if ('--params' %in% .args0)
  .args0[which(.args0 == '--params') + 1] else 'published'
stopifnot(PARAMS %in% c('simcore2026', 'published'))
Sys.setenv(PMSIM_PARAMS = PARAMS)
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')
stopifnot(param_set() == PARAMS)

args <- .args0
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
CHECK    <- '--check' %in% args
SMOKE    <- '--smoke' %in% args
TIME     <- '--time' %in% args
MAX_REPS <- as.integer(argv('--max-reps', NA))
CORES    <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
DATA_DIR <- 'analysis/data/carryover-sensitivity'
SET_DIR  <- if (PARAMS == 'published') file.path(DATA_DIR, 'published') else
  DATA_DIR
OUT_DIR  <- file.path(SET_DIR, if (SMOKE) 'a4-baseline-smoke' else
  'a4-baseline')
ck_dir <- file.path(OUT_DIR, 'trials')
GUARD_TOL <- 1e-4

## ---------------- the baseline row ----------------
add_bl_rows <- function(dat, dat_long) {
  pts <- dat_long |> distinct(ptID, path, bm)
  bl <- dat |> transmute(ptID = paste(path, ptID, sep = '_'), Sx = BL)
  bl_rows <- pts |> inner_join(bl, by = 'ptID') |>
    mutate(timepoint_name = 'BL', t = 0, Db = 0, Dbc = 0, L = 0,
           tsd = 0, De = 0)
  stopifnot(nrow(bl_rows) == nrow(pts), min(dat_long$t) > 0)
  ## Sorted by patient and time: cr2_extract() clusters by row position.
  bind_rows(bl_rows, dat_long) |> arrange(ptID, t)
}

METHODS <- c('E1cr2', 'E2cr2', 'E7null_cr2')
LABEL <- c(E1cr2 = 'A2', E2cr2 = 'A3', E7null_cr2 = 'A4')
fit_set <- function(dl, suffix) {
  map_dfr(METHODS, function(m) fit_spec_s7(dl, m) |>
            transmute(method = paste0(LABEL[m], suffix), estimate,
                      std_error, p_value, best_t_half, converged))
}
fit_all <- function(tr) {
  if (is.null(tr))
    return(tibble(method = c(paste0(LABEL, ''), paste0(LABEL, 'bl')),
                  estimate = NA_real_, std_error = NA_real_,
                  p_value = NA_real_, best_t_half = NA_real_,
                  converged = FALSE))
  bind_rows(fit_set(tr$dat_long, ''),
            fit_set(add_bl_rows(tr$dat, tr$dat_long), 'bl'))
}

## ---------------- regeneration, as script 54 ----------------
gen_trial <- function(cell) {
  design_set <- build_design_set(cell$design)
  rho <- default_val(cell[['rho']], 0.7)
  model_param <- list(
    N = cell$N, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
    c.tv = rho, c.pb = rho, c.br = rho, c.cf1t = 0.2, c.cfct = 0.2)
  model_param$c.bm.bl <- default_val(cell[['c_bl']], 0)
  model_param$c.bm.tv <- default_val(cell[['c_tv']], 0)
  model_param$b.bm.tv <- default_val(cell[['b_tv']], 0)
  model_param$b.bm.tv.trend <- default_val(cell[['b_tv_trend']], 0)
  model_param$sd_slope <- default_val(cell[['sd_slope']], 0)
  model_param$c.bm.pb <- default_val(cell[['c_pb']], 0)
  dat <- tryCatch(
    generate_data_multi_path(model_param, default_resp_param(),
                             default_baseline_param(), design_set,
                             stop_on_nonpd = TRUE),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat_long <- prepare_long_data(
    dat, design_set, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form, weibull_shape = cell$weibull_shape)
  list(dat = dat, dat_long = dat_long)
}

t1_grid <- bind_rows(
  expand_grid(design = c('CO', 'Hybrid', 'OLBDC'), N = c(35, 70),
              c_bm = c(0, 0.30, 0.45)) |> mutate(t1half = 1.0, panel = 'A'),
  expand_grid(design = c('CO', 'Hybrid', 'OLBDC'), N = c(35, 70),
              t1half = c(0, 0.5)) |> mutate(c_bm = 0.45, panel = 'B')
) |>
  mutate(carryover_form = 'exponential', weibull_shape = 1.0,
         dgp_arch = 'mvn',
         cell_id = sprintf('panel%s-design%s-N%d-cbm%.2f-t%.2f',
                           panel, design, N, c_bm, t1half))
seed_for <- function(r) 24102600L + r
cell_of <- function(i) as.list(t1_grid[i, ])

## ---------------- --check ----------------
if (CHECK) {
  ok <- TRUE
  for (d in c('CO', 'Hybrid', 'OLBDC')) {
    i <- which(t1_grid$cell_id ==
                 sprintf('panelA-design%s-N35-cbm0.45-t1.00', d))
    set.seed(seed_for(1L))
    tr <- gen_trial(cell_of(i))
    dl <- tr$dat_long
    da <- add_bl_rows(tr$dat, dl)
    n_pt <- n_distinct(dl$ptID)
    bl <- filter(da, timepoint_name == 'BL')
    for (hl in c(0.25, 0.5, 1, 2))
      stopifnot(!anyNA(exposure_weight(da$Db, da$tsd, hl)),
                all(exposure_weight(bl$Db, bl$tsd, hl) == 0))
    dup_t <- da |> count(ptID, t) |> filter(n > 1) |> nrow()
    f_sorted <- fit_spec_s7(da, 'E2cr2')
    f_unsorted <- fit_spec_s7(
      bind_rows(filter(da, timepoint_name == 'BL'),
                filter(da, timepoint_name != 'BL')), 'E2cr2')
    same <- isTRUE(all.equal(f_sorted$p_value, f_unsorted$p_value,
                             tolerance = 1e-8))
    bl_sx <- all(bl$Sx == tr$dat$BL[match(bl$ptID,
                   paste(tr$dat$path, tr$dat$ptID, sep = '_'))])
    cat(sprintf(paste0('%-6s rows %d -> %d (+%d patients); dup t %d; ',
                       'BL Sx matches %s; CR2 order-invariant %s ',
                       '(p %.6g vs %.6g)\n'),
                d, nrow(dl), nrow(da), n_pt, dup_t, bl_sx, same,
                f_sorted$p_value, f_unsorted$p_value))
    ok <- ok && nrow(da) == nrow(dl) + n_pt && dup_t == 0 && bl_sx
    print(as.data.frame(fit_all(tr)), digits = 4)
  }
  cat(if (ok) '\nCHECK PASSED\n' else '\nCHECK FAILED\n')
  quit(save = 'no', status = if (ok) 0 else 1)
}

## ---------------- --time ----------------
if (TIME) {
  probe <- t1_grid |> group_by(design, N) |> slice(1) |> ungroup()
  secs <- map_dbl(seq_len(nrow(probe)), function(i) {
    set.seed(seed_for(1L)); t0 <- Sys.time()
    for (k in 1:2) fit_all(gen_trial(as.list(probe[i, ])))
    as.numeric(Sys.time() - t0, units = 'secs') / 2
  })
  load <- probe |> transmute(design, N, secs) |>
    mutate(core_h = secs * 500 * 5 / 3600)
  print(as.data.frame(load), digits = 3)
  cat(sprintf('Projected %.1f core-hours; %.1f h wall on %d cores\n',
              sum(load$core_h), sum(load$core_h) / CORES, CORES))
  quit(save = 'no')
}

## ---------------- run ----------------
if (PARAMS == 'published') guard_out_dir(OUT_DIR)
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
reps <- if (SMOKE) 2L else if (!is.na(MAX_REPS)) MAX_REPS else 500L
jobs <- expand_grid(i = seq_len(nrow(t1_grid)), rep = seq_len(reps))
cat(sprintf('A4 with the baseline row: %d cells, %d (cell, rep) jobs\n',
            nrow(t1_grid), nrow(jobs)))

one_job <- function(k) {
  i <- jobs$i[k]; r <- jobs$rep[k]; cid <- t1_grid$cell_id[i]
  f <- file.path(ck_dir, sprintf('t1--%s-rep-%04d.rds', cid, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(seed_for(r))
  res <- tryCatch(fit_all(gen_trial(cell_of(i))), error = function(e) NULL)
  if (!is.null(res)) saveRDS(mutate(res, cell_id = cid, rep = r, .before = 1), f)
  invisible(NULL)
}
t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job, mc.cores = CORES,
                   mc.preschedule = FALSE))
cat(sprintf('Fits done in %.1f min\n',
            as.numeric(Sys.time() - t0, units = 'mins')))

files <- list.files(ck_dir, full.names = TRUE, pattern = '^t1--')
files <- files[as.integer(sub('.*-rep-(\\d+)\\.rds$', '\\1', files)) <= reps]
read_ck <- function(f) tryCatch(readRDS(f), error = function(e) {
  file.remove(f); NULL })
res <- bind_rows(lapply(files, read_ck)) |>
  left_join(t1_grid |> select(cell_id, panel, design, N, c_bm, t1half),
            by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))

## ---------------- guard: plain refits against job 27 ----------------
stored <- readRDS(file.path(SET_DIR, 'tier1-g-vs-f-e0/all.rds')) |>
  filter(method %in% METHODS) |>
  transmute(cell_id, rep, method = LABEL[method], est_stored = estimate,
            h_stored = h)
guard <- res |> filter(method %in% LABEL) |>
  select(cell_id, rep, method, est_new = estimate, h_new = best_t_half) |>
  inner_join(stored, by = c('cell_id', 'rep', 'method')) |>
  mutate(abs_diff = abs(est_new - est_stored),
         h_match = method != 'A4' | (is.na(h_new) & is.na(h_stored)) |
           (!is.na(h_new) & !is.na(h_stored) & h_new == h_stored),
         match = h_match & ((is.na(est_new) & is.na(est_stored)) |
           (!is.na(abs_diff) & abs_diff < GUARD_TOL)))
guard_rep <- guard |> group_by(cell_id, rep) |>
  summarise(match = all(match), max_diff = max(abs_diff, na.rm = TRUE),
            h_match = all(h_match), .groups = 'drop')
cat(sprintf(paste0('\nGuard: %d of %d replicates reproduce stored A2-A4 ',
                   '(max diff %.2g; A4 half-life mismatches %d)\n'),
            sum(guard_rep$match), nrow(guard_rep),
            max(guard_rep$max_diff), sum(!guard_rep$h_match)))
saveRDS(guard_rep, file.path(OUT_DIR, 'guard.rds'))

## ---------------- summary ----------------
ok <- guard_rep |> filter(match) |> select(cell_id, rep)
summ <- res |> semi_join(ok, by = c('cell_id', 'rep')) |>
  group_by(cell_id, panel, design, N, c_bm, t1half, method) |>
  summarise(n = sum(converged), reject = mean(p_value[converged] < 0.05),
            mcse = sqrt(reject * (1 - reject) / n),
            nonconv = sum(!converged),
            mean_est = mean(estimate[converged]),
            mean_se = mean(std_error[converged]),
            emp_se = sd(estimate[converged]), .groups = 'drop')
saveRDS(summ, file.path(OUT_DIR, 'summary.rds'))
saveRDS(list(date = Sys.time(), param_set = param_set(), reps = reps,
             guard_tol = GUARD_TOL, guard_pass = all(guard_rep$match),
             r_version = R.version.string,
             nlme = as.character(packageVersion('nlme')),
             clubSandwich = as.character(packageVersion('clubSandwich'))),
        file.path(OUT_DIR, 'meta.rds'))
cat('\nNull rejection by method, pooled over designs:\n')
print(as.data.frame(summ |> filter(c_bm == 0) |> group_by(method) |>
  summarise(reject = sum(reject * n) / sum(n), n = sum(n), .groups = 'drop')))
