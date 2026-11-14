## analysis/scripts/carryover-sensitivity/46-run-sample-size.R
##
## Power curves over the total sample size N for paper 02's recommended
## analysis A4 (stored E7null_cr2, G10 in the companion papers), with
## A2 (E1cr2, the binary indicator with CR2) as the reference, so that
## the N needed for 80 percent power can be read for each design,
## carryover half-life and effect size. A0 is not fitted: its test is
## liberal, so a "required N" for it would be meaningless.
##
## Grid: design in {CO, OLBDC, Hybrid} x true t1/2 in {0, 0.5, 1} weeks
## x c_bm in {0.2, 0.3, 0.4} x five values of N chosen per c_bm to
## bracket the N for 80 percent power (multiples of 4, so every path of
## every design gets the same count). 135 cells. Exponential decay,
## covariance moderation (mvn), predictive-only biomarker, rho = 0.7.
## 47-sample-size-curves.R fits a probit curve of power on sqrt(N) per
## (design, t1/2, c_bm, analysis) and inverts it at 0.8.
##
## Seeds: each (cell, replicate) seeds from label_seed() of its own
## label, as in 36-run-referee-variants.R, so cells are independent;
## A2 and A4 are paired within a replicate (same simulated trial).
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/46-run-sample-size.R \
##     [--params simcore2026|published] [--design CO|OLBDC|Hybrid|all] \
##     [--cbm 0.2|0.3|0.4|all] [--reps R] [--cores K] [--out-dir PATH] \
##     [--preflight] [--smoke] [--time]
##   --null runs the size check instead: c_bm = 0, t1/2 in {0, 1},
##   N in {152, 300} (simcore2026) or {300, 600} (published), so that
##   the size of A4 and A2 is known at the N of the sample-size table
##   (referee report 2026-11-14, M1); default 1,000 replicates;
##   --preflight checks every cell's covariance matrix and exits;
##   --smoke runs two replicates of every cell; --time runs three
##   replicates of the smallest and largest N on one worker and
##   reports the per-replicate wall time.
## Completed (cell, replicate) results are checkpointed under
## <out-dir>/trials/, so an interrupted run resumes.
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')
## The fitting code catches errors per fit, so a missing package would
## surface only as 100 percent non-convergence. Fail loudly instead.
for (p in c('nlme', 'clubSandwich', 'corpcor')) {
  if (!requireNamespace(p, quietly = TRUE)) stop('package missing: ', p)
}

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
DESIGN  <- argv('--design', 'all')
CBM     <- argv('--cbm', 'all')
NULL_RUN <- '--null' %in% args
REPS    <- as.integer(argv('--reps', if (NULL_RUN) '1000' else '500'))
CORES   <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
PREFLIGHT <- '--preflight' %in% args
SMOKE     <- '--smoke' %in% args
TIME      <- '--time' %in% args
PARAMS    <- Sys.getenv('PMSIM_PARAMS')
## One directory per (parameter set, design), so that separate host
## jobs never write the same result file.
OUT_DIR <- argv('--out-dir', file.path(
  'analysis/data/carryover-sensitivity',
  paste0(if (SMOKE) 'sample-size-smoke' else 'sample-size',
         if (NULL_RUN) '-null' else '', '-', PARAMS,
         if (DESIGN != 'all') paste0('-', DESIGN) else '')))
guard_out_dir(OUT_DIR)
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

SPECS <- c('E1cr2', 'E7null_cr2')
## The published set has larger component variances and needs about
## twice the N for the same power (local pilot 2026-11-11: Hybrid,
## c_bm = 0.4, t1/2 = 1, A4 power 0.86 at N = 100 under simcore2026
## against 0.61 under published), so its grid is roughly doubled.
N_GRID <- if (PARAMS == 'published') {
  list('0.4' = c(48, 80, 128, 200, 280),
       '0.3' = c(80, 140, 220, 340, 480),
       '0.2' = c(180, 300, 480, 740, 1120))
} else {
  list('0.4' = c(32, 48, 72, 100, 140),
       '0.3' = c(48, 80, 120, 172, 240),
       '0.2' = c(100, 172, 260, 380, 560))
}
stopifnot(all(unlist(N_GRID) %% 4 == 0))

N_NULL <- if (PARAMS == 'published') c(300, 600) else c(152, 300)
stopifnot(all(N_NULL %% 4 == 0))

grid <- if (NULL_RUN) {
  expand_grid(design = c('CO', 'OLBDC', 'Hybrid'), t1half = c(0, 1),
              c_bm = 0, N = N_NULL)
} else {
  expand_grid(design = c('CO', 'OLBDC', 'Hybrid'),
              t1half = c(0, 0.5, 1),
              c_bm = c(0.2, 0.3, 0.4)) |>
    mutate(N = map(c_bm, ~ N_GRID[[sprintf('%.1f', .x)]])) |>
    unnest(N)
}
grid <- grid |>
  mutate(cell_id = sprintf('%s-%s-t%.1f-cbm%.2f-N%03d',
                           if (NULL_RUN) 'ssnull' else 'ss', design, t1half,
                           c_bm, N),
         cell = pmap(list(design, t1half, c_bm, N), function(d, t, cb, n)
           list(design = d, N = n, t1half = t, c_bm = cb,
                carryover_form = 'exponential', weibull_shape = 1,
                dgp_arch = 'mvn')))
if (DESIGN != 'all') grid <- grid |> filter(design == DESIGN)
if (CBM != 'all') grid <- grid |> filter(abs(c_bm - as.numeric(CBM)) < 1e-9)
stopifnot(nrow(grid) > 0, !anyDuplicated(grid$cell_id))
grid$reps <- if (SMOKE) 2L else REPS

## ---------------- pre-flight: every cell must be PD on every path ----------------
preflight_cell <- function(cell) {
  design_set <- build_design_set(cell$design)
  mp <- list(N = cell$N, carryover_t1half = cell$t1half,
             carryover_form = cell$carryover_form,
             weibull_shape = cell$weibull_shape,
             dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
             c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0, c.bm.pb = 0)
  min(map_dbl(design_set, function(ts) {
    b <- build_sigma_matrix(mp, default_resp_param(), default_baseline_param(),
                            ts, stop_on_nonpd = FALSE)
    raw <- outer(b$standard_deviations, b$standard_deviations) * b$correlations
    pos <- b$standard_deviations > 0
    min(eigen(raw[pos, pos], symmetric = TRUE, only.values = TRUE)$values)
  }))
}
## The covariance does not depend on N, so one check per
## (design, t1/2, c_bm) suffices.
pre <- grid |>
  distinct(design, t1half, c_bm, .keep_all = TRUE) |>
  mutate(min_eigenvalue = map_dbl(cell, preflight_cell),
         feasible = min_eigenvalue > 0) |>
  select(design, t1half, c_bm, min_eigenvalue, feasible)
saveRDS(pre, file.path(OUT_DIR, 'preflight.rds'))
cat(sprintf('Parameter set %s. Pre-flight: %d settings, %d feasible, min eigenvalue %.4f\n',
            PARAMS, nrow(pre), sum(pre$feasible), min(pre$min_eigenvalue)))
if (!all(pre$feasible)) {
  print(as.data.frame(pre |> filter(!feasible)), row.names = FALSE)
  stop('Infeasible cells present; adjust the grid before running.')
}
if (PREFLIGHT) quit(save = 'no')

run_one <- function(cell_id, cell, r) {
  set.seed(label_seed(sprintf('%s-rep-%d', cell_id, r)))
  simulate_cell_g10(cell, n_reps = 1, specs = SPECS)
}

## ---------------- --time: per-replicate cost at the extremes of N ----------------
if (TIME) {
  for (id in c(grid$cell_id[which.min(grid$N)], grid$cell_id[which.max(grid$N)])) {
    k <- match(id, grid$cell_id)
    t0 <- Sys.time()
    for (r in 1:3) run_one(id, grid$cell[[k]], r)
    secs <- as.numeric(Sys.time() - t0, units = 'secs') / 3
    cat(sprintf('%s: %.2f s per replicate on one worker\n', id, secs))
  }
  quit(save = 'no')
}

jobs <- grid |> select(cell_id, reps) |>
  mutate(rep = map(reps, seq_len)) |> select(-reps) |> unnest(rep) |>
  ## Largest cells first, so the run does not end on a few slow jobs.
  left_join(grid |> select(cell_id, N), by = 'cell_id') |>
  arrange(desc(N)) |> select(-N)
cat(sprintf('Sample-size grid: %d cells, %d (cell, rep) jobs\n',
            nrow(grid), nrow(jobs)))

one_job <- function(i) {
  cell_id <- jobs$cell_id[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  res <- tryCatch(run_one(cell_id, grid$cell[[match(cell_id, grid$cell_id)]], r),
                  error = function(e) NULL)
  if (!is.null(res)) saveRDS(mutate(res, rep = r, cell_id = cell_id,
                                    .before = 1), f)
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job,
                   mc.cores = CORES, mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

files <- list.files(ck_dir, full.names = TRUE)
files <- files[sub('-rep-\\d+\\.rds$', '', basename(files)) %in% grid$cell_id]
res <- bind_rows(lapply(files, readRDS)) |>
  left_join(grid |> select(cell_id, design, t1half, c_bm, N), by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))

summary_tab <- res |>
  group_by(cell_id, design, t1half, c_bm, N, method) |>
  summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
            n = sum(!is.na(p_value)),
            mcse = sqrt(reject * (1 - reject) / n),
            nonconv = mean(!converged),
            mean_se = mean(std_error, na.rm = TRUE),
            emp_se = sd(estimate, na.rm = TRUE), .groups = 'drop') |>
  mutate(kappa = mean_se / emp_se)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))
saveRDS(list(params = PARAMS, reps = REPS, specs = SPECS,
             n_grid = if (NULL_RUN) N_NULL else N_GRID, null_run = NULL_RUN,
             date = Sys.time(), elapsed_min = elapsed_min),
        file.path(OUT_DIR, 'meta.rds'))
cat(sprintf('%d of %d (cell, rep) jobs present; %.1f min this run\n',
            length(unique(paste(res$cell_id, res$rep))), nrow(jobs),
            elapsed_min))
