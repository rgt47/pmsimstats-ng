## 19-g9-planning-and-fig4b.R
##
## Power of analysis G9 (exposure-weighted Dbc, half-life chosen by
## ML-AIC over {0.25, 0.5, 1, 2} weeks with the interaction in the
## selection model, REML refit, CR2 cluster-robust SE with Satterthwaite
## df, random intercept plus corCAR1; stored as `E7cr2`) under the
## extended DGP, for two inputs of paper 01. The fit is
## fit_spec_s7(dat_long, 'E7cr2') of
## analysis/scripts/carryover-sensitivity/simulation-core.R, unchanged.
##
##   P  planning column. Designs Hybrid, OL+BDC, CO; c_bm in
##      {0.25, 0.35, 0.45}; t1/2 in {0, 0.2, 1} weeks; separable AR(1)
##      rho 0.7, c_1 = c_x = 0.2, graded coupling, anchored carryover,
##      stop on non-PD, fixed allocation to paths. Each cell runs at
##      N = 70 and at the total N that the closed-form E9 needs for 80
##      percent power (planning-n80.csv of 17-rho-stage.R, capped at
##      N = 400; a capped cell runs N = 400 and is flagged). Alternatives
##      at --reps-alt (1,000), plus the null (c_bm = 0) at N = 70 for
##      each design and half-life at --reps-null (2,000).
##   F  Figure 4B panel. The original study's grid under the extended
##      DGP and G9: OL, OL+BDC, CO, N-of-1 (= Hybrid), N in {35, 70},
##      t1/2 in {0, 0.1, 0.2}, the four original censoring sets plus no
##      censoring, c_bm = 0.25, --reps-fig (250) replicates per cell.
##      Each replicate draws one data set and applies the five
##      censoring settings to it, as generateSimulatedResults() does in
##      58b32a9. Censoring uses censor_data(), the tidyverse port of
##      censordata(), on the full data set with the first path's design
##      (20-check-g9-planning.R compares it with the 58b32a9 function).
##      The open-label design has no off-drug contrast, so Dbc is
##      constant at 1 and bm:Dbc is aliased with bm: G9 is not
##      estimable. The original analysis fits bm*t in that case
##      (lme_analysis.R, varInDb == FALSE); G9 has no such fallback and
##      the OL cells are recorded as not estimable, not simulated.
##
## Seeds: every (cell, replicate) job seeds from
## label_seed(sprintf('%s-rep-%d', cell_id, r)), so no cell shares a
## random-number stream with another cell or with any stored run.
##
## --specs E7null_cr2 fits G10 instead (half-life chosen by ML-AIC on
## the null model Sx ~ bm + t + Dbc, REML refit with bm:Dbc, CR2) on
## the same simulated trials, and writes to g10-planning/. With
## --guard-reps R the first R replicates of every cell also refit G9
## (E7cr2); the summaries compare that refit with the stored G9 runs
## (block P 2026-10-28; block F zeus job 56, 2026-10-29) to confirm
## the trials are the same ones.
##
## Usage (repository root):
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/19-g9-planning-and-fig4b.R \
##     [--block P|F|all] [--reps-alt R] [--reps-null R] [--reps-fig R] \
##     [--cores K] [--out-dir PATH] [--smoke] [--preflight] [--time] \
##     [--cross-check [--reps-cross R]] [--specs S] [--guard-reps R]
##   --preflight checks every covariance matrix and exits; --smoke runs
##   2 replicates per cell; --time reports per-replicate cost on one
##   worker and a projected host wall time; --cross-check runs the
##   Hybrid, c_bm 0.45, t1/2 0.5 week, N = 70 cell of paper 02 block A
##   (default 200 replicates) and compares it with the stored value.
##
## Outputs in --out-dir: all.rds (per replicate: estimate, std_error,
## p_value, converged, h = selected half-life), summary.rds (list with
## elements P and F), meta.rds, preflight.rds, planning-g9.csv (block
## P), fig4b-g9-cbm<C>.csv (block F), trials/ (gitignored checkpoints).
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel); library(clubSandwich)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
BLOCK      <- argv('--block', 'all')
REPS_ALT   <- as.integer(argv('--reps-alt', '1000'))
REPS_NULL  <- as.integer(argv('--reps-null', '2000'))
REPS_FIG   <- as.integer(argv('--reps-fig', '250'))
REPS_CROSS <- as.integer(argv('--reps-cross', '200'))
FIG_CBM    <- as.numeric(argv('--cbm-fig', '0.25'))
CORES      <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
PREFLIGHT  <- '--preflight' %in% args
SMOKE      <- '--smoke' %in% args
TIME       <- '--time' %in% args
CROSS      <- '--cross-check' %in% args
SPEC       <- argv('--specs', 'E7cr2')
GUARD_REPS <- as.integer(argv('--guard-reps', '0'))
G_LABEL <- c(E7cr2 = 'G9', E7null_cr2 = 'G10')
stopifnot(SPEC %in% names(G_LABEL))
LABEL <- unname(G_LABEL[SPEC])
G9_DIR <- 'analysis/data/quick-sim/carryover-closed-form/g9-planning'
BASE_DIR <- file.path(dirname(G9_DIR), sprintf('%s-planning', tolower(LABEL)))
OUT_DIR <- argv('--out-dir', if (SMOKE) file.path(BASE_DIR, 'smoke') else BASE_DIR)
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(BLOCK %in% c('P', 'F', 'all'))

PLAN_FILE <- 'analysis/data/quick-sim/carryover-closed-form/rho-stage/planning-n80.csv'
N_CAP <- 400L
RHO <- 0.7
HOST_CORES <- 8
HOST_SPEEDUP <- 2

## Parameters of the four original censoring sets, from
## results_core$parameterselections$censorparams in 58b32a9; the check
## script compares them with the stored object.
CENSOR <- tibble(
  names = c('balanced', 'more of flat', 'more of biased', 'high dropout'),
  beta0 = c(0.05, 0.05, 0.05, 0.15),
  beta1 = c(0.5, 0.2, 0.8, 0.5),
  eb1   = c(2, 2, 2, 2))
CENSOR_LABELS <- c('No censoring', CENSOR$names)

PRESET <- c(Hybrid = 'Hybrid', `OL+BDC` = 'OLBDC', CO = 'CO', OL = 'OL',
            `N-of-1` = 'Hybrid')
DESIGN_SETS <- map(set_names(unique(PRESET)), build_design_set)

## ---------------- grids ----------------
plan <- read.csv(PLAN_FILE) |> as_tibble() |> mutate(N80 = as.integer(N80))
stopifnot(nrow(plan) == 27, all(plan$design %in% names(PRESET)))

p_common <- function(d) d |> mutate(block = 'P', simulate = TRUE)
alt_70 <- plan |>
  transmute(design, c_bm, t_half, N = 70L, N_label = '70', N80, capped = FALSE,
            e9_power = power_N70, reps = REPS_ALT)
alt_80 <- plan |>
  transmute(design, c_bm, t_half,
            N = as.integer(pmin(coalesce(N80, N_CAP), N_CAP)),
            N_label = 'N80_E9', N80,
            capped = is.na(N80) | N80 > N_CAP,
            e9_power = if_else(is.na(N80) | N80 > N_CAP, NA_real_, 0.8),
            reps = REPS_ALT)
null_70 <- plan |>
  distinct(design, t_half) |>
  mutate(c_bm = 0, N = 70L, N_label = '70', N80 = NA_integer_,
         capped = FALSE, e9_power = NA_real_, reps = REPS_NULL)
grid_p <- bind_rows(alt_70, alt_80, null_70) |>
  p_common() |>
  mutate(preset = unname(PRESET[design]),
         cell_id = sprintf('P-%s-cbm%.2f-t%g-%s', design, c_bm, t_half, N_label))

grid_f <- expand_grid(design = c('OL', 'OL+BDC', 'CO', 'N-of-1'),
                      N = c(35L, 70L), t_half = c(0, 0.1, 0.2)) |>
  mutate(block = 'F', preset = unname(PRESET[design]), c_bm = FIG_CBM,
         N_label = as.character(N), N80 = NA_integer_, capped = FALSE,
         e9_power = NA_real_, reps = REPS_FIG, simulate = design != 'OL',
         cell_id = sprintf('F-%s-N%d-t%g-cbm%.2f', design, N, t_half, c_bm))

grid_x <- tibble(block = 'X', design = 'Hybrid', preset = 'Hybrid', N = 70L,
                 t_half = 0.5, c_bm = 0.45, N_label = '70',
                 N80 = NA_integer_, capped = FALSE, e9_power = NA_real_,
                 reps = REPS_CROSS, simulate = TRUE,
                 cell_id = 'X-Hybrid-cbm0.45-t0.5-70')

grid <- if (CROSS) grid_x else bind_rows(grid_p, grid_f)
if (!CROSS && BLOCK != 'all') grid <- filter(grid, block == BLOCK)
stopifnot(nrow(grid) > 0, !anyDuplicated(grid$cell_id))
if (SMOKE) grid <- mutate(grid, reps = 2L)

## ---------------- model parameters, data, fit ----------------
make_mp <- function(N, t_half, c_bm) {
  list(N = N, carryover_t1half = t_half, carryover_form = 'exponential',
       weibull_shape = 1, dgp_architecture = 'mvn', c.bm = c_bm,
       c.tv = RHO, c.pb = RHO, c.br = RHO, c.cf1t = 0.2, c.cfct = 0.2,
       c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0, b.bm.tv.trend = 0,
       sd_slope = 0, c.bm.pb = 0)
}

NA_FIT <- tibble(estimate = NA_real_, std_error = NA_real_,
                 p_value = NA_real_, h = NA_real_, converged = FALSE)

fit_one <- function(dat_long, spec) {
  r <- tryCatch(fit_spec_s7(dat_long, spec), error = function(e) NULL)
  if (is.null(r)) return(NA_FIT)
  transmute(r, estimate, std_error, p_value, h = best_t_half, converged)
}

one_rep <- function(cell, r) {
  set.seed(label_seed(sprintf('%s-rep-%d', cell$cell_id, r)))
  ds <- DESIGN_SETS[[cell$preset]]
  sets <- if (cell$block == 'F') 0:4 else NA_integer_
  specs <- if (r <= GUARD_REPS) union(SPEC, 'E7cr2') else SPEC
  dat <- tryCatch(
    generate_data_multi_path(make_mp(cell$N, cell$t_half, cell$c_bm),
                             default_resp_param(), default_baseline_param(),
                             ds, stop_on_nonpd = TRUE),
    error = function(e) NULL)
  map_dfr(sets, function(k) {
    out <- if (is.null(dat) || nrow(dat) == 0) {
      map_dfr(specs, ~ mutate(NA_FIT, method = .x))
    } else {
      d <- if (is.na(k) || k == 0) dat else censor_data(dat, ds[[1]], CENSOR[k, ])
      dl <- prepare_long_data(d, ds, carryover_t1half = cell$t_half)
      map_dfr(specs, ~ mutate(fit_one(dl, .x), method = .x))
    }
    mutate(out, censorparamset = k)
  }) |>
    mutate(cell_id = cell$cell_id, rep = r, .before = 1)
}

## ---------------- pre-flight: every cell PD on every path ----------------
preflight_cell <- function(preset, c_bm, t_half) {
  ds <- DESIGN_SETS[[preset]]
  mp <- make_mp(70L, t_half, c_bm)
  map_dfr(seq_along(ds), function(i) {
    b <- build_sigma_matrix(mp, default_resp_param(), default_baseline_param(),
                            ds[[i]], stop_on_nonpd = FALSE)
    raw <- outer(b$standard_deviations, b$standard_deviations) * b$correlations
    pos <- b$standard_deviations > 0
    ev <- eigen(raw[pos, pos], symmetric = TRUE, only.values = TRUE)$values
    tibble(path = i, min_eigenvalue = min(ev), feasible = min(ev) > 0)
  })
}
pre <- grid |>
  distinct(block, design, preset, c_bm, t_half) |>
  mutate(pf = pmap(list(preset, c_bm, t_half), preflight_cell)) |>
  unnest(pf) |>
  group_by(block, design, c_bm, t_half) |>
  summarise(min_eigenvalue = min(min_eigenvalue),
            feasible = all(feasible), .groups = 'drop')
saveRDS(pre, file.path(OUT_DIR, 'preflight.rds'))
cat(sprintf('Pre-flight: %d (design, c_bm, t1/2) cells, %d feasible, min eigenvalue %.4f\n',
            nrow(pre), sum(pre$feasible), min(pre$min_eigenvalue)))
if (!all(pre$feasible)) {
  print(as.data.frame(filter(pre, !feasible)), row.names = FALSE)
  stop('Infeasible cells present; adjust the grid before running.')
}
if (PREFLIGHT) {
  print(as.data.frame(pre |> mutate(min_eigenvalue = round(min_eigenvalue, 4))),
        row.names = FALSE)
  quit(save = 'no')
}

## ---------------- --time ----------------
## The per-replicate cost is close to a power law in N (the lme and
## clubSandwich steps), so it is measured at N = 70 and at N = 400 and
## interpolated; the projection ignores design differences.
if (TIME) {
  time_cell <- function(cell, n) {
    t0 <- Sys.time()
    for (r in seq_len(n)) one_rep(cell, r)
    as.numeric(Sys.time() - t0, units = 'secs') / n
  }
  sim_p <- filter(grid_p, simulate)
  heavy_p <- filter(sim_p, design == 'Hybrid', N == max(N))[1, ]
  ref_p <- filter(sim_p, design == 'Hybrid', c_bm == 0.45, t_half == 0,
                  N_label == '70')[1, ]
  heavy_f <- filter(grid_f, design == 'N-of-1', N == 70, t_half == 0.2)[1, ]
  s_p_heavy <- time_cell(heavy_p, 3)
  s_p_ref <- time_cell(ref_p, 5)
  s_f_heavy <- time_cell(heavy_f, 3)
  cat(sprintf('Block P heaviest cell %s (N = %d): %.2f s per replicate, one worker\n',
              heavy_p$cell_id, heavy_p$N, s_p_heavy))
  cat(sprintf('Block P reference cell %s (N = %d): %.2f s per replicate\n',
              ref_p$cell_id, ref_p$N, s_p_ref))
  cat(sprintf('Block F heaviest cell %s (N = %d, five censoring analyses per replicate): %.2f s per replicate\n',
              heavy_f$cell_id, heavy_f$N, s_f_heavy))
  b <- log(s_p_heavy / s_p_ref) / log(heavy_p$N / ref_p$N)
  cost_p <- function(N) s_p_ref * (N / ref_p$N)^b
  cost_f <- function(N) s_f_heavy * (N / heavy_f$N)^b
  proj <- bind_rows(
    sim_p |> transmute(block, sec = reps * cost_p(N)),
    filter(grid_f, simulate) |> transmute(block, sec = reps * cost_f(N)))
  tot <- proj |> group_by(block) |> summarise(local_core_h = sum(sec) / 3600)
  cat(sprintf('Fitted exponent of cost in N: %.2f\n', b))
  for (i in seq_len(nrow(tot))) {
    cat(sprintf('Block %s: %.1f local core-hours; host wall time on %d cores at %gx per-core speed: %.1f h\n',
                tot$block[i], tot$local_core_h[i], HOST_CORES, HOST_SPEEDUP,
                tot$local_core_h[i] / HOST_SPEEDUP / HOST_CORES))
  }
  quit(save = 'no')
}

cat(sprintf('%s planning and Figure 4B, block %s: %d cells, %d (cell, rep) jobs\n',
            LABEL, if (CROSS) 'cross-check' else BLOCK, nrow(grid),
            sum(grid$reps[grid$simulate])))

jobs <- grid |>
  filter(simulate) |>
  arrange(desc(N)) |>
  select(cell_id, reps) |>
  mutate(rep = map(reps, seq_len)) |>
  select(-reps) |>
  unnest(rep)

one_job <- function(i) {
  cell_id <- jobs$cell_id[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%05d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  res <- tryCatch(one_rep(grid[match(cell_id, grid$cell_id), ], r),
                  error = function(e) NULL)
  if (!is.null(res)) saveRDS(res, f)
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job,
                   mc.cores = CORES, mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

files <- list.files(ck_dir, full.names = TRUE)
files <- files[sub('-rep-\\d+\\.rds$', '', basename(files)) %in% grid$cell_id]
meta_cols <- c('cell_id', 'block', 'design', 'N', 't_half', 'c_bm', 'N_label')
res <- bind_rows(lapply(files, readRDS)) |>
  left_join(grid[meta_cols], by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
cat(sprintf('\n%d of %d (cell, rep) jobs present in %.1f min (this run)\n',
            n_distinct(paste(res$cell_id, res$rep)), nrow(jobs), elapsed_min))

## The G9 refit on the first --guard-reps replicates must reproduce the
## stored G9 run; block P and block F were stored in separate runs.
## g9-planning/block-F/all.rds predates the censordata() sign fix
## (24f187bc6); the block F that paper 01 reports is zeus job 56, run
## after the fix, so the guard reads that one.
G9_F_ALL <- file.path('analysis/remote_runs/job-56-24f187b/results',
                      G9_DIR, 'block-F', 'all.rds')
guard <- NULL
if (SPEC != 'E7cr2' && GUARD_REPS > 0 && !SMOKE && !CROSS) {
  stored_g9 <- bind_rows(
    readRDS(file.path(G9_DIR, 'all.rds')),
    readRDS(G9_F_ALL)) |>
    select(cell_id, rep, censorparamset, p_stored = p_value, h_stored = h)
  guard <- res |>
    filter(method == 'E7cr2') |>
    select(cell_id, block, rep, censorparamset, p_value, h) |>
    left_join(stored_g9, by = c('cell_id', 'rep', 'censorparamset')) |>
    mutate(abs_diff = abs(p_value - p_stored),
           match = !is.na(abs_diff) & abs_diff < 1e-6 &
             abs(h - h_stored) < 1e-9,
           flip = (p_value < 0.05) != (p_stored < 0.05))
  saveRDS(guard, file.path(OUT_DIR, 'guard.rds'))
  cat(sprintf('Guard (G9 refit against stored G9): %d of %d fits match within 1e-6; %d missing in the stored run; %d rejection flips\n',
              sum(guard$match), nrow(guard), sum(is.na(guard$p_stored)),
              sum(guard$flip, na.rm = TRUE)))
}
res <- filter(res, method == SPEC)

## ---------------- summaries ----------------
summ_power <- function(d) {
  d |>
    summarise(n_reps = n_distinct(rep),
              n_p = sum(!is.na(p_value)),
              power = mean(p_value < 0.05, na.rm = TRUE),
              mean_estimate = mean(estimate, na.rm = TRUE),
              mean_se = mean(std_error, na.rm = TRUE),
              mean_selected_half_life = mean(h[converged], na.rm = TRUE),
              nonconv = mean(!converged), .groups = 'drop') |>
    mutate(mcse = sqrt(power * (1 - power) / n_p))
}

sum_p <- res |>
  filter(block %in% c('P', 'X')) |>
  group_by(cell_id) |>
  summ_power() |>
  left_join(grid |> select(cell_id, design, c_bm, t_half, N_label, N, N80,
                           capped, e9_power, reps), by = 'cell_id') |>
  mutate(size = if_else(c_bm == 0, power, NA_real_)) |>
  rename(N80_E9 = N80, n_reps_requested = reps) |>
  arrange(match(design, c('Hybrid', 'OL+BDC', 'CO')), c_bm, t_half,
          desc(N_label == '70'))

fig_rows <- res |>
  filter(block == 'F') |>
  group_by(cell_id, design, N, t_half, censorparamset) |>
  summ_power() |>
  transmute(censorparamset, power, mean_beta = mean_estimate, n = n_reps,
            version = paste0('extended-', LABEL), design, N, t_half,
            repairs = 0L,
            censoring = CENSOR_LABELS[censorparamset + 1],
            estimable = TRUE, n_conv = n_p, nonconv,
            mean_selected_half_life) |>
  ungroup()
ol_rows <- grid |>
  filter(block == 'F', !simulate) |>
  select(design, N, t_half) |>
  crossing(censorparamset = 0:4) |>
  mutate(power = NA_real_, mean_beta = NA_real_, n = 0L,
         version = paste0('extended-', LABEL), repairs = 0L,
         censoring = CENSOR_LABELS[censorparamset + 1], estimable = FALSE,
         n_conv = 0L, nonconv = NA_real_, mean_selected_half_life = NA_real_)
sum_f <- bind_rows(fig_rows, ol_rows) |>
  select(censorparamset, power, mean_beta, n, version, design, N, t_half,
         repairs, censoring, estimable, n_conv, nonconv,
         mean_selected_half_life) |>
  arrange(match(design, c('OL', 'OL+BDC', 'CO', 'N-of-1')), N, t_half,
          censorparamset)

if (nrow(sum_p) > 0 && !CROSS) {
  write.csv(sum_p |>
              select(design, c_bm, t_half, N_label, N, N80_E9, capped,
                     n_reps, n_reps_requested, power, mcse, size,
                     mean_selected_half_life, nonconv, e9_power),
            file.path(OUT_DIR, sprintf('planning-%s.csv', tolower(LABEL))),
            row.names = FALSE)
  cat('\nBlock P (power; size where c_bm = 0):\n')
  print(as.data.frame(sum_p |>
    select(design, c_bm, t_half, N_label, N, capped, n_reps, power, mcse,
           mean_selected_half_life, nonconv) |>
    mutate(across(c(power, mcse, mean_selected_half_life, nonconv),
                  ~ round(.x, 3)))), max = 4000)
}
if (nrow(sum_f) > 0 && !CROSS && any(grid$block == 'F')) {
  write.csv(sum_f, file.path(OUT_DIR, sprintf('fig4b-%s-cbm%.2f.csv',
                                              tolower(LABEL), FIG_CBM)),
            row.names = FALSE, na = 'NA')
  cat('\nBlock F (power; head):\n')
  print(as.data.frame(head(sum_f |> mutate(across(where(is.double), ~ round(.x, 3))), 20)),
        row.names = FALSE)
}
saveRDS(list(P = sum_p, F = sum_f), file.path(OUT_DIR, 'summary.rds'))

if (CROSS) {
  stored <- readRDS('analysis/data/carryover-sensitivity/referee-variants/block-A/summary.rds') |>
    filter(cell_id == 'alt-t0.5-N70-cbm0.45', method == 'E7cr2')
  cat(sprintf('\nCross-check, Hybrid, c_bm 0.45, t1/2 0.5, N 70:\n  this driver: power %.3f (MCSE %.3f, %d replicates, non-convergence %.3f)\n  stored paper 02 block A: power %.3f (MCSE %.3f, %d replicates)\n',
              sum_p$power, sum_p$mcse, sum_p$n_p, sum_p$nonconv,
              stored$reject, stored$mcse, stored$n))
}

meta <- list(
  script = '19-g9-planning-and-fig4b.R', block = BLOCK, date = Sys.time(),
  reps_alt = REPS_ALT, reps_null = REPS_NULL, reps_fig = REPS_FIG,
  fig_cbm = FIG_CBM, n_cap = N_CAP, smoke = SMOKE, cross_check = CROSS,
  spec = SPEC, guard_reps = GUARD_REPS,
  elapsed_min_this_run = elapsed_min, r_version = R.version.string,
  git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                     error = function(e) NA_character_)
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
cat(sprintf('\nWritten: %s\n', OUT_DIR))
