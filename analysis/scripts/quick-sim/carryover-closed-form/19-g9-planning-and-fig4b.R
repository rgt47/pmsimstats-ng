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
##   N  G10 sample-size search (referee report 2026-11-01, M5; block N
##      runs only when asked for, not under --block all). The Hybrid
##      cells of block P (c_bm in {0.25, 0.35, 0.45}, t1/2 in
##      {0, 0.2, 1}) at N = m * N80(E9) for m in --n80-mult (default
##      0.35, 0.5, 0.65, 0.8, 1, 1.25), --reps-n80 (500) replicates each. For
##      each cell and analysis a probit regression of rejection on
##      log N gives the N at which power is 0.80 (n80-<g>.csv; the
##      per-N rejection rates are in n80-grid-<g>.csv). A cell whose
##      0.80 crossing lies outside the simulated N range is flagged
##      as extrapolated.
##
## --specs E7null_cr2 fits G10 (half-life chosen by ML-AIC on the null
## model Sx ~ bm + t + Dbc, REML refit with bm:Dbc, CR2);
## --specs E7null_cr2,E7cr2 fits G10 and G9 on the same simulated
## trials and writes the planning and Figure 4B files of each.
##
## Parameter set (simulation-core.R): --params published|simcore2026.
## Under 'published' the default --out-dir is published-planning/;
## under 'simcore2026' it is <g>-planning/ for a single analysis. The
## stored runs (g9-planning/ 2026-10-28, zeus job 56 for block F, and
## g10-planning/) are simcore2026. --guard-reps R refits G9 on the
## first R replicates and compares it with those stored G9 runs, so it
## is refused under any other set; --cross-check compares with paper
## 02's block A of the same set (the published one once docs/71's
## p02-ref-A job is applied).
##
## Usage (repository root):
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/19-g9-planning-and-fig4b.R \
##     [--block P|F|N|all] [--reps-alt R] [--reps-null R] [--reps-fig R] \
##     [--reps-n80 R] [--n80-mult M1,M2,...] \
##     [--cores K] [--out-dir PATH] [--smoke] [--preflight] [--time] \
##     [--cross-check [--reps-cross R]] [--specs S1[,S2]] [--guard-reps R] \
##     [--params published|simcore2026] [--rho R --plan-file PATH]
##   --rho and --plan-file (block P only) run the planning cells at
##   another serial correlation, with the grid of a planning file from
##   24-rho-sweep.R; cells it flags infeasible are skipped. Output goes
##   to published-planning/rho<R>/ by default.
##   --preflight checks every covariance matrix and exits; --smoke runs
##   2 replicates per cell; --time reports per-replicate cost on one
##   worker and a projected host wall time; --cross-check runs the
##   Hybrid, c_bm 0.45, t1/2 0.5 week, N = 70 cell of paper 02 block A
##   (default 200 replicates) and compares it with the stored value.
##
## Outputs in --out-dir: all.rds (per replicate and analysis: estimate,
## std_error, p_value, converged, h = selected half-life, method),
## summary.rds (list with elements P, F and N, or one such list per
## analysis when several are fitted), meta.rds, preflight.rds,
## planning-<g>.csv (block P), fig4b-<g>-cbm<C>.csv (block F),
## n80-<g>.csv and n80-grid-<g>.csv (block N), PARAM_SET, trials/
## (gitignored checkpoints); <g> is g9 or g10.
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
REPS_N80   <- as.integer(argv('--reps-n80', '500'))
N80_MULT   <- as.numeric(strsplit(argv('--n80-mult', '0.35,0.5,0.65,0.8,1,1.25'),
                                  ',', fixed = TRUE)[[1]])
FIG_CBM    <- as.numeric(argv('--cbm-fig', '0.25'))
CORES      <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
PREFLIGHT  <- '--preflight' %in% args
SMOKE      <- '--smoke' %in% args
TIME       <- '--time' %in% args
CROSS      <- '--cross-check' %in% args
SPECS      <- strsplit(argv('--specs', 'E7cr2'), ',', fixed = TRUE)[[1]]
GUARD_REPS <- as.integer(argv('--guard-reps', '0'))
RHO        <- as.numeric(argv('--rho', '0.7'))
G_LABEL <- c(E7cr2 = 'G9', E7null_cr2 = 'G10')
stopifnot(length(SPECS) >= 1, all(SPECS %in% names(G_LABEL)),
          !anyDuplicated(SPECS), !anyNA(N80_MULT), all(N80_MULT > 0))
LABELS <- unname(G_LABEL[SPECS])
G9_DIR <- 'analysis/data/quick-sim/carryover-closed-form/g9-planning'
## The stored G9 and G10 runs are simcore2026; a published-set run gets
## its own directory so neither can resume from the other's checkpoints.
## Runs at another rho reuse the cell ids, and so the seeds, of the
## rho = 0.7 run; their own subdirectory keeps the checkpoints apart.
BASE_DIR <- if (param_set() == 'published' && RHO != 0.7) {
  file.path(dirname(G9_DIR), 'published-planning', sprintf('rho%g', RHO))
} else if (param_set() == 'published') {
  file.path(dirname(G9_DIR), 'published-planning')
} else if (length(SPECS) == 1) {
  file.path(dirname(G9_DIR), sprintf('%s-planning', tolower(LABELS)))
} else NA_character_
## Block N runs as its own job; its own subdirectory keeps it from
## rewriting the all.rds and meta.rds of the block P and F run.
if (BLOCK == 'N' && !is.na(BASE_DIR)) BASE_DIR <- file.path(BASE_DIR, 'block-N')
OUT_DIR <- argv('--out-dir',
  if (SMOKE && !is.na(BASE_DIR)) file.path(BASE_DIR, 'smoke') else BASE_DIR)
if (is.na(OUT_DIR))
  stop('several --specs under ', param_set(), ' need an explicit --out-dir')
if (GUARD_REPS > 0 && param_set() != 'simcore2026')
  stop('--guard-reps compares with the stored simcore2026 G9 runs; ',
       'it does not apply under --params ', param_set())
guard_out_dir(OUT_DIR)
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(BLOCK %in% c('P', 'F', 'N', 'all'))

PLAN_FILE <- argv('--plan-file',
  'analysis/data/quick-sim/carryover-closed-form/rho-stage/planning-n80.csv')
N_CAP <- 400L
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
## Planning files for rho other than 0.7 (24-rho-sweep.R) flag cells
## above the feasibility bound; those are not simulated.
if ('feasible' %in% names(plan)) plan <- filter(plan, feasible)
stopifnot(nrow(plan) > 0, all(plan$design %in% names(PRESET)))
if (RHO != 0.7 && BLOCK != 'P')
  stop('--rho other than 0.7 is set up for block P only')

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

## Block N: the Hybrid alternatives of block P at multiples of the E9
## sample size; cells with no E9 value (above 800) are left out.
grid_n <- plan |>
  filter(design == 'Hybrid', c_bm > 0, !is.na(N80)) |>
  select(design, c_bm, t_half, N80) |>
  crossing(mult = N80_MULT) |>
  mutate(N = pmax(20L, as.integer(round(N80 * mult))),
         block = 'N', preset = unname(PRESET[design]),
         N_label = sprintf('x%g', mult), capped = FALSE,
         e9_power = NA_real_, reps = REPS_N80, simulate = TRUE,
         cell_id = sprintf('N-%s-cbm%.2f-t%g-N%d', design, c_bm, t_half, N)) |>
  distinct(cell_id, .keep_all = TRUE) |>
  select(-mult)

grid <- if (CROSS) grid_x else
  if (BLOCK == 'N') grid_n else bind_rows(grid_p, grid_f)
if (!CROSS && BLOCK %in% c('P', 'F')) grid <- filter(grid, block == BLOCK)
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
  specs <- if (r <= GUARD_REPS) union(SPECS, 'E7cr2') else SPECS
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
  ref_p <- filter(sim_p, design == 'Hybrid', c_bm == max(c_bm), t_half == 0,
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
  ## one_rep() fits every analysis in --specs, so the timings already
  ## include all of them.
  proj <- bind_rows(
    sim_p |> transmute(block, sec = reps * cost_p(N)),
    filter(grid_f, simulate) |> transmute(block, sec = reps * cost_f(N)),
    grid_n |> transmute(block, sec = reps * cost_p(N)))
  tot <- proj |> group_by(block) |> summarise(local_core_h = sum(sec) / 3600)
  cat(sprintf('Fitted exponent of cost in N: %.2f\n', b))
  for (i in seq_len(nrow(tot))) {
    cat(sprintf('Block %s: %.1f local core-hours; host wall time on %d cores at %gx per-core speed: %.1f h\n',
                tot$block[i], tot$local_core_h[i], HOST_CORES, HOST_SPEEDUP,
                tot$local_core_h[i] / HOST_SPEEDUP / HOST_CORES))
  }
  quit(save = 'no')
}

cat(sprintf('%s planning and Figure 4B, block %s, parameter set %s: %d cells, %d (cell, rep) jobs\n',
            paste(LABELS, collapse = ' and '),
            if (CROSS) 'cross-check' else BLOCK, param_set(), nrow(grid),
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
res_all <- bind_rows(lapply(files, readRDS)) |>
  left_join(grid[meta_cols], by = 'cell_id')
## A checkpoint written by a run with other --specs lacks this run's
## analyses; refuse rather than summarize a partial set.
missing_specs <- setdiff(SPECS, unique(res_all$method))
if (length(missing_specs))
  stop('checkpoints in ', ck_dir, ' lack ', paste(missing_specs, collapse = ', '),
       '; they were written with other --specs. Use another --out-dir.')
saveRDS(res_all, file.path(OUT_DIR, 'all.rds'))
cat(sprintf('\n%d of %d (cell, rep) jobs present in %.1f min (this run)\n',
            n_distinct(paste(res_all$cell_id, res_all$rep)), nrow(jobs),
            elapsed_min))

## The G9 refit on the first --guard-reps replicates must reproduce the
## stored G9 run; block P and block F were stored in separate runs.
## g9-planning/block-F/all.rds predates the censordata() sign fix
## (24f187bc6); the block F that paper 01 reports is zeus job 56, run
## after the fix, so the guard reads that one. Both are simcore2026
## (refused above under any other set).
G9_F_ALL <- file.path('analysis/remote_runs/job-56-24f187b/results',
                      G9_DIR, 'block-F', 'all.rds')
guard <- NULL
if (!'E7cr2' %in% SPECS && GUARD_REPS > 0 && !SMOKE && !CROSS) {
  stored_g9 <- bind_rows(
    readRDS(file.path(G9_DIR, 'all.rds')),
    readRDS(G9_F_ALL)) |>
    select(cell_id, rep, censorparamset, p_stored = p_value, h_stored = h)
  guard <- res_all |>
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

## ---------------- summaries ----------------
summ_power <- function(d) {
  d |>
    summarise(n_reps = n_distinct(rep),
              n_p = sum(!is.na(p_value)),
              n_reject = sum(p_value < 0.05, na.rm = TRUE),
              power = mean(p_value < 0.05, na.rm = TRUE),
              mean_estimate = mean(estimate, na.rm = TRUE),
              mean_se = mean(std_error, na.rm = TRUE),
              mean_selected_half_life = mean(h[converged], na.rm = TRUE),
              nonconv = mean(!converged), .groups = 'drop') |>
    mutate(mcse = sqrt(power * (1 - power) / n_p))
}

## N at which power reaches `target`, from a probit regression of
## rejection on log N over the simulated N of one cell. The crossing
## is flagged when it lies outside the simulated range.
fit_n80 <- function(d, target = 0.8) {
  fit <- tryCatch(
    suppressWarnings(glm(cbind(n_reject, n_p - n_reject) ~ log(N),
                         family = binomial(link = 'probit'), data = d)),
    error = function(e) NULL)
  ok <- !is.null(fit) && all(is.finite(coef(fit))) && coef(fit)[2] > 0
  n_hat <- if (ok) exp((qnorm(target) - coef(fit)[1]) / coef(fit)[2]) else NA_real_
  tibble(N80_G = unname(n_hat), N_min = min(d$N), N_max = max(d$N),
         n_points = nrow(d),
         extrapolated = is.na(n_hat) | n_hat < min(d$N) | n_hat > max(d$N),
         power_at_min = d$power[which.min(d$N)],
         power_at_max = d$power[which.max(d$N)])
}

summarize_spec <- function(spec) {
  label <- unname(G_LABEL[spec])
  res <- filter(res_all, method == spec)

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
              version = paste0('extended-', label), design, N, t_half,
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
           version = paste0('extended-', label), repairs = 0L,
           censoring = CENSOR_LABELS[censorparamset + 1], estimable = FALSE,
           n_conv = 0L, nonconv = NA_real_, mean_selected_half_life = NA_real_)
  sum_f <- bind_rows(fig_rows, ol_rows) |>
    select(censorparamset, power, mean_beta, n, version, design, N, t_half,
           repairs, censoring, estimable, n_conv, nonconv,
           mean_selected_half_life) |>
    arrange(match(design, c('OL', 'OL+BDC', 'CO', 'N-of-1')), N, t_half,
            censorparamset)

  grid_n_sum <- res |>
    filter(block == 'N') |>
    group_by(cell_id) |>
    summ_power() |>
    left_join(grid |> select(cell_id, design, c_bm, t_half, N, N80),
              by = 'cell_id') |>
    rename(N80_E9 = N80) |>
    arrange(c_bm, t_half, N)
  sum_n <- if (nrow(grid_n_sum) == 0) tibble() else
    grid_n_sum |>
      group_by(design, c_bm, t_half, N80_E9) |>
      group_modify(~ fit_n80(.x)) |>
      ungroup() |>
      mutate(ratio_to_E9 = N80_G / N80_E9)

  if (nrow(sum_p) > 0 && !CROSS) {
    write.csv(sum_p |>
                select(design, c_bm, t_half, N_label, N, N80_E9, capped,
                       n_reps, n_reps_requested, power, mcse, size,
                       mean_selected_half_life, nonconv, e9_power),
              file.path(OUT_DIR, sprintf('planning-%s.csv', tolower(label))),
              row.names = FALSE)
    cat(sprintf('\n%s, block P (power; size where c_bm = 0):\n', label))
    print(as.data.frame(sum_p |>
      select(design, c_bm, t_half, N_label, N, capped, n_reps, power, mcse,
             mean_selected_half_life, nonconv) |>
      mutate(across(c(power, mcse, mean_selected_half_life, nonconv),
                    ~ round(.x, 3)))), max = 4000)
  }
  if (nrow(sum_f) > 0 && !CROSS && any(grid$block == 'F')) {
    write.csv(sum_f, file.path(OUT_DIR, sprintf('fig4b-%s-cbm%.2f.csv',
                                                tolower(label), FIG_CBM)),
              row.names = FALSE, na = 'NA')
    cat(sprintf('\n%s, block F (power; head):\n', label))
    print(as.data.frame(head(sum_f |> mutate(across(where(is.double), ~ round(.x, 3))), 20)),
          row.names = FALSE)
  }
  if (nrow(sum_n) > 0) {
    write.csv(grid_n_sum |>
                select(design, c_bm, t_half, N80_E9, N, n_reps, n_p,
                       n_reject, power, mcse, nonconv),
              file.path(OUT_DIR, sprintf('n80-grid-%s.csv', tolower(label))),
              row.names = FALSE)
    write.csv(sum_n, file.path(OUT_DIR, sprintf('n80-%s.csv', tolower(label))),
              row.names = FALSE)
    cat(sprintf('\n%s, block N (N for 80 percent power):\n', label))
    print(as.data.frame(sum_n |> mutate(across(where(is.double), ~ round(.x, 3)))),
          row.names = FALSE)
  }
  list(P = sum_p, F = sum_f, N = sum_n, N_grid = grid_n_sum)
}

summaries <- set_names(map(SPECS, summarize_spec), LABELS)
saveRDS(if (length(SPECS) == 1) summaries[[1]] else summaries,
        file.path(OUT_DIR, 'summary.rds'))

## The cross-check compares G9 with paper 02's block A cell of the same
## parameter set (docs/71 puts the published one under published/).
if (CROSS && 'E7cr2' %in% SPECS) {
  cross_file <- if (param_set() == 'simcore2026')
    'analysis/data/carryover-sensitivity/referee-variants/block-A/summary.rds' else
    file.path('analysis/data/carryover-sensitivity', param_set(),
              'referee-variants/block-A/summary.rds')
  if (!file.exists(cross_file)) stop('no stored block A under ', param_set(),
                                    ': ', cross_file)
  stored <- readRDS(cross_file) |>
    filter(cell_id == 'alt-t0.5-N70-cbm0.45', method == 'E7cr2')
  sx <- summaries[['G9']]$P
  cat(sprintf('\nCross-check, Hybrid, c_bm 0.45, t1/2 0.5, N 70:\n  this driver: power %.3f (MCSE %.3f, %d replicates, non-convergence %.3f)\n  stored paper 02 block A: power %.3f (MCSE %.3f, %d replicates)\n',
              sx$power, sx$mcse, sx$n_p, sx$nonconv,
              stored$reject, stored$mcse, stored$n))
}

meta <- list(
  script = '19-g9-planning-and-fig4b.R', block = BLOCK, date = Sys.time(),
  params = param_meta(),
  reps_alt = REPS_ALT, reps_null = REPS_NULL, reps_fig = REPS_FIG,
  reps_n80 = REPS_N80, n80_mult = N80_MULT,
  fig_cbm = FIG_CBM, n_cap = N_CAP, smoke = SMOKE, cross_check = CROSS,
  spec = SPECS, guard_reps = GUARD_REPS, rho = RHO, plan_file = PLAN_FILE,
  elapsed_min_this_run = elapsed_min, r_version = R.version.string,
  git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                     error = function(e) NA_character_)
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
cat(sprintf('\nWritten: %s\n', OUT_DIR))
