## analysis/scripts/prognostic-biomarker/01-run-bmt-block.R
##
## Paper 15, biomarker-by-time block (split plan 2026-10-29, Section
## 3.1) and the baseline-row contrast of paper 15 Section
## sec-reconcile. Data generation is that of block T of
## 36-run-referee-variants.R (simulation-core.R, unchanged; it is
## sourced, not edited). The analyses are defined here.
##
## Cells (Hybrid, t1/2 = 0.5 week, rho = 0.7, c1 = 0.2, fixed
## allocation):
##   T0-null  no trend, c_bm = 0,     N = 70    REPS_ALT
##   T0-alt   no trend, c_bm = 0.45,  N = 70    REPS_ALT
##   T3-null  trend 0.3, c_bm = 0,    N = 70    REPS_NULL
##   T3-alt   trend 0.3, c_bm = 0.45, N = 70    REPS_ALT
##   T3-null-N280  trend 0.3, c_bm = 0, N = 280 REPS_BIG
## The trend is b_tv_trend (block T): a natural-history mean shift of
## b * z_i * sigma_TV * g(w_t), g the TV mean course.
##
## Analyses, all fitted to the same trial in each replicate:
##   base  G8   (E2cr2): bm:Dbc at the true half-life, CR2
##         G9   (E7cr2): half-life by ML AIC on the interaction model
##         G10  (E7null_cr2): half-life by ML AIC on the null model
##         F1A0: F-family (De, bm:De), half-life by ML on the null
##               model over F_GRID
##   term  none | bmt (+ bm:t) | bmt2 (+ t2 + bm:t + bm:t2, t2 = t^2)
##         added to both the selection and the refit formulas
##   rows  bl = 0: first row is the first post-baseline visit
##               (the paper 02 analysis data)
##         bl = 1: plus a week-0 row holding BL, with Db = Dbc = De =
##               tsd = 0 (the second stress test's analysis data,
##               driver 11)
## The selected half-life of every selecting analysis is stored in h.
## Each replicate is seeded from its own label (label_seed()), as in
## driver 36.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/prognostic-biomarker/01-run-bmt-block.R \
##     --params published [--reps-null R] [--reps-alt R] \
##     [--reps-big R] [--cores K] [--out-dir PATH] [--smoke] [--time]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
REPS_NULL <- as.integer(argv('--reps-null', '5000'))
REPS_ALT  <- as.integer(argv('--reps-alt', '1000'))
REPS_BIG  <- as.integer(argv('--reps-big', '1000'))
CORES     <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
SMOKE     <- '--smoke' %in% args
TIME      <- '--time' %in% args
stopifnot(param_set() == 'published')
OUT_DIR <- argv('--out-dir',
  if (SMOKE) 'analysis/data/carryover-sensitivity/published/bmt-block-smoke'
  else 'analysis/data/carryover-sensitivity/published/bmt-block')
guard_out_dir(OUT_DIR)
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

T1HALF <- 0.5
G_GRID <- c(0.25, 0.5, 1, 2)

base_cell <- list(design = 'Hybrid', t1half = T1HALF,
                  carryover_form = 'exponential', weibull_shape = 1,
                  dgp_arch = 'mvn')
grid <- tribble(
  ~cell_id,         ~b_tv_trend, ~c_bm, ~N,  ~reps,
  'T0-null',        0,           0,     70,  REPS_ALT,
  'T0-alt',         0,           0.45,  70,  REPS_ALT,
  'T3-null',        0.3,         0,     70,  REPS_NULL,
  'T3-alt',         0.3,         0.45,  70,  REPS_ALT,
  'T3-null-N280',   0.3,         0,     280, REPS_BIG
)
if (SMOKE) grid$reps <- 2L

## ---------------- analyses ----------------
TERMS <- list(none = '', bmt = ' + bm:t', bmt2 = ' + t2 + bm:t + bm:t2')
G_FULL <- 'Sx ~ bm + t + Dbc + bm:Dbc'
G_NULL <- 'Sx ~ bm + t + Dbc'
F_FULL <- 'Sx ~ bm + t + Dbc + bm:Dbc + De + bm:De'
F_NULL <- 'Sx ~ bm + t + Dbc + De + bm:De'

lme_fit <- function(form, dl, method) tryCatch(
  nlme::lme(form, random = ~1 | ptID,
            correlation = nlme::corCAR1(form = ~t | ptID),
            data = dl, method = method,
            control = nlme::lmeControl(opt = 'optim', maxIter = 200,
                                       msMaxIter = 200)),
  error = function(e) NULL)

with_h <- function(dl, h) {
  dl$Dbc <- exposure_weight(dl$Db, dl$tsd, h)
  dl
}

## One analysis: choose the half-life on `grid_h` by the ML
## log-likelihood of `sel` (all candidates share one parameter count,
## so this is the AIC choice), refit `full` by REML there, test bm:Dbc
## with CR2. A single-value grid is a fixed half-life.
fit_one <- function(dl, full, sel, grid_h) {
  h <- grid_h
  if (length(grid_h) > 1) {
    ll <- map_dbl(grid_h, function(g) {
      f <- lme_fit(sel, with_h(dl, g), 'ML')
      if (is.null(f)) NA_real_ else as.numeric(logLik(f))
    })
    if (all(is.na(ll))) return(NULL)
    h <- grid_h[which.max(ll)]
  }
  dh <- with_h(dl, h)
  fit <- lme_fit(full, dh, 'REML')
  if (is.null(fit)) return(NULL)
  r <- cr2_extract(fit, dh, c('bm:Dbc', 'Dbc:bm'))
  if (is.null(r)) return(NULL)
  list(h = h, estimate = r$estimate, std_error = r$std_error,
       p_value = r$p_value)
}

ANALYSES <- list(
  G8   = list(full = G_FULL, sel = G_FULL, grid = T1HALF),
  G9   = list(full = G_FULL, sel = G_FULL, grid = G_GRID),
  G10  = list(full = G_FULL, sel = G_NULL, grid = G_GRID),
  F1A0 = list(full = F_FULL, sel = F_NULL, grid = F_GRID)
)

add_baseline_row <- function(dl, dat) {
  b <- dat |>
    transmute(ptID = paste(path, ptID, sep = '_'), path, bm, Sx = BL)
  b$bm <- dl$bm[match(b$ptID, dl$ptID)]
  b <- mutate(b, timepoint_name = 'BL', t = 0, Db = 0, Dbc = 0, L = 0,
              tsd = 0, De = 0)
  bind_rows(b[, names(dl)], dl) |> arrange(ptID, t)
}

fit_all <- function(dat, design_set) {
  dl0 <- prepare_long_data(dat, design_set, carryover_t1half = T1HALF)
  dl1 <- add_baseline_row(dl0, dat)
  out <- list()
  for (bl in 0:1) {
    dl <- if (bl == 0) dl0 else dl1
    dl$t2 <- dl$t^2
    for (tm in names(TERMS)) for (an in names(ANALYSES)) {
      a <- ANALYSES[[an]]
      r <- fit_one(dl, as.formula(paste0(a$full, TERMS[[tm]])),
                   as.formula(paste0(a$sel, TERMS[[tm]])), a$grid)
      out[[length(out) + 1]] <- tibble(
        analysis = an, term = tm, bl = bl,
        h = if (is.null(r)) NA_real_ else r$h,
        estimate = if (is.null(r)) NA_real_ else r$estimate,
        std_error = if (is.null(r)) NA_real_ else r$std_error,
        p_value = if (is.null(r)) NA_real_ else r$p_value,
        converged = !is.null(r))
    }
  }
  bind_rows(out)
}

simulate_one <- function(row) {
  cell <- modifyList(base_cell, list(N = row$N, c_bm = row$c_bm))
  design_set <- build_design_set(cell$design)
  mp <- list(N = cell$N, carryover_t1half = cell$t1half,
             carryover_form = cell$carryover_form,
             weibull_shape = cell$weibull_shape,
             dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2,
             c.cfct = 0.2, c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0,
             b.bm.tv.trend = row$b_tv_trend, sd_slope = 0, c.bm.pb = 0)
  dat <- generate_data_multi_path(mp, default_resp_param(),
                                  default_baseline_param(), design_set,
                                  stop_on_nonpd = TRUE)
  fit_all(dat, design_set)
}

## ---------------- --time ----------------
if (TIME) {
  for (k in c(1, 5)) {
    set.seed(label_seed(sprintf('%s-rep-%d', grid$cell_id[k], 1)))
    t0 <- Sys.time()
    res <- simulate_one(grid[k, ])
    secs <- as.numeric(Sys.time() - t0, units = 'secs')
    cat(sprintf('%s: %d rows, %d converged, %.1f s per replicate\n',
                grid$cell_id[k], nrow(res), sum(res$converged), secs))
  }
  print(as.data.frame(res))
  quit(save = 'no')
}

jobs <- grid |> mutate(rep = map(reps, seq_len)) |> unnest(rep)
cat(sprintf('bmt block: %d cells, %d (cell, rep) jobs, %d workers\n',
            nrow(grid), nrow(jobs), CORES))

one_job <- function(i) {
  j <- jobs[i, ]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', j$cell_id, j$rep))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(label_seed(sprintf('%s-rep-%d', j$cell_id, j$rep)))
  res <- tryCatch(simulate_one(j), error = function(e) NULL)
  if (!is.null(res)) {
    saveRDS(mutate(res, cell_id = j$cell_id, rep = j$rep, .before = 1), f)
  }
  invisible(NULL)
}

## Long jobs first, so the N = 280 cell does not finish last alone.
ord <- order(-jobs$N, seq_len(nrow(jobs)))
t0 <- Sys.time()
invisible(mclapply(ord, one_job, mc.cores = CORES,
                   mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

files <- list.files(ck_dir, full.names = TRUE)
res <- bind_rows(lapply(files, readRDS)) |>
  left_join(grid |> select(-reps), by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
cat(sprintf('%d of %d jobs present; %.1f min this run\n',
            length(unique(paste(res$cell_id, res$rep))), nrow(jobs),
            elapsed_min))

summary_tab <- res |>
  group_by(cell_id, b_tv_trend, c_bm, N, analysis, term, bl) |>
  summarise(n = sum(converged),
            reject = mean(p_value < 0.05, na.rm = TRUE),
            mcse = sqrt(reject * (1 - reject) / n),
            mean_est = mean(estimate, na.rm = TRUE),
            mcse_est = sd(estimate, na.rm = TRUE) / sqrt(n),
            emp_se = sd(estimate, na.rm = TRUE),
            mean_se = mean(std_error, na.rm = TRUE),
            kappa = mean_se / emp_se,
            nonconv = mean(!converged), .groups = 'drop')

## Paired contrasts on the same trials: each variant against the same
## analysis with no term and no baseline row.
ref <- res |> filter(term == 'none', bl == 0) |>
  select(cell_id, rep, analysis, p_ref = p_value, e_ref = estimate)
paired_tab <- res |>
  inner_join(ref, by = c('cell_id', 'rep', 'analysis')) |>
  mutate(d = (p_value < 0.05) - (p_ref < 0.05), de = estimate - e_ref) |>
  group_by(cell_id, analysis, term, bl) |>
  summarise(diff_reject = mean(d, na.rm = TRUE),
            se_diff_reject = sd(d, na.rm = TRUE) / sqrt(sum(!is.na(d))),
            diff_est = mean(de, na.rm = TRUE),
            se_diff_est = sd(de, na.rm = TRUE) / sqrt(sum(!is.na(de))),
            .groups = 'drop')

selection_tab <- res |>
  filter(analysis != 'G8', converged) |>
  count(cell_id, analysis, term, bl, h) |>
  group_by(cell_id, analysis, term, bl) |>
  mutate(share = n / sum(n)) |> ungroup()

saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))
saveRDS(paired_tab, file.path(OUT_DIR, 'paired.rds'))
saveRDS(selection_tab, file.path(OUT_DIR, 'selection.rds'))
saveRDS(list(script = '01-run-bmt-block.R', date = Sys.time(),
             params = param_meta(), grid = grid, terms = TERMS,
             analyses = ANALYSES, smoke = SMOKE,
             elapsed_min_this_run = elapsed_min,
             r_version = R.version.string,
             git_sha = tryCatch(system('git rev-parse --short HEAD',
                                       intern = TRUE),
                                error = function(e) NA_character_)),
        file.path(OUT_DIR, 'meta.rds'))
print(as.data.frame(summary_tab |>
  mutate(across(where(is.double), ~ round(.x, 3)))), max = 2000)
cat(sprintf('Written: %s\n', OUT_DIR))
