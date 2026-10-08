## analysis/scripts/carryover-sensitivity/26-co-dbc-fix-paired.R
##
## Paired check of the two defects fixed in simulation-core.R on
## 2026-10-06 (docs/45, Section 8):
##
##   1. Dbc coded off-drug occasions before any exposure as
##      carryover_decay(0, h) = 1 instead of 0. In CO this makes the
##      placebo-first path fully "exposed" at every visit.
##   2. E7/E7cr2 (G4/G9) chose the analysis half-life by comparing REML
##      AICs across candidates whose fixed-effect designs differ (the
##      Dbc column), which REML likelihoods do not support.
##
## Each simulated trial is fitted three ways, so differences are paired:
##   original : defective coding, REML AIC   (the manuscript's results)
##   coding   : corrected coding, REML AIC   (fix 1 only)
##   fixed    : corrected coding, ML AIC, REML refit  (both fixes)
## E2/E2cr2 (G3/G8) depend only on the coding; E1/E1cr2 (G1/G6) on
## neither and are the reference.
##
## Cells: the CO cells of 19-run-tier1-hendrickson-g9.R (N in {35, 70};
## Panel A c_bm in {0, 0.30, 0.45} at t1/2 = 1; Panel B t1/2 in
## {0, 0.5} at c_bm = 0.45), plus the Hybrid reference cell (N = 70,
## c_bm = 0.45, t1/2 = 1) and its null as a control. Data generation
## is that of simulate_cell_s7() (rho = 0.7, no dropout).
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/26-co-dbc-fix-paired.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
n_reps <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 500L

seed <- 20261006L
set.seed(seed)
half_lives <- c(0.25, 0.5, 1.0, 2.0)

## The coding used by manuscript 02 before the fix.
exposure_weight_original <- function(Db, tsd, halflife) {
  ifelse(Db == 1, 1, carryover_decay(tsd, halflife))
}

fit_lme <- function(form, dl, method = 'REML') tryCatch(
  nlme::lme(form, random = ~1 | ptID,
            correlation = nlme::corCAR1(form = ~t | ptID),
            data = dl, method = method,
            control = nlme::lmeControl(
              opt = 'optim', maxIter = 200, msMaxIter = 200)),
  error = function(e) NULL)

f_db  <- Sx ~ bm + t + Db + bm:Db
f_dbc <- Sx ~ bm + t + Dbc + bm:Dbc

## Model-based and CR2 results for one fitted model.
extract <- function(fit, dl, target) {
  if (is.null(fit)) return(c(est = NA_real_, p_model = NA_real_, p_cr2 = NA_real_))
  tt <- summary(fit)$tTable
  r <- cr2_extract(fit, dl, target)
  c(est = tt[target[1], 'Value'], p_model = tt[target[1], 'p-value'],
    p_cr2 = if (is.null(r)) NA_real_ else r$p_value)
}

with_dbc <- function(dat_long, w) { dl <- dat_long; dl$Dbc <- w; dl }

fit_replicate <- function(dat_long, analysis_t1half) {
  Db <- dat_long$Db; tsd <- dat_long$tsd
  rows <- list()
  add <- function(spec, variant, v, h = NA_real_)
    rows[[length(rows) + 1]] <<- tibble(spec = spec, variant = variant,
                                        estimate = v[['est']],
                                        p_model = v[['p_model']],
                                        p_cr2 = v[['p_cr2']], best_t_half = h)

  ## G1/G6: binary indicator, unaffected by either defect.
  add('G1/G6', 'reference', extract(fit_lme(f_db, dat_long), dat_long, c('bm:Db', 'Db:bm')))

  ## G3/G8: fixed analysis half-life (the DGP value).
  for (v in c('original', 'fixed')) {
    w <- if (v == 'original') exposure_weight_original(Db, tsd, analysis_t1half)
         else exposure_weight(Db, tsd, analysis_t1half)
    dl <- with_dbc(dat_long, w)
    add('G3/G8', v, extract(fit_lme(f_dbc, dl), dl, c('bm:Dbc', 'Dbc:bm')))
  }

  ## G4/G9: half-life selected by AIC over the candidate grid.
  cand <- function(wfun, method) map(half_lives, function(hl) {
    dl <- with_dbc(dat_long, wfun(Db, tsd, hl))
    f <- fit_lme(f_dbc, dl, method)
    list(hl = hl, dl = dl, fit = f, aic = if (is.null(f)) NA_real_ else AIC(f))
  })
  pick <- function(cs) {
    a <- map_dbl(cs, ~ .x$aic)
    if (all(is.na(a))) NULL else cs[[which.min(a)]]
  }
  orig <- pick(cand(exposure_weight_original, 'REML'))
  reml_new <- cand(exposure_weight, 'REML')
  cod <- pick(reml_new)
  ml_best <- pick(cand(exposure_weight, 'ML'))
  ## REML refit of the ML-selected half-life is the REML candidate fit.
  fix <- if (is.null(ml_best)) NULL else reml_new[[match(ml_best$hl, half_lives)]]
  for (v in list(list('original', orig), list('coding', cod), list('fixed', fix))) {
    b <- v[[2]]
    if (is.null(b)) add('G4/G9', v[[1]], c(est = NA_real_, p_model = NA_real_, p_cr2 = NA_real_))
    else add('G4/G9', v[[1]], extract(b$fit, b$dl, c('bm:Dbc', 'Dbc:bm')), b$hl)
  }
  bind_rows(rows)
}

run_cell <- function(cell) {
  design_set <- build_design_set(cell$design)
  model_param <- list(
    N = cell$N, carryover_t1half = cell$t1half,
    carryover_form = 'exponential', weibull_shape = 1,
    dgp_architecture = 'mvn', c.bm = cell$c_bm,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.1)
  furrr::future_map_dfr(seq_len(n_reps), function(rep_i) {
    dat <- tryCatch(generate_data_multi_path(model_param, default_resp_param(),
                                             default_baseline_param(), design_set),
                    error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0) return(tibble(rep = rep_i))
    dat_long <- prepare_long_data(dat, design_set, carryover_t1half = cell$t1half)
    fit_replicate(dat_long, cell$t1half) |> mutate(rep = rep_i, .before = 1)
  }, .options = furrr::furrr_options(seed = TRUE))
}

co <- bind_rows(
  expand_grid(N = c(35, 70), c_bm = c(0, 0.30, 0.45)) |> mutate(t1half = 1.0, panel = 'A'),
  expand_grid(N = c(35, 70), t1half = c(0, 0.5)) |> mutate(c_bm = 0.45, panel = 'B')) |>
  mutate(design = 'CO')
hyb <- tibble(design = 'Hybrid', N = 70, c_bm = c(0.45, 0), t1half = 1.0, panel = 'control')
grid <- bind_rows(co, hyb)

cat(sprintf('%d cells, n_reps = %d\n', nrow(grid), n_reps))
plan(multicore, workers = max(1, parallel::detectCores() - 1))
t_start <- Sys.time()
results <- map_dfr(seq_len(nrow(grid)), function(i) {
  out <- run_cell(grid[i, ])
  cat(sprintf('%s cell %d/%d done\n', format(Sys.time(), '%H:%M:%S'), i, nrow(grid)))
  bind_cols(grid[rep(i, nrow(out)), ], out)
})
elapsed <- as.numeric(Sys.time() - t_start, units = 'secs')

summary_tab <- results |>
  filter(!is.na(spec)) |>
  group_by(design, N, c_bm, t1half, panel, spec, variant) |>
  summarise(n = n(), fails = sum(is.na(p_model)),
            mean_est = mean(estimate, na.rm = TRUE),
            emp_sd = sd(estimate, na.rm = TRUE),
            reject_model = mean(p_model < 0.05, na.rm = TRUE),
            reject_cr2 = mean(p_cr2 < 0.05, na.rm = TRUE),
            h_025 = mean(best_t_half == 0.25, na.rm = TRUE),
            h_05 = mean(best_t_half == 0.5, na.rm = TRUE),
            h_1 = mean(best_t_half == 1, na.rm = TRUE),
            h_2 = mean(best_t_half == 2, na.rm = TRUE),
            .groups = 'drop')

out_dir <- file.path(repo_root, 'analysis/data/carryover-sensitivity')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
meta <- list(script = '26-co-dbc-fix-paired.R', date = Sys.time(), seed = seed,
             n_reps = n_reps, elapsed_secs = elapsed, r_version = R.version.string,
             git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                                error = function(e) NA_character_))
saveRDS(list(results = results, summary = summary_tab, grid = grid, meta = meta),
        file.path(out_dir, sprintf('26-co-dbc-fix-paired-reps%d.rds', n_reps)))
readr::write_csv(summary_tab,
                 file.path(out_dir, sprintf('26-co-dbc-fix-paired-reps%d.csv', n_reps)))
print(as.data.frame(summary_tab |> select(design, N, c_bm, t1half, spec, variant,
                                          reject_model, reject_cr2, mean_est, h_2)),
      digits = 3)
cat(sprintf('elapsed %.1f min\n', elapsed / 60))
