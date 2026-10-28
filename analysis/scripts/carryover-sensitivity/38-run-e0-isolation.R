## analysis/scripts/carryover-sensitivity/38-run-e0-isolation.R
##
## Isolation of the first rung of paper 02's analysis ladder. E0 (the
## original published analysis) and G1 (E1) differ in three respects at
## once: software and degrees-of-freedom method (lme4 with lmerTest
## Satterthwaite versus nlme with containment df), the residual
## correlation (none versus corCAR1), and the model terms and sample
## (expectancy De and exclusion of never-on-drug patients versus
## neither). This driver fits a chain that changes one thing at a time,
## all to the same simulated trial each replicate:
##
##   E0    fit_spec_e0(): lmerTest::lmer, random intercept, De when
##         expectancy varies, never-on-drug patients dropped, binary Db,
##         bm:Db, Satterthwaite df
##   E0n   E0's formula and exclusion rule, nlme::lme, no residual
##         correlation (isolates software and df method)
##   E0nc  E0n plus corCAR1(form = ~ t | ptID) (isolates the residual
##         correlation)
##   E1nc  G1's formula (no De, all patients, binary Db) with nlme::lme
##         and no residual correlation (compare E0n against E1nc to
##         isolate the expectancy term and patient dropping; compare
##         E1nc against E1 to isolate corCAR1 under G1's formula)
##   E1    fit_spec_s7(., 'E1'): G1 as paper 02 fits it, model-based
##
## The tested coefficient is the biomarker-by-drug interaction (bm:Db;
## bm:t if Db never varies within a patient, as fit_spec_e0 does) with
## the p-value from each fit's own reference distribution. Where the
## design's expectancy is constant in every path (the CO preset) E0, E0n
## and E0nc contain no De; the per-row flag `de_used` records this.
## `df` is NA for E0 and E1 because fit_spec_e0() and fit_spec_s7() do
## not return it.
##
## Grid: the 30 cells of 27-run-tier1-g-vs-f-e0.R (Panel A: c_bm in
## {0, 0.30, 0.45} at t1half = 1; Panel B: t1half in {0, 0.5} at c_bm =
## 0.45; designs CO, Hybrid, OLBDC; N in {35, 70}; exponential decay),
## same cell ids. DGP as simulation-core.R sets it (separable AR(1)
## rho = 0.7, c1 = 0.2, anchored carryover, stop on non-PD, fixed path
## allocation); pmsimstats.allow_repair is not set.
##
## Seeds: each (cell, replicate) seeds from label_seed() of its own
## label, so no cell shares a stream with another cell or with job 27.
## Comparisons across methods within a cell are paired.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/38-run-e0-isolation.R \
##     [--reps R] [--cores K] [--out-dir PATH] \
##     [--preflight] [--smoke] [--time]
##   --preflight checks every cell's covariance matrix and exits;
##   --smoke runs two replicates of every cell; --time runs 5 replicates
##   of the Hybrid N = 70 cell on one worker and reports the cost.
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

METHODS <- c('E0', 'E0n', 'E0nc', 'E1nc', 'E1')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
REPS      <- as.integer(argv('--reps', '1000'))
CORES     <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
PREFLIGHT <- '--preflight' %in% args
SMOKE     <- '--smoke' %in% args
TIME      <- '--time' %in% args
OUT_DIR <- argv('--out-dir',
  if (SMOKE) 'analysis/data/carryover-sensitivity/e0-isolation-smoke'
  else 'analysis/data/carryover-sensitivity/e0-isolation')

grid <- bind_rows(
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
stopifnot(nrow(grid) == 30, !anyDuplicated(grid$cell_id))

## The data-generating calls mirror simulate_cell_referee().
draw_trial <- function(cell) {
  design_set <- build_design_set(cell$design)
  rho <- default_val(cell[['rho']], 0.7)
  model_param <- list(
    N = cell$N, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
    c.tv = rho, c.pb = rho, c.br = rho, c.cf1t = 0.2, c.cfct = 0.2
  )
  dat <- generate_data_multi_path(model_param, default_resp_param(),
                                  default_baseline_param(), design_set,
                                  stop_on_nonpd = TRUE)
  prepare_long_data(dat, design_set, carryover_t1half = cell$t1half,
                    carryover_form = cell$carryover_form,
                    weibull_shape = cell$weibull_shape)
}

## Exclusion rule and formula of fit_spec_e0(), reproduced so that E0n
## and E0nc use exactly the data and terms E0 uses.
e0_setup <- function(dat_long) {
  dd <- dat_long
  db_tab <- dd |> filter(t > 0) |> group_by(ptID) |>
    summarise(mean_db = mean(Db), .groups = 'drop')
  var_in_db <- any(db_tab$mean_db != 0 & db_tab$mean_db != 1)
  if (!var_in_db) {
    dd <- filter(dd, ptID %in% db_tab$ptID[db_tab$mean_db == 1])
  }
  var_in_exp <- length(unique(dd$De[dd$t > 0])) > 1
  form <- if (var_in_db) {
    if (var_in_exp) Sx ~ bm + De + Db + t + bm:Db else Sx ~ bm + Db + t + bm:Db
  } else {
    if (var_in_exp) Sx ~ bm + De + t + bm:t else Sx ~ bm + t + bm:t
  }
  list(dd = dd, form = form, var_in_exp = var_in_exp,
       targets = if (var_in_db) c('bm:Db', 'Db:bm') else c('bm:t', 't:bm'))
}

row_na <- function(method, de_used = NA, n_pat = NA_integer_) {
  tibble(method = method, estimate = NA_real_, std_error = NA_real_,
         statistic = NA_real_, df = NA_real_, p_value = NA_real_,
         converged = FALSE, de_used = de_used, n_pat = n_pat)
}

fit_nlme_row <- function(method, form, dd, targets, corr, de_used) {
  n_pat <- dplyr::n_distinct(dd$ptID)
  fit <- tryCatch(
    if (corr) {
      nlme::lme(form, random = ~ 1 | ptID,
                correlation = nlme::corCAR1(form = ~ t | ptID), data = dd,
                control = nlme::lmeControl(opt = 'optim', maxIter = 200,
                                           msMaxIter = 200))
    } else {
      nlme::lme(form, random = ~ 1 | ptID, data = dd,
                control = nlme::lmeControl(opt = 'optim', maxIter = 200,
                                           msMaxIter = 200))
    },
    error = function(e) NULL)
  if (is.null(fit)) return(row_na(method, de_used, n_pat))
  cc <- summary(fit)$tTable
  tgt <- intersect(targets, rownames(cc))
  if (length(tgt) == 0) return(row_na(method, de_used, n_pat))
  tgt <- tgt[1]
  tibble(method = method, estimate = cc[tgt, 'Value'],
         std_error = cc[tgt, 'Std.Error'], statistic = cc[tgt, 't-value'],
         df = cc[tgt, 'DF'], p_value = cc[tgt, 'p-value'], converged = TRUE,
         de_used = de_used, n_pat = n_pat)
}

E1_FORM <- Sx ~ bm + t + Db + bm:Db

fit_chain <- function(dat_long) {
  s <- e0_setup(dat_long)
  e0 <- fit_spec_e0(dat_long)
  e0_row <- tibble(method = 'E0', estimate = e0$estimate,
                   std_error = e0$std_error,
                   statistic = e0$estimate / e0$std_error, df = NA_real_,
                   p_value = e0$p_value, converged = e0$converged,
                   de_used = s$var_in_exp,
                   n_pat = dplyr::n_distinct(s$dd$ptID))
  e1 <- fit_spec_s7(dat_long, 'E1')
  e1_row <- tibble(method = 'E1', estimate = e1$estimate,
                   std_error = e1$std_error,
                   statistic = e1$estimate / e1$std_error, df = NA_real_,
                   p_value = e1$p_value, converged = e1$converged,
                   de_used = FALSE,
                   n_pat = dplyr::n_distinct(dat_long$ptID))
  bind_rows(
    e0_row,
    fit_nlme_row('E0n', s$form, s$dd, s$targets, FALSE, s$var_in_exp),
    fit_nlme_row('E0nc', s$form, s$dd, s$targets, TRUE, s$var_in_exp),
    fit_nlme_row('E1nc', E1_FORM, dat_long, c('bm:Db', 'Db:bm'), FALSE, FALSE),
    e1_row
  )
}

one_trial <- function(cell_list, cell_id, r) {
  set.seed(label_seed(sprintf('%s-rep-%d', cell_id, r)))
  dat_long <- tryCatch(draw_trial(cell_list), error = function(e) NULL)
  if (is.null(dat_long) || nrow(dat_long) == 0) {
    return(map_dfr(METHODS, row_na) |> mutate(rep = r, cell_id = cell_id,
                                              .before = 1))
  }
  fit_chain(dat_long) |> mutate(rep = r, cell_id = cell_id, .before = 1)
}

if (!exists('E0_ISOLATION_LIB_ONLY')) {
  if (SMOKE) REPS <- 2L

  ## ---- pre-flight: every cell must be PD on every path ----
  preflight_cell <- function(cell) {
    design_set <- build_design_set(cell$design)
    rho <- default_val(cell[['rho']], 0.7)
    mp <- list(N = cell$N, carryover_t1half = cell$t1half,
               carryover_form = cell$carryover_form,
               weibull_shape = cell$weibull_shape,
               dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
               c.tv = rho, c.pb = rho, c.br = rho, c.cf1t = 0.2,
               c.cfct = 0.2)
    map_dfr(seq_along(design_set), function(i) {
      b <- build_sigma_matrix(mp, default_resp_param(),
                              default_baseline_param(), design_set[[i]],
                              stop_on_nonpd = FALSE)
      raw <- outer(b$standard_deviations, b$standard_deviations) *
        b$correlations
      pos <- b$standard_deviations > 0
      ev <- eigen(raw[pos, pos], symmetric = TRUE, only.values = TRUE)$values
      tibble(path = i, min_eigenvalue = min(ev), feasible = min(ev) > 0)
    })
  }
  pre <- grid |>
    mutate(pf = map(seq_len(n()), ~ preflight_cell(as.list(grid[.x, ])))) |>
    select(cell_id, pf) |>
    unnest(pf) |>
    group_by(cell_id) |>
    summarise(min_eigenvalue = min(min_eigenvalue),
              feasible = all(feasible), .groups = 'drop')
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(pre, file.path(OUT_DIR, 'preflight.rds'))
  cat(sprintf('Pre-flight: %d cells, %d feasible, min eigenvalue %.4f\n',
              nrow(pre), sum(pre$feasible), min(pre$min_eigenvalue)))
  if (!all(pre$feasible)) {
    print(as.data.frame(filter(pre, !feasible)), row.names = FALSE)
    stop('Infeasible cells present; adjust the grid before running.')
  }
  if (PREFLIGHT) {
    print(as.data.frame(mutate(pre, min_eigenvalue = round(min_eigenvalue, 4))),
          row.names = FALSE)
    quit(save = 'no')
  }

  ## ---- --time ----
  if (TIME) {
    k <- match('panelA-designHybrid-N70-cbm0.45-t1.00', grid$cell_id)
    cell <- as.list(grid[k, ])
    t0 <- Sys.time()
    res <- map_dfr(1:5, ~ one_trial(cell, grid$cell_id[k], .x))
    secs <- as.numeric(Sys.time() - t0, units = 'secs')
    cat(sprintf('%s: 5 replicates, %d method rows, %.2f s per replicate\n',
                grid$cell_id[k], nrow(res), secs / 5))
    quit(save = 'no')
  }

  ck_dir <- file.path(OUT_DIR, 'trials')
  dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf('E0 isolation: %d cells, %d reps, %d (cell, rep) jobs\n',
              nrow(grid), REPS, nrow(grid) * REPS))

  jobs <- grid |> select(cell_id) |>
    mutate(rep = list(seq_len(REPS))) |> unnest(rep)

  one_job <- function(i) {
    cell_id <- jobs$cell_id[i]; r <- jobs$rep[i]
    f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
    if (file.exists(f)) return(invisible(NULL))
    cell <- as.list(grid[match(cell_id, grid$cell_id), ])
    res <- tryCatch(one_trial(cell, cell_id, r), error = function(e) NULL)
    if (!is.null(res)) saveRDS(res, f)
    invisible(NULL)
  }

  t0 <- Sys.time()
  invisible(mclapply(seq_len(nrow(jobs)), one_job,
                     mc.cores = CORES, mc.preschedule = FALSE))
  elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

  files <- list.files(ck_dir, full.names = TRUE)
  files <- files[sub('-rep-\\d+\\.rds$', '', basename(files)) %in%
                   grid$cell_id]
  cell_cov <- grid |> select(cell_id, design, N, c_bm, t1half)
  res <- bind_rows(lapply(files, readRDS)) |>
    left_join(cell_cov, by = 'cell_id') |>
    select(cell_id, rep, method, estimate, std_error, statistic, df,
           p_value, converged, design, N, c_bm, t1half, de_used, n_pat)
  saveRDS(res, file.path(OUT_DIR, 'all.rds'))
  n_done <- length(unique(paste(res$cell_id, res$rep)))
  cat(sprintf('\n%d of %d (cell, rep) jobs present in %.1f min (this run)\n',
              n_done, nrow(jobs), elapsed_min))

  summary_tab <- res |>
    group_by(cell_id, design, N, c_bm, t1half, method) |>
    summarise(n = sum(!is.na(p_value)),
              reject = mean(p_value < 0.05, na.rm = TRUE),
              mcse = sqrt(reject * (1 - reject) / n),
              nonconv = mean(!converged),
              mean_est = mean(estimate, na.rm = TRUE),
              mean_se = mean(std_error, na.rm = TRUE),
              emp_sd = sd(estimate, na.rm = TRUE),
              de_used = any(de_used, na.rm = TRUE),
              .groups = 'drop') |>
    mutate(kappa = mean_se / emp_sd,
           method = factor(method, levels = METHODS)) |>
    arrange(cell_id, method)
  cat('\nRejection rate per (cell, method); kappa = mean SE / empirical SD:\n')
  print(as.data.frame(summary_tab |>
    select(cell_id, method, reject, mcse, n, nonconv, kappa) |>
    mutate(across(c(reject, mcse, nonconv, kappa), ~ round(.x, 3)))),
    max = 4000)
  saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))

  meta <- list(
    script = '38-run-e0-isolation.R', date = Sys.time(), reps = REPS,
    smoke = SMOKE, methods = METHODS, elapsed_min_this_run = elapsed_min,
    note_de = paste('de_used is FALSE for E0, E0n and E0nc in the CO',
                    'preset, where expectancy is constant in every path.'),
    note_df = 'df is NA for E0 and E1 (fit_spec_e0 and fit_spec_s7 omit it).',
    r_version = R.version.string,
    git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                       error = function(e) NA_character_)
  )
  saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
  cat(sprintf('\nWritten: %s\n', OUT_DIR))
}
