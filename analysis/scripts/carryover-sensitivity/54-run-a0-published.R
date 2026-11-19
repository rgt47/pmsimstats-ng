## analysis/scripts/carryover-sensitivity/54-run-a0-published.R
##
## Refit of paper 02's rung 0 (A0) as Hendrickson et al. published it.
##
## The stored E0 (fit_spec_e0() in simulation-core.R) differs from the
## published analysis in two respects. It adds the expectancy main
## effect De when expectancy varies, copying lme_analysis()'s default
## useDE = TRUE, whereas the publication driver
## (vignettes/Produce_Publication_Results_1_generate_data.Rmd at commit
## 58b32a9) passes useDE = FALSE. It also has no baseline (BL, t = 0)
## outcome row, which lme_analysis() prepends. This driver fits the
## published analysis, labeled E0p:
##
##   lmerTest::lmer(Sx ~ bm + Db + t + bm:Db + (1 | ptID)), Satterthwaite
##   df, with a BL row per participant (Sx = BL, t = 0, Db = 0), no De,
##   and the published exclusion rule (participants never on drug are
##   dropped only when Db never varies within a participant).
##
## Stage 1 covers source t1 only: the 30 cells of job 27 (script 27,
## tier1-g-vs-f-e0), 500 replicates each. Each (cell, replicate) is
## regenerated with job 27's seed (24102600 + r) by the same calls as
## simulate_cell_g10(), and fitted with E0 (stored formula, the pairing
## guard) and E0p. The guard requires the refitted E0 estimate to match
## the stored E0 estimate to GUARD_TOL in every replicate.
##
## --validate fits one trial per design (and one dropout trial) with
## E0p and with the vendored lme_analysis() of commit 58b32a9
## (analysis/data/quick-sim/hendrickson-58b32a9, useDE = FALSE), and
## requires equal p-values. The two codes may put the biomarker on
## different scales (under the published set they differ by the
## biomarker SD, 15.36), so the estimates must agree up to one constant
## factor: est_mine / est_orig is the same in every probe and equals 1
## or the biomarker SD of the active set.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/54-run-a0-published.R \
##     [--validate] [--smoke] [--time] [--max-reps R] [--cores K] \
##     [--out-dir PATH] [--params simcore2026|published]
## --params names the parameter set of the job-27 run to refit
## (default simcore2026). With published, the stored E0 comes from
## published/tier1-g-vs-f-e0 and the output goes to
## published/a0-published, which carries a PARAM_SET marker.
.args0 <- commandArgs(trailingOnly = TRUE)
PARAMS <- if ('--params' %in% .args0)
  .args0[which(.args0 == '--params') + 1] else 'simcore2026'
stopifnot(PARAMS %in% c('simcore2026', 'published'))
Sys.setenv(PMSIM_PARAMS = PARAMS)
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')
stopifnot(param_set() == PARAMS)

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
VALIDATE <- '--validate' %in% args
SMOKE    <- '--smoke' %in% args
TIME     <- '--time' %in% args
MAX_REPS <- as.integer(argv('--max-reps', NA))
CORES    <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
DATA_DIR <- 'analysis/data/carryover-sensitivity'
SET_DIR  <- if (PARAMS == 'published') file.path(DATA_DIR, 'published') else
  DATA_DIR
OUT_DIR  <- argv('--out-dir',
  file.path(SET_DIR, if (SMOKE) 'a0-published-smoke' else 'a0-published'))
ck_dir <- file.path(OUT_DIR, 'trials')
GUARD_TOL <- 1e-6

## ---------------- the published analysis ----------------
## dat: the wide trial from generate_data_multi_path() (with BL);
## dat_long: prepare_long_data() output, whose bm and ptID are reused.
add_bl_rows <- function(dat, dat_long) {
  pts <- dat_long |> distinct(ptID, path, bm)
  bl <- dat |>
    transmute(ptID = paste(path, ptID, sep = '_'), Sx = BL)
  bl_rows <- pts |> inner_join(bl, by = 'ptID') |>
    mutate(timepoint_name = 'BL', t = 0, Db = 0)
  stopifnot(nrow(bl_rows) == nrow(pts))
  bind_rows(bl_rows,
            dat_long |> select(ptID, path, bm, timepoint_name, Sx, t, Db))
}

## Data, fixed-effect formula and target of the published analysis.
e0p_setup <- function(dat, dat_long) {
  dd <- add_bl_rows(dat, dat_long)
  db_tab <- dd |> filter(t > 0) |> group_by(ptID) |>
    summarise(mean_db = mean(Db), .groups = 'drop')
  var_in_db <- any(db_tab$mean_db != 0 & db_tab$mean_db != 1)
  if (!var_in_db) {
    ever_on <- db_tab$ptID[db_tab$mean_db == 1]
    dd <- filter(dd, ptID %in% ever_on)
  }
  list(dd = dd,
       form = if (var_in_db) Sx ~ bm + Db + t + bm:Db else Sx ~ bm + t + bm:t,
       targets = if (var_in_db) c('bm:Db', 'Db:bm') else c('bm:t', 't:bm'))
}

na_row <- function(spec) tibble(spec = spec, estimate = NA_real_,
                                std_error = NA_real_, p_value = NA_real_,
                                converged = FALSE)

fit_spec_e0p <- function(dat, dat_long, s = e0p_setup(dat, dat_long)) {
  form <- update(s$form, . ~ . + (1 | ptID))
  fit <- tryCatch(suppressMessages(lmerTest::lmer(form, data = s$dd)),
                  error = function(e) NULL)
  if (is.null(fit)) return(na_row('E0p'))
  cc <- summary(fit)$coefficients
  tgt <- intersect(s$targets, rownames(cc))
  if (length(tgt) == 0) return(na_row('E0p'))
  tibble(spec = 'E0p', estimate = cc[tgt[1], 'Estimate'],
         std_error = cc[tgt[1], 'Std. Error'],
         p_value = cc[tgt[1], 'Pr(>|t|)'], converged = TRUE)
}

## Isolation chain from the published analysis to G1, one change per
## step, on the same data:
##   E0p    lmerTest::lmer, Satterthwaite df (the published analysis)
##   E0pn   the same data and terms, nlme::lme with containment df
##          (isolates software and df method)
##   E0pnc  E0pn plus corCAR1(form = ~ t | ptID) (isolates the residual
##          correlation)
##   E1     G1 as paper 02 fits it (fit_spec_s7): no BL row, all
##          patients kept (isolates the BL row and the exclusion rule)
## The nlme settings are those of 38-run-e0-isolation.R.
fit_e0p_nlme <- function(spec, s, corr) {
  ctl <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)
  fit <- tryCatch(
    if (corr) nlme::lme(s$form, random = ~ 1 | ptID, data = s$dd,
                        correlation = nlme::corCAR1(form = ~ t | ptID),
                        control = ctl)
    else nlme::lme(s$form, random = ~ 1 | ptID, data = s$dd, control = ctl),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row(spec))
  cc <- summary(fit)$tTable
  tgt <- intersect(s$targets, rownames(cc))
  if (length(tgt) == 0) return(na_row(spec))
  tibble(spec = spec, estimate = cc[tgt[1], 'Value'],
         std_error = cc[tgt[1], 'Std.Error'],
         p_value = cc[tgt[1], 'p-value'], converged = TRUE)
}

## ---------------- regeneration, as simulate_cell_g10() ----------------
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
  rate <- default_val(cell[['dropout_rate']], 0)
  if (rate > 0)
    dat <- apply_dropout(dat, design_set, rate = rate,
                         mechanism = default_val(cell[['dropout_mech']], 'MCAR'),
                         per = default_val(cell[['dropout_per']], 'occasion'))
  dat_long <- prepare_long_data(
    dat, design_set, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form, weibull_shape = cell$weibull_shape,
    analysis_t1half = default_val(cell[['analysis_t1half']], cell$t1half),
    analysis_form = default_val(cell[['analysis_form']], cell$carryover_form),
    analysis_shape = default_val(cell[['analysis_shape']], cell$weibull_shape))
  list(dat = dat, dat_long = dat_long, design_set = design_set)
}

METHODS <- c('E0', 'E0p', 'E0pn', 'E0pnc', 'E1')
fit_both <- function(tr) {
  if (is.null(tr))
    return(tibble(method = METHODS, estimate = NA_real_,
                  std_error = NA_real_, p_value = NA_real_,
                  converged = FALSE))
  s <- e0p_setup(tr$dat, tr$dat_long)
  keep <- function(x) select(x, spec, estimate, std_error, p_value, converged)
  bind_rows(keep(fit_spec_e0(tr$dat_long)),
            fit_spec_e0p(tr$dat, tr$dat_long, s),
            fit_e0p_nlme('E0pn', s, FALSE),
            fit_e0p_nlme('E0pnc', s, TRUE),
            keep(fit_spec_s7(tr$dat_long, 'E1'))) |>
    rename(method = spec)
}

## ---------------- job 27 grid, as script 41 copies it ----------------
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

## ---------------- --validate: against lme_analysis() at 58b32a9 ----------------
if (VALIDATE) {
  suppressMessages({library(data.table); library(lme4); library(lmerTest)})
  hd <- 'analysis/data/quick-sim/hendrickson-58b32a9/R'
  he <- new.env()
  for (f in c('utilities.R', 'lme_analysis.R')) sys.source(file.path(hd, f), he)
  to_orig <- function(tr) {
    tds <- lapply(tr$design_set, function(td)
      data.table(timeptnames = td$timepoint_name, t_wk = td$t_wk,
                 e = td$e, tod = td$tod, tsd = td$tsd, tpb = td$tpb))
    tp <- unique(unlist(lapply(tds, `[[`, 'timeptnames')))
    d <- as.data.table(tr$dat)
    d <- d[, c('ptID', 'path', 'bm', 'BL', intersect(tp, names(d))),
           with = FALSE]
    list(tds = tds, dat = d)
  }
  probes <- list(
    list(cell = cell_of(which(t1_grid$cell_id ==
           'panelA-designHybrid-N70-cbm0.30-t1.00')), r = 1L),
    list(cell = cell_of(which(t1_grid$cell_id ==
           'panelA-designCO-N70-cbm0.30-t1.00')), r = 1L),
    list(cell = cell_of(which(t1_grid$cell_id ==
           'panelA-designOLBDC-N35-cbm0.45-t1.00')), r = 2L),
    list(cell = modifyList(cell_of(1), list(design = 'Hybrid', N = 70,
           c_bm = 0.45, t1half = 0.5, dropout_rate = 0.3,
           dropout_mech = 'MAR', dropout_per = 'patient')), r = 3L))
  out <- map_dfr(probes, function(p) {
    set.seed(seed_for(p$r))
    tr <- gen_trial(p$cell)
    mine <- fit_spec_e0p(tr$dat, tr$dat_long)
    o <- to_orig(tr)
    theirs <- he$lme_analysis(o$tds, o$dat,
                              list(useDE = FALSE, t_random_slope = FALSE,
                                   full_model_out = FALSE))
    tibble(design = p$cell$design, N = p$cell$N,
           dropout = default_val(p$cell[['dropout_rate']], 0),
           est_mine = mine$estimate, est_orig = theirs$beta,
           p_mine = mine$p_value, p_orig = theirs$p)
  })
  bm_sd <- default_baseline_param() |> filter(cat == 'bm') |> pull(sd)
  out <- mutate(out, ratio = est_mine / est_orig,
                d_p = abs(p_mine - p_orig))
  print(as.data.frame(out), digits = 8)
  k <- out$ratio[1]
  ok <- all(out$d_p < 1e-6) && all(abs(out$ratio / k - 1) < 1e-6) &&
    (abs(k - 1) < 1e-6 || abs(k / bm_sd - 1) < 0.01)
  cat(sprintf('\nEstimate ratio %.6g; biomarker SD of this set %.6g\n',
              k, bm_sd))
  cat(if (ok) '\nVALIDATION PASSED\n' else '\nVALIDATION FAILED\n')
  quit(save = 'no', status = if (ok) 0 else 1)
}

## ---------------- --time ----------------
if (TIME) {
  probe <- t1_grid |> group_by(design, N) |> slice(1) |> ungroup()
  secs <- map_dbl(seq_len(nrow(probe)), function(i) {
    set.seed(seed_for(1L)); t0 <- Sys.time()
    for (k in 1:3) fit_both(gen_trial(as.list(probe[i, ])))
    as.numeric(Sys.time() - t0, units = 'secs') / 3
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
cat(sprintf('A0 as published, source t1: %d cells, %d (cell, rep) jobs\n',
            nrow(t1_grid), nrow(jobs)))

one_job <- function(k) {
  i <- jobs$i[k]; r <- jobs$rep[k]; cid <- t1_grid$cell_id[i]
  f <- file.path(ck_dir, sprintf('t1--%s-rep-%04d.rds', cid, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(seed_for(r))
  res <- tryCatch(fit_both(gen_trial(cell_of(i))), error = function(e) NULL)
  if (!is.null(res)) saveRDS(mutate(res, cell_id = cid, rep = r, .before = 1), f)
  invisible(NULL)
}
t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job, mc.cores = CORES,
                   mc.preschedule = FALSE))
cat(sprintf('Fits done in %.1f min\n',
            as.numeric(Sys.time() - t0, units = 'mins')))

## Only replicates 1..reps count, so a run with a smaller --max-reps
## ignores checkpoints left by a larger earlier run; an unreadable
## checkpoint (a worker stopped mid-write) is deleted and refitted on
## the next run.
files <- list.files(ck_dir, full.names = TRUE, pattern = '^t1--')
files <- files[as.integer(sub('.*-rep-(\\d+)\\.rds$', '\\1', files)) <= reps]
read_ck <- function(f) tryCatch(readRDS(f), error = function(e) {
  file.remove(f); NULL })
res <- bind_rows(lapply(files, read_ck)) |>
  left_join(t1_grid |> select(cell_id, panel, design, N, c_bm, t1half),
            by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))

## ---------------- guard: refitted E0 against stored E0 ----------------
stored <- readRDS(file.path(SET_DIR, 'tier1-g-vs-f-e0/all.rds')) |>
  filter(method == 'E0') |>
  select(cell_id, rep, est_stored = estimate, p_stored = p_value)
guard <- res |> filter(method == 'E0') |>
  select(cell_id, rep, est_new = estimate) |>
  inner_join(stored, by = c('cell_id', 'rep')) |>
  mutate(abs_diff = abs(est_new - est_stored),
         match = (is.na(est_new) & is.na(est_stored)) |
           (!is.na(abs_diff) & abs_diff < GUARD_TOL))
stopifnot(nrow(guard) == nrow(distinct(res, cell_id, rep)))
cat(sprintf('\nGuard: %d of %d replicates reproduce stored E0 (max diff %.2g)\n',
            sum(guard$match), nrow(guard), max(guard$abs_diff, na.rm = TRUE)))
## The refitted G1 (E1) closes the chain; it must also match its stored
## fit, so the chain ends on the stored A1.
e1_chk <- res |> filter(method == 'E1') |>
  select(cell_id, rep, est_new = estimate) |>
  inner_join(readRDS(file.path(SET_DIR, 'tier1-g-vs-f-e0/all.rds')) |>
               filter(method == 'E1') |>
               select(cell_id, rep, est_stored = estimate),
             by = c('cell_id', 'rep')) |>
  mutate(abs_diff = abs(est_new - est_stored))
guard <- guard |>
  left_join(e1_chk |> select(cell_id, rep, e1_diff = abs_diff),
            by = c('cell_id', 'rep'))
cat(sprintf('E1 check: max |refit - stored| %.2g over %d replicates\n',
            max(e1_chk$abs_diff, na.rm = TRUE), nrow(e1_chk)))
saveRDS(guard, file.path(OUT_DIR, 'guard.rds'))

## ---------------- summary ----------------
ok <- guard |> filter(match) |> select(cell_id, rep)
summ <- res |> semi_join(ok, by = c('cell_id', 'rep')) |>
  group_by(cell_id, panel, design, N, c_bm, t1half, method) |>
  summarise(n = sum(converged), reject = mean(p_value[converged] < 0.05),
            mcse = sqrt(reject * (1 - reject) / n),
            nonconv = sum(!converged),
            mean_est = mean(estimate[converged]),
            mean_se = mean(std_error[converged]),
            emp_se = sd(estimate[converged]), .groups = 'drop')
saveRDS(summ, file.path(OUT_DIR, 'summary.rds'))
pooled <- summ |> filter(c_bm == 0) |> group_by(method) |>
  summarise(reject = sum(reject * n) / sum(n), n = sum(n), .groups = 'drop')
by_design <- summ |> filter(c_bm == 0) |> group_by(method, design) |>
  summarise(reject = sum(reject * n) / sum(n), n = sum(n), .groups = 'drop')
cat('\nNull rejection, pooled over designs:\n'); print(as.data.frame(pooled))
cat('\nNull rejection by design:\n'); print(as.data.frame(by_design))
saveRDS(list(date = Sys.time(), param_set = param_set(), reps = reps,
             guard_tol = GUARD_TOL, guard_pass = all(guard$match),
             r_version = R.version.string,
             lme4 = as.character(packageVersion('lme4')),
             lmerTest = as.character(packageVersion('lmerTest'))),
        file.path(OUT_DIR, 'meta.rds'))
