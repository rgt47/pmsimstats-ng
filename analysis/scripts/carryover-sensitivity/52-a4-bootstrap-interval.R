## analysis/scripts/carryover-sensitivity/52-a4-bootstrap-interval.R
##
## A selection-aware interval for A4 (E7null_cr2, G10 in the companion
## papers) at the reference cell of paper 02 (Hybrid, N = 70,
## c_bm = 0.45, t1/2 = 1 week; job 27, panel A), on the same trials as
## the stored A4 fits (referee reports 2026-11-14 M2 and 2026-11-15,
## suggestion 2). The CR2 interval of A4 covers at 0.891 there (Table 7
## of the paper).
##
## For each Monte Carlo replicate r:
##   1. regenerate the trial from job 27's seed (24102600 + r) and fit
##      A4: half-life chosen by ML AIC over {0.25, 0.5, 1, 2} on the
##      model without the interaction, interaction model refitted by
##      REML, CR2 test. The estimate must reproduce the stored backfill
##      estimate to GUARD_TOL (pairing guard);
##   2. parametric bootstrap: draw B datasets from the fitted
##      interaction model at the selected half-life (fixed effects,
##      random intercept and corCAR1 residuals at their REML values;
##      the same patients, visits, biomarker values and paths) and
##      rerun the whole A4 procedure, selection included, on each;
##   3. store the B bootstrap estimates and standard errors.
## The intervals (percentile, basic, bootstrap-SE Wald and bootstrap-t)
## are formed and scored against theta = -c_bm * sigma_BR = -3.6 in
## 53-a4-bootstrap-summary.R (or the --summarise step below).
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/52-a4-bootstrap-interval.R \
##     [--reps-from A] [--reps-to B] [--boot B]
##     [--cores K] [--out-dir PATH] [--time] [--smoke]
##   --time runs one replicate with 10 resamples on one worker and
##   reports the cost; --smoke runs replicates 1 and 2 with 5 resamples.
## Completed replicates are checkpointed under <out-dir>/trials/, so a
## run can be split over hosts by replicate range and resumed.
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')
for (p in c('nlme', 'clubSandwich')) {
  if (!requireNamespace(p, quietly = TRUE)) stop('package missing: ', p)
}

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
SMOKE <- '--smoke' %in% args
TIME  <- '--time' %in% args
R_FROM <- as.integer(argv('--reps-from', '1'))
R_TO   <- as.integer(argv('--reps-to', if (SMOKE) '2' else '500'))
B      <- as.integer(argv('--boot', if (SMOKE) '5' else '199'))
CORES  <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
PARAMS <- Sys.getenv('PMSIM_PARAMS')
## PMSIM_PARAMS selects the set. Under published, the stored A4 fits
## come from that set's own job 27 (A4 was fitted in the run) and the
## output goes to a4-bootstrap-published/.
stopifnot(PARAMS %in% c('simcore2026', 'published'))
PUB <- PARAMS == 'published'
DATA_DIR <- 'analysis/data/carryover-sensitivity'
OUT_DIR <- argv('--out-dir', file.path(DATA_DIR, paste0(
  if (SMOKE) 'a4-bootstrap-smoke' else 'a4-bootstrap',
  if (PUB) '-published' else '')))
if (PUB) guard_out_dir(OUT_DIR)
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
GUARD_TOL <- 1e-6
CELL_ID <- 'panelA-designHybrid-N70-cbm0.45-t1.00'
CELL <- list(design = 'Hybrid', N = 70, c_bm = 0.45, t1half = 1,
             carryover_form = 'exponential', weibull_shape = 1,
             dgp_arch = 'mvn')
FULL <- Sx ~ bm + t + Dbc + bm:Dbc

stored <- (if (PUB)
  readRDS(file.path(DATA_DIR, 'published/tier1-g-vs-f-e0/all.rds'))
  else readRDS(file.path(DATA_DIR, 'g10-backfill/all.rds')) |>
    filter(source == 't1')) |>
  filter(cell_id == CELL_ID, method == 'E7null_cr2')
stopifnot(nrow(stored) == 500)

## ---------------- trial regeneration (as simulate_cell_g10) ----------------
regen <- function(r) {
  set.seed(24102600L + r)
  design_set <- build_design_set(CELL$design)
  mp <- list(N = CELL$N, carryover_t1half = CELL$t1half,
             carryover_form = CELL$carryover_form,
             weibull_shape = CELL$weibull_shape,
             dgp_architecture = CELL$dgp_arch, c.bm = CELL$c_bm,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2,
             c.cfct = 0.2, c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0,
             b.bm.tv.trend = 0, sd_slope = 0, c.bm.pb = 0)
  dat <- generate_data_multi_path(mp, default_resp_param(),
                                  default_baseline_param(), design_set,
                                  stop_on_nonpd = TRUE)
  prepare_long_data(dat, design_set, carryover_t1half = CELL$t1half,
                    carryover_form = CELL$carryover_form,
                    weibull_shape = CELL$weibull_shape,
                    analysis_t1half = CELL$t1half,
                    analysis_form = CELL$carryover_form,
                    analysis_shape = CELL$weibull_shape)
}

## A4 on one dataset: the stored procedure (select_half_life_aic with
## the null selection model, then CR2). Returns the fit for resampling.
a4 <- function(dl) {
  best <- select_half_life_aic(dl, selection_formula = E7_NULL_SELECTION)
  if (is.null(best)) return(NULL)
  r <- cr2_extract(best$fit, best$dl, c('bm:Dbc', 'Dbc:bm'))
  if (is.null(r) || !is.finite(r$std_error)) return(NULL)
  list(est = r$estimate, se = r$std_error, p = r$p_value, h = best$hl,
       fit = best$fit, dl = best$dl)
}

## One parametric resample from the fitted interaction model.
sim_full <- function(fit, dd) {
  mu <- drop(model.matrix(FULL[-2], dd) %*% nlme::fixef(fit))
  tau <- as.numeric(nlme::VarCorr(fit)['(Intercept)', 'StdDev'])
  phi <- as.numeric(coef(fit$modelStruct$corStruct,
                         unconstrained = FALSE))
  sig <- fit$sigma
  y <- mu
  for (idx in split(seq_len(nrow(dd)), dd$ptID)) {
    tt <- dd$t[idx]
    S <- tau^2 + sig^2 * phi^abs(outer(tt, tt, '-'))
    y[idx] <- mu[idx] + drop(crossprod(chol(S), rnorm(length(idx))))
  }
  out <- dd
  out$Sx <- y
  out
}

one_rep <- function(r, n_boot = B) {
  f <- file.path(ck_dir, sprintf('rep-%04d.rds', r))
  if (file.exists(f)) return(invisible(NULL))
  dl <- regen(r)
  obs <- a4(dl)
  st <- stored[stored$rep == r, ]
  guard <- !is.null(obs) && nrow(st) == 1 && isTRUE(st$converged) &&
    abs(obs$est - st$estimate) < GUARD_TOL
  res <- list(rep = r, guard = guard,
              est = if (is.null(obs)) NA_real_ else obs$est,
              se = if (is.null(obs)) NA_real_ else obs$se,
              p = if (is.null(obs)) NA_real_ else obs$p,
              h = if (is.null(obs)) NA_real_ else obs$h,
              boot = NULL)
  if (!is.null(obs)) {
    ## The resampling seed depends only on r, so a split run gives the
    ## same resamples as a single one.
    set.seed(label_seed(sprintf('a4boot-%s-rep-%d', CELL_ID, r)))
    bt <- map_dfr(seq_len(n_boot), function(b) {
      ds <- sim_full(obs$fit, obs$dl)
      ds$Dbc <- NULL
      x <- a4(ds)
      if (is.null(x)) return(tibble(b = b, est = NA_real_, se = NA_real_,
                                    h = NA_real_))
      tibble(b = b, est = x$est, se = x$se, h = x$h)
    })
    res$boot <- bt
  }
  saveRDS(res, f)
  invisible(NULL)
}

if (TIME) {
  t0 <- Sys.time()
  one_rep(R_FROM, n_boot = 10)
  secs <- as.numeric(Sys.time() - t0, units = 'secs')
  cat(sprintf('one replicate with 10 resamples: %.1f s; projected %.1f s per replicate at B = %d\n',
              secs, secs * (B + 1) / 11, B))
  unlink(file.path(ck_dir, sprintf('rep-%04d.rds', R_FROM)))
  quit(save = 'no')
}

t0 <- Sys.time()
invisible(mclapply(R_FROM:R_TO, one_rep, mc.cores = CORES,
                   mc.preschedule = FALSE))
cat(sprintf('replicates %d to %d, B = %d: %.1f min\n', R_FROM, R_TO, B,
            as.numeric(Sys.time() - t0, units = 'mins')))
saveRDS(list(params = PARAMS, cell_id = CELL_ID, boot = B,
             reps = c(R_FROM, R_TO), date = Sys.time()),
        file.path(OUT_DIR, sprintf('meta-%04d-%04d.rds', R_FROM, R_TO)))
