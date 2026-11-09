## analysis/scripts/carryover-sensitivity/44-g10-random-slope.R
##
## G10 with a random slope on time in the two patient-slope cells of
## referee-variants block S (paper 02 referee report 2026-11-08, m17).
## Each stored trial is regenerated with its block-S seed and fitted by
##   G10    E7null_cr2 with a random intercept, the pairing guard: its
##          estimate must reproduce the stored G10 of
##          g10-backfill/all.rds (source rvS) to GUARD_TOL;
##   G10rs  the same selection and refit with random = ~1 + t | ptID.
## The output also records, for each fit, whether the ML selection
## step failed or only the REML refit did (m18).
##
## Usage (repository root):
##   Rscript analysis/scripts/carryover-sensitivity/44-g10-random-slope.R \
##     [--max-reps R] [--cores K] [--time]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
MAX_REPS <- as.integer(argv('--max-reps', 1000))
CORES    <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
TIME     <- '--time' %in% args
DATA_DIR <- 'analysis/data/carryover-sensitivity'
OUT_DIR  <- file.path(DATA_DIR, if (TIME) 'g10-random-slope-time'
                      else 'g10-random-slope')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
GUARD_TOL <- 1e-6

## Block S cells, as in 41-run-g10-backfill.R.
rv_base <- list(design = 'Hybrid', N = 70, t1half = 0.5,
                carryover_form = 'exponential', weibull_shape = 1,
                dgp_arch = 'mvn', c_bm = 0, sd_slope = 0.3)
cells <- map(c(0, 0.45), ~ modifyList(rv_base, list(c_bm = .x)))
names(cells) <- sprintf('slope0.3-cbm%.2f-N70', c(0, 0.45))

## One G10 fit with a given random-effects structure; `stage` records
## where a failure happened.
fit_g10 <- function(dat_long, random, label) {
  out <- function(stage, est = NA_real_, se = NA_real_, p = NA_real_,
                  h = NA_real_) {
    tibble(method = label, h = h, estimate = est, std_error = se,
           p_value = p, converged = stage == 'ok', stage = stage)
  }
  fit_ml <- function(dl) tryCatch(
    nlme::lme(E7_NULL_SELECTION, random = random,
              correlation = nlme::corCAR1(form = ~t | ptID), data = dl,
              method = 'ML', control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  cands <- map(c(0.25, 0.5, 1.0, 2.0), function(hl) {
    dl <- dat_long
    dl$Dbc <- exposure_weight(dl$Db, dl$tsd, hl)
    f <- fit_ml(dl)
    list(hl = hl, dl = dl, aic = if (is.null(f)) NA_real_ else AIC(f))
  })
  aics <- map_dbl(cands, 'aic')
  if (all(is.na(aics))) return(out('selection'))
  best <- cands[[which.min(aics)]]
  fit <- tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc, random = random,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = best$dl, method = 'REML', control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(out('refit', h = best$hl))
  r <- cr2_extract(fit, best$dl, c('bm:Dbc', 'Dbc:bm'))
  if (is.null(r)) return(out('cr2', h = best$hl))
  out('ok', r$estimate, r$std_error, r$p_value, best$hl)
}

## Data generation as in simulate_cell_g10() (no dropout, no analysis
## overrides), so the trial matches the stored block-S trial.
one_trial <- function(cell) {
  design_set <- build_design_set(cell$design)
  model_param <- list(
    N = cell$N, carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
    c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0, b.bm.tv.trend = 0,
    sd_slope = cell$sd_slope, c.bm.pb = 0)
  dat <- generate_data_multi_path(model_param, default_resp_param(),
                                  default_baseline_param(), design_set,
                                  stop_on_nonpd = TRUE)
  prepare_long_data(dat, design_set, carryover_t1half = cell$t1half,
                    carryover_form = cell$carryover_form,
                    weibull_shape = cell$weibull_shape,
                    analysis_t1half = cell$t1half,
                    analysis_form = cell$carryover_form,
                    analysis_shape = cell$weibull_shape)
}

jobs <- expand.grid(cell_id = names(cells), rep = seq_len(MAX_REPS),
                    stringsAsFactors = FALSE)
one_job <- function(i) {
  id <- jobs$cell_id[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', id, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(label_seed(sprintf('%s-rep-%d', id, r)))
  dl <- tryCatch(one_trial(cells[[id]]), error = function(e) NULL)
  if (is.null(dl)) return(invisible(NULL))
  res <- bind_rows(fit_g10(dl, ~1 | ptID, 'G10'),
                   fit_g10(dl, ~1 + t | ptID, 'G10rs')) |>
    mutate(cell_id = id, rep = r, c_bm = cells[[id]]$c_bm, .before = 1)
  saveRDS(res, f)
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job, mc.cores = CORES,
                   mc.preschedule = FALSE))
cat(sprintf('%.1f min for %d jobs\n',
            as.numeric(Sys.time() - t0, units = 'mins'), nrow(jobs)))

res <- bind_rows(lapply(list.files(ck_dir, full.names = TRUE), readRDS))
if (!TIME) saveRDS(res, file.path(OUT_DIR, 'all.rds'))

## Pairing guard against the stored G10 of the backfill. The backfill
## itself matched its G9 reference in all but 11 of 55,000 replicates,
## so isolated numerical misses are recorded, not fatal; a pairing
## failure shows as a large share of misses.
bf <- readRDS(file.path(DATA_DIR, 'g10-backfill/all.rds')) |>
  filter(source == 'rvS', method == 'E7null_cr2') |>
  select(cell_id, rep, est_ref = estimate)
chk <- res |> filter(method == 'G10', converged) |>
  inner_join(bf, by = c('cell_id', 'rep')) |>
  mutate(gap = abs(estimate - est_ref))
n_miss <- sum(!(chk$gap < GUARD_TOL))
cat(sprintf('guard: %d paired fits, %d beyond %.0e, max |diff| %.2g\n',
            nrow(chk), n_miss, GUARD_TOL, max(chk$gap, na.rm = TRUE)))
if (!TIME) saveRDS(chk |> select(cell_id, rep, gap),
                   file.path(OUT_DIR, 'guard.rds'))
stopifnot(nrow(chk) > 0, n_miss <= 0.01 * nrow(chk))

summ <- res |>
  group_by(cell_id, c_bm, method) |>
  summarise(n = n(), reject = mean(p_value[converged] < 0.05),
            mcse = sqrt(reject * (1 - reject) / sum(converged)),
            fail = mean(!converged), fail_selection = mean(stage == 'selection'),
            fail_refit = mean(stage == 'refit'), fail_cr2 = mean(stage == 'cr2'),
            kappa = mean(std_error, na.rm = TRUE) / sd(estimate, na.rm = TRUE),
            .groups = 'drop')
print(as.data.frame(summ))
if (!TIME) saveRDS(summ, file.path(OUT_DIR, 'summary.rds'))
