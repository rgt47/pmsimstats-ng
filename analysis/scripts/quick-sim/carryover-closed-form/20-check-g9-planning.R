## 20-check-g9-planning.R
##
## Checks for 19-g9-planning-and-fig4b.R.
##   (default)     one replicate per design: fitted formula, selected
##                 half-life and the Dbc values on one path (the
##                 placebo-first path for CO); G9 on the open-label
##                 design; the censoring port against the 58b32a9
##                 censordata(); the stored censoring parameters.
##   --smoke       driver with 2 replicates per cell on 4 workers.
##   --preflight   driver pre-flight.
##   --time        driver timing and host projection.
##   --cross-check driver cross-check against paper 02 block A.
## The driver flags run the driver as a child process with --out-dir
## under OUT (default: the session scratchpad given by --out).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/20-check-g9-planning.R \
##     [--smoke] [--preflight] [--time] [--cross-check] [--out DIR]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(clubSandwich)
})
args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
DRIVER <- 'analysis/scripts/quick-sim/carryover-closed-form/19-g9-planning-and-fig4b.R'
OUT <- argv('--out', file.path(tempdir(), 'g9-check'))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

run_driver <- function(extra) {
  status <- system2('Rscript', c(DRIVER, extra))
  stopifnot(status == 0)
}
driver_flags <- c('--smoke', '--preflight', '--time', '--cross-check') %in% args
if (any(driver_flags)) {
  if ('--preflight' %in% args) run_driver(c('--preflight', '--out-dir', file.path(OUT, 'pf')))
  if ('--smoke' %in% args) run_driver(c('--smoke', '--cores', '4', '--out-dir', file.path(OUT, 'g9-smoke')))
  if ('--time' %in% args) run_driver(c('--time', '--out-dir', file.path(OUT, 'time')))
  if ('--cross-check' %in% args) run_driver(c('--cross-check', '--cores', '4', '--out-dir', file.path(OUT, 'cross')))
  quit(save = 'no')
}

source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

mp <- function(N, t_half, c_bm) list(
  N = N, carryover_t1half = t_half, carryover_form = 'exponential',
  weibull_shape = 1, dgp_architecture = 'mvn', c.bm = c_bm, c.tv = 0.7,
  c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2, c.bm.bl = 0,
  c.bm.tv = 0, b.bm.tv = 0, b.bm.tv.trend = 0, sd_slope = 0, c.bm.pb = 0)

cat('--- 1. Formula and Dbc, one replicate per design (t1/2 = 0.5, N = 70, c_bm = 0.25) ---\n')
for (nm in c('OL', 'OLBDC', 'CO', 'Hybrid')) {
  set.seed(label_seed(paste0('check-', nm)))
  ds <- build_design_set(nm)
  dat <- generate_data_multi_path(mp(70, 0.5, 0.25), default_resp_param(),
                                  default_baseline_param(), ds, stop_on_nonpd = TRUE)
  dl <- prepare_long_data(dat, ds, carryover_t1half = 0.5)
  best <- select_half_life_aic(dl)
  cat(sprintf('\n[%s] %d paths, %d long rows\n', nm, length(ds), nrow(dl)))
  if (is.null(best)) {
    cat('  select_half_life_aic returned NULL: no estimable fit (G9 not estimable)\n')
    cat('  Dbc at h = 0.5 over the design (constant => bm:Dbc aliased with bm):',
        paste(unique(exposure_weight(dl$Db, dl$tsd, 0.5)), collapse = ', '), '\n')
    next
  }
  cat('  formula:', deparse(formula(best$fit)), '\n')
  cat('  selected half-life:', best$hl, '\n')
  r <- cr2_extract(best$fit, best$dl, c('bm:Dbc', 'Dbc:bm'))
  cat(sprintf('  bm:Dbc estimate %.3f, CR2 SE %.3f, p %.4f\n',
              r$estimate, r$std_error, r$p_value))
  path_show <- if (nm == 'CO') 2L else 1L
  one <- best$dl |>
    filter(path == path_show) |>
    distinct(t, Db, tsd, Dbc) |>
    arrange(t)
  cat(sprintf('  Dbc on path %d%s:\n', path_show,
              if (nm == 'CO') ' (placebo first: 0 before first exposure)' else ''))
  print(as.data.frame(one), row.names = FALSE)
}

## Since 2026-10-29 censor_data() weights the biased term by
## (100 - D_)^eb1, the change in symptoms since baseline as the
## published Methods describe; 58b32a9 uses (100 + D_)^eb1. The 58b32a9
## function therefore reproduces the port when its D_ columns are
## negated (identical whenever no improvement exceeds 100 points).
cat('\n--- 2. Censoring port against the 58b32a9 censordata() ---\n')
suppressMessages({library(data.table); library(MASS)})
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
orig <- new.env()
sys.source(file.path(src, 'R', 'utilities.R'), envir = orig)
sys.source(file.path(src, 'R', 'censordata.R'), envir = orig)
pub <- new.env(); load(file.path(src, 'data', 'results_core.rda'), envir = pub)
cp <- pub$results_core$parameterselections$censorparams
ds <- build_design_set('Hybrid')
set.seed(1)
dat <- generate_data_multi_path(mp(70, 0.2, 0.25), default_resp_param(),
                                default_baseline_param(), ds)
td <- ds[[1]]
all_equal <- TRUE
for (k in seq_len(nrow(cp))) {
  set.seed(100 + k)
  a <- censor_data(dat, td, cp[k, ])
  set.seed(100 + k)
  td_dt <- data.table(timeptname = td$timepoint_name, t_wk = td$t_wk)
  dneg <- as.data.table(dat)
  dcols <- paste0('D_', td$timepoint_name)
  dneg[, (dcols) := lapply(.SD, `-`), .SDcols = dcols]
  b <- orig$censordata(dneg, td_dt, cp[k, ])
  same <- identical(is.na(as.matrix(a[td$timepoint_name])),
                    is.na(as.matrix(as.data.frame(b)[td$timepoint_name])))
  all_equal <- all_equal && same
  cat(sprintf('  %-15s NA pattern identical: %s (%d of %d cells censored)\n',
              cp$names[k], same, sum(is.na(a[td$timepoint_name])),
              prod(dim(a[td$timepoint_name]))))
}
cat('  Censoring port reproduces censordata() with D_ negated:', all_equal, '\n')

cat('\n--- 3. Stored censoring parameters against the driver constants ---\n')
drv <- readLines(DRIVER)
cat('  stored:', paste(cp$names, collapse = ', '), '\n')
cat('  beta0:', cp$beta0, ' beta1:', cp$beta1, ' eb1:', cp$eb1, '\n')
cat('  driver constants hard-coded: beta0 0.05 0.05 0.05 0.15; beta1 0.5 0.2 0.8 0.5; eb1 2\n')
stopifnot(all.equal(cp$beta0, c(0.05, 0.05, 0.05, 0.15)),
          all.equal(cp$beta1, c(0.5, 0.2, 0.8, 0.5)), all(cp$eb1 == 2),
          identical(as.character(cp$names),
                    c('balanced', 'more of flat', 'more of biased', 'high dropout')))
cat('  match: TRUE\n')
