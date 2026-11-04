## analysis/scripts/carryover-sensitivity/42-check-param-sets.R
##
## Checks of the response/baseline parameter sets in simulation-core.R
## (docs/70, docs/71), run before any production job under the
## published set:
##
##   1. published_*() equal Hendrickson's values at commit 58b32a9
##      (data/extracted_rp.rda, data/extracted_bp.rda), and the
##      simcore2026 and fixture sets keep their stored values.
##   2. --params / PMSIM_PARAMS dispatch: default_resp_param() and
##      default_baseline_param() follow the active set.
##   3. Moments under the published set: the empirical SD of change
##      from baseline at two Hybrid timepoints matches the value
##      implied by the covariance matrix (within 2 percent).
##   4. prepare_long_data() puts bm on the SD scale (mean near 0, SD
##      near 1) under the published set.
##   5. E9 power under the published set reproduces paper 01's closed
##      form (planning-n80.csv) in the Hybrid, OL+BDC and crossover
##      designs (N = 70, c_bm = 0.25, t1/2 = 0.2 week), within 2 MCSE.
##   6. guard_out_dir() refuses a directory written under another set.
##
## Usage (repository root):
##   Rscript analysis/scripts/carryover-sensitivity/42-check-param-sets.R \
##     [--reps R]   (E9 replicates, default 1000; 0 skips check 5)
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
REPS <- if ('--reps' %in% args)
  as.integer(args[which(args == '--reps') + 1]) else 1000L

ok <- function(cond, msg) {
  cat(sprintf('[%s] %s\n', if (isTRUE(cond)) 'PASS' else 'FAIL', msg))
  if (!isTRUE(cond)) stop('check failed: ', msg)
}

## 1. stored values -------------------------------------------------------
rp <- published_resp_param(); bp <- published_baseline_param()
ok(identical(rp$cat, c('tv', 'pb', 'br')) &&
     isTRUE(all.equal(rp$max, c(6.50647, 6.50647, 10.98604),
                      tolerance = 1e-6)) &&
     isTRUE(all.equal(rp$disp, c(5, 5, 5))) &&
     isTRUE(all.equal(rp$rate, c(0.35, 0.35, 0.42))) &&
     isTRUE(all.equal(rp$sd, c(10, 10, 8))),
   'published response set equals 58b32a9 extracted_rp')
ok(identical(bp$cat, c('bm', 'BL')) &&
     isTRUE(all.equal(bp$m, c(124.32759, 83.06897), tolerance = 1e-6)) &&
     isTRUE(all.equal(bp$sd, c(15.36159, 18.48267), tolerance = 1e-6)),
   'published baseline set equals 58b32a9 extracted_bp')
ok(isTRUE(all.equal(simcore2026_resp_param()$sd, c(6, 2, 8))) &&
     isTRUE(all.equal(simcore2026_baseline_param()$m, c(0, 70))),
   'simcore2026 set unchanged')
ok(isTRUE(all.equal(fixture_resp_param()$sd, c(5, 2, 5))),
   'fixture set equals inst/tinytest')

## 2. dispatch ------------------------------------------------------------
old <- Sys.getenv('PMSIM_PARAMS')
Sys.setenv(PMSIM_PARAMS = 'simcore2026')
ok(isTRUE(all.equal(default_resp_param()$sd, c(6, 2, 8))),
   'dispatch follows PMSIM_PARAMS = simcore2026')
Sys.setenv(PMSIM_PARAMS = 'published')
ok(isTRUE(all.equal(default_resp_param()$sd, c(10, 10, 8))) &&
     isTRUE(all.equal(default_baseline_param()$sd[1], 15.36159,
                      tolerance = 1e-6)),
   'dispatch follows PMSIM_PARAMS = published')
ok(inherits(try(Sys.setenv(PMSIM_PARAMS = 'nonesuch') && param_set(),
                silent = TRUE), 'try-error'),
   'unknown set is refused')
Sys.setenv(PMSIM_PARAMS = 'published')

## 3. SD of change from baseline under the published set ----------------
mp <- list(N = 20000, carryover_t1half = 0.2,
           carryover_form = 'exponential', weibull_shape = 1,
           dgp_architecture = 'mvn', c.bm = 0.25,
           c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
           c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0, b.bm.tv.trend = 0,
           sd_slope = 0, c.bm.pb = 0)
ds1 <- build_design_set('Hybrid')[1]
set.seed(4201L)
dat <- generate_data_multi_path(mp, default_resp_param(),
                                default_baseline_param(), ds1)
b <- build_sigma_matrix(mp, default_resp_param(), default_baseline_param(),
                        ds1[[1]], stop_on_nonpd = TRUE)
tp <- ds1[[1]]$timepoint_name
for (nm in c(tp[2], tp[length(tp)])) {
  idx <- match(paste(nm, c('tv', 'pb', 'br'), sep = '.'), b$labels)
  sd_th <- sqrt(sum(b$sigma[idx, idx]))
  sd_emp <- sd(dat[[nm]] - dat$BL)
  ok(abs(sd_emp / sd_th - 1) < 0.02,
     sprintf('SD of change at %s: empirical %.2f, implied %.2f', nm,
             sd_emp, sd_th))
}

## 4. biomarker scale in the analysis data -------------------------------
dl <- prepare_long_data(dat, ds1, carryover_t1half = 0.2)
bm_pt <- dl |> distinct(ptID, bm) |> pull(bm)
ok(abs(mean(bm_pt)) < 0.05 && abs(sd(bm_pt) - 1) < 0.05,
   sprintf('analysis bm on the SD scale: mean %.3f, SD %.3f',
           mean(bm_pt), sd(bm_pt)))

## 5. E9 power against the closed form ----------------------------------
## One cell per design (c_bm 0.25, t1/2 0.2 week, N = 70), against the
## closed-form power of planning-n80.csv (17-rho-stage.R), which paper
## 01's planning table reports (referee report 2026-11-01, M1).
if (REPS > 0) {
  cf <- read.csv(here::here('analysis/data/quick-sim/carryover-closed-form',
                            'rho-stage/planning-n80.csv'))
  mp$N <- 70
  for (des in c('Hybrid', 'OL+BDC', 'CO')) {
    ds <- build_design_set(c(Hybrid = 'Hybrid', `OL+BDC` = 'OLBDC',
                             CO = 'CO')[[des]])
    one <- function(r) {
      set.seed(1000L + r)
      d <- generate_data_multi_path(mp, default_resp_param(),
                                    default_baseline_param(), ds)
      l <- prepare_long_data(d, ds, carryover_t1half = 0.2) |>
        filter(t > 0) |>
        group_by(ptID, bm) |>
        summarise(delta = mean(Sx[Db == 1]) - mean(Sx[Db == 0]),
                  .groups = 'drop')
      summary(lm(delta ~ bm, data = l))$coefficients['bm', 'Pr(>|t|)'] < 0.05
    }
    p <- mean(vapply(seq_len(REPS), one, logical(1)))
    p0 <- cf$power_N70[cf$design == des & cf$c_bm == 0.25 & cf$t_half == 0.2]
    stopifnot(length(p0) == 1)
    mcse <- sqrt(p0 * (1 - p0) / REPS)
    ok(abs(p - p0) < 2 * mcse,
       sprintf('E9 power, %s: %.3f vs closed form %.3f (2 MCSE = %.3f, %d reps)',
               des, p, p0, 2 * mcse, REPS))
  }
}

## 6. output-directory guard ----------------------------------------------
td <- tempfile('guard-')
dir.create(td)
writeLines('x', file.path(td, 'old-output.rds'))
ok(inherits(try(guard_out_dir(td), silent = TRUE), 'try-error'),
   'guard refuses an unmarked directory with files under published')
Sys.setenv(PMSIM_PARAMS = 'simcore2026')
guard_out_dir(td)
ok(readLines(file.path(td, 'PARAM_SET')) == 'simcore2026',
   'guard marks a directory with its set')
Sys.setenv(PMSIM_PARAMS = 'published')
ok(inherits(try(guard_out_dir(td), silent = TRUE), 'try-error'),
   'guard refuses a marked directory under another set')
unlink(td, recursive = TRUE)
if (nzchar(old)) Sys.setenv(PMSIM_PARAMS = old)

cat('\nAll parameter-set checks passed.\n')
