## analysis/scripts/carryover-sensitivity/diagnose-null-meanmisfit-referee.R
##
## Paper 10: the paired mean-misfit check of the 2026-10-28 referee
## report (Major comment 1), kept here for provenance; the referee's
## original was p10-meanmisfit-sim.R in a session scratchpad. Its
## output is output/diag-null-meanmisfit-referee-2026-10-28.rds
## (1,487 usable of 1,500 replicates, seed 20261028).
##
## Each replicate of the primary null cell is analyzed twice with M3
## and E1 on the same data: once on Sx and once on Sx minus its true
## path-by-visit mean. Removing the mean misfit moves kappa toward 1.
##
## DGP note. The run used c.cfct = 0.1 (the pre-correction constant,
## docs/61) with the current package's carryover, which decays from
## the last on-drug mean. It is therefore neither the 2026-05-20
## stored DGP (compounding carryover) nor the refined DGP
## (c.cfct = 0.2). The corrected counterpart is the 'oracle' arm of
## diagnose-null-type1.R.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-meanmisfit-referee.R [cfct] [reps]

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(furrr); library(nlme)
})
repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))
rp <- default_resp_param()
bp <- default_baseline_param()
ds <- build_design_set('Hybrid')
args <- commandArgs(trailingOnly = TRUE)
CFCT <- as.numeric(if (length(args)) args[1] else 0.1)
NREP <- as.integer(if (length(args) > 1) args[2] else 1500)
mp <- list(N = 70, c.bm = 0, carryover_t1half = 1,
  carryover_form = 'exponential', weibull_shape = 1,
  dgp_architecture = 'mvn', c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.2, c.cfct = CFCT)
mu_tab <- map_dfr(seq_along(ds), function(i) {
  s <- build_sigma_matrix(mp, rp, bp, ds[[i]], dgp_architecture = 'mvn')
  tpn <- ds[[i]]$timepoint_name
  lab <- s$labels
  m <- s$means
  mu <- sapply(tpn, function(tp) m[lab == 'BL'] -
    sum(m[lab %in% paste(tp, c('tv', 'pb', 'br'), sep = '.')]))
  tibble(path = i, timepoint_name = tpn, mu = mu)
})
fit1 <- function(dl, y) {
  dl$y <- dl[[y]]
  f <- tryCatch(lme(y ~ bm + t + Db + bm:Db, random = ~1 | ptID,
    correlation = corCAR1(form = ~t | ptID), data = dl,
    control = lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA, NA))
  tt <- summary(f)$tTable
  c(tt['bm:Db', 'Value'], tt['bm:Db', 'Std.Error'], tt['bm:Db', 'p-value'])
}
one <- function(i) {
  dat <- generate_data_multi_path(mp, rp, bp, ds)
  dl <- prepare_long_data(dat, ds, carryover_t1half = 1)
  dl$path_i <- as.integer(sub('_.*', '', dl$ptID))
  dl <- dl |> left_join(mu_tab, by = c(path_i = 'path', 'timepoint_name'))
  dl$Sx_c <- dl$Sx - dl$mu
  a <- fit1(dl, 'Sx')
  b <- fit1(dl, 'Sx_c')
  tibble(rep = i, est_raw = a[1], se_raw = a[2], p_raw = a[3],
         est_c = b[1], se_c = b[2], p_c = b[3])
}
plan(multicore, workers = max(1, parallelly::availableCores() - 1))
res <- future_map_dfr(seq_len(NREP), one,
  .options = furrr_options(seed = 20261028L)) |> na.omit()
cat(sprintf('cfct=%.1f reps=%d kappa raw=%.4f mean-corrected=%.4f\n',
  CFCT, nrow(res), mean(res$se_raw) / sd(res$est_raw),
  mean(res$se_c) / sd(res$est_c)))
saveRDS(res, file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output',
  sprintf('diag-null-meanmisfit-referee-cfct%.1f.rds', CFCT)))
