## analysis/scripts/latent-class-mixture-application/13-check-referee-rerun.R
##
## Verification companion to 12-referee-rerun.R (referee report
## 2026-10-28, M1 and M11). Runs one null replicate at N = 35 on the
## Hybrid design with the driver's settings and prints:
##   1. the per-path allocation of N;
##   2. the fixed-effects formula lme_analysis() actually fits, and
##      whether the corCAR1 residual correlation was retained;
##   3. the Dbc, De, Db and tsd values on path A, visit by visit, so the
##      reader can confirm Dbc is 1 on drug, (1/2)^(tsd / 0.2) off drug
##      after first exposure, and 0 before exposure;
##   4. the parameter names of the hlme ng = 1 and ng = 2 fits, to show
##      the mixture arm carries the identical fixed-effects terms.
## It then inspects a finished driver output directory, if one is given,
## and reports that every analysis returned a value in the smoke run.
##
## Usage (from the repo root):
##   Rscript analysis/scripts/latent-class-mixture-application/13-check-referee-rerun.R \
##     [--out-dir PATH]

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(lcmm)
})
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
OUT_DIR <- argv('--out-dir', NA_character_)

RHO <- 0.7; C_CF <- 0.2; T_HALF_DGP <- 0.5; T_HALF_ANALYSIS <- 0.2
N_CHECK <- 35L

spec <- design_preset('Hybrid')
td <- buildtrialdesign(
  name_longform = 'Hybrid', name_shortform = 'Hybrid',
  timepoints = spec$timepoints, timeptnames = spec$timeptnames,
  expectancies = spec$expectancies, ondrug = spec$ondrug
)
design_set <- td$trialpaths
n_paths <- length(design_set)
n_path <- N_CHECK %/% n_paths +
  c(rep(1, N_CHECK %% n_paths), rep(0, n_paths - N_CHECK %% n_paths))
cat(sprintf('Design: Hybrid, %d paths; N = %d allocated as %s\n',
            n_paths, N_CHECK, paste(n_path, collapse = '/')))
cat('Schedule (path A):\n')
print(as.data.frame(design_set[[1]]), row.names = FALSE)

rp <- as.data.table(default_resp_param())
bp <- as.data.table(default_baseline_param())
set.seed(label_seed('null-N35-rep-1'))
dat <- rbindlist(lapply(seq_len(n_paths), function(g) {
  mp <- data.table(N = n_path[g], c.bm = 0, carryover_t1half = T_HALF_DGP,
                   c.tv = RHO, c.pb = RHO, c.br = RHO,
                   c.cf1t = C_CF, c.cfct = C_CF)
  d <- generateData(mp, rp, bp, design_set[[g]], empirical = FALSE,
                    makePositiveDefinite = TRUE, dgp_architecture = 'mvn')
  d[, path := g]
  d
}))
op <- list(useDE = TRUE, t_random_slope = FALSE, full_model_out = TRUE,
           carryover_t1half = T_HALF_ANALYSIS, simplecarryover = FALSE,
           carryover_scalefactor = 1)
prep <- lme_analysis(design_set, dat, op)

cat('\nFORMULA FIT BY lme_analysis:', paste(deparse(prep$form), collapse = ''), '\n')
cat('Residual correlation structure:',
    paste(class(prep$fit$modelStruct$corStruct), collapse = '/'), '\n')
cat(sprintf('Analysis half-life for Dbc: %.1f week\n', T_HALF_ANALYSIS))
cat('\nDbc on path A (participant 1), ordered by time:\n')
dbc <- prep$datamerged[ptID == 1, .(timeptnames, t, De, Db, tsd, Dbc)]
setorder(dbc, t)
print(as.data.frame(dbc), row.names = FALSE, digits = 6)
cat('\nlme tTable:\n')
print(round(summary(prep$fit)$tTable, 4))

df <- as.data.frame(prep$datamerged)
m1 <- hlme(fixed = prep$form, random = ~1, subject = 'ptID', ng = 1,
           data = df, verbose = FALSE, maxiter = 100)
m2 <- gridsearch(
  m = hlme(fixed = prep$form, mixture = ~Dbc, random = ~1, subject = 'ptID',
           ng = 2, classmb = ~bm, data = df, verbose = FALSE, maxiter = 100),
  rep = 5, maxiter = 30, minit = m1)
cat('\nhlme ng = 1 parameters:', paste(names(m1$best), collapse = ', '), '\n')
cat('hlme ng = 2 parameters:', paste(names(m2$best), collapse = ', '), '\n')
cat(sprintf('hlme conv (ng = 1, ng = 2): %d, %d; niter: %d, %d; LRT = %.3f\n',
            m1$conv, m2$conv, m1$niter, m2$niter,
            2 * (m2$loglik - m1$loglik)))

if (!is.na(OUT_DIR) && file.exists(file.path(OUT_DIR, 'all.rds'))) {
  all <- readRDS(file.path(OUT_DIR, 'all.rds'))
  summ <- readRDS(file.path(OUT_DIR, 'summary.rds'))
  cat(sprintf('\nOutput directory %s: %d rows in all.rds, %d cells, %d analyses\n',
              OUT_DIR, nrow(all), uniqueN(all$cell_id), uniqueN(all$analysis)))
  cat('Columns of all.rds:', paste(names(all), collapse = ', '), '\n')
  cat('Columns of summary.rds:', paste(names(summ), collapse = ', '), '\n')
  cat('\nValues returned per (cell, analysis) [finite estimate / rows]:\n')
  print(as.data.frame(all[, .(rows = .N, finite_estimate = sum(is.finite(estimate)),
                              finite_p = sum(is.finite(p_value)),
                              boot_lrt = sum(is.finite(boot_lrt_stat))),
                          by = .(cell_id, analysis)]), row.names = FALSE, max = 5000)
  cat('\nDistinct lme formulas in the run:\n')
  print(unique(all$lme_formula))
  cat('\nsummary.rds head:\n')
  print(head(as.data.frame(summ), 12))
}
