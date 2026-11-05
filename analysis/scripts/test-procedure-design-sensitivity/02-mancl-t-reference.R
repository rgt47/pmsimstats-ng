## 02-mancl-t-reference.R
##
## Referee report 2026-10-28, M2: the stored Mancl-DeRouen arm used a
## normal reference for the Wald statistic. This script recomputes the
## arm with a t reference on K - p degrees of freedom (K = N clusters,
## one per participant; p = 5 fixed-effect parameters: intercept, bm,
## t, Dbc, bm:Dbc) from the stored per-replicate beta and se. No model
## is refitted. It also runs an exact binomial test of every null cell
## of every arm against 0.05 and tabulates the calibration ratio kappa
## (mean model SE over the empirical SD of the estimates) for every
## arm and cell.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/test-procedure-design-sensitivity/02-mancl-t-reference.R \
##     [--prefix PREFIX | --params NAME]
##
## The prefix is the driver's OUT_PREFIX. It comes from `--prefix`,
## else the OUT_PREFIX environment variable, else the driver's default
## for the parameter set named by `--params` or PMSIM_PARAMS (default
## 'fixture'): analysis/data/quick-sim/08-test under fixture, the
## stored 2026-10-08 run; analysis/data/quick-sim/08-test-<set>/08-test
## under any other set.
##
## Inputs:  <PREFIX>-replicates.rds
##          <PREFIX>-betatrue.txt
## Outputs: <PREFIX>-mancl-tref.txt
##          <PREFIX>-null-binomial.txt
##          <PREFIX>-kappa.txt
##
## The production driver (01-study1-test-procedure.R) now uses the t
## reference natively; this script exists so that the manuscript's
## corrected numbers trace to the stored 2026-10-08 run without a rerun.

suppressPackageStartupMessages(library(data.table))

arg_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == flag)
  if (length(i) && i[1] < length(args)) args[i[1] + 1] else NA_character_
}

pset <- arg_value('--params')
if (is.na(pset)) pset <- Sys.getenv('PMSIM_PARAMS', 'fixture')
default_prefix <- if (pset == 'fixture') 'analysis/data/quick-sim/08-test' else
  sprintf('analysis/data/quick-sim/08-test-%s/08-test', pset)
prefix <- arg_value('--prefix')
if (is.na(prefix)) prefix <- Sys.getenv('OUT_PREFIX', default_prefix)
cat(sprintf('Prefix: %s\n', prefix))

## A run of the converted driver records its set; the stored run
## predates the metadata file.
meta_file <- paste0(prefix, '-meta.rds')
if (file.exists(meta_file)) {
  cat(sprintf('Parameter set of the run: %s\n',
              readRDS(meta_file)$params$param_set))
}

p_fixed <- 5L
rep_file <- paste0(prefix, '-replicates.rds')
bt_file  <- paste0(prefix, '-betatrue.txt')

r <- readRDS(rep_file)
bt_tab <- read.delim(bt_file, header = FALSE)
beta_true <- as.numeric(bt_tab$V2[bt_tab$V1 == 'beta_true_mc'])
stopifnot(is.finite(beta_true))

md <- r[procedure == 'gee_mancl' & converged & is.finite(beta) &
          is.finite(se)]
md[, df := N - p_fixed]
## The z reference is recomputed rather than read from p_interaction:
## the converted driver stores the t-reference p-value there, so the
## stored column is the z one only for the 2026-10-08 run.
md[, p_z := 2 * pnorm(-abs(beta / se))]
md[, lo_z := beta - qnorm(0.975) * se]
md[, hi_z := beta + qnorm(0.975) * se]
md[, p_t := 2 * pt(-abs(beta / se), df)]
md[, tc := qt(0.975, df)]
md[, lo_t := beta - tc * se]
md[, hi_t := beta + tc * se]

tref <- md[, {
  bt_cell <- if (c.bm[1] == 0) 0 else beta_true
  rej_t <- mean(p_t < 0.05)
  rej_z <- mean(p_z < 0.05)
  list(n_reps = .N, df = df[1],
       rejection_rate_z = rej_z,
       mcse_rate_z = sqrt(rej_z * (1 - rej_z) / .N),
       rejection_rate_t = rej_t,
       mcse_rate_t = sqrt(rej_t * (1 - rej_t) / .N),
       coverage_z = mean(lo_z <= bt_cell & bt_cell <= hi_z),
       coverage_t = mean(lo_t <= bt_cell & bt_cell <= hi_t),
       binom_p_z = binom.test(sum(p_z < 0.05), .N, 0.05)$p.value,
       binom_p_t = binom.test(sum(p_t < 0.05), .N, 0.05)$p.value)
}, by = .(c.bm, N)]
setorder(tref, c.bm, N)
write.table(tref, paste0(prefix, '-mancl-tref.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('Mancl-DeRouen arm, z versus t(K - p) reference:\n')
print(tref, digits = 4)

nulls <- r[c.bm == 0 & converged & is.finite(p_interaction), {
  k <- sum(p_interaction < 0.05)
  list(n_reps = .N, rate = k / .N,
       z_vs_nominal = (k / .N - 0.05) / sqrt(0.05 * 0.95 / .N),
       binom_p = binom.test(k, .N, 0.05)$p.value)
}, by = .(procedure, N)]
setorder(nulls, procedure, N)
write.table(nulls, paste0(prefix, '-null-binomial.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('\nExact binomial tests of the null cells (stored reference):\n')
print(nulls, digits = 4)

kappa <- r[procedure != 'rmanova' & converged & is.finite(beta),
           .(kappa = mean(se) / sd(beta)), by = .(procedure, c.bm, N)]
kappa_wide <- dcast(kappa, procedure + c.bm ~ N, value.var = 'kappa')
write.table(kappa_wide, paste0(prefix, '-kappa.txt'),
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('\nCalibration ratio kappa by arm and cell:\n')
print(kappa_wide, digits = 3)
invisible(NULL)
