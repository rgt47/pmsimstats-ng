## Paper 01, Stage 2 (referee comment M6, 2026-11-03). Reruns the
## original code (commit 58b32a9, as published) at t1/2 = 0.1 week, Hybrid,
## N = 70, c.bm = 0.25, no censoring, with more replicates and fresh seeds.
## Mirrors 14-graded-patch-58b32a9.R except: published version only,
## censorparams = NA (only the uncensored analysis is run), and seeds
## seed0 + chunk. '--seed0 1401 --chunks 1 --reps 500' reproduces the
## stored 0.126 of graded-patch/graded-patch-reps500.csv.
##
## Requires 00-fetch-58b32a9.sh. Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/17-published-th01-more-reps.R \
##     --reps 250 --chunks 16 --cores 8 \
##     --out analysis/data/quick-sim/hendrickson-problems/graded-patch/published-th0.1-reps4000.csv
suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(corpcor)
  library(MASS); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.numeric(args[which(args == k) + 1]) else d
reps <- arg('--reps', 250)
chunks <- arg('--chunks', 16)
cores <- arg('--cores', 8)
th <- arg('--thalf', 0.1)
seed0 <- as.integer(arg('--seed0', 7000))
out_file <- if ('--out' %in% args) args[which(args == '--out') + 1] else NA
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
env <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R', 'censordata.R', 'lme_analysis.R',
            'generateSimulatedResults.R', 'generateData.R'))
  sys.source(file.path(src, 'R', f), envir = env)
env$tic <- function(...) invisible(NULL)
env$toc <- function(...) invisible(NULL)
env$.repairs <- 0L
env$make.positive.definite <- function(m, tol) {
  env$.repairs <- env$.repairs + 1L
  corpcor::make.positive.definite(m, tol = tol)
}
pub <- new.env(); load(file.path(src, 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
ap <- expand.grid(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE)
td <- env$buildtrialdesign(
  name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0), pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                pathC = c(1, 1, 1, 0, 0, 0, 1, 0), pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))
mp0 <- ps$modelparams[2, , drop = FALSE]
stopifnot(mp0$N == 70)
run <- function(k) {
  mp <- mp0; mp$c.bm <- 0.25; mp$carryover_t1half <- th
  set.seed(seed0 + k)
  res <- env$generateSimulatedResults(
    trialdesigns = list(td), respparamsets = ps$respparamsets[1],
    blparamsets = ps$blparamsets[1], censorparams = NA,
    modelparams = mp, simparam = list(Nreps = reps, progressiveSave = FALSE,
                                      saveunit2start = 1),
    analysisparams = ap, rawdataout = FALSE)
  r <- as.data.table(res$results)[censorparamset == 0]
  data.table(chunk = k, reps = nrow(r), rej = sum(r$p < 0.05), mean_beta = mean(r$beta),
             repairs = env$.repairs)
}
t0 <- Sys.time()
out <- rbindlist(mclapply(seq_len(chunks), run, mc.cores = cores))
print(out)
R <- sum(out$reps); p <- sum(out$rej) / R
cat(sprintf('t_half %.2f: R = %d, power = %.4f, MCSE = %.4f, mean beta = %.4f, repairs = %d, %.1f min\n',
            th, R, p, sqrt(p * (1 - p) / R), weighted.mean(out$mean_beta, out$reps),
            sum(out$repairs), as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
if (!is.na(out_file)) fwrite(cbind(t_half = th, out), out_file)
