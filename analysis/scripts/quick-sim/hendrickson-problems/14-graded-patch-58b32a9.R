## Check of the coupling change proposed in docs/37, Appendix A, run on
## the 58b32a9 ('orig') code itself. The published generateData() is
## loaded twice: as published, and with only the biomarker-coupling lines
## (generateData.R lines 134-137) replaced by graded coupling. Both feed
## the published generateSimulatedResults() and lme_analysis(), Hybrid
## design, N = 70, c.bm = 0.25 (below the compound-symmetry ceiling of
## 0.256, so no matrix is repaired), t1/2 in {0, 0.1, 0.2, 1}, no
## censoring, published analysis options.
##
## Requires 00-fetch-58b32a9.sh. Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/14-graded-patch-58b32a9.R \
##     [--reps R] [--cores K]
suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(corpcor)
  library(MASS); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
reps <- arg('--reps', 500L)
cores <- arg('--cores', 8L)
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
## As in 10-strawman-58b32a9.R: censordata.R may live in H58_EXTRA.
extra <- Sys.getenv('H58_EXTRA', file.path(src, 'R'))
find_src <- function(f) {
  p <- c(file.path(src, 'R', f), file.path(extra, f))
  p <- p[file.exists(p)]
  if (!length(p)) stop('missing 58b32a9 source ', f)
  p[1]
}
out_dir <- 'analysis/data/quick-sim/hendrickson-problems/graded-patch'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## The replacement shown in docs/37, Appendix A (graded coupling: full on
## drug, decaying with the drug effect off drug, zero before exposure).
patch <- c(
  '      if(d[p]$tod>0){',
  '        g<-1',
  '      }else if(d[p]$tsd>0 && modelparam$carryover_t1half>0){',
  '        g<-(1/2)^(d[p]$tsd/modelparam$carryover_t1half)',
  '      }else{',
  '        g<-0',
  '      }',
  '      correlations[n1,\'bm\']<-modelparam$c.bm*g',
  '      correlations[\'bm\',n1]<-modelparam$c.bm*g')
gd <- readLines(file.path(src, 'R', 'generateData.R'))
i <- grep('if(means[which(n1==labels)]!=0){', gd, fixed = TRUE)
stopifnot(length(i) == 1, grepl('^\\s*}\\s*$', gd[i + 3]))
gd_graded <- c(gd[seq_len(i - 1)], patch, gd[(i + 4):length(gd)])
writeLines(gd_graded, file.path(out_dir, 'generateData-graded.R'))

make_env <- function(gen_file) {
  env <- new.env()
  for (f in c('utilities.R', 'buildtrialdesign.R', 'censordata.R', 'lme_analysis.R',
              'generateSimulatedResults.R'))
    sys.source(find_src(f), envir = env)
  sys.source(gen_file, envir = env)
  env$tic <- function(...) invisible(NULL)
  env$toc <- function(...) invisible(NULL)
  env$.repairs <- 0L
  env$make.positive.definite <- function(m, tol) {
    env$.repairs <- env$.repairs + 1L
    corpcor::make.positive.definite(m, tol = tol)
  }
  env
}
envs <- list(published = make_env(file.path(src, 'R', 'generateData.R')),
             graded = make_env(file.path(out_dir, 'generateData-graded.R')))

pub <- new.env(); load(file.path(src, 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
ap <- expand.grid(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE)
td <- envs$published$buildtrialdesign(
  name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0), pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                pathC = c(1, 1, 1, 0, 0, 0, 1, 0), pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))
mp0 <- ps$modelparams[2, , drop = FALSE]   # N = 70, published correlations
stopifnot(mp0$N == 70)

cells <- CJ(version = names(envs), t_half = c(0, 0.1, 0.2, 1))
cells[, seed := 1400L + .GRP, by = t_half]   # same seed for both versions at a half-life
run <- function(k) {
  cl <- cells[k]
  env <- envs[[cl$version]]
  env$.repairs <- 0L
  mp <- mp0; mp$c.bm <- 0.25; mp$carryover_t1half <- cl$t_half
  set.seed(cl$seed)
  res <- env$generateSimulatedResults(
    trialdesigns = list(td), respparamsets = ps$respparamsets[1],
    ## The published censoring sets are passed (the loop expects them);
    ## only the uncensored results (censorparamset 0) are used below.
    blparamsets = ps$blparamsets[1], censorparams = ps$censorparams,
    modelparams = mp, simparam = list(Nreps = reps, progressiveSave = FALSE,
                                      saveunit2start = 1),
    analysisparams = ap, rawdataout = FALSE)
  r <- as.data.table(res$results)[censorparamset == 0]
  data.table(version = cl$version, t_half = cl$t_half, reps = nrow(r),
             mean_beta = mean(r$beta), power = mean(r$p < 0.05),
             mcse = sqrt(mean(r$p < 0.05) * (1 - mean(r$p < 0.05)) / nrow(r)),
             repairs = env$.repairs)
}
out <- mclapply(seq_len(nrow(cells)), function(k) tryCatch(run(k), error = function(e)
  conditionMessage(e)), mc.cores = cores)
bad <- which(!vapply(out, is.data.frame, logical(1)))
if (length(bad)) stop('cells failed: ', paste(unique(unlist(out[bad])), collapse = ' | '))
out <- rbindlist(out)
fwrite(out, file.path(out_dir, sprintf('graded-patch-reps%d.csv', reps)))
print(dcast(melt(out, id.vars = c('version', 't_half'),
                 measure.vars = c('mean_beta', 'power', 'repairs')),
            variable + version ~ t_half), digits = 3)
