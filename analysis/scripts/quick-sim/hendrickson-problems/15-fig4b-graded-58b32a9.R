## Figure 4B of Hendrickson et al. (2020) rerun on their own code (commit
## 58b32a9), as published and with graded coupling (docs/37, Appendix A;
## paper 01, Appendix B.1), on the Figure 4B grid: the four published
## trial designs x N {35, 70} x t1/2 {0, 0.1, 0.2} x the five censoring
## sets, published analysis options. c.bm = 0.25, below the
## compound-symmetry ceiling of 0.256, so neither version repairs a
## matrix and the coupling rule is the only difference. The published
## Figure 4B used c.bm = 0.6, above every covariance ceiling, which
## would confound the comparison with the repair. Both versions use the
## same seed for each (design, N, t1/2) cell.
##
## Requires 00-fetch-58b32a9.sh (censordata.R may come from H58_EXTRA).
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/15-fig4b-graded-58b32a9.R \
##     [--reps R] [--cores K] [--cbm C]
suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(corpcor)
  library(MASS); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.numeric(args[which(args == k) + 1]) else d
reps <- as.integer(arg('--reps', 200))
cores <- as.integer(arg('--cores', 8))
cbm <- arg('--cbm', 0.25)
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
extra <- Sys.getenv('H58_EXTRA', file.path(src, 'R'))
find_src <- function(f) {
  p <- c(file.path(src, 'R', f), file.path(extra, f))
  p <- p[file.exists(p)]
  if (!length(p)) stop('missing 58b32a9 source ', f)
  p[1]
}
out_dir <- 'analysis/data/quick-sim/hendrickson-problems/fig4b-graded'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## The graded-coupling replacement, as in 14-graded-patch-58b32a9.R.
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
graded_file <- file.path(out_dir, 'generateData-graded.R')
writeLines(c(gd[seq_len(i - 1)], patch, gd[(i + 4):length(gd)]), graded_file)

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
gen_files <- c(published = file.path(src, 'R', 'generateData.R'), graded = graded_file)

pub <- new.env(); load(file.path(src, 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
ap <- expand.grid(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE)
design_names <- c('OL', 'OL+BDC', 'CO', 'N-of-1')
## The four Figure 4 designs, rebuilt with the 58b32a9 buildtrialdesign()
## from the definitions in its publication vignette (the saved objects in
## results_core do not match the structure this code version expects).
btd <- make_env(gen_files[['published']])$buildtrialdesign
tds <- list(
  btd(name_longform = 'open label', name_shortform = 'OL',
      timepoints = cumsum(rep(2.5, 8)), timeptnames = paste0('OL', 1:8),
      expectancies = rep(1, 8), ondrug = list(pathA = rep(1, 8))),
  btd(name_longform = 'open label+blinded discontinuation', name_shortform = 'OL+BDC',
      timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
      timeptnames = c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2', 'BD3', 'BD4'),
      expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
      ondrug = list(pathA = c(1, 1, 1, 1, 1, 1, 0, 0), pathB = c(1, 1, 1, 1, 1, 0, 0, 0))),
  btd(name_longform = 'traditional crossover', name_shortform = 'CO',
      timepoints = cumsum(rep(2.5, 8)), timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
      expectancies = rep(.5, 8),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0), pathB = c(0, 0, 0, 0, 1, 1, 1, 1))),
  btd(name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
      timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
      timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
      expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0), pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                    pathC = c(1, 1, 1, 0, 0, 0, 1, 0), pathD = c(1, 1, 1, 0, 0, 0, 0, 1))))
cens_names <- c('No censoring', as.character(ps$censorparams$names))

cells <- CJ(version = names(gen_files), design = 1:4, N = c(35, 70), t_half = c(0, 0.1, 0.2))
cells[, seed := 1500L + .GRP, by = .(design, N, t_half)]
run <- function(k) {
  cl <- cells[k]
  env <- make_env(gen_files[[cl$version]])
  mp <- ps$modelparams[ps$modelparams$N == cl$N, , drop = FALSE][1, , drop = FALSE]
  mp$c.bm <- cbm
  mp$carryover_t1half <- cl$t_half
  set.seed(cl$seed)
  res <- env$generateSimulatedResults(
    trialdesigns = tds[cl$design], respparamsets = ps$respparamsets[1],
    blparamsets = ps$blparamsets[1], censorparams = ps$censorparams,
    modelparams = mp, simparam = list(Nreps = reps, progressiveSave = FALSE,
                                      saveunit2start = 1),
    analysisparams = ap, rawdataout = FALSE)
  r <- as.data.table(res$results)
  r[, .(power = mean(p < 0.05, na.rm = TRUE), mean_beta = mean(beta, na.rm = TRUE),
        n = .N), by = censorparamset][
    , `:=`(version = cl$version, design = design_names[cl$design], N = cl$N,
           t_half = cl$t_half, repairs = env$.repairs)]
}
out <- mclapply(seq_len(nrow(cells)), function(k) tryCatch(run(k), error = function(e)
  conditionMessage(e)), mc.cores = cores, mc.preschedule = FALSE)
bad <- which(!vapply(out, is.data.frame, logical(1)))
if (length(bad)) stop('cells failed: ', paste(unique(unlist(out[bad])), collapse = ' | '))
out <- rbindlist(out)
out[, censoring := cens_names[censorparamset + 1]]
fwrite(out, file.path(out_dir, sprintf('fig4b-graded-cbm%.2f-reps%d.csv', cbm, reps)))
print(dcast(out[censorparamset == 0], design + N ~ version + t_half, value.var = 'power'),
      digits = 2)
cat('repairs by version:\n'); print(out[, .(repairs = sum(unique(repairs))), by = version])
