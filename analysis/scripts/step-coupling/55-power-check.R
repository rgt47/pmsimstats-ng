## docs/55, Section 7.3 and 8.3: power of the published analysis on data
## from three versions of generateData(): the 58b32a9 code as published
## (step coupling); the same file with only lines 134-137 replaced by
## graded coupling (report 01, Appendix B.1); and the 2024 upstream file
## (ratio coupling, scalefactor = 2, its default). Each feeds the 58b32a9
## lme_analysis() (lmer, useDE = FALSE). Hybrid design, N = 70
## (18/18/17/17), c.bm = 0.25, no censoring. The data are generated as
## the published generateSimulatedResults() generates them (all
## replicates of a path in one generateData() call, makePositiveDefinite
## = TRUE), and each replicate is analyzed with full_model_out = TRUE so
## that the drug coefficient can be recorded as well as the interaction.
## The same seed is used for every version at a half-life.
##
## Run from the repository root:
##   Rscript analysis/scripts/step-coupling/55-power-check.R \
##     [--reps R] [--cores K]
suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(corpcor)
  library(MASS); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) {
  if (k %in% args) as.integer(args[which(args == k) + 1]) else d
}
reps <- arg('--reps', 500L)
cores <- arg('--cores', 8L)
t0 <- Sys.time()
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
head_file <- paste0('analysis/scripts/quick-sim/',
  'hendrickson-original-comparison/vendored-hendrickson-generateData.R')
out_dir <- 'analysis/data/quick-sim/step-coupling'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

graded_patch <- c(
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
stopifnot(length(i) == 1, i == 134, grepl('^\\s*}\\s*$', gd[i + 3]))
graded_file <- file.path(out_dir, 'generateData-graded.R')
writeLines(c(gd[seq_len(i - 1)], graded_patch, gd[(i + 4):length(gd)]),
  graded_file)

make_env <- function(gen_file) {
  env <- new.env()
  for (f in c('utilities.R', 'buildtrialdesign.R', 'lme_analysis.R')) {
    sys.source(file.path(src, 'R', f), envir = env)
  }
  sys.source(gen_file, envir = env)
  env$.repairs <- 0L
  env$make.positive.definite <- function(m, tol) {
    env$.repairs <- env$.repairs + 1L
    corpcor::make.positive.definite(m, tol = tol)
  }
  env
}
files <- c(published = file.path(src, 'R', 'generateData.R'),
  graded = graded_file, upstream = head_file)
envs <- lapply(files, make_env)

pub <- new.env()
load(file.path(src, 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
op <- data.frame(useDE = FALSE, t_random_slope = FALSE, full_model_out = TRUE)
td <- envs$published$buildtrialdesign(
  name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1), pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))$trialpaths
rp <- ps$respparamsets[[1]]$param
bp <- ps$blparamsets[[1]]$param
mp0 <- ps$modelparams[2, , drop = FALSE]
stopifnot(mp0$N == 70)

cells <- CJ(version = names(envs), t_half = c(0, 0.001, 0.01, 0.1, 1))
cells[, seed := 5500L + .GRP, by = t_half]
run <- function(k) {
  cl <- cells[k]
  env <- envs[[cl$version]]
  env$.repairs <- 0L
  mp <- mp0
  mp$c.bm <- 0.25
  mp$carryover_t1half <- cl$t_half
  set.seed(cl$seed)
  nP <- length(td)
  Ns <- mp$N %/% nP + c(rep(1, mp$N %% nP), rep(0, nP - mp$N %% nP))
  dat <- rbindlist(lapply(seq_len(nP), function(iP) {
    m <- mp
    m$N <- Ns[iP] * reps
    d <- env$generateData(m, rp, bp, td[[iP]], FALSE, TRUE)
    d[, path := iP]
    d[, replicate := rep(1:reps, Ns[iP])]
    d
  }))
  out <- rbindlist(lapply(seq_len(reps), function(iS) {
    dd <- dat[replicate == iS]
    a <- env$lme_analysis(td, dd, op)
    cf <- summary(a$fit)$coefficients
    data.table(version = cl$version, t_half = cl$t_half, rep = iS,
      beta = cf['bm:DbTRUE', 'Estimate'], p = cf['bm:DbTRUE', 'Pr(>|t|)'],
      b_D = cf['DbTRUE', 'Estimate'],
      b_D_centered = cf['DbTRUE', 'Estimate'] +
        cf['bm:DbTRUE', 'Estimate'] * mean(dd$bm),
      p_D = cf['DbTRUE', 'Pr(>|t|)'])
  }))
  out[, repairs := env$.repairs]
  out
}
res <- mclapply(seq_len(nrow(cells)), function(k)
  tryCatch(run(k), error = function(e) conditionMessage(e)), mc.cores = cores)
bad <- which(!vapply(res, is.data.frame, logical(1)))
if (length(bad)) {
  stop('cells failed: ', paste(unique(unlist(res[bad])), collapse = ' | '))
}
res <- rbindlist(res)
fwrite(res, file.path(out_dir, sprintf('power-check-reps%d-raw.csv', reps)))
summ <- res[, .(reps = .N, mean_beta = mean(beta), power = mean(p < 0.05),
  mcse = sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / .N),
  mean_b_D = mean(b_D), mean_b_D_centered = mean(b_D_centered),
  sd_b_D_centered = sd(b_D_centered), power_D = mean(p_D < 0.05),
  repairs = repairs[1]), by = .(version, t_half)]
fwrite(summ, file.path(out_dir, sprintf('power-check-reps%d.csv', reps)))
print(summ, digits = 4)
el <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))
cat(sprintf('elapsed %.1f minutes on %d cores, %d replicates per cell\n',
  el, cores, reps))
writeLines(sprintf('elapsed_minutes=%.1f cores=%d reps=%d', el, cores, reps),
  file.path(out_dir, sprintf('power-check-reps%d-runtime.txt', reps)))
