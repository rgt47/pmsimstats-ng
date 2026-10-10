## docs/53, Section 8: does the anchored fix alone change the published
## analysis? Adapted from analysis/scripts/quick-sim/hendrickson-problems/
## 14-graded-patch-58b32a9.R. Four versions of the 58b32a9 generateData():
## published; anchored (lines 83-92 replaced); graded coupling (the
## docs/37 Appendix A patch); graded plus anchored. Each feeds the
## published generateSimulatedResults() and lme_analysis(), Hybrid
## design, N = 70, c.bm = 0.25, no censoring, published analysis options.
## The same seed is used for all versions at a half-life, and the
## anchored change alters only the mean vector, so paired versions share
## identical noise.
##
## Run from the repository root:
##   Rscript <path>/53-power-check.R [--reps R] [--cores K] [--out DIR]
suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(corpcor)
  library(MASS); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) {
  if (k %in% args) args[which(args == k) + 1] else d
}
reps <- as.integer(arg('--reps', 500L))
cores <- as.integer(arg('--cores', 8L))
out_dir <- arg('--out', file.path(tempdir(), 'docs53'))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
graded_file <- paste0('analysis/data/quick-sim/hendrickson-problems/',
  'graded-patch/generateData-graded.R')

anchor <- function(lines) {
  i <- grep('brmeans[p]<-brmeans[p]+brmeans[p-1]*', lines, fixed = TRUE)
  j <- grep('brmeans<-modgompertz(d$tod', lines, fixed = TRUE)
  stopifnot(length(i) == 1, length(j) == 1, i - j == 5)
  new <- c(
    '      brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)',
    '      if(nP>1){',
    '        last_on<-0',
    '        for(p in 1:nP){',
    '          if(d[p]$onDrug){',
    '            last_on<-brmeans[p]',
    '          }else if(d[p]$tsd>0){',
    paste0('            brmeans[p]<-brmeans[p]+last_on*',
      '(1/2)^(d$tsd[p]/modelparam$carryover_t1half)'),
    '          }',
    '        }',
    '      }')
  c(lines[seq_len(j - 1)], new, lines[(j + 10):length(lines)])
}
pub_lines <- readLines(file.path(src, 'R', 'generateData.R'))
gr_lines <- readLines(graded_file)
files <- c(published = file.path(src, 'R', 'generateData.R'),
  anchored = file.path(out_dir, 'gd-anchored.R'), graded = graded_file,
  graded_anchored = file.path(out_dir, 'gd-graded-anchored.R'))
writeLines(anchor(pub_lines), files[['anchored']])
writeLines(anchor(gr_lines), files[['graded_anchored']])

make_env <- function(gen_file) {
  env <- new.env()
  for (f in c('utilities.R', 'buildtrialdesign.R', 'censordata.R',
    'lme_analysis.R', 'generateSimulatedResults.R')) {
    sys.source(file.path(src, 'R', f), envir = env)
  }
  sys.source(gen_file, envir = env)
  env$tic <- function(...) invisible(NULL)
  env$toc <- function(...) invisible(NULL)
  env
}
envs <- lapply(files, make_env)

pub <- new.env()
load(file.path(src, 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
ap <- expand.grid(useDE = FALSE, t_random_slope = FALSE,
  full_model_out = FALSE)
td <- envs$published$buildtrialdesign(
  name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1), pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))
mp0 <- ps$modelparams[2, , drop = FALSE]
stopifnot(mp0$N == 70)

cells <- CJ(version = names(envs), t_half = c(0, 0.1, 0.2, 1, 2))
cells[, seed := 5300L + .GRP, by = t_half]
run <- function(k) {
  cl <- cells[k]
  env <- envs[[cl$version]]
  mp <- mp0
  mp$c.bm <- 0.25
  mp$carryover_t1half <- cl$t_half
  set.seed(cl$seed)
  res <- env$generateSimulatedResults(
    trialdesigns = list(td), respparamsets = ps$respparamsets[1],
    blparamsets = ps$blparamsets[1], censorparams = ps$censorparams,
    modelparams = mp, simparam = list(Nreps = reps, progressiveSave = FALSE,
      saveunit2start = 1),
    analysisparams = ap, rawdataout = FALSE)
  r <- as.data.table(res$results)[censorparamset == 0]
  data.table(version = cl$version, t_half = cl$t_half, rep = seq_len(nrow(r)),
    beta = r$beta, p = r$p)
}
out <- mclapply(seq_len(nrow(cells)), function(k) {
  tryCatch(run(k), error = function(e) conditionMessage(e))
}, mc.cores = cores)
bad <- which(!vapply(out, is.data.frame, logical(1)))
if (length(bad)) {
  stop('cells failed: ', paste(unique(unlist(out[bad])), collapse = ' | '))
}
out <- rbindlist(out)
fwrite(out, file.path(out_dir, sprintf('power-check-reps%d-raw.csv', reps)))

summ <- out[, .(reps = .N, mean_beta = mean(beta), power = mean(p < 0.05),
  mcse = sqrt(mean(p < 0.05) * (1 - mean(p < 0.05)) / .N)),
  by = .(version, t_half)]
fwrite(summ, file.path(out_dir, sprintf('power-check-reps%d.csv', reps)))
print(dcast(summ, t_half ~ version, value.var = 'power'), digits = 3)
print(dcast(summ, t_half ~ version, value.var = 'mean_beta'), digits = 4)

## Paired comparison: same seed, so rejection decisions can be compared
## replicate by replicate.
pair <- function(a, b) {
  w <- dcast(out[version %in% c(a, b)], t_half + rep ~ version,
    value.var = c('beta', 'p'))
  w[, .(pair = paste(b, 'vs', a), mean_dbeta = mean(get(paste0('beta_', b)) -
      get(paste0('beta_', a))),
    max_abs_dbeta = max(abs(get(paste0('beta_', b)) -
      get(paste0('beta_', a)))),
    discordant = sum((get(paste0('p_', b)) < 0.05) !=
      (get(paste0('p_', a)) < 0.05)),
    dpower = mean(get(paste0('p_', b)) < 0.05) -
      mean(get(paste0('p_', a)) < 0.05)), by = t_half]
}
pairs <- rbind(pair('published', 'anchored'),
  pair('graded', 'graded_anchored'))
fwrite(pairs, file.path(out_dir, sprintf('power-check-reps%d-pairs.csv', reps)))
print(pairs, digits = 3)
