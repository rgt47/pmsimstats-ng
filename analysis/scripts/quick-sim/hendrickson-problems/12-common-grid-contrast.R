## Common-grid contrast: the strawman against the candidate fixes, all run
## through Hendrickson's own 58b32a9 pipeline, so that they differ only in
## the data-generating construct.
##
## Configurations:
##   S  strawman: their generateData() untouched, repair switched on
##   B  AR(1) within factor, decaying cross-factor, step coupling
##   F  AR(1), decaying (covar) cross-factor, graded coupling
##   C  AR(1), separable cross-factor, graded coupling
##   D  CS as published, interaction in the mean (mean moderation)
##   E  AR(1), separable, mean moderation
## B-F use the generator of 04-power-simulation.R (definitions read from
## that file unchanged). Their generateData() still assembles the means,
## names and outcome; only its call to mvrnorm() is intercepted: the
## matrix is rebuilt for the same path with the 04 generator, the mean
## vector is checked against theirs, and for D and E the on-drug BR is
## shifted by c.bm * sd_BR * z. Their repair is off for B-F: a set whose
## matrix fails corpcor::is.positive.definite() is not run and is
## recorded as infeasible (that count is the Problem 1 measure; see also
## 11-pd-repair-count.R).
##
## Everything else is theirs: trial designs (rebuilt as in their vignette),
## allocation of subjects to paths, the four censoring patterns
## (censordata), the analysis (lme_analysis) and the replicate layout.
## --analysis car1 keeps their lme_analysis() but sends its model fit to
## nlme::lme with corCAR1(form = ~ t | ptID), same fixed effects and
## random intercept.
##
## Grid: OL, OL+BDC, CO, N-of-1; N = 70; c.bm {0, .25, .3};
## t1/2 {0, .1, .2, .5, 1}; no censoring + 4 censoring patterns.
## Seeds depend on the grid cell only, so every configuration and both
## analyses use the same seed for a cell (the data still differ between
## configurations, since the matrices differ).
##
## Each set is saved as it finishes; a rerun resumes. Run from the
## repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/12-common-grid-contrast.R \
##     [--reps R] [--analysis published|car1] [--configs SBFCDE] [--cores K]
suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(corpcor)
  library(MASS); library(nlme); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
reps <- as.integer(arg('--reps', '100'))
analysis <- arg('--analysis', 'published')
stopifnot(analysis %in% c('published', 'car1'))
configs <- strsplit(arg('--configs', 'SBFCDE'), '')[[1]]
stopifnot(all(configs %in% c('S', 'B', 'F', 'C', 'D', 'E')))
cores <- as.integer(arg('--cores', '8'))

here <- 'analysis/scripts/quick-sim/hendrickson-problems'
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
extra <- Sys.getenv('H58_EXTRA', file.path(src, 'R'))
find_src <- function(f) {
  p <- c(file.path(src, 'R', f), file.path(extra, f))
  p <- p[file.exists(p)]
  if (!length(p)) stop('missing 58b32a9 source ', f)
  p[1]
}
out_root <- 'analysis/data/quick-sim/hendrickson-problems/contrast'
set_dir <- file.path(out_root, sprintf('sets-reps%d-%s', reps, analysis))
dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------- their pipeline ----------------
env <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R', 'generateData.R',
            'censordata.R', 'lme_analysis.R', 'generateSimulatedResults.R'))
  sys.source(find_src(f), envir = env)
env$tic <- function(...) invisible(NULL)
env$toc <- function(...) invisible(NULL)

## ---------------- the 04 generator, read verbatim ----------------
keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
dgp <- new.env()
for (e in parse(file.path(here, '04-power-simulation.R'))) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))

## ---------------- configuration hooks ----------------
## The current configuration and path, set by the generateData() wrapper
## and read by the mvrnorm() intercept.
cur <- new.env()
cur$config <- 'S'
cur$repairs <- 0L
orig_generateData <- env$generateData
env$generateData <- function(modelparam, respparam, blparam, trialdesign,
                             empirical = FALSE, makePositiveDefinite = TRUE) {
  cur$td <- trialdesign
  cur$mp <- modelparam
  orig_generateData(modelparam, respparam, blparam, trialdesign, empirical,
                    makePositiveDefinite = (cur$config == 'S'))
}
env$make.positive.definite <- function(m, tol) {
  cur$repairs <- cur$repairs + 1L
  corpcor::make.positive.definite(m, tol = tol)
}
path_of <- function(td) {
  list(w = cumsum(td$t_wk), e = td$e, on = td$tod > 0, tod = td$tod,
       tsd = td$tsd, tpb = td$tpb)
}
env$mvrnorm <- function(n, mu, Sigma, empirical = FALSE) {
  if (cur$config == 'S') return(MASS::mvrnorm(n, mu, Sigma, empirical = empirical))
  pi <- path_of(cur$td)
  b <- dgp$build(pi, cur$config, cur$mp$c.bm, cur$mp$carryover_t1half)
  ## 04 hardcodes the published parameters to 5 decimals (means agree to
  ## ~5e-6), so check the mapping at 1e-5 and draw with their exact means
  ## and SDs; only the correlation matrix comes from 04.
  sds <- sqrt(diag(Sigma))
  if (!isTRUE(all.equal(unname(b$mu), unname(mu), tolerance = 1e-5)) ||
      !isTRUE(all.equal(unname(sqrt(diag(b$Sigma))), unname(sds), tolerance = 1e-5)))
    stop('04 generator means or SDs differ from 58b32a9 generateData()')
  S <- outer(sds, sds) * b$R
  if (!is.positive.definite(S)) stop('infeasible matrix reached the draw')
  X <- MASS::mvrnorm(n, mu, S)
  if (b$mean_mod) {
    nv <- b$n
    br <- (3 + 2 * nv):(2 + 3 * nv)
    z <- (X[, 1] - mu[1]) / sds[1]
    sd_br <- sds[br]
    X[, br[pi$on]] <- X[, br[pi$on]] + outer(cur$mp$c.bm * z, sd_br[pi$on])
  }
  colnames(X) <- colnames(Sigma)
  X
}

## ---------------- corCAR1 fit inside their lme_analysis() ----------------
if (analysis == 'car1') {
  env$lmer <- function(formula, data, ...) {
    fixed <- lme4::nobars(formula)
    d <- as.data.frame(data)
    fit_with <- function(ctrl) nlme::lme(
      fixed, random = ~ 1 | ptID, data = d, na.action = na.omit,
      correlation = nlme::corCAR1(form = ~ t | ptID), control = ctrl)
    fit <- tryCatch(fit_with(nlme::lmeControl()),
                    error = function(e) fit_with(nlme::lmeControl(opt = 'optim')))
    structure(list(fit = fit), class = 'car1fit')
  }
  ## Found by UseMethod() from inside lme_analysis(), whose enclosure is env.
  env$summary.car1fit <- function(object, ...) {
    tt <- summary(object$fit)$tTable[, c('Value', 'Std.Error', 't-value', 'p-value'),
                                     drop = FALSE]
    colnames(tt) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
    list(coefficients = tt, fitMsgs = character(0),
         optinfo = list(conv = list(lme4 = list(messages = list()))))
  }
}

## ---------------- designs and parameters (as in 10) ----------------
pub <- new.env(); load(file.path(src, 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
btd <- env$buildtrialdesign
tds <- list(
  btd(name_longform = 'open label', name_shortform = 'OL',
      timepoints = env$cumulative(rep(2.5, 8)), timeptnames = paste0('OL', 1:8),
      expectancies = rep(1, 8), ondrug = list(pathA = rep(1, 8))),
  btd(name_longform = 'open label+blinded discontinuation', name_shortform = 'OL+BDC',
      timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
      timeptnames = c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2', 'BD3', 'BD4'),
      expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
      ondrug = list(pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
                    pathB = c(1, 1, 1, 1, 1, 0, 0, 0))),
  btd(name_longform = 'traditional crossover', name_shortform = 'CO',
      timepoints = env$cumulative(rep(2.5, 8)),
      timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
      expectancies = rep(.5, 8),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
                    pathB = c(0, 0, 0, 0, 1, 1, 1, 1))),
  btd(name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
      timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
      timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
      expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0), pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                    pathC = c(1, 1, 1, 0, 0, 0, 1, 0), pathD = c(1, 1, 1, 0, 0, 0, 0, 1))))
ap <- expand.grid(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE)
mp_grid <- expand.grid(N = 70, c.bm = c(0, .25, .3),
                       carryover_t1half = c(0, .1, .2, .5, 1),
                       c.tv = .8, c.pb = .8, c.br = .8, c.cf1t = .2, c.cfct = .1)
stopifnot(dgp$RHO == .8, dgp$C1 == .2, dgp$CX == .1)
cells <- CJ(td = seq_along(tds), mp = seq_len(nrow(mp_grid)))
cells[, cell := .I]

feasible <- function(config, td, mp) {
  if (config == 'S') return(TRUE)
  all(vapply(tds[[td]]$trialpaths, function(p) is.positive.definite(
    dgp$build(path_of(p), config, mp$c.bm, mp$carryover_t1half)$Sigma), logical(1)))
}

run_set <- function(config, k) {
  g <- cells[k]
  f <- file.path(set_dir, sprintf('%s-%03d.rds', config, g$cell))
  if (file.exists(f)) return(readRDS(f)$status)
  mp <- mp_grid[g$mp, , drop = FALSE]
  info <- data.table(config = config, cell = g$cell, trialdesign = g$td,
                     design = tds[[g$td]]$metadata$name_shortform, mp[, 1:3])
  if (!feasible(config, g$td, mp)) {
    saveRDS(list(status = cbind(info, feasible = FALSE, repairs = NA_integer_,
                                seconds = 0), results = NULL), f)
    return(readRDS(f)$status)
  }
  cur$config <- config
  cur$repairs <- 0L
  set.seed(31000 + g$cell)
  t0 <- proc.time()[['elapsed']]
  res <- env$generateSimulatedResults(
    trialdesigns = tds[g$td], respparamsets = ps$respparamsets[1],
    blparamsets = ps$blparamsets[1], censorparams = ps$censorparams,
    modelparams = mp,
    simparam = list(Nreps = reps, progressiveSave = FALSE, saveunit2start = 1),
    analysisparams = ap, rawdataout = FALSE)
  r <- res$results
  r[, `:=`(config = config, cell = g$cell, trialdesign = g$td, modelparamset = g$mp,
           irep = (irep - 1L) %% reps + 1L)]
  r[, warning := NULL]
  status <- cbind(info, feasible = TRUE, repairs = cur$repairs,
                  seconds = proc.time()[['elapsed']] - t0)
  saveRDS(list(status = status, results = r), f)
  status
}

log_f <- file.path(set_dir, 'progress.log')
cat(sprintf('%s start: analysis %s, configs %s, %d cells, %d reps, %d cores\n',
            format(Sys.time()), analysis, paste(configs, collapse = ''),
            nrow(cells), reps, cores), file = log_f, append = TRUE)
for (cf in configs) {
  t0 <- proc.time()[['elapsed']]
  st <- mclapply(cells$cell, function(k) tryCatch(run_set(cf, k), error = function(e)
    data.table(config = cf, cell = k, error = conditionMessage(e))),
    mc.cores = cores, mc.preschedule = FALSE)
  st <- rbindlist(st, fill = TRUE)
  n_err <- if ('error' %in% names(st)) sum(!is.na(st$error)) else 0L
  msg <- sprintf('%s config %s: %d sets, %d run, %d infeasible, %d errors, %.1f min\n',
                 format(Sys.time()), cf, nrow(st), sum(st$feasible %in% TRUE),
                 sum(st$feasible %in% FALSE), n_err,
                 (proc.time()[['elapsed']] - t0) / 60)
  cat(msg); cat(msg, file = log_f, append = TRUE)
  if (n_err) print(unique(st[!is.na(error), .(cell, error)])[1:min(5, n_err)])
}
