## Verification of the orig-to-refined patches against pmsimstats at commit
## 58b32a9 ('orig') and against the stress-test driver. Called by verify.sh,
## which fetches orig, builds the patched trees and pins the driver:
##   bash analysis/patches/orig-to-refined/verify.sh [workdir]
## or, with an existing work directory, from the repository root:
##   Rscript analysis/patches/orig-to-refined/verify.R <workdir> [--fast]
## <workdir> holds the trees orig, core (orig + core patch), options (core +
## options patch) and full (core + options + analysis patches), and ref/ with
## the driver 11-strawman-stress-test.R and its helper 04-power-simulation.R as
## committed in 3238916 (the fixed reference). --fast uses fewer replicates
## in the Monte Carlo and analysis checks. One PASS/FAIL line per check.
suppressMessages({
  library(data.table); library(MASS); library(corpcor); library(nlme)
})
args <- commandArgs(trailingOnly = TRUE)
work <- args[1]
fast <- '--fast' %in% args
stopifnot(!is.na(work), dir.exists(file.path(work, 'full')))
n_fail <- 0L
report <- function(ok, label, detail = '') {
  cat(sprintf('%s %s%s\n', if (isTRUE(ok)) 'PASS' else 'FAIL', label,
              if (nzchar(detail)) paste0(': ', detail) else ''))
  if (!isTRUE(ok)) n_fail <<- n_fail + 1L
}
fmt <- function(x) formatC(x, format = 'g', digits = 3)

## ---------------- 0. dependencies and sourcing ----------------
need <- c('data.table', 'MASS', 'corpcor', 'lme4', 'lmerTest', 'nlme',
          'clubSandwich', 'mmrm')
have <- vapply(need, requireNamespace, logical(1), quietly = TRUE)
report(all(have), 'dependencies installed',
       if (all(have)) paste(need, collapse = ', ')
       else paste('missing:', paste(need[!have], collapse = ', ')))
if (!all(have[c('data.table', 'MASS', 'corpcor')])) stop('cannot continue')
## orig's lme_analysis() calls lmer() unqualified and reads lmerTest p-values
suppressMessages(library(lmerTest))

make_env <- function(tree) {
  env <- new.env()
  for (f in c('utilities.R', 'buildtrialdesign.R', 'generateData.R',
              'censordata.R', 'lme_analysis.R', 'generateSimulatedResults.R'))
    sys.source(file.path(work, tree, 'R', f), envir = env)
  env$tic <- function(...) invisible(NULL)
  env$toc <- function(...) invisible(NULL)
  env$.repairs <- 0L
  env$make.positive.definite <- function(m, tol) {
    env$.repairs <- env$.repairs + 1L
    corpcor::make.positive.definite(m, tol = tol)
  }
  ## Records the mean and covariance handed to the sampler; the draw itself
  ## is MASS::mvrnorm unchanged, so the random stream is not altered.
  env$mvrnorm <- function(n, mu, Sigma, ...) {
    env$.last <- list(mu = mu, Sigma = Sigma)
    MASS::mvrnorm(n, mu, Sigma, ...)
  }
  env
}
envs <- tryCatch(lapply(c(orig = 'orig', core = 'core', options = 'options',
                          full = 'full'), make_env),
                 error = function(e) e)
report(!inherits(envs, 'error'), 'source orig and the patched trees',
       if (inherits(envs, 'error')) conditionMessage(envs) else 'orig, core, options, full')
if (inherits(envs, 'error')) quit(status = 1)
orig <- envs$orig; core <- envs$core; opts <- envs$options; full <- envs$full

## The driver, as committed: only its definitions are evaluated.
load_driver <- function(f11, f04) {
  env <- new.env()
  env$here04 <- f04
  wanted <- c('keep', 'dgp', 'D', 'PIS', 'TV_MAX', 'base', 'cell', 'decay',
              'build', 'draw', 'dropout_visit', 'TH_A', 'GRID', 'dbc', 'na3',
              'lctrl', 'fit_lme', 'cr2', 'mb', 'fit_mmrm', 'mm', 'F_S', 'F_F1',
              'F_F2', 'F_E', 'F_A', 'F_EA', 'F_R', 'with_x', 'fit_all',
              'fit_extra', 'fit_extra2', 'make_trial')
  for (e in parse(f11)) {
    if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% wanted) eval(e, env)
    ## the loop that loads the 04-power-simulation.R definitions into dgp
    if (is.call(e) && identical(e[[1]], as.name('for')) &&
        identical(e[[2]], as.name('e'))) eval(e, env)
  }
  stopifnot(all(wanted %in% ls(env)), all(env$keep %in% ls(env$dgp)))
  env
}
suppressMessages(library(mmrm))
drv <- load_driver(file.path(work, 'ref', '11-strawman-stress-test.R'),
                   file.path(work, 'ref', '04-power-simulation.R'))
drv_dropout_visit <- drv$dropout_visit
dgp <- drv$dgp

## ---------------- shared inputs ----------------
pub <- new.env()
load(file.path(work, 'orig', 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
## The driver hard-codes the published parameters rounded to the digits
## below; the comparisons use the same rounded values. The published
## reproduction checks (Section 4) use the stored values themselves.
rp <- copy(ps$respparamsets[[1]]$param)
bp <- copy(ps$blparamsets[[1]]$param)
rp[cat %in% c('tv', 'pb'), max := 6.50647]
rp[cat == 'br', max := 10.98604]
bp[cat == 'BL', `:=`(m = dgp$bl_m, sd = dgp$bl_sd)]
bp[cat == 'bm', `:=`(m = dgp$bm_m, sd = dgp$bm_sd)]
stopifnot(max(abs(rp$max - ps$respparamsets[[1]]$param$max)) < 1e-4,
          max(abs(bp$sd - ps$blparamsets[[1]]$param$sd)) < 1e-4)
td_hyb <- orig$buildtrialdesign(
  name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0), pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                pathC = c(1, 1, 1, 0, 0, 0, 1, 0), pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))
td_co <- orig$buildtrialdesign(
  name_longform = 'traditional crossover', name_shortform = 'CO',
  timepoints = orig$cumulative(rep(2.5, 8)),
  timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
  expectancies = rep(.5, 8),
  ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0), pathB = c(0, 0, 0, 0, 1, 1, 1, 1)))
designs <- list(Hybrid = list(td = td_hyb, drv = dgp$designs$Hybrid),
                CO = list(td = td_co, drv = dgp$designs$CO))
for (nm in names(designs)) {
  x <- designs[[nm]]
  stopifnot(length(x$td$trialpaths) == length(x$drv$paths))
  for (k in seq_along(x$drv$paths))
    stopifnot(all.equal(x$td$trialpaths[[k]]$tod > 0, x$drv$paths[[k]] == 1))
}

mp <- function(cbm, th, N = 2L, ...) {
  m <- ps$modelparams[2, , drop = FALSE]
  m$N <- N; m$c.bm <- cbm; m$carryover_t1half <- th
  x <- list(...)
  for (k in names(x)) m[[k]] <- x[[k]]
  m
}
mp_refined <- function(cbm, th, N = 2L, ...) mp(cbm, th, N, rho = 0.7, c1 = 0.2, ...)
capture <- function(env, td_path, m, ...) {
  env$generateData(m, rp, bp, td_path, FALSE, TRUE, ...)
  env$.last
}
idx <- function(n) list(tv = 3:(2 + n), pb = (3 + n):(2 + 2 * n), br = (3 + 2 * n):(2 + 3 * n))
## Closed forms, computed from the trial path alone
g_graded <- function(tp, th, k = 1) {
  on <- tp$tod > 0
  w <- if (th > 0) {
    if (k == 1) 0.5^(tp$tsd / th) else exp(-(log(2)^(1 / k) * tp$tsd / th)^k)
  } else rep(0, nrow(tp))
  ifelse(on, 1, ifelse(tp$tsd > 0, w, 0))
}
br_anchored <- function(tp, th) {
  br <- orig$modgompertz(tp$tod, rp[cat == 'br']$max, rp[cat == 'br']$disp, rp[cat == 'br']$rate)
  last <- 0
  for (p in seq_len(nrow(tp))) {
    if (tp$tod[p] > 0) last <- br[p]
    else if (tp$tsd[p] > 0 && th > 0) br[p] <- last * 0.5^(tp$tsd[p] / th)
  }
  br
}

## ---------------- 1. refined DGP against the driver (base cell) ----------------
cat('\n## Refined DGP, stress-test base cell (rho 0.7, c1 0.2, c_bm 0.3, t1/2 0.5 week)\n')
cl0 <- drv$cell('curve-C07', 0.30)
stopifnot(cl0$rho == 0.7, cl0$t_half == 0.5, cl0$N == 70)
res1 <- list()
for (nm in names(designs)) {
  x <- designs[[nm]]
  for (k in seq_along(x$drv$paths)) {
    tp <- x$td$trialpaths[[k]]
    pi <- dgp$path_info(x$drv$t_wk, x$drv$e, x$drv$paths[[k]])
    b <- drv$build(pi, cl0)
    a <- capture(core, tp, mp_refined(0.3, 0.5))
    r <- capture(core, tp, mp_refined(0.3, 0.5), carryover = 'recursive')
    n <- nrow(tp); ix <- idx(n)
    Ra <- cov2cor(a$Sigma)
    res1[[length(res1) + 1]] <- data.table(
      design = nm, path = k,
      dR = max(abs(Ra - b$R)), dSigma = max(abs(a$Sigma - b$Sigma)),
      dcoup_drv = max(abs(Ra[1, ix$br] - b$R[1, ix$br])),
      dcoup_closed = max(abs(Ra[1, ix$br] - 0.3 * g_graded(tp, 0.5))),
      pre_exposure_coup = if (any(cumsum(tp$tod > 0) == 0))
        max(abs(Ra[1, ix$br][cumsum(tp$tod > 0) == 0])) else NA_real_,
      dmu_tvpb = max(abs(a$mu[c(ix$tv, ix$pb)] - b$mu[c(ix$tv, ix$pb)])),
      dbr_anchored = max(abs(a$mu[ix$br] - br_anchored(tp, 0.5))),
      dmu_recursive_drv = max(abs(r$mu - b$mu)),
      dSigma_recursive_drv = max(abs(r$Sigma - b$Sigma)),
      anchored_vs_recursive = max(abs(a$mu[ix$br] - r$mu[ix$br])))
  }
}
res1 <- rbindlist(res1)
print(res1[, lapply(.SD, fmt), by = .(design, path)], row.names = FALSE)
report(max(res1$dR) < 1e-12, 'correlation matrix equals the driver (Hybrid and CO, all paths)',
       sprintf('max abs difference %s', fmt(max(res1$dR))))
report(max(res1$dSigma) < 1e-10, 'covariance matrix equals the driver',
       sprintf('max abs difference %s (entries up to %.0f)', fmt(max(res1$dSigma)),
               max(abs(drv$build(dgp$path_info(designs$Hybrid$drv$t_wk, designs$Hybrid$drv$e,
                                               designs$Hybrid$drv$paths[[1]]), cl0)$Sigma))))
report(max(res1$dcoup_closed) < 1e-12 && max(res1$dcoup_drv) < 1e-12,
       'graded coupling vector = c_bm * g_t (closed form and driver)',
       sprintf('max abs difference %s (closed form), %s (driver)',
               fmt(max(res1$dcoup_closed)), fmt(max(res1$dcoup_drv))))
report(isTRUE(all(res1[design == 'CO' & path == 2]$pre_exposure_coup == 0)),
       'coupling is 0 before first exposure (CO, placebo-first path)',
       sprintf('max abs coupling before exposure %s',
               fmt(max(res1$pre_exposure_coup, na.rm = TRUE))))
report(max(res1$dmu_tvpb) < 1e-12, 'tv and pb means equal the driver',
       sprintf('max abs difference %s', fmt(max(res1$dmu_tvpb))))
report(max(res1$dbr_anchored) < 1e-12,
       'anchored br means equal the closed form mu_last * 2^(-t_sd/t_half)',
       sprintf('max abs difference %s', fmt(max(res1$dbr_anchored))))
report(max(res1$dmu_recursive_drv) < 1e-12 && max(res1$dSigma_recursive_drv) < 1e-10,
       "with carryover = 'recursive' the means and covariance equal the driver's",
       sprintf('max abs difference %s (means), %s (covariance)',
               fmt(max(res1$dmu_recursive_drv)), fmt(max(res1$dSigma_recursive_drv))))
cat(sprintf('  note: anchored minus recursive br mean, max abs over paths: %s (Hybrid), %s (CO)\n',
            fmt(max(res1[design == 'Hybrid']$anchored_vs_recursive)),
            fmt(max(res1[design == 'CO']$anchored_vs_recursive))))

## The driver's compound-symmetry curve (construct CS, rho 0.8) is the
## published structure with graded coupling
cl_cs <- drv$cell('curve-CS', 0.2, construct = 'CS', rho = 0.8)
d_cs <- max(sapply(seq_along(drv$PIS), function(k) {
  b <- drv$build(drv$PIS[[k]], cl_cs)
  a <- capture(core, td_hyb$trialpaths[[k]], mp(0.2, 0.5), corstructure = 'published',
               carryover = 'recursive')
  max(abs(cov2cor(a$Sigma) - b$R), abs(a$mu - b$mu))
}))
report(d_cs < 1e-12,
       "corstructure = 'published' with graded coupling equals the driver's CS curve (rho 0.8)",
       sprintf('max abs difference %s', fmt(d_cs)))

## Ceiling and the positive-definite stop
## Ceiling (u' M^-1 u)^(-1/2), with u the coupling profile g_t read off a
## feasible matrix (c_bm = 0.1) and M the rest of the correlation matrix
ceil_path <- function(R, cbm) {
  u <- R[1, -1] / cbm
  M <- R[-1, -1]
  1 / sqrt(drop(t(u) %*% solve(M, u)))
}
ceil <- min(sapply(td_hyb$trialpaths, function(tp)
  ceil_path(cov2cor(capture(core, tp, mp_refined(0.1, 0.5))$Sigma), 0.1)))
stops <- function(cbm) inherits(tryCatch({
  for (tp in td_hyb$trialpaths) core$generateData(mp_refined(cbm, 0.5), rp, bp, tp, FALSE, TRUE)
  NULL
}, error = function(e) e), 'error')
drv_min_eig <- function(cbm) {
  cl <- drv$cell('curve-C07', cbm)
  min(sapply(drv$PIS, function(pi) min(eigen(drv$build(pi, cl)$R, TRUE, TRUE)$values)))
}
below <- ceil - 1e-3; above <- ceil + 1e-3
report(abs(ceil - 0.518) < 5e-4, 'Hybrid ceiling of the base cell is 0.518',
       sprintf('ceiling %.4f', ceil))
report(!stops(below) && stops(above) && drv_min_eig(below) > 0 && drv_min_eig(above) <= 0,
       'generator stops above the ceiling and not below it, as the driver',
       sprintf('c_bm %.4f: stop %s, driver min eigenvalue %s; c_bm %.4f: stop %s, driver min eigenvalue %s',
               below, stops(below), fmt(drv_min_eig(below)), above, stops(above),
               fmt(drv_min_eig(above))))
msg <- tryCatch({ core$generateData(mp_refined(0.9, 0.5), rp, bp, td_hyb$trialpaths[[1]],
                                    FALSE, TRUE); '' }, error = conditionMessage)
cat(sprintf('  stop message at c_bm 0.9: %s\n', msg))
core$.repairs <- 0L
invisible(capture(core, td_hyb$trialpaths[[1]], mp_refined(0.9, 0.5), nonPD = 'repair'))
report(nzchar(msg) && core$.repairs == 1L,
       "c_bm 0.9 stops by default; nonPD = 'repair' restores the published repair",
       sprintf('repairs %d', core$.repairs))

## ---------------- 2. Monte Carlo check of the patched sampler ----------------
cat('\n## Large-sample draw from the patched generator (base cell, anchored carryover)\n')
n_mc <- if (fast) 5000L else 20000L
set.seed(20261008)
mc <- rbindlist(lapply(seq_along(td_hyb$trialpaths), function(k) {
  tp <- td_hyb$trialpaths[[k]]
  dat <- core$generateData(mp_refined(0.3, 0.5, N = n_mc), rp, bp, tp, FALSE, TRUE)
  brn <- paste(tp$timeptnames, 'br', sep = '.')
  r <- sapply(brn, function(v) cor(dat$bm, dat[[v]]))
  target <- 0.3 * g_graded(tp, 0.5)
  m <- sapply(brn, function(v) mean(dat[[v]]))
  data.table(path = k, visit = tp$timeptnames, r = r, target = target,
             z_cor = (r - target) / ((1 - target^2) / sqrt(n_mc)),
             mean = m, anchored = br_anchored(tp, 0.5),
             z_mean = (m - br_anchored(tp, 0.5)) / (rp[cat == 'br']$sd / sqrt(n_mc)))
}))
print(mc[, .(path, visit, r = round(r, 4), target = round(target, 4), z_cor = round(z_cor, 2),
             mean = round(mean, 3), anchored = round(anchored, 3), z_mean = round(z_mean, 2))],
      row.names = FALSE)
report(max(abs(mc$z_cor)) < 4, sprintf('empirical Cor(bm, BR_t) = c_bm g_t (N = %d per path)', n_mc),
       sprintf('max |z| %.2f over %d visits', max(abs(mc$z_cor)), nrow(mc)))
report(max(abs(mc$z_mean)) < 4, 'empirical BR means = anchored closed form',
       sprintf('max |z| %.2f over %d visits', max(abs(mc$z_mean)), nrow(mc)))

## ---------------- 3. randomized allocation ----------------
cat('\n## Path allocation in generateSimulatedResults()\n')
alloc_run <- function(allocation, nreps) {
  env <- make_env('core')
  counts <- list()
  env$lme_analysis <- function(td, dat, op) {
    if (uniqueN(dat$replicate) == 1)
      counts[[length(counts) + 1]] <<- tabulate(dat$path, length(td))
    data.table(beta = 0, betaSE = 1, p = 0.5, issingular = FALSE, warning = NA_character_)
  }
  env$generateSimulatedResults(
    trialdesigns = list(td_hyb), respparamsets = ps$respparamsets[1],
    blparamsets = ps$blparamsets[1], censorparams = NA,
    modelparams = mp_refined(0.3, 0.5, N = 70L),
    simparam = list(Nreps = nreps, progressiveSave = FALSE, saveunit2start = 1),
    analysisparams = expand.grid(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE),
    allocation = allocation)
  do.call(rbind, counts)
}
n_alloc <- if (fast) 300L else 1000L
set.seed(1)
cr <- alloc_run('randomized', n_alloc)
cf <- alloc_run('fixed', 5L)
share <- colSums(cr) / sum(cr)
v <- apply(cr, 2, var)
report(nrow(cr) == n_alloc && all(rowSums(cr) == 70L) && all(cr > 0),
       'randomized: every replicate has N = 70 and no empty path',
       sprintf('%d replicates', nrow(cr)))
report(max(abs(share - 0.25)) < 4 * sqrt(0.25 * 0.75 / sum(cr)),
       'randomized: path shares are 1/4',
       sprintf('shares %s (SE %.4f)', paste(round(share, 4), collapse = ', '),
               sqrt(0.25 * 0.75 / sum(cr))))
report(all(abs(v / (70 * 0.25 * 0.75) - 1) < 4 * sqrt(2 / (n_alloc - 1))),
       'randomized: per-replicate counts vary as Binomial(70, 1/4)',
       sprintf('variance %s against %.3f; distinct count vectors %d',
               paste(round(v, 2), collapse = ', '), 70 * 0.25 * 0.75, nrow(unique(cr))))
report(all(apply(cf, 1, identical, c(18L, 18L, 17L, 17L))),
       "fixed: the published split 18, 18, 17, 17 in every replicate")

## ---------------- 4. published behavior is reproduced ----------------
cat('\n## Published settings: patched trees against unpatched orig\n')
pub_opts <- list(coupling = 'step', corstructure = 'published', carryover = 'recursive',
                 nonPD = 'repair')
rp_pub <- ps$respparamsets[[1]]$param
bp_pub <- ps$blparamsets[[1]]$param
gd_pub <- function(env, m, tp) do.call(env$generateData,
                                       c(list(m, rp_pub, bp_pub, tp, FALSE, TRUE), pub_opts))
cells_pub <- data.table(cbm = c(0, 0.25, 0.3, 0.6, 0.3), th = c(0, 0, 0, 0.1, 0.2))
same <- TRUE; reps <- integer(0)
for (i in seq_len(nrow(cells_pub))) for (k in 1:4) {
  m <- mp(cells_pub$cbm[i], cells_pub$th[i], N = 50L)
  tp <- td_hyb$trialpaths[[k]]
  orig$.repairs <- 0L
  set.seed(100 * i + k); a <- orig$generateData(m, rp_pub, bp_pub, tp, FALSE, TRUE)
  reps <- c(reps, orig$.repairs)
  for (env in list(core, opts, full)) {
    set.seed(100 * i + k); b <- gd_pub(env, m, tp)
    same <- same && identical(a, b)
  }
}
report(same, 'generateData(): published options reproduce orig bit for bit (core, options, full)',
       sprintf('%d cells x 4 paths', nrow(cells_pub)))
rep_by <- cells_pub[, .(cbm, th, repairs_over_paths = colSums(matrix(reps, nrow = 4)))]
print(rep_by, row.names = FALSE)
report(rep_by[cbm == 0.25 & th == 0]$repairs_over_paths == 0 &&
         rep_by[cbm == 0.3 & th == 0]$repairs_over_paths > 0,
       'orig repairs above the compound-symmetry ceiling 0.256 (t1/2 = 0) and not below')

gsr <- function(env, extra = list(), ap = expand.grid(useDE = FALSE, t_random_slope = FALSE,
                                                      full_model_out = FALSE)) {
  set.seed(4242)
  do.call(env$generateSimulatedResults, c(list(
    trialdesigns = list(td_hyb), respparamsets = ps$respparamsets[1],
    blparamsets = ps$blparamsets[1], censorparams = ps$censorparams,
    modelparams = mp(0.3, 0.1, N = 35L),
    simparam = list(Nreps = 3L, progressiveSave = FALSE, saveunit2start = 1),
    analysisparams = ap), extra))$results
}
ap_orig <- data.frame(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE,
                      analysis = 'orig')
g0 <- suppressMessages(suppressWarnings(gsr(orig)))
g1 <- suppressMessages(suppressWarnings(gsr(full, list(dgpoptions = pub_opts, allocation = 'fixed'),
                                             ap_orig)))
g2 <- suppressMessages(suppressWarnings(gsr(core, list(dgpoptions = pub_opts, allocation = 'fixed'))))
cols <- c('censorparamset', 'irep', 'frac_NA', 'ETbeta', 'ETbetaSE', 'beta', 'betaSE', 'p',
          'issingular')
report(identical(g0[, ..cols], g1[, ..cols]) && identical(g0[, ..cols], g2[, ..cols]),
       'generateSimulatedResults(): published DGP, fixed allocation and analysis orig reproduce orig',
       sprintf('%d result rows (3 replicates x 5 censoring sets), beta and p identical', nrow(g0)))

g3 <- suppressMessages(suppressWarnings(gsr(full, list(), data.frame(
  useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE,
  analysis = c('F1', 'S', 'orig')))))
report(nrow(g3) == 45L && !anyNA(g3$beta) && setequal(g3$analysis, c('F1', 'S', 'orig')) &&
         all(is.na(g3[analysis == 'orig']$df)) && !anyNA(g3[analysis != 'orig']$df),
       'generateSimulatedResults() runs end to end with the refined defaults (smoke test)',
       sprintf('%d rows; analyses %s; failed fits %d; mean df F1 %.1f', nrow(g3),
               paste(unique(g3$analysis), collapse = ', '), sum(is.na(g3$beta)),
               mean(g3[analysis == 'F1']$df)))

## ---------------- 4b. the pooled (whole-population) fit ----------------
cat('\n## Pooled fit (ETbeta): population_fit\n')
et_orig <- g1[analysis == 'orig']
report(identical(g0$ETbeta, et_orig$ETbeta) && identical(g0$ETbetaSE, et_orig$ETbetaSE) &&
         all(et_orig$ETbetaSE_type == 'full'),
       "analysis 'orig' with the default population_fit = 'full' gives orig's ETbeta and ETbetaSE",
       sprintf('ETbeta %s', fmt(unique(g0$ETbeta))))
et_f1 <- g3[analysis == 'F1']; et_o <- g3[analysis == 'orig']
report(!anyNA(et_f1$ETbeta) && all(et_f1$ETbetaSE_type == 'model') &&
         all(et_o$ETbetaSE_type == 'full'),
       "refined analyses default to population_fit = 'model'; orig to 'full' (end-to-end run)",
       sprintf('ETbeta F1 %s (SE %s), S %s, orig %s', fmt(unique(et_f1$ETbeta)),
               fmt(unique(et_f1$ETbetaSE)), fmt(unique(g3[analysis == 'S']$ETbeta)),
               fmt(unique(et_o$ETbeta))))
g4 <- suppressMessages(suppressWarnings(gsr(full, list(), data.frame(
  useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE, analysis = c('F1', 'R'),
  population_fit = c(FALSE, NA)))))
report(nrow(g4) == 30L && all(is.na(g4$ETbeta)) && all(is.na(g4$ETbetaSE)) &&
         all(is.na(g4$ETbetaSE_type)) && !anyNA(g4$beta),
       'population_fit = FALSE skips the pooled fit; the mmrm analysis R defaults to FALSE',
       sprintf('%d rows (F1 with FALSE, R with the default)', nrow(g4)))

make_pool <- function(nreps, seed) {
  set.seed(seed)
  rbindlist(lapply(1:4, function(k) {
    x <- full$generateData(mp_refined(0.3, 0.5, N = 18L * nreps), rp, bp, td_hyb$trialpaths[[k]],
                           FALSE, TRUE)
    x[, `:=`(path = k, replicate = rep(seq_len(nreps), 18L))]
  }))
}
pool <- make_pool(20L, 31)
pf <- rbindlist(lapply(c('F1', 'C1', 'DS', 'A', 'R'), function(a) {
  t0 <- Sys.time()
  m <- full$lme_analysis(td_hyb$trialpaths, pool, list(analysis = a, inference = 'model'))
  t1 <- Sys.time()
  f <- full$lme_analysis(td_hyb$trialpaths, pool, list(analysis = a))
  t2 <- Sys.time()
  data.table(analysis = a, beta_model = m$beta, beta_full = f$beta, d = abs(m$beta - f$beta),
             se_model = m$betaSE, se_full = f$betaSE, h_model = m$t_half_a, h_full = f$t_half_a,
             sec_model = as.numeric(difftime(t1, t0, units = 'secs')),
             sec_full = as.numeric(difftime(t2, t1, units = 'secs')))
}))
print(pf[, .(analysis, beta_model = signif(beta_model, 8), d = fmt(d), se_model = signif(se_model, 4),
             se_full = signif(se_full, 4), h_model, h_full, sec_model = round(sec_model, 1),
             sec_full = round(sec_full, 1))], row.names = FALSE)
for (i in seq_len(nrow(pf)))
  report(!is.na(pf$d[i]) && pf$d[i] < 1e-8 && pf$h_model[i] == pf$h_full[i],
         sprintf("pooled fit, analysis %s: 'model' point estimate equals 'full' (1,440 participants)",
                 pf$analysis[i]),
         sprintf('beta %s, abs difference %s; SE model %s, full %s', fmt(pf$beta_model[i]),
                 fmt(pf$d[i]), fmt(pf$se_model[i]), fmt(pf$se_full[i])))
pool100 <- make_pool(if (fast) 20L else 100L, 32)
t0 <- Sys.time()
m100 <- full$lme_analysis(td_hyb$trialpaths, pool100, list(analysis = 'F1', inference = 'model'))
sec100 <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
report(!is.na(m100$beta) && sec100 < 300,
       sprintf("pooled 'model' fit of F1 on %d participants finishes within 5 minutes",
               nrow(pool100)),
       sprintf('%.1f s; beta %s', sec100, fmt(m100$beta)))

## ---------------- 5. analyses against the driver's code ----------------
cat('\n## Analyses: patched lme_analysis() against the driver code on the same data\n')
an_codes <- c('S', 'DS', 'C1', 'C2', 'F1', 'F1A', 'F2', 'ES', 'G1', 'G2', 'A', 'DA', 'EA', 'R')
tpn <- td_hyb$trialpaths[[1]]$timeptnames
make_dat <- function(seed, dropout = NA_character_) {
  set.seed(seed)
  Ns <- c(18L, 18L, 17L, 17L)
  dat <- rbindlist(lapply(1:4, function(k) {
    x <- full$generateData(mp_refined(0.3, 0.5, N = Ns[k]), rp, bp, td_hyb$trialpaths[[k]],
                           FALSE, TRUE)
    x[, `:=`(path = k, replicate = 1L)]
  }))
  if (!is.na(dropout))
    dat <- full$censordata(dat, td_hyb$trialpaths[[1]],
                           data.frame(names = dropout, beta0 = 0.2, beta1 = 0, eb1 = 2,
                                      type = dropout))
  dat
}
## The driver's make_trial() assembles its long data from draw() and
## dropout_visit(); here they hand it the same participants instead, and the
## dropout visit is read off the missing values.
driver_long <- function(dat) {
  fake <- split(dat, by = 'path', sorted = TRUE)
  k <- 0L
  drv$draw <- function(b, pi, cl, n_draw) {
    k <<- k + 1L
    x <- fake[[k]]
    stopifnot(nrow(x) == n_draw)
    list(bm = x$bm, BL = x$BL, Y = as.matrix(x[, ..tpn]))
  }
  drv$dropout_visit <- function(Y, type) apply(Y, 1, function(y) {
    j <- which(is.na(y)); if (length(j)) j[1] else Inf })
  cl <- drv$cell('curve-C07', 0.3)
  drv$make_trial(cl, vector('list', 4), as.integer(table(dat$path)), 1)
}
datasets <- list(list(seed = 11, dropout = NA_character_), list(seed = 12, dropout = NA_character_),
                 list(seed = 13, dropout = 'MAR'))
if (!fast) datasets <- c(datasets, list(list(seed = 14, dropout = 'MCAR'),
                                        list(seed = 15, dropout = NA_character_)))
cmp <- list(); orig_same <- TRUE
for (ds in datasets) {
  dat <- make_dat(ds$seed, ds$dropout)
  dl <- driver_long(dat)
  dr <- rbind(drv$fit_all(dl), drv$fit_extra(dl), drv$fit_extra2(dl))
  dr <- dr[!duplicated(analysis)]
  for (a in an_codes) {
    pt <- full$lme_analysis(td_hyb$trialpaths, dat, list(analysis = a, full_model_out = FALSE))
    d1 <- dr[analysis == a]
    cmp[[length(cmp) + 1]] <- data.table(
      seed = ds$seed, dropout = ds$dropout, analysis = a,
      est = pt$beta, se = pt$betaSE, p = pt$p, h = pt$t_half_a, df = pt$df,
      d_est = abs(pt$beta - d1$est), d_se = abs(pt$betaSE - d1$se), d_p = abs(pt$p - d1$p),
      d_h = abs(pt$t_half_a - d1$h), both_na = is.na(pt$beta) && is.na(d1$est))
  }
  for (de in c(FALSE, TRUE)) {
    op <- list(useDE = de, t_random_slope = FALSE, full_model_out = FALSE)
    a <- suppressMessages(orig$lme_analysis(td_hyb$trialpaths, dat, op))
    b <- suppressMessages(full$lme_analysis(td_hyb$trialpaths, dat, c(op, analysis = 'orig')))
    orig_same <- orig_same && identical(a, b)
  }
}
cmp <- rbindlist(cmp)
summ <- cmp[, .(datasets = .N, fails = sum(is.na(est)),
                max_d_est = max(d_est, na.rm = TRUE), max_d_se = max(d_se, na.rm = TRUE),
                max_d_p = max(d_p, na.rm = TRUE), max_d_h = max(d_h, na.rm = TRUE),
                mean_df = mean(df, na.rm = TRUE)), by = analysis]
print(summ[, .(analysis, datasets, fails, max_d_est = fmt(max_d_est), max_d_se = fmt(max_d_se),
               max_d_p = fmt(max_d_p), max_d_h = fmt(max_d_h), mean_df = round(mean_df, 1))],
      row.names = FALSE)
for (a in an_codes) {
  s <- summ[analysis == a]; x <- cmp[analysis == a]
  ok <- all(x$both_na | (!is.na(x$d_est) & x$d_est < 1e-8 & x$d_se < 1e-8 & x$d_p < 1e-8 &
                           x$d_h == 0))
  report(ok, sprintf('analysis %s equals the driver on %d datasets', a, s$datasets),
         sprintf('max abs difference est %s, SE %s, p %s', fmt(s$max_d_est), fmt(s$max_d_se),
                 fmt(s$max_d_p)))
}
report(orig_same, "analysis 'orig' is identical to unpatched lme_analysis() (useDE FALSE and TRUE)",
       sprintf('%d datasets', length(datasets)))

## Pre-exposure guard of the exposure-decayed indicator
tp2 <- td_co$trialpaths[[2]]
w_guard <- full$exposuredecay(tp2$tod > 0, tp2$tsd, 0.5)
w_naive <- ifelse(tp2$tod > 0, 1, 0.5^(tp2$tsd / 0.5))
set.seed(3)
dco <- rbindlist(lapply(1:2, function(k) full$generateData(mp_refined(0.3, 0.5, N = 35L), rp, bp,
                                                          td_co$trialpaths[[k]], FALSE, TRUE)[
  , `:=`(path = k, replicate = 1L)]))
fm <- full$lme_analysis(td_co$trialpaths, dco, list(analysis = 'F1', full_model_out = TRUE))
dm <- fm$datamerged
pre <- dm[as.integer(as.character(ptID)) > 35 & t > 0 & t <= 10]
report(all(w_guard[1:4] == 0) && all(pre$Dbc == 0) && all(w_naive[1:4] == 1),
       'Dbc = 0 before first exposure (CO, placebo-first path)',
       sprintf('Dbc by visit %s; unguarded 0.5^(tsd/h) would give %s; %d long-data rows checked',
               paste(round(w_guard, 3), collapse = ' '), paste(round(w_naive, 3), collapse = ' '),
               nrow(pre)))

## ---------------- 6. options against the driver ----------------
cat('\n## Options patch\n')
opt_cells <- list(
  list(label = 'prog-TV-const', drv = list(c_tv = 0.3), pat = list(c.bm.tv = 0.3)),
  list(label = 'prog-PB', drv = list(c_pb = 0.3), pat = list(c.bm.pb = 0.3)),
  list(label = 'prog-BL', drv = list(c_bl = 0.3), pat = list(c.bm.bl = 0.3)),
  list(label = 'weibull0.5', drv = list(decay_k = 0.5), pat = list(decay.k = 0.5)),
  list(label = 'weibull2', drv = list(decay_k = 2), pat = list(decay.k = 2)),
  list(label = 'thalf0', drv = list(t_half = 0), pat = list(), th = 0),
  list(label = 'thalf1', drv = list(t_half = 1), pat = list(), th = 1))
res6 <- rbindlist(lapply(opt_cells, function(oc) {
  cl <- do.call(drv$cell, c(list(oc$label, 0.3), oc$drv))
  th <- if (is.null(oc$th)) 0.5 else oc$th
  rbindlist(lapply(seq_along(drv$PIS), function(k) {
    b <- drv$build(drv$PIS[[k]], cl)
    m <- do.call(mp_refined, c(list(0.3, th), oc$pat))
    a <- capture(opts, td_hyb$trialpaths[[k]], m, carryover = 'recursive')
    data.table(cell = oc$label, path = k, dR = max(abs(cov2cor(a$Sigma) - b$R)),
               dSigma = max(abs(a$Sigma - b$Sigma)), dmu = max(abs(a$mu - b$mu)))
  }))
}))
summ6 <- res6[, .(dR = max(dR), dSigma = max(dSigma), dmu = max(dmu)), by = cell]
print(summ6[, .(cell, dR = fmt(dR), dSigma = fmt(dSigma), dmu = fmt(dmu))], row.names = FALSE)
for (i in seq_len(nrow(summ6)))
  report(summ6$dR[i] < 1e-12 && summ6$dSigma[i] < 1e-10 && summ6$dmu[i] < 1e-12,
         sprintf("option cell %s: matrix and means (carryover = 'recursive') equal the driver",
                 summ6$cell[i]),
         sprintf('max abs difference R %s, Sigma %s, mu %s', fmt(summ6$dR[i]),
                 fmt(summ6$dSigma[i]), fmt(summ6$dmu[i])))

## Mean moderation of tv: same draw with and without, so the difference is the shift
tp1 <- td_hyb$trialpaths[[1]]
set.seed(77); x0 <- opts$generateData(mp_refined(0, 0.5, N = 2000L), rp, bp, tp1, FALSE, TRUE)
set.seed(77); x1 <- opts$generateData(mp_refined(0, 0.5, N = 2000L, b.bm.tv = 0.3), rp, bp, tp1,
                                      FALSE, TRUE)
tvn <- paste(tp1$timeptnames, 'tv', sep = '.')
z <- (x0$bm - dgp$bm_m) / dgp$bm_sd
shift_drv <- 0.3 * 10 * outer(z, dgp$gomp(drv$PIS[[1]]$w, 6.50647, 5, 0.35) / drv$TV_MAX)
d_tv <- max(abs(as.matrix(x1[, ..tvn]) - as.matrix(x0[, ..tvn]) - shift_drv))
d_y <- max(abs(as.matrix(x1[, ..tpn]) - as.matrix(x0[, ..tpn]) + shift_drv))
report(d_tv < 1e-10 && d_y < 1e-10, "option b.bm.tv: tv shift equals the driver's formula",
       sprintf('max abs difference %s (tv), %s (observed score)', fmt(d_tv), fmt(d_y)))

## Measurement error: same draw with and without
set.seed(78); y0 <- opts$generateData(mp_refined(0.3, 0.5, N = n_mc), rp, bp, tp1, FALSE, TRUE)
set.seed(78); y1 <- opts$generateData(mp_refined(0.3, 0.5, N = n_mc, me.ol = 3, me.bl = 6), rp, bp,
                                      tp1, FALSE, TRUE)
s_t <- sapply(tpn, function(v) sd(y1[[v]] - y0[[v]]))
s_target <- ifelse(tp1$e == 1, 3, 6)
z_me <- (s_t - s_target) / (s_target / sqrt(2 * n_mc))
other <- setdiff(names(y0), tpn)
comp_same <- identical(y0[, ..other], y1[, ..other])
report(max(abs(z_me)) < 4 && comp_same,
       'option me.ol/me.bl: error SD 3 open label, 6 blinded, on the observed score only',
       sprintf('SDs %s; max |z| %.2f; components and D_ columns unchanged %s',
               paste(round(s_t, 3), collapse = ' '), max(abs(z_me)), comp_same))

## Dropout: one path of one trial, same seed, against the driver's dropout_visit()
first_na <- function(dat) apply(as.matrix(dat[, ..tpn]), 1, function(y) {
  j <- which(is.na(y)); if (length(j)) j[1] else Inf })
same_drop <- TRUE; mono <- TRUE
for (ty in c('MCAR', 'MAR')) for (s in 1:5) {
  set.seed(500 + s)
  x <- opts$generateData(mp_refined(0.3, 0.5, N = 70L), rp, bp, tp1, FALSE, TRUE)
  x[, `:=`(path = 1L, replicate = 1L)]
  Y <- as.matrix(x[, ..tpn])
  set.seed(900 + s); dv <- drv_dropout_visit(Y, ty)
  set.seed(900 + s)
  xc <- opts$censordata(x, tp1, data.frame(names = ty, beta0 = 0.2, beta1 = 0, eb1 = 2, type = ty))
  fn <- first_na(xc)
  same_drop <- same_drop && identical(as.numeric(fn), as.numeric(dv))
  nas <- is.na(as.matrix(xc[, ..tpn]))
  mono <- mono && all(apply(nas, 1, function(r) all(diff(r) >= 0))) && !any(nas[, 1:2])
}
report(same_drop, "option type = 'MCAR'/'MAR': dropout visits identical to the driver's (same seed)",
       '5 seeds each, N = 70')
report(mono, 'dropout is monotone and starts at the third timepoint')
nr <- if (fast) 50L else 200L
set.seed(6)
big <- rbindlist(lapply(1:4, function(k) {
  x <- opts$generateData(mp_refined(0.3, 0.5, N = 18L * nr), rp, bp, td_hyb$trialpaths[[k]],
                         FALSE, TRUE)
  x[, `:=`(path = k, replicate = rep(seq_len(nr), 18L))]
}))
rates <- sapply(c('MCAR', 'MAR'), function(ty) {
  xc <- opts$censordata(big, td_hyb$trialpaths[[1]],
                        data.frame(names = ty, beta0 = 0.2, beta1 = 0, eb1 = 2, type = ty))
  mean(is.finite(first_na(xc)))
})
se_rate <- sqrt(0.2 * 0.8 / nrow(big))
report(all(abs(rates - 0.2) < 4 * se_rate), 'dropout rates are 20%',
       sprintf('MCAR %.4f, MAR %.4f (SE %.4f, %d participants)', rates[1], rates[2], se_rate,
               nrow(big)))

cat(sprintf('\n%s: %d check(s) failed\n', if (n_fail == 0) 'ALL PASS' else 'FAILURES', n_fail))
quit(status = if (n_fail == 0) 0 else 1)
