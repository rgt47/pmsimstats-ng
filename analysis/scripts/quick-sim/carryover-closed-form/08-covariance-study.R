## Monte Carlo study: does the working covariance of the analysis, and
## its interplay with the time adjustment, determine the size and power
## of the biomarker-by-drug interaction test? Companion to docs/40,
## Sections 8 and 9.
##
## Data (04-power-simulation.R construct, published means and SDs,
## published allocation 18/18/17/17 and 35/35, N = 70):
##   CS   published compound symmetry (rho = 0.8 at every lag, cross-
##        factor 0.2 / 0.1), graded coupling, strength 0.25
##   C07  separable AR(1) (configuration C), rho = 0.7 per week,
##        graded coupling, strength 0.45 (ceiling 0.474 or more in every
##        cell: Hybrid 0.480, CO 0.616 at t1/2 = 0)
## each with a null twin (strength 0), designs Hybrid and CO,
## t1/2 = 0 and 1 week.
##
## Analyses, all mmrm on the same replicate (biomarker centered):
##   mean   linear t | ns(t, 3) | visit factor, each with bmc + Db + bmc:Db
##   cov    cs | csh | us (visit factor) | sp_exp (continuous time)
##   test   Satterthwaite, asymptotic vcov; plus Kenward-Roger
##          (Kenward-Roger-Linear vcov) for us with linear t and with
##          visit factor
##   coding baseline as a ninth response row (published), and baseline
##          as a covariate on the 8 post-baseline rows (cs and us, each
##          mean form; docs/41)
## cs with linear t, baseline as a row, is the published random-intercept
## model.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/08-covariance-study.R \
##     [--reps-alt R] [--reps-null R] [--cores K] [--chunk C] [--arms main|followup]
## Replicates are run in chunks of C (default 50), each saved on
## completion; a relaunch with the same arguments resumes.
suppressMessages({
  library(data.table); library(MASS); library(mmrm); library(splines); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
reps_alt <- arg('--reps-alt', 500L)
reps_null <- arg('--reps-null', 1000L)
cores <- arg('--cores', 8L)
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/covariance-study'
arm_arg <- if ('--arms' %in% args) args[which(args == '--arms') + 1] else 'main'
arm_tag <- if (arm_arg == 'main') '' else paste0('-', arm_arg)
run_tag <- sprintf('alt%d-null%d%s', reps_alt, reps_null, arm_tag)
cell_dir <- file.path(out_dir, paste0('cells-', run_tag))
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build', 'draw_path')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))
N_TOTAL <- dgp$N_TOTAL

## Graded coupling under CS is not an 04 arm: the arm A matrix with the
## graded biomarker row (as in 06 and 07). Configuration C is arm C.
build_cell <- function(pi, construct, cbm, th) {
  if (construct == 'C07') return(dgp$build(pi, 'C', cbm, th))
  b <- dgp$build(pi, 'A', 0, th)
  n <- b$n
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, n)
  phi[pi$on | pi$tsd <= 0] <- 0
  idx <- (3 + 2 * n):(2 + 3 * n)
  b$R[1, idx] <- b$R[idx, 1] <- cbm * ifelse(pi$on, 1, phi)
  sds <- sqrt(diag(b$Sigma))
  b$Sigma <- outer(sds, sds) * b$R
  b
}

MEANS <- c(linear = 'Sx ~ bmc + Db + t + bmc:Db',
           spline = 'Sx ~ bmc + Db + ns(t, df = 3) + bmc:Db',
           visit = 'Sx ~ bmc + Db + visit + bmc:Db')
COVS <- c(cs = 'cs(visit | ptID)', csh = 'csh(visit | ptID)', us = 'us(visit | ptID)',
          sp_exp = 'sp_exp(t | ptID)')
## coding 'row': baseline is a ninth row of the response (published).
## coding 'covariate': post-baseline rows only, baseline as covariate BL
## (docs/41), with cs and us covariance over the 8 post-baseline visits.
## --arms main: the first run's 20 analyses. --arms followup: Kenward-
## Roger and CR2 arms on the same replicates (same cells and seeds), plus
## one check of whether the row coding's power under AR(1) comes from
## the biomarker main effect it shares with the baseline row ('sepbm'
## frees that effect at baseline).
ARMS <- if ('--arms' %in% args) args[which(args == '--arms') + 1] else 'main'
stopifnot(ARMS %in% c('main', 'followup', 'dbc'))
## --arms dbc: the exposure-weighted drug indicator Dbc (1 on drug,
## 2^(-tsd / h) off drug after exposure, 0 before any exposure and at
## baseline) at analysis half-lives h = 0.5, 1, 2 weeks, under three
## inferences: us with Kenward-Roger, us with CR2 (mmrm), and the
## random intercept with corCAR1 and CR2 (nlme + clubSandwich; manuscript
## 02's G8). Binary Db under the last inference is the reference (02 G6).
## Baseline stays a response row throughout.
DBC_H <- c(0.5, 1, 2)
ANALYSES <- if (ARMS == 'dbc') {
  rbind(
    CJ(coding = 'row', mean = c('linear', 'visit'), cov = 'us',
       test = c('Kenward-Roger', 'CR2'), expo = 'Dbc', th_a = DBC_H, engine = 'mmrm',
       sorted = FALSE),
    CJ(coding = 'row', mean = c('linear', 'visit'), cov = 'car1', test = 'CR2',
       expo = 'Dbc', th_a = DBC_H, engine = 'lme', sorted = FALSE),
    CJ(coding = 'row', mean = c('linear', 'visit'), cov = 'car1', test = 'CR2',
       expo = 'Db', th_a = NA_real_, engine = 'lme', sorted = FALSE))[, extra := 'none']
} else if (ARMS == 'main') {
  rbind(
    CJ(coding = 'row', mean = names(MEANS), cov = names(COVS), test = 'Satterthwaite',
       sorted = FALSE),
    data.table(coding = 'row', mean = c('linear', 'visit'), cov = 'us', test = 'Kenward-Roger'),
    CJ(coding = 'covariate', mean = names(MEANS), cov = c('cs', 'us'), test = 'Satterthwaite',
       sorted = FALSE))[, extra := 'none']
} else {
  rbind(
    data.table(coding = 'row', mean = 'spline', cov = 'us', test = 'Kenward-Roger'),
    CJ(coding = 'row', mean = c('linear', 'visit'), cov = c('cs', 'csh'),
       test = 'Kenward-Roger', sorted = FALSE),
    CJ(coding = 'covariate', mean = names(MEANS), cov = 'us', test = 'Kenward-Roger',
       sorted = FALSE),
    CJ(coding = 'covariate', mean = c('linear', 'visit'), cov = 'cs', test = 'Kenward-Roger',
       sorted = FALSE),
    CJ(coding = 'row', mean = c('linear', 'visit'), cov = names(COVS), test = 'CR2',
       sorted = FALSE),
    CJ(coding = 'covariate', mean = c('linear', 'visit'), cov = c('cs', 'us'), test = 'CR2',
       sorted = FALSE))[, extra := 'none'][
    , rbind(.SD, data.table(coding = 'row', mean = 'visit', cov = 'us',
                            test = 'Kenward-Roger', extra = 'sepbm'))]
}
if (!'expo' %in% names(ANALYSES)) ANALYSES[, `:=`(expo = 'Db', th_a = NA_real_, engine = 'mmrm')]
ANALYSES[, label := paste(coding, cov, mean, test, extra, expo, th_a, engine, sep = '|')]

## Exposure-weighted indicator at analysis half-life h. Off-drug visits
## before any exposure (and baseline) have tsd = 0 and are coded 0, not
## exp(0) = 1.
make_dbc <- function(Db, tsd, h) ifelse(Db, 1, ifelse(tsd > 0, 0.5^(tsd / h), 0))

fit_lme_cr2 <- function(f, long, target) {
  fit <- tryCatch(
    nlme::lme(f, random = ~ 1 | ptID, correlation = nlme::corCAR1(form = ~ t | ptID),
              data = long, control = nlme::lmeControl(opt = 'optim', maxIter = 200,
                                                      msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(c(est = NA_real_, se = NA_real_, p = NA_real_))
  ct <- clubSandwich::coef_test(fit, vcov = 'CR2', test = 'Satterthwaite')
  k <- which(ct$Coef == target)
  c(est = ct$beta[k], se = ct$SE[k], p = ct$p_Satt[k])
}

fit_one <- function(long, a) {
  mean_f <- MEANS[[a$mean]]
  if (a$coding == 'covariate') {
    long <- long[t > 0]
    long[, visit := droplevels(visit)]
    mean_f <- sub('Sx ~ ', 'Sx ~ BL + ', mean_f, fixed = TRUE)
  }
  if (a$extra == 'sepbm') mean_f <- paste(mean_f, '+ bmc:bl0')
  target <- 'bmc:DbTRUE'
  if (a$expo == 'Dbc') {
    long <- copy(long)
    long[, Dbc := make_dbc(Db, tsd, a$th_a)]
    mean_f <- gsub('\\bDb\\b', 'Dbc', mean_f)
    target <- 'bmc:Dbc'
  }
  if (a$engine == 'lme') return(fit_lme_cr2(as.formula(mean_f), long, target))
  f <- as.formula(paste(mean_f, '+', COVS[[a$cov]]))
  tryCatch({
    fit <- suppressWarnings(switch(a$test,
      'Kenward-Roger' = mmrm(f, data = long, method = 'Kenward-Roger',
                             vcov = if (a$cov == 'us') 'Kenward-Roger-Linear' else 'Kenward-Roger'),
      'CR2' = mmrm(f, data = long, method = 'Satterthwaite', vcov = 'Empirical-Bias-Reduced'),
      mmrm(f, data = long)))
    cf <- summary(fit)$coefficients
    c(est = cf[target, 'Estimate'], se = cf[target, 'Std. Error'],
      p = cf[target, 'Pr(>|t|)'])
  }, error = function(e) c(est = NA_real_, se = NA_real_, p = NA_real_))
}

## Work is split into chunks of CHUNK replicates, each with its own seed
## and saved as it completes, so an interrupted run loses at most one
## chunk per worker and a relaunch resumes.
CHUNK <- arg('--chunk', 50L)
run_chunk <- function(job) {
  set.seed(job$seed)
  assign('RHO', job$rho, envir = dgp)
  d <- dgp$designs[[job$design]]
  pis <- lapply(d$paths, function(od) dgp$path_info(d$t_wk, d$e, od))
  nP <- length(pis)
  Ns <- N_TOTAL %/% nP + c(rep(1L, N_TOTAL %% nP), rep(0L, nP - N_TOTAL %% nP))
  bs <- lapply(pis, build_cell, construct = job$construct, cbm = job$cbm, th = job$t_half)
  min_eig <- min(sapply(bs, function(b) min(eigen(b$R, TRUE, TRUE)$values)))
  stopifnot(min_eig > 0)
  reps <- job$n_reps
  draws <- lapply(seq_len(nP), function(k)
    dgp$draw_path(bs[[k]], pis[[k]], job$cbm, Ns[k] * reps))
  vis_levels <- as.character(c(0, cumsum(d$t_wk)))
  res <- rbindlist(lapply(seq_len(reps), function(r) {
    long <- rbindlist(lapply(seq_len(nP), function(k) {
      rows <- ((r - 1) * Ns[k] + 1):(r * Ns[k])
      dr <- draws[[k]]; pi <- pis[[k]]
      id <- (k - 1) * 1000 + seq_len(Ns[k])
      rbind(data.table(ptID = id, bm = dr$bm[rows], Sx = dr$BL[rows], t = 0, Db = FALSE,
                       tsd = 0),
            data.table(ptID = rep(id, times = length(pi$w)),
                       bm = rep(dr$bm[rows], times = length(pi$w)),
                       Sx = as.vector(dr$Y[rows, , drop = FALSE]),
                       t = rep(pi$w, each = Ns[k]), Db = rep(pi$on, each = Ns[k]),
                       tsd = rep(pi$tsd, each = Ns[k])))
    }))
    long[, `:=`(ptID = factor(ptID), visit = factor(as.character(t), levels = vis_levels),
                bmc = bm - mean(bm[t == 0]))]
    long[, BL := Sx[t == 0], by = ptID]
    long[, bl0 := as.numeric(t == 0)]
    rbindlist(lapply(seq_len(nrow(ANALYSES)), function(j) {
      a <- ANALYSES[j]
      as.data.table(as.list(fit_one(long, a)))[, `:=`(label = a$label, rep = r)]
    }))
  }))
  res[, `:=`(rep = (job$chunk - 1L) * CHUNK + rep, cell = job$id, min_eig = min_eig)]
  merge(res, ANALYSES, by = 'label')
}

run_chunk_checkpointed <- function(job) {
  f <- file.path(cell_dir, sprintf('cell-%03d-chunk-%02d.rds', job$id, job$chunk))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  out <- run_chunk(job)
  saveRDS(out, f)
  cat(sprintf('%s cell %2d chunk %2d %-4s %-6s %-5s t1/2=%.0f  %.1f min\n',
              format(Sys.time(), '%H:%M:%S'), job$id, job$chunk, job$construct, job$design,
              job$hyp, job$t_half, as.numeric(difftime(Sys.time(), t0, units = 'mins'))),
      file = file.path(cell_dir, 'progress.log'), append = TRUE)
  out
}

cells <- CJ(construct = c('CS', 'C07'), design = c('Hybrid', 'CO'), hyp = c('alt', 'null'),
            t_half = c(0, 1), sorted = FALSE)
cells[, `:=`(rho = ifelse(construct == 'CS', 0.8, 0.7),
             cbm = ifelse(hyp == 'null', 0, ifelse(construct == 'CS', 0.25, 0.45)),
             reps = ifelse(hyp == 'null', reps_null, reps_alt))]
cells[, id := .I]
stopifnot(all(cells$reps %% CHUNK == 0))
## one job per chunk; the chunk seed depends only on the cell and chunk
jobs <- cells[, .(chunk = seq_len(reps / CHUNK)), by = id][cells, on = 'id']
jobs[, `:=`(n_reps = CHUNK, seed = 20261006L + 1000L * id + chunk)]
cat(sprintf('cells: %d, chunks: %d of %d reps, analyses per replicate: %d, reps alt/null: %d/%d, cores: %d\n',
            nrow(cells), nrow(jobs), CHUNK, nrow(ANALYSES), reps_alt, reps_null, cores))
t0 <- Sys.time()
out <- mclapply(split(jobs, seq_len(nrow(jobs))), run_chunk_checkpointed,
                mc.cores = cores, mc.preschedule = FALSE)
bad <- which(sapply(out, inherits, 'try-error'))
if (length(bad)) stop('chunks failed: ', paste(bad, collapse = ', '))
reps_all <- rbindlist(out)
saveRDS(reps_all, file.path(out_dir, paste0('replicates-', run_tag, '.rds')))
summ <- reps_all[, .(reject = mean(p < 0.05, na.rm = TRUE), mean_est = mean(est, na.rm = TRUE),
                     emp_sd = sd(est, na.rm = TRUE), mean_se = mean(se, na.rm = TRUE),
                     fails = sum(is.na(p)), fits = .N, min_eig = min(min_eig)),
                 by = .(cell, coding, mean, cov, test, extra, expo, th_a, engine)]
summ[, kappa := mean_se / emp_sd]
summ <- merge(cells[, .(cell = id, construct, design, hyp, t_half, rho, cbm)], summ, by = 'cell')
summ[, mcse := sqrt(reject * (1 - reject) / (fits - fails))]
fwrite(summ, file.path(out_dir, paste0('summary-', run_tag, '.csv')))
cat(sprintf('elapsed: %.1f min; failed fits: %d of %d\n',
            as.numeric(difftime(Sys.time(), t0, units = 'mins')),
            sum(summ$fails), sum(summ$fits)))
summ[, analysis := paste0(ifelse(coding == 'covariate', 'BL+', ''), cov,
                          c(Satterthwaite = '', `Kenward-Roger` = '-KR', CR2 = '-CR2')[test],
                          ifelse(extra == 'sepbm', '-sepbm', ''),
                          ifelse(expo == 'Dbc', paste0('-Dbc', th_a), ''))]
for (h in c('null', 'alt')) {
  cat(sprintf('\nRejection rate (%s), rows: construct x design x t1/2 x mean; cols: covariance\n', h))
  print(dcast(summ[hyp == h], construct + design + t_half + mean ~ analysis,
              value.var = 'reject'), digits = 3)
}
cat('\nSE ratio kappa = mean model SE / empirical SD (alternative cells)\n')
print(dcast(summ[hyp == 'alt'], construct + design + t_half + mean ~ analysis,
            value.var = 'kappa'), digits = 3)
