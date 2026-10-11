## Check of the prognostic-baseline (prog-BL) inflation of the warp-speed
## bootstrap arms in the stress-test extension (docs/57, Section 9;
## RERUN-ext3.md). Two parts.
##
## 1. Stored replicates (no fitting): for every null cell of
##    replicates-full-extra3-anchored.rds, the warp-speed size of the B
##    arms with a Monte Carlo standard error that resamples trials, so
##    that it carries both the binomial error of the observed statistics
##    and the error of the pooled resampled quantile.
## 2. Fresh trials (fitting): cells 2 (curve-C07, null) and 40 (prog-BL,
##    null) of the driver, refined process, new seeds. Per trial, the
##    observed F1 and F1A0 statistics and one resample from each of three
##    null models, rerunning F1A0's selection on the resample:
##      P0  the stored parametric resample (FM0 at the F1A0 choice, REML);
##      P1  the same with a baseline-specific biomarker slope, bmc:bl0;
##      W   a wild cluster bootstrap: FM0's fitted mean plus each
##          participant's residual vector times a Rademacher sign.
##    The warp-speed size of each variant is reported with a trial-
##    resampling MCSE.
##
## Run from the repository root:
##   VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript \
##     analysis/scripts/data-chosen-half-life/57-progbl-check.R \
##     [--reps R] [--cores K] [--out PATH]
## --reps 0 runs part 1 only. Writes PATH.rds (per-trial rows) and prints
## the tables.
suppressMessages(library(data.table))
source('analysis/scripts/data-chosen-half-life/57-common.R')
args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
REPS <- as.integer(argv('--reps', '400'))
CORES <- as.integer(argv('--cores', '6'))
OUT <- argv('--out', 'analysis/data/quick-sim/carryover-closed-form/strawman-stress/progbl-check')
STORED <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress/replicates-full-extra3-anchored.rds'
NB <- 2000L

## Warp-speed size: reject when |t| exceeds the cell's pooled resampled
## statistics in the p* sense, p* = (1 + #{T* >= |t|}) / (n + 1) < 0.05.
ws_size <- function(a, s) {
  s <- sort(s[is.finite(s)]); a <- a[is.finite(a)]; n <- length(s)
  mean((1 + n - findInterval(a, s, left.open = TRUE)) / (n + 1) < 0.05)
}
## MCSE by resampling trials (pairs of observed and resampled statistic).
ws_mcse <- function(a, s, nb = NB) {
  k <- length(a)
  sd(replicate(nb, { i <- sample.int(k, k, replace = TRUE); ws_size(a[i], s[i]) }))
}

## ---------------- part 1: stored replicates ----------------
set.seed(57101)
x <- readRDS(STORED)
BOOT_OF <- c(F1 = 'F1B', F1h02 = 'F1h02B', F1A = 'F1AB', F1A0 = 'F1A0B', F1G = 'F1GB',
             MA = 'MAB', MA0 = 'MA0B', MaxTBonf = 'MaxTB')
st <- x[hyp == 'null' & analysis %in% names(BOOT_OF),
        .(size = ws_size(abs(t), tstar), mcse = ws_mcse(abs(t), tstar),
          q95_obs = quantile(abs(t), 0.95, na.rm = TRUE),
          q95_star = quantile(tstar, 0.95, na.rm = TRUE)),
        by = .(cell, label, analysis)]
cat('Stored replicates: warp-speed size by null cell, with trial-resampling MCSE\n')
print(dcast(st, cell + label ~ analysis, value.var = 'size'), digits = 3, nrows = 30)
cat('\nMCSE of the warp-speed size (trial resampling), by analysis:\n')
print(st[, .(mean_mcse = mean(mcse), binom = sqrt(0.05 * 0.95 / 1000)), by = analysis], digits = 3)
cat('\nSpread of the warp-speed size over the 22 null cells:\n')
print(st[, .(mean = mean(size), sd = sd(size), min = min(size), max = max(size),
             z_progBL = (size[cell == 40] - mean(size)) / mean(mcse)), by = analysis], digits = 3)
cat('\nq95 of the resampled statistic, prog-BL against the other 21 null cells:\n')
print(st[, .(q95_star_40 = q95_star[cell == 40], mean_other = mean(q95_star[cell != 40]),
             sd_other = sd(q95_star[cell != 40]),
             rank40 = rank(q95_star)[cell == 40]), by = analysis], digits = 3)
if (REPS == 0) quit(save = 'no')

## ---------------- part 2: fresh trials ----------------
load_driver('anchored', 'random')
library(parallel)
FM0b <- Sx ~ bmc + t + X + De + bmc:De + bmc:bl0
setups <- list()
for (cid in c(2L, 40L)) setups[[cid]] <- cell_setup(cells[id == cid])
stopifnot(cells[id == 2L, label] == 'curve-C07', cells[id == 40L, label] == 'prog-BL',
          cells[id == 2L, cbm] == 0, cells[id == 40L, cbm] == 0)

sim_from <- function(fit, f0, d, h0) {
  dd <- with_x(d, h0)
  mu <- drop(model.matrix(f0[-2], dd) %*% nlme::fixef(fit))
  tau <- as.numeric(nlme::VarCorr(fit)['(Intercept)', 'StdDev'])
  phi <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  sig <- fit$sigma
  y <- mu
  for (idx in split(seq_len(nrow(dd)), dd$ptID)) {
    tt <- dd$t[idx]
    S <- tau^2 + sig^2 * phi^abs(outer(tt, tt, '-'))
    y[idx] <- mu[idx] + drop(crossprod(chol(S), rnorm(length(idx))))
  }
  ds <- copy(d); ds[, Sx := y]; ds
}
sim_wild <- function(fit, d, h0) {
  dd <- with_x(d, h0)
  mu <- drop(model.matrix(FM0[-2], dd) %*% nlme::fixef(fit))
  r <- dd$Sx - mu
  v <- sample(c(-1, 1), nlevels(dd$ptID), replace = TRUE)[as.integer(dd$ptID)]
  ds <- copy(d); ds[, Sx := mu + v * r]; ds
}
## |t| of F1 (0.5 week) and of F1A0 (ML over the grid without the
## interaction, REML refit with it), and F1A0's choice.
f1_f1a0 <- function(d) {
  ll0 <- vapply(G5, function(h) ll_or_na(fit_lme(FM0, with_x(d, h), 'ML')), numeric(1))
  h0 <- if (all(is.na(ll0))) NA_real_ else G5[which.max(ll0)]
  t1 <- cr2t(fit_lme(FM1, with_x(d, 0.5)))
  t0 <- if (is.na(h0)) na5 else if (h0 == 0.5) t1 else cr2t(fit_lme(FM1, with_x(d, h0)))
  c(F1 = abs(t1[['t']]), F1A0 = abs(t0[['t']]), h0 = h0, p_F1 = t1[['p']], p_F1A0 = t0[['p']])
}
one_trial <- function(cid, r) {
  t_start <- proc.time()[['elapsed']]
  cl <- cells[id == cid]; s <- setups[[cid]]
  set.seed(57990000L + 10000L * cid + r)
  d <- make_trial(cl, s$bs, alloc(cl, s), r)
  obs <- f1_f1a0(d)
  out <- data.table(cell = cid, label = cl$label, rep = r, h0 = obs[['h0']],
                    F1 = obs[['F1']], F1A0 = obs[['F1A0']],
                    p_F1 = obs[['p_F1']], p_F1A0 = obs[['p_F1A0']])
  h0 <- obs[['h0']]
  for (v in c('P0', 'P1', 'W')) {
    st <- c(F1 = NA_real_, F1A0 = NA_real_)
    if (is.finite(h0)) {
      ds <- tryCatch(switch(v,
        P0 = { f <- null_fit(d, h0); if (is.null(f)) NULL else sim_from(f, FM0, d, h0) },
        P1 = { f <- fit_lme(FM0b, with_x(d, h0)); if (is.null(f)) NULL else sim_from(f, FM0b, d, h0) },
        W  = { f <- null_fit(d, h0); if (is.null(f)) NULL else sim_wild(f, d, h0) }),
        error = function(e) NULL)
      if (!is.null(ds)) st <- f1_f1a0(ds)[c('F1', 'F1A0')]
    }
    out[, paste0(c('F1_', 'F1A0_'), v) := as.list(st)]
  }
  out[, secs := proc.time()[['elapsed']] - t_start]
  out
}
## Trials alternate between the two cells and each is checkpointed, so
## an interrupted run keeps a balanced partial result and resumes.
jobs <- CJ(rep = seq_len(REPS), cell = c(2L, 40L))
ck_dir <- paste0(OUT, '-trials')
dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), function(i) {
  f <- file.path(ck_dir, sprintf('cell-%03d-rep-%04d.rds', jobs$cell[i], jobs$rep[i]))
  if (!file.exists(f)) {
    r <- tryCatch(one_trial(jobs$cell[i], jobs$rep[i]), error = function(e) NULL)
    if (!is.null(r)) saveRDS(r, f)
  }
  NULL
}, mc.cores = CORES, mc.preschedule = FALSE))
res <- rbindlist(lapply(list.files(ck_dir, full.names = TRUE), readRDS), fill = TRUE)
res <- res[rep <= REPS]
saveRDS(res, paste0(OUT, '.rds'))
cat(sprintf('\nFresh trials: %d of %d done in %.1f min; median %.1f s per trial\n',
            nrow(res), nrow(jobs), as.numeric(difftime(Sys.time(), t0, units = 'mins')),
            median(res$secs)))

set.seed(57102)
long <- rbindlist(lapply(c('P0', 'P1', 'W'), function(v) rbindlist(lapply(c('F1', 'F1A0'),
  function(s) res[, .(cell, label, variant = v, stat = s, a = get(s),
                      tstar = get(paste0(s, '_', v)))]))))
tab <- long[, .(trials = sum(is.finite(a)), resamples = sum(is.finite(tstar)),
                nominal_q95 = quantile(a, 0.95, na.rm = TRUE),
                q95_star = quantile(tstar, 0.95, na.rm = TRUE),
                m2_ratio = mean(a^2, na.rm = TRUE) / mean(tstar^2, na.rm = TRUE),
                size = ws_size(a, tstar), mcse = ws_mcse(a, tstar)),
            by = .(cell, label, stat, variant)]
cat('\nFresh trials: warp-speed size by cell, statistic and null model\n')
print(tab, digits = 3)
## Fresh resampled quantiles applied to the stored observed statistics of
## the same cell (independent trials, same estimand).
x40 <- x[hyp == 'null' & cell %in% c(2L, 40L) & analysis %in% c('F1', 'F1A0')]
cross <- rbindlist(lapply(c('P0', 'P1', 'W'), function(v) x40[, {
  s <- long[cell == .BY$cell & stat == .BY$analysis & variant == v, tstar]
  list(variant = v, size_stored_obs = ws_size(abs(t), s))
}, by = .(cell, analysis)]))
cat('\nStored observed |t| against the fresh pooled resampled statistics:\n')
print(cross, digits = 3)
cat('\nNominal (CR2/Satterthwaite) size on the fresh trials, binomial MCSE:\n')
print(res[, .(F1 = mean(p_F1 < 0.05, na.rm = TRUE), F1A0 = mean(p_F1A0 < 0.05, na.rm = TRUE),
              mcse = sqrt(0.05 * 0.95 / .N)), by = .(cell, label)], digits = 3)
cat('\nF1A0 choice on the fresh trials:\n')
print(dcast(res[, .N, by = .(cell, h0)], cell ~ h0, value.var = 'N'))
