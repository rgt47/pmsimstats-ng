## Stress test of the provisional (strawman) analysis for the biomarker-
## by-drug interaction in the Hybrid design, as specified in the working
## protocol docs/46-strawman-stress-test-protocol.md. Cells, analyses and
## settings follow that document; any change to the protocol is to be
## reported in docs/47.
##
## Data: the 04-power-simulation.R construct (published Gompertz means
## and SDs, Hybrid design, allocation per --allocation below), graded
## coupling of the biomarker to BR, extended with prognostic couplings
## to TV and PB, mean moderation of TV, phase-specific measurement
## error, Weibull carryover decay and monotone dropout (docs/46,
## Section 5). No positive-definiteness repair: a cell whose matrix is
## not positive definite stops the run.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R \
##     [--pilot R] [--cores K] [--chunk C] [--check]
##     [--carryover anchored|recursive] [--allocation random|fixed]
##     [--out-dir PATH]
## --check runs the pre-run checks only. --pilot R runs R replicates per
## cell into a separate directory. Chunks are checkpointed; a relaunch
## with the same arguments resumes. --cores defaults to SIM_CORES (set
## by simrun) or 8.
##
## Refined process (2026-10-08; paper 01, tab:spec). --carryover
## anchored (default) decays the BR mean from its last on-drug value,
## mu_last * w(tsd), as R/generateData.R does; recursive is the earlier
## br[p] + br[p - 1] * w(tsd), which compounds the decay at the second
## and later off-drug visits because tsd is cumulative. --allocation
## random (default) assigns each participant to a path uniformly at
## random in every trial; fixed is the earlier near-equal split
## (18/18/17/17 at N = 70). '--carryover recursive --allocation fixed'
## reproduces the stored results and their untagged file names; the
## defaults write files tagged '-anchored', other combinations
## '-anchored-fixed' or '-recursive-random'.
##
## --extra3 (2026-10-09; RERUN-ext3.md): data-chosen analysis half-lives
## for the interaction (docs/57, Section 7.3) on fresh trials: a new seed
## base (57261009 in place of 20261007), the 41 cells plus 17 new ones
## (--ext3-cells docs57 keeps only the 12 of docs/57, Section 7.3), and
## the definitions of analysis/scripts/data-chosen-half-life/57-common.R.
## Writes files tagged '-extra3'. Without --extra3 nothing changes.
suppressMessages({
  library(data.table); library(MASS); library(mmrm); library(nlme); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
sarg <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
cores <- arg('--cores', as.integer(Sys.getenv('SIM_CORES', '8')))
CARRYOVER <- sarg('--carryover', 'anchored')
ALLOCATION <- sarg('--allocation', 'random')
stopifnot(CARRYOVER %in% c('anchored', 'recursive'), ALLOCATION %in% c('random', 'fixed'))
dgp_tag <- if (CARRYOVER == 'recursive' && ALLOCATION == 'fixed') '' else
  paste0('-', CARRYOVER, if (ALLOCATION == 'fixed') '-fixed' else
    if (CARRYOVER == 'recursive') '-random' else '')
CHUNK <- arg('--chunk', 25L)
pilot <- arg('--pilot', 0L)
check_only <- '--check' %in% args
## --extra: post hoc analyses only (not specified in docs/46), on the
## same cells and seeds, so the trials are identical to the main run. S is
## refitted as a pairing check against the main run.
extra <- '--extra' %in% args
## --extra2: post hoc F1 with AIC-selected half-life (F1A), and S again.
extra2 <- '--extra2' %in% args
extra3 <- '--extra3' %in% args
EXT3_CELLS <- sarg('--ext3-cells', 'all')
stopifnot(sum(extra, extra2, extra3) <= 1, EXT3_CELLS %in% c('all', 'docs57'))
SEED_BASE <- if (extra3) 57261009L else 20261007L
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
out_dir <- sarg('--out-dir', 'analysis/data/quick-sim/carryover-closed-form/strawman-stress')
run_tag <- paste0(if (pilot > 0) sprintf('pilot%d', pilot) else 'full',
                  if (extra) '-extra' else if (extra2) '-extra2' else if (extra3) '-extra3' else '',
                  if (extra3 && EXT3_CELLS == 'docs57') '-docs57' else '', dgp_tag)
cell_dir <- file.path(out_dir, paste0('cells-', run_tag))
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('C1', 'CX', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd', 'gomp', 'designs', 'path_info')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))
D <- dgp$designs$Hybrid
PIS <- lapply(D$paths, function(od) dgp$path_info(D$t_wk, D$e, od))
TV_MAX <- max(dgp$gomp(PIS[[1]]$w, 6.50647, 5, 0.35))

## ---------------- cells (docs/46, Section 5) ----------------
base <- list(construct = 'C07', rho = 0.7, N = 70L, t_half = 0.5, decay_k = 1,
             dropout = 'none', c_tv = 0, b_tv = 0, c_pb = 0, c_bl = 0, meas = FALSE)
cell <- function(label, cbm, ...) {
  x <- modifyList(base, list(...))
  as.data.table(c(list(label = label, cbm = cbm), x))
}
cells <- rbindlist(c(
  lapply(c(-0.30, 0, 0.075, 0.15, 0.225, 0.30, 0.375, 0.45),
         function(v) cell('curve-C07', v)),
  lapply(c(-0.20, 0, 0.05, 0.10, 0.15, 0.20, 0.25),
         function(v) cell('curve-CS', v, construct = 'CS', rho = 0.8)),
  unlist(lapply(c(0, 0.30), function(v) list(
    cell('meas-error', v, meas = TRUE),
    cell('N35', v, N = 35L),
    cell('N150', v, N = 150L),
    cell('thalf0', v, t_half = 0),
    cell('thalf1', v, t_half = 1),
    cell('weibull0.5', v, decay_k = 0.5),
    cell('weibull2', v, decay_k = 2),
    cell('dropout-MCAR', v, dropout = 'MCAR'),
    cell('dropout-MAR', v, dropout = 'MAR'),
    cell('prog-TV-const', v, c_tv = 0.3),
    cell('prog-TV-grow', v, b_tv = 0.3),
    cell('prog-PB', v, c_pb = 0.3))), recursive = FALSE),
  ## Added after the main run (2026-10-07), appended so that the ids and
  ## seeds of cells 1-39 are unchanged: the biomarker correlated with
  ## baseline severity.
  lapply(c(0, 0.30), function(v) cell('prog-BL', v, c_bl = 0.3))))
## --extra3 only, appended after cell 41 so that ids 1-41 are unchanged.
## Cells 42-53 (docs/57, Section 7.3): each prognostic pathway of cells
## 34-39 crossed with a true half-life the fixed analysis value (0.5 week)
## gets wrong, null and c_bm = 0.30. Cells 54-58 (omitted with
## --ext3-cells docs57): alternatives at a wrong fixed half-life beyond
## cells 25 and 27, and the null twin of the new t1/2 = 2 label.
if (extra3) {
  ext3 <- list()
  for (th in c(0, 1)) for (v in c(0, 0.30)) ext3 <- c(ext3, list(
    cell(sprintf('prog-TV-const-thalf%g', th), v, c_tv = 0.3, t_half = th),
    cell(sprintf('prog-TV-grow-thalf%g', th), v, b_tv = 0.3, t_half = th),
    cell(sprintf('prog-PB-thalf%g', th), v, c_pb = 0.3, t_half = th)))
  if (EXT3_CELLS == 'all') ext3 <- c(ext3, list(
    cell('thalf2', 0, t_half = 2), cell('thalf2', 0.30, t_half = 2),
    cell('thalf0', 0.45, t_half = 0), cell('thalf1', 0.45, t_half = 1),
    cell('thalf2', 0.45, t_half = 2)))
  cells <- rbindlist(c(list(cells), ext3))
}
cells[, `:=`(id = .I, hyp = fifelse(cbm == 0, 'null', 'alt'))]
cells[, reps := if (pilot > 0) pilot else fifelse(hyp == 'null', 1000L, 500L)]
if (pilot == 0) stopifnot(all(cells$reps %% CHUNK == 0))

## ---------------- data-generating process ----------------
## Carryover decay of the drug effect: exponential (k = 1) or Weibull
## with the same half-life; zero weight before any exposure.
decay <- function(tsd, th, k) {
  if (th <= 0) return(rep(0, length(tsd)))
  lw <- log(2)^(1 / k) / th
  ifelse(tsd > 0, exp(-(lw * tsd)^k), 0)
}

build <- function(pi, cl) {
  n <- length(pi$w)
  br <- dgp$gomp(pi$tod, 10.98604, 5, 0.42)
  w <- decay(pi$tsd, cl$t_half, cl$decay_k)
  if (CARRYOVER == 'recursive') {
    for (p in 2:n) if (!pi$on[p] && pi$tsd[p] > 0) br[p] <- br[p] + br[p - 1] * w[p]
  } else {
    ## tsd is cumulative, so w is applied once to the mean at
    ## discontinuation (R/generateData.R); recursing on br[p - 1] would
    ## count the elapsed time twice.
    last_on <- 0
    for (p in seq_len(n)) {
      if (pi$on[p]) last_on <- br[p]
      else if (pi$tsd[p] > 0) br[p] <- br[p] + last_on * w[p]
    }
  }
  mu <- c(dgp$bm_m, dgp$bl_m, dgp$gomp(pi$w, 6.50647, 5, 0.35),
          dgp$gomp(pi$tpb, 6.50647, 5, 0.35) * pi$e, br)
  sds <- c(dgp$bm_sd, dgp$bl_sd, rep(10, n), 10 * pi$e, rep(8, n))
  g <- abs(outer(pi$w, pi$w, '-'))
  if (cl$construct == 'CS') {
    A <- matrix(cl$rho, n, n); diag(A) <- 1
    C <- matrix(dgp$CX, n, n); diag(C) <- dgp$C1
  } else {
    A <- cl$rho^g
    C <- dgp$C1 * A
  }
  R <- diag(2 + 3 * n)
  R[3:(2 + 3 * n), 3:(2 + 3 * n)] <- rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
  i_tv <- 3:(2 + n); i_pb <- (3 + n):(2 + 2 * n); i_br <- (3 + 2 * n):(2 + 3 * n)
  coup <- cl$cbm * ifelse(pi$on, 1, w)
  R[1, i_br] <- R[i_br, 1] <- coup
  R[1, i_tv] <- R[i_tv, 1] <- cl$c_tv
  R[1, i_pb] <- R[i_pb, 1] <- cl$c_pb
  R[1, 2] <- R[2, 1] <- cl$c_bl
  list(mu = mu, Sigma = outer(sds, sds) * R, R = R, n = n)
}

draw <- function(b, pi, cl, n_draw) {
  X <- mvrnorm(n_draw, b$mu, b$Sigma)
  n <- b$n
  bm <- X[, 1]; BL <- X[, 2]
  tv <- X[, 3:(2 + n), drop = FALSE]
  pb <- X[, (3 + n):(2 + 2 * n), drop = FALSE]
  br <- X[, (3 + 2 * n):(2 + 3 * n), drop = FALSE]
  if (cl$b_tv != 0) {
    z <- (bm - dgp$bm_m) / dgp$bm_sd
    shape <- dgp$gomp(pi$w, 6.50647, 5, 0.35) / TV_MAX
    tv <- tv + cl$b_tv * 10 * outer(z, shape)
  }
  Y <- BL - (tv + pb + br)
  if (cl$meas) {
    s <- ifelse(pi$e == 1, 3, 6)
    Y <- Y + matrix(rnorm(length(Y)), nrow(Y)) * rep(s, each = nrow(Y))
  }
  list(bm = bm, BL = BL, Y = Y)
}

## Monotone dropout: returns, per participant, the first post-baseline
## visit index that is missing (Inf if complete). Dropout can start at
## visit 3. MAR hazard depends on the standardized previous Y, with the
## intercept solved so the expected overall dropout is 20%.
dropout_visit <- function(Y, type) {
  n <- nrow(Y); nv <- ncol(Y)
  out <- rep(Inf, n)
  if (type == 'MCAR') {
    hit <- runif(n) < 0.2
    out[hit] <- sample(3:nv, sum(hit), replace = TRUE)
  } else if (type == 'MAR') {
    zp <- scale(Y[, 2:(nv - 1), drop = FALSE])
    p_drop <- function(a) mean(1 - apply(1 - plogis(a + zp), 1, prod))
    a <- uniroot(function(a) p_drop(a) - 0.2, c(-15, 5))$root
    h <- plogis(a + zp)
    u <- matrix(runif(length(h)), nrow(h))
    for (i in seq_len(n)) {
      j <- which(u[i, ] < h[i, ])
      if (length(j)) out[i] <- j[1] + 2
    }
  }
  out
}

## ---------------- analyses (docs/46, Section 4) ----------------
TH_A <- 0.5
GRID <- c(0, 0.25, 0.5, 1, 2)
dbc <- function(Db, tsd, h) {
  if (h <= 0) return(as.numeric(Db))
  ifelse(Db, 1, ifelse(tsd > 0, 0.5^(tsd / h), 0))
}
na3 <- c(est = NA_real_, se = NA_real_, p = NA_real_)
lctrl <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)

fit_lme <- function(f, d, method = 'REML', hetero = FALSE) tryCatch(
  nlme::lme(f, random = ~ 1 | ptID, correlation = nlme::corCAR1(form = ~ t | ptID),
            weights = if (hetero) nlme::varIdent(form = ~ 1 | phase) else NULL,
            data = d, method = method, control = lctrl),
  error = function(e) NULL)
cr2 <- function(fit, target = 'bmc:X') {
  if (is.null(fit)) return(na3)
  ct <- tryCatch(clubSandwich::coef_test(fit, vcov = 'CR2', test = 'Satterthwaite'),
                 error = function(e) NULL)
  if (is.null(ct)) return(na3)
  k <- which(ct$Coef == target)
  c(est = ct$beta[k], se = ct$SE[k], p = ct$p_Satt[k])
}
mb <- function(fit, target = 'bmc:X') {
  if (is.null(fit)) return(na3)
  tt <- summary(fit)$tTable
  c(est = tt[target, 'Value'], se = tt[target, 'Std.Error'], p = tt[target, 'p-value'])
}
fit_mmrm <- function(f, d, kind) tryCatch(suppressWarnings(switch(kind,
  KR = mmrm(f, data = d, method = 'Kenward-Roger', vcov = 'Kenward-Roger-Linear'),
  ML = mmrm(f, data = d, reml = FALSE),
  Satt = mmrm(f, data = d))), error = function(e) NULL)
mm <- function(fit, target = 'bmc:X') {
  if (is.null(fit)) return(na3)
  cf <- summary(fit)$coefficients
  c(est = cf[target, 'Estimate'], se = cf[target, 'Std. Error'], p = cf[target, 'Pr(>|t|)'])
}

F_S  <- Sx ~ bmc + t + X + bmc:X
F_F1 <- Sx ~ bmc + t + X + bmc:X + De + bmc:De
F_F2 <- Sx ~ bmc + t + X + bmc:X + bmc:t
F_E  <- Sx ~ bmc + t + X + bmc:X + bmc:bl0
F_A  <- Sx ~ bmc + t + X + bmc:X + us(visit | ptID)
F_EA <- Sx ~ bmc + t + X + bmc:X + bmc:bl0 + us(visit | ptID)
F_R  <- Sx ~ bmc + t + X + bmc:X + cs(visit | ptID)

with_x <- function(d, h) { d <- copy(d); d[, X := dbc(Db, tsd, h)]; d }

fit_all <- function(d) {
  d05 <- with_x(d, TH_A)
  res <- list()
  s_fit <- fit_lme(F_S, d05)
  res$S  <- c(cr2(s_fit), h = TH_A)
  c_fit <- fit_lme(F_S, d05, hetero = TRUE)
  res$C1 <- c(cr2(c_fit), h = TH_A)
  res$C2 <- c(mb(c_fit), h = TH_A)
  res$A  <- c(mm(fit_mmrm(F_A, d05, 'KR')), h = TH_A)
  res$ES <- c(cr2(fit_lme(F_E, d05)), h = TH_A)
  res$EA <- c(mm(fit_mmrm(F_EA, d05, 'KR')), h = TH_A)
  res$F1 <- c(cr2(fit_lme(F_F1, d05)), h = TH_A)
  res$F2 <- c(cr2(fit_lme(F_F2, d05)), h = TH_A)
  res$R  <- c(mm(fit_mmrm(F_R, with_x(d, 0), 'Satt')), h = 0)
  ## AIC selection of the analysis half-life: ML for the comparison,
  ## REML refit of the selected half-life for inference.
  ds <- lapply(GRID, function(h) with_x(d, h))
  aic_s <- sapply(ds, function(dd) { f <- fit_lme(F_S, dd, 'ML')
                                     if (is.null(f)) NA_real_ else AIC(f) })
  if (all(is.na(aic_s))) res$DS <- c(na3, h = NA_real_) else {
    k <- which.min(aic_s)
    res$DS <- c(if (GRID[k] == TH_A) cr2(s_fit) else cr2(fit_lme(F_S, ds[[k]])), h = GRID[k])
  }
  aic_a <- sapply(ds, function(dd) { f <- fit_mmrm(F_A, dd, 'ML')
                                     if (is.null(f)) NA_real_ else AIC(f) })
  if (all(is.na(aic_a))) res$DA <- c(na3, h = NA_real_) else {
    k <- which.min(aic_a)
    res$DA <- c(if (GRID[k] == TH_A) res$A[1:3] else mm(fit_mmrm(F_A, ds[[k]], 'KR')),
                h = GRID[k])
  }
  rbindlist(lapply(names(res), function(nm) as.data.table(as.list(res[[nm]]))[, analysis := nm]))
}

## Post hoc (docs/48, Section 5): F1's mean model with phase-specific
## residual variances, CR2 (G1) and model-based (G2) inference; and S again.
fit_extra <- function(d) {
  d05 <- with_x(d, TH_A)
  g_fit <- fit_lme(F_F1, d05, hetero = TRUE)
  res <- list(S = c(cr2(fit_lme(F_S, d05)), h = TH_A),
              G1 = c(cr2(g_fit), h = TH_A), G2 = c(mb(g_fit), h = TH_A))
  rbindlist(lapply(names(res), function(nm) as.data.table(as.list(res[[nm]]))[, analysis := nm]))
}
## Post hoc: F1 with the analysis half-life chosen by AIC over GRID (ML for
## the comparison, REML refit of the selected half-life, CR2), as DS does
## for S; and S again as the pairing check.
fit_extra2 <- function(d) {
  d05 <- with_x(d, TH_A)
  ds <- lapply(GRID, function(h) with_x(d, h))
  aic <- sapply(ds, function(dd) { f <- fit_lme(F_F1, dd, 'ML')
                                   if (is.null(f)) NA_real_ else AIC(f) })
  f1a <- if (all(is.na(aic))) c(na3, h = NA_real_) else {
    k <- which.min(aic)
    c(cr2(fit_lme(F_F1, ds[[k]])), h = GRID[k])
  }
  res <- list(S = c(cr2(fit_lme(F_S, d05)), h = TH_A), F1A = f1a)
  rbindlist(lapply(names(res), function(nm) as.data.table(as.list(res[[nm]]))[, analysis := nm]))
}
if (extra) fit_all <- fit_extra
if (extra2) fit_all <- fit_extra2
## Data-chosen analysis half-lives (docs/57; RERUN-ext3.md). The
## grid fits, F1, F1A, F1A0, AIC-weight averaging (MA), the Bonferroni
## maximum and the null-model bootstrap draw are the definitions of
## 57-common.R, evaluated here; F1G and MA0 (post hoc, docs/57, Section
## 4.2 and 6.2) follow 57-posthoc-rules.R. Added: F1 at the
## pharmacological 0.2 week (F1h02), and the maximum |t| over the grid
## calibrated by a multivariate normal (MaxTMVN) or t with the smallest
## Satterthwaite df rounded down (MaxTMVT), with the correlation of the T(h) from the
## participants' model-based influence contributions (docs/58, Section
## 9.2). One null-model resample per trial (warp-speed method, Giacomini,
## Politis and White 2013) gives tstar; the B arms are formed after the
## run from the pooled tstar of each cell.
if (extra3) {
  keep57 <- c('G5', 'FM1', 'FM0', 'TARGET', 'na5', 'cr2t', 'll_or_na', 'grid_fits',
              'pick', 'wmax', 'grid_methods', 'null_fit', 'sim_null')
  for (e in parse('analysis/scripts/data-chosen-half-life/57-common.R')) {
    if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% keep57) eval(e, globalenv())
  }
  stopifnot(all(keep57 %in% ls(globalenv())), identical(G5, GRID))
  suppressMessages(library(mvtnorm))
  TH_P <- 0.2
  BOOT_OF <- c(F1 = 'F1B', F1h02 = 'F1h02B', F1A = 'F1AB', F1A0 = 'F1A0B', F1G = 'F1GB',
               MA = 'MAB', MA0 = 'MA0B', MaxTBonf = 'MaxTB')
  POSTHOC <- c('F1G', 'MA0', 'F1GB', 'MA0B')

  ## Participant-level influence of the interaction estimate in an lme
  ## fit with random intercept and CAR(1) (infl_X of 58-common.R, for
  ## the target bmc:X).
  infl_int <- function(fit, d) {
    su2 <- as.numeric(nlme::VarCorr(fit)['(Intercept)', 'Variance'])
    phi <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
    s2 <- fit$sigma^2
    Xm <- model.matrix(formula(fit)[-2], d)
    r <- d$Sx - drop(Xm %*% nlme::fixef(fit))
    parts <- lapply(split(seq_len(nrow(d)), d$ptID), function(ii) {
      tt <- d$t[ii]; Vi <- solve(su2 + s2 * phi^abs(outer(tt, tt, '-')))
      list(B = crossprod(Xm[ii, , drop = FALSE], Vi %*% Xm[ii, , drop = FALSE]),
           u = drop(crossprod(Xm[ii, , drop = FALSE], Vi %*% r[ii])))
    })
    Bi <- solve(Reduce('+', lapply(parts, `[[`, 'B')))
    k <- which(colnames(Xm) == TARGET)
    vapply(parts, function(p) drop(Bi %*% p$u)[k], numeric(1))
  }
  ## grid_fits of 57-common.R, also returning the influence vectors.
  grid_fits3 <- function(d, grid = G5, infl = TRUE) {
    U <- list()
    g <- rbindlist(lapply(grid, function(h) {
      dd <- with_x(d, h)
      fr <- fit_lme(FM1, dd)
      if (infl && !is.null(fr)) U[[as.character(h)]] <<-
        tryCatch(infl_int(fr, dd), error = function(e) NULL)
      c(list(h = h, ll1 = ll_or_na(fit_lme(FM1, dd, 'ML')),
             ll0 = ll_or_na(fit_lme(FM0, dd, 'ML'))), as.list(cr2t(fr)))
    }))
    list(g = g, U = U)
  }
  ## F1G and MA0 exactly as in 57-posthoc-rules.R.
  posthoc_rules <- function(g) {
    dev05 <- 2 * (max(g$ll0, na.rm = TRUE) - g$ll0[g$h == 0.5])
    sw <- isTRUE(dev05 > qchisq(0.95, 1))
    f1g <- pick(g, if (sw) which.max(g$ll0) else which(g$h == 0.5), 'F1G')[, switched := sw]
    ok <- is.finite(g$ll0) & is.finite(g$est)
    if (!any(ok)) return(rbind(f1g, pick(g, NA, 'MA0'), fill = TRUE))
    w <- exp(g$ll0[ok] - max(g$ll0[ok])); w <- w / sum(w)
    eb <- sum(w * g$est[ok]); sb <- sum(w * sqrt(g$se[ok]^2 + (g$est[ok] - eb)^2))
    dfb <- sum(w * g$df[ok])
    rbind(f1g, data.table(method = 'MA0', h = sum(w * g$h[ok]), est = eb, se = sb,
                          t = eb / sb, df = dfb, p = 2 * pt(-abs(eb / sb), dfb)), fill = TRUE)
  }
  ## Maximum |T(h)| calibrated by the multivariate normal and t.
  maxt_cal <- function(g, U) {
    ok <- is.finite(g$t) & as.character(g$h) %in% names(U)
    if (sum(ok) < 2) return(rbind(pick(g, NA, 'MaxTMVN'), pick(g, NA, 'MaxTMVT')))
    gk <- g[ok]
    km <- which.max(abs(gk$t)); tm <- abs(gk$t[km]); k <- nrow(gk)
    R <- cov2cor(crossprod(do.call(cbind, U[as.character(gk$h)])))
    pn <- tryCatch(1 - as.numeric(pmvnorm(lower = rep(-tm, k), upper = rep(tm, k), corr = R,
                   algorithm = GenzBretz(abseps = 1e-4))), error = function(e) NA_real_)
    ## pmvt needs an integer df: the smallest Satterthwaite df, rounded down.
    dft <- max(1L, as.integer(floor(min(gk$df))))
    pv <- tryCatch(1 - as.numeric(pmvt(lower = rep(-tm, k), upper = rep(tm, k), corr = R,
                   df = dft, algorithm = GenzBretz(abseps = 1e-4))),
                   error = function(e) NA_real_)
    a <- pick(gk, km, 'MaxTMVN')[, `:=`(t = tm, p = pn)]
    b <- pick(gk, km, 'MaxTMVT')[, `:=`(t = tm, df = dft, p = pv)]
    rbind(a, b)
  }
  ## Every nominal procedure from one trial's grid.
  methods3 <- function(d, G, cal = TRUE) {
    g <- G$g
    m <- rbind(grid_methods(g), posthoc_rules(g), fill = TRUE)
    m[method == 'MaxT', method := 'MaxTBonf']
    t02 <- cr2t(fit_lme(FM1, with_x(d, TH_P)))
    rbind(m, if (cal) maxt_cal(g, G$U),
          data.table(method = 'F1h02', h = TH_P, est = t02[['est']], se = t02[['se']],
                     t = t02[['t']], df = t02[['df']], p = t02[['p']]), fill = TRUE)
  }
  fit_extra3 <- function(d, aux_seed = 1L) {
    t0 <- proc.time()[['elapsed']]
    old <- .Random.seed
    on.exit(assign('.Random.seed', old, globalenv()))
    set.seed(aux_seed)
    m <- methods3(d, grid_fits3(d))
    m[, tstar := NA_real_]
    h0 <- m[method == 'F1A0', h]
    nf <- if (is.finite(h0)) null_fit(d, h0) else NULL
    if (!is.null(nf)) {
      ds <- tryCatch(sim_null(nf, d, h0), error = function(e) NULL)
      if (!is.null(ds)) {
        mb <- methods3(ds, grid_fits3(ds, infl = FALSE), cal = FALSE)
        m[method %in% names(BOOT_OF), tstar := abs(mb$t[match(method, mb$method)])]
      }
    }
    setnames(m, 'method', 'analysis')
    m[, `:=`(posthoc = analysis %in% POSTHOC, secs = proc.time()[['elapsed']] - t0)]
    m[]
  }
  ## Warp-speed arms: p* = (1 + #{tstar >= |t|}) / (n + 1) over the n
  ## resampled statistics of the cell; crit_ws is their 95% quantile.
  warp_speed <- function(x) {
    w <- x[analysis %in% names(BOOT_OF)]
    w[, c('p', 'crit_ws') := {
      s <- sort(tstar[is.finite(tstar)]); n <- length(s); a <- abs(t)
      list(fifelse(is.finite(a) & n > 0, (1 + n - findInterval(a, s, left.open = TRUE)) / (n + 1),
                   NA_real_), if (n > 0) unname(quantile(s, 0.95)) else NA_real_)
    }, by = .(cell, analysis)]
    w[, `:=`(analysis = BOOT_OF[analysis], tstar = NA_real_)]
    rbind(x, w, fill = TRUE)
  }
  fit_all <- fit_extra3
}

## ---------------- one trial ----------------
make_trial <- function(cl, bs, Ns, r) {
  long <- rbindlist(lapply(seq_along(PIS), function(k) {
    pi <- PIS[[k]]
    dr <- draw(bs[[k]], pi, cl, Ns[k])
    drop_at <- dropout_visit(dr$Y, cl$dropout)
    id <- (k - 1) * 1000 + seq_len(Ns[k])
    nv <- length(pi$w)
    post <- data.table(ptID = rep(id, times = nv), bm = rep(dr$bm, times = nv),
                       Sx = as.vector(dr$Y), vis = rep(seq_len(nv), each = Ns[k]),
                       t = rep(pi$w, each = Ns[k]), Db = rep(pi$on, each = Ns[k]),
                       tsd = rep(pi$tsd, each = Ns[k]), De = rep(pi$e, each = Ns[k]),
                       drop_at = rep(drop_at, times = nv))
    post <- post[vis < drop_at]
    rbind(data.table(ptID = id, bm = dr$bm, Sx = dr$BL, vis = 0L, t = 0, Db = FALSE,
                     tsd = 0, De = 0, drop_at = drop_at), post)
  }))
  long[, `:=`(ptID = factor(ptID),
              visit = factor(as.character(t), levels = as.character(c(0, PIS[[1]]$w))),
              bmc = bm - mean(bm[t == 0]), bl0 = as.numeric(t == 0),
              phase = factor(fifelse(t == 0, 'baseline', fifelse(De == 1, 'OL', 'blinded'))))]
  long
}

cell_setup <- function(cl) {
  bs <- lapply(PIS, build, cl = cl)
  min_eig <- min(sapply(bs, function(b) min(eigen(b$R, TRUE, TRUE)$values)))
  nP <- length(PIS)
  Ns <- cl$N %/% nP + c(rep(1L, cl$N %% nP), rep(0L, nP - cl$N %% nP))
  list(bs = bs, Ns = Ns, min_eig = min_eig)
}

## Per-trial path counts. Random allocation is redrawn when a path gets
## fewer than 2 participants, because the MAR dropout model standardizes
## Y within a path; the redraw probability is about 0.002 at N = 35 and
## negligible at N = 70 or 150.
alloc <- function(cl, st) {
  if (ALLOCATION == 'fixed') return(st$Ns)
  nP <- length(PIS)
  repeat {
    Ns <- tabulate(sample.int(nP, cl$N, replace = TRUE), nP)
    if (all(Ns >= 2L)) return(Ns)
  }
}

## ---------------- pre-run checks (docs/46, Section 9) ----------------
setups <- lapply(seq_len(nrow(cells)), function(i) cell_setup(cells[i]))
cells[, min_eig := sapply(setups, `[[`, 'min_eig')]
print(cells[, .(id, label, cbm, hyp, reps, min_eig = signif(min_eig, 3))], nrows = 100)
if (any(cells$min_eig <= 0)) stop('cells without a positive-definite matrix: ',
                                  paste(cells[min_eig <= 0, id], collapse = ', '))
if (check_only) {
  set.seed(1)
  cat('\nprognostic couplings in a large draw (target 0.3):\n')
  cl <- cells[label == 'prog-TV-const' & cbm == 0][1]
  b <- build(PIS[[1]], cl); X <- mvrnorm(2e5, b$mu, b$Sigma)
  cat(sprintf('  cor(bm, TV at week 20) = %.3f\n', cor(X[, 1], X[, 2 + length(PIS[[1]]$w)])))
  cl <- cells[label == 'prog-PB' & cbm == 0][1]
  b <- build(PIS[[1]], cl); X <- mvrnorm(2e5, b$mu, b$Sigma)
  n <- length(PIS[[1]]$w)
  cat(sprintf('  cor(bm, PB at week 20) = %.3f\n', cor(X[, 1], X[, 2 + 2 * n])))
  cat('dropout rates (target 0.20):\n')
  for (ty in c('MCAR', 'MAR')) {
    cl <- cells[label == paste0('dropout-', ty) & cbm == 0][1]
    dr <- draw(build(PIS[[1]], cl), PIS[[1]], cl, 2e4)
    cat(sprintf('  %s: %.3f\n', ty, mean(is.finite(dropout_visit(dr$Y, ty)))))
  }
  cat('AIC rule against the fixed half-life on one trial:\n')
  cl <- cells[label == 'curve-C07' & cbm == 0.30][1]
  d <- make_trial(cl, setups[[cl$id]]$bs, alloc(cl, setups[[cl$id]]), 1)
  print(fit_all(d))
  quit(save = 'no')
}

## ---------------- run ----------------
chunk_n <- if (pilot > 0) pilot else CHUNK
jobs <- cells[, .(chunk = seq_len(ceiling(reps / chunk_n))), by = id][cells, on = 'id']
jobs[, `:=`(n_reps = pmin(chunk_n, reps - (chunk - 1L) * chunk_n),
            seed = SEED_BASE + 1000L * id + chunk)]
run_chunk <- function(job) {
  f <- file.path(cell_dir, sprintf('cell-%03d-chunk-%02d.rds', job$id, job$chunk))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  set.seed(job$seed)
  st <- setups[[job$id]]
  out <- rbindlist(lapply(seq_len(job$n_reps), function(r) {
    d <- make_trial(job, st$bs, alloc(job, st), r)
    rr <- (job$chunk - 1L) * chunk_n + r
    ## --extra3: the bootstrap draw and the multivariate-normal integrals
    ## use their own stream (seed 571000000 + 10000 * cell + replicate),
    ## and the trial stream is restored afterwards.
    (if (extra3) fit_extra3(d, 571000000L + 10000L * job$id + rr) else fit_all(d))[, rep := rr]
  }))
  out[, cell := job$id]
  saveRDS(out, f)
  cat(sprintf('%s cell %2d chunk %2d %-14s cbm=%6.3f  %.1f min\n',
              format(Sys.time(), '%m-%d %H:%M:%S'), job$id, job$chunk, job$label, job$cbm,
              as.numeric(difftime(Sys.time(), t0, units = 'mins'))),
      file = file.path(cell_dir, 'progress.log'), append = TRUE)
  out
}
cat(sprintf('cells: %d, chunks: %d, replicates: %d, cores: %d\n',
            nrow(cells), nrow(jobs), sum(cells$reps), cores))
t0 <- Sys.time()
out <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j) tryCatch(run_chunk(j),
  error = function(e) paste('chunk', j$id, j$chunk, conditionMessage(e))),
  mc.cores = cores, mc.preschedule = FALSE)
bad <- which(!vapply(out, is.data.frame, logical(1)))
if (length(bad)) stop('chunks failed: ', paste(unlist(out[bad]), collapse = ' | '))
reps_all <- merge(rbindlist(out), cells, by.x = 'cell', by.y = 'id')
if (extra3) reps_all <- warp_speed(reps_all)
saveRDS(reps_all, file.path(out_dir, sprintf('replicates-%s.rds', run_tag)))
summ <- reps_all[, .(reject = mean(p < 0.05, na.rm = TRUE), fails = sum(is.na(p)), fits = .N,
                     mean_est = mean(est, na.rm = TRUE), emp_sd = sd(est, na.rm = TRUE),
                     mean_se = mean(se, na.rm = TRUE)),
                 by = .(cell, label, cbm, hyp, analysis)]
summ[, `:=`(kappa = mean_se / emp_sd, mcse = sqrt(reject * (1 - reject) / (fits - fails)))]
fwrite(summ, file.path(out_dir, sprintf('summary-%s.csv', run_tag)))
cat(sprintf('elapsed %.1f min; failed fits %d of %d\n',
            as.numeric(difftime(Sys.time(), t0, units = 'mins')), sum(summ$fails), sum(summ$fits)))
print(dcast(summ, label + cbm ~ analysis, value.var = 'reject'), digits = 3, nrows = 100)
