## Stress test of the provisional (strawman) analysis for the biomarker-
## by-drug interaction in the Hybrid design, as pre-registered in
## docs/46-strawman-stress-test-preregistration.md. Cells, analyses and
## settings follow that document; any departure is a deviation to be
## reported in docs/47.
##
## Data: the 04-power-simulation.R construct (published Gompertz means
## and SDs, Hybrid design, published allocation scaled to N), graded
## coupling of the biomarker to BR, extended with prognostic couplings
## to TV and PB, mean moderation of TV, phase-specific measurement
## error, Weibull carryover decay and monotone dropout (docs/46,
## Section 5). No positive-definiteness repair: a cell whose matrix is
## not positive definite stops the run.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R \
##     [--pilot R] [--cores K] [--chunk C] [--check]
## --check runs the pre-run checks only. --pilot R runs R replicates per
## cell into a separate directory. Chunks are checkpointed; a relaunch
## with the same arguments resumes.
suppressMessages({
  library(data.table); library(MASS); library(mmrm); library(nlme); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
cores <- arg('--cores', 8L)
CHUNK <- arg('--chunk', 25L)
pilot <- arg('--pilot', 0L)
check_only <- '--check' %in% args
## --extra: post hoc analyses only (not pre-registered in docs/46), on the
## same cells and seeds, so the trials are identical to the main run. S is
## refitted as a pairing check against the main run.
extra <- '--extra' %in% args
## --extra2: post hoc F1 with AIC-selected half-life (F1A), and S again.
extra2 <- '--extra2' %in% args
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'
run_tag <- paste0(if (pilot > 0) sprintf('pilot%d', pilot) else 'full',
                  if (extra) '-extra' else if (extra2) '-extra2' else '')
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
  for (p in 2:n) if (!pi$on[p] && pi$tsd[p] > 0) br[p] <- br[p] + br[p - 1] * w[p]
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
  d <- make_trial(cl, setups[[cl$id]]$bs, setups[[cl$id]]$Ns, 1)
  print(fit_all(d))
  quit(save = 'no')
}

## ---------------- run ----------------
chunk_n <- if (pilot > 0) pilot else CHUNK
jobs <- cells[, .(chunk = seq_len(ceiling(reps / chunk_n))), by = id][cells, on = 'id']
jobs[, `:=`(n_reps = pmin(chunk_n, reps - (chunk - 1L) * chunk_n),
            seed = 20261007L + 1000L * id + chunk)]
run_chunk <- function(job) {
  f <- file.path(cell_dir, sprintf('cell-%03d-chunk-%02d.rds', job$id, job$chunk))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  set.seed(job$seed)
  st <- setups[[job$id]]
  out <- rbindlist(lapply(seq_len(job$n_reps), function(r) {
    d <- make_trial(job, st$bs, st$Ns, r)
    fit_all(d)[, rep := (job$chunk - 1L) * chunk_n + r]
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
