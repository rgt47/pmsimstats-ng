## Power under the three adjustments, at Hendrickson's published
## parameters and with the published analysis model.
##
## Arms (one change at a time from the published construct):
##   A  published   : CS within, constant cross, step coupling
##   B  adj 1       : AR(1) within, decaying cross, step coupling
##   C  adj 1 + 2   : separable AR(1), graded coupling
##   D  adj 3       : CS within, constant cross, mean moderation
##   E  adj 1+2+3   : separable AR(1), mean moderation
##   F  adj 1 + 2, current covar form: AR(1) within, decaying cross,
##                   graded coupling; with C it isolates separability
## Everything else is held at the published values: rho = 0.8,
## c1 = 0.2, cx = 0.1, Gompertz means, SDs, recursive carryover on the
## BR mean, 8-visit designs. No positive-definiteness repair: a cell
## whose matrix is not PD is skipped and recorded.
##
## Analysis (--analysis):
##   published  the published lme_analysis() model, lmerTest::lmer
##              Sx ~ bm + Db + t + bm*Db + (1|ptID), Satterthwaite p
##   car1       the same fixed effects and random intercept, plus
##              continuous-time AR(1) residuals, nlme::lme with
##              corCAR1(form = ~ t | ptID)
## Both use the baseline as a t = 0 row and the binary Db, so the only
## difference is the residual serial correlation. Seeds are per cell,
## so both analyses see identical data.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R \
##     [--reps R] [--analysis published|car1] [--arms ABCDEF]
## --arms restricts the run to the listed arms; cell ids (and so seeds)
## do not depend on it, and a restricted run writes its tables with the
## arm letters appended, for 05-summarize-power.R to merge.
## Whitepaper: docs/36-hendrickson-pd-and-carryover.md, Section 9.
suppressMessages({
  library(data.table); library(MASS); library(lmerTest); library(parallel)
  library(nlme)
})
args <- commandArgs(trailingOnly = TRUE)
reps <- if ('--reps' %in% args) as.integer(args[which(args == '--reps') + 1]) else 500L
analysis <- if ('--analysis' %in% args) args[which(args == '--analysis') + 1] else 'published'
stopifnot(analysis %in% c('published', 'car1'))
tag <- if (analysis == 'published') '' else paste0('-', analysis)
FIT_LIMIT <- if ('--fit-limit' %in% args) as.integer(args[which(args == '--fit-limit') + 1]) else 30L
ALL_ARMS <- c('A', 'B', 'C', 'D', 'E', 'F')
arms <- if ('--arms' %in% args) strsplit(args[which(args == '--arms') + 1], '')[[1]] else ALL_ARMS[1:5]
stopifnot(all(arms %in% ALL_ARMS))
## The original run covered A-E; any other selection gets its own files.
arm_tag <- if (identical(sort(arms), ALL_ARMS[1:5])) '' else paste0('-', paste(sort(arms), collapse = ''))
out_dir <- 'analysis/data/quick-sim/hendrickson-problems/power-sim'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------- published parameters ----------------
RHO <- 0.8; C1 <- 0.2; CX <- 0.1; N_TOTAL <- 70L
bm_m <- 124.32759; bm_sd <- 15.36159; bl_m <- 83.06897; bl_sd <- 18.48267
gomp <- function(t, maxr, disp, rate) {
  y <- maxr * exp(-disp * exp(-rate * t))
  off <- maxr * exp(-disp)
  (y - off) * maxr / (maxr - off)
}

designs <- list(
  CO = list(t_wk = rep(2.5, 8), e = rep(0.5, 8),
            paths = list(c(1,1,1,1,0,0,0,0), c(0,0,0,0,1,1,1,1))),
  Hybrid = list(t_wk = diff(c(0, 4, 8, 9, 10, 11, 12, 16, 20)),
                e = c(1, 1, .5, .5, .5, .5, .5, .5),
                paths = list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1),
                             c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1))),
  `OL+BDC` = list(t_wk = diff(c(0, 4, 8, 12, 16, 17, 18, 19, 20)),
                  e = c(1, 1, 1, 1, .5, .5, .5, .5),
                  paths = list(c(1,1,1,1,1,1,0,0), c(1,1,1,1,1,0,0,0))))

## Published buildtrialdesign() semantics: tod and tsd accumulate within
## a run and reset when the state changes; tsd is zero before first
## exposure; tpb accumulates while expectancy is positive.
path_info <- function(t_wk, e, od) {
  n <- length(od)
  tod <- t_wk * od; tsd <- t_wk * (od != 1); tpb <- t_wk * (e > 0)
  for (i in 2:n) {
    if (tod[i] > 0) tod[i] <- tod[i] + tod[i - 1] else tsd[i] <- tsd[i] + tsd[i - 1]
    if (tpb[i] > 0) tpb[i] <- tpb[i] + tpb[i - 1]
  }
  tsd <- (cumsum(od) > 0) * tsd
  list(w = cumsum(t_wk), e = e, on = od == 1, tod = tod, tsd = tsd, tpb = tpb)
}

## ---------------- DGP ----------------
build <- function(pi, arm, cbm, t_half) {
  n <- length(pi$w)
  br <- gomp(pi$tod, 10.98604, 5, 0.42)
  for (p in 2:n) if (!pi$on[p] && pi$tsd[p] > 0)
    br[p] <- br[p] + br[p - 1] * 0.5^(pi$tsd[p] / t_half)
  mu <- c(bm_m, bl_m, gomp(pi$w, 6.50647, 5, 0.35),
          gomp(pi$tpb, 6.50647, 5, 0.35) * pi$e, br)
  sds <- c(bm_sd, bl_sd, rep(10, n), 10 * pi$e, rep(8, n))

  g <- abs(outer(pi$w, pi$w, '-'))
  within <- if (arm %in% c('A', 'D')) 'CS' else 'AR1'
  A <- if (within == 'CS') { m <- matrix(RHO, n, n); diag(m) <- 1; m } else RHO^g
  C <- if (arm %in% c('A', 'D')) { m <- matrix(CX, n, n); diag(m) <- C1; m }
       else if (arm %in% c('B', 'F')) { m <- CX * RHO^g; diag(m) <- C1; m }
       else C1 * A
  M <- rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
  R <- diag(2 + 3 * n)
  R[3:(2 + 3 * n), 3:(2 + 3 * n)] <- M

  coup <- switch(arm,
    A = , B = as.numeric(br != 0),
    C = , F = {
      phi <- if (t_half > 0) 0.5^(pi$tsd / t_half) else rep(0, n)
      phi[pi$on | pi$tsd <= 0] <- 0
      ifelse(pi$on, 1, phi)
    },
    D = , E = rep(0, n))
  idx <- (3 + 2 * n):(2 + 3 * n)
  R[1, idx] <- R[idx, 1] <- cbm * coup
  list(mu = mu, Sigma = outer(sds, sds) * R, R = R, n = n,
       mean_mod = arm %in% c('D', 'E'))
}

draw_path <- function(b, pi, cbm, n_draw) {
  X <- mvrnorm(n_draw, b$mu, b$Sigma)
  n <- b$n
  bm <- X[, 1]; BL <- X[, 2]
  tv <- X[, 3:(2 + n), drop = FALSE]
  pb <- X[, (3 + n):(2 + 2 * n), drop = FALSE]
  br <- X[, (3 + 2 * n):(2 + 3 * n), drop = FALSE]
  if (b$mean_mod) {
    z <- (bm - bm_m) / bm_sd
    br[, pi$on] <- br[, pi$on] + cbm * 8 * z
  }
  Y <- BL - (tv + pb + br)
  list(bm = bm, BL = BL, Y = Y)
}

## ---------------- analysis: published model ----------------
analyze_published <- function(long) {
  fit <- suppressMessages(suppressWarnings(
    lmerTest::lmer(Sx ~ bm + Db + t + bm * Db + (1 | ptID), data = long)))
  cf <- summary(fit)$coefficients
  data.table(beta = cf['bm:DbTRUE', 'Estimate'],
             se = cf['bm:DbTRUE', 'Std. Error'],
             p = cf['bm:DbTRUE', 'Pr(>|t|)'],
             singular = isSingular(fit))
}

## nlme's default optimizer occasionally fails on these designs; the
## project's carryover-sensitivity core falls back to optim the same way.
analyze_car1 <- function(long) {
  fit_with <- function(ctrl) nlme::lme(
    Sx ~ bm + Db + t + bm:Db, random = ~ 1 | ptID,
    correlation = nlme::corCAR1(form = ~ t | ptID),
    data = long, control = ctrl)
  fit <- tryCatch(fit_with(nlme::lmeControl()),
                  error = function(e) fit_with(nlme::lmeControl(opt = 'optim')))
  tt <- summary(fit)$tTable
  data.table(beta = tt['bm:DbTRUE', 'Value'],
             se = tt['bm:DbTRUE', 'Std.Error'],
             p = tt['bm:DbTRUE', 'p-value'],
             singular = NA)
}
analyze <- if (analysis == 'published') analyze_published else analyze_car1

run_cell <- function(cell) {
  set.seed(cell$seed)
  d <- designs[[cell$design]]
  pis <- lapply(d$paths, function(od) path_info(d$t_wk, d$e, od))
  bs <- lapply(pis, function(pi) build(pi, cell$arm, cell$c.bm, cell$t_half))
  min_eig <- min(sapply(bs, function(b) min(eigen(b$R, TRUE, TRUE)$values)))
  if (min_eig <= 0)
    return(list(summary = cbind(cell, feasible = FALSE, min_eig = min_eig,
                                power = NA_real_, mean_beta = NA_real_,
                                emp_sd = NA_real_, mean_se = NA_real_,
                                singular = NA_real_, fits = 0L),
                reps = NULL))
  nP <- length(pis)
  Ns <- N_TOTAL %/% nP + c(rep(1, N_TOTAL %% nP), rep(0, nP - N_TOTAL %% nP))
  draws <- lapply(seq_len(nP), function(k)
    draw_path(bs[[k]], pis[[k]], cell$c.bm, Ns[k] * reps))
  res <- vector('list', reps)
  for (r in seq_len(reps)) {
    parts <- lapply(seq_len(nP), function(k) {
      rows <- ((r - 1) * Ns[k] + 1):(r * Ns[k])
      dr <- draws[[k]]; pi <- pis[[k]]
      id <- (k - 1) * 1000 + seq_len(Ns[k])
      base <- data.table(ptID = id, bm = dr$bm[rows], Sx = dr$BL[rows],
                         t = 0, Db = FALSE)
      post <- data.table(ptID = rep(id, times = length(pi$w)),
                         bm = rep(dr$bm[rows], times = length(pi$w)),
                         Sx = as.vector(dr$Y[rows, , drop = FALSE]),
                         t = rep(pi$w, each = Ns[k]),
                         Db = rep(pi$on, each = Ns[k]))
      rbind(base, post)
    })
    long <- rbindlist(parts)
    t_fit <- Sys.time()
    ## A fit that has not finished within FIT_LIMIT seconds is recorded
    ## as a timeout rather than allowed to stall the whole cell.
    res[[r]] <- tryCatch({
      setTimeLimit(elapsed = FIT_LIMIT, transient = TRUE)
      cbind(analyze(long), status = 'ok')
    }, error = function(e) data.table(
      beta = NA_real_, se = NA_real_, p = NA_real_, singular = NA,
      status = if (grepl('time limit', conditionMessage(e))) 'timeout' else 'error'))
    setTimeLimit(elapsed = Inf, transient = FALSE)
    res[[r]][, fit_sec := as.numeric(difftime(Sys.time(), t_fit, units = 'secs'))]
  }
  rr <- rbindlist(res, fill = TRUE)[, rep := seq_len(.N)]
  rr <- cbind(cell[rep(1, nrow(rr))], rr)
  list(summary = cbind(cell, feasible = TRUE, min_eig = min_eig,
                       power = mean(rr$p < 0.05, na.rm = TRUE),
                       mean_beta = mean(rr$beta, na.rm = TRUE),
                       emp_sd = sd(rr$beta, na.rm = TRUE),
                       mean_se = mean(rr$se, na.rm = TRUE),
                       singular = mean(rr$singular, na.rm = TRUE),
                       fits = sum(!is.na(rr$p)),
                       timeouts = sum(rr$status == 'timeout'),
                       errors = sum(rr$status == 'error'),
                       max_fit_sec = max(rr$fit_sec),
                       total_fit_sec = sum(rr$fit_sec)),
       reps = rr)
}

## Checkpointed wrapper: each finished cell is saved to its own file and
## logged, and a restart skips cells already saved.
run_cell_checkpointed <- function(cell) {
  f <- file.path(cell_dir, sprintf('cell-%03d.rds', cell$id))
  if (file.exists(f)) return(readRDS(f))
  t_cell <- Sys.time()
  out <- run_cell(cell)
  saveRDS(out, f)
  cat(sprintf('%s cell %3d %s %-6s c_bm=%.2f t1/2=%.1f  %.1f min\n',
              format(Sys.time(), '%H:%M:%S'), cell$id, cell$arm, cell$design,
              cell$c.bm, cell$t_half,
              as.numeric(difftime(Sys.time(), t_cell, units = 'mins'))),
      file = progress_file, append = TRUE)
  out
}

## The grid is defined over all arms so that a cell's id and seed never
## depend on which arms a run selects (CJ sorts by arm, so A-E keep ids
## 1-300 and F takes 301-360).
cells <- CJ(arm = ALL_ARMS,
            design = names(designs),
            c.bm = c(0, 0.25, 0.3, 0.6),
            t_half = c(0, 0.1, 0.2, 0.5, 1))
cells[, seed := 20261001L + .I]
cells[, id := .I]
cells <- cells[arm %in% arms]
cell_dir <- file.path(out_dir, sprintf('cells-reps%d%s', reps, tag))
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
progress_file <- file.path(cell_dir, 'progress.log')
done <- length(list.files(cell_dir, pattern = '^cell-.*\\.rds$'))
cat(sprintf('analysis: %s, cells: %d (%d already done), reps per cell: %d, cores: %d, fit limit %ds\n',
            analysis, nrow(cells), done, reps, detectCores(), FIT_LIMIT))
t0 <- Sys.time()
out <- mclapply(split(cells, seq_len(nrow(cells))), run_cell_checkpointed,
                mc.cores = 8, mc.preschedule = FALSE)
bad <- which(sapply(out, inherits, 'try-error'))
if (length(bad)) stop('cells failed: ', paste(bad, collapse = ', '))
summ <- rbindlist(lapply(out, `[[`, 'summary'), fill = TRUE)
cat(sprintf('fits: %d ok, %d timed out, %d errors; slowest fit %.1f s\n',
            sum(summ$fits), sum(summ$timeouts, na.rm = TRUE),
            sum(summ$errors, na.rm = TRUE), max(summ$max_fit_sec, na.rm = TRUE)))
fwrite(summ, file.path(out_dir, sprintf('summary-reps%d%s%s.csv', reps, tag, arm_tag)))
saveRDS(rbindlist(lapply(out, `[[`, 'reps'), fill = TRUE),
        file.path(out_dir, sprintf('replicates-reps%d%s%s.rds', reps, tag, arm_tag)))
cat(sprintf('elapsed: %.1f min\n',
            as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
cat(sprintf('feasible cells: %d of %d\n', sum(summ$feasible), nrow(summ)))
print(dcast(summ, arm + design + c.bm ~ t_half, value.var = 'power'),
      digits = 2)
