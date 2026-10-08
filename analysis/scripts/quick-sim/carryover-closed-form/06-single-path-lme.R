## Monte Carlo check: power of the interaction test when all N = 70
## patients follow one path, against the published allocation over all
## paths (Hybrid 18/18/17/17, CO 35/35). Companion to the closed-form
## comparison of the E9 statistic (docs/37, single-path discussion).
##
## Data: the 04-power-simulation.R DGP, CS construct (rho = 0.8,
## c1 = 0.2, cx = 0.1), with three interaction rules at strength 0.25:
##   step    published coupling (arm A)
##   graded  coupling 1 on drug, 2^(-tsd / t1/2) off drug, CS construct
##   mean    mean moderation (arm D)
## Each replicate is analyzed three ways on the same data:
##   published  lmer Sx ~ bm + Db + t + bm*Db + (1|ptID), Satterthwaite
##   path       the same plus a fixed path (sequence) effect, mixed
##              allocations only
##   e9         OLS of the on-minus-off contrast on bm
## and the E9 closed form of 01-e9-closed-form.R is attached per cell.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/06-single-path-lme.R \
##     [--reps R] [--cores K]
suppressMessages({
  library(data.table); library(MASS); library(lmerTest); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
reps <- arg('--reps', 1000L)
cores <- arg('--cores', 8L)
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
here01 <- 'analysis/scripts/quick-sim/carryover-closed-form/01-e9-closed-form.R'
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/single-path-lme'
cell_dir <- file.path(out_dir, sprintf('cells-reps%d', reps))
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build', 'draw_path')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))
CBM <- 0.25

## The graded rule under CS is not one of the 04 arms: take the CS
## matrix of arm A and replace its biomarker row by the graded pattern.
build_rule <- function(pi, rule, th) {
  if (rule == 'step') return(dgp$build(pi, 'A', CBM, th))
  if (rule == 'mean') return(dgp$build(pi, 'D', CBM, th))
  b <- dgp$build(pi, 'A', 0, th)
  n <- b$n
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, n)
  phi[pi$on | pi$tsd <= 0] <- 0
  idx <- (3 + 2 * n):(2 + 3 * n)
  b$R[1, idx] <- b$R[idx, 1] <- CBM * ifelse(pi$on, 1, phi)
  sds <- sqrt(diag(b$Sigma))
  b$Sigma <- outer(sds, sds) * b$R
  b
}

fit_lmer <- function(long, with_path) {
  f <- if (with_path) Sx ~ bm + Db + t + bm * Db + factor(path) + (1 | ptID)
       else Sx ~ bm + Db + t + bm * Db + (1 | ptID)
  out <- tryCatch({
    fit <- suppressMessages(suppressWarnings(lmerTest::lmer(f, data = long)))
    cf <- summary(fit)$coefficients
    c(beta = cf['bm:DbTRUE', 'Estimate'], p = cf['bm:DbTRUE', 'Pr(>|t|)'])
  }, error = function(e) c(beta = NA_real_, p = NA_real_))
  out
}

## E9: on-drug mean minus off-drug mean over the post-baseline visits.
fit_e9 <- function(bm, Y, on) {
  d <- rowMeans(Y[, on, drop = FALSE]) - rowMeans(Y[, !on, drop = FALSE])
  cf <- summary(lm(d ~ bm))$coefficients
  c(beta = cf['bm', 'Estimate'], p = cf['bm', 'Pr(>|t|)'])
}

## Allocations: 'mix' is the published split; a letter puts all N on
## that path.
alloc_of <- function(design, a) {
  nP <- length(dgp$designs[[design]]$paths)
  if (a == 'mix') return(N_TOTAL %/% nP + c(rep(1L, N_TOTAL %% nP), rep(0L, nP - N_TOTAL %% nP)))
  k <- match(a, LETTERS)
  replace(integer(nP), k, N_TOTAL)
}
N_TOTAL <- dgp$N_TOTAL

run_cell <- function(cell) {
  set.seed(cell$seed)
  d <- dgp$designs[[cell$design]]
  pis <- lapply(d$paths, function(od) dgp$path_info(d$t_wk, d$e, od))
  Ns <- alloc_of(cell$design, cell$alloc)
  use <- which(Ns > 0)
  bs <- lapply(pis, build_rule, rule = cell$rule, th = cell$t_half)
  min_eig <- min(sapply(bs[use], function(b) min(eigen(b$R, TRUE, TRUE)$values)))
  stopifnot(min_eig > 0)
  draws <- lapply(seq_along(pis), function(k)
    if (Ns[k] > 0) dgp$draw_path(bs[[k]], pis[[k]], CBM, Ns[k] * reps) else NULL)
  res <- matrix(NA_real_, reps, 6,
                dimnames = list(NULL, c('beta_pub', 'p_pub', 'beta_path', 'p_path',
                                        'beta_e9', 'p_e9')))
  for (r in seq_len(reps)) {
    parts <- lapply(use, function(k) {
      rows <- ((r - 1) * Ns[k] + 1):(r * Ns[k])
      dr <- draws[[k]]; pi <- pis[[k]]
      id <- (k - 1) * 1000 + seq_len(Ns[k])
      Y <- dr$Y[rows, , drop = FALSE]
      list(long = rbind(
        data.table(ptID = id, path = k, bm = dr$bm[rows], Sx = dr$BL[rows],
                   t = 0, Db = FALSE),
        data.table(ptID = rep(id, times = length(pi$w)), path = k,
                   bm = rep(dr$bm[rows], times = length(pi$w)),
                   Sx = as.vector(Y), t = rep(pi$w, each = Ns[k]),
                   Db = rep(pi$on, each = Ns[k]))),
        bm = dr$bm[rows], Y = Y, on = pi$on)
    })
    long <- rbindlist(lapply(parts, `[[`, 'long'))
    res[r, 1:2] <- fit_lmer(long, FALSE)
    if (length(use) > 1) res[r, 3:4] <- fit_lmer(long, TRUE)
    ## E9 pooled over paths with a common intercept, as in the closed form.
    dd <- rbindlist(lapply(parts, function(q) data.table(
      bm = q$bm, d = rowMeans(q$Y[, q$on, drop = FALSE]) -
        rowMeans(q$Y[, !q$on, drop = FALSE]))))
    cf <- summary(lm(d ~ bm, data = dd))$coefficients
    res[r, 5:6] <- c(cf['bm', 'Estimate'], cf['bm', 'Pr(>|t|)'])
  }
  rr <- as.data.table(res)
  cbind(cell, min_eig = min_eig,
        power_pub = mean(rr$p_pub < 0.05, na.rm = TRUE),
        power_path = if (length(use) > 1) mean(rr$p_path < 0.05, na.rm = TRUE) else NA_real_,
        power_e9 = mean(rr$p_e9 < 0.05, na.rm = TRUE),
        mean_beta_pub = mean(rr$beta_pub, na.rm = TRUE),
        sd_beta_pub = sd(rr$beta_pub, na.rm = TRUE),
        fails_pub = sum(is.na(rr$p_pub)), fits = reps)
}

run_cell_checkpointed <- function(cell) {
  f <- file.path(cell_dir, sprintf('cell-%03d.rds', cell$id))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  out <- run_cell(cell)
  saveRDS(out, f)
  cat(sprintf('%s cell %3d %-6s %-4s %-6s t1/2=%.1f  %.1f min\n',
              format(Sys.time(), '%H:%M:%S'), cell$id, cell$design, cell$alloc,
              cell$rule, cell$t_half,
              as.numeric(difftime(Sys.time(), t0, units = 'mins'))),
      file = file.path(cell_dir, 'progress.log'), append = TRUE)
  out
}

cells <- rbind(CJ(design = 'Hybrid', alloc = c('mix', 'A', 'B', 'C', 'D'),
                  rule = c('step', 'graded', 'mean'), t_half = c(0, 1, 2), sorted = FALSE),
               CJ(design = 'CO', alloc = c('mix', 'A', 'B'),
                  rule = c('step', 'graded', 'mean'), t_half = c(0, 1, 2), sorted = FALSE))
cells[, id := .I]
cells[, seed := 20261003L + id]
cat(sprintf('cells: %d, reps: %d, cores: %d\n', nrow(cells), reps, cores))
t0 <- Sys.time()
out <- mclapply(split(cells, seq_len(nrow(cells))), run_cell_checkpointed,
                mc.cores = cores, mc.preschedule = FALSE)
bad <- which(sapply(out, inherits, 'try-error'))
if (length(bad)) stop('cells failed: ', paste(bad, collapse = ', '))
summ <- rbindlist(out)

## Closed form of the E9 statistic for each cell (01-e9-closed-form.R
## definitions, evaluated in their own environment).
cf <- new.env()
## use_design() sets these with <<-; they must exist in cf beforehand or
## the assignment lands in the global environment.
for (v in c('nP', 'Ns', 'paths')) assign(v, NULL, envir = cf)
ex <- parse(here01)
stop_at <- which(vapply(ex, function(e) is.call(e) && identical(e[[1]], as.name('<-')) &&
                          identical(e[[2]], as.name('pp')), logical(1)))[1]
for (e in ex[seq_len(stop_at - 1)]) eval(e, cf)
summ[, power_e9_closed := mapply(function(dn, a, rule, th) {
  cf$use_design(dn)
  m <- rbindlist(lapply(cf$paths, cf$path_moments, rule = rule, th = th,
                        variant = 'E9', construct = 'CS'))
  assign('Ns', alloc_of(dn, a), envir = cf)
  cf$pooled(m)$power
}, design, alloc, rule, t_half)]
summ[, mcse_pub := sqrt(power_pub * (1 - power_pub) / fits)]
summ[, z_e9 := (power_e9 - power_e9_closed) /
         sqrt(pmax(power_e9_closed * (1 - power_e9_closed), 1e-12) / fits)]

fwrite(summ, file.path(out_dir, sprintf('summary-reps%d.csv', reps)))
cat(sprintf('elapsed: %.1f min\n', as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
cat(sprintf('lmer failures: %d of %d fits; E9 MC vs closed form: max |z| %.2f\n',
            sum(summ$fails_pub), sum(summ$fits), max(abs(summ$z_e9))))
for (v in c('power_pub', 'power_path', 'power_e9', 'power_e9_closed')) {
  cat(sprintf('\n%s\n', v))
  print(dcast(summ, design + rule + t_half ~ alloc, value.var = v)[
    , c('design', 'rule', 't_half', 'mix', 'A', 'B', 'C', 'D'), with = FALSE], digits = 3)
}
