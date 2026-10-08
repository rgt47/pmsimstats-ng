## Check of the pre-exposure coding of Dbc used by manuscript 02
## (analysis/scripts/carryover-sensitivity/simulation-core.R, line 179):
## Dbc = 1 on drug, carryover_decay(tsd, h) off drug. Because tsd is zero
## before first exposure (R/buildtrialdesign.R, line 105) and
## carryover_decay(0, h) = 1 for h > 0, off-drug visits before any
## exposure are coded Dbc = 1. In CO this makes Dbc = 1 at every visit of
## the placebo-first path. The corrected coding sets those visits to 0.
##
## Data as in 08-covariance-study.R (CS construct, c_bm 0.25; configuration
## C at rho 0.7, c_bm 0.45; graded coupling; N = 70), alternative only,
## t1/2 = 1, analysis half-life h = 1. As in manuscript 02 there is no
## baseline row. Analyses: binary Db, Dbc corrected, Dbc as in 02; each
## with the random intercept and corCAR1, model-based (02 G1/G3) and CR2
## (02 G6/G8) standard errors. Hybrid has no off-drug visit before first
## exposure after baseline, so its two Dbc codings must coincide.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/10-dbc-preexposure-check.R \
##     [--reps R] [--cores K]
suppressMessages({
  library(data.table); library(MASS); library(nlme); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
reps <- arg('--reps', 200L)
cores <- arg('--cores', 8L)
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/dbc-preexposure'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build', 'draw_path')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))

## As build_cell() in 08-covariance-study.R.
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

dbc_corrected <- function(Db, tsd, h) ifelse(Db, 1, ifelse(tsd > 0, 0.5^(tsd / h), 0))
dbc_02 <- function(Db, tsd, h) ifelse(Db, 1, 0.5^(tsd / h))

fit <- function(long, x) {
  long <- copy(long)
  long[, X := x]
  m <- tryCatch(lme(Sx ~ bmc + t + X + bmc:X, random = ~ 1 | ptID,
                    correlation = corCAR1(form = ~ t | ptID), data = long,
                    control = lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)),
                error = function(e) NULL)
  if (is.null(m)) return(rep(NA_real_, 4))
  tt <- summary(m)$tTable['bmc:X', ]
  ct <- clubSandwich::coef_test(m, vcov = 'CR2', test = 'Satterthwaite')
  k <- which(ct$Coef == 'bmc:X')
  c(tt[['Value']], tt[['p-value']], ct$SE[k], ct$p_Satt[k])
}

cells <- CJ(construct = c('CS', 'C07'), design = c('Hybrid', 'CO'), sorted = FALSE)
cells[, `:=`(rho = ifelse(construct == 'CS', 0.8, 0.7),
             cbm = ifelse(construct == 'CS', 0.25, 0.45), t_half = 1, h = 1,
             seed = 20261007L + .I)]

run_cell <- function(k) {
  cl <- cells[k]
  set.seed(cl$seed)
  assign('RHO', cl$rho, envir = dgp)
  d <- dgp$designs[[cl$design]]
  pis <- lapply(d$paths, function(od) dgp$path_info(d$t_wk, d$e, od))
  nP <- length(pis)
  N <- dgp$N_TOTAL
  Ns <- N %/% nP + c(rep(1L, N %% nP), rep(0L, nP - N %% nP))
  bs <- lapply(pis, build_cell, construct = cl$construct, cbm = cl$cbm, th = cl$t_half)
  draws <- lapply(seq_len(nP), function(j)
    dgp$draw_path(bs[[j]], pis[[j]], cl$cbm, Ns[j] * reps))
  rbindlist(lapply(seq_len(reps), function(r) {
    long <- rbindlist(lapply(seq_len(nP), function(j) {
      rows <- ((r - 1) * Ns[j] + 1):(r * Ns[j])
      pi <- pis[[j]]; dr <- draws[[j]]
      data.table(ptID = rep((j - 1) * 1000 + seq_len(Ns[j]), times = length(pi$w)),
                 bm = rep(dr$bm[rows], times = length(pi$w)),
                 Sx = as.vector(dr$Y[rows, , drop = FALSE]),
                 t = rep(pi$w, each = Ns[j]), Db = rep(pi$on, each = Ns[j]),
                 tsd = rep(pi$tsd, each = Ns[j]))
    }))
    long[, `:=`(ptID = factor(ptID), bmc = bm - mean(bm))]
    codings <- list(Db = as.numeric(long$Db),
                    Dbc_corrected = dbc_corrected(long$Db, long$tsd, cl$h),
                    Dbc_02 = dbc_02(long$Db, long$tsd, cl$h))
    rbindlist(lapply(names(codings), function(nm) {
      v <- fit(long, codings[[nm]])
      data.table(construct = cl$construct, design = cl$design, rep = r, coding = nm,
                 est = v[1], p_model = v[2], se_cr2 = v[3], p_cr2 = v[4],
                 n_const = sum(tapply(codings[[nm]], long$ptID, function(z) diff(range(z)) == 0)))
    }))
  }))
}
out <- rbindlist(mclapply(seq_len(nrow(cells)), run_cell, mc.cores = min(cores, nrow(cells))))
fwrite(out, file.path(out_dir, sprintf('replicates-reps%d.csv', reps)))
summ <- out[, .(reps = .N, fails = sum(is.na(est)), mean_est = mean(est, na.rm = TRUE),
                emp_sd = sd(est, na.rm = TRUE),
                power_model = mean(p_model < 0.05, na.rm = TRUE),
                power_cr2 = mean(p_cr2 < 0.05, na.rm = TRUE),
                constant_participants = mean(n_const)),
            by = .(construct, design, coding)]
fwrite(summ, file.path(out_dir, sprintf('summary-reps%d.csv', reps)))
print(summ, digits = 3)
