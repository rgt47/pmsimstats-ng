## 24-rho-sweep.R
##
## Closed-form sweep of the serial correlation rho for paper 01's
## planning tables (referee report 2026-11-05, M1). Under the extended
## DGP (separable AR(1), c_1 = 0.2, graded coupling, anchored carryover)
## it computes, for rho from 0.4 to 0.95:
##
##   1. the feasibility bound on c_bm by design and carryover half-life;
##   2. the total N for 80% power of E9 and E9 power at N = 70, by
##      design, c_bm and half-life, with each cell flagged feasible or
##      not against its own bound.
##
## The functions are those of 17-rho-stage.R. At rho = 0.7 the N80
## values must reproduce rho-stage/planning-n80.csv; the script stops
## otherwise. Paper 01's inputs are not overwritten.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/24-rho-sweep.R
suppressMessages(library(data.table))
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
in_dir <- 'analysis/data/quick-sim/carryover-closed-form/rho-stage'
out_dir <- file.path(in_dir, 'rho-sweep')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))

ALPHA <- 0.05
sd_B <- dgp$bm_sd
sd_BR <- 8
ARM <- c(CS = 'A', two = 'F', sep = 'C')
DESIGNS <- c('Hybrid', 'OL+BDC', 'CO')
RHOS <- c(0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
CBMS <- c(0.25, 0.35, 0.45)
HALVES <- c(0, 0.2, 1)
N_MAX <- 2000L

design_paths <- function(name) {
  d <- dgp$designs[[name]]
  lapply(d$paths, function(od) dgp$path_info(d$t_wk, d$e, od))
}
alloc <- function(N, nP) N %/% nP + c(rep(1L, N %% nP), rep(0L, nP - N %% nP))
graded <- function(pi, th) {
  n <- length(pi$w)
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, n)
  phi[pi$on | pi$tsd <= 0] <- 0
  ifelse(pi$on, 1, phi)
}
response <- function(pi, construct, rho, th) {
  dgp$RHO <- rho
  b <- dgp$build(pi, ARM[[construct]], 0, th)
  n <- length(pi$w)
  idx <- 3:(2 + 3 * n)
  list(M = b$R[idx, idx], Sigma = b$Sigma[idx, idx], mu = b$mu[idx], n = n)
}
ceiling_of <- function(design, construct, rho, th = 0) {
  min(vapply(design_paths(design), function(pi) {
    r <- response(pi, construct, rho, th)
    u <- c(rep(0, 2 * r$n), graded(pi, th))
    if (!any(u != 0)) return(Inf)
    1 / sqrt(drop(crossprod(u, solve(r$M, u))))
  }, numeric(1)))
}

## The response moments do not depend on N, so they are computed once
## per cell and power is then evaluated cheaply over N.
cell_moments <- function(design, construct, rho, th, cbm) {
  rbindlist(lapply(design_paths(design), function(pi) {
    r <- response(pi, construct, rho, th)
    n <- r$n
    L <- cbind(diag(n), diag(n), diag(n))
    C_S <- L %*% r$Sigma %*% t(L)
    mu_S <- drop(L %*% r$mu)
    a <- ifelse(pi$on, 1 / sum(pi$on), -1 / sum(!pi$on))
    gamma <- -cbm * sd_BR / sd_B * sum(a * graded(pi, th))
    var_D <- drop(t(a) %*% C_S %*% a)
    data.table(gamma = gamma, tau2 = var_D - gamma^2 * sd_B^2,
               mu = -sum(a * mu_S))
  }))
}
power_at <- function(m, N) {
  w <- alloc(N, nrow(m)) / N
  gbar <- sum(w * m$gamma); mubar <- sum(w * m$mu)
  resid <- sum(w * (m$tau2 + (m$mu - mubar)^2 + sd_B^2 * (m$gamma - gbar)^2))
  chi_q <- qchisq((seq_len(400) - 0.5) / 400, N - 1)
  tc <- qt(1 - ALPHA / 2, N - 2)
  ncp <- gbar * sd_B * sqrt(chi_q) / sqrt(resid)
  mean(pt(-tc, N - 2, ncp = ncp) + pt(tc, N - 2, ncp = ncp, lower.tail = FALSE))
}
## Smallest even N from 20 with power >= 0.8, as in 17-rho-stage.R.
n80_of <- function(m) {
  Ns <- seq(20L, N_MAX, by = 2L)
  if (power_at(m, N_MAX) < 0.8) return(NA_integer_)
  lo <- 1L; hi <- length(Ns)
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (power_at(m, Ns[mid]) >= 0.8) hi <- mid else lo <- mid + 1L
  }
  Ns[lo]
}

## ---------------- 1. feasibility bounds ----------------
ceil <- CJ(rho = RHOS, design = DESIGNS, t_half = HALVES)
ceil[, ceiling := mapply(ceiling_of, design, 'sep', rho, t_half)]
stopifnot(abs(ceil[rho == 0.8 & t_half == 0, min(ceiling)] - 0.440) < 1e-3,
          abs(ceil[rho == 0.7 & t_half == 0, min(ceiling)] - 0.4741) < 1e-3)
fwrite(ceil, file.path(out_dir, 'ceilings-sep.csv'))

## ---------------- 2. N80 and power at N = 70 ----------------
plan <- CJ(rho = RHOS, design = DESIGNS, c_bm = CBMS, t_half = HALVES)
plan <- merge(plan, ceil, by = c('rho', 'design', 't_half'))
plan[, feasible := c_bm < ceiling]
plan[, c('N80', 'power_N70') := {
  r <- lapply(seq_len(.N), function(i) {
    if (!feasible[i]) return(list(NA_integer_, NA_real_))
    m <- cell_moments(design[i], 'sep', rho[i], t_half[i], c_bm[i])
    list(n80_of(m), power_at(m, 70L))
  })
  list(vapply(r, `[[`, integer(1), 1), vapply(r, `[[`, numeric(1), 2))
}]

## Reproduce the stored rho = 0.7 planning table (N80 capped at 800 there).
old <- fread(file.path(in_dir, 'planning-n80.csv'))
chk <- merge(old, plan[rho == 0.7], by = c('design', 'c_bm', 't_half'),
             suffixes = c('.old', '.new'))
stopifnot(nrow(chk) == 27,
          all(abs(chk$power_N70.old - chk$power_N70.new) < 1e-8),
          all(is.na(chk$N80.old) == (is.na(chk$N80.new) | chk$N80.new > 800)),
          all(chk$N80.old == chk$N80.new, na.rm = TRUE))
cat('rho = 0.7 reproduces planning-n80.csv (27 cells).\n')
fwrite(plan, file.path(out_dir, 'planning-n80-by-rho.csv'))

## ---------------- 3. planning files for driver 19 ----------------
## Paper 01 extends its planning tables to rho 0.5 and 0.9 (2026-11-06).
## At 0.9 the grid is moved inside the bound (0.20, 0.25, 0.30); at 0.5
## the reference grid is kept and infeasible cells are flagged, so that
## driver 19 (--plan-file, --rho) skips them.
PLAN_GRIDS <- list(`0.5` = c(0.25, 0.35, 0.45), `0.9` = c(0.20, 0.25, 0.30))
for (r in names(PLAN_GRIDS)) {
  rr <- as.numeric(r)
  p <- CJ(design = DESIGNS, c_bm = PLAN_GRIDS[[r]], t_half = HALVES)
  p <- merge(p, ceil[rho == rr, .(design, t_half, ceiling)],
             by = c('design', 't_half'))
  p[, feasible := c_bm < ceiling]
  p[, c('N80', 'power_N70') := {
    r2 <- lapply(seq_len(.N), function(i) {
      if (!feasible[i]) return(list(NA_integer_, NA_real_))
      m <- cell_moments(design[i], 'sep', rr, t_half[i], c_bm[i])
      n <- n80_of(m)
      ## Same cap as planning-n80.csv, which driver 19 expects.
      list(if (!is.na(n) && n > 800L) NA_integer_ else n, power_at(m, 70L))
    })
    list(vapply(r2, `[[`, integer(1), 1), vapply(r2, `[[`, numeric(1), 2))
  }]
  p[, slope_per_sd := -c_bm * sd_BR]
  setcolorder(p, c('design', 'c_bm', 't_half', 'N80', 'power_N70',
                   'slope_per_sd', 'ceiling', 'feasible'))
  f <- file.path(in_dir, sprintf('planning-n80-rho%s.csv', r))
  fwrite(p[order(design, c_bm, t_half)], f)
  cat(sprintf('\nWrote %s (%d cells, %d feasible)\n', f, nrow(p),
              sum(p$feasible)))
  print(dcast(p, design + c_bm ~ t_half, value.var = 'N80'))
}

## ---------------- printed summaries ----------------
cat('\nFeasibility bound on c_bm, separable AR(1), by design and half-life\n')
print(dcast(ceil, rho ~ design + t_half, value.var = 'ceiling'), digits = 3)
cat('\nMinimum over designs\n')
print(ceil[, .(min = min(ceiling), binding = design[which.min(ceiling)]),
           by = .(rho, t_half)][order(t_half, rho)], digits = 3)
fmt <- function(x, f) ifelse(!f, 'infeas', ifelse(is.na(x), '>2000',
                                                    as.character(x)))
plan[, N80_txt := fmt(N80, feasible)]
for (th in HALVES) {
  cat(sprintf('\nN for 80%% power (E9), t1/2 = %g week\n', th))
  print(dcast(plan[t_half == th], design + c_bm ~ rho, value.var = 'N80_txt'))
}
cat('\nE9 power at N = 70, t1/2 = 0.2 week\n')
print(dcast(plan[t_half == 0.2], design + c_bm ~ rho,
            value.var = 'power_N70'), digits = 2)
