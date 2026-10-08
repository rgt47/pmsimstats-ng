## Which schedule features move the c_bm ceiling?
##
## Synthetic schedules vary one feature at a time: number of
## occasions, spacing, heterogeneity of spacing, number of on/off
## switches, and on-drug fraction. Both within-factor structures are
## evaluated with their own cross-factor form: AR(1) with decaying
## cross-factor entries (covar) and CS with constant ones (orig). The
## coupling is the pure on/off pattern (t1/2 = 0), where the graded
## and step rules coincide, so only structure and schedule vary.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/cbm-ceiling/02-ceiling-determinants.R
## Whitepaper: docs/34-cbm-feasible-range.md

library(data.table)
out_dir <- 'analysis/data/quick-sim/cbm-ceiling'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#' Response block M for weeks w under a within/cross structure.
build_M <- function(w, rho, c1, cx, structure) {
  n <- length(w)
  gap <- abs(outer(w, w, '-'))
  if (structure == 'ar1') {
    A <- rho^gap
    C <- cx * rho^gap
  } else {
    A <- matrix(rho, n, n)
    diag(A) <- 1
    C <- matrix(cx, n, n)
  }
  diag(C) <- c1
  rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
}

ceiling_of <- function(w, u, rho = 0.7, c1 = 0.2, cx = 0.1,
                       structure) {
  M <- build_M(w, rho, c1, cx, structure)
  if (min(eigen(M, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    return(NA_real_)
  }
  v <- c(rep(0, 2 * length(w)), u)
  1 / sqrt(drop(t(v) %*% solve(M, v)))
}

#' Closed form under CS: depends on u only through n and sum(u), sum(u^2).
ceiling_cs_closed <- function(u, rho = 0.7, c1 = 0.2, cx = 0.1) {
  n <- length(u)
  J3 <- matrix(1, 3, 3)
  K1 <- (1 - rho) * diag(3) + (c1 - cx) * (J3 - diag(3))
  K2 <- rho * diag(3) + cx * (J3 - diag(3))
  a <- solve(K1)[3, 3]
  b <- solve(K1 + n * K2)[3, 3]
  q <- a * (sum(u^2) - sum(u)^2 / n) + b * sum(u)^2 / n
  1 / sqrt(q)
}

both <- function(w, u, ...) {
  c(ar1 = ceiling_of(w, u, structure = 'ar1', ...),
    cs = ceiling_of(w, u, structure = 'cs', ...))
}

blocks <- function(n, n_switch) {
  len <- n / (n_switch + 1)
  rep(rep(c(1, 0), length.out = n_switch + 1), each = len)
}

res <- list()
add <- function(experiment, label, w, u) {
  cc <- both(w, u)
  res[[length(res) + 1]] <<- data.table(
    experiment = experiment, label = label, n = length(w),
    n_on = sum(u), n_switch = sum(abs(diff(u))),
    span = diff(range(w)), ar1 = cc[['ar1']], cs = cc[['cs']],
    cs_closed = ceiling_cs_closed(u))
}

## E1. Number of occasions: equal spacing, half on, one switch.
for (delta in c(1, 2.5)) {
  for (n in c(4, 8, 16, 32)) {
    add('E1 occasions', sprintf('spacing %.1f wk', delta),
        seq(0, by = delta, length.out = n), blocks(n, 1))
  }
}

## E2. Number of on/off switches: n = 16, half on, equal spacing.
for (delta in c(1, 2.5)) {
  for (s in c(1, 3, 7, 15)) {
    add('E2 switches', sprintf('spacing %.1f wk', delta),
        seq(0, by = delta, length.out = 16), blocks(16, s))
  }
}

## E3. Uniform spacing: n = 8, pattern 11110000.
for (delta in c(0.5, 1, 2.5, 5, 10)) {
  add('E3 spacing', sprintf('spacing %.1f wk', delta),
      seq(0, by = delta, length.out = 8), blocks(8, 1))
}

## E4. Heterogeneity at fixed span (17.5 wk, n = 8, 11110000):
## vary the gap straddling the switch, other six gaps equal.
for (g in c(0.5, 1, 2.5, 5, 8)) {
  others <- (17.5 - g) / 6
  gaps <- c(rep(others, 3), g, rep(others, 3))
  add('E4 switch gap', sprintf('switch gap %.1f wk', g),
      c(0, cumsum(gaps)), blocks(8, 1))
}

## E5. Heterogeneity away from the switch (switch gap fixed 2.5 wk):
## on-block gaps short, off-block gaps long, and the reverse.
het <- list(
  'equal 2.5/2.5' = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5),
  'on short 1, off long 4' = c(1, 1, 1, 2.5, 4, 4, 4),
  'on long 4, off short 1' = c(4, 4, 4, 2.5, 1, 1, 1)
)
for (nm in names(het)) {
  add('E5 within-block gaps', nm, c(0, cumsum(het[[nm]])), blocks(8, 1))
}

## E6. On-drug fraction: n = 8, spacing 2.5, one switch.
for (k in 1:7) {
  add('E6 on fraction', sprintf('%d of 8 on', k),
      seq(0, by = 2.5, length.out = 8), c(rep(1, k), rep(0, 8 - k)))
}

## E7. Real designs against equal spacing over the same span.
real <- list(
  'Hybrid path C' = list(w = c(4, 8, 9, 10, 11, 12, 16, 20),
                         u = c(1, 1, 1, 0, 0, 0, 1, 0)),
  'OL+BDC path A' = list(w = c(4, 8, 12, 16, 17, 18, 19, 20),
                         u = c(1, 1, 1, 1, 1, 1, 0, 0)),
  'CO path A' = list(w = seq(2.5, 20, by = 2.5),
                     u = c(1, 1, 1, 1, 0, 0, 0, 0))
)
for (nm in names(real)) {
  r <- real[[nm]]
  add('E7 real vs equal', paste(nm, 'actual'), r$w, r$u)
  add('E7 real vs equal', paste(nm, 'equalized'),
      seq(min(r$w), max(r$w), length.out = length(r$w)), r$u)
}

## E8. Arrangement at fixed n, k and switches (CS invariance check).
arr <- list(
  '11110000' = c(1, 1, 1, 1, 0, 0, 0, 0),
  '00001111' = c(0, 0, 0, 0, 1, 1, 1, 1),
  '11000011' = c(1, 1, 0, 0, 0, 0, 1, 1)
)
for (nm in names(arr)) {
  add('E8 arrangement', paste(nm, 'Hybrid weeks'),
      c(4, 8, 9, 10, 11, 12, 16, 20), arr[[nm]])
}

out <- rbindlist(res)
cat(sprintf('CS closed form vs numeric, max abs diff: %.2e\n\n',
            max(abs(out$cs - out$cs_closed), na.rm = TRUE)))
print(out[, .(experiment, label, n, n_on, n_switch, span,
              ar1 = round(ar1, 3), cs = round(cs, 3))], nrows = 200)

## E9. Other parameters, on the Hybrid path C schedule.
w <- c(4, 8, 9, 10, 11, 12, 16, 20)
u <- c(1, 1, 1, 0, 0, 0, 1, 0)
par <- rbindlist(lapply(c(0.1, 0.2, 0.3), function(c1)
  rbindlist(lapply(c(0, 0.1, 0.2), function(cx)
    data.table(c1 = c1, cx = cx,
      ar1 = ceiling_of(w, u, 0.7, c1, cx, 'ar1'),
      cs = ceiling_of(w, u, 0.7, c1, cx, 'cs'))))))
cat('\nE9. Cross-factor parameters (rho = 0.7, Hybrid path C)\n')
print(par[cx <= c1], digits = 3)

## E10. Design lookup grid: visit spacing x number of switches,
## 16 occasions, half on drug, both structures.
grid <- CJ(spacing = c(1, 1.5, 2, 2.5, 3.5, 5), n_switch = c(1, 3, 7, 15))
grid[, c('ar1', 'cs') := as.list(both(
  seq(0, by = spacing, length.out = 16), blocks(16, n_switch))),
  by = .(spacing, n_switch)]
cat('\nE10. Lookup grid (16 occasions, 8 on drug)\n')
print(grid, digits = 3)

## E11. CS lookup: the closed form needs only n and the number of
## coupled occasions k, so this grid covers every CS design.
cs_grid <- CJ(n = c(8, 16, 24, 32, 48), frac = (1:7) / 8)
cs_grid[, k := n * frac]
cs_grid[, cs := ceiling_cs_closed(c(rep(1, k), rep(0, n - k))),
        by = .(n, k)]
cat('\nE11. CS lookup (closed form)\n')
print(dcast(cs_grid, n ~ frac, value.var = 'cs'), digits = 3)
fwrite(cs_grid, file.path(out_dir, 'determinants-cs-grid.csv'))

fwrite(out, file.path(out_dir, 'determinants.csv'))
fwrite(par, file.path(out_dir, 'determinants-crossfactor.csv'))
fwrite(grid, file.path(out_dir, 'determinants-grid.csv'))
