## Figures reconciled when the 2026-10-08 inventory of docs/09a to docs/52
## was merged into analysis/report/compendium-mathematics.tex.
## Run from the repository root: Rscript analysis/scripts/compendium-mathematics/check-reconciliations.R

## 1. Drug-expectancy collinearity in the Hybrid long data
## (Section "The published Hybrid construct"): r, r^2, VIF with and
## without linear time, with and without the baseline row.
w <- c(0, 4, 8, 9, 10, 11, 12, 16, 20)
eta <- c(0, 1, 1, rep(0.5, 6))
D <- rbind(
  c(0, 1, 1, 1, 1, 0, 0, 1, 0),
  c(0, 1, 1, 1, 1, 0, 0, 0, 1),
  c(0, 1, 1, 1, 0, 0, 0, 1, 0),
  c(0, 1, 1, 1, 0, 0, 0, 0, 1)
)
n_p <- c(18, 18, 17, 17)
dat <- do.call(rbind, lapply(1:4, function(p) {
  data.frame(p = p, w = w, eta = eta, D = D[p, ], wt = n_p[p])
}))
r2w <- function(f, d) {
  m <- lm(f, data = d, weights = wt)
  y <- model.response(model.frame(m))
  yb <- weighted.mean(y, d$wt)
  1 - sum(d$wt * resid(m)^2) / sum(d$wt * (y - yb)^2)
}
wcor <- function(x, y, wt) {
  mx <- weighted.mean(x, wt)
  my <- weighted.mean(y, wt)
  sum(wt * (x - mx) * (y - my)) /
    sqrt(sum(wt * (x - mx)^2) * sum(wt * (y - my)^2))
}
show <- function(lab, d) {
  r <- wcor(d$D, d$eta, d$wt)
  a <- r2w(eta ~ D + w, d)
  b <- r2w(D ~ eta + w, d)
  cat(sprintf(
    '%s: r = %.4f, r^2 = %.4f, 1/(1 - r^2) = %.3f; R2(eta | D, w) = %.4f, VIF %.3f; R2(D | eta, w) = %.4f, VIF %.3f\n',
    lab, r, r^2, 1 / (1 - r^2), a, 1 / (1 - a), b, 1 / (1 - b)))
  r
}
r9 <- show('nine rows', dat)
r8 <- show('post-baseline rows', dat[dat$w > 0, ])
stopifnot(abs(r9 - 0.588) < 5e-4, abs(r8 - 0.507) < 5e-4)

## Equal path weights, with D_bc at t_half = 1 week, baseline row added.
dbc <- function(on, wk, th) {
  out <- numeric(length(on)); last <- NA
  for (j in seq_along(on)) {
    if (on[j] == 1) { out[j] <- 1; last <- wk[j] }
    else if (!is.na(last)) out[j] <- 2^(-(wk[j] - last) / th)
  }
  out
}
xb <- unlist(lapply(1:4, function(p) dbc(D[p, ], w, 1)))
eb <- rep(eta, 4)
cat(sprintf('equal weights, baseline row included: binary %.4f; D_bc(t1/2 = 1) %.4f\n',
  cor(as.vector(t(D)), eb), cor(xb, eb)))
keep <- rep(w > 0, 4)
cat(sprintf('equal weights: D_bc(t1/2 = 1) without baseline %.4f\n',
  cor(xb[keep], eb[keep])))

## 2. Week-12 share of the TV and PB curves (Table of Gompertz landmarks).
g <- function(t, d, r) (exp(-d * exp(-r * t)) - exp(-d)) / (1 - exp(-d))
cat(sprintf('share at week 12: BR %.4f, TV/PB %.4f\n', g(12, 5, 0.42),
  g(12, 5, 0.35)))
stopifnot(abs(g(12, 5, 0.35) - 0.927) < 5e-4)

## 3. Rounding of the Hybrid covar (F) ceiling at rho = 0.5, t_half = 0.
csv <- 'analysis/data/quick-sim/carryover-closed-form/ceilings-graded-by-rho.csv'
if (file.exists(csv)) {
  cg <- read.csv(csv)
  v <- cg$ceiling[cg$rho == 0.5 & cg$design == 'Hybrid' &
    cg$t_half == 0 & grepl('covar', cg$construct)]
  cat(sprintf('Hybrid F ceiling, rho 0.5: %.5f (rounds to %.3f)\n', v,
    round(v, 3)))
  stopifnot(round(v, 3) == 0.489)
}
cat('done\n')
