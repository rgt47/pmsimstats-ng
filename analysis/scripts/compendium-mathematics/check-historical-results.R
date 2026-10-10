## Copied 2026-10-08 from the review scratch script inv_a_check.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-historical-results.R
# Checks for inventory part a (docs/09a, 19, 21, 22, 22, 23).
# Self-contained; runnable from the repo root as Rscript analysis/scripts/compendium-mathematics/check-historical-results.R
set.seed(20261008)
say <- function(...) cat(sprintf(...), '\n')

say('--- A. docs/09a: estimand and retention table')
a <- 0.5
beta_bm <- 0.3
say('beta_bm:D = -a beta_bm = %.3f', -a * beta_bm)
h_grid <- c(0.1, 0.2, 0.5, 1, 2)
tsd <- c(0, 1, 3, 7, 11)
ret <- outer(tsd, h_grid, function(s, h) 2^(-s / h))
dimnames(ret) <- list(paste0('tsd=', tsd), paste0('h=', h_grid))
print(round(ret, 3))

say('--- A2. docs/09a: integrated squared misspecification V(h)')
v_fixed <- function(h, T, A = 1) {
  k <- log(2) / h
  A^2 * (T - 2 * (1 - exp(-k * T)) / k + (1 - exp(-2 * k * T)) / (2 * k))
}
v_ls <- function(h, T, A = 1) {
  k <- log(2) / h
  i1 <- A * (1 - exp(-k * T)) / k
  i2 <- A^2 * (1 - exp(-2 * k * T)) / (2 * k)
  i2 - i1^2 / T
}
T <- 11
num_fixed <- function(h) integrate(function(s) (2^(-s / h) - 1)^2, 0, T)$value
num_ls <- function(h) {
  m <- integrate(function(s) 2^(-s / h), 0, T)$value / T
  integrate(function(s) (2^(-s / h) - m)^2, 0, T)$value
}
for (h in c(0.1, 0.5, 2, 10)) {
  say('h=%5.2f  fixed closed %.6f numeric %.6f | ls closed %.6f numeric %.6f',
      h, v_fixed(h, T), num_fixed(h), v_ls(h, T), num_ls(h))
}
hh <- exp(seq(log(1e-3), log(1e3), length.out = 4000))
vf <- v_fixed(hh, 1)
say('fixed case monotone decreasing in h: %s; limits %.4f, %.2e',
    all(diff(vf) < 0), v_fixed(1e-6, 1), v_fixed(1e6, 1))
opt <- optimize(function(h) -v_ls(h, 1), c(1e-4, 100))
say('LS case: argmax h/T = %.4f, max V/(A^2 T) = %.4f; limits %.2e, %.2e',
    opt$minimum, -opt$objective, v_ls(1e-6, 1), v_ls(1e6, 1))

say('--- B. two-period crossover estimator with carryover')
mu <- 3; p1 <- 0.4; p2 <- -0.2
tA1 <- -2; tB1 <- 0; tA2 <- -1.5; tB2 <- 0.1
chi1 <- 0.7; chi0 <- 0.2
y11 <- mu + p1 + tA1; y12 <- mu + p2 + tB2 + chi1
y21 <- mu + p1 + tB1; y22 <- mu + p2 + tA2 + chi0
est <- 0.5 * ((y11 - y12) - (y21 - y22))
say('estimator %.4f; formula 0.5(b1+b2-chi1+chi0) %.4f; source form %.4f',
    est, 0.5 * ((tA1 - tB1) + (tA2 - tB2) - chi1 + chi0),
    0.5 * ((tA1 - tB1) + (tA2 - tB2) - chi1 - chi0))

say('--- C. docs/19 failed formulations')
expo <- 2^(-c(2.5, 5, 7.5, 10) / 1)
say('attempt 2 exposures: %s', paste(signif(expo, 3), collapse = ', '))
say('attempt 2 shifts per SD: %s',
    paste(signif(0.45 * expo * 8, 3), collapse = ', '))
mu_bm <- 124.33; sd_bm <- 15.36
say('attempt 3 at tod=10: total %.1f, biomarker-free part %.1f, per-SD %.2f',
    124 * 10 * 0.45, 0.45 * mu_bm * 10, 0.45 * sd_bm * 10)
say('per-SD shift / sigma_BR = %.2f; / calibrated 3.6 = %.2f',
    0.45 * sd_bm * 10 / 8, 0.45 * sd_bm * 10 / 3.6)

say('--- D. docs/19 calibration: regression of mean on-drug BR on b')
rho <- 0.7; sig <- 8; cbm <- 0.45
w <- c(2.5, 5, 7.5, 10)
Abr <- rho^abs(outer(w, w, '-'))
var_mean <- sig^2 * sum(Abr) / 16
say('Var(mean of 4 on-drug BR) = %.3f, SD %.3f', var_mean, sqrt(var_mean))
var_b <- var_mean - (cbm * sig)^2
say('covariance moderation residual var %.3f, SD %.3f', var_b, sqrt(var_b))
sim_once <- function(arch, n = 500) {
  b <- rnorm(n)
  L <- chol(Abr)
  if (arch == 'mean') {
    X <- matrix(rnorm(n * 4), n) %*% L * sig
    br <- rowMeans(X) + cbm * sig * b
  } else {
    S <- rbind(c(1, rep(cbm, 4)), cbind(rep(cbm, 4), Abr))
    Z <- matrix(rnorm(n * 5), n) %*% chol(S)
    b <- Z[, 1]
    br <- rowMeans(Z[, -1]) * sig
  }
  f <- lm(br ~ b)
  c(coef(f)[2], sigma(f))
}
for (arch in c('mean', 'cov')) {
  r <- t(replicate(2000, sim_once(arch)))
  say('%s: slope mean %.3f sd %.3f; resid SD mean %.3f sd %.3f',
      arch, mean(r[, 1]), sd(r[, 1]), mean(r[, 2]), sd(r[, 2]))
}
say('single-visit residual SDs: %.3f, %.3f', sig, sig * sqrt(1 - cbm^2))

say('--- E. docs/23 4.3 residual correlation after a random intercept')
sT2 <- 10; su2 <- 3; se2 <- sT2 - su2
for (hl in c(0.25, 1, 4)) {
  direct <- (sT2 * rho^hl - su2) / se2
  formula <- rho^hl - (1 - rho^hl) * su2 / se2
  say('lag %.2f: direct %.5f formula %.5f; sign change at rho^h = %.3f',
      hl, direct, formula, su2 / sT2)
}

say('--- F. central 95%% binomial range, 200 replicates, size 0.05')
say('[%.3f, %.3f]', qbinom(0.025, 200, 0.05) / 200,
    qbinom(0.975, 200, 0.05) / 200)

say('--- G. docs/23 5.5-5.7 conditional means')
mu_on <- 10.19; mu_off <- mu_on / 2
say('on: %.2f %.2f; off graded: %.2f %.2f; ratios %.3f %.3f',
    mu_on + 0.5 * 8, mu_on - 0.5 * 8, mu_off + 0.25 * 8, mu_off - 0.25 * 8,
    8 / mu_on, 4 / mu_off)
say('off slow decay lambda=1: r %.3f, %.2f %.2f', 0.5 * exp(-1),
    mu_off + 0.5 * exp(-1) * 8, mu_off - 0.5 * exp(-1) * 8)
say('off full: %.2f %.2f', mu_off + 4, mu_off - 4)
m01 <- mu_on * 2^(-10)
say('t1/2=0.1: mean %.4f; step %.2f %.2f; graded %.4f %.4f; B at +-1SD %.1f %.1f',
    m01, m01 - 0.6 * 8, m01 + 0.6 * 8, m01 - 0.6 * 2^-10 * 8,
    m01 + 0.6 * 2^-10 * 8, mu_bm - sd_bm, mu_bm + sd_bm)

say('--- H. within-participant sum of squares of a binary exposure')
ss <- function(T, non) T * (non / T) * (1 - non / T)
say('OL+BDC 8/6 %.3f; CO 8/4 %.3f; N-of-1 12/7 %.3f; Hybrid 8/5 %.3f, 8/4 %.3f',
    ss(8, 6), ss(8, 4), ss(12, 7), ss(8, 5), ss(8, 4))
say('identity check (1/p+1/q)^-1 at 5,3: %.3f', 1 / (1 / 5 + 1 / 3))

say('--- I. dichotomizing a normal biomarker at the median')
say('sqrt(2/pi) = %.4f, 2/pi = %.4f, pi/2 = %.4f', sqrt(2 / pi), 2 / pi, pi / 2)
tr <- replicate(4000, {
  n <- 200; b <- rnorm(n); d <- 0.15 * b + rnorm(n)
  x <- as.numeric(b > median(b))
  c(summary(lm(d ~ b))$coef[2, 3], summary(lm(d ~ x))$coef[2, 3])
})
say('mean t continuous %.3f, dichotomized %.3f, ratio %.3f',
    mean(tr[1, ]), mean(tr[2, ]), mean(tr[2, ]) / mean(tr[1, ]))

say('--- K. OLS design effect for a between-participant effect')
deff <- 1 + 7 * 0.3
say('deff %.2f; realized size %.3f', deff, 2 * pnorm(-1.96 / sqrt(deff)))

say('--- L. proportional rule under the compounding recursion')
h <- 1
tsd_run <- c(1, 2, 3)
mu_c <- numeric(4); mu_c[1] <- 10.19
for (j in 1:3) mu_c[j + 1] <- mu_c[j] * 2^(-tsd_run[j] / h)
say('ratios %s; graded g %s', paste(round(mu_c[-1] / mu_c[-4], 4), collapse = ' '),
    paste(round(2^(-tsd_run / h), 4), collapse = ' '))
