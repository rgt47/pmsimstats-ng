repo <- '.'

g <- function(t, m = 10.99, d = 5, r = 0.42) {
  m * (exp(-d * exp(-r * t)) - exp(-d)) / (1 - exp(-d))
}
cat('inflection ln(d)/r:', log(5) / 0.42, '\n')
cat('half-max t:', uniroot(function(t) g(t) - 10.99 / 2, c(0, 20))$root, '\n')
cat('90% t:', uniroot(function(t) g(t) - 0.9 * 10.99, c(0, 30))$root, '\n')
cat('84% at week 8:', g(8) / 10.99, '\n')

z <- qnorm(0.975)
cat('d alpha/d kappa at 1:', -2 * z * dnorm(z), '\n')
cat('elasticity:', -2 * z * dnorm(z) / 0.05, '\n')
cat('alpha at kappa 1.05, 1.10:', 2 * pnorm(-1.05 * z), 2 * pnorm(-1.10 * z), '\n')
cat('n for MCSE 0.02 at 0.5:', 0.25 / 0.02^2, '\n')

b <- readRDS(file.path(repo, 'analysis/data/derived_data/component-decomposition',
  'study-b-balanced/study-b-balanced-summary.rds'))
cat('study B balanced N values:', sort(unique(b$N)), '\n')

rm_sim <- function(n = 30, T = 3, sig_s = 2, sig_e = 1, reps = 2000) {
  out <- replicate(reps, {
    s <- rnorm(2 * n, 0, sig_s)
    y1 <- matrix(rnorm(2 * n * T, 0, sig_e), 2 * n) + s
    y2 <- matrix(rnorm(2 * n * T, 0, sig_e), 2 * n) + s
    d <- rowMeans(y1) - rowMeans(y2)
    grp <- rep(1:2, each = n)
    L <- mean(d[grp == 1]) - mean(d[grp == 2])
    sp2 <- (var(d[grp == 1]) + var(d[grp == 2])) / 2
    c(L = L, sd_check = sp2, t = L / sqrt(sp2 * 2 / n))
  })
  out
}
set.seed(1)
o <- rm_sim()
cat('Var(L) sim:', var(o['L', ]), ' theory 4 sig_e^2/(nT) = ', 4 / (30 * 3), '\n')
cat('pooled var of d: mean', mean(o['sd_check', ]), ' theory 2 sig_e^2/T =',
  2 / 3, '\n')
cat('rejection rate of pooled t, df 2(n-1):',
  mean(abs(o['t', ]) > qt(0.975, 58)), '\n')
