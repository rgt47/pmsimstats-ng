set.seed(20261008)
cat('== Gompertz landmarks (m = 10.99, d = 5, r = 0.42) ==\n')
gomp <- function(t, m = 10.99, d = 5, r = 0.42) {
  m * (exp(-d * exp(-r * t)) - exp(-d)) / (1 - exp(-d))
}
cat('inflection ln(d)/r =', log(5) / 0.42, '\n')
cat('t at 50% =', uniroot(\(t) gomp(t) / 10.99 - 0.5, c(0, 20))$root, '\n')
cat('t at 90% =', uniroot(\(t) gomp(t) / 10.99 - 0.9, c(0, 30))$root, '\n')
cat('G(4) =', gomp(4), ' share =', gomp(4) / 10.99, '\n')

cat('\n== Turnover model and PD sensitivity moments ==\n')
t_half_r <- 1
k_out <- log(2) / t_half_r
visits <- 1:10
on <- visits <= 6
f_at <- function(tt) {
  ifelse(tt <= 6, 1 - 2^(-tt / t_half_r),
    (1 - 2^(-6 / t_half_r)) * 2^(-(tt - 6) / t_half_r))
}
ode_f <- function(tt, dt = 1e-4) {
  grid <- seq(0, max(tt), by = dt)
  f <- numeric(length(grid))
  for (j in 2:length(grid)) {
    target <- as.numeric(grid[j - 1] < 6)
    f[j] <- f[j - 1] + dt * k_out * (target - f[j - 1])
  }
  approx(grid, f, tt)$y
}
cat('max |closed form - Euler ODE| =',
  max(abs(f_at(visits) - ode_f(visits))), '\n')
f <- f_at(visits)
a <- ifelse(on, 1 / sum(on), -1 / sum(!on))
F_val <- sum(a * f)
cat('F =', F_val, '\n')
n <- 4e5
omega_bar <- 10
zeta <- 0.3
sigma_v <- 2
sigma_xi <- 3
rho_xi <- 0.5
sigma_bm <- 2
b <- rnorm(n)
B <- 5 + sigma_bm * b
v <- rnorm(n, 0, sigma_v)
omega <- omega_bar * (1 + zeta * b) + v
R <- rho_xi^abs(outer(visits, visits, '-'))
xi <- matrix(rnorm(n * 10), n) %*% chol(R) * sigma_xi
br <- outer(omega, f) + xi
other <- matrix(rnorm(n * 10, 0, 4), n)
y <- 20 - (br + other)
delta <- as.vector(y %*% a)
H <- omega_bar^2 * zeta^2 + sigma_v^2
V0 <- sigma_xi^2 * drop(t(a) %*% R %*% a) + 16 * sum(a^2)
cat('beta sim =', coef(lm(delta ~ B))[2],
  ' formula =', -omega_bar * zeta * F_val / sigma_bm, '\n')
cat('Var(Delta) sim =', var(delta), ' formula =', H * F_val^2 + V0, '\n')
cat('rho sim =', cor(delta, B), ' formula =',
  -omega_bar * zeta * F_val / sqrt(H * F_val^2 + V0), '\n')
cat('Cor(B, BR_6) sim =', cor(B, br[, 6]), ' formula =',
  omega_bar * zeta * f[6] / sqrt(H * f[6]^2 + sigma_xi^2), '\n')
cat('cap omega_bar zeta / sqrt(H) =', omega_bar * zeta / sqrt(H),
  ' attainable at |F| = 1:', omega_bar * zeta / sqrt(H + V0), '\n')
cat('recursion check: mu2/mu_stop =', 2^(-1) * 2^(-2), ' = 2^-3\n')

cat('\n== Threshold moderation projection ==\n')
blin <- integrate(\(x) x * dnorm(x), 0, Inf)$value
cat('beta_lin/Delta =', blin, ' phi(0) =', dnorm(0), '\n')
cat('R2 = phi(0)^2/(1/4) =', dnorm(0)^2 / 0.25, ' 2/pi =', 2 / pi, '\n')
cat('ncp ratio sqrt(R2) =', sqrt(dnorm(0)^2 / 0.25), '\n')
reps <- 4000
m <- 200
tl <- ti <- numeric(reps)
for (r in seq_len(reps)) {
  bb <- rnorm(m)
  yy <- 0.5 * (bb > 0) + rnorm(m, 0, 2)
  tl[r] <- summary(lm(yy ~ bb))$coefficients[2, 3]
  ti[r] <- summary(lm(yy ~ I(bb > 0)))$coefficients[2, 3]
}
cat('mean t ratio linear/indicator (sim) =', mean(tl) / mean(ti), '\n')

cat('\n== Reliability attenuation ==\n')
rb <- 0.8
m <- 400
tt_true <- tt_obs <- numeric(reps)
for (r in seq_len(reps)) {
  bt <- rnorm(m)
  bo <- bt + rnorm(m, 0, sqrt((1 - rb) / rb))
  yy <- 0.25 * bt + rnorm(m)
  tt_true[r] <- summary(lm(yy ~ bt))$coefficients[2, 3]
  tt_obs[r] <- summary(lm(yy ~ bo))$coefficients[2, 3]
}
cat('ncp ratio sim =', mean(tt_obs) / mean(tt_true),
  ' sqrt(R_B) =', sqrt(rb), '\n')
rho0 <- 0.25 / sqrt(1 + 0.25^2)
cat('exact ratio sqrt(R)*sqrt((1-rho^2)/(1-R rho^2)) =',
  sqrt(rb) * sqrt((1 - rho0^2) / (1 - rb * rho0^2)), '\n')
bo_big <- rnorm(1e6)
bt_big <- rb * bo_big + rnorm(1e6, 0, sqrt(rb * (1 - rb)))
cat('raw-scale slope of B on B* =', coef(lm(bt_big ~ bo_big))[2], '\n')

cat('\n== Range restriction ==\n')
mills <- function(z) dnorm(z) / (1 - pnorm(z))
for_r <- c(0.3, 0.5, 0.8)
zz <- 1
s <- rnorm(2e6)
for (rr in for_r) {
  bb <- rr * s + sqrt(1 - rr^2) * rnorm(2e6)
  sel <- s > zz
  cz <- mills(zz)
  cat('r =', rr, ' sim var =', var(bb[sel]), ' formula =',
    1 - rr^2 * cz * (cz - zz), '\n')
}

cat('\n== Regression to the mean (Barnett et al. 2005) ==\n')
sb <- 1
sw <- 0.7
tau <- 1
bl <- rnorm(2e6, 0, sb)
scr <- bl + rnorm(2e6, 0, sw)
later <- bl + rnorm(2e6, 0, sw)
sel <- scr > tau
sd_tot <- sqrt(sb^2 + sw^2)
cat('sim fall =', mean(scr[sel] - later[sel]), ' formula =',
  sw^2 / sd_tot * mills(tau / sd_tot), '\n')

cat('\n== Mixture of exponentials ==\n')
sh <- 0.5
th <- 1
lam <- log(2) / th
gbar <- function(t) {
  integrate(\(z) 2^(-t / (th * exp(sh * z))) * dnorm(z), -Inf, Inf)$value
}
for (tsd in c(0.1, 0.5, 1, 2, 4, 8)) {
  cat('t_sd =', tsd, ' gbar =', round(gbar(tsd), 4), ' median =',
    round(2^(-tsd / th), 4), '\n')
}
cat('initial slope -gbar\'(0) =', (1 - gbar(1e-4)) / 1e-4,
  ' lambda exp(sh^2/2) =', lam * exp(sh^2 / 2), '\n')
rate_at <- function(t, h = 1e-3) -(log(gbar(t + h)) - log(gbar(t - h))) / (2 * h)
cat('local decay rate at 0.5, 2, 8:', rate_at(0.5), rate_at(2), rate_at(8),
  ' (median rate', lam, ')\n')

cat('\n== Unblinding bias ==\n')
chi <- 0.6
pbar <- 3
spb <- 1.2
bb <- rnorm(4e5)
P <- pbar + spb * bb
eta_on <- 0.5 + chi * 0.5
eta_off <- 0.5 - chi * 0.5
y_on <- 20 - eta_on * P + rnorm(4e5)
y_off <- 20 - eta_off * P + rnorm(4e5)
dd <- y_on - y_off
cat('mean contrast =', mean(dd), ' -chi*pbar =', -chi * pbar, '\n')
cat('slope on b =', coef(lm(dd ~ bb))[2], ' -chi*s_pb =', -chi * spb, '\n')

cat('\n== Adherence attenuation ==\n')
abar <- 0.8
ad <- rbeta(4e5, 8, 2)
bb <- rnorm(4e5)
eff <- ad * 10 * (1 + 0.3 * bb)
cat('mean adherence =', mean(ad), ' slope =', coef(lm(eff ~ bb))[2],
  ' abar*3 =', mean(ad) * 3, '\n')
