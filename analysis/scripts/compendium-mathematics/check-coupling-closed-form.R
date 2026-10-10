## Independent check of docs/37 closed-form identities (fragment cg).
repo <- '.'
dd <- file.path(repo, 'analysis/data/quick-sim/carryover-closed-form')

s_tv <- 10; s_pb <- 10; s_br <- 8; s_bm <- 15.36
rho <- 0.8; c1 <- 0.2; cx <- 0.1
eta <- c(1, 1, rep(0.5, 6))
on <- list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1),
           c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1))
tsd <- list(c(0,0,0,0,1,2,0,4), c(0,0,0,0,1,2,6,0),
            c(0,0,0,1,2,3,0,4), c(0,0,0,1,2,3,7,0))
n_p <- c(18, 18, 17, 17); N <- 70; q <- n_p / N
cbm <- 0.25; b0 <- -cbm * s_br / s_bm

## full response covariance (TV, PB, BR; 8 visits each), CS construct
resp_cov <- function() {
  sd <- list(rep(s_tv, 8), s_pb * eta, rep(s_br, 8))
  S <- matrix(0, 24, 24)
  for (f in 1:3) for (g in 1:3) for (t in 1:8) for (s in 1:8) {
    r <- if (f == g) (if (t == s) 1 else rho) else (if (t == s) c1 else cx)
    S[(f-1)*8 + t, (g-1)*8 + s] <- r * sd[[f]][t] * sd[[g]][s]
  }
  S
}
S <- resp_cov()

weights <- function(p, variant) {
  o <- on[[p]]; non <- sum(o); noff <- sum(1 - o)
  dn <- if (variant == 'E9') noff else noff + 1
  ifelse(o == 1, 1 / non, -1 / dn)
}
var_formula <- function(a) {
  m1 <- sum(a); m2 <- sum(a^2); me <- sum(a * eta)
  mee <- sum(a^2 * eta^2); m2e <- sum(a^2 * eta)
  (s_tv^2 + s_br^2) * (rho * m1^2 + (1 - rho) * m2) +
    s_pb^2 * (rho * me^2 + (1 - rho) * mee) +
    2 * s_tv * s_br * (cx * m1^2 + (c1 - cx) * m2) +
    2 * (s_tv + s_br) * s_pb * (cx * m1 * me + (c1 - cx) * m2e)
}
cat('--- Var(Delta): formula vs full matrix ---\n')
for (v in c('E9', 'E9b')) for (p in 1:4) {
  a <- weights(p, v); A <- c(a, a, a)
  cat(v, p, round(var_formula(a), 4), round(drop(t(A) %*% S %*% A), 4),
      '\n')
}

coupling <- function(p, rule, th) {
  o <- on[[p]]; ts <- tsd[[p]]
  g <- ifelse(o == 1, 1, 0)
  if (th > 0) {
    off_exp <- o == 0 & ts > 0
    g[off_exp] <- if (rule == 'step') 1 else 2^(-ts[off_exp] / th)
  }
  g
}
gbar_off <- function(p, rule, th) mean(coupling(p, rule, th)[on[[p]] == 0])

## slope by direct covariance vs Proposition 2
cat('--- slope: Cov(Delta,B)/var(B) vs Prop 2 (graded, th = 1) ---\n')
for (v in c('E9', 'E9b')) for (p in 1:4) {
  a <- weights(p, v); g <- coupling(p, 'graded', 1)
  direct <- -cbm * s_bm * s_br * sum(a * g) / s_bm^2
  noff <- sum(1 - on[[p]]); gb <- gbar_off(p, 'graded', 1)
  prop <- if (v == 'E9') b0 * (1 - gb) else b0 * (1 - noff / (noff + 1) * gb)
  cat(v, p, signif(direct, 8), signif(prop, 8), '\n')
}
pp <- read.csv(file.path(dd, 'e9-per-path.csv'))
sub <- pp[pp$design == 'Hybrid' & pp$construct == 'CS' &
          pp$rule == 'graded' & abs(pp$t_half - 1) < 1e-9, ]
cat('stored per-path graded th=1:\n'); print(sub[, c('variant', 'path',
  'gamma', 'var_D', 'tau2', 'gbar_off')])

## conditional identity: covariance binding = mean moderation with h = g
## plus rank-one smaller residual covariance
cat('--- conditional identity, path 1, graded, th = 1 ---\n')
g <- coupling(1, 'graded', 1)
cvec <- c(rep(0, 16), cbm * g * s_br * s_bm)
Sfull <- rbind(c(s_bm^2, cvec), cbind(cvec, S))
Sc <- S - outer(cvec, cvec) / s_bm^2
br_idx <- 17:24
rank1 <- S[br_idx, br_idx] - cbm^2 * s_br^2 * outer(g, g)
cat('max |Cond cov BR - (Sigma_BR - c^2 s^2 g g^T)|:',
    max(abs(Sc[br_idx, br_idx] - rank1)), '\n')
cat('min eigen full joint:', min(eigen(Sfull)$values), '\n')
a <- weights(1, 'E9'); A <- c(a, a, a)
tau_cov <- drop(t(A) %*% Sc %*% A)
beta_p <- b0 * (1 - gbar_off(1, 'graded', 1))
cat('tau2 cov:', tau_cov, ' Var - beta^2 s^2:',
    var_formula(a) - beta_p^2 * s_bm^2, ' tau2 decayed MM:', var_formula(a),
    '\n')

## pooled power (docs/37 Section 4.4)
pooled_power <- function(beta, tau2, mu, alpha = 0.05, K = 400) {
  bb <- sum(q * beta); mb <- sum(q * mu)
  s2 <- sum(q * (tau2 + (mu - mb)^2 + s_bm^2 * (beta - bb)^2))
  ss <- s_bm^2 * qchisq(((1:K) - 0.5) / K, N - 1)
  ncp <- bb * sqrt(ss) / sqrt(s2); cr <- qt(1 - alpha / 2, N - 2)
  mean(pt(-cr, N - 2, ncp) + 1 - pt(cr, N - 2, ncp))
}
mu_lookup <- function(rule, th, v) {
  s <- pp[pp$design == 'Hybrid' & pp$construct == 'CS' & pp$rule == rule &
          abs(pp$t_half - th) < 1e-9 & pp$variant == v, ]
  if (nrow(s) == 4) s$mu[order(s$path)] else rep(0, 4)
}
cat('rules in per-path csv:', unique(pp$rule), '\n')
cat('t_half in per-path csv:', sort(unique(pp$t_half)), '\n')
cat('--- pooled power by rule (mu from csv when present) ---\n')
for (rule in c('step', 'graded', 'meanmod', 'decayed')) {
  out <- c()
  for (th in c(0, 0.1, 0.2, 0.5, 1, 2)) {
    for (v in c('E9', 'E9b')) {
      beta <- tau2 <- numeric(4)
      for (p in 1:4) {
        a <- weights(p, v); noff <- sum(1 - on[[p]])
        r <- if (rule %in% c('step', 'graded')) rule else
          if (rule == 'decayed') 'graded' else 'none'
        gb <- if (r == 'none') 0 else gbar_off(p, r, th)
        k <- if (v == 'E9') 1 else noff / (noff + 1)
        beta[p] <- b0 * (1 - k * gb)
        tau2[p] <- if (rule %in% c('step', 'graded'))
          var_formula(a) - beta[p]^2 * s_bm^2 else var_formula(a)
      }
      csvrule <- switch(rule, meanmod = 'mean',
                        decayed = 'mean_decayed', rule)
      mu <- mu_lookup(csvrule, th, v)
      out <- c(out, sprintf('%s:%.3f', v, pooled_power(beta, tau2, mu)))
    }
  }
  cat(rule, out, '\n')
}

## coupling gap, patient weighted
cat('--- graded gap (E9, patient weighted) ---\n')
for (th in c(0.2, 0.5, 1, 2)) {
  gb <- sapply(1:4, function(p) gbar_off(p, 'graded', th))
  cat(th, round(1 - sum(q * gb), 4), '\n')
}
cat('bound check th=0.1, max gbar:',
    max(sapply(1:4, function(p) gbar_off(p, 'graded', 0.1))),
    ' 2^-10 =', 2^-10, '\n')

## derivative formula vs finite difference, path 3, th = 1
p <- 3; th <- 1; h <- 1e-6
fd <- (gbar_off(p, 'graded', th + h) - gbar_off(p, 'graded', th - h)) / (2*h)
o <- on[[p]]; ts <- tsd[[p]][o == 0]
an <- log(2) / (length(ts) * th^2) * sum(ts * 2^(-ts / th))
cat('dgbar/dth fd:', fd, ' analytic:', an, '\n')

## power vs gap (uniform off coupling 1 - gap), E9
cat('--- power vs gap ---\n')
for (gap in c(0, 0.25, 0.5, 0.75, 1)) {
  r <- c()
  for (v in c('E9', 'E9b')) {
    beta <- tau2 <- numeric(4)
    for (p in 1:4) {
      a <- weights(p, v); noff <- sum(1 - on[[p]])
      k <- if (v == 'E9') 1 else noff / (noff + 1)
      beta[p] <- b0 * (1 - k * (1 - gap))
      tau2[p] <- var_formula(a) - beta[p]^2 * s_bm^2
    }
    r <- c(r, sprintf('%s:%.3f', v, pooled_power(beta, tau2, rep(0, 4))))
  }
  cat(gap, r, '\n')
}
gv <- read.csv(file.path(dd, 'e9-power-vs-gap.csv'))
print(gv[abs(gv$gap - round(gv$gap * 4) / 4) < 1e-9, ])

## sample-size scaling, single path, exact: n - 1 = L^2 (1/rho^2 - 1)
a <- weights(1, 'E9'); vD <- var_formula(a)
for (gap in c(1, 0.5)) {
  rr <- b0 * gap * s_bm / sqrt(vD)
  cat('gap', gap, 'rho_DB', round(rr, 4), ' (1/rho^2 - 1):',
      round(1 / rr^2 - 1, 2), '\n')
}

## compounding vs anchored mean, docs/40 example (path 1, th = 1)
m_last <- 10.19; ts1 <- c(1, 2)
anch <- m_last * 2^(-ts1 / 1)
comp <- m_last * 2^(-cumsum(ts1) / 1)
cat('anchored:', anch, ' compounding:', comp, '\n')
