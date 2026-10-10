## Copied 2026-10-08 from the review scratch script inv_b_check.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-gompertz-rmanova-cycles.R
## Inventory part b: checks for docs/24-28 results.
## Self-contained (base R + nlme); run as Rscript analysis/scripts/compendium-mathematics/check-gompertz-rmanova-cycles.R.

options(digits = 6)
hdr <- function(s) cat('\n==== ', s, ' ====\n', sep = '')

## ---------------------------------------------------------------
hdr('1. Gompertz landmarks and conditioning (docs/25)')
mg <- function(t, m, d, r) {
  y <- m * exp(-d * exp(-r * t))
  v <- m * exp(-d)
  (y - v) * (m / (m - v))
}
cat('inflection BR ln5/0.42 =', log(5) / 0.42, '\n')
cat('inflection TV/PB ln5/0.35 =', log(5) / 0.35, '\n')
cat('BR at week 12 / maxr =', mg(12, 10.99, 5, 0.42) / 10.99, '\n')
frac_inf <- (exp(-1) - exp(-5)) / (1 - exp(-5))
cat('shifted fraction at inflection, d=5 =', frac_inf, '\n')
tt <- seq(0, 20, by = 1e-4)
yy <- mg(tt, 1, 5, 0.42)
sl <- diff(yy) / diff(tt)
i_max <- which.max(sl)
cat('numeric inflection =', tt[i_max], ' value =', yy[i_max], '\n')
for (d in c(2, 1, 0.5)) {
  ti <- log(d) / 0.42
  f <- (exp(-1) - exp(-d)) / (1 - exp(-d))
  cat(sprintf('d=%.2f t*=%.3f fraction=%.4f\n', d, ti, f))
}
cat('rescale factor 1/(1-e^-d) at d=5:', 1 / (1 - exp(-5)),
    ' d=0.01:', 1 / (1 - exp(-0.01)), '\n')
for (d in c(1e-4, 1e-8, 1e-12, 1e-14)) {
  num <- mg(5, 1, d, 0.42)
  ex <- (exp(-d * exp(-0.42 * 5)) - exp(-d)) / (-expm1(-d))
  stable <- -expm1(-d * exp(-0.42 * 5)) * 0 +
    (expm1(-d * exp(-0.42 * 5)) - expm1(-d)) / (-expm1(-d))
  cat(sprintf('d=%g naive=%.12f stable=%.12f limit=%.12f relerr=%.2e\n',
              d, num, stable, 1 - exp(-0.42 * 5),
              abs(num - stable) / stable))
}
cat('disp = 0 gives', mg(5, 1, 0, 0.42), '\n')

## ---------------------------------------------------------------
hdr('2. Hill / Emax inflection and slope at zero (docs/25)')
hill <- function(t, E, c, h) E * t^h / (c^h + t^h)
for (h in c(1.5, 2, 3)) {
  cc <- 4
  tt <- seq(1e-6, 20, by = 1e-4)
  sl <- diff(hill(tt, 1, cc, h)) / diff(tt)
  cat(sprintf('h=%.1f numeric t*=%.4f formula=%.4f  y(t*)=%.4f (h-1)/(2h)=%.4f\n',
              h, tt[which.max(sl)], cc * ((h - 1) / (h + 1))^(1 / h),
              hill(tt[which.max(sl)], 1, cc, h), (h - 1) / (2 * h)))
}
for (h in c(0.5, 1, 2)) {
  e <- 1e-8
  cat(sprintf('h=%.1f slope near 0 = %.4g\n', h, hill(e, 1, 4, h) / e))
}
cat('Hill as logistic in ln t: y(c e^s) = 1/(1+e^{-hs}) check:',
    hill(4 * exp(0.3), 1, 4, 2), 1 / (1 + exp(-2 * 0.3)), '\n')

## ---------------------------------------------------------------
hdr('3. RM-ANOVA sums of squares (docs/26)')
set.seed(20261008)
sim_split <- function(a, b, n, Tt, su = 1, se = 1, eff = 0) {
  d <- expand.grid(time = seq_len(Tt), phase = seq_len(b),
                   id = seq_len(n), g = seq_len(a))
  d$subj <- interaction(d$g, d$id)
  u <- rnorm(a * n, 0, su)
  d$y <- 10 + 0.5 * d$g + 0.3 * d$phase + 0.1 * d$time +
    eff * (d$g == 1) * (d$phase == 1) +
    u[as.integer(d$subj)] + rnorm(nrow(d), 0, se)
  d$g <- factor(d$g); d$phase <- factor(d$phase)
  d$time <- factor(d$time)
  d
}
ss_manual <- function(d, a, b, n, Tt) {
  gm <- mean(d$y)
  mg_ <- tapply(d$y, d$g, mean)
  mj <- tapply(d$y, d$phase, mean)
  mgj <- tapply(d$y, list(d$g, d$phase), mean)
  mgi <- tapply(d$y, d$subj, mean)
  mgji <- tapply(d$y, list(d$subj, d$phase), mean)
  g_of_subj <- tapply(as.integer(d$g), d$subj, function(z) z[1])
  SSA <- b * n * Tt * sum((mg_ - gm)^2)
  SSSA <- b * Tt * sum((mgi - mg_[g_of_subj])^2)
  SSB <- a * n * Tt * sum((mj - gm)^2)
  int <- outer(mg_, mj, '+')
  SSAB <- n * Tt * sum((mgj - int + gm)^2)
  res <- mgji - matrix(mgi, ncol = b, nrow = length(mgi)) -
    mgj[g_of_subj, ] + matrix(mg_[g_of_subj], ncol = b,
                              nrow = length(mgi))
  SSBSA <- Tt * sum(res^2)
  c(SSA = SSA, SSSA = SSSA, SSB = SSB, SSAB = SSAB, SSBSA = SSBSA)
}
for (cfg in list(c(2, 2, 10, 3), c(3, 3, 8, 4))) {
  a <- cfg[1]; b <- cfg[2]; n <- cfg[3]; Tt <- cfg[4]
  d <- sim_split(a, b, n, Tt)
  s <- ss_manual(d, a, b, n, Tt)
  fit <- aov(y ~ g * phase * time + Error(subj / (phase * time)),
             data = d)
  sm <- summary(fit)
  cat(sprintf('a=%d b=%d n=%d T=%d\n', a, b, n, Tt))
  print(round(s, 6))
  print(sm[['Error: subj']])
  print(sm[['Error: subj:phase']])
  print(sm[['Error: subj:phase:time']])
  FAB <- (s['SSAB'] / ((a - 1) * (b - 1))) /
    (s['SSBSA'] / (a * (n - 1) * (b - 1)))
  cat('manual F_AB =', FAB, '\n')
  if (a == 2 && b == 2) {
    mgj <- tapply(d$y, list(d$g, d$phase), mean)
    L <- (mgj[1, 1] - mgj[1, 2]) - (mgj[2, 1] - mgj[2, 2])
    xi_phi_11 <- mgj[1, 1] - mean(mgj[1, ]) - mean(mgj[, 1]) +
      mean(mgj)
    cat('xi.phi_11 =', xi_phi_11, ' L/4 =', L / 4, '\n')
    cat('SS_AB =', s['SSAB'], ' nT L^2/4 =', n * Tt * L^2 / 4, '\n')
  }
}

## ---------------------------------------------------------------
hdr('4. Box epsilon (docs/26)')
box_eps <- function(S) {
  m <- nrow(S)
  C <- qr.Q(qr(cbind(1, diag(m))[, 1:m]))[, -1]
  M <- t(C) %*% S %*% C
  sum(diag(M))^2 / ((m - 1) * sum(diag(M %*% M)))
}
for (m in c(2, 3, 4, 8)) {
  for (rho in c(0.5, 0.8)) {
    S <- rho^abs(outer(1:m, 1:m, '-'))
    cat(sprintf('AR(1) m=%d rho=%.1f eps=%.4f  bound=%.4f\n', m, rho,
                box_eps(S), 1 / (m - 1)))
  }
  Scs <- 0.5 + 0.5 * diag(m)
  cat(sprintf('CS m=%d eps=%.4f\n', m, box_eps(Scs)))
}
m <- 4
v <- c(3, -1, -1, -1)
cat('rank-one contrast covariance eps =',
    box_eps(diag(m) * 1e-12 + outer(v, v)), ' 1/(m-1) =', 1 / (m - 1),
    '\n')

## ---------------------------------------------------------------
hdr('5. RM-ANOVA sample size by iteration (docs/26 eq. 19)')
pow_n <- function(n, Tt, es, alpha = 0.05) {
  ncp <- n * Tt * es^2 / 4
  df2 <- 2 * (n - 1)
  1 - pf(qf(1 - alpha, 1, df2), 1, df2, ncp = ncp)
}
for (es in c(1, 0.5)) {
  for (Tt in c(1, 4)) {
    n <- 2
    while (pow_n(n, Tt, es) < 0.8) n <- n + 1
    cat(sprintf('L1/sigma=%.1f T=%d: n per stratum=%d power=%.4f Lambda=%.3f\n',
                es, Tt, n, pow_n(n, Tt, es), n * Tt * es^2 / 4))
  }
}

## ---------------------------------------------------------------
hdr('6. Phase-augmented model: aliasing (docs/24, docs/28)')
build_long <- function(tp, eta, ondrug, t12, nper = 30) {
  out <- list()
  dt <- c(tp[1], diff(tp))
  for (p in seq_along(ondrug)) {
    od <- ondrug[[p]]
    tsd <- dt * (od != 1)
    for (j in 2:length(tp)) if (od[j] != 1) tsd[j] <- tsd[j] + tsd[j - 1]
    tsd <- (cumsum(od) > 0) * tsd
    Dbc <- ifelse(od == 1, 1, ifelse(tsd > 0, 0.5^(tsd / t12), 0))
    for (i in seq_len(nper)) {
      out[[length(out) + 1]] <- data.frame(
        ptID = paste(p, i), t = tp, eta = eta, D = od, Dbc = Dbc)
    }
  }
  d <- do.call(rbind, out)
  bm <- rnorm(length(unique(d$ptID)))
  names(bm) <- unique(d$ptID)
  d$bm <- bm[d$ptID]
  d
}
tp_full <- seq(1, 20, by = 1.25)
eta_full <- c(rep(1, 8), rep(0.5, 8))
od_full <- list(AB = c(rep(1, 8), rep(0, 4), 1, 1, 0, 0),
                BA = c(rep(1, 8), rep(0, 4), 0, 0, 1, 1))
d <- build_long(tp_full, eta_full, od_full, 1)
d$phase <- factor(ifelse(d$eta >= 0.99, 'OL', 'BDC'),
                  levels = c('OL', 'BDC'))
f_full <- ~ bm + t + Dbc + phase + Dbc:phase + bm:Dbc + bm:Dbc:phase
f_red <- ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase
X <- model.matrix(f_full, d)
cat('full formula: ncol =', ncol(X), ' rank =', qr(X)$rank, '\n')
X2 <- model.matrix(f_red, d)
cat('fallback formula: ncol =', ncol(X2), ' rank =', qr(X2)$rank, '\n')
d$y <- rnorm(nrow(d)) + rnorm(length(unique(d$ptID)))[as.integer(factor(d$ptID))]
cf <- coef(lm(update(f_full, y ~ .), data = d))
cat('lm NA coefficients:', names(cf)[is.na(cf)], '\n')
r <- tryCatch(nlme::lme(update(f_full, y ~ .), random = ~ 1 | ptID,
                        correlation = nlme::corCAR1(form = ~ t | ptID),
                        data = d),
              error = function(e) conditionMessage(e))
cat('nlme::lme on full formula:', if (is.character(r)) r else 'fitted',
    '\n')
lhs <- X[, 'Dbc']
rhs <- X[, '(Intercept)'] - X[, 'phaseBDC'] + X[, 'Dbc:phaseBDC']
cat('max |Dbc - (1 - phaseBDC + Dbc:phaseBDC)| =', max(abs(lhs - rhs)),
    '\n')

## docs/24 three-phase Hybrid (8 visits, 4 paths)
tp8 <- c(4, 8, 9, 10, 11, 12, 16, 20)
eta8 <- c(1, 1, rep(0.5, 6))
od8 <- list(c(1, 1, 1, 1, 0, 0, 1, 0), c(1, 1, 1, 1, 0, 0, 0, 1),
            c(1, 1, 1, 0, 0, 0, 1, 0), c(1, 1, 1, 0, 0, 0, 0, 1))
d8 <- build_long(tp8, eta8, od8, 0.5)
d8$phase <- factor(ifelse(d8$t <= 8, 'OL', ifelse(d8$t <= 12, 'BDC',
                                                  'CO')),
                   levels = c('OL', 'BDC', 'CO'))
X8 <- model.matrix(f_full, d8)
X8b <- model.matrix(~ bm + t + Dbc + bm:Dbc, d8)
cat('3-phase Hybrid: added columns =', ncol(X8) - ncol(X8b),
    ' rank gain =', qr(X8)$rank - qr(X8b)$rank, '\n')

## ---------------------------------------------------------------
hdr('7. Corr(D_bc, eta) by design (paper 06, docs/28)')
cor_de <- function(tp, eta, ondrug, t12) {
  d <- build_long(tp, eta, ondrug, t12, nper = 1)
  x <- if (is.na(t12)) d$D else d$Dbc
  cor(x, d$eta)
}
cat('full hybrid (2 paths, 16 visits): binary',
    cor_de(tp_full, eta_full, od_full, NA), ' t1/2=1',
    cor_de(tp_full, eta_full, od_full, 1), '\n')
cat('8-visit Hybrid (4 paths): binary', cor_de(tp8, eta8, od8, NA),
    ' t1/2=0.5', cor_de(tp8, eta8, od8, 0.5), ' t1/2=1',
    cor_de(tp8, eta8, od8, 1), '\n')
tp_bp <- seq(1.5, 22.5, by = 1.5)
eta_bp <- c(rep(1, 3), rep(0, 6), rep(1, 6))
od_bp <- list(c(rep(1, 9), rep(0, 6)))
cat('balanced placebo: binary', cor_de(tp_bp, eta_bp, od_bp, NA),
    ' (exact -2/3 =', -2 / 3, ') t1/2=1',
    cor_de(tp_bp, eta_bp, od_bp, 1), '\n')
cor_de_bl <- function(tp, eta, ondrug, t12) {
  d <- build_long(tp, eta, ondrug, t12, nper = 1)
  x <- c(if (is.na(t12)) d$D else d$Dbc, rep(0, length(ondrug)))
  e <- c(d$eta, rep(0, length(ondrug)))
  cor(x, e)
}
cat('with baseline row (D=0, eta=0) per path:\n')
cat('  full hybrid binary', cor_de_bl(tp_full, eta_full, od_full, NA),
    ' t1/2=1', cor_de_bl(tp_full, eta_full, od_full, 1), '\n')
cat('  8-visit Hybrid binary', cor_de_bl(tp8, eta8, od8, NA),
    ' t1/2=1', cor_de_bl(tp8, eta8, od8, 1), '\n')
cat('  balanced binary', cor_de_bl(tp_bp, eta_bp, od_bp, NA),
    ' t1/2=1', cor_de_bl(tp_bp, eta_bp, od_bp, 1), '\n')
## closed form for a single path with indicator vectors
n1 <- 3; n2 <- 6; n3 <- 6
nn <- n1 + n2 + n3
pD <- (n1 + n2) / nn; pE <- (n1 + n3) / nn; pDE <- n1 / nn
cat('closed form balanced:', (pDE - pD * pE) /
      sqrt(pD * (1 - pD) * pE * (1 - pE)), '\n')

## ---------------------------------------------------------------
hdr('8. Cycle-period information (docs/27)')
cyc <- function(Tol, k, Ton, Toff) {
  D <- c(rep(1, Tol), rep(c(rep(1, Ton), rep(0, Toff)), k))
  eta <- c(rep(1, Tol), rep(0.5, k * (Ton + Toff)))
  t <- seq_along(D)
  list(D = D, eta = eta, t = t)
}
## Information on beta_bm:D per unit sum((b - bbar)^2), sigma^2 = 1,
## from the full design matrix of N participants on one schedule.
## Random intercept with variance su2 handled by GLS.
set.seed(7)
bN <- rnorm(40); bN <- (bN - mean(bN)) / sd(bN)
info <- function(z, adj = c('none', 't', 'eta', 't_eta', 't_bmeta'),
                 su2 = 1) {
  adj <- match.arg(adj)
  n <- length(z$D)
  rows <- lapply(seq_along(bN), function(i) {
    b <- bN[i]
    base <- cbind(1, b, z$D, b * z$D)
    extra <- switch(adj, none = NULL, t = cbind(z$t),
                    eta = cbind(z$eta), t_eta = cbind(z$t, z$eta),
                    t_bmeta = cbind(z$t, z$eta, b * z$eta))
    cbind(base, extra)
  })
  Vi <- solve(diag(n) + su2 * matrix(1, n, n))
  XtVX <- Reduce('+', lapply(rows, function(X) t(X) %*% Vi %*% X))
  1 / solve(XtVX)[4, 4] / sum((bN - mean(bN))^2)
}
lead <- function(k, Ton, Toff) k * Ton * Toff / (Ton + Toff)
olf <- function(Tol, k, Ton, Toff) {
  non <- Tol + k * Ton; noff <- k * Toff
  non * noff / (non + noff)
}
cat('k Ton Toff | lead | noOL | OL8 | OL8 formula | OL8+t | OL8+t+eta |',
    'OL8+t+eta+bm:eta\n')
for (k in 1:4) for (Ton in c(2, 4, 6)) for (Toff in c(2, 4, 6)) {
  if (k * (Ton + Toff) > 24) next
  z0 <- cyc(0, k, Ton, Toff); z8 <- cyc(8, k, Ton, Toff)
  cat(sprintf('%d %d %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f\n',
              k, Ton, Toff, lead(k, Ton, Toff), info(z0), info(z8),
              olf(8, k, Ton, Toff), info(z8, 't'), info(z8, 't_eta'),
              info(z8, 't_bmeta')))
}
cat('fixed S=8, k=2: Ton | OL8 (+t) | OL8+t+eta+bm:eta\n')
for (Ton in 2:6) {
  z8 <- cyc(8, 2, Ton, 8 - Ton)
  cat(sprintf('%d | %.4f | %.4f\n', Ton, info(z8, 't'),
              info(z8, 't_bmeta')))
}
cat('optimum formula S/2 - Tol/(2k) =', 8 / 2 - 8 / (2 * 2), '\n')
cat('independence (su2=0) check, k=2 Ton=Toff=4, OL8+t:',
    info(cyc(8, 2, 4, 4), 't', su2 = 0), '\n')
cnt <- 0
for (k in 1:4) for (Ton in c(2, 4, 6)) for (Toff in c(2, 4, 6))
  if (k * (Ton + Toff) <= 24) cnt <- cnt + 1
cat('feasible cycle structures =', cnt, ' cells =', cnt * 12,
    ' with 5 c_bm values =', cnt * 2 * 3 * 5, '\n')

## ---------------------------------------------------------------
hdr('9. Worked example of docs/24 against plim identity')
A <- c(BR = 4, PB = 6, TV = 1); B <- c(BR = 8, PB = 1, TV = 1)
for (p in list(A, B)) cat('one-component contrast =',
                          p['BR'] + 0.5 * p['PB'], '\n')
cat('mean =', mean(c(A['BR'] + 0.5 * A['PB'], B['BR'] + 0.5 * B['PB'])),
    '\n')
