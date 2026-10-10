## Copied 2026-10-08 from the review scratch script inv_e_check.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-closed-form-extensions.R
## Inventory part e: independent checks of the results added from
## docs/32, 34, 35, 36, 37, 38. Run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-closed-form-extensions.R
## Covariances and ceilings are built here from the stated definitions.
## Only the path means mu_p (which enter the path-mean term of the pooled
## E9 power) are taken from the source's own mean builder
## (analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R).
suppressMessages(library(MASS))
f3 <- function(x) ifelse(is.na(x), 'NA', formatC(as.numeric(x), format = 'f', digits = 3))
hdr <- function(s) cat('\n=====', s, '=====\n')

## ---------------- designs ----------------
mk_path <- function(w, on, e) {
  gaps <- c(w[1], diff(w)); n <- length(on)
  tsd <- numeric(n); ever <- FALSE
  for (t in seq_len(n)) {
    if (on[t] == 1) { tsd[t] <- 0; ever <- TRUE }
    else tsd[t] <- if (ever) (if (t > 1) tsd[t - 1] else 0) + gaps[t] else 0
  }
  list(w = w, on = on == 1, tsd = tsd, e = e)
}
W_H <- c(4, 8, 9, 10, 11, 12, 16, 20)
W_OB <- c(4, 8, 12, 16, 17, 18, 19, 20)
W_CO <- cumsum(rep(2.5, 8))
E_H <- c(1, 1, rep(0.5, 6)); E_OB <- c(1, 1, 1, 1, rep(0.5, 4))
E_CO <- rep(0.5, 8)
D <- list(
  Hybrid = list(paths = lapply(list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1),
                                    c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1)),
                               mk_path, w = W_H, e = E_H),
                n_p = c(18, 18, 17, 17)),
  `OL+BDC` = list(paths = lapply(list(c(1,1,1,1,1,1,0,0), c(1,1,1,1,1,0,0,0)),
                                 mk_path, w = W_OB, e = E_OB),
                  n_p = c(35, 35)),
  CO = list(paths = lapply(list(c(1,1,1,1,0,0,0,0), c(0,0,0,0,1,1,1,1)),
                           mk_path, w = W_CO, e = E_CO),
            n_p = c(35, 35)))

## ---------------- response correlation block ----------------
## within: 'CS' or 'AR'; cross (t != s): 'const' c_x, 'decay' c_x rho^gap,
## 'sep' c1 rho^gap (with AR within this is K (x) A)
build_M <- function(w, rho, c1, cx, within, cross) {
  n <- length(w); g <- abs(outer(w, w, '-'))
  A <- if (within == 'AR') rho^g else { m <- matrix(rho, n, n); diag(m) <- 1; m }
  C <- switch(cross,
    const = { m <- matrix(cx, n, n); diag(m) <- c1; m },
    decay = { m <- cx * rho^g; diag(m) <- c1; m },
    sep = c1 * A)
  rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
}
STRUCT <- list(CS = c('CS', 'const'), covar = c('AR', 'decay'),
               sep = c('AR', 'sep'), ARconst = c('AR', 'const'),
               CSdecay = c('CS', 'decay'))
Mof <- function(p, s, rho = 0.8, c1 = 0.2, cx = 0.1)
  build_M(p$w, rho, c1, cx, STRUCT[[s]][1], STRUCT[[s]][2])

decay_fn <- function(tsd, th, k = 1) {
  if (k == 1) 0.5^(tsd / th) else exp(-((log(2))^(1 / k) / th * tsd)^k)
}
coupling <- function(p, rule, th, k = 1) {
  g <- as.numeric(p$on)
  if (th > 0) {
    oe <- !p$on & p$tsd > 0
    g[oe] <- if (rule == 'step') as.numeric(decay_fn(p$tsd[oe], th, k) > 0)
             else decay_fn(p$tsd[oe], th, k)
  }
  g
}
ceil_path <- function(M, g) {
  n <- length(g); u <- c(rep(0, 2 * n), g)
  1 / sqrt(drop(t(u) %*% solve(M, u)))
}
ceil_design <- function(d, s, rule, th, rho = 0.8, c1 = 0.2, cx = 0.1, k = 1) {
  M1 <- lapply(d$paths, Mof, s = s, rho = rho, c1 = c1, cx = cx)
  if (any(sapply(M1, function(M) min(eigen(M, TRUE, TRUE)$values)) <= 0))
    return(NA)
  min(mapply(function(p, M) ceil_path(M, coupling(p, rule, th, k)),
             d$paths, M1))
}

## ---------------- E9 moments ----------------
s_tv <- 10; s_pb <- 10; s_br <- 8; s_bm <- 15.36159
weights <- function(p, v) {
  non <- sum(p$on); noff <- sum(!p$on) + (v == 'E9b')
  ifelse(p$on, 1 / non, -1 / noff)
}
resp_cov <- function(p, s, rho = 0.8) {
  M <- Mof(p, s, rho); sd <- c(rep(s_tv, 8), s_pb * p$e, rep(s_br, 8))
  outer(sd, sd) * M
}
var_D <- function(p, s, v, rho = 0.8) {
  a <- weights(p, v); A <- c(a, a, a)
  drop(t(A) %*% resp_cov(p, s, rho) %*% A)
}

## path means from the source's builder (means only)
keep <- c('RHO', 'C1', 'CX', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd', 'gomp',
          'designs', 'path_info', 'build')
dgp <- new.env()
f04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
for (e in parse(f04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
mu_path <- function(dn, j, th, v) {
  d <- dgp$designs[[dn]]
  pi <- dgp$path_info(d$t_wk, d$e, d$paths[[j]])
  b <- dgp$build(pi, 'A', 0, th); n <- 8
  muS <- b$mu[3:(2 + n)] + b$mu[(3 + n):(2 + 2 * n)] + b$mu[(3 + 2 * n):(2 + 3 * n)]
  -sum(weights(D[[dn]]$paths[[j]], v) * muS)
}
## mean structure must match: compare the source's tsd with ours
for (dn in names(D)) for (j in seq_along(D[[dn]]$paths)) {
  d <- dgp$designs[[dn]]
  pi <- dgp$path_info(d$t_wk, d$e, d$paths[[j]])
  stopifnot(max(abs(pi$tsd - D[[dn]]$paths[[j]]$tsd)) < 1e-12,
            max(abs(pi$e - D[[dn]]$paths[[j]]$e)) < 1e-12)
}

N <- 70; KQ <- 400
moments <- function(dn, s, rule, th, v, cbm = 0.25, g_override = NULL,
                    drop_mu = FALSE) {
  d <- D[[dn]]; P <- length(d$paths)
  out <- t(sapply(seq_len(P), function(j) {
    p <- d$paths[[j]]; a <- weights(p, v); vd <- var_D(p, s, v)
    if (rule %in% c('step', 'graded') || !is.null(g_override)) {
      g <- if (is.null(g_override)) coupling(p, rule, th) else g_override[[j]]
      beta <- -cbm * s_br / s_bm * sum(a * g); tau2 <- vd - beta^2 * s_bm^2
    } else {
      h <- if (rule == 'mean') as.numeric(p$on) else coupling(p, 'graded', th)
      beta <- -cbm * s_br / s_bm * sum(a * h); tau2 <- vd
    }
    mu <- if (drop_mu) 0 else mu_path(dn, j, th, v)
    c(beta = beta, tau2 = tau2, mu = mu, vd = vd)
  }))
  as.data.frame(out)
}
pooled <- function(m, n_p, path_term = TRUE) {
  q <- n_p / sum(n_p); bb <- sum(q * m$beta); mb <- sum(q * m$mu)
  pm <- if (path_term) sum(q * (m$mu - mb)^2) else 0
  s2 <- sum(q * (m$tau2 + s_bm^2 * (m$beta - bb)^2)) + pm
  ss <- s_bm^2 * qchisq(((1:KQ) - 0.5) / KQ, N - 1)
  ncp <- bb * sqrt(ss) / sqrt(s2); cr <- qt(0.975, N - 2)
  V <- (sum(q * m$tau2) + pm) / ((N - 3) * s_bm^2) +
    2 * sum(q * (m$beta - bb)^2) / N
  list(power = mean(pt(-cr, N - 2, ncp) + 1 - pt(cr, N - 2, ncp)),
       mean = bb, sd = sqrt(V), s2 = s2, pm = pm)
}
pw <- function(dn, s, rule, th, v = 'E9', cbm = 0.25, path_term = TRUE)
  pooled(moments(dn, s, rule, th, v, cbm), D[[dn]]$n_p, path_term)

## ======================================================================
hdr('docs/37 5: Hybrid CS mean slope and E9b mean-moderation power')
TH <- c(0, 0.1, 0.2, 0.5, 1, 2)
for (rule in c('step', 'graded', 'mean', 'decayed'))
  for (v in c('E9', 'E9b')) {
    r <- if (rule == 'decayed') 'decayed' else rule
    res <- sapply(TH, function(th) {
      x <- pw('Hybrid', 'CS', r, th, v); c(x$power, x$mean)
    })
    cat(sprintf('%-8s %-4s power: %s\n', rule, v, paste(f3(res[1, ]), collapse = ' ')))
    cat(sprintf('%-8s %-4s slope: %s\n', rule, v, paste(f3(res[2, ]), collapse = ' ')))
  }

cat('graded E9b slope t=1, 4 dp:', sprintf('%.4f', pw('Hybrid', 'CS', 'graded', 1, 'E9b')$mean), '\n')
hdr('docs/37 5.4: marginal power along the gap (Hybrid CS E9)')
gap_power <- function(gap, v = 'E9') {
  g <- lapply(D$Hybrid$paths, function(p) ifelse(p$on, 1, 1 - gap))
  m <- moments('Hybrid', 'CS', 'graded', 0, v, g_override = g, drop_mu = TRUE)
  pooled(m, D$Hybrid$n_p)$power
}
cat('gap 1, 0.9, 0.1, 0:', f3(sapply(c(1, 0.9, 0.1, 0), gap_power)), '\n')

hdr('docs/37 5.5/5.6: power at strengths 0.2, 0.3 (Hybrid CS)')
for (rule in c('graded', 'mean', 'decayed', 'step')) for (cb in c(0.2, 0.3)) {
  for (v in c('E9', 'E9b')) {
    res <- sapply(TH, function(th)
      pooled(moments('Hybrid', 'CS', rule, th, v, cb), D$Hybrid$n_p)$power)
    cat(sprintf('%-8s c=%.1f %-4s %s\n', rule, cb, v, paste(f3(res), collapse = ' ')))
  }
}
cat('graded CS Hybrid ceilings:', f3(sapply(TH, function(th)
  ceil_design(D$Hybrid, 'CS', 'graded', th))), '\n')

hdr('docs/37 5.7: OL+BDC and CO, CS, E9')
cat('OL+BDC graded:', f3(sapply(TH, function(th) pw('OL+BDC', 'CS', 'graded', th)$power)), '\n')
cat('OL+BDC step  :', f3(sapply(TH, function(th) pw('OL+BDC', 'CS', 'step', th)$power)), '\n')
for (v in c('E9', 'E9b'))
  cat('CO step pooled slope', v, f3(pw('CO', 'CS', 'step', 1, v)$mean), '\n')
TC <- c(0, 0.1, 0.2, 1, 2)
for (rule in c('step', 'graded', 'mean')) for (pt_ in c(TRUE, FALSE))
  cat(sprintf('CO %-7s path term %-5s %s\n', rule, pt_,
              paste(f3(sapply(TC, function(th)
                pw('CO', 'CS', rule, th, path_term = pt_)$power)), collapse = ' ')))
for (th in c(0, 2)) {
  m <- moments('CO', 'CS', 'mean', th, 'E9')
  x <- pooled(m, D$CO$n_p)
  cat('CO mean th', th, 'mu_p:', f3(m$mu), ' path term:', f3(x$pm),
      ' resid var (tau2 part):', f3(sum(0.5 * m$tau2)), '\n')
}
for (dn in c('Hybrid', 'OL+BDC')) {
  d <- max(abs(sapply(TH, function(th) pw(dn, 'CS', 'graded', th)$power -
                           pw(dn, 'CS', 'graded', th, path_term = FALSE)$power)))
  cat(dn, 'max change from removing path term (graded):', signif(d, 3), '\n')
}

hdr('docs/37 5.8: factorized variance, mean Var(Delta), power, ceilings')
## factorized form for the separable structure with constant eta (CO):
p <- D$CO$paths[[1]]; a <- weights(p, 'E9')
Ksep <- (1 - 0.2) * diag(3) + 0.2 * matrix(1, 3, 3)
A <- 0.8^abs(outer(p$w, p$w, '-')); vv <- c(s_tv, s_pb * 0.5, s_br)
cat('CO sep path 1: full', f3(var_D(p, 'sep', 'E9')), ' (v^T K v)(a^T A a):',
    f3(drop(t(vv) %*% Ksep %*% vv) * drop(t(a) %*% A %*% a)), '\n')
for (dn in names(D)) cat(dn, 'mean Var(Delta) CS/covar/sep:',
  sapply(c('CS', 'covar', 'sep'), function(s)
    f3(mean(sapply(D[[dn]]$paths, var_D, s = s, v = 'E9')))), '\n')
for (dn in names(D)) for (s in c('CS', 'covar', 'sep')) {
  st <- c(pw(dn, s, 'step', 0)$power, pw(dn, s, 'step', 1)$power)
  gr <- sapply(c(0, 0.5, 1, 2), function(th) pw(dn, s, 'graded', th)$power)
  mm <- sapply(c(0, 2), function(th) pw(dn, s, 'mean', th)$power)
  cat(sprintf('%-7s %-6s step %s | graded %s | mean %s\n', dn, s,
              paste(f3(st), collapse = '/'), paste(f3(gr), collapse = '/'),
              paste(f3(mm), collapse = '-')))
}
for (dn in names(D)) for (s in c('CS', 'covar', 'sep')) {
  cat(sprintf('%-7s %-6s ceiling step 0/1: %s  graded 0/0.5/1/2: %s\n', dn, s,
    paste(f3(c(ceil_design(D[[dn]], s, 'step', 0),
               ceil_design(D[[dn]], s, 'step', 1))), collapse = '/'),
    paste(f3(sapply(c(0, 0.5, 1, 2), function(th)
      ceil_design(D[[dn]], s, 'graded', th))), collapse = '/')))
}
## factored ceiling vs direct, separable
mx <- 0
for (dn in names(D)) for (th in c(0, 1)) for (rule in c('step', 'graded')) {
  d <- D[[dn]]
  fac <- (solve(Ksep)[3, 3] * max(sapply(d$paths, function(p) {
    A <- 0.8^abs(outer(p$w, p$w, '-')); g <- coupling(p, rule, th)
    drop(t(g) %*% solve(A, g)) })))^(-1/2)
  mx <- max(mx, abs(fac - ceil_design(d, 'sep', rule, th)))
}
cat('factored vs direct separable ceiling, max diff:', signif(mx, 3), '\n')

hdr('docs/37 4.4: pooled-slope variance, Monte Carlo of the stated formula')
## Delta_i = mu_p + beta_p (B - mean) + e, fixed allocation, OLS slope
set.seed(37)
n_p <- c(18, 18, 17, 17); mu <- c(0, 1.5, -1, 2); be <- c(-0.2, -0.1, -0.05, 0)
ta <- c(40, 42, 45, 47)
sl <- replicate(40000, {
  pid <- rep(1:4, n_p); B <- rnorm(70, 0, s_bm)
  y <- mu[pid] + be[pid] * B + rnorm(70, 0, sqrt(ta[pid]))
  coef(lm.fit(cbind(1, B), y))[2]
})
q <- n_p / 70; bb <- sum(q * be); mb <- sum(q * mu)
Vf <- sum(q * (ta + (mu - mb)^2)) / (67 * s_bm^2) + 2 * sum(q * (be - bb)^2) / 70
V1 <- sum(q * (ta + (mu - mb)^2 + s_bm^2 * (be - bb)^2)) / (67 * s_bm^2)
cat('MC var', signif(var(sl), 4), ' formula (2/N):', signif(Vf, 4),
    ' homoscedastic sigma_res^2 version:', signif(V1, 4), '\n')

## ======================================================================
hdr('docs/32: evenly spaced 8 visits 2 weeks apart, rho = 0.7')
pe <- list(paths = list(mk_path(seq(2, 16, 2), c(1,1,1,1,0,0,0,0), rep(1, 8))))
cat('orig t=0:', f3(ceil_design(pe, 'CS', 'step', 0, rho = 0.7)),
    ' orig t>0:', f3(ceil_design(pe, 'CS', 'step', 1, rho = 0.7)),
    ' covar graded t=1:', f3(ceil_design(pe, 'covar', 'graded', 1, rho = 0.7)), '\n')
cat('graded t=1: CS const', f3(ceil_design(pe, 'CS', 'graded', 1, rho = 0.7)),
    ' AR const', f3(ceil_design(pe, 'ARconst', 'graded', 1, rho = 0.7)),
    ' covar', f3(ceil_design(pe, 'covar', 'graded', 1, rho = 0.7)), '\n')
cat('rho=0.8 orig on this schedule:', f3(ceil_design(pe, 'CS', 'step', 0)),
    f3(ceil_design(pe, 'CS', 'step', 1)), '\n')

hdr('docs/34 4.1: crossed ceilings, rho = 0.7, t = 0/0.5/1')
for (s in c('covar', 'ARconst', 'CS', 'CSdecay')) for (rule in c('graded', 'step')) {
  cat(sprintf('%-8s %-6s', s, rule))
  for (dn in c('CO', 'Hybrid', 'OL+BDC'))
    cat(' |', paste(f3(sapply(c(0, 0.5, 1), function(th)
      ceil_design(D[[dn]], s, rule, th, rho = 0.7))), collapse = '/'))
  cat('\n')
}
hdr('docs/34 4.2: smallest eigenvalue of M, rho = 0.7 (and docs/36 at 0.8)')
for (rr in c(0.7, 0.8)) for (s in c('covar', 'ARconst', 'CS', 'CSdecay', 'sep')) {
  cat(sprintf('rho %.1f %-8s', rr, s))
  for (dn in c('CO', 'Hybrid', 'OL+BDC'))
    cat(' ', f3(min(sapply(D[[dn]]$paths, function(p)
      min(eigen(Mof(p, s, rho = rr), TRUE, TRUE)$values)))))
  cat('\n')
}
hdr('docs/34 4.4: rho sweep at t = 0, CO/Hybrid/OL+BDC')
for (rr in c(0.3, 0.5, 0.7, 0.8, 0.9)) cat('rho', rr, 'CS:',
  sapply(c('CO', 'Hybrid', 'OL+BDC'), function(dn)
    f3(ceil_design(D[[dn]], 'CS', 'graded', 0, rho = rr))), ' covar:',
  sapply(c('CO', 'Hybrid', 'OL+BDC'), function(dn)
    f3(ceil_design(D[[dn]], 'covar', 'graded', 0, rho = rr))), '\n')

hdr('docs/34 5.1: cross-factor parameters, Hybrid path C, rho = 0.7')
pC <- D$Hybrid$paths[[3]]
for (cc in list(c(.1, 0), c(.2, .1), c(.3, .2), c(.1, .1), c(.2, .2), c(.2, 0),
                c(.3, .1))) {
  g <- coupling(pC, 'graded', 0)
  Mc <- Mof(pC, 'CS', 0.7, cc[1], cc[2]); Ma <- Mof(pC, 'covar', 0.7, cc[1], cc[2])
  ca <- if (min(eigen(Ma, TRUE, TRUE)$values) > 0) f3(ceil_path(Ma, g)) else 'notPD'
  cat(cc, ' CS', f3(ceil_path(Mc, g)), ' AR', ca, '\n')
}
hdr('docs/34 5.2: synthetic schedules under AR(1) covar, rho = 0.7')
syn <- function(w, on) {
  p <- mk_path(w, on, rep(1, length(w)))
  M <- Mof(p, 'covar', 0.7)
  if (min(eigen(M, TRUE, TRUE)$values) <= 0) return(NA)
  ceil_path(M, as.numeric(on))
}
cat('visits 4..32, 1-wk, half on:', f3(sapply(c(4, 8, 16, 32), function(n)
  syn(seq_len(n), rep(1:0, each = n / 2)))), '\n')
cat('visits 4..32, 2.5-wk:', f3(sapply(c(4, 8, 16, 32), function(n)
  syn(2.5 * seq_len(n), rep(1:0, each = n / 2)))), '\n')
cat('CS visits 4..32 half on (rho .7):', f3(sapply(c(4, 8, 16, 32), function(n) {
  p <- mk_path(seq_len(n), rep(1:0, each = n / 2), rep(1, n))
  ceil_path(Mof(p, 'CS', 0.7), as.numeric(p$on)) })), '\n')
sw <- function(k, sp) {
  ## 16 visits, 8 on, k switches: blocks of equal length 16/(k+1) where possible
  pats <- list(`1` = rep(1:0, each = 8), `3` = rep(rep(1:0, each = 4), 2),
               `7` = rep(rep(1:0, each = 2), 4), `15` = rep(1:0, 8))
  syn(sp * (1:16), pats[[as.character(k)]])
}
cat('switches 1,3,7,15 at 1 wk:', f3(sapply(c(1, 3, 7, 15), sw, sp = 1)),
    ' at 2.5 wk:', f3(sapply(c(1, 3, 7, 15), sw, sp = 2.5)), '\n')
cat('uniform spacing 0.5,1,2.5,10 (8 visits, 4 on):', f3(sapply(c(0.5, 1, 2.5, 10),
  function(sp) syn(sp * (1:8), rep(1:0, each = 4)))), '\n')
cat('on-count 1..7 of 8, 2.5 wk:', f3(sapply(1:7, function(k)
  syn(2.5 * (1:8), c(rep(1, k), rep(0, 8 - k))))), '\n')
## switch gap at fixed study length 17.5 wk (8 visits, 4 on, one switch)
gapw <- function(d) {
  ## visits 1..4 then gap d then 5..8, total length 17.5 from first visit
  sp <- (17.5 - d) / 6
  w <- c(sp * (0:3), sp * 3 + d + sp * (0:3)) + 1
  syn(w, rep(1:0, each = 4))
}
cat('switch gap 0.5, 1, 2.5, 8 (fixed length):', f3(sapply(c(0.5, 1, 2.5, 8), gapw)), '\n')

hdr('docs/34 5.5: |r_t - r_s| bound, numeric check')
## max feasible r_t with r_s = 0 for a single pair at correlation rho:
for (rr in c(0.7, 0.9)) {
  f <- function(r) min(eigen(matrix(c(1, r, 0, r, 1, rr, 0, rr, 1), 3), TRUE, TRUE)$values)
  rmax <- uniroot(f, c(0, 0.999))$root
  cat('rho', rr, 'max |r_t - r_s| with r_s=0:', f3(rmax),
      ' sqrt(2(1-rho)):', f3(sqrt(2 * (1 - rr))), ' sqrt(1-rho^2):',
      f3(sqrt(1 - rr^2)), '\n')
}

hdr('docs/34 5.6: Weibull decay shapes, covar, rho = 0.7')
D$Multi <- list(paths = list(mk_path(1:16, rep(c(1, 1, 0, 0), 4), rep(1, 16))))
for (dn in c('Hybrid', 'OL+BDC', 'Multi', 'CO')) for (th in c(0.5, 1, 2)) {
  cat(sprintf('%-7s th %.1f:', dn, th), f3(sapply(c(0.25, 0.5, 1, 2, 4), function(k)
    ceil_design(D[[dn]], 'covar', 'graded', th, rho = 0.7, k = k))),
    ' none:', f3(ceil_design(D[[dn]], 'covar', 'graded', 0, rho = 0.7)), '\n')
}
cat('step k=4 th=0.5 Hybrid/OL+BDC:', f3(ceil_design(D$Hybrid, 'CS', 'step', 0.5, 0.7, k = 4)),
    f3(ceil_design(D$`OL+BDC`, 'CS', 'step', 0.5, 0.7, k = 4)),
    ' th=1 Hybrid:', f3(ceil_design(D$Hybrid, 'CS', 'step', 1, 0.7, k = 4)), '\n')
cat('step k<=2 th=0.5 Hybrid:', f3(ceil_design(D$Hybrid, 'CS', 'step', 0.5, 0.7, k = 2)), '\n')

## ======================================================================
hdr('docs/35 4.4, 4.5: covar eigenvalues, Toeplitz thresholds, orig levels')
for (cc in list(c(.2, .1), c(.3, .2), c(.1, .05))) {
  k <- (cc[1] - cc[2]) / (1 - cc[2]); phi <- (1 - k) / (1 + k)
  cat(cc, ' threshold k', f3(k), ' min gap rho .7:', f3(log(phi) / log(0.7)),
      ' rho .8:', f3(log(phi) / log(0.8)), '\n')
}
## explicit eigenvalues of covar
p <- D$Hybrid$paths[[1]]; A <- 0.8^abs(outer(p$w, p$w, '-'))
om <- eigen(A, TRUE, TRUE)$values
pred <- sort(c(om * (1 + 2 * 0.1) + 2 * 0.1, rep(om * 0.9 - 0.1, 2)))
cat('covar eigenvalue formula vs eigen, max diff:',
    signif(max(abs(pred - sort(eigen(Mof(p, 'covar'), TRUE, TRUE)$values))), 3), '\n')
## half-week spacing failures, rho = 0.7
for (cc in list(c(.2, .1), c(.3, .2), c(.1, .05))) for (sp in c(0.25, 0.5, 1)) {
  w <- sp * (1:16); A <- 0.7^abs(outer(w, w, '-'))
  cat(cc, 'spacing', sp, ' omega_min', f3(min(eigen(A, TRUE, TRUE)$values)),
      ' threshold', f3((cc[1] - cc[2]) / (1 - cc[2])), '\n')
}
for (rr in c(0.7, 0.8)) cat('rho', rr, 'orig person/occasion cross cor:',
  f3(0.1 / rr), f3(0.1 / (1 - rr)), '\n')
cat('2-component [K^-1]_22 at c1=.2:', f3(solve(matrix(c(1, .2, .2, 1), 2))[2, 2]),
    ' 1/(1-c1^2):', f3(1 / (1 - 0.04)), '\n')
hdr('docs/35 6: covar vs separable ceilings, rho = 0.7, t = 0')
cp <- function(p, s) { M <- Mof(p, s, 0.7)
  if (min(eigen(M, TRUE, TRUE)$values) <= 0) 'notPD' else f3(ceil_path(M, as.numeric(p$on))) }
for (nm in c('CO A', 'Hybrid A', 'Hybrid C', 'OL+BDC A')) {
  p <- switch(nm, `CO A` = D$CO$paths[[1]], `Hybrid A` = D$Hybrid$paths[[1]],
              `Hybrid C` = D$Hybrid$paths[[3]], `OL+BDC A` = D$`OL+BDC`$paths[[1]])
  cat(nm, cp(p, 'covar'), cp(p, 'sep'), '\n')
}
cat('Multi-cycle', cp(D$Multi$paths[[1]], 'covar'), cp(D$Multi$paths[[1]], 'sep'), '\n')
for (nn in c(8, 16)) for (pat in list('half', 'alt')) {
  on <- if (pat == 'half') rep(1:0, each = nn / 2) else rep(c(1, 1, 0, 0), nn / 4)
  pd <- mk_path(0.5 * (1:nn), on, rep(1, nn))
  cat('dense 0.5 n', nn, pat, cp(pd, 'covar'), cp(pd, 'sep'), '\n')
}

## ======================================================================
hdr('docs/36 6.5: ceilings A/B/F/C at rho = 0.8, OL schedule candidates')
for (wn in c('OB', 'H')) {
  w <- if (wn == 'OB') W_OB else W_H
  pol <- list(paths = list(mk_path(w, rep(1, 8), rep(1, 8))))
  cat('OL on', wn, 'weeks: A', f3(ceil_design(pol, 'CS', 'step', 0)),
      ' B/F', f3(ceil_design(pol, 'covar', 'step', 0)),
      ' C', f3(ceil_design(pol, 'sep', 'step', 0)), '\n')
}
TH36 <- c(0, 0.1, 0.2, 0.5, 1)
for (cf in list(c('A', 'CS', 'step'), c('B', 'covar', 'step'), c('F', 'covar', 'graded'),
                c('C', 'sep', 'graded'), c('adj2', 'CS', 'graded')))
  for (dn in c('OL+BDC', 'CO', 'Hybrid'))
    cat(cf[1], dn, f3(sapply(TH36, function(th) ceil_design(D[[dn]], cf[2], cf[3], th))), '\n')
hdr('docs/36 6.4: covar threshold at published values')
cat('(c1-cx)/(1-cx) =', f3(0.1 / 0.9), ' min gap at rho .8 =', f3(log(0.8) / log(0.8)), '\n')
hdr('docs/36 repair counts from ceilings')
D$OL <- list(paths = list(mk_path(W_CO, rep(1, 8), rep(1, 8))))
count_fail <- function(s, rule, cbs, ths) {
  tot <- 0; fl <- 0
  for (dn in c('OL', 'OL+BDC', 'CO', 'Hybrid')) for (p in D[[dn]]$paths)
    for (cb in cbs) for (th in ths) {
      tot <- tot + 1
      if (rule == 'mean') next
      M <- Mof(p, s); c_ <- ceil_path(M, coupling(p, rule, th))
      if (cb >= c_) fl <- fl + 1
    }
  c(fl, tot)
}
for (cf in list(c('A', 'CS', 'step'), c('B', 'covar', 'step'), c('F', 'covar', 'graded'),
                c('C', 'sep', 'graded'), c('D', 'CS', 'mean'))) {
  a1 <- count_fail(cf[2], cf[3], c(0.3, 0.6), c(0, 0.1, 0.2))
  a2 <- count_fail(cf[2], cf[3], c(0.25, 0.3), c(0, 0.1, 0.2, 0.5, 1))
  cat(cf[1], ' published grid', a1[1], '/', a1[2], ' common grid', a2[1], '/', a2[2], '\n')
}
hdr('docs/36 A.5/A.7: four-visit inventory part (weeks 10, 11, 12, 16)')
fw <- c(10, 11, 12, 16); A <- 0.8^abs(outer(fw, fw, '-'))
for (s in c('CS', 'covar', 'sep')) {
  M <- build_M(fw, 0.8, 0.2, 0.1, STRUCT[[s]][1], STRUCT[[s]][2])
  cat(s, 'min eig', f3(min(eigen(M, TRUE, TRUE)$values)), '\n')
}
cat('eig(A):', f3(eigen(A, TRUE, TRUE)$values), '\n')
phi <- 0.8^c(1, 1, 4)
steps <- function(u) c(u[1]^2, (u[-1] - phi * u[-4])^2 / (1 - phi^2))
Msep <- build_M(fw, 0.8, 0.2, 0.1, 'AR', 'sep')
for (u in list(c(1, 0, 0, 1), c(1, 0.5, 0.25, 1), c(1, 1, 1, 1))) {
  st <- steps(u); qf <- sum(st)
  cat('u', u, ' steps', f3(st), ' sum', f3(qf), ' ceiling', sprintf('%.4f',
      (solve(Ksep)[3, 3] * qf)^(-1/2)), ' full', sprintf('%.4f', ceil_path(Msep, u)), '\n')
}
u <- c(1, 0.5, 0.25, 1); Kinv <- solve(Ksep)
Wm <- 0.3 * outer(Kinv[, 3], solve(A, u))
cat('W rows:\n'); print(round(Wm, 4))
cat('1 - c^2 q:', sprintf('%.4f', 1 - 0.09 * Kinv[3, 3] * drop(t(u) %*% solve(A, u))), '\n')
pA <- D$Hybrid$paths[[1]]
cat('Hybrid paths 1-4 graded t=0.5 sep ceiling (rho .8):', sprintf('%.6f',
    sapply(D$Hybrid$paths, function(p) ceil_path(Mof(p, 'sep'), coupling(p, 'graded', 0.5)))), '\n')
cat('same at rho .7:', sprintf('%.6f',
    sapply(D$Hybrid$paths, function(p) ceil_path(Mof(p, 'sep', 0.7), coupling(p, 'graded', 0.5)))), '\n')
cat('OL on CO weeks (constant g): B/F', f3(ceil_path(Mof(D$CO$paths[[1]], 'covar'), rep(1, 8))),
    ' C', f3(ceil_path(Mof(D$CO$paths[[1]], 'sep'), rep(1, 8))), '\n')

## ======================================================================
hdr('docs/38 / sensitivity binding: R^2 of B on the BR trajectory')
set.seed(38)
fvec <- c(0.3, 0.7, 0.95, 1, 0.5, 0.25, 0.9, 0.2)
G <- 10 * 0.3; sv <- 2; sx <- 3; H <- G^2 + sv^2
Rx <- 0.5^abs(outer(1:8, 1:8, '-')); Sxi <- sx^2 * Rx
Sbr <- H * outer(fvec, fvec) + Sxi
qv <- drop(t(fvec) %*% solve(Sxi, fvec))
R2_direct <- G^2 * drop(t(fvec) %*% solve(Sbr, fvec))
cat('R2 direct', signif(R2_direct, 6), ' G^2 q/(1+Hq)', signif(G^2 * qv / (1 + H * qv), 6),
    ' bound G^2/H', signif(G^2 / H, 6), '\n')
## joint correlation matrix of (B, BR) and its smallest eigenvalue
Cbb <- c(1, G * fvec)
Sj <- rbind(Cbb, cbind(G * fvec, Sbr))
cat('min eigen joint covariance (B standardized):', signif(min(eigen(Sj)$values), 4), '\n')
## conditional moments given b (simulate)
nn <- 4e5; b <- rnorm(nn); v <- rnorm(nn, 0, sv)
xi <- mvrnorm(nn, rep(0, 8), Sxi)
om <- 10 * (1 + 0.3 * b) + v
BR <- outer(om, fvec) + xi
fit <- lm(BR[, 4] ~ b)
cat('E(BR_4 | b) slope', signif(coef(fit)[2], 4), ' expected G f_4', G * fvec[4], '\n')
res <- BR - outer(10 + G * b, fvec)
cat('max |Cov(BR | b) - (sv^2 f f^T + Sxi)|:',
    signif(max(abs(cov(res) - (sv^2 * outer(fvec, fvec) + Sxi))), 3), '\n')
## implied coupling relative to the peak
cr <- G * fvec / sqrt(H * fvec^2 + sx^2)
cat('Cor(B,BR_t):', f3(cr), '\n')
cat('coupling ratio Cor/Cor at f=1:', f3(cr / (G / sqrt(H + sx^2))), ' vs f:', f3(fvec), '\n')

hdr('docs/36 6.7 check vector: source used 0.0625 at visit 8 (tsd = 4)')
u_src <- c(1, 1, 1, 1, 0.25, 0.0625, 1, 0.0625)
cat('ceiling with source vector', sprintf('%.6f', ceil_path(Mof(pA, 'sep'), u_src)),
    ' correct vector', sprintf('%.6f', ceil_path(Mof(pA, 'sep'), coupling(pA, 'graded', 0.5))),
    ' correct g:', signif(coupling(pA, 'graded', 0.5), 4), '\n')

hdr('docs/37 6: validation summary recomputed from stored output')
dd <- 'analysis/data/quick-sim/carryover-closed-form'
val <- read.csv(file.path(dd, 'e9-closed-form-vs-mc.csv'))
cat('rows:', nrow(val), ' distinct cells:',
    nrow(unique(val[, c('design', 'construct', 'rule', 'c.bm', 't_half')])), '\n')
for (bk in unique(val$block)) {
  s <- val[val$block == bk, ]
  cat(sprintf('%-50s n=%3d max|z|=%.2f mean z=%.2f beyond2=%d maxdiff=%.4f\n', bk, nrow(s),
      max(abs(s$power_z)), mean(s$power_z), sum(abs(s$power_z) > 2),
      max(abs(s$mean - s$mc_mean))))
}
cat('overall beyond 2:', sum(abs(val$power_z) > 2), 'of', nrow(val), '\n')
h <- val[val$design == 'Hybrid' & val$construct == 'CS' & val$c.bm == 0.25 & val$variant == 'E9', ]
h <- h[order(h$rule, h$t_half), c('rule', 't_half', 'mean', 'mc_mean', 'sd', 'mc_sd', 'power', 'mc_power')]
print(format(h, digits = 4), row.names = FALSE)
lm_ <- read.csv(file.path(dd, 'e9-vs-lmer.csv'))
print(names(lm_)); print(format(lm_, digits = 3), row.names = FALSE)

hdr('docs/34 5.3, 5.4: lookup grid corners and equalized paper 01 schedules')
for (sp in c(1, 1.5, 2.5, 5)) cat('16 visits spacing', sp, ' 1 switch', f3(sw(1, sp)),
                                 ' 15 switches', f3(sw(15, sp)), '\n')
eq <- function(p) { w <- seq(min(p$w), max(p$w), length.out = length(p$w)); syn(w, as.numeric(p$on)) }
cat('Hybrid C as scheduled / equalized:', f3(syn(D$Hybrid$paths[[3]]$w, as.numeric(D$Hybrid$paths[[3]]$on))),
    f3(eq(D$Hybrid$paths[[3]])), '\n')
cat('OL+BDC A as scheduled / equalized:', f3(syn(D$`OL+BDC`$paths[[1]]$w, as.numeric(D$`OL+BDC`$paths[[1]]$on))),
    f3(eq(D$`OL+BDC`$paths[[1]])), '\n')

hdr('docs/37 5.8: step power at t = 0.1, 1, 2 and mean-rule range over t')
for (dn in c('Hybrid', 'OL+BDC', 'CO')) for (s in c('CS', 'covar', 'sep'))
  cat(sprintf('%-7s %-6s step 0.1/1/2: %s  mean 0..2: %s\n', dn, s,
    paste(f3(sapply(c(0.1, 1, 2), function(th) pw(dn, s, 'step', th)$power)), collapse = '/'),
    paste(f3(sapply(TH, function(th) pw(dn, s, 'mean', th)$power)), collapse = ' ')))
