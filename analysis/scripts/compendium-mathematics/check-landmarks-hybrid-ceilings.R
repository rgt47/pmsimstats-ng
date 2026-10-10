## Copied 2026-10-08 from the review scratch script inv_c_check.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-landmarks-hybrid-ceilings.R
# Inventory part c: checks for docs/29, 30, 31, 33, 40.
# Run from the repository root: Rscript analysis/scripts/compendium-mathematics/check-landmarks-hybrid-ceilings.R
# Self-contained; reads no repository data.

gomp <- function(t, m, d, r) {
  m * (exp(-d * exp(-r * t)) - exp(-d)) / (1 - exp(-d))
}
landmark <- function(frac, d, r) {
  uniroot(function(t) gomp(t, 1, d, r) - frac, c(0, 100), tol = 1e-12)$root
}

cat('== 1. Gompertz landmarks (docs/31 s2; docs/40 s4) ==\n')
for (pr in list(c(10.99, 5, 0.42), c(6.51, 5, 0.35))) {
  m <- pr[1]; d <- pr[2]; r <- pr[3]
  cat(sprintf(paste0('m %.2f d %g r %.2f: inflection %.3f, share at ',
    'inflection %.4f, half %.3f, 90%% %.3f, share wk1/4/8 %.3f %.3f %.3f\n'),
    m, d, r, log(d) / r, (exp(-1) - exp(-d)) / (1 - exp(-d)),
    landmark(0.5, d, r), landmark(0.9, d, r),
    gomp(1, 1, d, r), gomp(4, 1, d, r), gomp(8, 1, d, r)))
}
cat(sprintf('half-max of BR in days %.1f\n', 7 * landmark(0.5, 5, 0.42)))
cat(sprintf('rate doubling: half-time at r = 0.84 %.4f vs half of 0.42 value %.4f\n',
  landmark(0.5, 5, 0.84), landmark(0.5, 5, 0.42) / 2))
cat(sprintf('d -> 0 limit at t = 3, r = 0.42: %.6f vs 1 - exp(-rt) %.6f\n',
  gomp(3, 1, 1e-7, 0.42), 1 - exp(-0.42 * 3)))

cat('\n== 2. Maximum slope of the rescaled Gompertz as a function of d ==\n')
max_slope_closed <- function(m, d, r) {
  if (d >= 1) m * r / (exp(1) * (1 - exp(-d))) else m * r * d / (exp(d) - 1)
}
max_slope_num <- function(m, d, r) {
  tt <- seq(0, 60, by = 1e-4)
  dg <- m * d * r * exp(-r * tt) * exp(-d * exp(-r * tt)) / (1 - exp(-d))
  max(dg)
}
for (d in c(0.25, 0.5, 1, 2, 5, 10, 20)) {
  cat(sprintf('d %5.2f: closed %.5f numeric %.5f (m = 1, r = 1)\n', d,
    max_slope_closed(1, d, 1), max_slope_num(1, d, 1)))
}
cat(sprintf('BR at d = 5: max slope %.4f points/wk at t* = %.3f\n',
  max_slope_closed(10.99, 5, 0.42), log(5) / 0.42))

cat('\n== 3. Hybrid component means by week (docs/40 s3) ==\n')
# Exact shipped maxima from data/extracted_rp.rda (check-extracted-parameters.R prints them).
m_tv <- 6.50647; m_br <- 10.98604
w <- c(4, 8, 9, 10, 11, 12, 16, 20)
eta <- c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
tv <- gomp(w, m_tv, 5, 0.35)
pb <- eta * gomp(w, m_tv, 5, 0.35)
paths <- list(
  c(1, 1, 1, 1, 0, 0, 1, 0), c(1, 1, 1, 1, 0, 0, 0, 1),
  c(1, 1, 1, 0, 0, 0, 1, 0), c(1, 1, 1, 0, 0, 0, 0, 1))
tod_of <- function(w, on) {
  out <- numeric(length(w)); start <- NA
  prev_on <- 0; prev_w <- 0
  for (k in seq_along(w)) {
    if (on[k] == 1) {
      if (prev_on == 0) start <- prev_w
      out[k] <- w[k] - start
    }
    prev_on <- on[k]; prev_w <- w[k]
  }
  out
}
tsd_of <- function(w, on) {
  out <- rep(NA_real_, length(w)); last_on <- NA
  for (k in seq_along(w)) {
    if (on[k] == 1) {
      out[k] <- 0; last_on <- w[k]
    } else if (!is.na(last_on)) {
      out[k] <- w[k] - last_on
    }
  }
  out
}
cat('TV ', sprintf('%.2f', tv), '\n')
cat('PB ', sprintf('%.2f', pb), '\n')
for (p in seq_along(paths)) {
  tod <- tod_of(w, paths[[p]])
  br <- ifelse(paths[[p]] == 1, gomp(tod, m_br, 5, 0.42), 0)
  cat(sprintf('path %d tod %s | BR %s\n', p, paste(tod, collapse = ' '),
    paste(sprintf('%.2f', br), collapse = ' ')))
}
cat(sprintf('total at week 8: %.2f; at week 9: %.2f\n',
  tv[2] + pb[2] + gomp(8, m_br, 5, 0.42), tv[3] + pb[3] + gomp(9, m_br, 5, 0.42)))
cat(sprintf('expectancy drop week 8 to 9 in PB: %.2f; 0.5*G(9) = %.2f\n',
  pb[2] - pb[3], 0.5 * tv[3]))

cat('\n== 4. Per-path E9 slopes (docs/40 s5) ==\n')
s_br <- 8; s_bm <- 15.36; cbm <- 0.25
for (th in c(0, 1)) {
  sl <- sapply(paths, function(on) {
    tsd <- tsd_of(w, on)
    g <- ifelse(on == 1, 1, if (th > 0) 2^(-tsd / th) else 0)
    -cbm * s_br / s_bm * (1 - mean(g[on == 0]))
  })
  cat(sprintf('graded t1/2 %g: %s\n', th, paste(sprintf('%.4f', sl), collapse = ' ')))
}
cat(sprintf('2^-10 = %.6f\n', 2^-10))

cat('\n== 5. Outcome covariance under the CS construct, path 1, c_bm = 0 ==\n')
n <- 8
s_tv <- rep(10, n); s_pb <- 10 * eta; s_brv <- rep(8, n); s_bl <- 18.48
rho <- 0.8; c1 <- 0.2; cx <- 0.1
Ablk <- matrix(rho, n, n); diag(Ablk) <- 1
Cblk <- matrix(cx, n, n); diag(Cblk) <- c1
Rr <- rbind(cbind(Ablk, Cblk, Cblk), cbind(Cblk, Ablk, Cblk),
  cbind(Cblk, Cblk, Ablk))
sv <- c(s_tv, s_pb, s_brv)
Sig <- diag(sv) %*% Rr %*% diag(sv)
L <- cbind(diag(n), diag(n), diag(n))
covY <- s_bl^2 + L %*% Sig %*% t(L)
cat(sprintf('Var(BL) %.1f\n', s_bl^2))
cat(sprintf('Var OL %.1f, Var blinded %.1f\n', covY[1, 1], covY[3, 3]))
cat(sprintf('Cov within OL %.1f, within blinded %.1f (range %.1f-%.1f), OL-blinded %.1f\n',
  covY[1, 2], covY[3, 4], min(covY[3:8, 3:8][upper.tri(covY[3:8, 3:8])]),
  max(covY[3:8, 3:8][upper.tri(covY[3:8, 3:8])]), covY[1, 3]))
st <- c(10, 10, 8); ss <- c(10, 5, 8)
Kss <- matrix(c1, 3, 3); diag(Kss) <- 1
Kts <- matrix(cx, 3, 3); diag(Kts) <- rho
cat(sprintf('closed: Var OL %.1f Var BD %.1f Cov OL-OL %.1f BD-BD %.1f OL-BD %.1f\n',
  s_bl^2 + t(st) %*% Kss %*% st, s_bl^2 + t(ss) %*% Kss %*% ss,
  s_bl^2 + t(st) %*% Kts %*% st, s_bl^2 + t(ss) %*% Kts %*% ss,
  s_bl^2 + t(st) %*% Kts %*% ss))

cat('\n== 6. Hybrid long data: Cor(De, Db), R^2 of De on Db and t ==\n')
alloc <- c(18, 18, 17, 17)
rows <- do.call(rbind, lapply(seq_along(paths), function(p) {
  data.frame(wt = alloc[p], t = c(0, w), Db = c(0, paths[[p]]), De = c(0, eta))
}))
wcor <- function(x, y, wt) {
  mx <- sum(wt * x) / sum(wt); my <- sum(wt * y) / sum(wt)
  sum(wt * (x - mx) * (y - my)) / sqrt(sum(wt * (x - mx)^2) * sum(wt * (y - my)^2))
}
fit <- lm(De ~ Db + t, data = rows, weights = wt)
r2 <- summary(fit)$r.squared
cat(sprintf('Cor(De, Db) %.4f; R^2 %.4f; VIF %.3f\n',
  wcor(rows$De, rows$Db, rows$wt), r2, 1 / (1 - r2)))
rows_nb <- rows[rows$t > 0, ]
fit2 <- lm(De ~ Db + t, data = rows_nb, weights = wt)
cat(sprintf('without baseline rows: Cor %.4f; R^2 %.4f\n',
  wcor(rows_nb$De, rows_nb$Db, rows_nb$wt), summary(fit2)$r.squared))

cat('\n== 7. K-component compound-symmetry spectrum (docs/31 s20) ==\n')
cs_MK <- function(K, n, rho, c1, cx) {
  A <- matrix(rho, n, n); diag(A) <- 1
  C <- matrix(cx, n, n); diag(C) <- c1
  M <- kronecker(diag(K), A) + kronecker(matrix(1, K, K) - diag(K), C)
  M
}
spec_closed <- function(K, n, rho, c1, cx) {
  ev <- c((1 + (n - 1) * rho) + (K - 1) * (c1 + (n - 1) * cx),
    (1 - rho) + (K - 1) * (c1 - cx))
  if (K >= 2) ev <- c(ev, (1 + (n - 1) * rho) - (c1 + (n - 1) * cx),
    (1 - rho) - (c1 - cx))
  sort(ev)
}
for (K in 1:4) {
  for (nn in c(4, 8, 12)) {
    ev <- eigen(cs_MK(K, nn, 0.7, 0.2, 0.1), symmetric = TRUE,
      only.values = TRUE)$values
    cl <- spec_closed(K, nn, 0.7, 0.2, 0.1)
    ok <- all(sapply(round(unique(round(ev, 10)), 8),
      function(e) any(abs(e - cl) < 1e-8)))
    cat(sprintf('K %d n %2d: min eig %.4f closed min %.4f distinct %s all-match %s\n',
      K, nn, min(ev), min(cl),
      paste(sprintf('%.2f', sort(unique(round(ev, 8)), decreasing = TRUE)),
        collapse = ' '), ok))
  }
}
cat('counterexample c_x > rho (rho 0.1, c1 0.2, cx 0.15, K 3, n 8):\n')
ev <- eigen(cs_MK(3, 8, 0.1, 0.2, 0.15), symmetric = TRUE, only.values = TRUE)$values
cat(sprintf('  min eig %.4f vs (1-rho)-(c1-cx) = %.4f\n', min(ev), 0.9 - 0.05))
cat('counterexample c1 < cx (rho 0.7, c1 0.05, cx 0.2, K 3, n 8):\n')
ev <- eigen(cs_MK(3, 8, 0.7, 0.05, 0.2), symmetric = TRUE, only.values = TRUE)$values
cat(sprintf('  min eig %.4f vs (1-rho)-(c1-cx) = %.4f, (1-rho)+(K-1)(c1-cx) = %.4f\n',
  min(ev), 0.3 + 0.15, 0.3 + 2 * (0.05 - 0.2)))

cat('\n== 8. Ceilings under the non-separable AR(1) construct, K = 3 vs 2 ==\n')
covar_M <- function(K, w, rho, c1, cx) {
  A <- rho^abs(outer(w, w, '-'))
  C <- cx * A; diag(C) <- c1
  kronecker(diag(K), A) + kronecker(matrix(1, K, K) - diag(K), C)
}
graded <- function(w, on, th) {
  tsd <- tsd_of(w, on)
  g <- numeric(length(w))
  g[on == 1] <- 1
  off_after <- on == 0 & !is.na(tsd)
  if (th > 0) g[off_after] <- 2^(-tsd[off_after] / th)
  g
}
ceil <- function(M, g, K) {
  u <- c(rep(0, (K - 1) * length(g)), g)
  1 / sqrt(drop(t(u) %*% solve(M, u)))
}
hyb_w <- w
olb_w <- c(4, 8, 12, 16, 17, 18, 19, 20)
co_w <- seq(2.5, 20, 2.5)
designs <- list(
  CO = list(co_w, list(c(1, 1, 1, 1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, 1, 1, 1))),
  Hybrid = list(hyb_w, paths),
  'OL+BDC' = list(olb_w, list(c(1, 1, 1, 1, 1, 1, 0, 0), c(1, 1, 1, 1, 1, 0, 0, 0))))
gains <- c()
for (dn in names(designs)) {
  dw <- designs[[dn]][[1]]
  for (th in c(0, 0.5, 1)) {
    cK <- sapply(c(3, 2), function(K) {
      M <- covar_M(K, dw, 0.7, 0.2, 0.1)
      min(sapply(designs[[dn]][[2]], function(on) ceil(M, graded(dw, on, th), K)))
    })
    gains <- c(gains, cK[2] - cK[1])
    cat(sprintf('%-7s t1/2 %.1f: K=3 %.3f K=2 %.3f gain %.4f\n', dn, th,
      cK[1], cK[2], cK[2] - cK[1]))
  }
}
cat(sprintf('gain range %.4f to %.4f\n', min(gains), max(gains)))

cat('\n== 9. Separable K-component ceiling ratio (derived here) ==\n')
r2K <- function(K, c1) (K - 1) * c1^2 / (1 + (K - 2) * c1)
for (c1v in c(0.1, 0.2, 0.3)) {
  for (K in 2:4) {
    Kc <- matrix(c1v, K, K); diag(Kc) <- 1
    cat(sprintf('c1 %.1f K %d: R2 closed %.5f vs 1 - 1/[Kinv]KK %.5f\n', c1v, K,
      r2K(K, c1v), 1 - 1 / solve(Kc)[K, K]))
  }
  cat(sprintf('  ratio c*(K=2)/c*(K=3) = %.5f\n',
    sqrt((1 - r2K(2, c1v)) / (1 - r2K(3, c1v)))))
}
sep_ceil <- function(K, w, g, rho, c1) {
  Kc <- matrix(c1, K, K); diag(Kc) <- 1
  M <- kronecker(Kc, rho^abs(outer(w, w, '-')))
  ceil(M, g, K)
}
for (dn in names(designs)) {
  dw <- designs[[dn]][[1]]
  c3 <- min(sapply(designs[[dn]][[2]], function(on) sep_ceil(3, dw, graded(dw, on, 0), 0.7, 0.2)))
  c2 <- min(sapply(designs[[dn]][[2]], function(on) sep_ceil(2, dw, graded(dw, on, 0), 0.7, 0.2)))
  cat(sprintf('separable rho 0.7 %-7s t1/2 0: K=3 %.4f K=2 %.4f ratio %.5f\n',
    dn, c3, c2, c2 / c3))
}

cat('\n== 10. CS versus separable ceilings under graded coupling, rho 0.8 (docs/33 s3.4) ==\n')
for (dn in names(designs)) {
  dw <- designs[[dn]][[1]]
  for (th in c(0, 0.5, 1, 2)) {
    ccs <- min(sapply(designs[[dn]][[2]], function(on)
      ceil(cs_MK(3, 8, 0.8, 0.2, 0.1), graded(dw, on, th), 3)))
    csep <- min(sapply(designs[[dn]][[2]], function(on)
      sep_ceil(3, dw, graded(dw, on, th), 0.8, 0.2)))
    cat(sprintf('%-7s t1/2 %.1f: CS %.4f separable %.4f\n', dn, th, ccs, csep))
  }
}

cat('\n== 11. First off-drug coupling by design (docs/33 s2.5) ==\n')
for (th in c(0.5, 1)) {
  cat(sprintf('t1/2 %.1f: CO (2.5 wk) %.4f ; Hybrid/OL+BDC (1 wk) %.4f\n',
    th, 2^(-2.5 / th), 2^(-1 / th)))
}

cat('\n== 12. Conditional SD 8 sqrt(1 - 0.45^2) (docs/33 s4) ==\n')
cat(sprintf('%.4f\n', 8 * sqrt(1 - 0.45^2)))

cat('\n== 13. Working-independence information loss attempts (docs/33 s2.6) ==\n')
dbc <- function(w, on, th) graded(w, on, th)
ss_path <- function(w, on, th, base) {
  x <- dbc(w, on, th)
  if (base) x <- c(0, x)
  sum((x - mean(x))^2)
}
for (base in c(FALSE, TRUE)) {
  for (th in c(0.5, 1)) {
    out <- sapply(names(designs), function(dn) {
      dw <- designs[[dn]][[1]]
      s0 <- sum(sapply(designs[[dn]][[2]], function(on) ss_path(dw, on, 0, base)))
      s1 <- sum(sapply(designs[[dn]][[2]], function(on) ss_path(dw, on, th, base)))
      1 - s1 / s0
    })
    cat(sprintf('baseline %s t1/2 %.1f: %s\n', base, th,
      paste(sprintf('%s %.1f%%', names(out), 100 * out), collapse = '; ')))
  }
}
