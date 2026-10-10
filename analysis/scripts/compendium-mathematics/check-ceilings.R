ar1 <- function(w, rho) rho^abs(outer(w, w, '-'))
equi <- function(d, o) {
  K <- matrix(o, 3, 3)
  diag(K) <- d
  K
}
cs_M <- function(n, rho, c1, cx) {
  A <- matrix(rho, n, n)
  diag(A) <- 1
  C <- matrix(cx, n, n)
  diag(C) <- c1
  rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
}
ceil_full <- function(M, g) {
  u <- c(rep(0, 2 * length(g)), g)
  1 / sqrt(drop(t(u) %*% solve(M, u)))
}
ceil_sep_closed <- function(w, g, rho, c1) {
  phi <- rho^diff(w)
  q <- g[1]^2 + sum((g[-1] - phi * g[-length(g)])^2 / (1 - phi^2))
  r2 <- 2 * c1^2 / (1 + c1)
  sqrt(1 - r2) / sqrt(q)
}
cs_closed <- function(n, n_on, rho, c1, cx) {
  K1 <- equi(1 - rho, c1 - cx)
  K2 <- equi(rho, cx)
  a <- solve(K1)[3, 3]
  an <- solve(K1 + n * K2)[3, 3]
  1 / sqrt(a * n_on * (n - n_on) / n + an * n_on^2 / n)
}
min_eig <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

hyb_w <- c(4, 8, 9, 10, 11, 12, 16, 20)
olb_w <- c(4, 8, 12, 16, 17, 18, 19, 20)
co_w <- seq(2.5, 20, 2.5)
paths <- list(
  'CO A' = list(co_w, c(1, 1, 1, 1, 0, 0, 0, 0)),
  'CO B' = list(co_w, c(0, 0, 0, 0, 1, 1, 1, 1)),
  'Hybrid A' = list(hyb_w, c(1, 1, 1, 1, 0, 0, 1, 0)),
  'Hybrid B' = list(hyb_w, c(1, 1, 1, 1, 0, 0, 0, 1)),
  'Hybrid C' = list(hyb_w, c(1, 1, 1, 0, 0, 0, 1, 0)),
  'Hybrid D' = list(hyb_w, c(1, 1, 1, 0, 0, 0, 0, 1)),
  'OL+BDC A' = list(olb_w, c(1, 1, 1, 1, 1, 1, 0, 0)),
  'OL+BDC B' = list(olb_w, c(1, 1, 1, 1, 1, 0, 0, 0))
)

cat('== CS (published) ceilings, no carryover, c1 = 0.2, cx = 0.1 ==\n')
for (rho in c(0.8, 0.7)) {
  for (nm in names(paths)) {
    g <- paths[[nm]][[2]]
    M <- cs_M(8, rho, c1 = 0.2, cx = 0.1)
    cat(sprintf('rho %.1f %-9s full %.4f closed %.4f minEigM %.4f\n', rho, nm,
      ceil_full(M, g), cs_closed(8, sum(g), rho, 0.2, 0.1), min_eig(M)))
  }
}

cat('\n== Separable AR(1) ceilings, no carryover, c1 = 0.2 ==\n')
K <- equi(1, 0.2)
for (rho in c(0.8, 0.7)) {
  for (nm in names(paths)) {
    w <- paths[[nm]][[1]]
    g <- paths[[nm]][[2]]
    M <- kronecker(K, ar1(w, rho))
    cat(sprintf('rho %.1f %-9s full %.4f closed %.4f\n', rho, nm,
      ceil_full(M, g), ceil_sep_closed(w, g, rho, 0.2)))
  }
}

cat('\n== R^2 = 2 c1^2 / (1 + c1) against [K^-1]_33 ==\n')
for (c1 in c(0.1, 0.2, 0.3)) {
  Kc <- equi(1, c1)
  cat(sprintf('c1 %.1f R2 %.4f 1 - 1/[Kinv]33 %.4f\n', c1,
    2 * c1^2 / (1 + c1), 1 - 1 / solve(Kc)[3, 3]))
}

cat('\n== Per-step costs, single step with gap 1, rho 0.7 ==\n')
phi <- 0.7
A2 <- ar1(c(0, 1), phi)
qf <- function(g) drop(t(g) %*% solve(A2, g))
cat(sprintf('on-on %.4f formula %.4f\n', qf(c(1, 1)) - 1, (1 - phi) / (1 + phi)))
cat(sprintf('on-off %.4f formula %.4f\n', qf(c(1, 0)) - 1, phi^2 / (1 - phi^2)))
cat(sprintf('off-on %.4f formula %.4f\n', qf(c(0, 1)), 1 / (1 - phi^2)))

cat('\n== Generalized ceiling, Hybrid path A, rho 0.7, c1 0.2 ==\n')
w <- hyb_w
g <- paths[['Hybrid A']][[2]]
A <- ar1(w, 0.7)
Ki <- solve(K)
Ai <- solve(A)
gen_lhs <- function(vtv, vpb, vbr, cbl) {
  V <- cbind(vtv, vpb, vbr)
  cbl^2 + sum(Ki * (t(V) %*% Ai %*% V))
}
z <- rep(0, 8)
one <- rep(1, 8)
lhs_both <- gen_lhs(0.3 * one, z, 0.3 * g, 0.3)
n <- 8
R <- diag(2 + 3 * n)
R[3:(2 + 3 * n), 3:(2 + 3 * n)] <- kronecker(K, A)
rv <- c(0.3 * one, z, 0.3 * g)
R[1, 3:(2 + 3 * n)] <- rv
R[3:(2 + 3 * n), 1] <- rv
R[1, 2] <- R[2, 1] <- 0.3
Mfull <- R[2:(2 + 3 * n), 2:(2 + 3 * n)]
schur <- drop(t(c(0.3, rv)) %*% solve(Mfull, c(0.3, rv)))
cat(sprintf('lhs cbm=cbmTV=cbl=0.3: %.4f full-matrix %.4f minEigR %.4f\n',
  lhs_both, schur, min_eig(R)))
br_only <- gen_lhs(z, z, 0.3 * g, 0)
tv_only <- gen_lhs(0.3 * one, z, z, 0)
both <- gen_lhs(0.3 * one, z, 0.3 * g, 0)
cat(sprintf('BR only %.4f TV only %.4f both %.4f cross %.4f\n',
  br_only, tv_only, both, both - br_only - tv_only))
cat(sprintf('Kinv offdiag %.5f formula %.5f\n', Ki[1, 3],
  -0.2 / ((1 - 0.2) * (1 + 0.4))))

cat('\n== CS ceiling versus number of visits, half on drug ==\n')
for (rho in c(0.7, 0.8)) {
  for (nv in c(4, 8, 16, 32)) {
    g <- c(rep(1, nv / 2), rep(0, nv / 2))
    cat(sprintf('rho %.1f n %2d closed %.4f full %.4f\n', rho, nv,
      cs_closed(nv, nv / 2, rho, 0.2, 0.1), ceil_full(cs_M(nv, rho, 0.2, 0.1), g)))
  }
}

cat('\n== CS minimum eigenvalue (1 - rho) - (c1 - cx), rho 0.8 ==\n')
for (nv in c(4, 8, 16, 32)) {
  cat(sprintf('n %2d minEig %.4f\n', nv, min_eig(cs_M(nv, 0.8, 0.2, 0.1))))
}

cat('\n== Separable AR(1) constant pattern: sqrt((1-R2)/n_eff) ==\n')
w <- hyb_w
A <- ar1(w, 0.7)
neff <- 1 + sum((1 - 0.7^diff(w)) / (1 + 0.7^diff(w)))
cat(sprintf('n_eff %.4f 1^T A^-1 1 %.4f ceiling %.4f full %.4f\n', neff,
  sum(solve(A, one)), sqrt((1 - 2 * 0.04 / 1.2) / neff),
  ceil_full(kronecker(K, A), one)))
