## Copied 2026-10-08 from the review scratch script inv_c_check3.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-crossover-binding-path.R
# Fragment c: which CO path binds the ceiling, by half-life.
covar_M <- function(K, w, rho, c1, cx) {
  A <- rho^abs(outer(w, w, '-'))
  C <- cx * A; diag(C) <- c1
  kronecker(diag(K), A) + kronecker(matrix(1, K, K) - diag(K), C)
}
sep_M <- function(K, w, rho, c1) {
  Kc <- matrix(c1, K, K); diag(Kc) <- 1
  kronecker(Kc, rho^abs(outer(w, w, '-')))
}
graded <- function(w, on, th) {
  g <- numeric(length(w)); last <- NA
  for (k in seq_along(w)) {
    if (on[k] == 1) { g[k] <- 1; last <- w[k] }
    else if (!is.na(last) && th > 0) g[k] <- 2^(-(w[k] - last) / th)
  }
  g
}
ceil <- function(M, g, K) {
  u <- c(rep(0, (K - 1) * length(g)), g)
  1 / sqrt(drop(t(u) %*% solve(M, u)))
}
co_w <- seq(2.5, 20, 2.5)
pa <- c(1, 1, 1, 1, 0, 0, 0, 0); pbp <- c(0, 0, 0, 0, 1, 1, 1, 1)
for (th in c(0, 0.5, 1, 2)) {
  M <- covar_M(3, co_w, 0.7, 0.2, 0.1)
  S <- sep_M(3, co_w, 0.8, 0.2)
  cat(sprintf('t1/2 %.1f covar rho0.7: A %.4f B %.4f | separable rho0.8: A %.4f B %.4f\n',
    th, ceil(M, graded(co_w, pa, th), 3), ceil(M, graded(co_w, pbp, th), 3),
    ceil(S, graded(co_w, pa, th), 3), ceil(S, graded(co_w, pbp, th), 3)))
}
