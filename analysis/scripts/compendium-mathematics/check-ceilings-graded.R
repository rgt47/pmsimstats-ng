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
graded <- function(w, on, th) {
  g <- numeric(length(w))
  last_on <- NA
  for (t in seq_along(w)) {
    if (on[t] == 1) {
      g[t] <- 1
      last_on <- w[t]
    } else if (!is.na(last_on) && th > 0) {
      g[t] <- 0.5^((w[t] - last_on) / th)
    }
  }
  g
}
hyb_w <- c(4, 8, 9, 10, 11, 12, 16, 20)
hp <- list(A = c(1, 1, 1, 1, 0, 0, 1, 0), B = c(1, 1, 1, 1, 0, 0, 0, 1),
  C = c(1, 1, 1, 0, 0, 0, 1, 0), D = c(1, 1, 1, 0, 0, 0, 0, 1))
M <- cs_M(8, 0.8, 0.2, 0.1)
for (th in c(0, 0.1, 0.2, 1, 2)) {
  v <- sapply(hp, function(on) ceil_full(M, graded(hyb_w, on, th)))
  cat(sprintf('CS graded Hybrid rho 0.8 t1/2 %.1f: %s min %.4f\n', th,
    paste(sprintf('%.4f', v), collapse = ' '), min(v)))
}
