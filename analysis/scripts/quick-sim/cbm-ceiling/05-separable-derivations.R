## Numerical checks for the derivations in
## docs/35-separable-response-covariance.md.
##
##   D1  covar  M = K_a (x) A + K_b (x) I, K_b indefinite; PD iff
##       lambda_min(A) (1 - cx) > c1 - cx  (common rho).
##   D2  orig   M = K_o (x) I + K_p (x) J  (occasion + person terms).
##   D3  separable  [K^-1]_33 = (1 + c1) / ((1 - c1)(1 + 2 c1)).
##   D4  separable  u'A^-1 u for 0/1 u = u_1 + sum over steps of
##       on-on (1-f)/(1+f), on-off f^2/(1-f^2), off-on 1/(1-f^2).
##   D5  orig is separable iff cx = c1 * rho.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/cbm-ceiling/05-separable-derivations.R

library(data.table)

ar1 <- function(w, rho) rho^abs(outer(w, w, '-'))
equi <- function(d, off) {
  m <- matrix(off, 3, 3)
  diag(m) <- d
  m
}
min_eig <- function(m) min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)

covar_M <- function(w, rho, c1, cx) {
  A <- ar1(w, rho)
  C <- cx * A
  diag(C) <- c1
  rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
}
orig_M <- function(n, rho, c1, cx) {
  A <- matrix(rho, n, n); diag(A) <- 1
  C <- matrix(cx, n, n); diag(C) <- c1
  rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
}

schedules <- list(
  'weekly x8'    = 0:7,
  'half-week x8' = seq(0, 3.5, 0.5),
  'quarter-week x8' = seq(0, 1.75, 0.25),
  'Hybrid'       = c(4, 8, 9, 10, 11, 12, 16, 20),
  'OL+BDC'       = c(4, 8, 12, 16, 17, 18, 19, 20),
  'weekly x16'   = 0:15
)
pars <- list(c(0.2, 0.1), c(0.1, 0.05), c(0.3, 0.2))

cat('== D1: covar decomposition and exact PD condition (rho = 0.7) ==\n')
d1 <- rbindlist(lapply(names(schedules), function(nm) {
  w <- schedules[[nm]]
  rbindlist(lapply(pars, function(p) {
    c1 <- p[1]; cx <- p[2]
    A <- ar1(w, 0.7); n <- length(w)
    Ka <- equi(1, cx); Kb <- equi(0, c1 - cx)
    M_dec <- kronecker(Ka, A) + kronecker(Kb, diag(n))
    M_dir <- covar_M(w, 0.7, c1, cx)
    data.table(schedule = nm, c1 = c1, cx = cx,
      recon_err = max(abs(M_dec - M_dir)),
      lambda_min_A = min_eig(A),
      threshold = (c1 - cx) / (1 - cx),
      predicted_PD = min_eig(A) * (1 - cx) > (c1 - cx),
      min_eig_M = min_eig(M_dir),
      pred_min_eig = min(min_eig(A) * (1 - cx) - (c1 - cx),
                         min_eig(A) * (1 + 2 * cx) + 2 * (c1 - cx)))
  }))
}))
print(d1, digits = 3)
cat(sprintf('max reconstruction error %.1e; max |min_eig - predicted| %.1e\n',
  max(d1$recon_err), max(abs(d1$min_eig_M - d1$pred_min_eig))))

cat('\n== D2: orig as occasion + person terms (rho = 0.7) ==\n')
d2 <- rbindlist(lapply(c(4, 8, 16, 32), function(n) {
  rbindlist(lapply(pars, function(p) {
    c1 <- p[1]; cx <- p[2]
    Ko <- equi(1 - 0.7, c1 - cx); Kp <- equi(0.7, cx)
    M_dec <- kronecker(Ko, diag(n)) + kronecker(Kp, matrix(1, n, n))
    data.table(n = n, c1 = c1, cx = cx,
      recon_err = max(abs(M_dec - orig_M(n, 0.7, c1, cx))),
      min_eig_Ko = min_eig(Ko), min_eig_Kp = min_eig(Kp),
      person_cross_cor = cx / 0.7,
      occasion_cross_cor = (c1 - cx) / (1 - 0.7),
      min_eig_M = min_eig(orig_M(n, 0.7, c1, cx)))
  }))
}))
print(d2, digits = 3)

cat('\n== D3: [K^-1]_33 closed form ==\n')
d3 <- rbindlist(lapply(c(0.05, 0.1, 0.2, 0.3), function(c1) {
  K <- equi(1, c1)
  data.table(c1 = c1, numeric = solve(K)[3, 3],
             closed = (1 + c1) / ((1 - c1) * (1 + 2 * c1)),
             R2 = 1 - 1 / solve(K)[3, 3])
}))
print(d3, digits = 4)

cat('\n== D4: switch-cost form of u\'A^-1 u (rho = 0.7) ==\n')
switch_cost <- function(w, u, rho) {
  f <- rho^diff(w)
  a <- u[-length(u)]; b <- u[-1]
  u[1]^2 + sum(ifelse(a == 1 & b == 1, (1 - f) / (1 + f),
               ifelse(a == 1 & b == 0, f^2 / (1 - f^2),
               ifelse(a == 0 & b == 1, 1 / (1 - f^2), 0))))
}
cases <- list(
  list('CO A', seq(2.5, 20, 2.5), c(1,1,1,1,0,0,0,0)),
  list('CO B', seq(2.5, 20, 2.5), c(0,0,0,0,1,1,1,1)),
  list('Hybrid A', c(4,8,9,10,11,12,16,20), c(1,1,1,1,0,0,1,0)),
  list('Hybrid C', c(4,8,9,10,11,12,16,20), c(1,1,1,0,0,0,1,0)),
  list('OL+BDC A', c(4,8,12,16,17,18,19,20), c(1,1,1,1,1,1,0,0)),
  list('Multi 16wk', 1:16, rep(c(1,1,0,0), 4)))
d4 <- rbindlist(lapply(cases, function(cs) {
  w <- cs[[2]]; u <- cs[[3]]
  q_mat <- drop(t(u) %*% solve(ar1(w, 0.7), u))
  q_sw <- switch_cost(w, u, 0.7)
  kfac <- (1 + 0.2) / ((1 - 0.2) * (1 + 2 * 0.2))
  data.table(path = cs[[1]], quad_matrix = q_mat, quad_switch = q_sw,
             ceiling = 1 / sqrt(kfac * q_sw))
}))
print(d4, digits = 4)
cat(sprintf('max |matrix - switch-cost| %.1e\n',
            max(abs(d4$quad_matrix - d4$quad_switch))))

cat('\n== D5: orig separable iff cx = c1 * rho ==\n')
n <- 6
for (cx in c(0.1, 0.2 * 0.7)) {
  M <- orig_M(n, 0.7, 0.2, cx)
  K <- equi(1, 0.2)
  A <- matrix(0.7, n, n); diag(A) <- 1
  cat(sprintf('cx = %.3f: max |M - K (x) A_cs| = %.3f\n', cx,
              max(abs(M - kronecker(K, A)))))
}

cat('\n== Long-lag within-factor correlation ==\n')
print(data.table(gap_weeks = c(1, 4, 8, 16),
  ar1 = 0.7^c(1, 4, 8, 16), orig_cs = 0.7), digits = 3)
