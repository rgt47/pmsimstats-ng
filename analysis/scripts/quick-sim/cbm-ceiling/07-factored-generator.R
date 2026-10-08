## Factored generator for the separable AR(1) covariance-moderation DGP
## (docs/35-separable-response-covariance.md), checked against the
## target correlation matrix.
##
## Target, for one path with n visits at weeks w:
##   response block  M = K (x) A,  K = (1 - c1) I3 + c1 J3,
##                   A_ts = rho^|w_t - w_s|
##   biomarker row   r = c_bm (e3 (x) u),  u the coupling pattern
## The generator never forms or factors the (2 + 3n) x (2 + 3n) matrix:
##   1. three independent AR(1) paths by recursion (bidiagonal factor)
##   2. mix them across factors with the 3 x 3 Cholesky factor of K
##   3. draw the biomarker given the response, using the closed-form
##      AR(1) precision; its residual variance is 1 - c_bm^2 q, and the
##      draw stops if that is not positive (the ceiling check).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/cbm-ceiling/07-factored-generator.R

## AR(1) paths on unequal gaps: X_1 ~ N(0,1), X_t = f_t X_{t-1} +
## sqrt(1 - f_t^2) e_t with f_t = rho^(w_t - w_{t-1}). Rows are draws.
ar1_paths <- function(n_draw, w, rho) {
  f <- rho^diff(w)
  Z <- matrix(rnorm(n_draw * length(w)), n_draw)
  X <- Z
  for (t in seq_along(w)[-1])
    X[, t] <- f[t - 1] * X[, t - 1] + sqrt(1 - f[t - 1]^2) * Z[, t]
  X
}

## A^{-1} u from the bidiagonal precision A^{-1} = L'L, without solve().
ar1_prec_times <- function(u, w, rho) {
  f <- rho^diff(w); n <- length(u)
  d <- c(1, 1 / sqrt(1 - f^2))                 # diag of L
  o <- -f / sqrt(1 - f^2)                      # sub-diagonal of L
  Lu <- d * u + c(0, o * u[-n])                # L u
  ## L' (L u): transpose of a lower bidiagonal is upper bidiagonal
  d * Lu + c(o * Lu[-1], 0)
}

factored_draw <- function(n_draw, w, rho, c1, cbm, u) {
  n <- length(w)
  K <- (1 - c1) * diag(3) + c1
  LK <- t(chol(K))                             # 3 x 3, K = LK LK'
  P <- lapply(1:3, function(j) ar1_paths(n_draw, w, rho))
  ## factor c = sum_j LK[c, j] * P_j   (n_draw x n each)
  Fct <- lapply(1:3, function(c) Reduce(`+`, lapply(1:3, function(j)
    LK[c, j] * P[[j]])))
  Kinv <- solve(K)                             # 3 x 3 only
  Au <- ar1_prec_times(u, w, rho)              # A^{-1} u, closed form
  q <- Kinv[3, 3] * sum(u * Au)                # u'A^{-1}u [K^{-1}]_33
  resid_var <- 1 - cbm^2 * q
  if (resid_var <= 0)
    stop(sprintf('c_bm = %.3f exceeds the ceiling %.3f', cbm, 1 / sqrt(q)))
  ## E[B | X] = r~' M^{-1} X = cbm * sum_c Kinv[c,3] * (Fct_c %*% Au)
  mean_B <- cbm * Reduce(`+`, lapply(1:3, function(c)
    Kinv[c, 3] * drop(Fct[[c]] %*% Au)))
  B <- mean_B + sqrt(resid_var) * rnorm(n_draw)
  BL <- rnorm(n_draw)
  cbind(B, BL, Fct[[1]], Fct[[2]], Fct[[3]])
}

target_R <- function(w, rho, c1, cbm, u) {
  n <- length(w)
  A <- rho^abs(outer(w, w, '-'))
  K <- (1 - c1) * diag(3) + c1
  R <- diag(2 + 3 * n)
  R[3:(2 + 3 * n), 3:(2 + 3 * n)] <- kronecker(K, A)
  idx <- (3 + 2 * n):(2 + 3 * n)
  R[1, idx] <- R[idx, 1] <- cbm * u
  R
}

set.seed(20261002)
w <- c(4, 8, 9, 10, 11, 12, 16, 20)            # Hybrid
u <- c(1, 1, 1, 1, 0.25, 0.0625, 1, 0.0625)    # path A, graded, t1/2 = 0.5
rho <- 0.8; c1 <- 0.2; cbm <- 0.3
X <- factored_draw(400000, w, rho, c1, cbm, u)
R_hat <- cor(X)
R_tgt <- target_R(w, rho, c1, cbm, u)
cat(sprintf('max |empirical - target| correlation over 26 x 26: %.4f\n',
            max(abs(R_hat - R_tgt))))
cat(sprintf('Monte Carlo SE of a correlation near 0 with 400k draws: %.4f\n',
            1 / sqrt(400000)))
idx <- (3 + 16):(2 + 24)
cat('biomarker-BR correlations, empirical vs target:\n')
print(rbind(empirical = round(R_hat[1, idx], 4), target = round(R_tgt[1, idx], 4)))

cat('\nceiling check:\n')
q <- (solve((1 - c1) * diag(3) + c1)[3, 3]) * sum(u * ar1_prec_times(u, w, rho))
v <- c(rep(0, 16), u)
M <- R_tgt[3:26, 3:26]
cat(sprintf('closed form %.6f  vs  full-matrix %.6f\n', 1 / sqrt(q),
            1 / sqrt(drop(t(v) %*% solve(M, v)))))
cat(sprintf('A^{-1}u closed form vs solve(): max diff %.2e\n',
            max(abs(ar1_prec_times(u, w, rho) -
                    solve(rho^abs(outer(w, w, '-')), u)))))
## Benchmark: the same number of draws from the full matrix, so the
## factored generator's deviations can be compared with sampling noise.
dev <- function(draws) max(abs(cor(draws) - R_tgt))
rep_dev <- t(sapply(1:5, function(s) {
  set.seed(100 + s)
  c(factored = dev(factored_draw(400000, w, rho, c1, cbm, u)),
    full_mvrnorm = dev(MASS::mvrnorm(400000, rep(0, 26), R_tgt)))
}))
cat('\nmax |empirical - target| over 5 seeds, 400k draws each:\n')
print(round(rep_dev, 4))
res <- tryCatch(factored_draw(10, w, rho, c1, 0.6, u),
                error = function(e) conditionMessage(e))
cat('draw at c_bm = 0.6:', res, '\n')
