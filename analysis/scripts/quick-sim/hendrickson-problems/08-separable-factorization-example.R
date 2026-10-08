## Worked example of the separable factorization for Appendix A of
## docs/36-hendrickson-pd-and-carryover.md, at the published parameters,
## on the four Hybrid path-A visits of Appendix A.5 (weeks 10, 11, 12,
## 16: on, off, off, on).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/08-separable-factorization-example.R
w <- c(10, 11, 12, 16)
rho <- 0.8; c1 <- 0.2; cx <- 0.1; cbm <- 0.3
n <- length(w)
f4 <- function(x) sprintf('%.4f', x)
tex <- function(m, d = 4) {
  rows <- apply(m, 1, function(r) paste(sprintf(paste0('%.', d, 'f'), r),
                                        collapse = ' & '))
  paste0('\\begin{pmatrix}\n', paste(rows, collapse = ' \\\\\n'), '\n\\end{pmatrix}')
}
say <- function(...) cat(sprintf(...), '\n', sep = '')

## ---- factor part K ----
K <- (1 - c1) * diag(3) + c1
Kinv <- solve(K); LK <- t(chol(K))
say('K eigenvalues: %s', paste(f4(eigen(K)$values), collapse = ' '))
cat('K^{-1}\n', tex(Kinv), '\n', sep = '')
say('[K^-1]_33 = %s (closed form %s); R^2 = %s',
    f4(Kinv[3, 3]), f4((1 + c1) / ((1 - c1) * (1 + 2 * c1))), f4(1 - 1 / Kinv[3, 3]))
cat('L_K (lower Cholesky of K)\n', tex(LK), '\n', sep = '')

## ---- time part A ----
g <- abs(outer(w, w, '-'))
A <- rho^g
phi <- rho^diff(w)
say('gaps %s; phi_t = rho^gap: %s; sqrt(1 - phi^2): %s',
    paste(diff(w), collapse = ' '), paste(f4(phi), collapse = ' '),
    paste(f4(sqrt(1 - phi^2)), collapse = ' '))
say('A eigenvalues: %s', paste(f4(eigen(A)$values), collapse = ' '))
Ainv <- solve(A)
cat('A^{-1} (tridiagonal)\n', tex(Ainv), '\n', sep = '')
## Bidiagonal precision factor: A^{-1} = L'L.
L <- diag(c(1, 1 / sqrt(1 - phi^2)))
for (t in 2:n) L[t, t - 1] <- -phi[t - 1] / sqrt(1 - phi[t - 1]^2)
cat('L (bidiagonal, A^{-1} = L\'L)\n', tex(L), '\n', sep = '')
say('max |L\'L - A^{-1}| = %.1e', max(abs(t(L) %*% L - Ainv)))
cat('chol factor of A (lower, A = L_A L_A\')\n', tex(t(chol(A))), '\n', sep = '')

## ---- Kronecker identities ----
M <- kronecker(K, A)
say('max |(K x A)^{-1} - K^{-1} x A^{-1}| = %.1e',
    max(abs(solve(M) - kronecker(Kinv, Ainv))))
say('max |chol(K x A) - chol(K) x chol(A)| = %.1e',
    max(abs(chol(M) - kronecker(chol(K), chol(A)))))
say('det(K x A) = %.6e; det(K)^n det(A)^3 = %.6e', det(M), det(K)^n * det(A)^3)

## ---- ceilings via switch costs ----
step_cost <- function(u) {
  a <- u[-n]; b <- u[-1]
  (b - phi * a)^2 / (1 - phi^2)
}
pats <- list(no_carryover = c(1, 0, 0, 1), graded_t1 = c(1, 0.5, 0.25, 1),
             step_any = c(1, 1, 1, 1))
for (nm in names(pats)) {
  u <- pats[[nm]]
  sc <- step_cost(u)
  q_t <- u[1]^2 + sum(sc)
  R <- diag(2 + 3 * n); R[3:(2 + 3 * n), 3:(2 + 3 * n)] <- M
  idx <- (3 + 2 * n):(2 + 3 * n); R[1, idx] <- R[idx, 1] <- u
  v <- c(rep(0, 2 * n), u)
  full <- 1 / sqrt(drop(t(v) %*% solve(M, v)))
  say('%-13s u = (%s): u1^2 = %s, step costs %s, u\'A^-1u = %s; ceiling %s (full matrix %s)',
      nm, paste(u, collapse = ', '), f4(u[1]^2), paste(f4(sc), collapse = ' '),
      f4(q_t), f4(1 / sqrt(Kinv[3, 3] * q_t)), f4(full))
}

## ---- the current covar block, decomposed ----
Ka <- (1 - cx) * diag(3) + cx
Kb <- (c1 - cx) * (matrix(1, 3, 3) - diag(3))
Cc <- cx * A; diag(Cc) <- c1
Mc <- rbind(cbind(A, Cc, Cc), cbind(Cc, A, Cc), cbind(Cc, Cc, A))
say('covar: max |K_a x A + K_b x I - M| = %.1e',
    max(abs(kronecker(Ka, A) + kronecker(Kb, diag(n)) - Mc)))
say('K_b eigenvalues: %s', paste(f4(eigen(Kb)$values), collapse = ' '))
lam <- min(eigen(A)$values)
say('lambda_min(A) = %s; threshold (c1-cx)/(1-cx) = %s; predicted min eig %s; actual %s',
    f4(lam), f4((c1 - cx) / (1 - cx)),
    f4(min(lam * (1 - cx) - (c1 - cx), lam * (1 + 2 * cx) + 2 * (c1 - cx))),
    f4(min(eigen(Mc)$values)))

## ---- the published block, decomposed ----
Kp <- rho * diag(3) + cx * (matrix(1, 3, 3) - diag(3))
Ko <- (1 - rho) * diag(3) + (c1 - cx) * (matrix(1, 3, 3) - diag(3))
Acs <- matrix(rho, n, n); diag(Acs) <- 1
Cs <- matrix(cx, n, n); diag(Cs) <- c1
Mo <- rbind(cbind(Acs, Cs, Cs), cbind(Cs, Acs, Cs), cbind(Cs, Cs, Acs))
say('orig: max |K_o x I + K_p x J - M| = %.1e',
    max(abs(kronecker(Ko, diag(n)) + kronecker(Kp, matrix(1, n, n)) - Mo)))
say('K_o eigenvalues %s; K_p eigenvalues %s; min eig M = %s',
    paste(f4(eigen(Ko)$values), collapse = ' '),
    paste(f4(eigen(Kp)$values), collapse = ' '), f4(min(eigen(Mo)$values)))

## ---- factored draw: biomarker given response, graded u at t1/2 = 1 ----
u <- pats$graded_t1
Au <- drop(Ainv %*% u)
wts <- cbm * outer(Kinv[, 3], Au)       # rows: TV, PB, BR; cols: visits
cat('biomarker regression weights c_bm [K^-1]_{c3} (A^-1 u)_t\n', tex(wts), '\n', sep = '')
q <- Kinv[3, 3] * sum(u * Au)
say('residual variance 1 - c_bm^2 q = %s (q = %s, ceiling %s)',
    f4(1 - cbm^2 * q), f4(q), f4(1 / sqrt(q)))
