## Numerical checks for the docs/35 review.
w <- c(0, 4, 8, 9, 10, 11, 12, 16, 20)[-1]   # Hybrid post-baseline weeks
n <- length(w)
## 1. orig at RH's published rho = 0.8 versus the doc's 0.7.
for (rho in c(0.7, 0.8)) {
  c1 <- 0.2; cx <- 0.1
  Ko <- (1 - rho) * diag(3) + (c1 - cx) * (matrix(1, 3, 3) - diag(3))
  Kp <- rho * diag(3) + cx * (matrix(1, 3, 3) - diag(3))
  M <- kronecker(Ko, diag(n)) + kronecker(Kp, matrix(1, n, n))
  cat(sprintf('rho %.1f: min eig orig %.4f; lambda_min(Ko) %.4f; trait cor %.3f; shock cor %.3f\n',
              rho, min(eigen(M, TRUE, TRUE)$values), min(eigen(Ko, TRUE, TRUE)$values),
              cx / rho, (c1 - cx) / (1 - rho)))
}
## 2. Generalized ceiling with couplings to several components and to BL,
## under separability: r' M^-1 r = sum_cc' [K^-1]_cc' v_c' A^-1 v_c' (+ c_bl^2).
rho <- 0.7; c1 <- 0.2
A <- rho^abs(outer(w, w, '-'))
K <- (1 - c1) * diag(3) + c1 * matrix(1, 3, 3)
M <- kronecker(K, A)
on <- c(1, 1, 1, 1, 0, 0, 1, 0)                 # Hybrid path A
v_tv <- rep(0.3, n); v_pb <- rep(0, n); v_br <- 0.3 * on
r <- c(v_tv, v_pb, v_br)
direct <- drop(t(r) %*% solve(M) %*% r)
Ki <- solve(K); Ai <- solve(A); V <- list(v_tv, v_pb, v_br)
formula <- sum(sapply(1:3, function(a) sapply(1:3, function(b)
  Ki[a, b] * drop(t(V[[a]]) %*% Ai %*% V[[b]]))))
cat(sprintf('generalized quadratic form: direct %.10f, formula %.10f\n', direct, formula))
## Full matrix with BL coupling 0.3: PD iff c_bl^2 + r'M^-1 r < 1.
c_bl <- 0.3
R <- diag(2 + 3 * n)
R[3:(2 + 3 * n), 3:(2 + 3 * n)] <- M
R[1, 3:(2 + 3 * n)] <- R[3:(2 + 3 * n), 1] <- r
R[1, 2] <- R[2, 1] <- c_bl
cat(sprintf('c_bl^2 + quad = %.4f (< 1 needed); min eig of R %.4f\n',
            c_bl^2 + direct, min(eigen(R, TRUE, TRUE)$values)))
## 3. Same-sign couplings to correlated components are cheaper: compare
## BR-only quadratic form with BR plus TV.
q_br <- Ki[3, 3] * drop(t(v_br) %*% Ai %*% v_br)
cat(sprintf('BR only %.4f; TV only %.4f; both %.4f; cross term %.4f\n', q_br,
            Ki[1, 1] * drop(t(v_tv) %*% Ai %*% v_tv), direct,
            2 * Ki[1, 3] * drop(t(v_tv) %*% Ai %*% v_br)))
## 4. Summed-response correlation when the PB SD varies by phase:
## correlation factor (s_t' K s_s) / sqrt(s_t' K s_t s_s' K s_s) at OL vs blinded.
s_ol <- c(10, 10, 8); s_bl <- c(10, 5, 8)
f <- drop(t(s_ol) %*% K %*% s_bl) / sqrt(drop(t(s_ol) %*% K %*% s_ol) * drop(t(s_bl) %*% K %*% s_bl))
cat(sprintf('summed-response correlation factor, open label vs blinded: %.4f\n', f))
