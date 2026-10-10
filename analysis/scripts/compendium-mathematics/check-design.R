## Checks for the bp fragment, section 2: contrast variances by design.
qa <- function(a, R) drop(t(a) %*% R %*% a)
wts <- function(on) ifelse(on == 1, 1 / sum(on), -1 / sum(on == 0))
ar1 <- function(w, rho) rho^abs(outer(w, w, '-'))

## Paper 08 mini: S_m, C_pq, V_CS, V_AR1 against the quadratic form.
S <- function(m, r) m * (1 + r) / (1 - r) - 2 * r * (1 - r^m) / (1 - r)^2
Cc <- function(p, q, r) r * (1 - r^p) * (1 - r^q) / (1 - r)^2
V_cs <- function(p, q, r) (1 - r) * (1 / p + 1 / q)
V_ar <- function(p, q, r) S(p, r) / p^2 + S(q, r) / q^2 - 2 * Cc(p, q, r) / (p * q)
for (pq in list(c(1, 1), c(4, 4), c(3, 5), c(8, 8))) {
  p <- pq[1]; q <- pq[2]; on <- c(rep(1, p), rep(0, q))
  Rcs <- matrix(0.5, p + q, p + q); diag(Rcs) <- 1
  cat(sprintf('p=%d q=%d  V_CS %.4f (qf %.4f)  V_AR1 %.4f (qf %.4f)  ratio %.3f\n',
              p, q, V_cs(p, q, .5), qa(wts(on), Rcs), V_ar(p, q, .5),
              qa(wts(on), ar1(seq_len(p + q), .5)), V_ar(p, q, .5) / V_cs(p, q, .5)))
}
## Toeplitz cross weight min(l, p, q, p+q-l) for a non-AR(1) kernel.
rk <- function(k) ifelse(k == 0, 1, 0.6 * 0.7^(k - 1))
for (pq in list(c(2, 6), c(3, 5), c(4, 4))) {
  p <- pq[1]; q <- pq[2]; Tt <- p + q
  R <- outer(1:Tt, 1:Tt, function(i, j) rk(abs(i - j)))
  l <- 1:(Tt - 1)
  scr <- sum(pmin(l, p, q, Tt - l) * rk(l))
  son <- p + 2 * sum((p - seq_len(p - 1)) * rk(seq_len(p - 1)))
  soff <- q + 2 * sum((q - seq_len(q - 1)) * rk(seq_len(q - 1)))
  cat(sprintf('Toeplitz p=%d q=%d block-sum %.6f quadratic %.6f\n', p, q,
              son / p^2 + soff / q^2 - 2 * scr / (p * q),
              qa(c(rep(1 / p, p), rep(-1 / q, q)), R)))
}
## Random slope adds tau^2 (tbar_on - tbar_off)^2.
w <- 1:8; on <- c(1, 1, 1, 1, 0, 0, 0, 0); a <- wts(on)
cat(sprintf('random slope term: qf %.4f vs formula %.4f\n',
            qa(a, outer(w, w)), (mean(w[on == 1]) - mean(w[on == 0]))^2))

## docs/39: time part a'Aa under AR(1), rho = 0.8 per week.
pat <- function(s) as.numeric(strsplit(s, '')[[1]])
for (s in c('11110000', '11001100', '10101010')) {
  cat(sprintf('pattern %s at 2.5 wk: a\'Aa = %.3f\n', s,
              qa(wts(pat(s)), ar1(2.5 * (1:8), 0.8))))
}
hyb_w <- c(4, 8, 9, 10, 11, 12, 16, 20)
hyb <- list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1),
            c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1))
cat('Hybrid a\'Aa:', round(sapply(hyb, function(o) qa(wts(o), ar1(hyb_w, 0.8))), 3), '\n')
cat('ABAB at Hybrid times:', round(qa(wts(pat('11001100')), ar1(hyb_w, 0.8)), 3), '\n')

## Large-sample alternating approximation and the optimum of f(d).
for (d in c(1, 2.5, 4)) {
  cat(sprintf('d=%.1f exact %.3f approx %.3f\n', d,
              qa(wts(pat('10101010')), ar1(d * (1:8), 0.8)),
              0.5 * (1 - 0.8^d) / (1 + 0.8^d)))
}
f <- function(d, th, r = 0.8) (1 - 2^(-d / th))^2 * (1 + r^d) / (1 - r^d)
cat('argmax f(d):', round(sapply(c(0.25, 0.5, 1, 2), function(th)
  optimize(f, c(0.01, 20), th = th, maximum = TRUE)$maximum), 2), '\n')

## CS components of Var(Delta), docs/39 Section 3 (A.1 of paper 01).
cs_var <- function(on, e, rho = 0.8, sTV = 10, sBR = 8, sPB = 10,
                   c1 = 0.2, cx = 0.1) {
  a <- wts(on); m1 <- sum(a); m2 <- sum(a^2); me <- sum(a * e)
  mee <- sum(a^2 * e^2); m2e <- sum(a^2 * e)
  bal <- (sTV^2 + sBR^2) * (rho * m1^2 + (1 - rho) * m2)
  pb <- sPB^2 * (rho * me^2 + (1 - rho) * mee)
  x1 <- 2 * sTV * sBR * (cx * m1^2 + (c1 - cx) * m2)
  x2 <- 2 * (sTV + sBR) * sPB * (cx * m1 * me + (c1 - cx) * m2e)
  c(balance = bal, total = bal + pb + x1 + x2)
}
eh <- c(1, 1, .5, .5, .5, .5, .5, .5)
cat('Hybrid path 1:', round(cs_var(hyb[[1]], eh), 2), '\n')
cat('Hybrid path 3:', round(cs_var(hyb[[3]], eh), 2), '\n')
cat('Hybrid path 1 blinded:', round(cs_var(hyb[[1]], rep(.5, 8)), 2), '\n')
cat('Blinded 4/4:', round(cs_var(pat('10101010'), rep(.5, 8)), 2), '\n')

## Paper 12: independence approximation versus exact a'Aa (rho 0.7).
for (o in list(pat('11110000'), hyb[[1]], hyb[[3]])) {
  ww <- if (length(unique(o)) && identical(o, pat('11110000'))) 2.5 * (1:8) else hyb_w
  R <- ar1(ww, 0.7); a <- wts(o)
  ion <- o == 1
  v_ind <- sum(R[ion, ion]) / sum(ion)^2 + sum(R[!ion, !ion]) / sum(!ion)^2
  cat(sprintf('v independence %.3f exact %.3f\n', v_ind, qa(a, R)))
}
