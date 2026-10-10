## Copied 2026-10-08 from the review scratch script inv_d_check_ol.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-open-label-ceiling.R
## OL constant-pattern ceiling on an equally spaced 2.5-week schedule
## (docs/44, Step 4: 0.893 under CS; 0.576 under AR(1) with the decaying
## cross-factor block, which is docs/36 configuration B (step) or F
## (graded); the constant published block with AR(1) is not one of the
## docs/36 configurations and gives 0.535; rho 0.8).
w <- seq(2.5, 20, 2.5)
n <- 8
lag <- abs(outer(w, w, '-'))
mk <- function(A, C) rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
ceil1 <- function(M) { u <- c(rep(0, 16), rep(1, 8)); 1 / sqrt(drop(crossprod(u, solve(M, u)))) }
rho <- 0.8
Acs <- matrix(rho, n, n); diag(Acs) <- 1
Cc <- matrix(0.1, n, n); diag(Cc) <- 0.2
Cd <- 0.1 * rho^lag; diag(Cd) <- 0.2
cat(sprintf('CS (A) %.4f  AR1+const cross %.4f  AR1+decay cross (B/F) %.4f  separable (C) %.4f\n',
            ceil1(mk(Acs, Cc)), ceil1(mk(rho^lag, Cc)), ceil1(mk(rho^lag, Cd)),
            ceil1(mk(rho^lag, 0.2 * rho^lag))))
## Hybrid, no carryover, AR(1) with the constant (published) cross block
## versus the decaying block, rho 0.8.
wh <- c(4, 8, 9, 10, 11, 12, 16, 20)
lh <- abs(outer(wh, wh, '-'))
Cch <- matrix(0.1, 8, 8); diag(Cch) <- 0.2
Cdh <- 0.1 * rho^lh; diag(Cdh) <- 0.2
pats <- list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1), c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1))
cg <- function(M, g) { u <- c(rep(0, 16), g); 1 / sqrt(drop(crossprod(u, solve(M, u)))) }
cat(sprintf('Hybrid t_half 0, rho 0.8: AR1+const %.4f  AR1+decay %.4f\n',
            min(sapply(pats, cg, M = mk(rho^lh, Cch))), min(sapply(pats, cg, M = mk(rho^lh, Cdh)))))
