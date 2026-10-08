## Worked-example matrices for the appendix of
## docs/36-hendrickson-pd-and-carryover.md, at the published parameters.
## Four Hybrid path-A visits at weeks 10, 11, 12, 16: on, off, off, on;
## time since discontinuation (0, 1, 2, 0).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/06-worked-example-matrices.R
w <- c(10, 11, 12, 16)
on <- c(1, 0, 0, 1)
tsd <- c(0, 1, 2, 0)
rho <- 0.8; c1 <- 0.2; cx <- 0.1; cbm <- 0.3
g <- abs(outer(w, w, '-'))

fmt <- function(m) {
  rows <- apply(m, 1, function(r) paste(sprintf('%.4f', r), collapse = ' & '))
  paste0('\\begin{pmatrix}\n', paste(rows, collapse = ' \\\\\n'),
         '\n\\end{pmatrix}')
}
cat('gaps\n'); print(g)

A_cs <- matrix(rho, 4, 4); diag(A_cs) <- 1
A_ar <- rho^g
cat('\nA_cs\n', fmt(A_cs), '\n\nA_ar1\n', fmt(A_ar), '\n', sep = '')

C_const <- matrix(cx, 4, 4); diag(C_const) <- c1
C_decay <- cx * rho^g; diag(C_decay) <- c1
C_sep <- c1 * A_ar
cat('\nC constant (published)\n', fmt(C_const), '\n', sep = '')
cat('\nC decaying (covar)\n', fmt(C_decay), '\n', sep = '')
cat('\nC separable\n', fmt(C_sep), '\n', sep = '')

## Coupling vectors. The step rule tests the recursive carryover mean,
## which is nonzero at both off-drug visits for any t_half > 0.
for (th in c(0, 0.1, 1)) {
  phi <- if (th > 0) 0.5^(tsd / th) else rep(0, 4)
  phi[on == 1] <- 0
  graded <- cbm * ifelse(on == 1, 1, phi)
  step <- cbm * ifelse(on == 1, 1, as.numeric(th > 0))
  cat(sprintf('\nt_half = %.1f\n  step   r = %s\n  graded r = %s\n', th,
              paste(sprintf('%.4f', step), collapse = ' '),
              paste(sprintf('%.6f', graded), collapse = ' ')))
}

## Smallest eigenvalue of each response block for this 4-visit example.
blk <- function(A, C) rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
me <- function(m) min(eigen(m, TRUE, TRUE)$values)
cat(sprintf('\nmin eigenvalue: published %.4f, covar %.4f, separable %.4f\n',
            me(blk(A_cs, C_const)), me(blk(A_ar, C_decay)), me(blk(A_ar, C_sep))))
