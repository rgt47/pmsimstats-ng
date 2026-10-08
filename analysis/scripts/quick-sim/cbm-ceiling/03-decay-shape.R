## Does the shape of carryover decay change the c_bm ceiling?
##
## Decay shape enters only through the off-drug entries of the
## coupling vector v; the response block M does not depend on it.
## Shapes are the Weibull family used in paper 02. Under the graded
## coupling the shape sets the off-drug values; under the step
## coupling only its zero pattern matters. Weibull decay never reaches
## zero analytically, but fast shapes underflow to exactly 0 in double
## precision at long times, which the step coupling detects.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/cbm-ceiling/03-decay-shape.R
## Output: analysis/data/quick-sim/cbm-ceiling/decay-shape.csv
## Whitepaper: docs/34-cbm-feasible-range.md, Section 5.6

suppressMessages(pkgload::load_all('.', quiet = TRUE))
library(data.table)
out_dir <- 'analysis/data/quick-sim/cbm-ceiling'

## Weibull family of paper 02 (k = 0.25, 0.5, 2, 4), with k = 1 the
## exponential. lambda_w = (ln 2)^(1/k) / t_half puts phi = 0.5 at
## t_half for every k. k < 1: heavier, slower tail; k > 1: faster
## washout.
shape_k <- c(weibull_0.25 = 0.25, weibull_0.5 = 0.5, exponential = 1,
             weibull_2 = 2, weibull_4 = 4)

decay <- function(tsd, t_half, shape) {
  if (t_half == 0) return(rep(0, length(tsd)))
  k <- shape_k[[shape]]
  exp(-((log(2))^(1 / k) / t_half * tsd)^k)
}

pattern <- function(tod, tsd, t_half, shape, coupling) {
  on <- tod > 0
  phi <- decay(tsd, t_half, shape)
  phi[on | tsd <= 0] <- 0
  off <- if (coupling == 'graded') phi else as.numeric(phi > 0)
  ifelse(on, 1, off)
}

ceiling_of <- function(w, u, structure, rho = 0.7, c1 = 0.2, cx = 0.1) {
  n <- length(w)
  gap <- abs(outer(w, w, '-'))
  if (structure == 'ar1') {
    A <- rho^gap
    C <- cx * rho^gap
  } else {
    A <- matrix(rho, n, n)
    diag(A) <- 1
    C <- matrix(cx, n, n)
  }
  diag(C) <- c1
  M <- rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
  v <- c(rep(0, 2 * n), u)
  1 / sqrt(drop(t(v) %*% solve(M, v)))
}

mk <- function(tp, od) {
  tp <- as.data.table(tp)
  list(w = cumulative(tp$t_wk), tod = tp$tod, tsd = tp$tsd)
}

designs <- list(
  CO = list(timepoints = cumulative(rep(2.5, 8)),
            ondrug = list(c(1, 1, 1, 1, 0, 0, 0, 0),
                          c(0, 0, 0, 0, 1, 1, 1, 1))),
  Hybrid = list(timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
                ondrug = list(c(1, 1, 1, 1, 0, 0, 1, 0),
                              c(1, 1, 1, 1, 0, 0, 0, 1),
                              c(1, 1, 1, 0, 0, 0, 1, 0),
                              c(1, 1, 1, 0, 0, 0, 0, 1))),
  `OL+BDC` = list(timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
                  ondrug = list(c(1, 1, 1, 1, 1, 1, 0, 0),
                                c(1, 1, 1, 1, 1, 0, 0, 0))),
  `Multi-cycle` = list(timepoints = 1:16,
                       ondrug = list(rep(c(1, 1, 0, 0), 4)))
)

paths <- lapply(designs, function(d) {
  td <- buildtrialdesign(
    name_longform = 'x', name_shortform = 'x',
    timepoints = d$timepoints,
    timeptnames = paste0('V', seq_along(d$timepoints)),
    expectancies = rep(1, length(d$timepoints)),
    ondrug = setNames(d$ondrug, paste0('p', seq_along(d$ondrug))))
  lapply(td$trialpaths, mk)
})

shapes <- names(shape_k)
grid <- CJ(design = names(designs), shape = shapes,
           t_half = c(0.5, 1, 2),
           cfg = c('ar1/graded', 'cs/step', 'ar1/step', 'cs/graded'))
grid[, c('structure', 'coupling') := tstrsplit(cfg, '/')]
grid[, c_star := min(sapply(paths[[design]], function(p)
  ceiling_of(p$w, pattern(p$tod, p$tsd, t_half, shape, coupling),
             structure))),
  by = .(design, shape, t_half, cfg)]

base <- CJ(design = names(designs), cfg = unique(grid$cfg))
base[, c('structure', 'coupling') := tstrsplit(cfg, '/')]
base[, c_star0 := min(sapply(paths[[design]], function(p)
  ceiling_of(p$w, pattern(p$tod, p$tsd, 0, 'exponential', coupling),
             structure))), by = .(design, cfg)]

cat('Ceiling at t1/2 = 0 (no carryover; shape irrelevant):\n')
print(dcast(base, design ~ cfg, value.var = 'c_star0'), digits = 3)

cat('\nHalf-life check: value of each shape at t = t1/2\n')
print(sapply(shapes, function(s) decay(1, 1, s)), digits = 3)

cat('\nResidual fraction at 1, 2, 4, 7 weeks, t1/2 = 0.5\n')
print(sapply(shapes, function(s) decay(c(1, 2, 4, 7), 0.5, s)),
      digits = 3)

for (cf in c('ar1/graded', 'cs/step')) {
  cat(sprintf('\n== %s ==\n', cf))
  print(dcast(grid[cfg == cf], design + t_half ~ shape,
              value.var = 'c_star')[, c('design', 't_half', shapes),
                                    with = FALSE], digits = 3)
}

fwrite(rbind(grid, base[, .(design, shape = 'none', t_half = 0, cfg,
                            structure, coupling, c_star = c_star0)],
             fill = TRUE),
       file.path(out_dir, 'decay-shape.csv'))
