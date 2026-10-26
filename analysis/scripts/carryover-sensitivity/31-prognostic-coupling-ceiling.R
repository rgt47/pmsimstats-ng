## analysis/scripts/carryover-sensitivity/31-prognostic-coupling-ceiling.R
##
## Exact feasibility ceilings for a prognostic coupling (c.bm.pb,
## c.bm.tv or c.bm.bl) alongside the predictive coupling c.bm, by the
## Schur-complement argument of the compendium (sec. "The joint
## correlation matrix and the exact ceiling", eq. ce-ceiling) extended
## to a biomarker row with more than one pattern.
##
## With bm removed, the rest of the correlation matrix R_0 does not
## depend on any biomarker coupling. The bm row is
##     r = c_bm * g + c_pb * p + c_tv * q + c_bl * b,
## where g, p, q, b are the unit-coupling patterns the DGP writes into
## that row (g: BR on drug, decaying off drug; p: PB at every
## timepoint; q: TV at every timepoint; b: BL). R is PD iff R_0 is PD
## and  r' R_0^{-1} r < 1.  That is the quadratic form
##     Q(c) = c' A c,   A_jk = x_j' R_0^{-1} x_k,   x in {g, p, q, b},
## so for fixed c_bm (and c_bl) the largest feasible second coupling
## is the positive root of a quadratic. The design ceiling is the
## minimum over its paths, since all paths share one parameter set.
##
## Uses the paper 02 driver's own Sigma (build_sigma_matrix with the
## simulation-core.R model_param conventions: Hybrid, t1/2 = 0.5,
## rho = 0.7, c_1 = 0.2), so the numbers apply to jobs 28/30 and any
## successor. Each closed-form ceiling is verified against the
## smallest eigenvalue of the full R just below and just above it.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/31-prognostic-coupling-ceiling.R
## Output: analysis/data/carryover-sensitivity/g9-prognostic-stress/coupling-ceilings.rds
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

OUT_DIR <- 'analysis/data/carryover-sensitivity/g9-prognostic-stress'

design_set <- build_design_set('Hybrid')
resp_param <- default_resp_param()
baseline_param <- default_baseline_param()

mk_param <- function(c_bm = 0, c_pb = 0, c_tv = 0, c_bl = 0) {
  list(N = 70, carryover_t1half = 0.5, carryover_form = 'exponential',
       weibull_shape = 1, dgp_architecture = 'mvn', c.bm = c_bm,
       c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
       c.bm.bl = c_bl, c.bm.tv = c_tv, b.bm.tv = 0, c.bm.pb = c_pb)
}

## Raw correlation matrix (never repaired: build_sigma_matrix returns
## `correlations` untouched; stop_on_nonpd = FALSE only so that
## unit-coupling probes, which are deliberately infeasible, build).
corr_of <- function(mp, td) {
  build_sigma_matrix(mp, resp_param, baseline_param, td,
                     stop_on_nonpd = FALSE)$correlations
}

## Per path: R_0 and the four unit patterns from the bm row.
path_pieces <- function(td) {
  R0_full <- corr_of(mk_param(), td)
  i_bm <- which(rownames(R0_full) == 'bm')
  R0 <- R0_full[-i_bm, -i_bm]
  pattern <- function(mp) corr_of(mp, td)[i_bm, -i_bm]
  X <- cbind(g = pattern(mk_param(c_bm = 1)),
             p = pattern(mk_param(c_pb = 1)),
             q = pattern(mk_param(c_tv = 1)),
             b = pattern(mk_param(c_bl = 1)))
  stopifnot(min(eigen(R0, symmetric = TRUE, only.values = TRUE)$values) > 0)
  list(R0 = R0, X = X, A = t(X) %*% solve(R0, X))
}
pieces <- map(design_set, path_pieces)

## Largest feasible value of coupling `y` given fixed values of the
## others, from Q(c) = 1, for one path.
ceiling_given <- function(A, y, fixed) {
  others <- setdiff(colnames(A), y)
  cf <- setNames(rep(0, length(others)), others)
  cf[names(fixed)] <- fixed
  a <- A[y, y]
  bq <- 2 * sum(A[y, others] * cf)
  cq <- drop(t(cf) %*% A[others, others] %*% cf) - 1
  disc <- bq^2 - 4 * a * cq
  if (disc < 0) return(NA_real_)
  (-bq + sqrt(disc)) / (2 * a)
}
design_ceiling <- function(y, fixed = c()) {
  min(map_dbl(pieces, ~ ceiling_given(.x$A, y, fixed)))
}

## Eigenvalue verification of a (c_bm, c_pb, c_tv, c_bl) point: the
## smallest eigenvalue of the full R, minimum over paths.
min_eig_full <- function(c_bm = 0, c_pb = 0, c_tv = 0, c_bl = 0) {
  min(map_dbl(design_set, function(td) {
    R <- corr_of(mk_param(c_bm, c_pb, c_tv, c_bl), td)
    min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  }))
}

label <- c(g = 'c_bm', p = 'c_pb', q = 'c_tv', b = 'c_bl')
cases <- tribble(
  ~y,  ~fixed,
  'g', c(),
  'p', c(),
  'q', c(),
  'b', c(),
  'p', c(g = 0.30),
  'q', c(g = 0.30),
  'b', c(g = 0.30),
  'p', c(g = 0.45),
  'q', c(g = 0.45),
  'b', c(g = 0.45),
  'g', c(p = 0.30),
  'g', c(q = 0.30),
  'g', c(p = 0.60),
  'g', c(q = 0.60)
)

results <- pmap_dfr(cases, function(y, fixed) {
  cstar <- design_ceiling(y, fixed)
  args <- function(val) {
    v <- c(g = 0, p = 0, q = 0, b = 0)
    v[names(fixed)] <- fixed
    v[y] <- val
    min_eig_full(v[['g']], v[['p']], v[['q']], v[['b']])
  }
  tibble(
    coupling = label[[y]],
    given = if (length(fixed)) paste(label[names(fixed)], '=', fixed,
                                     collapse = ', ') else 'alone',
    ceiling = cstar,
    min_eig_below = if (is.na(cstar)) NA_real_ else args(0.999 * cstar),
    min_eig_above = if (is.na(cstar)) NA_real_ else args(1.001 * cstar)
  )
})

stopifnot(all(is.na(results$ceiling) |
              (results$min_eig_below > 0 & results$min_eig_above < 0)))

cat('Exact feasibility ceilings, Hybrid (min over 4 paths), t1/2 = 0.5,\n',
    'rho = 0.7, c_1 = 0.2. Verified: smallest eigenvalue of the full R\n',
    'is positive at 0.999 x ceiling and negative at 1.001 x ceiling.\n\n')
print(as.data.frame(results |>
  mutate(ceiling = round(ceiling, 4),
         min_eig_below = signif(min_eig_below, 2),
         min_eig_above = signif(min_eig_above, 2))), row.names = FALSE)

cat('\nQuadratic-form coefficient matrices A (per path, patterns g p q b):\n')
walk(seq_along(pieces), function(i) {
  cat(sprintf('path %d\n', i)); print(round(pieces[[i]]$A, 3))
})

saveRDS(list(results = results, A = map(pieces, 'A'),
             meta = list(date = Sys.time(), design = 'Hybrid',
                         t1half = 0.5, rho = 0.7, c1 = 0.2)),
        file.path(OUT_DIR, 'coupling-ceilings.rds'))
cat(sprintf('\nWritten: %s\n', file.path(OUT_DIR, 'coupling-ceilings.rds')))
