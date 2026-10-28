## analysis/scripts/quick-sim/carryover-closed-form/18-prognostic-joint-ceiling.R
##
## Joint feasibility ceilings for a prognostic coupling (c.bm.pb,
## c.bm.tv or c.bm.bl) alongside the predictive coupling c.bm, under
## paper 01's refined-process configuration (separable AR(1), rho = 0.7,
## c_1 = 0.2) for the Hybrid, OL+BDC and crossover designs and carryover
## half-lives of 0, 0.2 and 0.5 week. The method is that of
## analysis/scripts/carryover-sensitivity/31-prognostic-coupling-ceiling.R
## (docs/65, Section 4): with the biomarker row written as
##     r = c_bm g + c_pb p + c_tv q + c_bl b,
## the correlation matrix is positive definite iff R_0 is and
## c' A c < 1 with A_jk = x_j' R_0^{-1} x_k. The ceiling on one coupling
## given the others is the positive root of a quadratic; a design's
## ceiling is the minimum over its paths. Every ceiling is verified
## against the smallest eigenvalue of the full R just below and above.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/18-prognostic-joint-ceiling.R
## Output:
##   analysis/data/quick-sim/carryover-closed-form/rho-stage/prognostic-joint-ceilings.csv
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

out_dir <- 'analysis/data/quick-sim/carryover-closed-form/rho-stage'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

resp_param <- default_resp_param()
baseline_param <- default_baseline_param()

mk_param <- function(t1half, c_bm = 0, c_pb = 0, c_tv = 0, c_bl = 0) {
  list(N = 70, carryover_t1half = t1half, carryover_form = 'exponential',
       weibull_shape = 1, dgp_architecture = 'mvn', c.bm = c_bm,
       c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
       c.bm.bl = c_bl, c.bm.tv = c_tv, b.bm.tv = 0, c.bm.pb = c_pb)
}

## stop_on_nonpd = FALSE only so that the unit-coupling probes, which
## are deliberately infeasible, build; `correlations` is never repaired.
corr_of <- function(mp, td) {
  build_sigma_matrix(mp, resp_param, baseline_param, td,
                     stop_on_nonpd = FALSE)$correlations
}

path_pieces <- function(td, t1half) {
  R0_full <- corr_of(mk_param(t1half), td)
  i_bm <- which(rownames(R0_full) == 'bm')
  R0 <- R0_full[-i_bm, -i_bm]
  pattern <- function(mp) corr_of(mp, td)[i_bm, -i_bm]
  X <- cbind(g = pattern(mk_param(t1half, c_bm = 1)),
             p = pattern(mk_param(t1half, c_pb = 1)),
             q = pattern(mk_param(t1half, c_tv = 1)),
             b = pattern(mk_param(t1half, c_bl = 1)))
  stopifnot(min(eigen(R0, symmetric = TRUE, only.values = TRUE)$values) > 0)
  list(A = t(X) %*% solve(R0, X))
}

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

label <- c(g = 'c_bm', p = 'c_pb', q = 'c_tv', b = 'c_bl')
cases <- tribble(
  ~y,  ~fixed,
  'g', c(),
  'p', c(),
  'b', c(),
  'p', c(g = 0.25),
  'b', c(g = 0.25),
  'p', c(g = 0.30),
  'b', c(g = 0.30),
  'p', c(g = 0.45),
  'b', c(g = 0.45),
  'g', c(p = 0.30),
  'g', c(p = 0.60)
)

one_config <- function(design, t1half) {
  design_set <- build_design_set(design)
  pieces <- map(design_set, path_pieces, t1half = t1half)
  min_eig_full <- function(c_bm, c_pb, c_tv, c_bl) {
    min(map_dbl(design_set, function(td) {
      R <- corr_of(mk_param(t1half, c_bm, c_pb, c_tv, c_bl), td)
      min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
    }))
  }
  pmap_dfr(cases, function(y, fixed) {
    cstar <- min(map_dbl(pieces, ~ ceiling_given(.x$A, y, fixed)))
    at <- function(val) {
      v <- c(g = 0, p = 0, q = 0, b = 0)
      v[names(fixed)] <- fixed
      v[y] <- val
      min_eig_full(v[['g']], v[['p']], v[['q']], v[['b']])
    }
    tibble(design = design, t1half = t1half, coupling = label[[y]],
           given = if (length(fixed))
             paste(label[names(fixed)], '=', fixed, collapse = ', ')
           else 'alone',
           ceiling = cstar,
           min_eig_below = if (is.na(cstar)) NA_real_ else at(0.999 * cstar),
           min_eig_above = if (is.na(cstar)) NA_real_ else at(1.001 * cstar))
  })
}

grid <- expand_grid(design = c('Hybrid', 'OLBDC', 'CO'),
                    t1half = c(0, 0.2, 0.5))
results <- pmap_dfr(grid, one_config)

stopifnot(all(is.na(results$ceiling) |
              (results$min_eig_below > 0 & results$min_eig_above < 0)))

## The PB and TV patterns are structurally identical (docs/65, 4.3), so
## c_tv ceilings equal c_pb ceilings; checked once here.
pb_tv_check <- one_config('Hybrid', 0.2)
hyb <- build_design_set('Hybrid')
pc <- map(hyb, path_pieces, t1half = 0.2)
stopifnot(abs(min(map_dbl(pc, ~ ceiling_given(.x$A, 'q', c()))) -
              min(map_dbl(pc, ~ ceiling_given(.x$A, 'p', c())))) < 1e-10)

## Tie to the single-coupling ceiling of paper 01 (rho-stage, sep, t 0.2):
## Hybrid 0.480 at t 0 reported in docs/34 used c_x = 0.1; here c_1 = c_x.
print(as.data.frame(results |>
  mutate(ceiling = round(ceiling, 3),
         min_eig_below = signif(min_eig_below, 2),
         min_eig_above = signif(min_eig_above, 2))), row.names = FALSE)

write.csv(results, file.path(out_dir, 'prognostic-joint-ceilings.csv'),
          row.names = FALSE)
cat(sprintf('\nWritten: %s\n',
            file.path(out_dir, 'prognostic-joint-ceilings.csv')))
