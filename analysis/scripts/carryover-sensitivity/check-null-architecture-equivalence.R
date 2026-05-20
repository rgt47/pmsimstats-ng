## analysis/scripts/carryover-sensitivity/check-null-architecture-equivalence.R
##
## Verifies the claim that at c_bm = 0 the two DGP architectures
## ('mvn' = Architecture B, 'mean_moderation' = Architecture A)
## collapse to one identical data-generating process, so the
## architecture axis carries no information in the null (c_bm = 0)
## cells of the carryover-sensitivity factorial.
##
## Why this should hold, read from implementations/tidyverse/R/
## functions.R:
##
##   1. build_correlation_matrix() writes the biomarker-BR
##      differential-correlation block only inside the
##      `dgp_architecture == "mvn"` branch (functions.R:256-281).
##      The values written are model_param$c.bm on-drug and
##      c.bm * decay off-drug. At c_bm = 0 every such cell is 0.
##      Under 'mean_moderation' the block is never written and the
##      cells keep their identity-initialised value of 0. Hence the
##      correlation matrix, and therefore sigma, is architecture-
##      invariant when c_bm = 0.
##
##   2. generate_data() applies the 'mean_moderation' additive shift
##      beta_bm * bm_z * br_sd with beta_bm = model_param$c.bm
##      (functions.R:422-435). At c_bm = 0 the shift is identically
##      zero.
##
## With a common RNG seed the two architectures should therefore
## emit bit-identical participant data when c_bm = 0, and should
## diverge when c_bm > 0. This script confirms both directions.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/check-null-architecture-equivalence.R

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

resp_param     <- default_resp_param()
baseline_param <- default_baseline_param()

make_model_param <- function(c_bm, t1half) {
  list(
    N                = 40,
    c.bm             = c_bm,
    carryover_t1half = t1half,
    carryover_form   = 'exponential',
    weibull_shape    = 1,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.2, c.cfct = 0.1
  )
}

## generate_data() consumes a single-path trial design; the
## architecture switch is path-independent, so path 1 suffices.
trial_path <- function(design) build_design_set(design)[[1]]

## Level 1: is the covariance structure architecture-invariant?
cmp_sigma <- function(c_bm, t1half, design) {
  td <- trial_path(design)
  mp <- make_model_param(c_bm, t1half)
  s_b <- build_sigma_matrix(mp, resp_param, baseline_param, td,
                            dgp_architecture = 'mvn')
  s_a <- build_sigma_matrix(mp, resp_param, baseline_param, td,
                            dgp_architecture = 'mean_moderation')
  list(
    cor_identical    = identical(s_a$correlations, s_b$correlations),
    max_abs_cor_diff = max(abs(s_a$correlations - s_b$correlations))
  )
}

## Level 2: is the generated participant data architecture-invariant
## under a common seed?
cmp_data <- function(c_bm, t1half, design, seed = 4242) {
  td <- trial_path(design)
  mp <- make_model_param(c_bm, t1half)
  d_b <- generate_data(mp, resp_param, baseline_param, td,
                       empirical = FALSE, make_positive_definite = TRUE,
                       seed = seed, dgp_architecture = 'mvn')
  d_a <- generate_data(mp, resp_param, baseline_param, td,
                       empirical = FALSE, make_positive_definite = TRUE,
                       seed = seed, dgp_architecture = 'mean_moderation')
  num_b <- as.matrix(d_b |> dplyr::select(where(is.numeric)))
  num_a <- as.matrix(d_a |> dplyr::select(where(is.numeric)))
  list(
    data_identical = isTRUE(all.equal(num_a, num_b)),
    max_abs_diff   = max(abs(num_a - num_b), na.rm = TRUE)
  )
}

grid <- expand.grid(
  design = c('CO', 'Hybrid'),
  t1half = c(0.0, 1.0),
  c_bm   = c(0.0, 0.45),
  stringsAsFactors = FALSE
)

cat(sprintf('%-8s %7s %6s | %-13s %-14s | %-13s %-14s\n',
            'design', 't1half', 'c_bm',
            'cor_identical', 'max_cor_diff',
            'data_identical', 'max_data_diff'))
cat(strrep('-', 82), '\n')

ok <- TRUE
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  s <- cmp_sigma(g$c_bm, g$t1half, g$design)
  d <- cmp_data(g$c_bm, g$t1half, g$design)
  cat(sprintf('%-8s %7.1f %6.2f | %-13s %-14.3e | %-13s %-14.3e\n',
              g$design, g$t1half, g$c_bm,
              s$cor_identical, s$max_abs_cor_diff,
              d$data_identical, d$max_abs_diff))
  expect_identical <- g$c_bm == 0
  row_ok <- (s$cor_identical == expect_identical) &&
            (d$data_identical == expect_identical)
  ok <- ok && row_ok
}

cat(strrep('-', 82), '\n')
cat('Expectation: c_bm = 0 rows identical (TRUE, 0 diff);',
    'c_bm = 0.45 rows divergent.\n')
cat(if (ok) 'RESULT: assertion CONFIRMED.\n'
    else    'RESULT: assertion NOT confirmed - inspect rows above.\n')
