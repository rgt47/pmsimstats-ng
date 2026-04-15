#!/usr/bin/env Rscript
# test-parity-extended-tidyverse.R
#
# Confirms that implementations/original-extended/ (data.table style,
# camelCase) and implementations/tidyverse/ (tidyverse style, snake_case)
# produce numerically identical results under a fixed random seed, for
# both dgp_architecture = "mvn" (Architecture B) and
# dgp_architecture = "mean_moderation" (Architecture A).
#
# Parity is checked at three levels:
#   1. modgompertz() vs mod_gompertz()
#   2. buildSigma() vs build_sigma_matrix()
#   3. generateData() vs generate_data() under both architectures
#
# Run from the project root:
#   Rscript implementations/test-parity-extended-tidyverse.R

suppressPackageStartupMessages({
  library(data.table)
  library(tibble)
  library(dplyr)
  library(purrr)
  library(nlme)
  library(MASS)
  library(Matrix)
  library(corpcor)
})

# ---- 1. Source both implementations into separate environments ------------

here <- function(...) {
  file.path('implementations', ...)
}

ext_env <- new.env()
for (f in list.files(here('original-extended', 'R'), full.names = TRUE)) {
  sys.source(f, envir = ext_env)
}

tid_env <- new.env()
sys.source(here('tidyverse', 'R', 'functions.R'), envir = tid_env)

# ---- 2. Shared test fixtures ---------------------------------------------

timepoints <- cumsum(rep(2.5, 8))

td_ext <- ext_env$buildtrialdesign(
  name_longform = 'open label',
  name_shortform = 'OL',
  timepoints = timepoints,
  timeptnames = paste0('OL', 1:8),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)
td_ext_path <- td_ext$trialpaths[[1]]

# Tidyverse path expects snake_case timepoint_name; same numeric fields.
td_tid_path <- td_ext_path
td_tid_path$timepoint_name <- td_tid_path$timeptnames

model_param <- data.table(
  N = 35, c.bm = 0.3, carryover_t1half = 1.0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.1, c.cfct = 0.05
)

resp_param_dt <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)
resp_param_tbl <- as_tibble(resp_param_dt)

baseline_param_dt <- data.table(
  cat = c('bm', 'BL'),
  m = c(0, 70),
  sd = c(1, 10)
)
baseline_param_tbl <- as_tibble(baseline_param_dt)

tol <- 1e-10

report <- function(label, result) {
  ok <- isTRUE(result)
  tag <- if (ok) 'PASS' else 'FAIL'
  cat(sprintf('  [%s] %s\n', tag, label))
  if (!ok) cat('       ', format(result), '\n')
  invisible(ok)
}

# ---- 3. Level 1: Gompertz -------------------------------------------------

cat('=== Level 1: Gompertz ===\n')
t_seq <- seq(0, 30, by = 0.5)
g_ext <- ext_env$modgompertz(t_seq, 10.98604, 5, 0.42)
g_tid <- tid_env$mod_gompertz(t_seq, 10.98604, 5, 0.42)
lev1 <- report('modgompertz == mod_gompertz',
               all.equal(g_ext, g_tid, tolerance = tol))

# ---- 4. Level 2: Sigma matrix under both architectures --------------------

cat('\n=== Level 2: Sigma matrix ===\n')
check_sigma <- function(arch) {
  s_ext <- ext_env$buildSigma(
    model_param, resp_param_dt, baseline_param_dt, td_ext_path,
    dgp_architecture = arch
  )
  s_tid <- tid_env$build_sigma_matrix(
    model_param, resp_param_tbl, baseline_param_tbl, td_tid_path,
    factor_types = c('tv', 'pb', 'br'),
    factor_abbreviations = c('tv', 'pb', 'br'),
    dgp_architecture = arch
  )
  list(
    labels = all(s_ext$labels == s_tid$labels),
    means  = all.equal(s_ext$means, s_tid$means, tolerance = tol),
    sigma  = all.equal(unname(s_ext$sigma), unname(s_tid$sigma),
                       tolerance = tol)
  )
}

lev2 <- c()
for (arch in c('mvn', 'mean_moderation')) {
  cat(sprintf('-- %s\n', arch))
  r <- check_sigma(arch)
  lev2 <- c(lev2,
    report(sprintf('labels match (%s)', arch), r$labels),
    report(sprintf('means match  (%s)', arch), r$means),
    report(sprintf('sigma match  (%s)', arch), r$sigma)
  )
}

# ---- 5. Level 3: Data generation under both architectures -----------------

cat('\n=== Level 3: generateData vs generate_data ===\n')

check_data <- function(arch) {
  set.seed(42)
  d_ext <- ext_env$generateData(
    model_param, resp_param_dt, baseline_param_dt, td_ext_path,
    empirical = FALSE, makePositiveDefinite = TRUE,
    dgp_architecture = arch
  )
  set.seed(42)
  d_tid <- tid_env$generate_data(
    model_param, resp_param_tbl, baseline_param_tbl, td_tid_path,
    empirical = FALSE, make_positive_definite = TRUE,
    dgp_architecture = arch
  )

  shared <- intersect(names(d_ext), names(d_tid))
  stopifnot(length(shared) > 0)

  results <- vapply(shared, function(col) {
    isTRUE(all.equal(d_ext[[col]], d_tid[[col]], tolerance = tol))
  }, logical(1))

  list(shared = shared, ok = results,
       nrow_ext = nrow(d_ext), nrow_tid = nrow(d_tid))
}

lev3 <- c()
for (arch in c('mvn', 'mean_moderation')) {
  cat(sprintf('-- %s\n', arch))
  r <- check_data(arch)
  cat(sprintf('   nrow(ext)=%d  nrow(tid)=%d  shared cols=%d\n',
              r$nrow_ext, r$nrow_tid, length(r$shared)))
  for (col in r$shared) {
    lev3 <- c(lev3,
              report(sprintf('col %-10s match (%s)', col, arch),
                     r$ok[[col]]))
  }
}

# ---- 6. Summary -----------------------------------------------------------

cat('\n=== Summary ===\n')
all_ok <- all(c(lev1, lev2, lev3))
cat(sprintf('Level 1 (gompertz)     : %s\n',
            if (isTRUE(lev1)) 'PASS' else 'FAIL'))
cat(sprintf('Level 2 (sigma, 2 arch): %d / %d PASS\n',
            sum(lev2), length(lev2)))
cat(sprintf('Level 3 (data, 2 arch) : %d / %d PASS\n',
            sum(lev3), length(lev3)))
cat('\n')
if (all_ok) {
  cat('ALL PARITY CHECKS PASSED: original-extended and tidyverse\n')
  cat('produce identical results under both dgp_architecture values.\n')
  quit(status = 0)
} else {
  cat('PARITY FAILURES DETECTED. See failed checks above.\n')
  quit(status = 1)
}
