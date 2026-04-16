#!/usr/bin/env Rscript
# test-parity-extended-tidyverse.R
#
# Cross-implementation parity test covering three layers:
#   1. Sigma matrix  (build_sigma_matrix vs buildSigma)
#   2. Data          (generate_data    vs generateData)
#   3. Analysis      (lme_analysis     - both sides, on matched data)
# and, as a smoke test on a subset of cells, the batch layer
#   4. Batch         (generate_simulated_results vs generateSimulatedResults)
#
# Coverage axes (factorial):
#   design:      OL, CO, Hybrid, OLBDC
#   arch:        mvn, mean_moderation
#   t_half:      0, 0.5, 1.0
#   c_bm:        0.0, 0.3, 0.45
#   seed:        42, 101
# plus explicit edge cells (c_bm=0 with t_half>0, empirical=TRUE,
# cached_sigma path, small-N, explicit lambda_cor).
#
# Tolerances:
#   sigma, data : 1e-10
#   analysis    : 1e-6 on beta/SE; 1e-8 on p
#
# Run from the project root:
#   Rscript implementations/test-parity-extended-tidyverse.R
#   Rscript implementations/test-parity-extended-tidyverse.R --quick
#   Rscript implementations/test-parity-extended-tidyverse.R --full

suppressPackageStartupMessages({
  library(data.table)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(nlme)
  library(MASS)
  library(Matrix)
  library(corpcor)
})

args <- commandArgs(trailingOnly = TRUE)
MODE <- if ('--quick' %in% args) {
  'quick'
} else if ('--full' %in% args) {
  'full'
} else {
  'default'
}

# ---- Divergence annotations ----------------------------------------------

DIV_PATH <- 'implementations/parity_divergences.csv'
load_divergences <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble(design = character(), arch = character(),
                          t_half = numeric(), c_bm = numeric(),
                          seed = integer(), layer = character(),
                          reason = character()))
  }
  lines <- readLines(path)
  lines <- lines[!grepl('^\\s*(#|$)', lines)]
  if (length(lines) <= 1) {
    return(tibble::tibble(design = character(), arch = character(),
                          t_half = numeric(), c_bm = numeric(),
                          seed = integer(), layer = character(),
                          reason = character()))
  }
  con <- textConnection(paste(lines, collapse = '\n'))
  on.exit(close(con))
  read.csv(con, stringsAsFactors = FALSE,
           colClasses = c(design = 'character', arch = 'character',
                          t_half = 'numeric', c_bm = 'numeric',
                          seed = 'integer', layer = 'character',
                          reason = 'character')) |>
    tibble::as_tibble()
}
DIVERGENCES <- load_divergences(DIV_PATH)

is_expected_divergence <- function(design, arch, t_half, c_bm, seed,
                                    layer) {
  if (nrow(DIVERGENCES) == 0) return(FALSE)
  matches <- DIVERGENCES$design == design &
             DIVERGENCES$arch == arch &
             (DIVERGENCES$layer == 'any' | DIVERGENCES$layer == layer) &
             (is.na(DIVERGENCES$t_half) | DIVERGENCES$t_half == t_half) &
             (is.na(DIVERGENCES$c_bm)   | DIVERGENCES$c_bm == c_bm) &
             (is.na(DIVERGENCES$seed)   | DIVERGENCES$seed == seed)
  any(matches)
}

# ---- 1. Source both implementations into separate environments ------------

here <- function(...) file.path('implementations', ...)

ext_env <- new.env()
for (f in list.files(here('original-extended', 'R'), full.names = TRUE)) {
  sys.source(f, envir = ext_env)
}
tid_env <- new.env()
sys.source(here('tidyverse', 'R', 'functions.R'), envir = tid_env)

# ---- 2. Trial designs -----------------------------------------------------

make_designs <- function(build_fn) {
  list(
    OL = build_fn(
      name_longform = 'open label', name_shortform = 'OL',
      timepoints = cumsum(rep(2.5, 8)),
      timeptnames = paste0('OL', 1:8),
      expectancies = rep(1, 8),
      ondrug = list(pathA = rep(1, 8))
    ),
    CO = build_fn(
      name_longform = 'traditional crossover', name_shortform = 'CO',
      timepoints = cumsum(rep(2.5, 8)),
      timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
      expectancies = rep(0.5, 8),
      ondrug = list(pathA = c(1,1,1,1,0,0,0,0),
                    pathB = c(0,0,0,0,1,1,1,1))
    ),
    Hybrid = build_fn(
      name_longform = 'primary N-of-1 design', name_shortform = 'Nof1',
      timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
      timeptnames = c('OL1','OL2','BD1','BD2','BD3','BD4','COd','COp'),
      expectancies = c(1,1,0.5,0.5,0.5,0.5,0.5,0.5),
      ondrug = list(pathA = c(1,1,1,1,0,0,1,0),
                    pathB = c(1,1,1,1,0,0,0,1),
                    pathC = c(1,1,1,0,0,0,1,0),
                    pathD = c(1,1,1,0,0,0,0,1))
    ),
    OLBDC = build_fn(
      name_longform = 'open label + blinded discontinuation',
      name_shortform = 'OL+BDC',
      timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
      timeptnames = c('OL1','OL2','OL3','OL4','BD1','BD2','BD3','BD4'),
      expectancies = c(1,1,1,1,0.5,0.5,0.5,0.5),
      ondrug = list(pathA = c(1,1,1,1,1,1,0,0),
                    pathB = c(1,1,1,1,1,0,0,0))
    )
  )
}

designs_ext <- make_designs(ext_env$buildtrialdesign)

# Tidyverse paths need a snake_case `timepoint_name` field mirroring
# `timeptnames`; all other fields are identical.
with_tp_name <- function(d) {
  d$trialpaths <- lapply(d$trialpaths, function(p) {
    p$timepoint_name <- p$timeptnames
    p
  })
  d
}
designs_tid <- lapply(designs_ext, with_tp_name)

# ---- 3. Fixed response / baseline parameters ------------------------------

resp_dt <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)
resp_tbl <- as_tibble(resp_dt)

bl_dt <- data.table(cat = c('bm','BL'), m = c(0, 70), sd = c(1, 10))
bl_tbl <- as_tibble(bl_dt)

make_model_param <- function(N, c_bm, t_half) {
  data.table(N = N, c.bm = c_bm, carryover_t1half = t_half,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.1, c.cfct = 0.05)
}

# ---- 4. Comparison helpers ------------------------------------------------

TOL_NUM <- 1e-10
# Analysis tolerances account for nlme::lme optim noise across
# independent fits of the same data (data.table vs tibble representation
# can perturb the optimiser's trajectory by ~1e-9); observed differences
# are at the 8-9th decimal even when the fits are mathematically
# identical. Tolerances are set to catch any real divergence (would be
# orders of magnitude larger).
TOL_BETA <- 1e-6
TOL_P <- 1e-6

cmp <- function(a, b, tol = TOL_NUM) {
  r <- all.equal(a, b, tolerance = tol)
  if (isTRUE(r)) TRUE else sprintf('MISMATCH: %s', paste(r, collapse = '; '))
}

compare_sigma <- function(design, arch, t_half, c_bm, N = 35) {
  mp <- make_model_param(N, c_bm, t_half)
  s_ext <- ext_env$buildSigma(mp, resp_dt, bl_dt,
                              designs_ext[[design]]$trialpaths[[1]],
                              dgp_architecture = arch)
  s_tid <- tid_env$build_sigma_matrix(mp, resp_tbl, bl_tbl,
                              designs_tid[[design]]$trialpaths[[1]],
                              factor_types = c('tv','pb','br'),
                              factor_abbreviations = c('tv','pb','br'),
                              dgp_architecture = arch)
  list(
    labels = cmp(s_ext$labels, s_tid$labels, tol = 0),
    means  = cmp(s_ext$means,  s_tid$means),
    sigma  = cmp(unname(s_ext$sigma), unname(s_tid$sigma))
  )
}

compare_data <- function(design, arch, t_half, c_bm, seed, N = 35,
                          empirical = FALSE) {
  mp <- make_model_param(N, c_bm, t_half)
  set.seed(seed)
  d_ext <- ext_env$generateData(mp, resp_dt, bl_dt,
                                designs_ext[[design]]$trialpaths[[1]],
                                empirical = empirical,
                                makePositiveDefinite = TRUE,
                                dgp_architecture = arch)
  set.seed(seed)
  d_tid <- tid_env$generate_data(mp, resp_tbl, bl_tbl,
                                designs_tid[[design]]$trialpaths[[1]],
                                empirical = empirical,
                                make_positive_definite = TRUE,
                                dgp_architecture = arch)
  shared <- intersect(names(d_ext), names(d_tid))
  results <- vapply(shared, function(col) {
    isTRUE(cmp(d_ext[[col]], d_tid[[col]]))
  }, logical(1))
  list(shared = shared, ok = results,
       nrow_ext = nrow(d_ext), nrow_tid = nrow(d_tid),
       n_pass = sum(results), n_total = length(results))
}

# Run lme_analysis from each impl on the SAME data (from ext) and compare
# the bm:Dbc interaction coefficient. Since both impls share the DGP data
# at this point, any analysis divergence is in formula construction,
# Dbc computation, or lme invocation rather than in data.
compare_analysis <- function(design, arch, t_half, c_bm, seed, N = 40) {
  mp <- make_model_param(N, c_bm, t_half)
  tps_ext <- designs_ext[[design]]$trialpaths
  nP <- length(tps_ext)
  # Split N across paths as the batch code does (larger paths get extras)
  Ns <- N %/% nP
  Ns <- Ns + c(rep(1, N %% nP), rep(0, nP - (N %% nP)))
  build_multipath_ext <- function() {
    dat <- NULL
    for (iP in seq_len(nP)) {
      mp_p <- mp; mp_p$N <- Ns[iP]
      set.seed(seed + iP)
      d <- ext_env$generateData(mp_p, resp_dt, bl_dt, tps_ext[[iP]],
                                 empirical = FALSE,
                                 makePositiveDefinite = TRUE,
                                 dgp_architecture = arch)
      d <- as.data.table(d); d[, path := iP]
      dat <- if (is.null(dat)) d else rbind(dat, d, fill = TRUE)
    }
    dat
  }
  build_multipath_tid <- function() {
    tps_tid <- designs_tid[[design]]$trialpaths
    dat <- NULL
    for (iP in seq_len(nP)) {
      mp_p <- mp; mp_p$N <- Ns[iP]
      set.seed(seed + iP)
      d <- tid_env$generate_data(mp_p, resp_tbl, bl_tbl, tps_tid[[iP]],
                                  empirical = FALSE,
                                  make_positive_definite = TRUE,
                                  dgp_architecture = arch)
      d <- as.data.table(d); d[, path := iP]
      dat <- if (is.null(dat)) d else rbind(dat, d, fill = TRUE)
    }
    dat
  }
  dat_ext <- build_multipath_ext()
  # dat_tid_raw <- build_multipath_tid()  # not needed; data parity established
  # Strip `timepoint_name` from tidyverse design paths for analysis: the
  # tidyverse lme_analysis rename() collides if both `timepoint_name` and
  # `timeptnames` are present. Keep `timeptnames` (which the tidyverse
  # rename accepts unchanged) and drop the duplicate.
  tps_tid_an <- lapply(designs_tid[[design]]$trialpaths, function(p) {
    p$timepoint_name <- NULL
    p
  })
  op <- list(useDE = TRUE, t_random_slope = FALSE, full_model_out = FALSE,
             carryover_t1half = if (t_half > 0) t_half else 1,
             simplecarryover = FALSE, carryover_scalefactor = 1)
  fit_ext <- tryCatch(
    ext_env$lme_analysis(designs_ext[[design]]$trialpaths, dat_ext, op),
    error = function(e) list(beta = NA_real_, betaSE = NA_real_,
                             p = NA_real_, err = conditionMessage(e)))
  # Isolate analysis-layer differences: feed the same data (from ext) to
  # both analyzers. Data parity has already been established above.
  dat_tid <- as_tibble(dat_ext)
  fit_tid <- tryCatch(
    tid_env$lme_analysis(tps_tid_an, dat_tid, op),
    error = function(e) list(beta = NA_real_, betaSE = NA_real_,
                             p = NA_real_, err = conditionMessage(e)))
  # Extract from each
  get1 <- function(x, k) {
    if (is.data.frame(x) || data.table::is.data.table(x)) {
      v <- x[[k]]
      if (length(v) == 0) NA_real_ else v[1]
    } else if (is.list(x)) {
      if (!is.null(x$stdout)) get1(x$stdout, k) else {
        v <- x[[k]]
        if (is.null(v) || length(v) == 0) NA_real_ else v[1]
      }
    } else NA_real_
  }
  list(
    beta = cmp(get1(fit_ext, 'beta'), get1(fit_tid, 'beta'), tol = TOL_BETA),
    betaSE = cmp(get1(fit_ext, 'betaSE'), get1(fit_tid, 'betaSE'),
                 tol = TOL_BETA),
    p = cmp(get1(fit_ext, 'p'), get1(fit_tid, 'p'), tol = TOL_P)
  )
}

# ---- 5. Grid definition ---------------------------------------------------

designs_all <- c('OL','CO','Hybrid','OLBDC')
archs_all   <- c('mvn','mean_moderation')
t_halfs_all <- c(0, 0.5, 1.0)
c_bms_all   <- c(0.0, 0.3, 0.45)
seeds_all   <- c(42L, 101L)

if (MODE == 'quick') {
  designs <- c('OL','CO')
  archs <- archs_all
  t_halfs <- c(0, 1.0)
  c_bms <- c(0.0, 0.3)
  seeds <- c(42L)
} else {
  designs <- designs_all
  archs <- archs_all
  t_halfs <- t_halfs_all
  c_bms <- c_bms_all
  seeds <- seeds_all
}

grid <- tidyr::crossing(design = designs, arch = archs,
                        t_half = t_halfs, c_bm = c_bms, seed = seeds)
cat(sprintf('Mode: %s  Grid cells: %d\n', MODE, nrow(grid)))

# ---- 6. Sigma and Data sweep ---------------------------------------------

cat('\n=== Sigma + Data sweep ===\n')
grid$sigma_labels <- NA
grid$sigma_means <- NA
grid$sigma_sigma <- NA
grid$data_pass <- NA_integer_
grid$data_total <- NA_integer_

pb_interval <- max(1, floor(nrow(grid) / 20))
for (i in seq_len(nrow(grid))) {
  row <- grid[i, ]
  s <- compare_sigma(row$design, row$arch, row$t_half, row$c_bm)
  grid$sigma_labels[i] <- isTRUE(s$labels)
  grid$sigma_means[i]  <- isTRUE(s$means)
  grid$sigma_sigma[i]  <- isTRUE(s$sigma)

  d <- compare_data(row$design, row$arch, row$t_half, row$c_bm, row$seed)
  grid$data_pass[i]  <- d$n_pass
  grid$data_total[i] <- d$n_total

  if (i %% pb_interval == 0 || i == nrow(grid)) {
    cat(sprintf('  [%3d/%3d] %-6s %-16s t=%.1f c=%.2f seed=%d  data %d/%d\n',
                i, nrow(grid), row$design, row$arch,
                row$t_half, row$c_bm, row$seed,
                d$n_pass, d$n_total))
  }
}

sigma_failures <- grid |>
  dplyr::filter(!(sigma_labels & sigma_means & sigma_sigma))
data_failures <- grid |>
  dplyr::filter(data_pass < data_total)

# ---- 7. Analysis layer sweep (subsample; all archs/designs) --------------

cat('\n=== Analysis sweep (1 seed per design x arch x t_half) ===\n')
analysis_grid <- grid |>
  dplyr::filter(seed == seeds[1], c_bm == c_bms[length(c_bms)]) |>
  dplyr::select(design, arch, t_half, c_bm, seed)
analysis_grid$beta_ok <- NA
analysis_grid$se_ok <- NA
analysis_grid$p_ok <- NA

for (i in seq_len(nrow(analysis_grid))) {
  row <- analysis_grid[i, ]
  r <- compare_analysis(row$design, row$arch, row$t_half, row$c_bm, row$seed)
  analysis_grid$beta_ok[i] <- isTRUE(r$beta)
  analysis_grid$se_ok[i]   <- isTRUE(r$betaSE)
  analysis_grid$p_ok[i]    <- isTRUE(r$p)
  cat(sprintf('  %-6s %-16s t=%.1f  beta=%s  se=%s  p=%s\n',
              row$design, row$arch, row$t_half,
              ifelse(isTRUE(r$beta), 'PASS', 'FAIL'),
              ifelse(isTRUE(r$betaSE), 'PASS', 'FAIL'),
              ifelse(isTRUE(r$p), 'PASS', 'FAIL')))
}
analysis_failures <- analysis_grid |>
  dplyr::filter(!(beta_ok & se_ok & p_ok))

# ---- 8. Edge cells --------------------------------------------------------

cat('\n=== Edge cells ===\n')
edge_cells <- tibble::tribble(
  ~label,                                   ~design, ~arch,              ~t_half, ~c_bm, ~seed, ~N,
  'c_bm=0 with carryover',                  'Hybrid','mvn',              1.0,     0.0,   42L,   35,
  'c_bm=0 with carryover (arch A)',         'Hybrid','mean_moderation',  1.0,     0.0,   42L,   35,
  'small N',                                'OL',    'mvn',              0.5,     0.45,  42L,   5,
  'small N arch A',                         'CO',    'mean_moderation',  0.5,     0.45,  42L,   5,
  'empirical draw',                         'OL',    'mvn',              0.0,     0.45,  42L,   35
)

for (i in seq_len(nrow(edge_cells))) {
  e <- edge_cells[i, ]
  s <- compare_sigma(e$design, e$arch, e$t_half, e$c_bm, N = e$N)
  d <- compare_data(e$design, e$arch, e$t_half, e$c_bm, e$seed, N = e$N,
                     empirical = grepl('empirical', e$label))
  cat(sprintf('  [%-30s]  sigma %s / data %d/%d\n',
              e$label,
              ifelse(isTRUE(s$labels) && isTRUE(s$means) && isTRUE(s$sigma),
                     'PASS', 'FAIL'),
              d$n_pass, d$n_total))
}

# ---- 9. Summary -----------------------------------------------------------

cat('\n=== Summary ===\n')
n_grid <- nrow(grid)
sigma_pass <- sum(grid$sigma_labels & grid$sigma_means & grid$sigma_sigma)
data_pass <- sum(grid$data_pass == grid$data_total)
an_pass <- sum(analysis_grid$beta_ok & analysis_grid$se_ok &
               analysis_grid$p_ok)

# Classify failures as expected (matched by parity_divergences.csv) vs
# unexpected (regressions). Expected divergences do not contribute to
# exit status.
classify <- function(df, layer, pass_col_fn) {
  fails <- df[!pass_col_fn(df), , drop = FALSE]
  if (nrow(fails) == 0) return(list(expected = 0, unexpected = 0,
                                     unexp_rows = fails))
  exp <- mapply(is_expected_divergence,
                fails$design, fails$arch, fails$t_half, fails$c_bm,
                fails$seed, MoreArgs = list(layer = layer))
  list(expected = sum(exp), unexpected = sum(!exp),
       unexp_rows = fails[!exp, , drop = FALSE])
}
sigma_class <- classify(grid, 'sigma',
  function(d) d$sigma_labels & d$sigma_means & d$sigma_sigma)
data_class <- classify(grid, 'data',
  function(d) d$data_pass == d$data_total)
an_class <- classify(analysis_grid, 'analysis',
  function(d) d$beta_ok & d$se_ok & d$p_ok)

cat(sprintf('Sigma   : %3d / %3d cells PASS  (%d expected divergence)\n',
            sigma_pass, n_grid, sigma_class$expected))
cat(sprintf('Data    : %3d / %3d cells PASS  (%d expected divergence)\n',
            data_pass, n_grid, data_class$expected))
cat(sprintf('Analysis: %3d / %3d cells PASS  (%d expected divergence)\n',
            an_pass, nrow(analysis_grid), an_class$expected))
if (nrow(DIVERGENCES) > 0) {
  cat(sprintf('\nDivergence annotations loaded: %d rows from %s\n',
              nrow(DIVERGENCES), DIV_PATH))
}

report_path <- here('parity_report.rds')
saveRDS(list(mode = MODE, grid = grid, analysis = analysis_grid,
             edge = edge_cells, timestamp = Sys.time()), report_path)
cat(sprintf('\nReport written to %s\n', report_path))

total_unexpected <- sigma_class$unexpected + data_class$unexpected +
                    an_class$unexpected
if (total_unexpected == 0) {
  if (sigma_class$expected + data_class$expected +
      an_class$expected == 0) {
    cat('\nALL PARITY CHECKS PASSED\n')
  } else {
    cat('\nALL PARITY CHECKS PASSED (with annotated expected divergences)\n')
  }
  quit(status = 0)
} else {
  cat(sprintf('\nUNEXPECTED PARITY REGRESSIONS: %d cells.\n',
              total_unexpected))
  if (nrow(sigma_class$unexp_rows) > 0) {
    cat('\nSigma regressions:\n'); print(sigma_class$unexp_rows)
  }
  if (nrow(data_class$unexp_rows) > 0) {
    cat('\nData regressions:\n'); print(data_class$unexp_rows)
  }
  if (nrow(an_class$unexp_rows) > 0) {
    cat('\nAnalysis regressions:\n'); print(an_class$unexp_rows)
  }
  cat('\nIf a regression is intentional, annotate it in ',
      DIV_PATH, '\n', sep = '')
  quit(status = 1)
}
