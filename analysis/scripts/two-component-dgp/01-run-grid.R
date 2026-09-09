## analysis/scripts/two-component-dgp/01-run-grid.R
##
## Driver for papers 13 and 14 (two-component decompositions).
##
## Paper 13 repeats paper 01's Architecture A versus Architecture B
## comparison under a reduced data-generating process that carries only
## the pharmacological (BR) and time-variant natural-history (TV)
## components, dropping the placebo-belief component (PB). Retaining a
## trend term and omitting a modelled belief component follows the
## convention of the N-of-1 and crossover literature. Two changes from
## paper 01's grid:
##
## NOTE: trialdesign$e (the expectancy weight) drives PB alone. With PB
## dropped it has no effect on the DGP, so the open-label versus
## blinded contrast is invisible to data generation and the designs
## differ only through their on-drug patterns. See Section 2.2.3.
##
##   1. components = c('tv','br') is passed to generateData(), so the
##      covariance matrix is (2 + 2*nP) rather than (2 + 3*nP): 18x18
##      instead of 26x26 at eight occasions.
##   2. The 'orig' (vendored Hendrickson) arm is dropped. Paper 13
##      compares 'mean' (Architecture A) against 'covar' (Architecture
##      B) only, both run through the package's own generateData() with
##      the corrected carryover recursion (see 02-deviations.md).
##
## Grid: 2 architectures x 3 designs x 2 sample sizes x 3 biomarker
## levels x 3 carryover half-lives = 108 cells.
##
## Usage:
##   Rscript analysis/scripts/two-component-dgp/01-run-grid.R [--dev] [--reps N]
##
##   --dev     : 20 reps/cell (smoke test); default 500 (production).
##   --reps N  : override the rep count.
##   --drop X  : 'pb' (paper 13, BR+TV) or 'tv' (paper 14, BR+PB).
##               Defaults to 'pb'.
##
## Writes: analysis/data/{13-drop-pb,14-drop-tv}-grid.rds

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

repo_root <- here::here()
suppressMessages(pkgload::load_all(repo_root, quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
dev_mode <- '--dev' %in% args
reps_idx <- which(args == '--reps')
n_reps <- if (length(reps_idx) && reps_idx < length(args)) {
  as.integer(args[reps_idx + 1])
} else if (dev_mode) 20L else 500L

SEED <- 20260908L

## Which component to drop. --drop pb  -> BR + TV (paper 13)
##                          --drop tv  -> BR + PB (paper 14)
drop_idx <- which(args == '--drop')
DROP <- if (length(drop_idx) && drop_idx < length(args)) {
  args[drop_idx + 1]
} else 'pb'
if (!DROP %in% c('pb', 'tv')) stop("--drop must be 'pb' or 'tv'")
COMPONENTS <- setdiff(c('tv', 'pb', 'br'), DROP)
PAPER <- if (DROP == 'pb') 13L else 14L
OUTFILE <- sprintf('%d-drop-%s-grid.rds', PAPER, DROP)

data(extracted_bp, package = 'pmsimstats', envir = environment())
data(extracted_rp, package = 'pmsimstats', envir = environment())

## ---------------------------------------------------------------
## Trial designs: identical to paper 01 / paper 02
## ---------------------------------------------------------------

designs <- list(
  CO = buildtrialdesign(
    name_shortform = 'CO', timepoints = cumulative(rep(2.5, 8)),
    timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
    expectancies = rep(.5, 8),
    ondrug = list(pathA = c(1,1,1,1,0,0,0,0),
                  pathB = c(0,0,0,0,1,1,1,1))),
  Hybrid = buildtrialdesign(
    name_shortform = 'Hybrid', timepoints = c(4,8,9,10,11,12,16,20),
    timeptnames = c('OL1','OL2','BD1','BD2','BD3','BD4','COd','COp'),
    expectancies = c(1,1,.5,.5,.5,.5,.5,.5),
    ondrug = list(pathA = c(1,1,1,1,0,0,1,0),
                  pathB = c(1,1,1,1,0,0,0,1),
                  pathC = c(1,1,1,0,0,0,1,0),
                  pathD = c(1,1,1,0,0,0,0,1))),
  `OL+BDC` = buildtrialdesign(
    name_shortform = 'OL+BDC', timepoints = c(4,8,12,16,17,18,19,20),
    timeptnames = c('OL1','OL2','OL3','OL4','BD1','BD2','BD3','BD4'),
    expectancies = c(1,1,1,1,.5,.5,.5,.5),
    ondrug = list(pathA = c(1,1,1,1,1,1,0,0),
                  pathB = c(1,1,1,1,1,0,0,0)))
)

## N is the TOTAL number of participants, allocated as evenly as
## possible across the design's randomization paths.
allocate <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rem <- n_total %% n_paths
  rep(base, n_paths) + c(rep(1L, rem), rep(0L, n_paths - rem))
}

## ---------------------------------------------------------------
## One cell: n_reps trials, both architectures fitted to the same
## simulated data within a replicate is not possible (the DGP differs
## by architecture), so each architecture gets its own replicate loop
## from the same seed stream.
## ---------------------------------------------------------------

simulate_cell <- function(cell, n_reps) {
  td <- designs[[cell$design]]
  n_paths <- length(td$trialpaths)
  alloc <- allocate(cell$N, n_paths)

  lam <- if (cell$t1half > 0) log(2) / cell$t1half else 0

  ## Pre-build sigma per path; it does not vary across replicates.
  cached <- lapply(seq_len(n_paths), function(g) {
    mp <- list(N = alloc[g], c.bm = cell$c_bm,
               carryover_t1half = cell$t1half,
               c.tv = .7, c.pb = .7, c.br = .7,
               c.cf1t = .2, c.cfct = .1)
    buildSigma(mp, extracted_rp, extracted_bp, td$trialpaths[[g]],
               makePositiveDefinite = TRUE, lambda_cor = lam,
               dgp_architecture = cell$arch, components = COMPONENTS)
  })

  op <- list(useDE = FALSE, t_random_slope = FALSE,
             full_model_out = FALSE,
             carryover_t1half = cell$t1half,
             carryover_scalefactor = 1)

  out <- vector('list', n_reps)
  for (r in seq_len(n_reps)) {
    dat <- rbindlist(lapply(seq_len(n_paths), function(g) {
      mp <- list(N = alloc[g], c.bm = cell$c_bm,
                 carryover_t1half = cell$t1half,
                 c.tv = .7, c.pb = .7, c.br = .7,
                 c.cf1t = .2, c.cfct = .1)
      d <- generateData(mp, extracted_rp, extracted_bp,
                        td$trialpaths[[g]], FALSE, TRUE,
                        lambda_cor = lam,
                        cached_sigma = cached[[g]],
                        dgp_architecture = cell$arch,
                        components = COMPONENTS)
      d[, path := g]
      d
    }), fill = TRUE)
    dat[, ptID := .I]

    fit <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                    error = function(e) NULL)
    out[[r]] <- if (is.null(fit)) {
      data.table(rep = r, beta = NA_real_, betaSE = NA_real_,
                 p = NA_real_, converged = FALSE)
    } else {
      data.table(rep = r, beta = fit$beta, betaSE = fit$betaSE,
                 p = fit$p, converged = !is.na(fit$beta))
    }
  }
  rbindlist(out)
}

## ---------------------------------------------------------------
## Grid
## ---------------------------------------------------------------

grid <- tidyr::expand_grid(
  arch    = c('mean_moderation', 'mvn'),
  design  = names(designs),
  N       = c(35L, 70L),
  c_bm    = c(0, 0.30, 0.45),
  t1half  = c(0, 0.5, 1.0)
)

cat(sprintf('Paper %d two-component grid: %d cells, n_reps = %d\n',
            PAPER, nrow(grid), n_reps))
cat(sprintf('Components: %s (%s dropped) -> %s\n',
            paste(toupper(COMPONENTS), collapse = ' + '),
            toupper(DROP), OUTFILE))

t0 <- Sys.time()
results <- vector('list', nrow(grid))
for (i in seq_len(nrow(grid))) {
  cell <- grid[i, ]
  set.seed(SEED + i)
  r <- simulate_cell(cell, n_reps)
  results[[i]] <- cbind(cell[rep(1, nrow(r)), ], r)
  cat(sprintf('[%3d/%d] %-16s %-7s N=%2d c_bm=%.2f t1/2=%.1f  power=%.3f\n',
              i, nrow(grid), cell$arch, cell$design, cell$N,
              cell$c_bm, cell$t1half,
              mean(r$p < 0.05, na.rm = TRUE)))
}
results <- rbindlist(results)

summary_tbl <- results |>
  group_by(arch, design, N, c_bm, t1half) |>
  summarise(
    n = dplyr::n(),
    n_converged = sum(converged),
    non_convergence_rate = mean(!converged),
    power = mean(p < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / sum(!is.na(p))),
    ## NOTE: empirical_se must be computed BEFORE `beta` is reassigned.
    ## dplyr evaluates summarise() expressions in order, so putting
    ## `beta = mean(beta)` first would leave sd() operating on a single
    ## scalar and returning NA.
    empirical_se = stats::sd(beta, na.rm = TRUE),
    beta = mean(beta, na.rm = TRUE),
    mean_model_se = mean(betaSE, na.rm = TRUE),
    .groups = 'drop')

out_dir <- file.path(repo_root, 'analysis/data')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(
  results = results,
  summary = summary_tbl,
  grid = grid,
  meta = list(
    script = '01-run-grid.R',
    date = Sys.time(),
    seed = SEED,
    n_reps = n_reps,
    dev_mode = dev_mode,
    components = COMPONENTS,
    dropped = DROP,
    paper = PAPER,
    architectures = c('mean_moderation', 'mvn'),
    elapsed_secs = as.numeric(difftime(Sys.time(), t0, units = 'secs')),
    r_version = R.version.string,
    git_sha = tryCatch(
      system('git rev-parse --short HEAD', intern = TRUE),
      error = function(e) NA_character_))
), file.path(out_dir, OUTFILE))

cat(sprintf('\nWrote analysis/data/%s (%.0f sec)\n', OUTFILE,
            as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
