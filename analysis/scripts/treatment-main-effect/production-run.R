# production-run.R
#
# Driver for the paper 04 linear-DGP pilot (multi-crossover versus
# parallel-group RCT, treatment main effect). Version 2, 2026-10-28.
# Replaces the 2026-05-19 driver, whose stored output is kept at
# analysis/data/derived_data/04-mainfx-production.rds.
#
# What it runs (all cells, 2,000 replicates each):
#   M1  effect sweep   beta_D in {0, -0.5, -1, -2}, t_half = 0.2 wk
#   S1  half-life sweep t_half in {0, 0.5, 1, 2} wk, beta_D = -1
#       (t_half = 0.2 is the M1 cell at beta_D = -1; the carryover
#       term is proportional to beta_D, so a null cell does not
#       depend on t_half)
#
# N-of-1 analyses (ALL, XS, DBC, DBCL, and the labeled oracle XO) are
# defined in pilot-core.R, each with nlme containment df and CR2
# cluster-robust SEs with Satterthwaite df. RCT comparators: one visit
# with an error-free covariate (the 2026-05 comparator), one visit with
# covariate reliability 0.8, and eight visits (matched burden).
#
# Random numbers: one L'Ecuyer chain, advanced with nextRNGStream(),
# gives every (cell, replicate) its own stream. Results do not depend
# on the number of cores.
#
# Replicate-level rows are written per cell to
#   analysis/data/derived_data/04-mainfx-pilot/rep-<cell>.rds
# as each cell finishes. A rerun skips cells whose file exists; set
# FORCE <- TRUE to recompute.
#
# Usage (project root):
#   Rscript analysis/scripts/treatment-main-effect/production-run.R
#   Rscript analysis/scripts/treatment-main-effect/production-run.R --smoke
#
# --smoke runs 20 replicates of two cells into a temporary directory
# and does not touch the stored output.

suppressPackageStartupMessages({
  library(nlme)
  library(clubSandwich)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
SMOKE <- '--smoke' %in% args
FORCE <- '--force' %in% args

source(file.path('analysis', 'scripts', 'treatment-main-effect',
                 'pilot-core.R'))

MASTER_SEED <- 20261028L
N_SIM <- if (SMOKE) 20L else 2000L
N_CORES <- max(1L, min(4L, detectCores() - 1L))

Sys.setenv(VECLIB_MAXIMUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1',
           OMP_NUM_THREADS = '1')

out_dir <- if (SMOKE) {
  file.path(tempdir(), '04-mainfx-pilot')
} else {
  file.path('analysis', 'data', 'derived_data', '04-mainfx-pilot')
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells <- rbind(
  data.frame(study = 'M1', true_effect = c(0, -0.5, -1, -2), t_half = 0.2),
  data.frame(study = 'S1', true_effect = -1, t_half = c(0, 0.5, 1, 2))
)
cells$cell <- sprintf('%s-b%s-h%s', cells$study,
                      format(cells$true_effect), format(cells$t_half))
if (SMOKE) cells <- cells[c(2L, 6L), ]

streams <- make_streams(nrow(cells) * N_SIM, MASTER_SEED)

t_total <- proc.time()
for (ci in seq_len(nrow(cells))) {
  cell_file <- file.path(out_dir, paste0('rep-', cells$cell[ci], '.rds'))
  if (file.exists(cell_file) && !FORCE) {
    cat(sprintf('Cell %s exists, skipped\n', cells$cell[ci]))
    next
  }
  offset <- (ci - 1L) * N_SIM
  cat(sprintf('Cell %s: %d replicates on %d cores\n',
              cells$cell[ci], N_SIM, N_CORES))
  t0 <- proc.time()
  reps <- mclapply(
    seq_len(N_SIM),
    function(i) {
      tryCatch(
        run_replicate(i, streams[[offset + i]], cells$true_effect[ci],
                      cells$t_half[ci]),
        error = function(e) NULL
      )
    },
    mc.cores = N_CORES, mc.set.seed = FALSE, mc.preschedule = TRUE
  )
  failed <- vapply(reps, is.null, logical(1))
  stopifnot(
    'more than 5% of replicates failed outright' = mean(failed) < 0.05
  )
  reps <- do.call(rbind, reps[!failed])
  reps$cell <- cells$cell[ci]
  saveRDS(reps, cell_file)
  cat(sprintf('  done in %.1f min; %d failed\n',
              (proc.time() - t0)[['elapsed']] / 60, sum(failed)))
}
cat(sprintf('Total %.1f min\n',
            (proc.time() - t_total)[['elapsed']] / 60))

# ---- Summaries ----------------------------------------------------

files <- file.path(out_dir, paste0('rep-', cells$cell, '.rds'))
stopifnot('a cell file is missing' = all(file.exists(files)))
summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ci) {
  reps <- readRDS(files[ci])
  s <- summarize_cell(reps, cells$true_effect[ci])
  cbind(study = cells$study[ci], true_effect = cells$true_effect[ci],
        t_half = cells$t_half[ci], s)
}))
rownames(summ) <- NULL

paired <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ci) {
  reps <- readRDS(files[ci])
  pairs <- list(c('XS', 'RCT1_r08'), c('XS', 'RCT8_r08'),
                c('DBC', 'RCT8_r08'), c('XS', 'DBC'))
  do.call(rbind, lapply(pairs, function(p) {
    r <- paired_diff(reps, p[1L], p[2L])
    data.frame(cell = cells$cell[ci], true_effect = cells$true_effect[ci],
               t_half = cells$t_half[ci], arm_a = p[1L], arm_b = p[2L],
               diff = r[['diff']], mcse = r[['mcse']],
               mcse_independent = r[['mcse_independent']])
  }))
}))
rownames(paired) <- NULL

results <- list(
  summary = summ,
  paired = paired,
  cells = cells,
  meta = list(
    n_sim = N_SIM, rng_kind = RNGkind()[1L], master_seed = MASTER_SEED,
    n_patients = PILOT$n_patients, n_periods = PILOT$n_periods,
    sigma_u = PILOT$sigma_u, sigma_w = PILOT$sigma_w, rho = PILOT$rho,
    t_half_a = PILOT$t_half_a, r_bl = PILOT$r_bl,
    n_cores = N_CORES, run_date = format(Sys.time(), '%Y-%m-%d %H:%M %Z'),
    r_version = as.character(getRversion()),
    platform = R.version$platform, smoke = SMOKE
  )
)

if (!SMOKE) {
  saveRDS(results, file.path('analysis', 'data', 'derived_data',
                             '04-mainfx-pilot-summary.rds'))
  cat('Summary saved to analysis/data/derived_data/',
      '04-mainfx-pilot-summary.rds\n', sep = '')
} else {
  print(summ[, c('study', 'true_effect', 't_half', 'arm', 'n_ok',
                 'power_c', 'power_r', 'bias', 'emp_se', 'mean_se_c',
                 'mean_se_r')], digits = 3)
}
