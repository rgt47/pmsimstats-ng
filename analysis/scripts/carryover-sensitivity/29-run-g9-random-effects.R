## analysis/scripts/carryover-sensitivity/29-run-g9-random-effects.R
##
## Can a richer random-effects structure recover the power G9 loses
## under a prognostic confound? Job 28 (28-run-g9-prognostic-
## stress.R) found that a strength-0.6 prog-PB or prog-TV-const
## confound costs every specification, G9 and the F-family alike,
## about 20 points of power at N = 35 and about 10 at N = 70, with
## Type I error at nominal throughout. The loss is therefore added
## patient-level variance correlated with the biomarker, not bias,
## and no fixed effect shared across patients (the F-family's De
## term included) can absorb it. This driver tests the variance-
## component alternative: G9 (AIC-selected t1/2 + CR2, random
## intercept) against the same model with
##
##   G9rs   a random slope on t          (patient-specific trend; TV)
##   G9re   a random coefficient on De   (patient-specific placebo
##                                        responsiveness; PB)
##   G9rse  both
##
## at the two affected pathways, strength 0.6, N in {35, 70},
## c_bm in {0, 0.45}, plus matching no-confound cells so that any
## power or size cost of the richer structures in the absence of a
## confound (and their non-convergence rate at N = 35) is measured
## on the same footing. Hybrid design, corrected DGP constants, same
## seeds as job 28, so the confounded cells' trials are identical to
## job 28's when run on the same platform. All four variants are fit
## to the same trial each replicate (simulate_cell_g9_random(),
## simulation-core.R), so comparisons are paired.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/29-run-g9-random-effects.R \
##     [--reps R] [--cores K] [--out-dir PATH]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
REPS <- as.integer(argv('--reps', '500'))
CORES <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
OUT_DIR <- argv('--out-dir',
  'analysis/data/carryover-sensitivity/g9-random-effects')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

base_cell <- list(design = 'Hybrid', t1half = 0.5,
                  carryover_form = 'exponential', weibull_shape = 1,
                  dgp_arch = 'mvn')

pathway_field <- c(none = NA_character_, progtvconst = 'c_tv',
                   progpb = 'c_pb')
STRENGTH <- 0.6

grid <- tidyr::expand_grid(
  pathway = names(pathway_field),
  c_bm    = c(0, 0.45),
  N       = c(35, 70)
) |>
  dplyr::mutate(
    strength = ifelse(pathway == 'none', 0, STRENGTH),
    cell_id = sprintf('%s-conf%.1f-cbm%.2f-N%d', pathway, strength, c_bm, N)
  )

build_cell <- function(i) {
  cl <- modifyList(base_cell, list(N = grid$N[i], c_bm = grid$c_bm[i]))
  fld <- pathway_field[[grid$pathway[i]]]
  if (!is.na(fld)) cl[[fld]] <- grid$strength[i]
  cl
}

cat(sprintf(
  'G9 random-effects variants: %d cells (%d pathways x %d c_bm x %d N), n_reps = %d\n',
  nrow(grid), length(pathway_field), 2, 2, REPS))

jobs <- tidyr::expand_grid(cell_i = seq_len(nrow(grid)), rep = seq_len(REPS))

one_job <- function(i) {
  cell_i <- jobs$cell_i[i]; r <- jobs$rep[i]
  cell_id <- grid$cell_id[cell_i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(25102600L + r)
  res <- tryCatch(
    simulate_cell_g9_random(build_cell(cell_i), n_reps = 1),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    res <- dplyr::mutate(res, rep = r, cell_id = cell_id, .before = 1)
    saveRDS(res, f)
  }
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job,
                   mc.cores = CORES, mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

res <- bind_rows(lapply(list.files(ck_dir, full.names = TRUE), readRDS))
res <- dplyr::left_join(res, grid, by = 'cell_id')
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
n_done <- length(unique(paste(res$cell_id, res$rep)))
cat(sprintf(
  '\n%d of %d (cell, rep) jobs present in %.1f min (this run; checkpointed jobs from earlier runs are not re-timed)\n',
  n_done, nrow(jobs), elapsed_min))

## ---------------- summary: rejection rate, non-convergence, SE per variant ----------------
summary_tab <- res |>
  group_by(cell_id, pathway, strength, c_bm, N, method) |>
  summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
            n = sum(!is.na(p_value)),
            nonconv = mean(!converged),
            mean_se = mean(std_error, na.rm = TRUE),
            emp_se = sd(estimate, na.rm = TRUE),
            .groups = 'drop') |>
  mutate(method = factor(method, levels = names(G9_RANDOM_VARIANTS))) |>
  arrange(pathway, c_bm, N, method)

cat('\nRejection rate (Type I error at c_bm=0 cells, power at c_bm=0.45 cells):\n')
print(as.data.frame(summary_tab), max = 2000)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))

meta <- list(
  script = '29-run-g9-random-effects.R',
  date = Sys.time(),
  n_reps = REPS,
  variants = G9_RANDOM_VARIANTS,
  elapsed_min_this_run = elapsed_min,
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
