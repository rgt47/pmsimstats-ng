## analysis/scripts/carryover-sensitivity/28-run-g9-prognostic-stress.R
##
## Stress-tests G9 (AIC-selected t1/2 + CR2) against prognostic-
## biomarker confounds beyond what has been checked so far. Every
## prior prognostic check (docs/64; job 21, 06-g-vs-f-comparison.R)
## used c_bm = 0 (confound only, no true interaction) at one
## confound strength (0.3). Two gaps follow directly from that and
## are this script's reason to exist:
##
##   1. Does a STRONGER confound (0.6, double the strength checked
##      so far) inflate G9's Type I error where 0.3 did not?
##   2. Does a confound, present ALONGSIDE a true interaction
##      (c_bm = 0.45), reduce G9's power relative to the no-confound
##      alt baseline? This co-occurring scenario has never been run
##      (docs/64's own caveat: "whether a prognostic confound and
##      the true interaction co-occur, is untested").
##
## All four of paper 01 Study 3's prognostic pathways (prog-BL,
## prog-TV-const, prog-TV-grow, prog-PB) are crossed with both
## confound strengths, both c_bm levels, and N in {35, 70}, Hybrid
## design (where the F-family is identifiable, for a direct
## comparison), under the corrected c.cf1t/c.cfct = 0.2 DGP
## constants (already the default in simulate_cell_combined()).
## Fits all nine G-family specifications, E0, and the F-family
## (F1/F1A/F1A0/F1A0B) to the same simulated trial each replicate,
## via simulate_cell_combined() (simulation-core.R) -- the same
## machinery as 06-g-vs-f-comparison.R and 27-run-tier1-g-vs-f-e0.R.
##
## The pathway x N x c_bm = 0 x confound = 0.3 x Hybrid cells
## duplicate four of job 21's cells; they are re-run here (not
## merged from job 21's output) to keep this driver self-contained.
## That is 4 of this grid's 32 cells (12.5%), an acceptable
## duplication cost for a standalone, directly reproducible run.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/28-run-g9-prognostic-stress.R \
##     [--reps R] [--n-boot B] [--cores K] [--out-dir PATH]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
REPS <- as.integer(argv('--reps', '500'))
N_BOOT <- as.integer(argv('--n-boot', '1'))
CORES <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
OUT_DIR <- argv('--out-dir',
  'analysis/data/carryover-sensitivity/g9-prognostic-stress')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

base_cell <- list(design = 'Hybrid', t1half = 0.5,
                  carryover_form = 'exponential', weibull_shape = 1,
                  dgp_arch = 'mvn')

## Each pathway sets exactly one of c_bl/c_tv/b_tv/c_pb (docs/57
## Section 7.3); the other three stay at their default of 0.
pathway_field <- c(progbl = 'c_bl', progtvconst = 'c_tv',
                   progtvgrow = 'b_tv', progpb = 'c_pb')

grid <- tidyr::expand_grid(
  pathway  = names(pathway_field),
  strength = c(0.3, 0.6),
  c_bm     = c(0, 0.45),
  N        = c(35, 70)
) |>
  dplyr::mutate(
    cell_id = sprintf('%s-conf%.1f-cbm%.2f-N%d', pathway, strength, c_bm, N)
  )

build_cell <- function(i) {
  cl <- modifyList(base_cell, list(N = grid$N[i], c_bm = grid$c_bm[i]))
  cl[[pathway_field[[grid$pathway[i]]]]] <- grid$strength[i]
  cl
}

cat(sprintf(
  'G9 prognostic stress: %d cells (%d pathways x %d strengths x %d c_bm x %d N), n_reps = %d\n',
  nrow(grid), length(pathway_field), 2, 2, 2, REPS))

jobs <- tidyr::expand_grid(cell_i = seq_len(nrow(grid)), rep = seq_len(REPS))

one_job <- function(i) {
  cell_i <- jobs$cell_i[i]; r <- jobs$rep[i]
  cell_id <- grid$cell_id[cell_i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(25102600L + r)
  res <- tryCatch(
    simulate_cell_combined(build_cell(cell_i), n_reps = 1, n_boot = N_BOOT),
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

## ---------------- summary: nominal rejection rate per method, plus warp-speed F1A0B ----------------
summarize_cell <- function(d) {
  nominal <- d |>
    filter(!(family == 'F' & method == 'F1A0')) |>
    group_by(family, method) |>
    summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
              n = sum(!is.na(p_value)), .groups = 'drop')
  boot_pool <- d |> filter(family == 'F', method == 'F1A0') |>
    distinct(rep, boot_t_f1a0) |> pull(boot_t_f1a0)
  obs_f1a0 <- d |> filter(family == 'F', method == 'F1A0') |> pull(t)
  f1a0_nominal <- d |> filter(family == 'F', method == 'F1A0') |>
    summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
              n = sum(!is.na(p_value)))
  bind_rows(
    nominal,
    tibble(family = 'F', method = 'F1A0', reject = f1a0_nominal$reject,
          n = f1a0_nominal$n),
    tibble(family = 'F', method = 'F1A0B',
          reject = ws_size(abs(obs_f1a0), boot_pool),
          n = sum(!is.na(obs_f1a0)))
  )
}

summary_tab <- res |>
  group_by(cell_id, pathway, strength, c_bm, N) |>
  group_modify(~ summarize_cell(.x)) |>
  ungroup()

cat('\nRejection rate (Type I error at c_bm=0 cells, power at c_bm=0.45 cells):\n')
print(as.data.frame(summary_tab), max = 2000)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))

meta <- list(
  script = '28-run-g9-prognostic-stress.R',
  date = Sys.time(),
  n_reps = REPS,
  n_boot = N_BOOT,
  elapsed_min_this_run = elapsed_min,
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
