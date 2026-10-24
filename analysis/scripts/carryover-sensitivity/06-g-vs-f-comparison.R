## G-family (G1-G9) vs F-family (F1/F1A/F1A0/F1A0B) head-to-head on
## shared cells, via simulate_cell_combined() (simulation-core.R):
## one DGP draw per replicate, both families fit to the *same*
## simulated trial. docs/62 Part 3 names this comparison as missing
## entirely and names criterion 7 (prognostic biomarker) as untested
## for every G-family specification. 2026-10-23.
##
## Cells reused from 05-f1a0b-evaluation.R's F-only evaluation
## (docs/64), restricted to the ones meaningful for both families:
##   hybrid, hybridlo, olbdc      reference comparison and the
##            de-saturated Hybrid power cell
##   progbl, progtvconst, progtvgrow, progpb   docs/57 Section 7.3
##            prognostic pathways, null only -- closes criterion 7
##            for the G-family, untested anywhere until now
##   t1hifhi, shapemis, meanmod   criteria 3, 4, 6, null and alt
## The F-only grid-density/small-N/off-grid stress cells
## (densegrid, vdensegrid, smalln, offgrid, smalln_offgrid_dense)
## are NOT included here: they test F1A0's own AIC grid, which has
## no G-family analog (G4/G9's own AIC grid is fixed at paper 02's
## 4-point {0.25, 0.5, 1, 2} and not parameterized by this change).
##
## Run from the repository root (SIM_CORES is set by simrun; a bare
## local run falls back to 6):
##   Rscript analysis/scripts/carryover-sensitivity/06-g-vs-f-comparison.R \
##     [--reps R] [--n-boot B] [--cores K] [--out-dir PATH]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
REPS <- as.integer(argv('--reps', '500'))
N_BOOT <- as.integer(argv('--n-boot', '1'))
CORES <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '6')))
OUT_DIR <- argv('--out-dir', 'analysis/data/carryover-sensitivity/g-vs-f-comparison')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)

base_cell <- list(N = 70, t1half = 0.5, carryover_form = 'exponential',
                  weibull_shape = 1, dgp_arch = 'mvn')

CELLS <- list(
  hybrid_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0)),
  hybrid_alt  = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0.45)),
  hybrid_alt_lo = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0.25)),
  olbdc_null  = modifyList(base_cell, list(design = 'OLBDC',  c_bm = 0)),
  olbdc_alt   = modifyList(base_cell, list(design = 'OLBDC',  c_bm = 0.45)),
  progbl_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0, c_bl = 0.3)),
  progtvconst_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0, c_tv = 0.3)),
  progtvgrow_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0, b_tv = 0.3)),
  progpb_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0, c_pb = 0.3)),
  t1hifhi_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0,
                                             t1half = 2)),
  t1hifhi_alt  = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0.25,
                                             t1half = 2)),
  shapemis_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0,
                                              carryover_form = 'weibull',
                                              weibull_shape = 2,
                                              analysis_form = 'exponential',
                                              analysis_shape = 1)),
  shapemis_alt  = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0.25,
                                              carryover_form = 'weibull',
                                              weibull_shape = 2,
                                              analysis_form = 'exponential',
                                              analysis_shape = 1)),
  meanmod_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0,
                                             dgp_arch = 'mean_moderation')),
  meanmod_alt  = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0.25,
                                             dgp_arch = 'mean_moderation'))
)

jobs <- expand.grid(cell = names(CELLS), rep = seq_len(REPS),
                    stringsAsFactors = FALSE)
one_job <- function(i) {
  cell_name <- jobs$cell[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_name, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(23102600L + r)
  res <- tryCatch(
    simulate_cell_combined(CELLS[[cell_name]], n_reps = 1, n_boot = N_BOOT),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    res <- mutate(res, rep = r, cell = cell_name, .before = 1)
    saveRDS(res, f)
  }
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job,
                   mc.cores = CORES, mc.preschedule = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

res <- bind_rows(lapply(list.files(ck_dir, full.names = TRUE), readRDS))
saveRDS(res, file.path(OUT_DIR, 'all.rds'))
n_done <- length(unique(paste(res$cell, res$rep)))
cat(sprintf('\n%d of %d (cell, rep) jobs present in %.1f min (this run; checkpointed jobs from earlier runs are not re-timed)\n',
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
  group_by(cell) |>
  group_modify(~ summarize_cell(.x)) |>
  ungroup()

cat('\nRejection rate (Type I error at null cells, power at alt cells):\n')
print(summary_tab, n = 60)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))
