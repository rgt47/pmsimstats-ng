## Evaluation of F1A0B (and its F1/F1A/F1A0 relatives) against the
## docs/62 "what best would mean" criteria, on the refined-process
## DGP, scoped to three cells (docs/62; this session, 2026-10-22):
##   hybrid   reference design, null (c_bm = 0) and alt (c_bm = 0.45)
##   hybridlo same as hybrid_alt but c_bm = 0.25, to de-saturate the
##            power comparison among F1/F1A/F1A0/F1A0B: at 0.45 all
##            four are near the power ceiling (0.94-0.96) and cannot
##            be told apart (docs/64). Added 2026-10-22.
##   olbdc    design-robustness arm (CO's constant De made it
##            unidentifiable for the F-family model; see
##            simulation-core.R's F-family block), null and alt
##   progbl   prognostic-biomarker-null cell (c_bm = 0, c_bl = 0.3;
##            docs/57 Section 7.3 prog-BL), null only -- this cell
##            tests Type I error under a confound unrelated to the
##            interaction, not power
##   progtvconst, progtvgrow, progpb  the other three docs/57
##            Section 7.3 prognostic pathways (c_tv = 0.3, b_tv =
##            0.3, c_pb = 0.3 respectively), null only, same
##            rationale as progbl. Added 2026-10-22 once c.bm.pb was
##            mirrored into implementations/tidyverse and plumbed
##            through this driver (previously package-only).
##
## Flattens all (cell, replicate) jobs into one list so the 8
## available cores are shared across cells rather than one cell
## finishing before the next starts. Each job is checkpointed to its
## own .rds file, so an interrupted run resumes without recomputing
## finished replicates.
##
## Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/05-f1a0b-evaluation.R \
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
CORES <- as.integer(argv('--cores', '6'))
OUT_DIR <- argv('--out-dir', 'analysis/data/carryover-sensitivity/f1a0b-evaluation')
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
  progpb_null = modifyList(base_cell, list(design = 'Hybrid', c_bm = 0, c_pb = 0.3))
)

jobs <- expand.grid(cell = names(CELLS), rep = seq_len(REPS),
                    stringsAsFactors = FALSE)
one_job <- function(i) {
  cell_name <- jobs$cell[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_name, r))
  if (file.exists(f)) return(invisible(NULL))
  set.seed(22102600L + r)
  res <- tryCatch(
    simulate_cell_f(CELLS[[cell_name]], n_reps = 1, n_boot = N_BOOT),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    ## simulate_cell_f()'s own `rep` column is its internal loop
    ## index (always 1 for a single-replicate call here); overwrite
    ## with this job's actual replicate number.
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

## ---------------- summary: nominal F1/F1A/F1A0, warp-speed F1A0B ----------------
summarize_cell <- function(d) {
  nominal <- d |>
    filter(method %in% c('F1', 'F1A', 'F1A0')) |>
    group_by(method) |>
    summarise(reject = mean(p_value < 0.05, na.rm = TRUE),
              n = sum(!is.na(p_value)), .groups = 'drop')
  boot_pool <- d |> filter(method == 'F1A0') |> distinct(rep, boot_t_f1a0) |>
    pull(boot_t_f1a0)
  obs_f1a0 <- d |> filter(method == 'F1A0') |> pull(t)
  bind_rows(
    nominal,
    tibble(method = 'F1A0B',
          reject = ws_size(abs(obs_f1a0), boot_pool),
          n = sum(!is.na(obs_f1a0)))
  )
}

summary_tab <- res |>
  group_by(cell) |>
  group_modify(~ summarize_cell(.x)) |>
  ungroup()

cat('\nRejection rate (Type I error at null cells, power at alt cells):\n')
print(summary_tab, n = 30)
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))
