#' 01-study1.R
#'
#' Paper 09 Study 1: power x design x dropout pattern x N x c_bm.
#'
#' Cell grid (per 00-ademp-pre-registration.md):
#'   designs:  OL, OLBDC, CO, Hybrid          (4)
#'   dropouts: none, balanced, more_of_flat,
#'             more_of_biased, high_dropout   (5)
#'   N:        35, 70, 100                    (3)
#'   c_bm:     0, 0.3, 0.6                    (3)
#'   total = 4 * 5 * 3 * 3 = 180 cells.
#'
#' Replicates: 1000 per alternative cell (c_bm > 0); 5000 per type-I
#' cell (c_bm = 0).
#'
#' Seed: study_seed = 42 + 9 * 100 + 1 = 943.

source(here::here('analysis', 'scripts',
                  'informative-dropout-by-design', '00-common.R'))

set.seed(943)
STUDY_ID   <- '1'
STUDY_SEED <- 943L

design_levels  <- c('OL', 'OLBDC', 'CO', 'Hybrid')
dropout_levels <- names(DROPOUT_PATTERNS)
N_levels       <- c(35L, 70L, 100L)
c_bm_levels    <- c(0, 0.3, 0.6)

cells <- CJ(design_name  = design_levels,
            dropout_name = dropout_levels,
            N            = N_levels,
            c_bm         = c_bm_levels,
            sorted = FALSE)
cells[, n_reps  := ifelse(c_bm == 0, 5000L, 1000L)]
cells[, regime  := ifelse(c_bm == 0, 'null', 'alt')]
cells[, cell_id := .I]

cat(sprintf('Study 1: %d cells (%d alt + %d null); ',
            nrow(cells),
            sum(cells$regime == 'alt'),
            sum(cells$regime == 'null')))
cat(sprintf('total target reps = %d.\n', sum(cells$n_reps)))

t_study <- Sys.time()
cell_summaries <- vector('list', nrow(cells))

for (i in seq_len(nrow(cells))) {
  cell <- list(design_name  = cells$design_name[i],
               dropout_name = cells$dropout_name[i],
               N            = cells$N[i],
               c_bm         = cells$c_bm[i])
  reps_dt <- run_cell_09(cell,
                         n_reps    = cells$n_reps[i],
                         study_seed = STUDY_SEED,
                         cell_id   = cells$cell_id[i])
  reps_dt[, regime := cells$regime[i]]
  save_cell_09(reps_dt, STUDY_ID, cells$cell_id[i])

  tb <- TRUE_BETA_09(c_bm = cells$c_bm[i])
  s  <- summarise_cell_09(reps_dt, tb)
  s[, `:=`(cell_id = cells$cell_id[i],
           design_name = cells$design_name[i],
           dropout_name = cells$dropout_name[i],
           N = cells$N[i], c_bm = cells$c_bm[i],
           regime = cells$regime[i],
           true_beta = tb)]
  cell_summaries[[i]] <- s

  cat(sprintf('Study 1 cell %d/%d done.\n', i, nrow(cells)))
}

study_summary <- rbindlist(cell_summaries, use.names = TRUE,
                           fill = TRUE)
setcolorder(study_summary,
            c('cell_id', 'regime', 'design_name', 'dropout_name',
              'N', 'c_bm', 'true_beta',
              'n_reps', 'n_converged', 'conv_rate',
              'power', 'mcse_power',
              'mean_beta', 'sd_beta', 'mcse_mean_beta',
              'bias', 'mcse_bias',
              'coverage', 'mcse_coverage',
              'mean_dropout'))
setorder(study_summary, regime, design_name, dropout_name, N, c_bm)
save_summary_09(study_summary, STUDY_ID)

cat(sprintf('Study 1 complete: %.1f min.\n',
            as.numeric(difftime(Sys.time(), t_study,
                                units = 'mins'))))
