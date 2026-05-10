#' 02-study-A.R
#'
#' Paper 06 Study A: bias of the one-component analysis under the
#' three-component DGP.
#'
#' Cell grid (per 00-ademp-pre-registration.md):
#'   m_PB     in {0, 1, 3, 6, 10}     (5 levels)
#'   m_TV     in {-1, 0, 1, 2}        (4 levels)
#'   N        in {35, 70, 100, 150}   (4 levels)
#'   analyses x cells = 240 (one_component, phase_augmented at all
#'                          cells; full_nonlinear at N >= 100 only)
#'
#' Replicates: 1000 per cell for the alternative cells; 5000 for the
#' c_bm = 0 type-I anchor row (one cell per N at the m_PB = 0,
#' m_TV = 0 reference). Type-I cells use c_bm = 0; alternative cells
#' use c_bm = 0.45 per the prazosin-PTSD baseline.
#'
#' Seed: study_seed = 42 + 6*100 + 1 = 642.

source(here::here('analysis', 'scripts',
                  'component-decomposition', '00-common.R'))

set.seed(642)
STUDY_ID   <- 'A'
STUDY_SEED <- 642L

m_PB_levels <- c(0, 1, 3, 6, 10)
m_TV_levels <- c(-1, 0, 1, 2)
N_levels    <- c(35L, 70L, 100L, 150L)

#=====================================================================
# Cell construction
#=====================================================================

alt_cells <- CJ(m_PB = m_PB_levels, m_TV = m_TV_levels,
                N = N_levels, sorted = FALSE)
alt_cells[, c_bm    := 0.45]
alt_cells[, regime  := 'alt']
alt_cells[, n_reps  := 1000L]

null_cells <- data.table(
  m_PB = 0, m_TV = 0, N = N_levels,
  c_bm = 0, regime = 'null', n_reps = 5000L)

cells <- rbindlist(list(alt_cells, null_cells), use.names = TRUE,
                   fill = TRUE)
cells[, cell_id := .I]
cells[, analyses := list(list(c('one_component', 'phase_augmented')))]
## Add full_nonlinear at N >= 100 (per pre-reg).
for (i in seq_len(nrow(cells))) {
  if (cells$N[i] >= 100) {
    cells$analyses[[i]] <- c(cells$analyses[[i]],
                             'full_nonlinear')
  }
}

cat(sprintf('Study A: %d cells (%d alt + %d null); ',
            nrow(cells), nrow(alt_cells), nrow(null_cells)))
cat(sprintf('total target reps = %d.\n', sum(cells$n_reps)))

#=====================================================================
# Execute
#=====================================================================

t_study <- Sys.time()
cell_summaries <- vector('list', nrow(cells))

for (i in seq_len(nrow(cells))) {
  cell <- list(m_PB = cells$m_PB[i], m_TV = cells$m_TV[i],
               N = cells$N[i], c_bm = cells$c_bm[i],
               analyses = cells$analyses[[i]],
               design = if (cells$N[i] >= 100) 'full' else 'simple')
  reps_dt <- run_cell_06(cell,
                         n_reps    = cells$n_reps[i],
                         study_seed = STUDY_SEED,
                         cell_id   = cells$cell_id[i])
  reps_dt[, regime := cells$regime[i]]
  save_cell_06(reps_dt, STUDY_ID, cells$cell_id[i])

  ## Per-analysis summary for this cell
  tb <- if (cells$c_bm[i] == 0) 0 else TRUE_BETA(c_bm = cells$c_bm[i])
  by_an <- reps_dt[, summarise_cell_06(.SD, tb), by = analysis]
  by_an[, `:=`(cell_id = cells$cell_id[i],
               m_PB = cells$m_PB[i], m_TV = cells$m_TV[i],
               N = cells$N[i], c_bm = cells$c_bm[i],
               regime = cells$regime[i],
               true_beta = tb)]
  cell_summaries[[i]] <- by_an

  cat(sprintf('Study A cell %d/%d done.\n', i, nrow(cells)))
}

study_summary <- rbindlist(cell_summaries, use.names = TRUE,
                           fill = TRUE)
setcolorder(study_summary,
            c('cell_id', 'regime', 'm_PB', 'm_TV', 'N', 'c_bm',
              'analysis', 'true_beta',
              'n_reps', 'n_converged', 'conv_rate',
              'power', 'mcse_power',
              'mean_beta', 'sd_beta', 'mcse_mean_beta',
              'bias', 'mcse_bias',
              'coverage', 'mcse_coverage'))
setorder(study_summary, regime, N, m_TV, m_PB, analysis)

save_summary_06(study_summary, STUDY_ID)

cat(sprintf('Study A complete: %.1f min.\n',
            as.numeric(difftime(Sys.time(), t_study,
                                units = 'mins'))))
