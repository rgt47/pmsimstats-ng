#' 03-study-B.R
#'
#' Paper 06 Study B: identifiability of the full nonlinear three-
#' component decomposition.
#'
#' Cell grid (per 00-ademp-pre-registration.md):
#'   m_PB in {1, 6}, m_TV in {0, 1}, N in {35, 70, 100, 150}
#'   = 16 cells.
#'   c_bm = 0.45 throughout.
#'
#' Method: full_nonlinear_fit only (currently a stub; see
#' 00-common.R for TODO). 1000 reps per cell.
#'
#' Performance: convergence rate, fit time, parameter-covariance
#' condition number, MSE of (m_BR_hat, m_PB_hat, m_TV_hat) vs the
#' participant's true m. Reported as a function of N to characterise
#' the sample-size threshold for usable identifiability.
#'
#' NOTE. Per-participant component estimates require an extension of
#' full_nonlinear_fit() to return participant-level estimates rather
#' than the population-level beta. The stub returns NA. Until the
#' nlme implementation lands, this script writes empty-output marker
#' files to keep the run-all orchestration intact and allows the user
#' to verify the script-level wiring.
#'
#' Seed: study_seed = 42 + 6*100 + 2 = 643.

source(here::here('analysis', 'scripts',
                  'component-decomposition', '00-common.R'))

set.seed(643)
STUDY_ID   <- 'B'
STUDY_SEED <- 643L

m_PB_levels <- c(1, 6)
m_TV_levels <- c(0, 1)
N_levels    <- c(35L, 70L, 100L, 150L)

cells <- CJ(m_PB = m_PB_levels, m_TV = m_TV_levels,
            N = N_levels, sorted = FALSE)
cells[, c_bm    := 0.45]
cells[, regime  := 'alt']
cells[, n_reps  := 1000L]
cells[, cell_id := .I]
cells[, analyses := list(list('full_nonlinear'))]

cat(sprintf('Study B: %d cells; total target reps = %d.\n',
            nrow(cells), sum(cells$n_reps)))

t_study <- Sys.time()
cell_summaries <- vector('list', nrow(cells))

for (i in seq_len(nrow(cells))) {
  cell <- list(m_PB = cells$m_PB[i], m_TV = cells$m_TV[i],
               N = cells$N[i], c_bm = cells$c_bm[i],
               analyses = cells$analyses[[i]],
               design = 'full')
  reps_dt <- run_cell_06(cell,
                         n_reps    = cells$n_reps[i],
                         study_seed = STUDY_SEED,
                         cell_id   = cells$cell_id[i])
  save_cell_06(reps_dt, STUDY_ID, cells$cell_id[i])

  tb <- TRUE_BETA(c_bm = cells$c_bm[i])
  by_an <- reps_dt[, summarise_cell_06(.SD, tb), by = analysis]
  by_an[, `:=`(cell_id = cells$cell_id[i],
               m_PB = cells$m_PB[i], m_TV = cells$m_TV[i],
               N = cells$N[i], c_bm = cells$c_bm[i],
               regime = cells$regime[i],
               true_beta = tb)]
  cell_summaries[[i]] <- by_an
  cat(sprintf('Study B cell %d/%d done.\n', i, nrow(cells)))
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
save_summary_06(study_summary, STUDY_ID)

cat(sprintf('Study B complete: %.1f min.\n',
            as.numeric(difftime(Sys.time(), t_study,
                                units = 'mins'))))
