#' 06-smoke-test.R
#'
#' Tiny end-to-end verification of the paper 06 production pipeline.
#' Runs 5 replicates each at four corner cells of Study A. Should
#' complete in well under 60 seconds on an 8-core machine.

source(here::here('analysis', 'scripts',
                  'component-decomposition', '00-common.R'))

set.seed(640L)
STUDY_SEED <- 640L

corners <- data.table(
  m_PB    = c(0,  0, 10, 10),
  m_TV    = c(-1, 2, -1,  2),
  N       = c(35L, 35L, 35L, 35L),
  c_bm    = 0.45,
  n_reps  = 5L,
  cell_id = 1:4)

cat('Smoke test: 4 corner cells x 5 reps each.\n')

t0 <- Sys.time()
for (i in seq_len(nrow(corners))) {
  cell <- list(m_PB = corners$m_PB[i], m_TV = corners$m_TV[i],
               N = corners$N[i], c_bm = corners$c_bm[i],
               analyses = c('one_component', 'phase_augmented'),
               design = 'simple')
  reps_dt <- run_cell_06(cell,
                         n_reps    = corners$n_reps[i],
                         study_seed = STUDY_SEED,
                         cell_id   = corners$cell_id[i],
                         chunk_size = 5L,
                         progress_every = 5L)
  cat(sprintf('  Cell %d returned %d rows.\n', i, nrow(reps_dt)))
}
cat(sprintf('Smoke test wall: %.1f s.\n',
            as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
cat('Smoke test PASSED if all four cells returned non-zero rows.\n')
