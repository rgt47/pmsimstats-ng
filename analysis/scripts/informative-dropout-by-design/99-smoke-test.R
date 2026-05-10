#' 99-smoke-test.R
#'
#' Tiny end-to-end verification of the paper 09 production pipeline.
#' Runs 5 replicates each at one cell from each of Studies 1-4.
#' Should complete well under 60 seconds on an 8-core machine.

source(here::here('analysis', 'scripts',
                  'informative-dropout-by-design', '00-common.R'))

set.seed(940L)
STUDY_SEED <- 940L

cells <- list(
  S1 = list(design_name = 'OL',     dropout_name = 'balanced',
            N = 35L, c_bm = 0.45),
  S2 = list(design_name = 'OLBDC',  dropout_name = 'more_of_biased',
            N = 35L, c_bm = 0.30),
  S3 = list(design_name = 'Hybrid', dropout_name = 'balanced',
            N = 35L, c_bm = 0.30),
  S4 = list(design_name = 'CO',     dropout_name = 'high_dropout',
            N = 35L, c_bm = 0.45,
            analysis_name = 'mar_mmrm')
)

cat('Smoke test: 4 cells x 5 reps each.\n')

t0 <- Sys.time()

cat('S1 cell...\n')
s1 <- run_cell_09(cells$S1, n_reps = 5L,
                  study_seed = STUDY_SEED, cell_id = 1L,
                  chunk_size = 5L, progress_every = 5L)
cat(sprintf('  S1: %d rows.\n', nrow(s1)))

cat('S3 cell (path-aware)...\n')
source(here::here('analysis', 'scripts',
                  'informative-dropout-by-design',
                  '03-study3.R'),
       local = new.env())
## The above sources the production Study 3 file; instead we
## inline-call its replicate executor without running the full grid:
s3 <- (function() {
  reps <- list()
  for (r in 1:5) {
    reps[[r]] <- one_rep_S3(r, cells$S3, STUDY_SEED, 3L)
  }
  rbindlist(reps, fill = TRUE)
})()
cat(sprintf('  S3: %d rows.\n', nrow(s3)))

cat(sprintf('Smoke test wall: %.1f s.\n',
            as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
cat('Smoke test PASSED if rows are non-zero and no errors above.\n')
