# test-pd-08-current.R
#
# Test positive definiteness with c.tv=c.pb=c.br=0.8
# using current HEAD code (Architecture A + AR(1))
#
# Usage:
#   Rscript analysis/figure4/test-pd-08-current.R

library(devtools)
load_all('.')
library(data.table)
library(corpcor)
library(tictoc)

output_dir <- 'analysis/figure4/output'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv('NREPS', unset = '50'))
n_cores <- as.integer(Sys.getenv('NCORES', unset = '-1'))
if (n_cores < 0) n_cores <- max(1, parallel::detectCores() - 1)

cat('=== PD Test: c.tv=c.pb=c.br=0.8 (Current HEAD) ===\n')
cat('Architecture: mean moderation (A)\n')
cat('Covariance: AR(1)\n')
cat('Nreps:', Nreps, '\nCores:', n_cores, '\n\n')

# Instrumentation
pd_log <- new.env(parent = emptyenv())
pd_log$checks <- 0L
pd_log$failures <- 0L
pd_log$fixes <- 0L
pd_log$eigenvalues_at_failure <- list()

original_is_pd <- corpcor::is.positive.definite
original_make_pd <- corpcor::make.positive.definite

instrumented_is_pd <- function(m, tol = 1e-8, method = c('eigen', 'chol')) {
  pd_log$checks <<- pd_log$checks + 1L
  result <- original_is_pd(m, tol = tol)

  if (!result) {
    pd_log$failures <<- pd_log$failures + 1L
    evals <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
    pd_log$eigenvalues_at_failure[[
      length(pd_log$eigenvalues_at_failure) + 1]] <<- evals
  }
  result
}

instrumented_make_pd <- function(m, tol = 1e-3) {
  pd_log$fixes <<- pd_log$fixes + 1L
  original_make_pd(m, tol = tol)
}

environment(instrumented_is_pd) <- asNamespace('corpcor')
environment(instrumented_make_pd) <- asNamespace('corpcor')
assignInNamespace('is.positive.definite', instrumented_is_pd, ns = 'corpcor')
assignInNamespace('make.positive.definite', instrumented_make_pd, ns = 'corpcor')

# Trial designs
tdOL <- buildtrialdesign(
  name_longform = 'open label',
  name_shortform = 'OL',
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste('OL', 1:8, sep = ''),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)

tdOLBDC <- buildtrialdesign(
  name_longform = 'open label+blinded discontinuation',
  name_shortform = 'OL+BDC',
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptname = c('OL1', 'OL2', 'OL3', 'OL4',
                 'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

tdCO <- buildtrialdesign(
  name_longform = 'traditional crossover',
  name_shortform = 'CO',
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = c(paste('COa', 1:4, sep = ''),
                 paste('COb', 1:4, sep = '')),
  expectancies = rep(.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)

tdNof1 <- buildtrialdesign(
  name_longform = 'primary N-of-1 design',
  name_shortform = 'N-of-1',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptname = c('OL1', 'OL2', 'BD1', 'BD2',
                 'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)

trialdesigns <- list(OL = tdOL, OLBDC = tdOLBDC,
                     CO = tdCO, Nof1 = tdNof1)

# Parameters
data(extracted_bp)
data(extracted_rp)
blparamsets <- list(list(name = 'TR', param = extracted_bp))
respparamsets <- list(list(name = 'TR', param = extracted_rp))

censorparams <- data.table(
  names = c('balanced', 'more of flat',
            'more of biased', 'high dropout'),
  beta0 = c(.05, .05, .05, .15),
  beta1 = c(.5, .2, .8, .5),
  eb1 = c(2, 2, 2, 2)
)

# Test with 0.8 correlations (Hendrickson values)
coremodelparams <- expand.grid(
  N = c(35, 70),
  c.bm = c(0, .25, .45),
  carryover_t1half = c(0, .5, 1.0),
  c.tv = .8, c.pb = .8, c.br = .8,
  c.cf1t = .2, c.cfct = .1
)

analysisparams <- expand.grid(
  useDE = FALSE,
  t_random_slope = FALSE,
  full_model_out = FALSE
)

cat('Starting simulation at', format(Sys.time()), '\n\n')
tic('Simulation')

simresults <- generateSimulatedResults(
  trialdesigns = list(trialdesigns$OL, trialdesigns$OLBDC,
                      trialdesigns$CO, trialdesigns$Nof1),
  respparamsets = respparamsets,
  blparamsets = blparamsets,
  censorparams = censorparams,
  modelparams = coremodelparams,
  simparam = list(
    Nreps = Nreps,
    progressiveSave = FALSE,
    basesavename = 'pd_test_08_HEAD',
    nRep2save = 5,
    saveunit2start = 1,
    savedir = output_dir
  ),
  analysisparams = analysisparams,
  rawdataout = FALSE,
  n_cores = n_cores,
  dgp_architecture = 'mean_moderation'
)

toc()

# Report
cat('\n')
cat(strrep('=', 65), '\n')
cat('POSITIVE DEFINITENESS DIAGNOSTICS: c.tv=c.pb=c.br=0.8\n')
cat(strrep('=', 65), '\n\n')

cat('Total is.positive.definite() checks:', pd_log$checks, '\n')
cat('Failures (non-PD matrices):         ',
    pd_log$failures, '\n')
cat('make.positive.definite() calls:     ', pd_log$fixes, '\n')
cat('Failure rate:                        ',
    sprintf('%.4f (%d / %d)\n',
            pd_log$failures / max(pd_log$checks, 1),
            pd_log$failures, pd_log$checks))

if (length(pd_log$eigenvalues_at_failure) > 0) {
  cat('\nEigenvalue analysis of non-PD matrices:\n')
  all_min_evals <- sapply(pd_log$eigenvalues_at_failure, min)
  cat('  Count of failures:', length(all_min_evals), '\n')
  cat('  Min of min eigenvalues:   ', sprintf('%.6f\n', min(all_min_evals)))
  cat('  Median of min eigenvalues:', sprintf('%.6f\n', median(all_min_evals)))
  cat('  Max of min eigenvalues:   ', sprintf('%.6f\n', max(all_min_evals)))
}

cat('\nSimulation complete at', format(Sys.time()), '\n')

# Save results
diag_result <- list(
  correlation_value = 0.8,
  architecture = 'mean_moderation',
  covariance_structure = 'AR(1)',
  nreps = Nreps,
  checks = pd_log$checks,
  failures = pd_log$failures,
  fixes = pd_log$fixes,
  failure_rate = pd_log$failures / max(pd_log$checks, 1),
  eigenvalues_at_failure = pd_log$eigenvalues_at_failure
)

outfile <- file.path(output_dir,
  paste0('pd_test_08_HEAD_n', Nreps, '.rds'))
saveRDS(diag_result, outfile)
cat('\nDiagnostics saved to:', outfile, '\n')

# Restore
assignInNamespace('is.positive.definite', original_is_pd, ns = 'corpcor')
assignInNamespace('make.positive.definite', original_make_pd, ns = 'corpcor')
