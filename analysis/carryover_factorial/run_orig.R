# run_orig.R
# Factorial comparison: DGP carryover x analysis model
# Uses orig R/ package functions via devtools::load_all()
#
# Design: For each DGP half-life, loop over 3 analysis models
# explicitly. Each (dgp_t1half, analysis_model) pair runs a
# separate simulation so the analysis_model label is unambiguous.
# The same data cannot be reused across analysis models in this
# architecture, but at Nreps >= 200 the Monte Carlo error is
# small enough that this does not matter.
#
# Usage:
#   NREPS=10 Rscript analysis/carryover_factorial/run_orig.R
#   NREPS=1000 Rscript analysis/carryover_factorial/run_orig.R

library(devtools)
library(tictoc)
load_all('.')

output_dir <- 'analysis/carryover_factorial/output'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv('NREPS', unset = '10'))
cat(sprintf('Carryover factorial [orig]: Nreps = %d\n\n', Nreps))

# -----------------------------------------------------------------
# Trial designs
# -----------------------------------------------------------------

tdOL <- buildtrialdesign(
  name_longform = 'open label', name_shortform = 'OL',
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste0('OL', 1:8),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)

tdCO <- buildtrialdesign(
  name_longform = 'crossover', name_shortform = 'CO',
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = c(paste0('COa', 1:4), paste0('COb', 1:4)),
  expectancies = rep(0.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)

tdNof1 <- buildtrialdesign(
  name_longform = 'hybrid n-of-1', name_shortform = 'Nof1',
  timepoints = cumulative(c(4, 4, 1, 1, 1, 1, 4, 4)),
  timeptname = c(paste0('HYa', 1:4), paste0('HYb', 1:4)),
  expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 0, 1, 1, 0, 1, 0),
    pathD = c(1, 1, 0, 1, 1, 0, 0, 1)
  )
)

tdOLBDC <- buildtrialdesign(
  name_longform = 'open label + BDC', name_shortform = 'OL+BDC',
  timepoints = cumulative(c(4, 4, 4, 4, 1, 1, 1, 1)),
  timeptname = c(paste0('OBa', 1:4), paste0('OBb', 1:4)),
  expectancies = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

trialdesigns <- list(tdOL, tdCO, tdNof1, tdOLBDC)
design_names <- c('OL', 'CO', 'Nof1', 'OL+BDC')

# -----------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------

respparam <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)

blparam <- data.table(
  cat = c('bm', 'BL'),
  m = c(0, 70),
  sd = c(1, 10)
)

respparamsets <- list(list(name = 'default', param = respparam))
blparamsets <- list(list(name = 'default', param = blparam))

simparam <- list(
  Nreps = Nreps,
  progressiveSave = FALSE,
  saveunit2start = 1
)

# -----------------------------------------------------------------
# Factorial grid
# -----------------------------------------------------------------

dgp_halflives <- c(0, 0.25, 0.5, 1.0, 2.0)
cbm_values <- c(0, 0.3, 0.45)

make_analysisparams <- function(label, simplecarryover,
                                 analysis_t1half) {
  data.table(
    useDE = TRUE,
    t_random_slope = FALSE,
    full_model_out = FALSE,
    simplecarryover = simplecarryover,
    carryover_t1half = analysis_t1half,
    carryover_scalefactor = 1
  )
}

all_results <- list()
run_idx <- 0
total_runs <- length(dgp_halflives) * 3
current_run <- 0

for (dgp_t1half in dgp_halflives) {

  modelparams <- data.table(
    N = 70,
    c.bm = cbm_values,
    carryover_t1half = dgp_t1half,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )

  models <- list(
    list(label = 'ignore',
         ap = make_analysisparams('ignore', FALSE, 0)),
    list(label = 'separate_tsd',
         ap = make_analysisparams('separate_tsd', TRUE, 0)),
    list(label = 'continuous_Dbc',
         ap = make_analysisparams('continuous_Dbc', FALSE,
                                   dgp_t1half))
  )

  for (m in models) {
    current_run <- current_run + 1
    cat(sprintf('[%d/%d] DGP t1half=%.2f, analysis=%s\n',
                current_run, total_runs, dgp_t1half, m$label))

    result <- tryCatch(
      generateSimulatedResults(
        trialdesigns = trialdesigns,
        respparamsets = respparamsets,
        blparamsets = blparamsets,
        censorparams = NA,
        modelparams = modelparams,
        simparam = simparam,
        analysisparams = m$ap,
        lambda_cor = NA,
        n_cores = 1
      ),
      error = function(e) {
        cat(sprintf('  ERROR: %s\n', e$message))
        NULL
      }
    )

    if (!is.null(result)) {
      dt <- result$results
      dt[, dgp_t1half := dgp_t1half]
      dt[, analysis_model := m$label]
      dt[, design := design_names[trialdesign]]
      run_idx <- run_idx + 1
      all_results[[run_idx]] <- dt
    }
  }
}

# -----------------------------------------------------------------
# Combine and save
# -----------------------------------------------------------------

results <- rbindlist(all_results, fill = TRUE)

outfile <- file.path(output_dir,
  sprintf('factorial_orig_n%d.rds', Nreps))
saveRDS(results, outfile)
cat(sprintf('\nSaved %d rows to %s\n', nrow(results), outfile))

# -----------------------------------------------------------------
# Summary tables
# -----------------------------------------------------------------

cat('\n=== Type I error (c.bm = 0) ===\n')
t1e <- results[c.bm == 0,
  .(type1 = mean(p < 0.05, na.rm = TRUE),
    n = .N,
    n_na = sum(is.na(p))),
  by = .(design, dgp_t1half, analysis_model)]
print(t1e[order(design, dgp_t1half, analysis_model)])

cat('\n=== Power (c.bm = 0.3) ===\n')
pwr <- results[c.bm == 0.3,
  .(power = mean(p < 0.05, na.rm = TRUE),
    mean_beta = mean(beta, na.rm = TRUE),
    n = .N),
  by = .(design, dgp_t1half, analysis_model)]
print(pwr[order(design, dgp_t1half, analysis_model)])

cat('\n=== Power (c.bm = 0.45) ===\n')
pwr2 <- results[c.bm == 0.45,
  .(power = mean(p < 0.05, na.rm = TRUE),
    mean_beta = mean(beta, na.rm = TRUE),
    n = .N),
  by = .(design, dgp_t1half, analysis_model)]
print(pwr2[order(design, dgp_t1half, analysis_model)])

cat('\nDone.\n')
