# 01-compare-dgp-architectures.R
#
# Compares Architecture A (mean moderation) vs Architecture B (MVN)
# across carryover conditions, with and without carryover in the
# analysis model.
#
# Usage:
#   NREPS=50 Rscript analysis/architecture_comparison/01-compare-dgp-architectures.R

library(devtools)
library(tictoc)
load_all('.')
library(data.table)

output_dir <- 'analysis/architecture_comparison/output'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv('NREPS', unset = '50'))
cat(sprintf('Architecture comparison: Nreps=%d\n\n', Nreps))

tdOL <- buildtrialdesign(
  'open label', 'OL',
  cumulative(rep(2.5, 8)),
  paste0('OL', 1:8),
  rep(1, 8),
  list(pathA = rep(1, 8)))

tdCO <- buildtrialdesign(
  'crossover', 'CO',
  cumulative(rep(2.5, 8)),
  c(paste0('COa', 1:4), paste0('COb', 1:4)),
  rep(0.5, 8),
  list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)))

tdNof1 <- buildtrialdesign(
  'hybrid n-of-1', 'Nof1',
  c(4, 8, 9, 10, 11, 12, 16, 20),
  c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
  c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))

tdOLBDC <- buildtrialdesign(
  'open label + BDC', 'OL+BDC',
  c(4, 8, 12, 16, 17, 18, 19, 20),
  c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2', 'BD3', 'BD4'),
  c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)))

trialdesigns <- list(tdOL, tdCO, tdNof1, tdOLBDC)
design_names <- c('OL', 'CO', 'Nof1', 'OL+BDC')

data(extracted_bp)
data(extracted_rp)
blparamsets <- list(list(name = 'TR', param = extracted_bp))
respparamsets <- list(list(name = 'TR', param = extracted_rp))

simparam <- list(
  Nreps = Nreps, progressiveSave = FALSE, saveunit2start = 1)

dgp_halflives <- c(0, 0.5, 1.0)
cbm_values <- c(0, 0.25, 0.45)
architectures <- c('mvn', 'mean_moderation')

make_ap <- function(carryover_t1half) {
  data.table(
    useDE = FALSE, t_random_slope = FALSE,
    full_model_out = FALSE, simplecarryover = FALSE,
    carryover_t1half = carryover_t1half,
    carryover_scalefactor = 1)
}

all_results <- list()
run_idx <- 0
total_runs <- length(architectures) * length(dgp_halflives) * 2
current_run <- 0

for (arch in architectures) {
  for (dgp_t1half in dgp_halflives) {

    modelparams <- data.table(
      N = 70, c.bm = cbm_values,
      carryover_t1half = dgp_t1half,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.1, c.cfct = 0.05)

    analysis_models <- list(
      list(label = 'ignore', ap = make_ap(0)),
      list(label = 'matched_Dbc', ap = make_ap(dgp_t1half)))

    for (m in analysis_models) {
      current_run <- current_run + 1
      cat(sprintf(
        '[%d/%d] arch=%s, dgp_t1half=%.1f, analysis=%s\n',
        current_run, total_runs, arch, dgp_t1half, m$label))

      tic(sprintf('  %s / t1half=%.1f / %s', arch,
        dgp_t1half, m$label))
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
          n_cores = 1,
          dgp_architecture = arch),
        error = function(e) {
          cat(sprintf('  ERROR: %s\n', e$message))
          NULL
        })
      toc()

      if (!is.null(result)) {
        dt <- result$results
        dt[, dgp_architecture := arch]
        dt[, dgp_t1half := dgp_t1half]
        dt[, analysis_model := m$label]
        dt[, design := design_names[trialdesign]]
        run_idx <- run_idx + 1
        all_results[[run_idx]] <- dt
      }
    }
  }
}

results <- rbindlist(all_results, fill = TRUE)
outfile <- file.path(output_dir,
  sprintf('architecture_comparison_n%d.rds', Nreps))
saveRDS(results, outfile)
cat(sprintf('\nSaved %d rows to %s\n\n', nrow(results), outfile))

cat('=== Power by architecture, design, carryover, analysis model ===\n')
cat('    (c.bm = 0.45)\n\n')
pwr <- results[c.bm == 0.45,
  .(power = mean(p < 0.05, na.rm = TRUE),
    mean_beta = mean(beta, na.rm = TRUE),
    n = .N, n_na = sum(is.na(p))),
  by = .(dgp_architecture, design, dgp_t1half, analysis_model)]
print(
  pwr[order(dgp_architecture, design, dgp_t1half, analysis_model)],
  n = 100)

cat('\n=== Type I error (c.bm = 0) ===\n\n')
t1e <- results[c.bm == 0,
  .(type1 = mean(p < 0.05, na.rm = TRUE),
    n = .N, n_na = sum(is.na(p))),
  by = .(dgp_architecture, design, dgp_t1half, analysis_model)]
print(
  t1e[order(dgp_architecture, design, dgp_t1half, analysis_model)],
  n = 100)

cat('\nDone.\n')
