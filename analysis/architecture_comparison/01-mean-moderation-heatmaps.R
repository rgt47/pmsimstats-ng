# 01-mean-moderation-heatmaps.R
#
# Architecture A (mean moderation) heatmaps: side by side comparison
# of analysis model with vs without carryover term (Dbc).
#
# Usage:
#   NREPS=50 Rscript analysis/architecture_comparison/01-mean-moderation-heatmaps.R

library(devtools)
library(tictoc)
load_all('.')
library(data.table)
library(gridExtra)

output_dir <- 'analysis/architecture_comparison/output'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv('NREPS', unset = '50'))
n_cores <- as.integer(Sys.getenv('NCORES', unset = '-1'))
if (n_cores < 0) n_cores <- max(1, parallel::detectCores() - 1)
cat(sprintf('Mean moderation heatmaps: Nreps=%d, cores=%d\n\n',
  Nreps, n_cores))

tdOL <- buildtrialdesign(
  'open label', 'OL',
  cumulative(rep(2.5, 8)),
  paste0('OL', 1:8), rep(1, 8),
  list(pathA = rep(1, 8)))

tdOLBDC <- buildtrialdesign(
  'open label + BDC', 'OL+BDC',
  c(4, 8, 12, 16, 17, 18, 19, 20),
  c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2', 'BD3', 'BD4'),
  c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)))

tdCO <- buildtrialdesign(
  'crossover', 'CO',
  cumulative(rep(2.5, 8)),
  c(paste0('COa', 1:4), paste0('COb', 1:4)),
  rep(0.5, 8),
  list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)))

tdNof1 <- buildtrialdesign(
  'N-of-1', 'N-of-1',
  c(4, 8, 9, 10, 11, 12, 16, 20),
  c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
  c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))

trialdesigns <- list(tdOL, tdOLBDC, tdCO, tdNof1)

data(extracted_bp)
data(extracted_rp)
blparamsets <- list(list(name = 'TR', param = extracted_bp))
respparamsets <- list(list(name = 'TR', param = extracted_rp))

censorparams <- data.table(
  names = c('balanced', 'more of flat', 'more of biased', 'high dropout'),
  beta0 = c(0.05, 0.05, 0.05, 0.15),
  beta1 = c(0.5, 0.2, 0.8, 0.5),
  eb1 = c(2, 2, 2, 2))

coremodelparams <- data.table(expand.grid(
  N = c(35, 70),
  c.bm = c(0, 0.25, 0.45),
  carryover_t1half = c(0, 0.5, 1),
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.2, c.cfct = 0.1))

simparam <- list(
  Nreps = Nreps, progressiveSave = FALSE,
  basesavename = 'archA', nRep2save = 5,
  saveunit2start = 1, savedir = output_dir)

ap_ignore <- data.table(
  useDE = FALSE, t_random_slope = FALSE,
  full_model_out = FALSE, simplecarryover = FALSE,
  carryover_t1half = 0, carryover_scalefactor = 1)

cat('=== Panel A: No carryover in analysis model ===\n')
tic('Panel A')
sim_a <- generateSimulatedResults(
  trialdesigns = trialdesigns,
  respparamsets = respparamsets,
  blparamsets = blparamsets,
  censorparams = censorparams,
  modelparams = coremodelparams,
  simparam = simparam,
  analysisparams = ap_ignore,
  rawdataout = FALSE,
  n_cores = n_cores,
  dgp_architecture = 'mean_moderation')
toc()

saveRDS(sim_a, file.path(output_dir,
  sprintf('archA_no_carryover_n%d.rds', Nreps)))

cat('\n=== Panel B: Matched Dbc carryover in analysis ===\n')
tic('Panel B')
sim_b_parts <- list()
for (dgp_th in c(0, 0.5, 1)) {
  cat(sprintf('  DGP t1half = %.1f\n', dgp_th))
  mp_sub <- coremodelparams[carryover_t1half == dgp_th]
  ap_dbc <- data.table(
    useDE = FALSE, t_random_slope = FALSE,
    full_model_out = FALSE, simplecarryover = FALSE,
    carryover_t1half = dgp_th, carryover_scalefactor = 1)
  sim_part <- generateSimulatedResults(
    trialdesigns = trialdesigns,
    respparamsets = respparamsets,
    blparamsets = blparamsets,
    censorparams = censorparams,
    modelparams = mp_sub,
    simparam = simparam,
    analysisparams = ap_dbc,
    rawdataout = FALSE,
    n_cores = n_cores,
    dgp_architecture = 'mean_moderation')
  sim_b_parts[[as.character(dgp_th)]] <- sim_part$results
}
toc()

sim_b <- list(
  results = rbindlist(sim_b_parts),
  parameterselections = sim_a$parameterselections)
saveRDS(sim_b, file.path(output_dir,
  sprintf('archA_matched_carryover_n%d.rds', Nreps)))

cat('\n=== Generating heatmaps ===\n')

provenance_a <- sprintf(
  'Arch A (mean moderation) | No carryover in analysis | %d reps | %s',
  Nreps, format(Sys.time(), '%Y-%m-%d'))
provenance_b <- sprintf(
  'Arch A (mean moderation) | Matched Dbc in analysis | %d reps | %s',
  Nreps, format(Sys.time(), '%Y-%m-%d'))

p1a <- PlotModelingResults(sim_a,
  param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'c.bm', 'censorparamset'),
  param2hold = data.table(
    blparamset = 1, respparamset = 1, carryover_t1half = 0),
  param2nothold = c('c.cfct', 'modelparamset')) +
  ggtitle('A1. No carryover in analysis (t1half=0)') +
  labs(subtitle = provenance_a) +
  theme(plot.title = element_text(hjust = 0, size = 9),
    plot.subtitle = element_text(size = 6, color = 'grey40'))

p1b <- PlotModelingResults(sim_b,
  param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'c.bm', 'censorparamset'),
  param2hold = data.table(
    blparamset = 1, respparamset = 1, carryover_t1half = 0),
  param2nothold = c('c.cfct', 'modelparamset')) +
  ggtitle('B1. Matched Dbc in analysis (t1half=0)') +
  labs(subtitle = provenance_b) +
  theme(plot.title = element_text(hjust = 0, size = 9),
    plot.subtitle = element_text(size = 6, color = 'grey40'))

p2a <- PlotModelingResults(sim_a,
  param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'carryover_t1half',
    'censorparamset'),
  param2hold = data.table(
    blparamset = 1, respparamset = 1, c.bm = 0.45),
  param2nothold = c('c.cfct', 'modelparamset')) +
  ggtitle('A2. No carryover in analysis (c.bm=0.45)') +
  labs(subtitle = provenance_a) +
  theme(plot.title = element_text(hjust = 0, size = 9),
    plot.subtitle = element_text(size = 6, color = 'grey40'))

p2b <- PlotModelingResults(sim_b,
  param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'carryover_t1half',
    'censorparamset'),
  param2hold = data.table(
    blparamset = 1, respparamset = 1, c.bm = 0.45),
  param2nothold = c('c.cfct', 'modelparamset')) +
  ggtitle('B2. Matched Dbc in analysis (c.bm=0.45)') +
  labs(subtitle = provenance_b) +
  theme(plot.title = element_text(hjust = 0, size = 9),
    plot.subtitle = element_text(size = 6, color = 'grey40'))

ml <- arrangeGrob(grobs = list(p1a, p1b, p2a, p2b),
  nrow = 2, ncol = 2)

timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
pdf_path <- file.path(output_dir,
  paste0('archA_heatmaps_', timestamp, '.pdf'))
ggsave(pdf_path, ml, width = 14, height = 12,
  units = 'in', dpi = 300)
cat('Saved:', pdf_path, '\n')

png_path <- file.path(output_dir,
  paste0('archA_heatmaps_', timestamp, '.png'))
ggsave(png_path, ml, width = 14, height = 12,
  units = 'in', dpi = 150)
cat('Saved:', png_path, '\n')

cat('\nDone.\n')
