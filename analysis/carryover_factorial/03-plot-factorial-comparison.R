# 03-plot-factorial-comparison.R
# Figure 4 heatmaps: (a) no carryover in analysis model vs
#                     (b) continuous Dbc carryover in analysis model
#
# Usage:
#   NREPS=50 Rscript analysis/carryover_factorial/03-plot-factorial-comparison.R

library(devtools)
library(tictoc)
load_all('.')
library(gridExtra)

output_dir <- 'analysis/carryover_factorial/output'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv('NREPS', unset = '50'))
n_cores <- as.integer(Sys.getenv('NCORES', unset = '-1'))
if (n_cores < 0) n_cores <- max(1, parallel::detectCores() - 1)

cat(sprintf('Figure 4 comparison: Nreps=%d, cores=%d\n\n', Nreps, n_cores))

# -----------------------------------------------------------------
# Trial designs (identical to 04-generate-power-head.R)
# -----------------------------------------------------------------

tdOL <- buildtrialdesign(
  name_longform = 'open label', name_shortform = 'OL',
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste0('OL', 1:8),
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
  name_longform = 'traditional crossover', name_shortform = 'CO',
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = c(paste0('COa', 1:4), paste0('COb', 1:4)),
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

trialdesigns <- list(tdOL, tdOLBDC, tdCO, tdNof1)

# -----------------------------------------------------------------
# Parameters (matching Figure 4 exactly)
# -----------------------------------------------------------------

data(extracted_bp)
blparamsets <- list(list(name = 'TR', param = extracted_bp))

data(extracted_rp)
respparamsets <- list(list(name = 'TR', param = extracted_rp))

censorparams <- data.table(
  names = c('balanced', 'more of flat',
            'more of biased', 'high dropout'),
  beta0 = c(.05, .05, .05, .15),
  beta1 = c(.5, .2, .8, .5),
  eb1 = c(2, 2, 2, 2)
)

coremodelparams <- data.table(expand.grid(
  N = c(35, 70),
  c.bm = c(0, .25, .45),
  carryover_t1half = c(0, .5, 1.0),
  c.tv = .7, c.pb = .7, c.br = .7,
  c.cf1t = .2, c.cfct = .1
))

simparam <- list(
  Nreps = Nreps,
  progressiveSave = FALSE,
  basesavename = 'fig4cmp',
  nRep2save = 5,
  saveunit2start = 1,
  savedir = output_dir
)

# -----------------------------------------------------------------
# Run A: No carryover in analysis model
# -----------------------------------------------------------------

cat('=== Model A: No carryover in analysis ===\n')
tic('Model A')

ap_ignore <- data.table(
  useDE = FALSE,
  t_random_slope = FALSE,
  full_model_out = FALSE,
  simplecarryover = FALSE,
  carryover_t1half = 0,
  carryover_scalefactor = 1
)

sim_a <- generateSimulatedResults(
  trialdesigns = trialdesigns,
  respparamsets = respparamsets,
  blparamsets = blparamsets,
  censorparams = censorparams,
  modelparams = coremodelparams,
  simparam = simparam,
  analysisparams = ap_ignore,
  rawdataout = FALSE,
  n_cores = n_cores
)
toc()

saveRDS(sim_a, file.path(output_dir,
  sprintf('fig4_no_carryover_n%d.rds', Nreps)))

# -----------------------------------------------------------------
# Run B: Continuous Dbc carryover in analysis model
# (matched to DGP half-life for each parameter set)
#
# Since the DGP half-life varies across modelparams rows,
# we need per-row matching. But analysisparams is a single
# table applied uniformly. Solution: use carryover_t1half = 1.0
# as the analysis model's half-life (the middle of the DGP range).
#
# Actually, the proper approach: the analysis Dbc uses the same
# t1half as the DGP. In the orig code, op$carryover_t1half in
# the analysis reads from analysisparams, not from modelparams.
# So we need to run 3 separate sims (one per DGP t1half) with
# matched analysis t1half.
# -----------------------------------------------------------------

cat('\n=== Model B: Continuous Dbc (matched t1half) ===\n')
tic('Model B')

sim_b_parts <- list()
for (dgp_th in c(0, 0.5, 1.0)) {
  cat(sprintf('  DGP t1half = %.1f\n', dgp_th))

  mp_sub <- coremodelparams[carryover_t1half == dgp_th]

  ap_dbc <- data.table(
    useDE = FALSE,
    t_random_slope = FALSE,
    full_model_out = FALSE,
    simplecarryover = FALSE,
    carryover_t1half = dgp_th,
    carryover_scalefactor = 1
  )

  sim_part <- generateSimulatedResults(
    trialdesigns = trialdesigns,
    respparamsets = respparamsets,
    blparamsets = blparamsets,
    censorparams = censorparams,
    modelparams = mp_sub,
    simparam = simparam,
    analysisparams = ap_dbc,
    rawdataout = FALSE,
    n_cores = n_cores
  )

  sim_b_parts[[as.character(dgp_th)]] <- sim_part$results
}
toc()

sim_b <- list(
  results = rbindlist(sim_b_parts),
  parameterselections = sim_a$parameterselections
)

saveRDS(sim_b, file.path(output_dir,
  sprintf('fig4_with_carryover_n%d.rds', Nreps)))

# -----------------------------------------------------------------
# Plot: Side-by-side Figure 4 heatmaps
# -----------------------------------------------------------------

cat('\n=== Generating heatmaps ===\n')

provenance_a <- sprintf(
  'No carryover in analysis | %d reps | %s',
  Nreps, format(Sys.time(), '%Y-%m-%d'))

provenance_b <- sprintf(
  'Matched Dbc carryover in analysis | %d reps | %s',
  Nreps, format(Sys.time(), '%Y-%m-%d'))

p1a <- PlotModelingResults(
  sim_a, param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'c.bm', 'censorparamset'),
  param2hold = data.table(blparamset = 1, respparamset = 1,
                          carryover_t1half = 0),
  param2nothold = c('c.cfct', 'modelparamset')
) + ggtitle('A1. No carryover in analysis (t1half=0)') +
  labs(subtitle = provenance_a) +
  theme(plot.title = element_text(hjust = 0, size = 9),
        plot.subtitle = element_text(size = 6, color = 'grey40'))

p1b <- PlotModelingResults(
  sim_b, param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'c.bm', 'censorparamset'),
  param2hold = data.table(blparamset = 1, respparamset = 1,
                          carryover_t1half = 0),
  param2nothold = c('c.cfct', 'modelparamset')
) + ggtitle('B1. Matched Dbc carryover (t1half=0)') +
  labs(subtitle = provenance_b) +
  theme(plot.title = element_text(hjust = 0, size = 9),
        plot.subtitle = element_text(size = 6, color = 'grey40'))

p2a <- PlotModelingResults(
  sim_a, param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'carryover_t1half',
                 'censorparamset'),
  param2hold = data.table(blparamset = 1, respparamset = 1,
                          c.bm = 0.45),
  param2nothold = c('c.cfct', 'modelparamset')
) + ggtitle('A2. No carryover in analysis (c.bm=0.45)') +
  labs(subtitle = provenance_a) +
  theme(plot.title = element_text(hjust = 0, size = 9),
        plot.subtitle = element_text(size = 6, color = 'grey40'))

p2b <- PlotModelingResults(
  sim_b, param2plot = 'power',
  param2vary = c('trialdesign', 'N', 'carryover_t1half',
                 'censorparamset'),
  param2hold = data.table(blparamset = 1, respparamset = 1,
                          c.bm = 0.45),
  param2nothold = c('c.cfct', 'modelparamset')
) + ggtitle('B2. Matched Dbc carryover (c.bm=0.45)') +
  labs(subtitle = provenance_b) +
  theme(plot.title = element_text(hjust = 0, size = 9),
        plot.subtitle = element_text(size = 6, color = 'grey40'))

ml <- arrangeGrob(grobs = list(p1a, p1b, p2a, p2b),
                  nrow = 2, ncol = 2)

timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')

pdf_path <- file.path(output_dir,
  paste0('figure4_comparison_', timestamp, '.pdf'))
ggsave(pdf_path, ml, width = 14, height = 12,
       units = 'in', dpi = 300)
cat('Saved:', pdf_path, '\n')

png_path <- file.path(output_dir,
  paste0('figure4_comparison_', timestamp, '.png'))
ggsave(png_path, ml, width = 14, height = 12,
       units = 'in', dpi = 150)
cat('Saved:', png_path, '\n')

cat('\nDone.\n')
