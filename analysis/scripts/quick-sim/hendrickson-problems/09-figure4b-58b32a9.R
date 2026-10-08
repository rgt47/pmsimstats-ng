## Reproduce Figure 4, panel B of Hendrickson et al. (2020) from the
## commit 58b32a9 ('orig'): its saved results_core.rda, plotted with its
## own PlotModelingResults() and the exact call in its vignette
## Produce_Publication_Results_2_Look_at_data.Rmd. Nothing is
## re-simulated.
##
## Requires 00-fetch-58b32a9.sh. If plottingfunctions.R is not in the
## cache, point H58_PLOT at a copy fetched from 58b32a9. Run from the
## repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/09-figure4b-58b32a9.R
suppressMessages({library(data.table); library(ggplot2)})
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
plot_src <- Sys.getenv('H58_PLOT', file.path(src, 'R', 'plottingfunctions.R'))
stopifnot(file.exists(plot_src))
out_dir <- 'analysis/figures/quick-sim/hendrickson-problems'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

env <- new.env()
sys.source(plot_src, envir = env)
load(file.path(src, 'data', 'results_core.rda'))
simresults <- results_core

## Verbatim from the 58b32a9 vignette.
param2vary <- c('trialdesign', 'N', 'carryover_t1half', 'censorparamset')
param2hold <- data.table(blparamset = 1, respparamset = 1, c.bm = .6)
param2nothold <- c('c.cfct', 'modelparamset')
p2 <- env$PlotModelingResults(simresults, 'power', param2vary, param2hold,
                              param2nothold)
p2 <- p2 + ggtitle('B. Effect of nonzero carryover term on power') +
  theme(plot.title = element_text(hjust = 0))

ggsave(file.path(out_dir, 'fig4b-58b32a9.png'), p2, width = 9, height = 4.5,
       units = 'in', dpi = 200)
ggsave(file.path(out_dir, 'fig4b-58b32a9.pdf'), p2, width = 9, height = 4.5,
       units = 'in')

## The plotted cell values, for checking against the published figure.
d <- as.data.table(p2$data)
fwrite(d, file.path(out_dir, 'fig4b-58b32a9-values.csv'))
print(head(d, 3))
cat(sprintf('reps per cell: %d\n',
            nrow(simresults$results) /
              nrow(unique(simresults$results[, .(trialdesign, N, c.bm,
                carryover_t1half, censorparamset)]))))
cat('Wrote', file.path(out_dir, 'fig4b-58b32a9.png'), '\n')
