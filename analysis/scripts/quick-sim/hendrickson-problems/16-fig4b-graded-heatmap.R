## Heatmap in the layout of Hendrickson et al. (2020) Figure 4B, for the
## output of 15-fig4b-graded-58b32a9.R: power by censoring set (rows) and
## carryover half-life (columns), faceted by trial design and N, for the
## published coupling rule and for graded coupling.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/16-fig4b-graded-heatmap.R \
##     [--reps R] [--cbm C]
suppressMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.numeric(args[which(args == k) + 1]) else d
reps <- as.integer(arg('--reps', 250))
cbm <- arg('--cbm', 0.25)
in_file <- sprintf('analysis/data/quick-sim/hendrickson-problems/fig4b-graded/fig4b-graded-cbm%.2f-reps%d.csv',
                   cbm, reps)
out_dir <- 'analysis/figures/quick-sim/hendrickson-problems'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
d <- fread(in_file)
d[, design := factor(design, levels = c('OL', 'OL+BDC', 'CO', 'N-of-1'))]
d[, censoring := factor(censoring, levels = c('high dropout', 'more of biased',
                                              'more of flat', 'balanced', 'No censoring'))]
d[, version := factor(fifelse(version == 'published', 'Published coupling (58b32a9)',
                              'Graded coupling'),
                      levels = c('Published coupling (58b32a9)', 'Graded coupling'))]
d[, N := factor(paste0('N = ', N), levels = c('N = 35', 'N = 70'))]
p <- ggplot(d, aes(factor(t_half), censoring, fill = power)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = sprintf('%.2f', power)), size = 2.6) +
  facet_grid(version + N ~ design) +
  scale_fill_distiller(palette = 'RdYlBu', limits = c(0, 1), name = 'Power') +
  labs(x = 'Carryover half-life (weeks)', y = 'Censoring parameters') +
  theme_grey(base_size = 10) +
  theme(strip.text.y = element_text(size = 8))
base <- file.path(out_dir, sprintf('fig4b-graded-cbm%.2f', cbm))
ggsave(paste0(base, '.png'), p, width = 9, height = 8.5, units = 'in', dpi = 200)
ggsave(paste0(base, '.pdf'), p, width = 9, height = 8.5, units = 'in')
cat('Wrote', paste0(base, '.png'), '\n')
