## 21-fig4b-g9-heatmap.R
##
## Figure 4B layout of 16-fig4b-graded-heatmap.R with a third version
## row, "Extended DGP, G9 analysis", from 19-g9-planning-and-fig4b.R
## (block F). Reads the stored two-version CSV of
## 15-fig4b-graded-58b32a9.R and the G9 CSV; runs on stored files only.
## The open-label design is not estimable under G9 (Dbc is constant, so
## bm:Dbc is aliased with bm); its tiles are grey and marked n/e. The
## existing figure is left unchanged and the new file is written beside
## it under analysis/figures. Version labels follow STYLE.md: the
## original rule (constant correlation) and the decaying correlation.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/21-fig4b-g9-heatmap.R \
##     [--cbm C] [--reps R] [--base FILE] [--g9 FILE] [--out PNG]
suppressMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
cbm <- as.numeric(argv('--cbm', '0.25'))
reps <- as.integer(argv('--reps', '250'))
base_file <- argv('--base', sprintf(
  'analysis/data/quick-sim/hendrickson-problems/fig4b-graded/fig4b-graded-cbm%.2f-reps%d.csv',
  cbm, reps))
g9_file <- argv('--g9', sprintf(
  'analysis/data/quick-sim/carryover-closed-form/g9-planning/fig4b-g9-cbm%.2f.csv', cbm))
out_file <- argv('--out', sprintf(
  'analysis/figures/quick-sim/hendrickson-problems/fig4b-graded-g9-cbm%.2f.png', cbm))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

cols <- c('censorparamset', 'power', 'version', 'design', 'N', 't_half', 'censoring')
d <- rbind(fread(base_file)[, ..cols], fread(g9_file)[, ..cols])
stopifnot(setequal(unique(d$version), c('published', 'graded', 'extended-G9')))
v_levels <- c('Original rule (58b32a9)\noriginal analysis',
              'Decaying correlation\noriginal analysis',
              'Extended DGP\nG9 analysis')
d[, design := factor(design, levels = c('OL', 'OL+BDC', 'CO', 'N-of-1'))]
d[, censoring := factor(censoring, levels = c('high dropout', 'more of biased',
                                              'more of flat', 'balanced', 'No censoring'))]
d[, version := factor(fcase(version == 'published', v_levels[1],
                            version == 'graded', v_levels[2],
                            default = v_levels[3]), levels = v_levels)]
d[, N := factor(paste0('N = ', N), levels = c('N = 35', 'N = 70'))]
d[, label := fifelse(is.na(power), 'n/e', sprintf('%.2f', power))]
p <- ggplot(d, aes(factor(t_half), censoring, fill = power)) +
  geom_tile(color = 'white') +
  geom_text(aes(label = label), size = 2.6) +
  facet_grid(version + N ~ design) +
  scale_fill_distiller(palette = 'RdYlBu', limits = c(0, 1), name = 'Power',
                       na.value = 'grey85') +
  labs(x = 'Carryover half-life (weeks)', y = 'Censoring parameters',
       caption = 'n/e: not estimable. Under G9 the open-label design has no off-drug contrast.') +
  theme_grey(base_size = 10) +
  theme(strip.text.y = element_text(size = 8))
ggsave(out_file, p, width = 9, height = 12, units = 'in', dpi = 200)
cat('Wrote', out_file, '\n')
