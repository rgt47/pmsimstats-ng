## Summarize the common-grid contrast (12-common-grid-contrast.R): the
## strawman against the candidate fixes, N = 70.
##
##   Problem 1  share of correlation matrices that the published code
##              would repair (strawman: repaired during the run; fixes:
##              fail corpcor::is.positive.definite(), so not run)
##   Problem 2  power across carryover half-life, no censoring
##   Size       Type I error (c.bm = 0)
## Reads whichever analyses have finished. Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/13-summarize-contrast.R [--reps R]
suppressMessages({library(data.table); library(ggplot2); library(corpcor)})
args <- commandArgs(trailingOnly = TRUE)
reps <- if ('--reps' %in% args) as.integer(args[which(args == '--reps') + 1]) else 100L
here <- 'analysis/scripts/quick-sim/hendrickson-problems'
root <- 'analysis/data/quick-sim/hendrickson-problems/contrast'
fig_dir <- 'docs/figures'

lab <- c(S = 'S strawman (58b32a9)', B = 'B AR(1)', F = 'F AR(1) + graded, covar',
         C = 'C AR(1) + graded, separable', D = 'D mean moderation',
         E = 'E all three')
order_cf <- c('S', 'B', 'F', 'C', 'D', 'E')
cens_lab <- c('none', 'balanced', 'more of flat', 'more of biased', 'high dropout')

read_analysis <- function(a) {
  d <- file.path(root, sprintf('sets-reps%d-%s', reps, a))
  fs <- list.files(d, pattern = '^[SBFCDE]-[0-9]+\\.rds$', full.names = TRUE)
  if (!length(fs)) return(NULL)
  x <- lapply(fs, readRDS)
  list(status = rbindlist(lapply(x, `[[`, 'status'))[, analysis := a],
       results = rbindlist(lapply(x, `[[`, 'results'), fill = TRUE)[, analysis := a])
}
got <- Filter(Negate(is.null), lapply(c('published', 'car1'), read_analysis))
if (!length(got)) stop('no contrast output for reps = ', reps)
status <- rbindlist(lapply(got, `[[`, 'status'), fill = TRUE)
res <- rbindlist(lapply(got, `[[`, 'results'), fill = TRUE)
cat('Sets found:\n')
print(dcast(status[, .N, by = .(analysis, config)], config ~ analysis, value.var = 'N'))

## ---------------- Problem 1 ----------------
## One matrix per design path x c.bm x t1/2. Strawman: the run's repair
## log (repairs per set = paths repaired). Fixes: matrices built with the
## 04 generator for every path, including sets not run.
keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
dgp <- new.env()
for (e in parse(file.path(here, '04-power-simulation.R'))) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
designs <- dgp$designs
designs$OL <- list(t_wk = rep(2.5, 8), e = rep(1, 8), paths = list(rep(1, 8)))
dkey <- c(OL = 'OL', `OL+BDC` = 'OL+BDC', CO = 'CO', `N-of-1` = 'Hybrid')
grid <- unique(status[analysis == status$analysis[1], .(design, c.bm, carryover_t1half)])
p1 <- list()
for (cf in setdiff(order_cf, 'S')) for (i in seq_len(nrow(grid))) {
  g <- grid[i]; d <- designs[[dkey[[g$design]]]]
  for (pth in d$paths) {
    b <- dgp$build(dgp$path_info(d$t_wk, d$e, pth), cf, g$c.bm, g$carryover_t1half)
    p1[[length(p1) + 1]] <- data.table(config = cf, design = g$design, c.bm = g$c.bm,
                                       repaired = !is.positive.definite(b$Sigma))
  }
}
npaths <- c(OL = 1L, `OL+BDC` = 2L, CO = 2L, `N-of-1` = 4L)
sm <- unique(status[config == 'S' & analysis == status$analysis[1],
                    .(design, c.bm, carryover_t1half, repairs)])
sm <- sm[, .(config = 'S', design, c.bm, built = npaths[design], repairs)]
p1 <- rbind(rbindlist(p1)[, .(built = .N, repairs = sum(repaired)), by = .(config, design, c.bm)],
            sm[, .(built = sum(built), repairs = sum(repairs)), by = .(config, design, c.bm)])
prob1 <- p1[c.bm > 0, .(matrices = sum(built), repaired = sum(repairs)), by = config]
prob1[, pct := round(100 * repaired / matrices, 1)]
prob1[, config := factor(lab[config], levels = lab[order_cf])]
cat('\nProblem 1: matrices the published code would repair (c.bm 0.25 and 0.3)\n')
print(prob1[order(config)])
cat('\n... by c.bm\n')
p1c <- p1[c.bm > 0, .(cell = sprintf('%d/%d', sum(repairs), sum(built))), by = .(config, c.bm)]
print(dcast(p1c, factor(lab[config], levels = lab[order_cf]) ~ c.bm, value.var = 'cell'))

## ---------------- power and size ----------------
pw <- res[, .(power = mean(p < .05, na.rm = TRUE), n = sum(!is.na(p)),
              mean_beta = mean(beta, na.rm = TRUE)),
          by = .(analysis, config, design = c('OL', 'OL+BDC', 'CO', 'N-of-1')[trialdesign],
                 c.bm, t_half = carryover_t1half, censor = censorparamset)]
pw[, censor := cens_lab[censor + 1]]
fwrite(pw, file.path(root, sprintf('contrast-power-reps%d.csv', reps)))
fwrite(prob1, file.path(root, sprintf('contrast-problem1-reps%d.csv', reps)))
show <- function(d, title) {
  cat('\n', title, '\n', sep = '')
  d[, config := factor(lab[config], levels = lab[order_cf])]
  print(dcast(d, analysis + config + design ~ t_half, value.var = 'power'))
}
show(pw[censor == 'none' & c.bm == 0], 'Type I error, no censoring (columns: t1/2)')
show(pw[censor == 'none' & c.bm == 0.25], 'Power, c.bm = 0.25, no censoring')
show(pw[censor == 'none' & c.bm == 0.3], 'Power, c.bm = 0.3, no censoring (S repaired where noted above)')

## ---------------- figure ----------------
pal <- setNames(c('#2a78d6', '#eb6834', '#008300', '#1baf7a', '#eda100', '#e87ba4'),
                lab[order_cf])
shp <- setNames(c(16, 17, 1, 15, 18, 25), lab[order_cf])
lty <- setNames(c('solid', 'solid', 'dashed', 'solid', 'solid', 'solid'), lab[order_cf])
ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'
fd <- pw[censor == 'none' & c.bm == 0.25]
fd[, config := factor(lab[config], levels = lab[order_cf])]
fd[, model := factor(ifelse(analysis == 'published', 'Published analysis', 'corCAR1 analysis'),
                     levels = c('Published analysis', 'corCAR1 analysis'))]
fd[, design := factor(design, levels = c('N-of-1', 'OL+BDC', 'CO', 'OL'))]
fd[, t_lab := factor(format(t_half, drop0trailing = TRUE), levels = c('0', '0.1', '0.2', '0.5', '1'))]
g <- ggplot(fd, aes(t_lab, power, colour = config, shape = config, linetype = config,
                    group = config)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.4, fill = 'white', stroke = 0.9) +
  scale_colour_manual(values = pal) + scale_shape_manual(values = shp) +
  scale_linetype_manual(values = lty) + scale_y_continuous(limits = c(0, 1)) +
  facet_grid(model ~ design) +
  labs(x = 'Carryover half-life (weeks)', y = 'Power (nominal, alpha = 0.05)',
       title = 'Strawman against the fixes, all run through the 58b32a9 pipeline',
       subtitle = sprintf('N = 70, c.bm = 0.25 (feasible in every configuration), no censoring, %d replicates per cell', reps),
       caption = paste('Steps on the x-axis are evenly spaced for legibility; the half-lives are not.',
                       sprintf('Monte Carlo SE about %.2f at power 0.5.', sqrt(.25 / reps)), sep = '\n')) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = surface, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = '#e6e5e1', linewidth = 0.3),
        axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
        strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
        plot.title = element_text(colour = ink, face = 'bold'),
        plot.subtitle = element_text(colour = ink2),
        plot.caption = element_text(colour = ink2, hjust = 0),
        plot.title.position = 'plot', plot.caption.position = 'plot',
        legend.position = 'top', legend.justification = 'left',
        legend.title = element_blank(), legend.text = element_text(colour = ink)) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2))
f <- file.path(fig_dir, sprintf('36-fig4-contrast-reps%d.png', reps))
ggsave(f, g, width = 10, height = 6.4, dpi = 200, bg = surface)
cat('\nWrote', f, '\n')
