## docs/57, Figure 1, drawn from the outputs of 57-example-trial.R
## (which sources this file at its end). Can be run on its own from the
## repository root.
suppressMessages({ library(data.table); library(ggplot2) })
out <- 'analysis/data/quick-sim/data-chosen-half-life'
prof <- fread(file.path(out, '57-example-profile.csv'))
pl <- fread(file.path(out, '57-example-profile-interval.csv'))
lmax0 <- pl$ll_max[1]; lmax1 <- pl$ll_max[2]
grid5 <- c(0, 0.25, 0.5, 1, 2)
pp <- prof[h <= 20]
dl <- rbind(pp[, .(h, value = 2 * (lmax1 - ll1), what = 'with the interaction (F1)')],
            pp[, .(h, value = 2 * (lmax0 - ll0), what = 'without the interaction (null model)')])
brk <- c(0, 0.25, 0.5, 1, 2, 5, 10, 20)
xs <- scale_x_continuous(trans = scales::transform_log1p(), breaks = brk,
                         labels = format(brk, drop0trailing = TRUE, trim = TRUE))
p1 <- ggplot(dl, aes(h, value, linetype = what)) +
  geom_hline(yintercept = qchisq(0.95, 1), colour = 'grey50', linewidth = 0.3) +
  geom_line() +
  geom_point(data = dl[h %in% grid5], size = 1.8) +
  xs +
  labs(x = 'analysis half-life (weeks, log(1 + h) scale)',
       y = 'ML deviance from the maximum', linetype = NULL,
       title = 'A. Profile of the analysis half-life') +
  theme_bw(base_size = 10) + theme(legend.position = 'bottom')
p2 <- ggplot(pp, aes(h, t)) +
  geom_hline(yintercept = c(-1.96, 0, 1.96), colour = 'grey50', linewidth = 0.3) +
  geom_line() + geom_point(data = pp[h %in% grid5], size = 1.8) +
  xs + coord_cartesian(ylim = c(-2.2, 2.2)) +
  labs(x = 'analysis half-life (weeks, log(1 + h) scale)',
       y = 'interaction t statistic (CR2)', title = 'B. Interaction test across the half-life') +
  theme_bw(base_size = 10)
png('docs/figures/57-fig1-example-trial.png', width = 9, height = 4, units = 'in', res = 200)
gridExtra::grid.arrange(p1, p2, ncol = 2)
invisible(dev.off())
