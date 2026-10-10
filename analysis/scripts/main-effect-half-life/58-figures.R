## docs/58, figures from the stored outputs of 58-plim-and-power.R and
## 58-randomization-example.R. Writes docs/figures/58-fig1-closed-form.png
## and docs/figures/58-fig2-example-trial.png.
## Run from the repository root:
##   Rscript analysis/scripts/main-effect-half-life/58-figures.R
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
OUT <- 'analysis/data/quick-sim/main-effect-half-life'
PAL <- c('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4')
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = 'bottom',
        plot.title = element_text(size = 10, face = 'bold'))

alt <- fread(file.path(OUT, '58-cf-alt.csv'))[design == 'Hybrid' & type %in% c('F1', 'VF')]
alt[, `:=`(ratio = plim_exp / exp, th = factor(t_half))]
pA <- ggplot(alt, aes(h, ratio, colour = th)) + geom_hline(yintercept = 1, colour = 'grey60') +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) + facet_wrap(~ type) +
  scale_colour_manual(values = PAL, name = expression(true~t[1/2]~(weeks))) +
  labs(x = 'analysis half-life h (weeks)', y = expression(plim~hat(beta)[D]^{exp}*(h)/beta[D]^{exp}),
       title = 'A. Model-implied exposure contrast relative to the estimand') + th
pB <- ggplot(alt, aes(h, size_or_power, colour = th)) + geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) + facet_wrap(~ type) +
  scale_colour_manual(values = PAL, name = expression(true~t[1/2]~(weeks))) +
  labs(x = 'analysis half-life h (weeks)', y = 'power (normal approximation)',
       title = expression(bold('B. Power at ')*bolditalic(f[m])*bold(' = 0.4, ')*bolditalic(N)*bold(' = 70'))) + th
nul <- fread(file.path(OUT, '58-cf-null.csv'))[cov == 'rc' & is.finite(h) & design != 'CO']
pC <- ggplot(nul, aes(h, size_or_power, colour = type)) + geom_hline(yintercept = 0.05, colour = 'grey60') +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) + facet_wrap(~ design, scales = 'free_y') +
  scale_colour_manual(values = PAL[c(1, 2, 3)], name = 'analysis') +
  labs(x = 'analysis half-life h (weeks)', y = 'size under the sharp null',
       title = expression(bold('C. Size from the null bias ')*bolditalic(b)[0]*bold('(h), ')*bolditalic(N)*bold(' = 70'))) + th
rp <- fread(file.path(OUT, '58-cf-randomization-power.csv'))[type %in% c('F1', 'VF')]
rpl <- melt(rp, id.vars = c('design', 'type', 'f_m'), measure.vars = c('power_model', 'power_rand'))
rpl[, test := fifelse(variable == 'power_model', 'model-based', 'randomization')]
pD <- ggplot(rpl, aes(f_m, value, colour = type, linetype = test)) + geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) + facet_wrap(~ design) +
  scale_colour_manual(values = PAL[c(2, 3)], name = 'analysis') +
  scale_linetype_manual(values = c('solid', 'dashed'), name = 'test') +
  labs(x = expression(drug~effect~multiplier~f[m]), y = 'power',
       title = 'D. Model-based and randomization tests, h = 0.2') + th
ggsave('docs/figures/58-fig1-closed-form.png', (pA / pB / pC / pD), width = 7.5, height = 12, dpi = 200)

## Example trial
prof <- fread(file.path(OUT, '58-ex-profile.csv'))[is.finite(h) & h <= 6]
cm <- fread(file.path(OUT, '58-ex-cmax.csv'))
pl <- melt(prof, id.vars = 'h', measure.vars = c('tF', 'tV'))
pl[, type := fifelse(variable == 'tF', 'F1', 'VF')]
e1 <- ggplot(pl, aes(h, value, colour = type)) +
  geom_hline(yintercept = c(-1, 1) * qnorm(0.975), colour = 'grey60') +
  geom_hline(yintercept = -cm[type == 'F1', c_max], colour = PAL[2], linetype = 'dashed') +
  geom_hline(yintercept = -cm[type == 'VF', c_max], colour = PAL[3], linetype = 'dashed') +
  geom_line(linewidth = 0.7) + scale_colour_manual(values = PAL[c(2, 3)], name = 'analysis') +
  labs(x = 'analysis half-life h (weeks)', y = 'CR2 t statistic of the drug coefficient',
       title = 'A. T(h) across the half-life') + th
e2 <- ggplot(prof, aes(h, devF)) + geom_hline(yintercept = c(2.71, 3.84), colour = 'grey60') +
  geom_line(linewidth = 0.7, colour = PAL[2]) + geom_point(size = 1.2, colour = PAL[2]) +
  labs(x = 'analysis half-life h (weeks)', y = 'ML deviance of F1', title = 'B. Profile of h') + th
sims <- fread(file.path(OUT, '58-ex-randomization-fast-sims.csv'))
obs <- fread(file.path(OUT, '58-ex-randomization-fast.csv'))
sl <- rbind(data.table(type = 'F1', v = sims$F1_max_t), data.table(type = 'VF', v = sims$VF_max_t))
ol <- data.table(type = c('F1', 'VF'), v = obs[match(c('F1_max_t', 'VF_max_t'), statistic), observed])
e3 <- ggplot(sl, aes(v)) + geom_histogram(bins = 40, fill = 'grey70', colour = 'white') +
  geom_vline(data = ol, aes(xintercept = v), colour = PAL[1], linewidth = 0.8) + facet_wrap(~ type) +
  labs(x = 'max |t| over the grid, relabeled trials', y = 'count',
       title = 'C. Randomization distribution (999 permutations)') + th
ggsave('docs/figures/58-fig2-example-trial.png', (e1 | e2) / e3, width = 7.5, height = 6.5, dpi = 200)
cat('figures written\n')
