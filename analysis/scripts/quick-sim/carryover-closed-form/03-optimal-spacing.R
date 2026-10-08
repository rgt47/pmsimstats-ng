## Is there an optimal time between visits for a given carryover
## half-life? Closed form of 01-e9-closed-form.R (definitions evaluated,
## nothing simulated).
##
## Designs: 8 equally spaced blinded visits (expectancy 0.5) at spacing d
## weeks, two counterbalanced sequences, three arrangements:
##   alternating  10101010 / 01010101
##   ABAB         11001100 / 00110011
##   CO           11110000 / 00001111
## Rule: decayed mean moderation (same slope as graded coupling, no c_bm
## ceiling, so feasibility does not distort the comparison). Structures:
## CS and separable AR(1) (rho = 0.8 per week). The criterion is the
## per-patient information, the patient-weighted squared contrast-
## biomarker correlation rho^2 without the path term (sequence or period
## differences removed, as an analysis with period terms would), and the
## corresponding power at N = 70, c_bm = 0.25.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/03-optimal-spacing.R
suppressMessages({library(data.table); library(ggplot2)})
src <- readLines('analysis/scripts/quick-sim/carryover-closed-form/01-e9-closed-form.R')
stop_at <- grep('^## ---------------- per-path table', src)[1] - 1
eval(parse(text = src[1:stop_at]))
out_dir <- 'analysis/data/quick-sim/carryover-closed-form'
fig_dir <- 'docs/figures'

arrangements <- list(
  alternating = list(c(1, 0, 1, 0, 1, 0, 1, 0), c(0, 1, 0, 1, 0, 1, 0, 1)),
  ABAB = list(c(1, 1, 0, 0, 1, 1, 0, 0), c(0, 0, 1, 1, 0, 0, 1, 1)),
  CO = list(c(1, 1, 1, 1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, 1, 1, 1)))
d_grid <- seq(0.25, 6, by = 0.25)
halves <- c(0.25, 0.5, 1, 2, 4)

## Separable AR(1) has no measurement-error (white-noise) component, so
## as visits move closer together they become perfectly correlated and
## the contrast's noise vanishes; the optimum is then pushed to the
## shortest spacing for reasons that are not pharmacological. The
## 'nugget' structures add a white-noise share nu to the time kernel,
## A'_ts = (1 - nu) rho^|w_t - w_s| + nu [t = s], keeping the factor part
## K of the separable form (docs/37, Section 5.8).
nug_moments <- function(pi, th, nu) {
  a <- weights(pi, 'E9')
  h <- mod_profile(pi, 'mean_decayed', th)
  s <- cbind(10, 10 * pi$e, sd_BR)
  K <- (1 - dgp$C1) * diag(3) + dgp$C1
  At <- (1 - nu) * dgp$RHO^abs(outer(pi$w, pi$w, '-')) + nu * diag(n)
  aCa <- sum(outer(a, a) * At * (s %*% K %*% t(s)))
  shift <- CBM * sd_BR * sum(a * h)
  data.table(gamma = -CBM * sd_BR / sd_B * sum(a * h), var_D = aCa + shift^2,
             tau2 = aCa, mu = 0)
}
eval_design <- function(arr, d, construct, th) {
  name <- sprintf('%s-d%.2f', arr, d)
  dgp$designs[[name]] <- list(t_wk = rep(d, 8), e = rep(0.5, 8),
                              paths = arrangements[[arr]])
  use_design(name)
  m <- if (startsWith(construct, 'nug')) {
    rbindlist(lapply(paths, nug_moments, th = th,
                     nu = as.numeric(sub('nug', '', construct))))
  } else {
    rbindlist(lapply(paths, path_moments, rule = 'mean_decayed', th = th,
                     variant = 'E9', construct = construct))
  }
  w <- Ns / N
  info <- sum(w * (m$gamma * sd_B)^2 / m$var_D)
  p <- pooled(copy(m)[, mu := 0])$power
  data.table(arrangement = arr, d = d, construct = construct, t_half = th,
             info = info, power = p)
}
grid <- CJ(arrangement = names(arrangements), d = d_grid,
           construct = c('CS', 'sep', 'nug0.2', 'nug0.5'), t_half = halves)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(k)
  eval_design(grid$arrangement[k], grid$d[k], grid$construct[k], grid$t_half[k])))
fwrite(res, file.path(out_dir, 'optimal-spacing.csv'))

## Optimum per arrangement, structure and half-life. An optimum at the
## top of the grid means wider spacing would be better still.
opt <- res[, .SD[which.max(info)], by = .(arrangement, construct, t_half)]
opt[, at_edge := d %in% range(d_grid)]
## The power curve is often a plateau; the shortest spacing reaching 95%
## of the best power is the practical choice (shortest trial).
near <- res[, .(d95 = min(d[power >= 0.95 * max(power)])),
            by = .(arrangement, construct, t_half)]
opt <- merge(opt, near, by = c('arrangement', 'construct', 't_half'))
fwrite(opt, file.path(out_dir, 'optimal-spacing-optimum.csv'))
cat('Optimal spacing (weeks) maximizing per-patient information',
    '(* = at the edge of the 0.25-6 week grid); d95 = shortest spacing within 95% of the best power\n')
print(dcast(opt[, .(arrangement, construct, t_half,
                    v = sprintf('d*=%.2f%s d95=%.2f p=%.3f', d, ifelse(at_edge, '*', ''),
                                d95, power))],
            arrangement + construct ~ t_half, value.var = 'v'))

## Large-n approximation for the alternating design under separable AR(1):
## the contrast's time part is about (1/2)(1 - phi)/(1 + phi) with
## phi = rho^d, and each off visit is d weeks after a discontinuation,
## so information is proportional to
##   f(d) = (1 - 2^(-d/t1/2))^2 (1 + rho^d) / (1 - rho^d).
approx_opt <- rbindlist(lapply(halves, function(th) {
  f <- function(d) (1 - 2^(-d / th))^2 * (1 + dgp$RHO^d) / (1 - dgp$RHO^d)
  o <- optimize(f, c(0.01, 30), maximum = TRUE)
  data.table(t_half = th, d_approx = o$maximum)
}))
cat('\nLarge-n approximation, alternating, separable AR(1)\n')
print(merge(approx_opt, opt[arrangement == 'alternating' & construct == 'sep',
                            .(t_half, d_exact_grid = d)], by = 't_half'))

## Figure: information against spacing.
ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'
## Half-life is an ordered magnitude: one hue, light (short) to dark
## (long).
pal <- setNames(c('#9ec5f0', '#6aa3e3', '#3d82d1', '#2563b0', '#17447e'),
                sprintf('t1/2 = %s wk', format(halves)))
fd <- copy(res)[, th_lab := factor(sprintf('t1/2 = %s wk', format(t_half)), levels = names(pal))]
con_lev <- c(CS = 'Compound symmetry', sep = 'AR(1), no measurement error',
             nug0.2 = 'AR(1) + 20% measurement error', nug0.5 = 'AR(1) + 50% measurement error')
fd[, construct := factor(con_lev[construct], levels = con_lev)]
fd[, arrangement := factor(arrangement, levels = c('alternating', 'ABAB', 'CO'),
                           labels = c('Alternating every visit', 'ABAB', 'CO (one switch)'))]
od <- copy(opt)[, th_lab := factor(sprintf('t1/2 = %s wk', format(t_half)), levels = names(pal))]
od[, construct := factor(con_lev[construct], levels = con_lev)]
od[, arrangement := factor(arrangement, levels = c('alternating', 'ABAB', 'CO'),
                           labels = levels(fd$arrangement))]
g <- ggplot(fd, aes(d, power, colour = th_lab)) +
  geom_line(linewidth = 0.85) +
  geom_point(data = od, size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(construct ~ arrangement) +
  labs(x = 'Time between visits (weeks); 8 visits, so trial length is 8 times this',
       y = 'Power (alpha = 0.05)',
       title = 'Power against visit spacing, by response half-life',
       subtitle = 'Decayed mean moderation, N = 70, c_bm = 0.25, E9 with period terms, closed form',
       caption = 'Points: spacing that maximizes power for that half-life (at the grid edge if wider is better).') +
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
        legend.title = element_blank(), legend.text = element_text(colour = ink))
ggsave(file.path(fig_dir, '39-fig1-optimal-spacing.png'), g, width = 10.5, height = 10,
       dpi = 200, bg = surface)
cat('\nWrote', file.path(fig_dir, '39-fig1-optimal-spacing.png'), '\n')
