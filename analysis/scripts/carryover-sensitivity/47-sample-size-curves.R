## analysis/scripts/carryover-sensitivity/47-sample-size-curves.R
##
## Total N for 80 percent power, from the power curves of
## 46-run-sample-size.R. Reads stored output only; no simulation.
##
## Model: within a (parameter set, design, true t1/2, analysis), power
## is a probit curve in x = c_bm * sqrt(N),
##   Phi^{-1}(power) = a + b x,
## fitted by binomial maximum likelihood over the 15 cells (three c_bm,
## five N). The noncentrality of the interaction test is close to
## proportional to c_bm sqrt(N), so one curve serves all c_bm; the
## residual deviance of each fit is reported as the check, and paper
## 01's runs found no gain from separate curves per c_bm (p = 0.43 to
## 0.74). N80 = ((Phi^{-1}(0.8) - a) / b / c_bm)^2. Its 95 percent
## interval comes from 4,000 draws of (a, b) from their estimated
## normal distribution (parametric bootstrap of the fit, Monte Carlo
## error only).
##
## Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/47-sample-size-curves.R
## Writes analysis/data/carryover-sensitivity/sample-size-n80.csv,
## sample-size-fit.csv and analysis/figures/02-sample-size-curves-<set>.{pdf,png}.
suppressMessages({ library(dplyr); library(purrr); library(tidyr)
  library(ggplot2) })
source('analysis/scripts/carryover-sensitivity/spec-labels.R')

dat_dir <- 'analysis/data/carryover-sensitivity'
dirs <- Sys.glob(file.path(dat_dir, 'sample-size-*-*'))
dirs <- dirs[!grepl('smoke', dirs)]
s <- map_dfr(dirs, function(d) {
  m <- readRDS(file.path(d, 'meta.rds'))
  readRDS(file.path(d, 'summary.rds')) |> mutate(params = m$params)
})
stopifnot(nrow(s) == 6 * 45 * 2, all(s$n >= 450))

s <- s |> mutate(k = round(reject * n), x = c_bm * sqrt(N),
                 analysis = to_a(method))
C_TAB <- c(0.2, 0.25, 0.3, 0.35, 0.4)
z80 <- qnorm(0.8)
set.seed(20261112)

fits <- s |>
  group_by(params, design, t1half, analysis) |>
  group_modify(function(d, key) {
    f <- glm(cbind(k, n - k) ~ x, family = binomial('probit'), data = d)
    b <- MASS::mvrnorm(4000, coef(f), vcov(f))
    x80_draw <- (z80 - b[, 1]) / b[, 2]
    tibble(a = coef(f)[1], b = coef(f)[2],
           x80 = (z80 - coef(f)[1]) / coef(f)[2],
           x80_lo = unname(quantile(x80_draw, 0.025)),
           x80_hi = unname(quantile(x80_draw, 0.975)),
           deviance = deviance(f), df = df.residual(f),
           gof_p = pchisq(deviance(f), df.residual(f), lower.tail = FALSE))
  }) |> ungroup()

n80 <- fits |>
  crossing(c_bm = C_TAB) |>
  mutate(N80 = (x80 / c_bm)^2, N80_lo = (x80_lo / c_bm)^2,
         N80_hi = (x80_hi / c_bm)^2) |>
  select(params, design, t1half, analysis, c_bm, N80, N80_lo, N80_hi)
write.csv(n80, file.path(dat_dir, 'sample-size-n80.csv'), row.names = FALSE)
write.csv(fits, file.path(dat_dir, 'sample-size-fit.csv'), row.names = FALSE)

## Check of the common-curve model: a separate probit on sqrt(N) for
## each c_bm, inverted at 0.8, against the common-curve N80.
sep <- s |>
  group_by(params, design, t1half, analysis, c_bm) |>
  group_modify(function(d, key) {
    f <- glm(cbind(k, n - k) ~ sqrt(N), family = binomial('probit'), data = d)
    tibble(N80_sep = ((z80 - coef(f)[1]) / coef(f)[2])^2)
  }) |> ungroup() |>
  inner_join(n80, by = c('params', 'design', 't1half', 'analysis', 'c_bm')) |>
  mutate(ratio = N80_sep / N80,
         inside = N80_sep >= N80_lo & N80_sep <= N80_hi)
write.csv(sep, file.path(dat_dir, 'sample-size-sepcheck.csv'), row.names = FALSE)
cat(sprintf('Separate per-c_bm N80 / common N80: %.3f to %.3f; %d of %d inside the common interval\n',
            min(sep$ratio), max(sep$ratio), sum(sep$inside), nrow(sep)))

cat('Goodness of fit of the probit curves (deviance on df, p):\n')
print(as.data.frame(fits |> mutate(across(c(deviance, gof_p), ~ round(.x, 3))) |>
  select(params, design, t1half, analysis, deviance, df, gof_p)), row.names = FALSE)

## Figures: observed power with the fitted curves, one per parameter set.
design_lab <- c(CO = 'CO', OLBDC = 'OL+BDC', Hybrid = 'Hybrid')
for (p in unique(s$params)) {
  sp <- s |> filter(params == p) |>
    mutate(mcse = sqrt(reject * (1 - reject) / n))
  grid <- fits |> filter(params == p) |>
    crossing(c_bm = c(0.2, 0.3, 0.4)) |>
    group_by(params, design, t1half, analysis, c_bm) |>
    reframe(N = seq(min(sp$N), max(sp$N), length.out = 200)) |>
    left_join(fits |> filter(params == p),
              by = c('params', 'design', 't1half', 'analysis')) |>
    mutate(power = pnorm(a + b * c_bm * sqrt(N)))
  lab <- function(d) d |>
    mutate(design = factor(design_lab[design], levels = design_lab),
           half = factor(sprintf('t1/2 = %g wk', t1half)),
           c_bm = factor(c_bm))
  g <- ggplot(lab(sp), aes(N, reject, colour = c_bm, shape = analysis)) +
    geom_hline(yintercept = 0.8, linewidth = 0.3, colour = 'grey40') +
    geom_line(data = lab(grid), aes(N, power, colour = c_bm,
                                    linetype = analysis),
              inherit.aes = FALSE, linewidth = 0.45) +
    geom_errorbar(aes(ymin = reject - 2 * mcse, ymax = reject + 2 * mcse),
                  width = 0, linewidth = 0.3) +
    geom_point(size = 1.3) +
    facet_grid(half ~ design) +
    scale_x_log10(breaks = c(30, 50, 100, 200, 500, 1000)) +
    scale_colour_manual(values = c('0.2' = '#D55E00', '0.3' = '#009E73',
                                   '0.4' = '#0072B2'),
                        name = expression(c[bm])) +
    scale_shape_manual(values = c(A2 = 1, A4 = 16), name = NULL) +
    scale_linetype_manual(values = c(A2 = 'dashed', A4 = 'solid'),
                          name = NULL) +
    labs(x = 'Total sample size N (log scale)', y = 'Power') +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = 'bottom',
          strip.background = element_rect(fill = 'grey92', colour = NA))
  for (ext in c('pdf', 'png')) {
    ggsave(file.path('analysis/figures',
                     sprintf('02-sample-size-curves-%s.%s', p, ext)),
           g, width = 6.5, height = 5.2, dpi = 300)
  }
}

cat('\nN for 80 percent power (A4), by parameter set:\n')
print(as.data.frame(n80 |> filter(analysis == 'A4', c_bm %in% c(0.2, 0.3, 0.4)) |>
  mutate(N80 = round(N80)) |>
  select(params, design, t1half, c_bm, N80) |>
  pivot_wider(names_from = c_bm, values_from = N80, names_prefix = 'cbm')),
  row.names = FALSE)
