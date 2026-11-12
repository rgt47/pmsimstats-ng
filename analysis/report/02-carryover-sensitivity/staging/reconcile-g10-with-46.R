## analysis/report/02-carryover-sensitivity/staging/reconcile-g10-with-46.R
##
## Reconcile paper 01's G10 sample-size runs (19-g9-planning-and-fig4b.R,
## block N, --params published) with paper 02's A4 runs
## (46-run-sample-size.R, --params published). A4 is G10 (stored code
## E7null_cr2).
##
## Same DGP (inspected 2026-11-11): both drivers call
## generate_data_multi_path() of simulation-core.R with build_design_set(),
## c.tv = c.pb = c.br = rho = 0.7, c.cf1t = c.cfct = 0.2, exponential
## decay, mvn, no prognostic terms, no censoring, stop on non-PD, and fit
## E7null_cr2 with fit_spec_s7(). They differ only in seeds and grid:
##   01: Hybrid, c_bm 0.25/0.35/0.45, t1/2 0/0.2/1, six N per cell;
##   46: CO/OLBDC/Hybrid, c_bm 0.2/0.3/0.4, t1/2 0/0.5/1, five N per cell.
## The grids share no cell, so the comparison is through a curve.
##
## Model (an approximation, stated): within a (design, t1/2), power is
## a probit curve in x = c_bm * sqrt(N), since the noncentrality of the
## interaction test is close to proportional to c_bm sqrt(N) for small
## c_bm. Two checks on the Hybrid design at t1/2 = 0 and 1:
##   (1) likelihood-ratio test of a study shift in intercept and slope;
##   (2) each study's own fit inverted at power 0.8 for c_bm 0.25 to
##       0.40, side by side.
## Agreement within Monte Carlo error is expected; a difference would
## point to a code or seed-handling divergence, not to the DGP.
##
## Usage (from the repository root, after the 46 published Hybrid job
## has been fetched):
##   Rscript analysis/report/02-carryover-sensitivity/staging/reconcile-g10-with-46.R

f01 <- 'analysis/data/quick-sim/carryover-closed-form/published-planning/block-N/n80-grid-g10.csv'
f46 <- 'analysis/data/carryover-sensitivity/sample-size-published-Hybrid/summary.rds'
if (!file.exists(f46)) stop('46 output not found: ', f46)
stopifnot(identical(readLines(file.path(dirname(f01), 'PARAM_SET'), n = 1L),
                    'published'))
m46 <- readRDS(file.path(dirname(f46), 'meta.rds'))
stopifnot(identical(m46$params, 'published'))

d01 <- read.csv(f01)
d01 <- data.frame(study = '01', t_half = d01$t_half, c_bm = d01$c_bm,
                  N = d01$N, k = d01$n_reject, n = d01$n_p)
s46 <- readRDS(f46)
s46 <- s46[s46$method == 'E7null_cr2' & s46$design == 'Hybrid', ]
d46 <- data.frame(study = '46', t_half = s46$t1half, c_bm = s46$c_bm,
                  N = s46$N, k = round(s46$reject * s46$n), n = s46$n)
d <- rbind(d01, d46)
d <- d[d$t_half %in% c(0, 1), ]
d$x <- d$c_bm * sqrt(d$N)
stopifnot(all(table(d$study, d$t_half) > 0))

n80 <- function(fit, cb) {
  b <- coef(fit)
  ((qnorm(0.8) - b[[1]]) / b[[2]] / cb)^2
}
for (th in c(0, 1)) {
  dd <- d[d$t_half == th, ]
  f0 <- glm(cbind(k, n - k) ~ x, binomial('probit'), data = dd)
  f1 <- glm(cbind(k, n - k) ~ x * study, binomial('probit'), data = dd)
  lr <- anova(f0, f1, test = 'Chisq')
  cat(sprintf('\nHybrid, t1/2 = %g week: LR test of a study shift, chi2 = %.2f on %d df, p = %.3f\n',
              th, lr$Deviance[2], lr$Df[2], lr$`Pr(>Chi)`[2]))
  cat(sprintf('  residual deviance of the common curve: %.1f on %d df\n',
              deviance(f0), df.residual(f0)))
  g01 <- glm(cbind(k, n - k) ~ x, binomial('probit'), data = dd[dd$study == '01', ])
  g46 <- glm(cbind(k, n - k) ~ x, binomial('probit'), data = dd[dd$study == '46', ])
  cb <- c(0.25, 0.30, 0.35, 0.40)
  print(data.frame(c_bm = cb, N80_paper01 = round(n80(g01, cb)),
                   N80_paper02 = round(n80(g46, cb)),
                   N80_pooled = round(n80(f0, cb))), row.names = FALSE)
}
