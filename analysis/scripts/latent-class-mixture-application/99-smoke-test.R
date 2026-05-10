## Smoke test: confirm lcmm_analysis() runs to completion on a
## single simulated trial under the existing covariance-only DGP,
## and reports diagnostics matching the documented contract.
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/99-smoke-test.R

library(devtools)
load_all('.')
library(data.table)

source('analysis/scripts/latent-class-mixture-application/01-lcmm-wrapper.R')
source('analysis/scripts/latent-class-mixture-application/02-heterogeneous-slopes-dgp.R')

set.seed(20260506)

td <- buildtrialdesign(
  name_longform = 'open label',
  name_shortform = 'OL',
  timepoints = cumulative(rep(2.5, 8)),
  timeptnames = paste0('OL', 1:8),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)

mp <- data.table(
  N = 70, c.bm = 0.45,
  carryover_t1half = 1.0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.1, c.cfct = 0.05
)
rp <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)
bp <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)

dat <- generateData(mp, rp, bp, td$trialpaths[[1]],
                    empirical = FALSE, makePositiveDefinite = TRUE)
dat[, path := 1]

op <- list(
  useDE = TRUE,
  t_random_slope = FALSE,
  full_model_out = FALSE,
  carryover_t1half = 1.0,
  simplecarryover = FALSE,
  carryover_scalefactor = 1
)

cat('=== lme_analysis baseline ===\n')
t_lme <- system.time(out_lme <- lme_analysis(td$trialpaths, dat, op))
print(out_lme)
cat('Elapsed:', t_lme['elapsed'], 's\n\n')

cat('=== lcmm_analysis (ng=2) on MVN DGP ===\n')
t_lcmm <- system.time(out_lcmm <- lcmm_analysis(td$trialpaths, dat, op))
print(out_lcmm)
cat('Elapsed:', t_lcmm['elapsed'], 's\n\n')

## --- Heterogeneous-random-slopes DGP smoke test --------------
cat('=== generate_random_slopes_data ===\n')
mp_rs <- data.table(
  N = 70, c.bm = 0,                       # ignored under random-slopes
  carryover_t1half = 1.0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.1, c.cfct = 0.05,
  beta_R_factor = 1.0, beta_NR_factor = 0.0,
  gating_intercept = 0.0, gating_slope = 1.5
)
set.seed(20260507)
dat_rs <- generate_random_slopes_data(
  modelparam = mp_rs, respparam = rp, blparam = bp,
  trialdesign = td$trialpaths[[1]],
  empirical = FALSE, makePositiveDefinite = TRUE
)
dat_rs[, path := 1]
cat('class prevalence (responder fraction):',
    mean(dat_rs$latent_class), '\n')
cat('marginal Cor(B, on-drug response):\n')
on_drug_cols <- paste0('OL', 1:8, '.br')
br_means <- rowMeans(as.matrix(dat_rs[, ..on_drug_cols]))
cat('  cor =', round(cor(dat_rs$bm, br_means), 3), '\n\n')

cat('=== lme_analysis on random-slopes DGP ===\n')
t_lme_rs <- system.time(out_lme_rs <- lme_analysis(td$trialpaths,
                                                    dat_rs, op))
print(out_lme_rs)
cat('Elapsed:', t_lme_rs['elapsed'], 's\n\n')

cat('=== lcmm_analysis on random-slopes DGP ===\n')
t_lcmm_rs <- system.time(out_lcmm_rs <- lcmm_analysis(td$trialpaths,
                                                       dat_rs, op))
print(out_lcmm_rs)
cat('Elapsed:', t_lcmm_rs['elapsed'], 's\n')
