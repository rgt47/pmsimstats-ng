# test_alignment.R
# Validates that 2025 pm_functions.R produces identical output
# to orig R/ functions for the same inputs.

library(devtools)
load_all('.')
library(tidyverse)
source('analysis/2025/pm_functions.R')

cat('=== Level 1: Gompertz ===\n')
t_seq <- seq(0, 30, by = 0.5)
orig_vals <- modgompertz(t_seq, 10.98604, 5, 0.42)
new_vals <- mod_gompertz(t_seq, 10.98604, 5, 0.42)
stopifnot(all.equal(orig_vals, new_vals, tolerance = 1e-12))
cat('PASS: Gompertz functions match\n\n')

cat('=== Level 2: Sigma matrix ===\n')
td_orig <- buildtrialdesign(
  name_longform = 'open label',
  name_shortform = 'OL',
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste0('OL', 1:8),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)

mp <- data.table(
  N = 35, c.bm = 0.3, carryover_t1half = 0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.1, c.cfct = 0.05
)
rp <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)
bp <- data.table(
  cat = c('bm', 'BL'),
  m = c(0, 70),
  sd = c(1, 10)
)

orig_sigma <- buildSigma(mp, rp, bp, td_orig$trialpaths[[1]])

rp_2025 <- tibble(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)
bp_2025 <- tibble(
  cat = c('bm', 'BL'),
  m = c(0, 70),
  sd = c(1, 10)
)
td_2025 <- td_orig$trialpaths[[1]]
td_2025$timepoint_name <- td_2025$timeptnames

new_sigma <- build_sigma_matrix(
  mp, rp_2025, bp_2025, td_2025,
  factor_types = c('tv', 'pb', 'br'),
  factor_abbreviations = c('tv', 'pb', 'br')
)

cat('Orig labels:', paste(orig_sigma$labels[1:5], collapse = ', '), '...\n')
cat('2025 labels:', paste(new_sigma$labels[1:5], collapse = ', '), '...\n')

labels_match <- all(orig_sigma$labels == new_sigma$labels)
cat('Labels match:', labels_match, '\n')

means_match <- all.equal(orig_sigma$means, new_sigma$means,
                         tolerance = 1e-10)
cat('Means match:', means_match, '\n')

sigma_match <- all.equal(orig_sigma$sigma, new_sigma$sigma,
                         tolerance = 1e-10)
cat('Sigma match:', sigma_match, '\n\n')

cat('=== Level 3: Data generation ===\n')
set.seed(42)
orig_dat <- generateData(mp, rp, bp, td_orig$trialpaths[[1]],
                         empirical = FALSE,
                         makePositiveDefinite = TRUE)

set.seed(42)
new_dat <- generate_data(mp, rp_2025, bp_2025, td_2025,
                         empirical = FALSE,
                         make_positive_definite = TRUE)

bm_match <- all.equal(orig_dat$bm, new_dat$bm, tolerance = 1e-10)
cat('bm match:', bm_match, '\n')
bl_match <- all.equal(orig_dat$BL, new_dat$BL, tolerance = 1e-10)
cat('BL match:', bl_match, '\n')
ol1_match <- all.equal(orig_dat$OL1, new_dat$OL1, tolerance = 1e-10)
cat('OL1 match:', ol1_match, '\n\n')

if (isTRUE(bm_match) && isTRUE(bl_match) && isTRUE(ol1_match)) {
  cat('=== ALL LEVELS PASS ===\n')
} else {
  cat('=== FAILURES DETECTED ===\n')
}
