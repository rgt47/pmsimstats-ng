## analysis/scripts/carryover-sensitivity/30-check-prognostic-confound-pd.R
##
## Diagnostic for the prognostic-confound cells of jobs 28 and 30
## (28-run-g9-prognostic-stress.R, 29-run-g9-random-effects.R).
## Those drivers run with stop_on_nonpd = FALSE, so a non-positive-
## definite Sigma is silently repaired by
## corpcor::make.positive.definite() inside build_sigma_matrix().
## This script asks, for every (pathway, strength, c_bm) combination
## those jobs used, two questions that need no simulation at all:
##
##   1. Is the raw Sigma positive definite on each Hybrid path?
##      (min eigenvalue of the raw sd * correlation * sd matrix, and
##      build_sigma_matrix()'s own was_pd_corrected flag).
##   2. What is the population on-drug minus off-drug slope contrast
##      of Sx on bm, the quantity bm:Dbc estimates, before and after
##      the repair? If the repair changes it, the "power loss" seen
##      at that cell is attenuation of the generated signal, not a
##      property of any analysis model.
##
## The population slope of Sx_t on bm at timepoint t is derived from
## Sigma directly: Sx_t = BL - tv_t - pb_t - br_t, so
## slope_t = (Cov(bm, BL) - Cov(bm, tv_t) - Cov(bm, pb_t)
##            - Cov(bm, br_t)) / Var(bm).
## The contrast averages slope_t over on-drug timepoints (tod > 0,
## tsd == 0) minus off-drug timepoints, per path, then over paths.
##
## Result (2026-10-26, verified): every strength-0.6 cell on the
## prog-PB and prog-TV-const pathways, and the prog-BL cell at
## c_bm = 0.45, is non-PD on all four Hybrid paths and is repaired;
## the repair shrinks the contrast from about -3.3 to -2.4 (prog-PB)
## and -1.8 (prog-TV-const), matching the attenuated point estimates
## jobs 28 and 30 report. Strength-0.3 cells are PD on every path.
## prog-TV-grow is a mean shift (b.bm.tv), not a covariance entry,
## and never affects PD. See docs/60 (2026-10-26 update) and docs/64.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/30-check-prognostic-confound-pd.R
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

OUT_DIR <- 'analysis/data/carryover-sensitivity/g9-prognostic-stress'

design_set <- build_design_set('Hybrid')
resp_param <- default_resp_param()
baseline_param <- default_baseline_param()

## Same model_param construction as simulate_cell_combined() /
## simulate_cell_g9_random() in simulation-core.R (Hybrid, t1/2 = 0.5,
## exponential carryover, rho = 0.7, c_1 = 0.2). N does not enter
## Sigma, so one value suffices.
pathway_field <- c(none = NA, progbl = 'c.bm.bl', progtvconst = 'c.bm.tv',
                   progtvgrow = 'b.bm.tv', progpb = 'c.bm.pb')

mk_param <- function(c_bm, field, strength) {
  mp <- list(N = 70, carryover_t1half = 0.5, carryover_form = 'exponential',
             weibull_shape = 1, dgp_architecture = 'mvn', c.bm = c_bm,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
             c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0, c.bm.pb = 0)
  if (!is.na(field)) mp[[field]] <- strength
  mp
}

slopes <- function(sigma, labels, td) {
  i_bm <- which(labels == 'bm')
  i_bl <- which(labels == 'BL')
  v_bm <- sigma[i_bm, i_bm]
  map_dbl(td$timepoint_name, function(tp) {
    idx <- which(labels %in% paste(tp, c('tv', 'pb', 'br'), sep = '.'))
    (sigma[i_bm, i_bl] - sum(sigma[i_bm, idx])) / v_bm
  })
}

grid <- expand_grid(pathway = names(pathway_field),
                    strength = c(0.3, 0.6),
                    c_bm = c(0, 0.45)) |>
  filter(!(pathway == 'none' & strength == 0.6)) |>
  mutate(strength = ifelse(pathway == 'none', 0, strength),
         field = unname(pathway_field[pathway]))

per_path <- pmap_dfr(grid, function(pathway, strength, c_bm, field) {
  mp <- mk_param(c_bm, field, strength)
  map_dfr(seq_along(design_set), function(i) {
    td <- design_set[[i]]
    b <- build_sigma_matrix(mp, resp_param, baseline_param, td)
    raw <- outer(b$standard_deviations, b$standard_deviations) *
      b$correlations
    ev <- eigen(raw, symmetric = TRUE, only.values = TRUE)$values
    on <- td$tod > 0 & td$tsd == 0
    s_raw <- slopes(raw, b$labels, td)
    s_rep <- slopes(b$sigma, b$labels, td)
    tibble(pathway, strength, c_bm, path = i,
           repaired = b$was_pd_corrected,
           min_eigenvalue = min(ev),
           contrast_raw = mean(s_raw[on]) - mean(s_raw[!on]),
           contrast_repaired = mean(s_rep[on]) - mean(s_rep[!on]))
  })
})

stopifnot(nrow(per_path) == nrow(grid) * length(design_set),
          !anyNA(per_path$contrast_raw))

summary_tbl <- per_path |>
  group_by(pathway, strength, c_bm) |>
  summarise(paths_repaired = sum(repaired),
            min_eigenvalue = min(min_eigenvalue),
            contrast_raw = mean(contrast_raw),
            contrast_repaired = mean(contrast_repaired),
            .groups = 'drop') |>
  mutate(attenuation = ifelse(abs(contrast_raw) < 1e-6, NA_real_,
                              contrast_repaired / contrast_raw))

cat('Prognostic-confound cells (Hybrid, 4 paths): PD status and\n',
    'population on/off slope contrast before and after PD repair\n\n')
print(as.data.frame(summary_tbl |>
  mutate(across(c(min_eigenvalue, contrast_raw, contrast_repaired,
                  attenuation), ~ round(.x, 3)))),
  row.names = FALSE)

saveRDS(list(per_path = per_path, summary = summary_tbl,
             meta = list(date = Sys.time(),
                         design = 'Hybrid', t1half = 0.5, rho = 0.7)),
        file.path(OUT_DIR, 'pd-check.rds'))
cat(sprintf('\nWritten: %s\n', file.path(OUT_DIR, 'pd-check.rds')))
