#' 00-designs.R
#'
#' Trial designs, baseline parameters and the realized estimand for
#' paper 05. Sourced by 00-common.R (sweeps) and by report.Rmd
#' (design table), so it holds no simulation or file-system side
#' effects.
#'
#' Design set (all use the same eight-visit, 20-week calendar):
#'   Hybrid    Published Hendrickson et al. (2020) design: visits at
#'             weeks 4, 8, 9, 10, 11, 12, 16, 20; four paths from two
#'             randomizations; open-label expectancy 1, blinded 0.5.
#'   Compact2  The two-path, eight-weekly-visit schedule used in the
#'             2026-05 runs of this study (retained for comparison
#'             only; it is not the published design).
#'   Alt       Alternating on/off at every visit, two paths.
#'   XO        Two-period crossover, one assessment per period.
#'   RCT       Parallel groups, two arms, blinded (expectancy 0.5).
#'   RCT-OL    Parallel groups with expectancy 1 at every visit (the
#'             assumption of the 2026-05 runs); quantifies the effect
#'             of the comparator blinding assumption.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
})

#=====================================================================
# Baseline parameters (prazosin-for-PTSD calibration)
#=====================================================================

#' Baseline sample size (total across paths or arms).
BASE_N <- 35L

#' Baseline effect: mean drug response on drug, in outcome units,
#' common to all designs in a sweep unless the sweep varies it.
BASE_DELTA <- -0.5

#' Visit weeks shared by every design in the primary comparison.
STD_WEEKS <- c(4, 8, 9, 10, 11, 12, 16, 20)

#' Baseline response parameters for the natural-history (tv),
#' placebo-response (pb) and drug-response (br) components.
baseline_resp_param <- function(br_max = 2.0) {
  tibble(
    cat = c('tv', 'pb', 'br'),
    max = c(0.5, 1.0, br_max),
    disp = c(1.0, 1.0, 1.0),
    rate = c(0.3, 0.3, 0.3),
    sd = c(0.5, 0.5, 1.0)
  )
}

#' Biomarker and baseline-symptom means and SDs.
baseline_bl_param <- function() {
  tibble(
    cat = c('bm', 'BL'),
    m = c(120, 8.0),
    sd = c(15, 1.5)
  )
}

#' Baseline model-parameter data frame (one row).
#'
#' `carryover_t1half` is in weeks (3 / 7 is the 3-day baseline).
#' `tau` is the SD of participant-level treatment-effect
#' heterogeneity (0 means none; see `generate_data()`).
baseline_model_param <- function(N = 35, c.bm = 0.3,
                                 carryover_t1half = 3 / 7,
                                 c.tv = 0.3, c.pb = 0.3,
                                 c.br = 0.3,
                                 c.cf1t = 0.1, c.cfct = 0.05,
                                 carryover_form = 'exponential',
                                 weibull_shape = 1, tau = 0) {
  data.frame(
    N = N, c.bm = c.bm, tau = tau,
    carryover_t1half = carryover_t1half,
    c.tv = c.tv, c.pb = c.pb, c.br = c.br,
    c.cf1t = c.cf1t, c.cfct = c.cfct,
    carryover_form = carryover_form,
    weibull_shape = weibull_shape,
    stringsAsFactors = FALSE
  )
}

#=====================================================================
# Designs
#=====================================================================

#' Published Hendrickson et al. (2020) Hybrid design.
#'
#' Open label at weeks 4 and 8 (expectancy 1), blinded
#' discontinuation at weeks 9 to 12 (expectancy 0.5), blinded
#' crossover at weeks 16 and 20 (expectancy 0.5). Four paths: the
#' drug stops after week 9 or after week 10, and the crossover is
#' drug first or placebo first (docs/40, Section 2).
design_hybrid <- function() {
  nof1power::build_trial_design(
    name_longform = 'Hendrickson (2020) Hybrid',
    name_shortform = 'Hybrid',
    timepoints = STD_WEEKS,
    timeptnames = c('OL1', 'OL2', 'BD1', 'BD2',
                    'BD3', 'BD4', 'COd', 'COp'),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    ondrug = list(
      c(1, 1, 1, 1, 0, 0, 1, 0),
      c(1, 1, 1, 1, 0, 0, 0, 1),
      c(1, 1, 1, 0, 0, 0, 1, 0),
      c(1, 1, 1, 0, 0, 0, 0, 1)
    )
  )
}

#' Alias kept for the drivers: the primary N-of-1 design.
design_nof1 <- function() design_hybrid()

#' The two-path, eight-weekly-visit schedule of the 2026-05 runs.
#'
#' Not the published Hybrid. Two paths, two open-label visits.
design_compact2 <- function() {
  nof1power::build_trial_design(
    name_longform = 'Compact two-path schedule',
    name_shortform = 'Compact2',
    timepoints = seq(1, 8, 1),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    ondrug = list(
      c(1, 1, 1, 1, 0, 0, 1, 1),
      c(1, 1, 1, 1, 1, 1, 0, 0)
    )
  )
}

#' Parallel-group design: two arms, drug for all visits or never.
design_rct <- function(weeks = STD_WEEKS, expectancy = 0.5,
                       short = 'RCT') {
  nof1power::build_trial_design(
    name_longform = 'Parallel RCT',
    name_shortform = short,
    timepoints = weeks,
    expectancies = rep(expectancy, length(weeks)),
    ondrug = list(
      rep(1, length(weeks)),
      rep(0, length(weeks))
    )
  )
}

#' Parallel-group design with open-label expectancy (2026-05 setting).
design_rct_openlabel <- function(weeks = STD_WEEKS) {
  design_rct(weeks, expectancy = 1, short = 'RCT-OL')
}

#' Aggregated N-of-1, alternating on/off at every visit.
design_nof1_alternating <- function(weeks = STD_WEEKS,
                                    expectancy = 0.5) {
  n <- length(weeks)
  nof1power::build_trial_design(
    name_longform = 'Aggregated N-of-1 (alternating)',
    name_shortform = 'Alt',
    timepoints = weeks,
    expectancies = rep(expectancy, n),
    ondrug = list(
      rep(c(1, 0), length.out = n),
      rep(c(0, 1), length.out = n)
    )
  )
}

#' Aggregated N-of-1 with `n_periods` equal periods and one visit per
#' week for `n_visits` weeks (fixed trial length and fixed number of
#' assessments; the number of switches varies).
design_nof1_periods <- function(n_periods, n_visits = 16,
                                expectancy = 0.5) {
  stopifnot(n_visits %% n_periods == 0)
  len <- n_visits / n_periods
  pattern <- rep(rep(c(1, 0), length.out = n_periods), each = len)
  nof1power::build_trial_design(
    name_longform = sprintf('Alternating, %d periods', n_periods),
    name_shortform = sprintf('Alt-P%d', n_periods),
    timepoints = seq_len(n_visits),
    expectancies = rep(expectancy, n_visits),
    ondrug = list(pattern, 1 - pattern)
  )
}

#' Parallel RCT with `n_visits` weekly visits.
design_rct_weekly <- function(n_visits = 16, expectancy = 0.5) {
  design_rct(seq_len(n_visits), expectancy = expectancy)
}

#' Two-period crossover with one assessment at the end of each
#' period.
design_crossover <- function(weeks = c(10, 20), expectancy = 0.5) {
  nof1power::build_trial_design(
    name_longform = 'Two-period crossover',
    name_shortform = 'XO',
    timepoints = weeks,
    expectancies = rep(expectancy, length(weeks)),
    ondrug = list(c(1, 0), c(0, 1))
  )
}

#=====================================================================
# Realized estimand
#=====================================================================

#' Realized estimand of a design: mean drug response on drug.
#'
#' The mean, over every on-drug follow-up visit of every path
#' (paths weighted equally, as when N is divisible by the number of
#' paths), of the mean drug-response component, with the sign of the
#' symptom change. It depends on the design and on the response
#' curve only, so it is common to all designs and is recorded in
#' every stored artifact as `delta_on`.
#'
#' @param td Design from nof1power::build_trial_design().
#' @param rp Response-parameter tibble.
#' @param bp Baseline-parameter tibble.
#' @param mp One-row model-parameter data frame.
design_delta_on <- function(td, rp = baseline_resp_param(),
                            bp = baseline_bl_param(),
                            mp = baseline_model_param()) {
  vals <- lapply(td$trialpaths, function(path) {
    sg <- nof1power::build_sigma_matrix(
      mp, rp, bp, path, stop_on_nonpd = FALSE
    )
    on <- path$tod > 0
    list(sum = sum(utils::tail(sg$means, nrow(path))[on]),
         n = sum(on))
  })
  -sum(vapply(vals, `[[`, 0, 'sum')) /
    sum(vapply(vals, `[[`, 0, 'n'))
}

#' Drug-response ceiling `m` that gives a target `delta_on`.
#'
#' `delta_on` is linear in `m`, so one evaluation at `m = 1` fixes it.
br_max_for_delta <- function(td, target_delta,
                             rp = baseline_resp_param(),
                             bp = baseline_bl_param(),
                             mp = baseline_model_param()) {
  rp$max[rp$cat == 'br'] <- 1
  unit <- design_delta_on(td, rp, bp, mp)
  abs(target_delta) / abs(unit)
}

#' Visit-by-visit description of a design for the manuscript.
describe_design <- function(td) {
  purrr::imap_dfr(td$trialpaths, function(tp, q) {
    tibble(
      path = q,
      visit = tp$timeptnames,
      week = cumsum(tp$t_wk),
      drug = ifelse(tp$tod > 0, 'on', 'off'),
      expectancy = tp$e,
      t_od = tp$tod,
      t_sd = tp$tsd
    )
  })
}
