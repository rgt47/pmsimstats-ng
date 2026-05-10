#' 00-common.R
#'
#' Shared fixtures and helpers for the sensitivity sweeps described
#' in docs/sensitivity-sweep-plan-2026-04-17.md. Sourced by each of
#' the 01-effect-size.R through 12-morris-null.R scripts and by
#' run-all.R.
#'
#' All baseline parameters are calibrated to the prazosin-for-PTSD
#' context that motivates the primary simulation; parameter sweeps
#' perturb one axis at a time while holding the rest at the baseline.
#'
#' Invariants:
#'   - Column-name convention is pmsimstats short-form:
#'     `cat`, `m`, `sd` for baseline; `cat`, `max`, `disp`, `rate`,
#'     `sd` for response; `bm`, `BL`, `ptID`, `timeptnames` in data.
#'   - Every script sets its own RNG seed via `set.seed(42 + id)`.

suppressPackageStartupMessages({
  library(nof1power)
  library(tibble)
  library(dplyr)
})

#=====================================================================
# Paths
#=====================================================================

SENS_DATA_DIR <- file.path(
  "analysis", "data", "derived_data", "sensitivity"
)
SENS_FIG_DIR <- file.path(
  "analysis", "figures", "sensitivity"
)
if (!dir.exists(SENS_DATA_DIR)) {
  dir.create(SENS_DATA_DIR, recursive = TRUE)
}
if (!dir.exists(SENS_FIG_DIR)) {
  dir.create(SENS_FIG_DIR, recursive = TRUE)
}

#=====================================================================
# Baseline parameters (prazosin-for-PTSD calibration)
#=====================================================================

#' Baseline response-parameter tibble
#'
#' Three latent factors (time-variant, pharm-biomarker, bio-response)
#' parameterized with modified-Gompertz trajectories plus an
#' independent Normal error SD.
baseline_resp_param <- function() {
  tibble(
    cat = c("tv", "pb", "br"),
    max = c(0.5, 1.0, 2.0),
    disp = c(1.0, 1.0, 1.0),
    rate = c(0.3, 0.3, 0.3),
    sd = c(0.5, 0.5, 1.0)
  )
}

#' Baseline baseline-parameter tibble
#'
#' Biomarker and baseline-symptom means and SDs.
baseline_bl_param <- function() {
  tibble(
    cat = c("bm", "BL"),
    m = c(120, 8.0),
    sd = c(15, 1.5)
  )
}

#' Baseline model-parameter data.frame (single-row)
#'
#' Correlation coefficients follow pmsimstats short-form:
#'   c.tv, c.pb, c.br: within-factor AR(1) over timepoints.
#'   c.cf1t: cross-factor correlation at matched timepoints.
#'   c.cfct: cross-factor correlation at lagged timepoints.
#'   c.bm: biomarker-bio-response correlation (Architecture B)
#'         or regression coefficient (Architecture A).
#'   carryover_t1half: carryover half-life in the same time units
#'         as the trial design (weeks in the prazosin baseline).
baseline_model_param <- function(N = 35, c.bm = 0.3,
                                  carryover_t1half = 3 / 7,
                                  c.tv = 0.3, c.pb = 0.3,
                                  c.br = 0.3,
                                  c.cf1t = 0.1, c.cfct = 0.05,
                                  carryover_form = "exponential",
                                  weibull_shape = 1) {
  data.frame(
    N = N, c.bm = c.bm,
    carryover_t1half = carryover_t1half,
    c.tv = c.tv, c.pb = c.pb, c.br = c.br,
    c.cf1t = c.cf1t, c.cfct = c.cfct,
    carryover_form = carryover_form,
    weibull_shape = weibull_shape,
    stringsAsFactors = FALSE
  )
}

#=====================================================================
# Trial designs
#=====================================================================

#' Parallel RCT represented as two always-on / always-off paths
design_rct <- function(weeks = seq(1, 8, 1)) {
  build_trial_design(
    name_longform = "Parallel RCT",
    name_shortform = "RCT",
    timepoints = weeks,
    expectancies = rep(1, length(weeks)),
    ondrug = list(
      rep(1, length(weeks)),
      rep(0, length(weeks))
    )
  )
}

#' Aggregated N-of-1, alternating A/B/A/B variant
#'
#' Kept as a sensitivity comparator; the Hendrickson 2020 hybrid
#' is the primary aggregated-N-of-1 design (see `design_nof1`).
design_nof1_alternating <- function(weeks = seq(1, 8, 1)) {
  n <- length(weeks)
  seq1 <- rep(c(1, 0), length.out = n)
  seq2 <- rep(c(0, 1), length.out = n)
  build_trial_design(
    name_longform = "Aggregated N-of-1 (alternating)",
    name_shortform = "Nof1-alt",
    timepoints = weeks,
    expectancies = rep(1, n),
    ondrug = list(seq1, seq2)
  )
}

#' Primary aggregated N-of-1 design (Hendrickson 2020 hybrid)
#'
#' The default aggregated-N-of-1 design used throughout the
#' sensitivity sweeps is the Hendrickson 2020 hybrid
#' (open-label titration + blinded discontinuation + brief
#' crossover). See `design_hybrid()` for the explicit
#' construction. The older alternating A/B/A/B variant is
#' available as `design_nof1_alternating()`.
design_nof1 <- function() {
  design_hybrid()
}

#' Two-period crossover
design_crossover <- function() {
  build_trial_design(
    name_longform = "Two-period crossover",
    name_shortform = "XO",
    timepoints = c(4, 8),
    expectancies = c(1, 1),
    ondrug = list(c(1, 0), c(0, 1))
  )
}

#' Hybrid OL + blinded-discontinuation + brief crossover
#' (Hendrickson 2020 design 4)
design_hybrid <- function() {
  build_trial_design(
    name_longform = "OL + BDC + crossover",
    name_shortform = "Hybrid",
    timepoints = seq(1, 8, 1),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    ondrug = list(
      c(1, 1, 1, 1, 0, 0, 1, 1),
      c(1, 1, 1, 1, 1, 1, 0, 0)
    )
  )
}

#=====================================================================
# Cycle-period design grammar (S13)
#=====================================================================

#' Build a cycle-period N-of-1 trial design at protocol level.
#'
#' Higher-level interface to `nof1power::build_trial_design()`
#' that exposes the trial-protocol parameters directly: number of
#' blinded on/off cycles, period lengths, optional pre-treatment
#' run-in, optional open-label titration, optional post-cycle
#' run-out. Returns the same object as `build_trial_design()`,
#' so the downstream simulation pipeline consumes the result
#' unchanged.
#'
#' Convention (per docs/27-cycle-period-design-sweep-plan.md):
#' - **Run-in**: pre-treatment baseline period, no drug.
#'   Observations have $\text{Db} = 0$, contributing directly to
#'   the off-drug stratum.
#' - **Open-label**: on-drug, unblinded titration period.
#'   Observations have $\text{Db} = 1$.
#' - **Cycles**: blinded alternation of on-drug and off-drug
#'   periods. The `order` parameter sets which phase begins each
#'   cycle.
#' - **Run-out**: post-cycle observation period (optionally
#'   off-drug).
#'
#' @param T_runin        weeks of pre-treatment baseline (no drug)
#' @param T_openlabel    weeks of on-drug, unblinded titration
#' @param k              number of blinded on/off cycles
#' @param T_on           blinded on-drug period (weeks); scalar
#'   or length-k vector for across-cycle asymmetry
#' @param T_off          blinded off-drug period (weeks); scalar
#'   or length-k
#' @param T_runout       weeks of post-cycle observation
#' @param order          'on_first' or 'off_first' (which phase
#'   begins each cycle)
#' @param density        observations per week (default 1)
#' @param expectancy_blinded  blinded-phase expectancy (default 0.5)
#' @return list(metadata, trialpaths) per `build_trial_design()`
build_cycle_design <- function(T_runin = 0,
                               T_openlabel = 0,
                               k,
                               T_on, T_off,
                               T_runout = 0,
                               order = c("on_first", "off_first"),
                               density = 1,
                               expectancy_blinded = 0.5) {
  order <- match.arg(order)
  if (length(T_on)  == 1) T_on  <- rep(T_on,  k)
  if (length(T_off) == 1) T_off <- rep(T_off, k)
  stopifnot(length(T_on) == k, length(T_off) == k)
  stopifnot(density > 0)

  ## Build the phase sequence as a list of (duration, ondrug,
  ## expectancy) triples; convert to per-week timepoints below.
  phases <- list()
  if (T_runin > 0) {
    phases[[length(phases) + 1L]] <- list(
      duration = T_runin, ondrug = 0L, expectancy = 0)
  }
  if (T_openlabel > 0) {
    phases[[length(phases) + 1L]] <- list(
      duration = T_openlabel, ondrug = 1L, expectancy = 1)
  }
  for (j in seq_len(k)) {
    if (order == "on_first") {
      phases[[length(phases) + 1L]] <- list(
        duration = T_on[j], ondrug = 1L,
        expectancy = expectancy_blinded)
      phases[[length(phases) + 1L]] <- list(
        duration = T_off[j], ondrug = 0L,
        expectancy = expectancy_blinded)
    } else {
      phases[[length(phases) + 1L]] <- list(
        duration = T_off[j], ondrug = 0L,
        expectancy = expectancy_blinded)
      phases[[length(phases) + 1L]] <- list(
        duration = T_on[j], ondrug = 1L,
        expectancy = expectancy_blinded)
    }
  }
  if (T_runout > 0) {
    phases[[length(phases) + 1L]] <- list(
      duration = T_runout, ondrug = 0L,
      expectancy = expectancy_blinded)
  }

  ## Discretise each phase into observations at the requested
  ## density. An observation is placed at the *end* of each
  ## sub-period so that, e.g., a 4-week phase at density 1
  ## yields observations at weeks 1, 2, 3, 4 of the phase.
  weeks <- numeric()
  ondrug_vec <- integer()
  expect_vec <- numeric()
  cum_t <- 0
  for (ph in phases) {
    n_obs <- max(1L, as.integer(round(ph$duration * density)))
    dt <- ph$duration / n_obs
    for (i in seq_len(n_obs)) {
      cum_t <- cum_t + dt
      weeks <- c(weeks, cum_t)
      ondrug_vec <- c(ondrug_vec, ph$ondrug)
      expect_vec <- c(expect_vec, ph$expectancy)
    }
  }

  build_trial_design(
    name_longform = sprintf(
      "cycle k=%d Ton=%s Toff=%s runin=%g OL=%g runout=%g",
      k,
      paste(T_on, collapse = ","),
      paste(T_off, collapse = ","),
      T_runin, T_openlabel, T_runout),
    name_shortform = sprintf("C%d", k),
    timepoints = weeks,
    expectancies = expect_vec,
    ondrug = list(ondrug_vec)
  )
}

#=====================================================================
# Convenience wrappers
#=====================================================================

#' Wrap a response-parameter tibble and a baseline-parameter tibble
#' in the list-of-lists shape expected by
#' `generate_simulated_results()`.
wrap_param_set <- function(param_tbl, name = "default") {
  list(list(name = name, param = param_tbl))
}

#' Build an `analysisparams` data.frame with sensible defaults.
#'
#' @param target_coef Either `"main"` (default; extract Dbc main-effect
#'   coefficient) or `"interaction"` (extract bm:Dbc biomarker-by-
#'   treatment interaction coefficient). The nof1_power sensitivity
#'   sweeps target the main treatment effect by default to align
#'   with the majority of the N-of-1 simulation-comparison
#'   literature; interaction-targeted sweeps are archived under
#'   `analysis/archive/interaction-2026-04-18/`.
default_analysisparams <- function(target_coef = "main") {
  data.frame(
    useDE = TRUE,
    t_random_slope = FALSE,
    full_model_out = FALSE,
    carryover_t1half = 0,
    simplecarryover = FALSE,
    carryover_scalefactor = 1,
    target_coef = target_coef,
    stringsAsFactors = FALSE
  )
}

#' Build a `simparam` list with progressive save disabled.
default_simparam <- function(n_reps = 500,
                              basesavename = "sensitivity",
                              savedir = tempdir()) {
  list(
    Nreps = n_reps,
    progressiveSave = FALSE,
    basesavename = basesavename,
    nRep2save = .Machine$integer.max,
    saveunit2start = 1,
    savedir = savedir
  )
}

#' Save a sweep result to the sensitivity directory.
save_sweep <- function(result, id, name) {
  out_path <- file.path(
    SENS_DATA_DIR,
    sprintf("sensitivity-%02d-%s.rds", id, name)
  )
  saveRDS(result, out_path)
  message("Saved: ", out_path)
  invisible(out_path)
}
