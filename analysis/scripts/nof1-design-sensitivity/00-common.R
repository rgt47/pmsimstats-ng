#' 00-common.R
#'
#' Shared helpers for the sensitivity sweeps (protocol:
#' 00-ademp-protocol.md). Sourced by each of the 01-effect-size.R
#' through 12-morris-null.R scripts and by run-all.R. Designs and
#' baseline parameters are in 00-designs.R.
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
  library(tibble)
  library(dplyr)
})

#' Engine: the installed nof1power must contain the 0.1.0 changes
#' (feasibility skipping, ANCOVA arm contrast, fitted-formula column).
#' Set NOF1POWER_DEV=1 to load the working tree instead.
if (nzchar(Sys.getenv('NOF1POWER_DEV'))) {
  pkgload::load_all(
    here::here('implementations', 'nof1power'),
    quiet = TRUE
  )
} else {
  library(nof1power)
  if (utils::packageVersion('nof1power') < '0.1.0') {
    stop('nof1power >= 0.1.0 is required. Run ',
         'devtools::install_local(\'implementations/nof1power\').')
  }
}

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-designs.R'
))

#=====================================================================
# Paths
#=====================================================================

SENS_DATA_DIR <- Sys.getenv(
  'SENS_OUT_DIR',
  file.path('analysis', 'data', 'derived_data', 'sensitivity')
)
SENS_FIG_DIR <- file.path(
  'analysis', 'figures', 'sensitivity'
)
if (!dir.exists(SENS_DATA_DIR)) {
  dir.create(SENS_DATA_DIR, recursive = TRUE)
}
if (!dir.exists(SENS_FIG_DIR)) {
  dir.create(SENS_FIG_DIR, recursive = TRUE)
}

#' Replicates per cell: the driver default, or the value of the
#' environment variable SENS_NREPS (used only for smoke tests, with
#' SENS_OUT_DIR pointing away from the artifact directory).
sens_reps <- function(default) {
  v <- Sys.getenv('SENS_NREPS')
  if (nzchar(v)) as.integer(v) else as.integer(default)
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
wrap_param_set <- function(param_tbl, name = 'default') {
  list(list(name = name, param = param_tbl))
}

#' Analysis-parameter rows, one per analysis half-life.
#'
#' @param target_coef `'main'` (the `Dbc` coefficient; for a parallel
#'   design, the ANCOVA arm contrast) or `'interaction'`.
#' @param analysis_t1half Analysis half-lives in weeks used to build
#'   `Dbc`. The project policy value when the half-life is fixed is
#'   0.2 week. 0 gives the binary on/off indicator.
default_analysisparams <- function(target_coef = 'main',
                                   analysis_t1half = 0.2) {
  data.frame(
    useDE = TRUE,
    t_random_slope = FALSE,
    full_model_out = FALSE,
    carryover_t1half = analysis_t1half,
    simplecarryover = FALSE,
    carryover_scalefactor = 1,
    target_coef = target_coef,
    parallel_analysis = 'ancova',
    stringsAsFactors = FALSE
  )
}

#' Build a `simparam` list with progressive save disabled.
default_simparam <- function(n_reps = 500,
                             basesavename = 'sensitivity',
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

#' Number of workers: all but one logical core.
sweep_cores <- function() {
  max(1, parallel::detectCores() - 1)
}

#=====================================================================
# Running a sweep
#=====================================================================

#' One call to the engine with infeasible cells skipped (never
#' repaired) and design labels attached.
#'
#' @param design_names Labels in the order of `designs`.
#' @param design_ids Integer identifiers written to `trialdesign`
#'   so that results of several calls can be combined.
sweep_call <- function(designs, design_names, respparamsets,
                       blparamsets, modelparams, n_reps,
                       analysisparams = default_analysisparams(),
                       dgp_architecture = 'mvn',
                       design_ids = seq_along(designs)) {
  res <- generate_simulated_results(
    trialdesigns = designs,
    respparamsets = respparamsets,
    blparamsets = blparamsets,
    censorparams = NA,
    modelparams = modelparams,
    simparam = default_simparam(n_reps = n_reps),
    analysisparams = analysisparams,
    rawdataout = FALSE,
    dgp_architecture = dgp_architecture,
    n_cores = sweep_cores(),
    skip_infeasible = TRUE
  )
  local_td <- res$results$trialdesign
  res$results$design_name <- design_names[local_td]
  res$results$trialdesign <- design_ids[local_td]
  res$results$dgp_architecture <- dgp_architecture
  res$cells$design_name <- design_names[res$cells$trialdesign]
  res$cells$trialdesign <- design_ids[res$cells$trialdesign]
  res$infeasible <- res$cells[!res$cells$feasible, , drop = FALSE]
  res$cells$dgp_architecture <- dgp_architecture
  res$infeasible$dgp_architecture <- dgp_architecture
  res
}

#' Combine the results of several `sweep_call()`s.
combine_sweeps <- function(parts, tag = NULL) {
  bind <- function(el) {
    dplyr::bind_rows(lapply(parts, `[[`, el), .id = tag)
  }
  list(results = bind('results'), cells = bind('cells'),
       infeasible = bind('infeasible'))
}

#' Stop unless the stored values of a column equal the grid.
#'
#' Guards against the S9 defect, where a driver built a grid but
#' passed only its first row, so every cell was run (and then
#' labeled) at one value.
#'
#' @param df Stored per-replicate results.
#' @param column Column that must hold the swept parameter as the
#'   engine stored it (not a label attached afterwards).
#' @param expected The grid, as intended.
#' @param among Optional logical vector selecting the rows to check.
assert_stored_grid <- function(df, column, expected, among = NULL) {
  if (!is.null(among)) df <- df[among, , drop = FALSE]
  if (is.character(expected)) {
    stored <- sort(unique(df[[column]]))
    want <- sort(unique(expected))
    if (!identical(stored, want)) {
      stop(sprintf(
        'Stored values of `%s` (%s) differ from the grid (%s).',
        column, paste(stored, collapse = ', '),
        paste(want, collapse = ', ')
      ), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  stored <- sort(unique(round(df[[column]], 8)))
  want <- sort(unique(round(expected, 8)))
  if (length(stored) != length(want) ||
      any(abs(stored - want) > 1e-7)) {
    stop(sprintf(
      paste0('Stored values of `%s` (%s) differ from the grid ',
             '(%s): %d distinct stored, %d in the grid.'),
      column, paste(stored, collapse = ', '),
      paste(want, collapse = ', '),
      length(stored), length(want)
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless run cells plus dropped cells equal the intended grid.
assert_cell_count <- function(res, expected_cells) {
  n_run <- nrow(dplyr::distinct(
    res$results, trialdesign, respparamset, blparamset,
    modelparamset, dgp_architecture
  ))
  n_drop <- nrow(res$infeasible)
  if (n_run + n_drop != expected_cells) {
    stop(sprintf(
      'Ran %d cells and dropped %d; the grid has %d.',
      n_run, n_drop, expected_cells
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Save a sweep result with run metadata.
#'
#' The metadata record is what report.Rmd checks before it reads the
#' artifact (schema, replicate count, seed, engine version).
#'
#' @param checks Optional function of the result that stops on a
#'   failed grid assertion (see `assert_stored_grid()`).
save_sweep <- function(result, id, name, n_reps = NA, seed = 42 + id,
                       notes = NULL, checks = NULL) {
  rev <- tryCatch(
    system2('git', c('rev-parse', '--short', 'HEAD'),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_, warning = function(w) NA_character_
  )
  result$run_meta <- list(
    schema = 2L,
    sweep_id = id,
    sweep_name = name,
    n_reps = n_reps,
    seed = seed,
    base_N = BASE_N,
    engine_version = as.character(utils::packageVersion('nof1power')),
    git_rev = if (length(rev)) rev[1] else NA_character_,
    r_version = R.version.string,
    date = format(Sys.time(), '%Y-%m-%d %H:%M %Z'),
    notes = notes
  )
  out_path <- file.path(
    SENS_DATA_DIR,
    sprintf('sensitivity-%02d-%s.rds', id, name)
  )
  if (!is.null(checks)) {
    # A failed assertion must not discard hours of simulation: the
    # unchecked artifact is kept under a name report.Rmd ignores.
    ok <- tryCatch({ checks(result); TRUE }, error = function(e) e)
    if (inherits(ok, 'error')) {
      bad <- sub('\\.rds$', '-UNCHECKED.rds', out_path)
      saveRDS(result, bad)
      stop('Grid assertion failed (', conditionMessage(ok),
           '). Results kept in ', bad, call. = FALSE)
    }
  }
  saveRDS(result, out_path)
  message('Saved: ', out_path)
  if (!is.null(result$infeasible) && nrow(result$infeasible) > 0) {
    message('Infeasible cells dropped: ', nrow(result$infeasible),
            ' (see $infeasible in the artifact)')
  }
  invisible(out_path)
}

#' Response-parameter wrapper matched on the realized estimand.
#'
#' The engine sets the drug-response ceiling per design so that
#' `delta_on` equals `target_delta` (see `generate_simulated_results`).
#'
#' @param target_delta Mean drug response on drug, in outcome units
#'   (negative lowers symptoms; 0 is the null).
resp_set <- function(rp = baseline_resp_param(),
                     target_delta = BASE_DELTA,
                     name = sprintf('delta=%g', target_delta)) {
  list(name = name, param = rp, target_delta = target_delta)
}
