#' results-helpers.R
#'
#' Reads the stored sweep artifacts and builds the Results text,
#' tables and figure paths for report.Rmd. Every number in Results
#' comes from these functions. A sweep whose artifact is missing, or
#' was written by an earlier version of the drivers (no schema-2
#' metadata), yields a RESULTS-PENDING note with the rerun command
#' instead of numbers, so the document renders either way.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(knitr)
})

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-summarize.R'
))
source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-designs.R'
))

DATA_DIR <- Sys.getenv(
  'SENS_OUT_DIR',
  here::here('analysis', 'data', 'derived_data', 'sensitivity')
)
FIG_DIR <- Sys.getenv(
  'SENS_FIG_DIR',
  here::here('analysis', 'figures', 'sensitivity')
)

#=====================================================================
# Sweep specification (planned grids; counts are checked against the
# artifacts when they exist)
#=====================================================================

SPEC <- tibble::tribble(
  ~id, ~name, ~label, ~reps, ~cells,
  1L, 'effect-size', 'Effect size', 1000L, 42L,
  2L, 'carryover-halflife', 'Carryover half-life', 500L, 36L,
  3L, 'decay-form', 'Decay form', 500L, 15L,
  4L, 'tau2', 'Heterogeneity and baseline SD', 500L, 30L,
  5L, 'sigma2', 'Within-participant SD', 500L, 15L,
  6L, 'cycles', 'Periods at fixed length', 500L, 5L,
  7L, 'patients', 'Sample size', 500L, 18L,
  8L, 'ar1', 'AR(1) correlation', 500L, 15L,
  9L, 'period', 'Period length by half-life', 500L, 15L,
  10L, 'biomarker', 'Biomarker correlation by encoding', 500L, 48L,
  11L, 'misspec', 'Misspecification factorial', 500L, 24L,
  12L, 'morris-null', 'Null calibration', 2000L, 16L
)

art_file <- function(id, name) {
  file.path(DATA_DIR, sprintf('sensitivity-%02d-%s.rds', id, name))
}

.cache <- new.env()

#' Load a schema-2 artifact, or NULL.
load_art <- function(id, name = SPEC$name[SPEC$id == id]) {
  key <- sprintf('%d', id)
  if (!is.null(.cache[[key]])) return(.cache[[key]])
  path <- art_file(id, name)
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (!identical(obj$run_meta$schema, 2L)) return(NULL)
  r <- obj$results
  r$halflife_days <- round(r$carryover_t1half * 7, 4)
  for (col in c('delta_on', 'c.bm', 'c.br', 'tau', 'br_sd', 'bl_sd',
                'analysis_t1half')) {
    if (col %in% names(r)) r[[col]] <- round(r[[col]], 6)
  }
  obj$results <- r
  .cache[[key]] <- obj
  obj
}

avail <- function(id) !is.null(load_art(id))

#' Summaries of a sweep by the given grouping columns.
smry <- function(id, by) summarize_cells(load_art(id)$results, by)

#=====================================================================
# Formatting
#=====================================================================

f2 <- function(x) sprintf('%.2f', x)
f3 <- function(x) sprintf('%.3f', x)
f1 <- function(x) sprintf('%.1f', x)
pm <- function(p, mc) sprintf('%.3f (MCSE %.3f)', p, mc)
join_txt <- function(x) {
  x <- x[nzchar(x)]
  if (length(x) <= 1) return(paste(x, collapse = ''))
  if (length(x) == 2) return(paste(x, collapse = ' and '))
  paste0(paste(x[-length(x)], collapse = ', '), ', and ', x[length(x)])
}
big <- function(x) format(x, big.mark = ',')

#' Note shown in place of results that await a rerun.
pending <- function(id, extra = NULL) {
  sp <- SPEC[SPEC$id == id, ]
  script <- sprintf(
    'analysis/scripts/nof1-design-sensitivity/%02d-%s.R', id, sp$name)
  paste0(
    '> **RESULTS-PENDING (S', id, ').** The stored artifact is absent ',
    'or predates the corrected drivers, so no number is reported for ',
    'this sweep. Rerun `Rscript ', script, '` (', big(sp$reps),
    ' replicates per cell, ', sp$cells, ' cells), then ',
    '`Rscript analysis/scripts/nof1-design-sensitivity/plot-all.R`, ',
    'then re-render. ', extra %||% ''
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Select the rows of a summary that match the given values.
pick <- function(sm, ...) {
  cond <- list(...)
  keep <- rep(TRUE, nrow(sm))
  for (n in names(cond)) {
    v <- cond[[n]]
    keep <- keep & if (is.numeric(v)) {
      abs(sm[[n]] - v) < 1e-6
    } else {
      sm[[n]] == v
    }
  }
  sm[keep, , drop = FALSE]
}

#' Power or rate with MCSE for one cell, as text.
cell_txt <- function(sm, ..., what = 'power') {
  r <- pick(sm, ...)
  if (nrow(r) != 1) return('not run')
  mc <- if (what == 'power') r$mcse_power else r$mcse_coverage
  pm(r[[what]], mc)
}

#' First grid value at which power reaches `target`, else NA.
first_reaching <- function(sm, x, target = 0.8) {
  sm <- sm[order(abs(sm[[x]])), ]
  hit <- sm[sm$power >= target, ]
  if (nrow(hit) == 0) NA else hit[[x]][1]
}

fig_path <- function(file) {
  p <- file.path(FIG_DIR, file)
  art <- NULL
  if (file.exists(p)) p else NULL
}

#' TRUE when the figure exists and is newer than its artifact.
fig_ok <- function(id, suffix = '') {
  sp <- SPEC[SPEC$id == id, ]
  f <- file.path(FIG_DIR, sprintf('sens-%02d%s-%s.pdf', id, suffix,
                                  sp$name))
  a <- art_file(id, sp$name)
  file.exists(f) && file.exists(a) && avail(id) &&
    file.mtime(f) >= file.mtime(a)
}

fig_file <- function(id, suffix = '') {
  sp <- SPEC[SPEC$id == id, ]
  file.path(FIG_DIR, sprintf('sens-%02d%s-%s.pdf', id, suffix, sp$name))
}

#' Wide table of rates (MCSE is stated in the caption).
wide_table <- function(sm, rows, cols, what = 'power') {
  labels <- c(delta_on = 'Effect', halflife_days = 'Half-life (days)',
              decay_form = 'Decay form', tau = 'tau',
              bl_sd = 'Baseline SD', sd_scale = 'Multiplier', N = 'N',
              c.br = 'AR(1) correlation', period_days = 'Period (days)',
              c.bm = 'c_bm', dgp = 'DGP condition',
              config = 'Configuration')
  sm$cell <- sprintf('%.3f', sm[[what]])
  out <- sm |>
    select(all_of(c(rows, cols, 'cell'))) |>
    tidyr::pivot_wider(names_from = all_of(cols), values_from = 'cell')
  out <- as.data.frame(out)
  if (rows %in% names(labels)) names(out)[1] <- labels[[rows]]
  out
}

#' Largest MCSE of the rates in a summary, for table captions.
mcse_note <- function(sm) {
  sprintf('Monte Carlo standard error at most %.3f.',
          max(sm$mcse_power, na.rm = TRUE))
}

#' Calibration table: bias (in MCSE), kappa and coverage.
calib_table <- function(sm, keys) {
  out <- sm |>
    transmute(across(all_of(keys)),
              z = sprintf('%.1f', bias / mcse_bias),
              bias = sprintf('%.3f', bias),
              emp_se = sprintf('%.3f', emp_se),
              model_se = sprintf('%.3f', model_se),
              kappa = sprintf('%.3f', kappa),
              coverage = sprintf('%.3f', coverage),
              nonconv = sprintf('%.3f', nonconv))
  as.data.frame(out)
}

#' Count of cells outside the 2-MCSE band around the nominal rate.
outside_band <- function(sm, nominal = 0.05) {
  ok <- within_nominal(sm$power, sm$n_fit, nominal)
  c(above = sum(!ok & sm$power > nominal),
    below = sum(!ok & sm$power < nominal), total = nrow(sm))
}

#=====================================================================
# Design table (does not need artifacts)
#=====================================================================

DESIGNS <- function() {
  list(
    Hybrid = design_hybrid(),
    Compact2 = design_compact2(),
    Alt = design_nof1_alternating(),
    XO = design_crossover(),
    RCT = design_rct(),
    `RCT-OL` = design_rct_openlabel()
  )
}

design_summary_table <- function() {
  ds <- DESIGNS()
  purrr::imap_dfr(ds, function(td, nm) {
    tibble(
      Design = nm,
      Paths = length(td$trialpaths),
      Visits = length(td$metadata$timepoints),
      Weeks = paste(td$metadata$timepoints, collapse = ', '),
      Expectancy = paste(td$metadata$expectancies, collapse = ', '),
      `Realized estimand` = design_delta_on(td, baseline_resp_param())
    )
  })
}

hybrid_path_table <- function() {
  describe_design(design_hybrid()) |>
    mutate(cell = paste0(ifelse(drug == 'on', 'on', 'off'))) |>
    select(path, week, cell) |>
    tidyr::pivot_wider(names_from = week, values_from = cell) |>
    as.data.frame()
}

#=====================================================================
# Results text, one function per sweep
#=====================================================================

reps_of <- function(id) load_art(id)$run_meta$n_reps

txt_s1 <- function() {
  if (!avail(1)) return(pending(1))
  sm <- smry(1, c('design_name', 'delta_on'))
  null <- pick(sm, delta_on = 0)
  base <- pick(sm, delta_on = BASE_DELTA)
  nm <- unique(sm$design_name)
  t1 <- join_txt(sprintf('%s %s', null$design_name, pm(
    null$power, null$mcse_power)))
  bias <- join_txt(sprintf(
    '%s %s (%.1f MCSE)', null$design_name, f3(null$bias),
    null$bias / null$mcse_bias))
  kap <- join_txt(sprintf('%s %s', null$design_name, f3(null$kappa)))
  reach <- join_txt(vapply(nm, function(d) {
    r <- first_reaching(sm[sm$design_name == d, ], 'delta_on')
    sprintf('%s %s', d, if (is.na(r)) 'not reached' else f2(r))
  }, ''))
  inband <- outside_band(null)
  paste0(
    'We varied the mean drug response on drug over ',
    length(unique(sm$delta_on)), ' values from 0 to ',
    f2(min(sm$delta_on)), ' outcome units, with ', big(reps_of(1)),
    ' replicates per cell. The asymptote of the drug-response curve was set for ',
    'each design so that every design had the same value. ',
    'At the null, Type I error was ', t1, '. ',
    'Of ', inband[['total']], ' designs, ', inband[['above']],
    ' were above and ', inband[['below']], ' below the band of two ',
    'MCSE around 0.05. The mean estimate minus its large-sample ',
    'value at the null was ', bias, '. The ratio of model to ',
    'empirical SE ($\\kappa$) at the null was ', kap, '. ',
    'The smallest grid effect with power of at least 0.80 was ',
    reach, '. At the baseline effect of ', f2(BASE_DELTA),
    ', power was ',
    join_txt(sprintf('%s %s', base$design_name, pm(
      base$power, base$mcse_power))), '.'
  )
}

tab_s1 <- function() {
  sm <- smry(1, c('design_name', 'delta_on'))
  wide_table(sm, 'delta_on', 'design_name')
}

txt_s2 <- function() {
  if (!avail(2)) return(pending(2))
  sm <- smry(2, c('design_name', 'halflife_days', 'analysis_t1half'))
  hmin <- min(sm$halflife_days)
  hmax <- max(sm$halflife_days)
  an <- sort(unique(sm$analysis_t1half))
  parts <- vapply(an, function(a) {
    s <- sm[abs(sm$analysis_t1half - a) < 1e-9, ]
    paste0(
      'With the analysis half-life at ', f2(a), ' week, power at ',
      hmin, ' and ', hmax, ' days was ',
      join_txt(vapply(unique(s$design_name), function(d) {
        sprintf('%s %s and %s', d,
                f3(pick(s, design_name = d, halflife_days = hmin)$power),
                f3(pick(s, design_name = d,
                        halflife_days = hmax)$power))
      }, '')), '.')
  }, '')
  paste0(
    'We varied the data-generating half-life over ',
    length(unique(sm$halflife_days)), ' values from ', hmin, ' to ',
    hmax, ' days. Each dataset was analyzed with the binary ',
    'on/off indicator (analysis half-life 0) and with the decaying ',
    'indicator at 0.2 week. ', paste(parts, collapse = ' '),
    ' Each cell has ', big(reps_of(2)), ' replicates.'
  )
}

tab_s2 <- function() {
  sm <- smry(2, c('design_name', 'halflife_days', 'analysis_t1half'))
  sm$design_analysis <- paste0(sm$design_name, ', ',
                               f1(sm$analysis_t1half), ' wk')
  wide_table(sm, 'halflife_days', 'design_analysis')
}

txt_s3 <- function() {
  if (!avail(3)) return(pending(3))
  r <- load_art(3)$results
  r$decay_form <- ifelse(r$carryover_form == 'weibull',
                         sprintf('weibull_k=%g', r$weibull_shape),
                         r$carryover_form)
  sm <- summarize_cells(r, c('design_name', 'decay_form',
                             'analysis_t1half'))
  parts <- vapply(sort(unique(sm$analysis_t1half)), function(a) {
    s <- sm[abs(sm$analysis_t1half - a) < 1e-9, ]
    paste0(
      'With analysis half-life ', f2(a), ' week, power ranged across ',
      'decay forms from ', f3(min(s$power)), ' to ', f3(max(s$power)),
      ' (',
      join_txt(vapply(unique(s$design_name), function(d) {
        z <- s[s$design_name == d, ]
        sprintf('%s %s to %s', d, f3(min(z$power)), f3(max(z$power)))
      }, '')), ').')
  }, '')
  paste0(
    'We varied the form of the data-generating carryover over ',
    length(unique(sm$decay_form)), ' forms at a half-life of ',
    '3 days. The analysis always built the decaying indicator from ',
    'an exponential, at three half-lives. ', paste(parts, collapse = ' '),
    ' Each cell has ', big(reps_of(3)), ' replicates.'
  )
}

tab_s3 <- function() {
  r <- load_art(3)$results
  r$decay_form <- ifelse(r$carryover_form == 'weibull',
                         sprintf('weibull_k=%g', r$weibull_shape),
                         r$carryover_form)
  sm <- summarize_cells(r, c('design_name', 'decay_form',
                             'analysis_t1half'))
  sm$design_analysis <- paste0(sm$design_name, ', ',
                               f2(sm$analysis_t1half), ' wk')
  wide_table(sm, 'decay_form', 'design_analysis')
}

txt_s4 <- function() {
  if (!avail(4)) return(pending(4))
  sm <- smry(4, c('axis', 'design_name', 'tau', 'bl_sd'))
  a <- sm[sm$axis == 'tau', ]
  b <- sm[sm$axis == 'bl_sd', ]
  rng <- function(s, x) {
    join_txt(vapply(unique(s$design_name), function(d) {
      z <- s[s$design_name == d, ]
      z <- z[order(z[[x]]), ]
      sprintf('%s %s at %s and %s at %s', d, f3(z$power[1]), z[[x]][1],
              f3(z$power[nrow(z)]), z[[x]][nrow(z)])
    }, ''))
  }
  paste0(
    'We varied the SD of participant-level heterogeneity of the drug ',
    'effect ($\\tau$) over ', length(unique(a$tau)), ' values and, ',
    'separately, the SD of the baseline score over ',
    length(unique(b$bl_sd)), ' values. Across $\\tau$, power was ',
    rng(a, 'tau'), '. Across the baseline SD, power was ',
    rng(b, 'bl_sd'), '. Each cell has ', big(reps_of(4)),
    ' replicates.'
  )
}

tab_s4 <- function() {
  sm <- smry(4, c('axis', 'design_name', 'tau', 'bl_sd'))
  a <- wide_table(sm[sm$axis == 'tau', ], 'tau', 'design_name')
  b <- wide_table(sm[sm$axis == 'bl_sd', ], 'bl_sd', 'design_name')
  list(tau = a, bl_sd = b)
}

txt_s5 <- function() {
  if (!avail(5)) return(pending(5))
  r <- load_art(5)$results
  r$sd_scale <- round(r$br_sd, 4)
  sm <- summarize_cells(r, c('design_name', 'sd_scale'))
  paste0(
    'We multiplied the SDs of all three components by a common ',
    'factor from ', min(sm$sd_scale), ' to ', max(sm$sd_scale),
    ' (', length(unique(sm$sd_scale)), ' values). Power at the ',
    'smallest and largest factor was ',
    join_txt(vapply(unique(sm$design_name), function(d) {
      z <- sm[sm$design_name == d, ]
      z <- z[order(z$sd_scale), ]
      sprintf('%s %s and %s', d, f3(z$power[1]), f3(z$power[nrow(z)]))
    }, '')), '. Each cell has ', big(reps_of(5)), ' replicates.'
  )
}

tab_s5 <- function() {
  r <- load_art(5)$results
  r$sd_scale <- round(r$br_sd, 4)
  wide_table(summarize_cells(r, c('design_name', 'sd_scale')),
             'sd_scale', 'design_name')
}

txt_s6 <- function() {
  if (!avail(6)) return(pending(6))
  r <- load_art(6)$results
  r$periods <- r$n_switch + 1
  sm <- summarize_cells(r, c('design_name', 'periods'))
  alt <- sm[sm$design_name != 'RCT', ]
  alt <- alt[order(alt$periods), ]
  rct <- sm[sm$design_name == 'RCT', ]
  paste0(
    'We varied the number of equal periods of the alternating design ',
    'over ', paste(alt$periods, collapse = ', '), ', with 16 weekly ',
    'visits throughout. Power of the alternating design was ',
    join_txt(sprintf('%s with %d periods', f3(alt$power),
                     as.integer(alt$periods))),
    '. Power of the parallel group with the same 16 visits was ',
    pm(rct$power, rct$mcse_power), '. The realized estimand was ',
    f3(mean(sm$delta_on_cell)), ' for all designs. Each cell has ',
    big(reps_of(6)), ' replicates.'
  )
}

tab_s6 <- function() {
  r <- load_art(6)$results
  r$periods <- r$n_switch + 1
  sm <- summarize_cells(r, c('design_name', 'periods'))
  sm$power_txt <- pm(sm$power, sm$mcse_power)
  as.data.frame(sm[, c('design_name', 'periods', 'power_txt', 'kappa')])
}

txt_s7 <- function() {
  if (!avail(7)) return(pending(7))
  sm <- smry(7, c('design_name', 'N'))
  reach <- join_txt(vapply(unique(sm$design_name), function(d) {
    z <- sm[sm$design_name == d, ]
    z <- z[order(z$N), ]
    hit <- z$N[z$power >= 0.8]
    sprintf('%s %s', d, if (length(hit)) hit[1] else 'not reached')
  }, ''))
  paste0(
    'We varied the total sample size over ',
    paste(sort(unique(sm$N)), collapse = ', '),
    '. The smallest N in the grid with power of at least 0.80 was ',
    reach, '. Each cell has ', big(reps_of(7)), ' replicates.'
  )
}

tab_s7 <- function() {
  wide_table(smry(7, c('design_name', 'N')), 'N', 'design_name')
}

txt_s8 <- function() {
  if (!avail(8)) return(pending(8))
  art <- load_art(8)
  sm <- summarize_cells(art$results, c('design_name', 'c.br'))
  inf <- art$infeasible
  paste0(
    'We varied the AR(1) correlation per week over ',
    paste(c(0, 0.3, 0.5, 0.7, 0.9), collapse = ', '), '. ',
    nrow(inf), ' of ', SPEC$cells[8],
    ' cells had a covariance matrix that was not positive definite ',
    'and were dropped, not repaired',
    if (nrow(inf) > 0) paste0(
      ' (', join_txt(sprintf('%s at %s, smallest eigenvalue %s',
                             inf$design_name,
                             c(0, 0.3, 0.5, 0.7, 0.9)[inf$modelparamset],
                             f3(inf$min_eig))), ')') else '',
    '. Power at the smallest and largest run correlation was ',
    join_txt(vapply(unique(sm$design_name), function(d) {
      z <- sm[sm$design_name == d, ]
      z <- z[order(z$c.br), ]
      sprintf('%s %s at %s and %s at %s', d, f3(z$power[1]), z$c.br[1],
              f3(z$power[nrow(z)]), z$c.br[nrow(z)])
    }, '')), '. Each cell has ', big(reps_of(8)), ' replicates.'
  )
}

tab_s8 <- function() {
  sm <- summarize_cells(load_art(8)$results, c('design_name', 'c.br'))
  wide_table(sm, 'c.br', 'design_name')
}

txt_s9 <- function() {
  if (!avail(9)) return(pending(9))
  r <- load_art(9)$results
  r$period_days <- as.numeric(sub('^Alt-(\\d+)d$', '\\1', r$design_name))
  sm <- summarize_cells(r, c('period_days', 'halflife_days'))
  paste0(
    'We crossed ', length(unique(sm$period_days)), ' period lengths (',
    paste(sort(unique(sm$period_days)), collapse = ', '),
    ' days) with ', length(unique(sm$halflife_days)), ' half-lives (',
    paste(sort(unique(sm$halflife_days)), collapse = ', '),
    ' days). The stored half-life takes ',
    length(unique(round(r$carryover_t1half, 8))), ' distinct values ',
    'in ', nrow(sm), ' cells. Power ranged from ', f3(min(sm$power)),
    ' to ', f3(max(sm$power)), '. At the longest period, power at the ',
    'shortest and longest half-life was ',
    f3(pick(sm, period_days = max(sm$period_days),
            halflife_days = min(sm$halflife_days))$power), ' and ',
    f3(pick(sm, period_days = max(sm$period_days),
            halflife_days = max(sm$halflife_days))$power),
    '. Each cell has ', big(reps_of(9)), ' replicates.'
  )
}

tab_s9 <- function() {
  r <- load_art(9)$results
  r$period_days <- as.numeric(sub('^Alt-(\\d+)d$', '\\1', r$design_name))
  wide_table(summarize_cells(r, c('period_days', 'halflife_days')),
             'period_days', 'halflife_days')
}

txt_s10 <- function() {
  if (!avail(10)) return(pending(10))
  art <- load_art(10)
  sm <- summarize_cells(art$results, c('design_name', 'dgp_architecture',
                                       'c.bm', 'halflife_days'))
  inf <- art$infeasible
  inf_tab <- if (nrow(inf) > 0) {
    cb <- c(0, 0.2, 0.4, 0.6)
    paste0(nrow(inf), ' of ', SPEC$cells[10], ' cells were not ',
           'positive definite and were dropped: ',
           join_txt(unique(sprintf(
             'c_bm of %s under the %s encoding',
             rep(cb, each = 2)[inf$modelparamset], inf$dgp_architecture
           ))), '.')
  } else 'No cell was dropped.'
  paste0(
    'The sweep stored the main-effect coefficient. We varied the ',
    'biomarker-response correlation $c_{bm}$ over 0, 0.2, 0.4 and ',
    '0.6, the carryover half-life over 0 and 3 days, and the ',
    'encoding of the interaction (the correlation encoding `mvn` and ',
    'the mean-shift encoding `mean_moderation`). ', inf_tab,
    ' Across the cells that ran, power ranged from ', f3(min(sm$power)),
    ' to ', f3(max(sm$power)), ' (',
    join_txt(vapply(unique(sm$design_name), function(d) {
      z <- sm[sm$design_name == d, ]
      sprintf('%s %s to %s', d, f3(min(z$power)), f3(max(z$power)))
    }, '')), '). Each cell has ', big(reps_of(10)), ' replicates.'
  )
}

tab_s10 <- function() {
  sm <- summarize_cells(load_art(10)$results,
                        c('design_name', 'dgp_architecture', 'c.bm',
                          'halflife_days'))
  sm$setting <- paste0(sm$design_name, ', ', sm$dgp_architecture,
                       ', ', sm$halflife_days, ' d')
  wide_table(sm, 'c.bm', 'setting')
}

txt_s11 <- function() {
  if (!avail(11)) return(pending(11))
  sm <- smry(11, c('design_name', 'delta_label', 'carryover', 'ar1',
                   'analysis_t1half'))
  nl <- sm[sm$delta_label == 'null', ]
  parts <- vapply(unique(nl$design_name), function(d) {
    z <- nl[nl$design_name == d, ]
    ob <- outside_band(z)
    sprintf('%s %s to %s (%d of %d cells outside the band)', d,
            f3(min(z$power)), f3(max(z$power)), ob[['above']] +
              ob[['below']], ob[['total']])
  }, '')
  al <- sm[sm$delta_label == 'alt', ]
  paste0(
    'We crossed data-generating carryover (off, 3 days) with ',
    'data-generating AR(1) (off, 0.7), at the null and at the ',
    'baseline effect, and analyzed each dataset with three analysis ',
    'half-lives. Type I error ranged over cells and analyses from ',
    join_txt(parts), '. Under the alternative, power ranged from ',
    f3(min(al$power)), ' to ', f3(max(al$power)), '. Each cell has ',
    big(reps_of(11)), ' replicates.'
  )
}

tab_s11 <- function() {
  sm <- smry(11, c('design_name', 'delta_label', 'carryover', 'ar1',
                   'analysis_t1half'))
  sm$setting <- paste0(sm$design_name, ', ', f2(sm$analysis_t1half),
                       ' wk')
  sm$dgp <- paste0(sm$delta_label, ': carryover ', sm$carryover,
                   ', AR(1) ', sm$ar1)
  wide_table(sm, 'dgp', 'setting')
}

txt_s12 <- function() {
  if (!avail(12)) return(pending(12))
  r <- load_art(12)$results
  sm <- summarize_cells(r, c('design_name', 'config', 'analysis_t1half'))
  pooled <- r |>
    group_by(design_name, analysis_t1half) |>
    summarise(rate = mean(p < 0.05, na.rm = TRUE),
              n = sum(!is.na(p)), .groups = 'drop') |>
    mutate(mcse = sqrt(rate * (1 - rate) / n),
           z = (rate - 0.05) / sqrt(0.05 * 0.95 / n))
  parts <- vapply(unique(pooled$design_name), function(d) {
    z <- pooled[pooled$design_name == d, ]
    paste0(d, ' ', join_txt(sprintf(
      '%s (MCSE %s, %.1f MCSE from 0.05) at analysis half-life %s',
      f3(z$rate), f3(z$mcse), z$z, f2(z$analysis_t1half))))
  }, '')
  paste0(
    'We simulated a null drug effect in ', length(unique(sm$config)),
    ' configurations of the baseline SD, the drug-response SD and ',
    'the half-life, with ', big(reps_of(12)), ' replicates per cell. ',
    'Pooled over configurations, Type I error was: ',
    paste(parts, collapse = '; '), '.'
  )
}

tab_s12 <- function() {
  sm <- smry(12, c('design_name', 'config', 'analysis_t1half'))
  sm$setting <- paste0(sm$design_name, ', ', f2(sm$analysis_t1half),
                       ' wk')
  wide_table(sm, 'config', 'setting')
}

txt_s12b <- function() {
  if (is.null(load_art(15, 'xo-corcar1'))) {
    return(paste0(
      '> **RESULTS-PENDING (S12b).** Rerun `Rscript analysis/scripts/',
      'nof1-design-sensitivity/15-xo-corcar1.R` (2,000 replicates, ',
      '2 designs, 2 analyses).'))
  }
  r <- load_art(15, 'xo-corcar1')$results
  r$corCAR1 <- c(TRUE, FALSE)[r$analysis_row]
  sm <- summarize_cells(r, c('design_name', 'corCAR1'))
  paste0(
    'Refitting the null data without the serial-correlation ',
    'structure gave Type I error of ',
    join_txt(sprintf('%s for %s (%s corCAR1)', f3(sm$power),
                     sm$design_name,
                     ifelse(sm$corCAR1, 'with', 'without'))),
    '. Each cell has ', big(load_art(15, 'xo-corcar1')$run_meta$n_reps),
    ' replicates.'
  )
}

#=====================================================================
# Parallel-group check (script 14) and calibration summaries
#=====================================================================

check_file <- function() {
  file.path(DATA_DIR, 'sensitivity-14-rct-null-check.rds')
}

load_check <- function() {
  p <- check_file()
  if (!file.exists(p)) return(NULL)
  obj <- readRDS(p)
  if (!identical(obj$run_meta$schema, 2L)) return(NULL)
  r <- obj$results
  r$analysis <- c('legacy', 'ancova')[r$analysis_row]
  r$rho <- r$c.br
  list(results = r, meta = obj$run_meta)
}

rct_check_table <- function() {
  ck <- load_check()
  if (is.null(ck)) return(NULL)
  tab <- summarize_cells(ck$results,
                         c('design_name', 'analysis', 'delta_on', 'rho'))
  tab$z <- tab$bias_estimand / tab$mcse_bias
  data.frame(
    Design = tab$design_name,
    Analysis = tab$analysis,
    Effect = f2(ifelse(abs(tab$delta_on) < 1e-9, 0, tab$delta_on)),
    rho = tab$rho,
    `Mean estimate` = f3(tab$mean_beta),
    `Bias vs estimand` = f3(tab$bias_estimand),
    `Bias/MCSE` = f1(tab$z),
    `Empirical SE` = f3(tab$emp_se),
    `Model SE` = f3(tab$model_se),
    kappa = f3(tab$kappa),
    `Rejection rate` = f3(tab$power),
    check.names = FALSE
  )
}

rct_check_text <- function() {
  ck <- load_check()
  if (is.null(ck)) {
    return(paste0(
      '> **RESULTS-PENDING (parallel-group check).** Run `Rscript ',
      'analysis/scripts/nof1-design-sensitivity/14-rct-null-check.R`.'))
  }
  tab <- summarize_cells(ck$results,
                         c('design_name', 'analysis', 'delta_on', 'rho'))
  nl <- tab[abs(tab$delta_on) < 1e-9, ]
  leg <- nl[nl$analysis == 'legacy', ]
  anc <- nl[nl$analysis == 'ancova', ]
  zz <- function(d) abs(d$bias_estimand / d$mcse_bias)
  paste0(
    'Under the null, the mean estimate of the mixed-model analysis of ',
    'the parallel group ranged from ', f3(min(leg$mean_beta)), ' to ',
    f3(max(leg$mean_beta)), ', and its Type I error from ',
    f3(min(leg$power)), ' to ', f3(max(leg$power)),
    '. The largest $|$mean$|$/MCSE was ', f1(max(zz(leg))),
    ' and $\\kappa$ ranged from ', f3(min(leg$kappa)), ' to ',
    f3(max(leg$kappa)), '. The ANCOVA arm contrast had mean ',
    'estimates from ', f3(min(anc$mean_beta)), ' to ',
    f3(max(anc$mean_beta)), ' (largest $|$mean$|$/MCSE ',
    f1(max(zz(anc))), '), $\\kappa$ from ', f3(min(anc$kappa)),
    ' to ', f3(max(anc$kappa)), ', and Type I error from ',
    f3(min(anc$power)), ' to ', f3(max(anc$power)), '. Each cell has ',
    big(ck$meta$n_reps), ' replicates.'
  )
}


#=====================================================================
# Abstract and synthesis text
#=====================================================================

abs_results <- function() {
  if (!(avail(1) && avail(12))) {
    return(paste0(
      'The results of the corrected sweeps are pending. They are ',
      'generated from the stored artifacts when the sweeps have been ',
      'rerun.'))
  }
  sm1 <- smry(1, c('design_name', 'delta_on'))
  base <- pick(sm1, delta_on = BASE_DELTA)
  r <- load_art(12)$results
  pooled <- r |>
    group_by(design_name) |>
    summarise(rate = mean(p < 0.05, na.rm = TRUE), .groups = 'drop')
  paste0(
    'At a mean drug response of ', f2(BASE_DELTA), ' outcome units ',
    'on drug and N = 35, power was ',
    join_txt(sprintf('%s %s', base$design_name, f2(base$power))),
    '. In the null calibration, pooled Type I error was ',
    join_txt(sprintf('%s %s', pooled$design_name, f3(pooled$rate))),
    '.')
}

synth_text <- function() {
  if (!avail(1)) return(pending(1))
  sm1 <- smry(1, c('design_name', 'delta_on'))
  base <- pick(sm1, delta_on = BASE_DELTA)
  base <- base[order(-base$power), ]
  out <- paste0(
    'At the baseline effect, designs ordered by power were ',
    paste(sprintf('%s (%s)', base$design_name, f2(base$power)),
          collapse = ', '), '. ')
  if (avail(7)) {
    sm7 <- smry(7, c('design_name', 'N'))
    out <- paste0(out, 'The smallest N in the grid with power of at ',
      'least 0.80 was ',
      join_txt(vapply(unique(sm7$design_name), function(d) {
        z <- sm7[sm7$design_name == d, ]
        z <- z[order(z$N), ]
        hit <- z$N[z$power >= 0.8]
        sprintf('%s %s', d, if (length(hit)) hit[1] else 'not reached')
      }, '')), '. ')
  }
  out
}

#=====================================================================
# Feasibility of the S8 and S10 grids (script 17, no simulation)
#=====================================================================

load_feas <- function() {
  p <- file.path(DATA_DIR, 'feasibility-s8-s10.rds')
  if (!file.exists(p)) return(NULL)
  obj <- readRDS(p)
  if (!identical(obj$run_meta$schema, 2L)) return(NULL)
  obj$table
}

feas_text <- function() {
  ft <- load_feas()
  if (is.null(ft)) {
    return(paste0(
      '> **RESULTS-PENDING (feasibility).** Run `Rscript analysis/',
      'scripts/nof1-design-sensitivity/17-feasibility.R`. It needs ',
      'no simulation.'))
  }
  bad8 <- ft[ft$sweep == 'S8' & !ft$feasible, ]
  bad10 <- ft[ft$sweep == 'S10' & !ft$feasible, ]
  paste0(
    'In S8, ', nrow(bad8), ' of ', sum(ft$sweep == 'S8'),
    ' design and correlation cells are infeasible: ',
    join_txt(sprintf('%s at $\\rho$ = %s (smallest eigenvalue %s)',
                     bad8$design, bad8$rho, f3(bad8$min_eig))),
    '. In S10 under the correlation encoding, ', nrow(bad10), ' of ',
    sum(ft$sweep == 'S10'), ' cells are infeasible: ',
    join_txt(unique(sprintf('$c_{bm}$ = %s', bad10$c_bm))),
    ' for ', join_txt(unique(bad10$design)), '. The smallest ',
    'eigenvalue of the feasible cells is ',
    f3(min(ft$min_eig[ft$feasible])), '. The same check runs inside ',
    'every sweep, and the artifacts list the dropped cells.')
}

feas_tables <- function() {
  ft <- load_feas()
  s8 <- ft[ft$sweep == 'S8', ]
  a <- tidyr::pivot_wider(
    s8[, c('design', 'rho', 'min_eig')],
    names_from = 'design', values_from = 'min_eig')
  s10 <- ft[ft$sweep == 'S10', ]
  s10$setting <- paste0(s10$design, ', ', s10$halflife_days, ' d')
  b <- tidyr::pivot_wider(
    s10[, c('c_bm', 'setting', 'min_eig')],
    names_from = 'setting', values_from = 'min_eig')
  list(s8 = as.data.frame(a), s10 = as.data.frame(b))
}
