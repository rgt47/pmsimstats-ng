#' 04-study-C.R
#'
#' Paper 06 Study C: subadditivity sensitivity (mechanistic vs
#' emergent).
#'
#' Two DGP families per 00-ademp-pre-registration.md:
#'   1. Mechanistic subadditivity: gamma in {-0.05, 0, 0.02, 0.05,
#'      0.10}; the DGP adds gamma * BR * PB to the additive sum.
#'   2. Latent-class with class-correlated PB: reuses paper 03's
#'      mixture infrastructure.
#'
#' Falsifiability test: under the latent-class DGP with no true
#' gamma, the saturated mean-structure analysis should still recover
#' a non-zero gamma-hat. This demonstrates that gamma-hat alone is
#' not evidence of mechanistic subadditivity.
#'
#' NOTE. The mechanistic-subadditivity DGP requires modifying
#' generateData to add gamma * BR * PB after the MVN draw. This is
#' implemented inline below as a post-hoc shift, mirroring the
#' Architecture A 'mean_moderation' pattern from R/generateData.R.
#'
#' The latent-class arm requires paper 03's mixture infrastructure
#' to be exposed as a callable from this script. The relevant DGP
#' lives at analysis/scripts/latent-class-mixture-application/
#' 02-heterogeneous-slopes-dgp.R and adjacent files. Until that arm
#' is wired in, the latent-class cells are stubbed and the
#' mechanistic arm runs to completion.
#'
#' Seed: study_seed = 42 + 6*100 + 3 = 644.

source(here::here('analysis', 'scripts',
                  'component-decomposition', '00-common.R'))

set.seed(644)
STUDY_ID   <- 'C'
STUDY_SEED <- 644L

gamma_levels <- c(-0.05, 0, 0.02, 0.05, 0.10)
N_levels     <- c(35L, 70L, 100L, 150L)

mechanistic_cells <- CJ(gamma = gamma_levels, N = N_levels,
                        sorted = FALSE)
mechanistic_cells[, dgp_family := 'mechanistic']
mechanistic_cells[, c_bm := 0.45]
mechanistic_cells[, m_PB := 6]   # fixed at the prazosin baseline
mechanistic_cells[, m_TV := 0]
mechanistic_cells[, n_reps := 1000L]

## Latent-class arm: 4 N x 1 gamma_true=0 reference cell. The
## class-correlation parameter is varied via the latent-class DGP
## driver (TODO); for now, mark as stub.
latent_cells <- data.table(gamma = 0, N = N_levels,
                           dgp_family = 'latent_class',
                           c_bm = 0.45, m_PB = 6, m_TV = 0,
                           n_reps = 1000L)

cells <- rbindlist(list(mechanistic_cells, latent_cells),
                   use.names = TRUE, fill = TRUE)
cells[, cell_id := .I]
cells[, analyses := list(list(c('phase_augmented')))]

cat(sprintf('Study C: %d cells (%d mechanistic + %d latent_class).\n',
            nrow(cells), nrow(mechanistic_cells),
            nrow(latent_cells)))

#=====================================================================
# Mechanistic-subadditivity DGP shim
#=====================================================================
#
# Wraps generateData with a post-hoc additive shift gamma * BR * PB.
# Returns a data.table with the same columns as generateData, with
# the Sx_* columns adjusted.

generate_data_with_subadditivity <- function(mp, rp, bp, td_path,
                                              gamma) {
  dat <- generateData(mp, rp, bp, td_path,
                      empirical = FALSE,
                      makePositiveDefinite = TRUE,
                      dgp_architecture = 'mvn')
  if (gamma == 0) return(dat)
  ## Find Sx columns and the matching BR / PB columns. In pmsimstats
  ## the data table holds Sx_<tp>, BR_<tp>, PB_<tp> for each
  ## timepoint label tp. We shift Sx by -gamma * BR * PB at each tp.
  sx_cols <- grep('^Sx_', names(dat), value = TRUE)
  for (sx in sx_cols) {
    tp <- sub('^Sx_', '', sx)
    br_col <- paste0('BR_', tp)
    pb_col <- paste0('PB_', tp)
    if (br_col %in% names(dat) && pb_col %in% names(dat)) {
      dat[, (sx) := get(sx) - gamma * get(br_col) * get(pb_col)]
    }
  }
  dat
}

#=====================================================================
# Per-replicate (overrides one_rep_06 for Study C)
#=====================================================================

one_rep_C <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- design_hybrid_simple()
  mp <- make_model_params(N = cell$N, c_bm = cell$c_bm,
                          t1half = 1)
  rp <- make_resp_params(m_PB = cell$m_PB, m_TV = cell$m_TV)
  bp <- make_bl_params()

  dat <- if (cell$dgp_family == 'mechanistic') {
    tryCatch(
      generate_data_with_subadditivity(mp, rp, bp,
                                       td$trialpaths[[1]],
                                       gamma = cell$gamma),
      error = function(e) NULL)
  } else {
    ## TODO: latent-class DGP via paper 03 infrastructure.
    NULL
  }
  if (is.null(dat)) {
    return(data.table(analysis = 'phase_augmented',
                      rep_idx = rep_idx,
                      beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = if (cell$dgp_family ==
                                            'latent_class')
                        'lc_stub' else 'gen_failed'))
  }
  dat[, path := 1]
  fit <- phase_augmented_fit(td, dat, cell$N)
  cbind(analysis = 'phase_augmented', rep_idx = rep_idx, fit)
}

run_cell_C <- function(cell, n_reps, study_seed, cell_id,
                       n_cores = max(1L,
                                     parallel::detectCores() - 1L)) {
  results <- parallel::mclapply(
    seq_len(n_reps),
    function(r) one_rep_C(r, cell, study_seed, cell_id),
    mc.cores = n_cores)
  good <- !vapply(results, is.null, logical(1))
  out <- rbindlist(results[good], fill = TRUE)
  out[, `:=`(gamma = cell$gamma, dgp_family = cell$dgp_family,
             N = cell$N, c_bm = cell$c_bm,
             m_PB = cell$m_PB, m_TV = cell$m_TV,
             cell_id = cell_id)]
  out
}

#=====================================================================
# Execute
#=====================================================================

t_study <- Sys.time()
cell_summaries <- vector('list', nrow(cells))

for (i in seq_len(nrow(cells))) {
  cell <- list(gamma = cells$gamma[i],
               dgp_family = cells$dgp_family[i],
               N = cells$N[i], c_bm = cells$c_bm[i],
               m_PB = cells$m_PB[i], m_TV = cells$m_TV[i])
  reps_dt <- run_cell_C(cell, n_reps = cells$n_reps[i],
                        study_seed = STUDY_SEED,
                        cell_id = cells$cell_id[i])
  save_cell_06(reps_dt, STUDY_ID, cells$cell_id[i])

  tb <- TRUE_BETA(c_bm = cells$c_bm[i])
  by_an <- reps_dt[, summarise_cell_06(.SD, tb), by = analysis]
  by_an[, `:=`(cell_id = cells$cell_id[i],
               gamma = cells$gamma[i],
               dgp_family = cells$dgp_family[i],
               N = cells$N[i], c_bm = cells$c_bm[i],
               true_beta = tb)]
  cell_summaries[[i]] <- by_an
  cat(sprintf('Study C cell %d/%d done.\n', i, nrow(cells)))
}

study_summary <- rbindlist(cell_summaries, use.names = TRUE,
                           fill = TRUE)
save_summary_06(study_summary, STUDY_ID)

cat(sprintf('Study C complete: %.1f min.\n',
            as.numeric(difftime(Sys.time(), t_study,
                                units = 'mins'))))
