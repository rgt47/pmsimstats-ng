#' 04-study4.R
#'
#' Paper 09 Study 4: sensitivity to dropout-pattern misspecification.
#'
#' Per 00-ademp-pre-registration.md: 4 designs x 5 patterns x 3
#' analyses x N = 35 = 60 cells. 1000 reps per cell.
#'
#' Three analysis variants:
#'   1. complete_case: drop participants on first missing.
#'   2. mar_mmrm:      Study 1 analysis (lme with corCAR1, ML).
#'   3. mi_mnar:       multiple imputation under correct MNAR
#'                     mechanism.
#'
#' NOTE on MI-MNAR. The pre-registered third analysis requires a
#' multiple-imputation procedure that knows the true dropout
#' mechanism. A defensible implementation uses pattern-mixture
#' modelling with the dropout indicator as a covariate (see
#' Carpenter & Kenward 2013). Implementing this in `mice` requires a
#' custom imputation method that conditions on the dropout-pattern
#' indicator. Until that is implemented, the third analysis is a
#' STUB: the same lme fit on the dropped data, marked as 'mi_stub'
#' in the analysis_name column. Replace the stub with a real MI fit
#' before drawing inferences from this study.
#'
#' Seed: study_seed = 42 + 9 * 100 + 4 = 946.

source(here::here('analysis', 'scripts',
                  'informative-dropout-by-design', '00-common.R'))

set.seed(946)
STUDY_ID   <- '4'
STUDY_SEED <- 946L

design_levels   <- c('OL', 'OLBDC', 'CO', 'Hybrid')
dropout_levels  <- names(DROPOUT_PATTERNS)
analysis_levels <- c('complete_case', 'mar_mmrm', 'mi_stub')
N_FIXED         <- 35L
c_bm_levels     <- c(0.3, 0.6)

cells <- CJ(design_name   = design_levels,
            dropout_name  = dropout_levels,
            analysis_name = analysis_levels,
            c_bm          = c_bm_levels,
            sorted = FALSE)
cells[, N       := N_FIXED]
cells[, n_reps  := 1000L]
cells[, cell_id := .I]

cat(sprintf('Study 4: %d cells; total target reps = %d.\n',
            nrow(cells), sum(cells$n_reps)))

#=====================================================================
# Analysis dispatch
#=====================================================================

analyze_S4 <- function(td, dat, analysis_name) {
  if (analysis_name == 'complete_case') {
    ## Drop participants with any missing Sx_ value.
    sx_cols <- grep('^Sx_', names(dat), value = TRUE)
    keep <- complete.cases(dat[, ..sx_cols])
    dat_cc <- dat[keep]
    if (nrow(dat_cc) < 4L) {
      return(data.table(beta = NA_real_, betaSE = NA_real_,
                        p = NA_real_, issingular = NA,
                        warning = 'too_few_complete'))
    }
    return(lme_fit_09(td, dat_cc))
  } else if (analysis_name == 'mar_mmrm') {
    return(lme_fit_09(td, dat))
  } else if (analysis_name == 'mi_stub') {
    ## TODO: replace with MI-MNAR via pattern-mixture mice routine.
    ## For now mirror mar_mmrm so the pipeline runs end-to-end.
    res <- lme_fit_09(td, dat)
    res[, warning := 'mi_stub_returns_mar_mmrm']
    return(res)
  } else {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, issingular = NA,
                      warning = 'unknown_analysis'))
  }
}

one_rep_S4 <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- ALL_DESIGNS()[[cell$design_name]]
  paths <- td$trialpaths
  mp <- make_model_params_09(N = cell$N, c_bm = cell$c_bm)
  rp <- make_resp_params_09()
  bp <- make_bl_params_09()
  drop <- DROPOUT_PATTERNS[[cell$dropout_name]]

  dat_list <- vector('list', length(paths))
  drop_fracs <- numeric(length(paths))
  for (g in seq_along(paths)) {
    di <- tryCatch(
      generateData(mp, rp, bp, paths[[g]],
                   empirical = FALSE, makePositiveDefinite = TRUE,
                   dgp_architecture = 'mvn'),
      error = function(e) NULL)
    if (is.null(di)) next
    di[, path := g]
    di_d <- apply_hazard_dropout(di, paths[[g]],
                                  beta0 = drop$beta0,
                                  beta1 = drop$beta1)
    drop_fracs[g] <- attr(di_d, 'dropout_fraction') %||% 0
    dat_list[[g]] <- di_d
  }
  good <- !vapply(dat_list, is.null, logical(1))
  if (!any(good)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, dropout_fraction = NA_real_,
                      rep_idx = rep_idx))
  }
  dat <- rbindlist(dat_list[good], fill = TRUE)
  res <- analyze_S4(td, dat, cell$analysis_name)
  res[, dropout_fraction := mean(drop_fracs[good])]
  res[, rep_idx := rep_idx]
  res
}

run_cell_S4 <- function(cell, n_reps, study_seed, cell_id,
                        n_cores = max(1L,
                                      parallel::detectCores() - 1L),
                        chunk_size = 50L) {
  results <- vector('list', 0L)
  done <- 0L
  while (done < n_reps) {
    take <- min(chunk_size, n_reps - done)
    chunk <- parallel::mclapply(
      seq.int(done + 1L, done + take),
      function(r) one_rep_S4(r, cell, study_seed, cell_id),
      mc.cores = n_cores)
    good <- !vapply(chunk, function(x)
      inherits(x, 'try-error') || is.null(x), logical(1))
    results[[length(results) + 1L]] <- rbindlist(chunk[good],
                                                 fill = TRUE)
    done <- done + take
  }
  out <- rbindlist(results, fill = TRUE)
  out[, `:=`(design_name = cell$design_name,
             dropout_name = cell$dropout_name,
             analysis_name = cell$analysis_name,
             N = cell$N, c_bm = cell$c_bm,
             cell_id = cell_id)]
  out
}

t_study <- Sys.time()
cell_summaries <- vector('list', nrow(cells))

for (i in seq_len(nrow(cells))) {
  cell <- list(design_name   = cells$design_name[i],
               dropout_name  = cells$dropout_name[i],
               analysis_name = cells$analysis_name[i],
               N             = cells$N[i],
               c_bm          = cells$c_bm[i])
  reps_dt <- run_cell_S4(cell,
                         n_reps    = cells$n_reps[i],
                         study_seed = STUDY_SEED,
                         cell_id   = cells$cell_id[i])
  save_cell_09(reps_dt, STUDY_ID, cells$cell_id[i])

  tb <- TRUE_BETA_09(c_bm = cells$c_bm[i])
  s  <- summarise_cell_09(reps_dt, tb)
  s[, `:=`(cell_id = cells$cell_id[i],
           design_name = cells$design_name[i],
           dropout_name = cells$dropout_name[i],
           analysis_name = cells$analysis_name[i],
           N = cells$N[i], c_bm = cells$c_bm[i],
           true_beta = tb)]
  cell_summaries[[i]] <- s
  cat(sprintf('Study 4 cell %d/%d done.\n', i, nrow(cells)))
}

study_summary <- rbindlist(cell_summaries, use.names = TRUE,
                           fill = TRUE)
setcolorder(study_summary,
            c('cell_id', 'design_name', 'dropout_name',
              'analysis_name', 'N', 'c_bm', 'true_beta',
              'n_reps', 'n_converged', 'conv_rate',
              'power', 'mcse_power',
              'mean_beta', 'sd_beta', 'mcse_mean_beta',
              'bias', 'mcse_bias',
              'coverage', 'mcse_coverage',
              'mean_dropout'))
save_summary_09(study_summary, STUDY_ID)

cat(sprintf('Study 4 complete: %.1f min.\n',
            as.numeric(difftime(Sys.time(), t_study,
                                units = 'mins'))))
