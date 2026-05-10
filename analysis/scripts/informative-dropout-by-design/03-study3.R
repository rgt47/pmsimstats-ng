#' 03-study3.R
#'
#' Paper 09 Study 3: randomization-path 'happy accident' decomposition.
#'
#' Per 00-ademp-pre-registration.md: 2 designs (OLBDC, Hybrid) x
#' 2 patterns (balanced, more_of_biased) x 3 N (35, 70, 100) =
#' 12 cells. c_bm = 0.3 throughout. 1000 reps per cell.
#'
#' For each cell we compute path-conditional power: among replicates
#' assigned to each randomization path, what is P(p < 0.05)? The
#' weighted average uses the empirical post-dropout path distribution
#' (i.e., the distribution that survives dropout); the unweighted
#' average uses the design's nominal path distribution. Their
#' difference is the 'happy accident' selection effect.
#'
#' Path is recorded per participant via the path index returned by
#' generateData (di[, path := g]). Within each replicate the path
#' distribution is the participant-level mode (each participant is
#' assigned exactly one path at simulation time). For OLBDC we use
#' the single-path 'pathA'; the path heterogeneity arises from
#' which BDC sub-block holds the participant when dropout occurs,
#' which we capture via the first-drop visit index.
#'
#' Seed: study_seed = 42 + 9 * 100 + 3 = 945.

source(here::here('analysis', 'scripts',
                  'informative-dropout-by-design', '00-common.R'))

set.seed(945)
STUDY_ID   <- '3'
STUDY_SEED <- 945L

design_levels  <- c('OLBDC', 'Hybrid')
dropout_levels <- c('balanced', 'more_of_biased')
N_levels       <- c(35L, 70L, 100L)

cells <- CJ(design_name  = design_levels,
            dropout_name = dropout_levels,
            N            = N_levels,
            sorted = FALSE)
cells[, c_bm    := 0.3]
cells[, n_reps  := 1000L]
cells[, cell_id := .I]

cat(sprintf('Study 3: %d cells; total target reps = %d.\n',
            nrow(cells), sum(cells$n_reps)))

#=====================================================================
# Path-aware replicate executor
#=====================================================================
#
# Per-replicate path tracking: returns one row per participant with
# the path index (for hybrid: AB or BA) and the first-drop visit
# index (for OLBDC: integer in 1..20 or NA). Plus the population-
# level lme estimate. Path-conditional power is then computed at
# cell level by stratifying replicates by the modal path index.

one_rep_S3 <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- ALL_DESIGNS()[[cell$design_name]]
  paths <- td$trialpaths
  mp <- make_model_params_09(N = cell$N, c_bm = cell$c_bm)
  rp <- make_resp_params_09()
  bp <- make_bl_params_09()
  drop <- DROPOUT_PATTERNS[[cell$dropout_name]]

  dat_list <- vector('list', length(paths))
  drop_fracs <- numeric(length(paths))
  first_drops_list <- vector('list', length(paths))
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
    first_drops_list[[g]] <- attr(di_d, 'first_drop')
    dat_list[[g]] <- di_d
  }
  good <- !vapply(dat_list, is.null, logical(1))
  if (!any(good)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, dropout_fraction = NA_real_,
                      modal_path = NA_integer_,
                      mean_first_drop = NA_real_,
                      rep_idx = rep_idx))
  }
  dat <- rbindlist(dat_list[good], fill = TRUE)
  res <- lme_fit_09(td, dat)

  ## Modal path: which of the design paths most participants follow.
  ## At simulation time pmsimstats assigns participants uniformly
  ## over paths; for two-path designs the modal path is essentially
  ## random per replicate. For OLBDC there is only one path; we
  ## record the mean first-drop visit instead.
  modal_path <- if (length(paths) > 1L) {
    pcounts <- vapply(seq_along(paths), function(g)
      if (is.null(dat_list[[g]])) 0L else nrow(dat_list[[g]]),
      integer(1))
    which.max(pcounts)
  } else {
    1L
  }
  all_first_drops <- unlist(first_drops_list[good])
  mean_first_drop <- mean(all_first_drops, na.rm = TRUE)

  data.table(
    beta = if ('beta' %in% names(res)) res$beta else NA_real_,
    betaSE = if ('betaSE' %in% names(res)) res$betaSE else NA_real_,
    p = if ('p' %in% names(res)) res$p else NA_real_,
    dropout_fraction = mean(drop_fracs[good]),
    modal_path = modal_path,
    mean_first_drop = mean_first_drop,
    rep_idx = rep_idx)
}

run_cell_S3 <- function(cell, n_reps, study_seed, cell_id,
                        n_cores = max(1L,
                                      parallel::detectCores() - 1L),
                        chunk_size = 50L) {
  results <- vector('list', 0L)
  done <- 0L
  while (done < n_reps) {
    take <- min(chunk_size, n_reps - done)
    chunk <- parallel::mclapply(
      seq.int(done + 1L, done + take),
      function(r) one_rep_S3(r, cell, study_seed, cell_id),
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
             N = cell$N, c_bm = cell$c_bm,
             cell_id = cell_id)]
  out
}

t_study <- Sys.time()
cell_summaries <- vector('list', nrow(cells))

for (i in seq_len(nrow(cells))) {
  cell <- list(design_name  = cells$design_name[i],
               dropout_name = cells$dropout_name[i],
               N            = cells$N[i],
               c_bm         = cells$c_bm[i])
  reps_dt <- run_cell_S3(cell,
                         n_reps    = cells$n_reps[i],
                         study_seed = STUDY_SEED,
                         cell_id   = cells$cell_id[i])
  save_cell_09(reps_dt, STUDY_ID, cells$cell_id[i])

  ## Path-conditional power.
  conv <- reps_dt[!is.na(p)]
  if (nrow(conv) == 0) next

  ## Empirical (post-dropout) path-distribution-weighted power
  by_path <- conv[, .(power = mean(p < 0.05),
                      n = .N,
                      mean_dropout = mean(dropout_fraction)),
                  by = modal_path]
  emp_power <- weighted.mean(by_path$power,
                             w = by_path$n / sum(by_path$n))

  ## Nominal (uniform across paths) weighted power
  if (nrow(by_path) > 1L) {
    nom_power <- mean(by_path$power)
  } else {
    nom_power <- by_path$power
  }
  happy_accident_diff <- emp_power - nom_power

  s <- data.table(
    cell_id = cells$cell_id[i],
    design_name = cells$design_name[i],
    dropout_name = cells$dropout_name[i],
    N = cells$N[i], c_bm = cells$c_bm[i],
    n_reps = nrow(reps_dt),
    n_converged = nrow(conv),
    emp_weighted_power = emp_power,
    nom_weighted_power = nom_power,
    happy_accident_diff = happy_accident_diff,
    mean_dropout = mean(conv$dropout_fraction, na.rm = TRUE),
    paths_observed = nrow(by_path))
  cell_summaries[[i]] <- s
  cat(sprintf('Study 3 cell %d/%d done.\n', i, nrow(cells)))
}

study_summary <- rbindlist(cell_summaries, use.names = TRUE,
                           fill = TRUE)
save_summary_09(study_summary, STUDY_ID)

cat(sprintf('Study 3 complete: %.1f min.\n',
            as.numeric(difftime(Sys.time(), t_study,
                                units = 'mins'))))
