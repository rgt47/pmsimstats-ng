#' 02-study2.R
#'
#' Paper 09 Study 2: bias quantification (two-mode).
#'
#' Per 00-ademp-pre-registration.md: subset of Study 1 at
#' c_bm in {0.3, 0.6}. 4 designs x 5 dropouts x 3 N x 2 c_bm = 120
#' cells. 1000 reps per cell.
#'
#' For each replicate, we compute Δβ in BOTH modes:
#'   (a) Δβ_a = β_hat^dropped - β_true
#'   (b) Δβ_b = β_hat^dropped - β_hat^gold
#' where β_hat^gold is the analysis estimate from the full
#' uncensored data (i.e., the same simulated dataset before dropout
#' is applied).
#'
#' This requires running the analysis twice per replicate, so the
#' per-replicate cost is roughly double Study 1's per-replicate cost.
#'
#' Seed: study_seed = 42 + 9 * 100 + 2 = 944.

source(here::here('analysis', 'scripts',
                  'informative-dropout-by-design', '00-common.R'))

set.seed(944)
STUDY_ID   <- '2'
STUDY_SEED <- 944L

design_levels  <- c('OL', 'OLBDC', 'CO', 'Hybrid')
dropout_levels <- names(DROPOUT_PATTERNS)
N_levels       <- c(35L, 70L, 100L)
c_bm_levels    <- c(0.3, 0.6)

cells <- CJ(design_name  = design_levels,
            dropout_name = dropout_levels,
            N            = N_levels,
            c_bm         = c_bm_levels,
            sorted = FALSE)
cells[, n_reps  := 1000L]
cells[, cell_id := .I]

cat(sprintf('Study 2: %d cells; total target reps = %d.\n',
            nrow(cells), sum(cells$n_reps)))

#=====================================================================
# Replicate executor for Study 2 (two-mode bias)
#=====================================================================

one_rep_S2 <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- ALL_DESIGNS()[[cell$design_name]]
  paths <- td$trialpaths
  mp <- make_model_params_09(N = cell$N, c_bm = cell$c_bm)
  rp <- make_resp_params_09()
  bp <- make_bl_params_09()
  drop <- DROPOUT_PATTERNS[[cell$dropout_name]]

  ## Generate the FULL data first (no dropout). Analyze it for
  ## the gold-standard estimate. Then apply dropout to a copy and
  ## analyze that for the dropped estimate.
  full_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
    di <- tryCatch(
      generateData(mp, rp, bp, paths[[g]],
                   empirical = FALSE, makePositiveDefinite = TRUE,
                   dgp_architecture = 'mvn'),
      error = function(e) NULL)
    if (is.null(di)) next
    di[, path := g]
    full_list[[g]] <- di
  }
  good <- !vapply(full_list, is.null, logical(1))
  if (!any(good)) {
    return(data.table(beta_full = NA_real_, beta_drop = NA_real_,
                      betaSE_full = NA_real_, betaSE_drop = NA_real_,
                      p_full = NA_real_, p_drop = NA_real_,
                      dropout_fraction = NA_real_,
                      rep_idx = rep_idx))
  }
  dat_full <- rbindlist(full_list[good], fill = TRUE)
  res_full <- lme_fit_09(td, dat_full)

  ## Apply dropout to a copy.
  dropped_list <- vector('list', length(paths))
  drop_fracs <- numeric(length(paths))
  for (g in seq_along(paths)) {
    if (is.null(full_list[[g]])) next
    dj <- copy(full_list[[g]])
    dj_d <- apply_hazard_dropout(dj, paths[[g]],
                                  beta0 = drop$beta0,
                                  beta1 = drop$beta1)
    drop_fracs[g] <- attr(dj_d, 'dropout_fraction') %||% 0
    dropped_list[[g]] <- dj_d
  }
  dat_drop <- rbindlist(dropped_list[good], fill = TRUE)
  res_drop <- lme_fit_09(td, dat_drop)

  data.table(
    beta_full   = if ('beta'   %in% names(res_full)) res_full$beta   else NA_real_,
    betaSE_full = if ('betaSE' %in% names(res_full)) res_full$betaSE else NA_real_,
    p_full      = if ('p'      %in% names(res_full)) res_full$p      else NA_real_,
    beta_drop   = if ('beta'   %in% names(res_drop)) res_drop$beta   else NA_real_,
    betaSE_drop = if ('betaSE' %in% names(res_drop)) res_drop$betaSE else NA_real_,
    p_drop      = if ('p'      %in% names(res_drop)) res_drop$p      else NA_real_,
    dropout_fraction = mean(drop_fracs[good]),
    rep_idx = rep_idx)
}

run_cell_S2 <- function(cell, n_reps, study_seed, cell_id,
                        n_cores = max(1L,
                                      parallel::detectCores() - 1L),
                        chunk_size = 50L) {
  t0 <- Sys.time()
  results <- vector('list', 0L)
  done <- 0L
  while (done < n_reps) {
    take <- min(chunk_size, n_reps - done)
    chunk <- parallel::mclapply(
      seq.int(done + 1L, done + take),
      function(r) one_rep_S2(r, cell, study_seed, cell_id),
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
  cat(sprintf('  Cell %d wall: %.1f s\n', cell_id,
              as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
  out
}

t_study <- Sys.time()
cell_summaries <- vector('list', nrow(cells))

for (i in seq_len(nrow(cells))) {
  cell <- list(design_name  = cells$design_name[i],
               dropout_name = cells$dropout_name[i],
               N            = cells$N[i],
               c_bm         = cells$c_bm[i])
  reps_dt <- run_cell_S2(cell,
                         n_reps    = cells$n_reps[i],
                         study_seed = STUDY_SEED,
                         cell_id   = cells$cell_id[i])
  save_cell_09(reps_dt, STUDY_ID, cells$cell_id[i])

  tb <- TRUE_BETA_09(c_bm = cells$c_bm[i])
  good <- reps_dt[!is.na(beta_drop) & !is.na(beta_full)]
  n <- nrow(good)
  if (n > 0) {
    delta_a <- good$beta_drop - tb
    delta_b <- good$beta_drop - good$beta_full
    s <- data.table(
      cell_id = cells$cell_id[i],
      design_name = cells$design_name[i],
      dropout_name = cells$dropout_name[i],
      N = cells$N[i], c_bm = cells$c_bm[i],
      true_beta = tb,
      n_reps = nrow(reps_dt),
      n_converged = n,
      conv_rate = n / nrow(reps_dt),
      mean_delta_a = mean(delta_a),
      mcse_delta_a = sd(delta_a) / sqrt(n),
      mean_delta_b = mean(delta_b),
      mcse_delta_b = sd(delta_b) / sqrt(n),
      prop_sign_flip_a = mean(sign(delta_a) != sign(tb), na.rm = TRUE),
      prop_sign_flip_b = mean(sign(delta_b) != sign(good$beta_full),
                              na.rm = TRUE),
      mean_dropout = mean(reps_dt$dropout_fraction, na.rm = TRUE))
  } else {
    s <- data.table(
      cell_id = cells$cell_id[i],
      design_name = cells$design_name[i],
      dropout_name = cells$dropout_name[i],
      N = cells$N[i], c_bm = cells$c_bm[i], true_beta = tb,
      n_reps = nrow(reps_dt), n_converged = 0L, conv_rate = 0,
      mean_delta_a = NA_real_, mcse_delta_a = NA_real_,
      mean_delta_b = NA_real_, mcse_delta_b = NA_real_,
      prop_sign_flip_a = NA_real_, prop_sign_flip_b = NA_real_,
      mean_dropout = NA_real_)
  }
  cell_summaries[[i]] <- s
  cat(sprintf('Study 2 cell %d/%d done.\n', i, nrow(cells)))
}

study_summary <- rbindlist(cell_summaries, use.names = TRUE,
                           fill = TRUE)
save_summary_09(study_summary, STUDY_ID)

cat(sprintf('Study 2 complete: %.1f min.\n',
            as.numeric(difftime(Sys.time(), t_study,
                                units = 'mins'))))
