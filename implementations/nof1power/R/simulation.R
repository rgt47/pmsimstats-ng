#' Orchestrate simulation over a parameter grid
#'
#' Iterates over the Cartesian product of trial designs, response
#' parameter sets, baseline parameter sets, and model parameter rows;
#' caches sigma matrices per (design, path) combination; and runs
#' [generate_data()] plus [lme_analysis()] for each cell, optionally
#' with and without censoring. Direct port of the pmsimstats-ng
#' tidyverse `generate_simulated_results`.
#'
#' @param trialdesigns List of trial-design objects, each the output
#'   of [build_trial_design()] (a `list(metadata, trialpaths)`).
#' @param respparamsets List of response-parameter wrappers, each a
#'   list with element `param` that is a tibble accepted by
#'   [generate_data()]. A wrapper may also carry `target_delta`, a
#'   number. The drug-response ceiling (`max` of the `br` row) is then
#'   set, separately for every trial design, so that the realized
#'   estimand `delta_on` (the mean drug response over the on-drug
#'   visits of the design, see `.realized_delta_on()`) equals
#'   `target_delta`. This puts designs with different exposure
#'   schedules on the same effect. `target_delta = 0` is the null.
#' @param blparamsets List of baseline-parameter wrappers, each a
#'   list with element `param` that is a tibble accepted by
#'   [generate_data()].
#' @param censorparams Either `NA` (no censoring) or a data frame of
#'   censoring parameters accepted by [censor_data()].
#' @param modelparams Data frame of model parameters, one row per
#'   parameter set.
#' @param simparam List with simulation controls:
#'   \itemize{
#'     \item `Nreps`: number of replicates per parameter set.
#'     \item `progressiveSave`: logical, save in chunks.
#'     \item `basesavename`: file prefix.
#'     \item `nRep2save`: chunk size.
#'     \item `saveunit2start`: chunk to start from.
#'     \item `savedir`: directory for saved files.
#'   }
#' @param analysisparams Data frame of analysis options (one row per
#'   analysis variant), each with columns matching the `op` list in
#'   [lme_analysis()].
#' @param rawdataout Logical; return raw simulated data.
#' @param lambda_cor Correlation decay rate; see [generate_data()].
#' @param n_cores Number of parallel workers; `< 1` auto-detects.
#'   Currently sequential; reserved for future `furrr` integration.
#' @param dgp_architecture Either `"mvn"` or `"mean_moderation"`.
#' @param skip_infeasible Logical. If `FALSE` (default), a
#'   non-positive-definite covariance matrix in any cell stops the
#'   run. If `TRUE`, such cells are dropped, recorded in the
#'   `infeasible` element of the return value with the reason and the
#'   smallest eigenvalue, and the remaining cells run. No matrix is
#'   ever repaired.
#' @return A list with elements `results` (tibble of per-replicate
#'   beta/SE/p, the fitted formula, the realized estimand
#'   `delta_on`, the feasibility margin `min_eig`, `pd_repaired`, and
#'   the analysis variant `analysis_row`, `analysis_t1half`),
#'   `parameterselections` (record of inputs), `cells` (one row per
#'   run cell with its feasibility margin) and `infeasible` (cells
#'   dropped under `skip_infeasible = TRUE`).
#' @export
generate_simulated_results <- function(
    trialdesigns, respparamsets, blparamsets,
    censorparams, modelparams, simparam,
    analysisparams, rawdataout = FALSE,
    lambda_cor = NA, n_cores = 1, dgp_architecture = "mvn",
    skip_infeasible = FALSE) {

  dgp_architecture <- match.arg(dgp_architecture,
                                c("mvn", "mean_moderation"))

  if (missing(analysisparams)) {
    analysisparams <- data.frame(
      useDE = TRUE, t_random_slope = FALSE,
      full_model_out = FALSE
    )
  }

  if (n_cores < 1) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }

  vpg_master <- expand.grid(
    trialdesign = seq_along(trialdesigns),
    respparamset = seq_along(respparamsets),
    blparamset = seq_along(blparamsets),
    modelparamset = seq_len(nrow(modelparams))
  )
  n_combos <- nrow(vpg_master)

  extract_paths <- function(td_obj) {
    if (is.list(td_obj) && "trialpaths" %in% names(td_obj)) {
      td_obj$trialpaths
    } else if (is.list(td_obj) && length(td_obj) >= 2) {
      td_obj[[2]]
    } else {
      list(td_obj)
    }
  }

  message("Caching sigma matrices...")
  sigma_cache <- list()
  rp_cache <- list()
  cell_info <- vector("list", n_combos)
  for (i_r in seq_len(n_combos)) {
    td_paths <- extract_paths(
      trialdesigns[[vpg_master[i_r, "trialdesign"]]]
    )
    bpp <- blparamsets[[vpg_master[i_r, "blparamset"]]]$param
    mpp <- modelparams[vpg_master[i_r, "modelparamset"], ]
    rp <- .resolve_resp(
      respparamsets[[vpg_master[i_r, "respparamset"]]],
      td_paths, bpp, mpp, lambda_cor, dgp_architecture
    )
    rp_cache[[paste(
      vpg_master[i_r, "trialdesign"], vpg_master[i_r, "respparamset"],
      vpg_master[i_r, "blparamset"], vpg_master[i_r, "modelparamset"],
      sep = "_")]] <- rp

    min_eig <- Inf
    reason <- NA_character_
    for (i_p in seq_along(td_paths)) {
      cache_key <- paste(
        vpg_master[i_r, "trialdesign"],
        vpg_master[i_r, "respparamset"],
        vpg_master[i_r, "blparamset"],
        vpg_master[i_r, "modelparamset"],
        i_p, sep = "_"
      )
      td_path <- ensure_timepoint_name(td_paths[[i_p]])
      built <- tryCatch(
        build_sigma_matrix(
          mpp, rp, bpp, td_path,
          lambda_cor = lambda_cor,
          dgp_architecture = dgp_architecture
        ),
        error = function(e) e
      )
      if (inherits(built, "error")) {
        reason <- conditionMessage(built)
        if (!skip_infeasible) stop(built)
        raw <- tryCatch(
          build_sigma_matrix(
            mpp, rp, bpp, td_path, lambda_cor = lambda_cor,
            dgp_architecture = dgp_architecture,
            stop_on_nonpd = FALSE
          ),
          error = function(e) NULL
        )
        min_eig <- if (is.null(raw)) NA_real_ else
          min(min_eig, .min_eigenvalue(raw))
      } else {
        sigma_cache[[cache_key]] <- built
        min_eig <- min(min_eig, .min_eigenvalue(built))
      }
    }
    cell_info[[i_r]] <- tibble::tibble(
      trialdesign = vpg_master[i_r, "trialdesign"],
      respparamset = vpg_master[i_r, "respparamset"],
      blparamset = vpg_master[i_r, "blparamset"],
      modelparamset = vpg_master[i_r, "modelparamset"],
      feasible = is.na(reason),
      min_eig = if (is.infinite(min_eig)) NA_real_ else min_eig,
      reason = reason
    )
  }
  cells <- dplyr::bind_rows(cell_info)
  infeasible <- cells[!cells$feasible, , drop = FALSE]
  if (nrow(infeasible) > 0) {
    message(sprintf("Dropped %d infeasible cell(s) of %d",
                    nrow(infeasible), n_combos))
    vpg_master <- vpg_master[cells$feasible, , drop = FALSE]
    n_combos <- nrow(vpg_master)
  }
  if (n_combos == 0) stop("Every cell is infeasible; nothing to run.")
  message(sprintf("Cached %d unique sigma matrices",
                  length(sigma_cache)))

  if (!simparam$progressiveSave) {
    n_large_loops <- 1
    ll_starts <- 1
    ll_stops <- n_combos
  } else {
    n_large_loops <- ceiling(n_combos / simparam$nRep2save)
    if (n_large_loops > 1) {
      ll_starts <- c(1, 1 + simparam$nRep2save *
                       (1:(n_large_loops - 1)))
      ll_stops <- c(ll_starts[2:n_large_loops] - 1, n_combos)
    } else {
      ll_starts <- 1
      ll_stops <- n_combos
    }
  }

  outpt <- NULL

  for (i_ll in simparam$saveunit2start:n_large_loops) {
    vpg <- vpg_master[ll_starts[i_ll]:ll_stops[i_ll], , drop = FALSE]
    i_paramset <- ll_starts[i_ll]
    n_r <- nrow(vpg)

    run_one <- function(i_r) {
      td_paths <- extract_paths(
        trialdesigns[[vpg[i_r, "trialdesign"]]]
      )
      bpp <- blparamsets[[vpg[i_r, "blparamset"]]]$param
      mpp <- modelparams[vpg[i_r, "modelparamset"], ]
      rp <- rp_cache[[paste(
        vpg[i_r, "trialdesign"], vpg[i_r, "respparamset"],
        vpg[i_r, "blparamset"], vpg[i_r, "modelparamset"],
        sep = "_")]]

      n_paths <- length(td_paths)
      Ns <- mpp$N %/% n_paths
      Ns <- Ns + c(
        rep(1, mpp$N %% n_paths),
        rep(0, n_paths - mpp$N %% n_paths)
      )
      NNs <- Ns * simparam$Nreps

      mpp_copy <- mpp
      mpp_copy$N <- NNs[[1]]
      ck <- paste(
        vpg[i_r, "trialdesign"],
        vpg[i_r, "respparamset"],
        vpg[i_r, "blparamset"],
        vpg[i_r, "modelparamset"],
        1, sep = "_"
      )

      td_path <- ensure_timepoint_name(td_paths[[1]])
      dat <- generate_data(
        mpp_copy, rp, bpp, td_path,
        empirical = FALSE, make_positive_definite = TRUE,
        cached_sigma = sigma_cache[[ck]],
        dgp_architecture = dgp_architecture
      )
      dat$path <- 1
      dat$replicate <- rep(seq_len(simparam$Nreps), Ns[1])

      if (n_paths > 1) {
        for (i_p in 2:n_paths) {
          mpp_copy$N <- NNs[[i_p]]
          ck <- paste(
            vpg[i_r, "trialdesign"],
            vpg[i_r, "respparamset"],
            vpg[i_r, "blparamset"],
            vpg[i_r, "modelparamset"],
            i_p, sep = "_"
          )
          td_p <- ensure_timepoint_name(td_paths[[i_p]])
          dat2 <- generate_data(
            mpp_copy, rp, bpp, td_p,
            empirical = FALSE, make_positive_definite = TRUE,
            cached_sigma = sigma_cache[[ck]],
            dgp_architecture = dgp_architecture
          )
          dat2$path <- i_p
          dat2$replicate <- rep(seq_len(simparam$Nreps), Ns[i_p])
          dat <- dplyr::bind_rows(dat, dat2)
        }
      }
      mpp_copy$N <- sum(Ns)

      cell_row <- cells[
        cells$trialdesign == vpg[i_r, "trialdesign"] &
          cells$respparamset == vpg[i_r, "respparamset"] &
          cells$blparamset == vpg[i_r, "blparamset"] &
          cells$modelparamset == vpg[i_r, "modelparamset"], ]
      delta_on <- .realized_delta_on(
        vpg[i_r, ], td_paths, sigma_cache
      )
      pd_repaired <- any(vapply(seq_len(n_paths), function(i_p) {
        ck_p <- paste(
          vpg[i_r, "trialdesign"], vpg[i_r, "respparamset"],
          vpg[i_r, "blparamset"], vpg[i_r, "modelparamset"],
          i_p, sep = "_")
        isTRUE(sigma_cache[[ck_p]]$was_pd_corrected)
      }, logical(1)))
      cell_cols <- tibble::tibble(
        br_max = rp$max[rp$cat == "br"],
        br_sd = rp$sd[rp$cat == "br"],
        pb_max = rp$max[rp$cat == "pb"],
        tv_max = rp$max[rp$cat == "tv"],
        bl_sd = bpp$sd[bpp$cat == "BL"],
        n_visits = nrow(td_paths[[1]]),
        n_paths = n_paths,
        n_switch = max(vapply(td_paths, function(tp) {
          on <- prepare_trial_data(ensure_timepoint_name(tp))$on_drug
          sum(diff(as.integer(on)) != 0)
        }, numeric(1))),
        delta_on = delta_on,
        min_eig = cell_row$min_eig[1],
        pd_repaired = pd_repaired
      )

      make_result_row <- function(ap_row, et_out, a_out, i_s,
                                  censorparamset, frac_na, i_ap) {
        dplyr::bind_cols(
          tibble::as_tibble(vpg[i_r, , drop = FALSE]),
          tibble::as_tibble(as.list(mpp_copy)),
          cell_cols,
          tibble::tibble(
            censorparamset = censorparamset,
            analysis_row = i_ap,
            analysis_t1half = ap_row$carryover_t1half %||% 0,
            use_DE = ap_row$useDE,
            t_random_slope = ap_row$t_random_slope,
            irep = (i_r - 1) * simparam$Nreps + i_s,
            frac_NA = frac_na,
            ETbeta = et_out$beta,
            ETbetaSE = et_out$betaSE,
            beta = a_out$beta,
            betaSE = a_out$betaSE,
            betaDF = a_out$betaDF,
            p = a_out$p,
            issingular = a_out$issingular,
            warning = a_out$warning,
            fit_formula = a_out$fit_formula,
            car1_dropped = a_out$car1_dropped
          )
        )
      }

      results_list <- list()
      ridx <- 0
      nocensoring <- length(censorparams) < 2 &&
        all(is.na(censorparams))

      for (i_ap in seq_len(nrow(analysisparams))) {
        ap_row <- as.list(analysisparams[i_ap, ])
        et_out <- lme_analysis(td_paths, dat, ap_row)
        for (i_s in seq_len(simparam$Nreps)) {
          dat_rep <- dat |>
            dplyr::filter(.data$replicate == i_s)
          a_out <- lme_analysis(td_paths, dat_rep, ap_row)
          ridx <- ridx + 1
          results_list[[ridx]] <- make_result_row(
            ap_row, et_out, a_out, i_s, 0, 0, i_ap
          )
        }

        if (!nocensoring) {
          for (i_c in seq_len(nrow(censorparams))) {
            datc <- censor_data(dat, td_paths[[1]],
                                censorparams[i_c, ])
            for (i_s in seq_len(simparam$Nreps)) {
              datc_rep <- datc |>
                dplyr::filter(.data$replicate == i_s)
              frac_na <- sum(is.na(datc_rep)) /
                (mpp_copy$N * nrow(td_paths[[1]]))
              a_out <- lme_analysis(td_paths, datc_rep, ap_row)
              ridx <- ridx + 1
              results_list[[ridx]] <- make_result_row(
                ap_row, et_out, a_out, i_s, i_c, frac_na, i_ap
              )
            }
          }
        }
      }
      dplyr::bind_rows(results_list)
    }

    use_parallel <- n_cores > 1 &&
      requireNamespace("future", quietly = TRUE) &&
      requireNamespace("furrr", quietly = TRUE)

    if (use_parallel) {
      prev_plan <- future::plan(
        future::multisession, workers = n_cores
      )
      on.exit(future::plan(prev_plan), add = TRUE)
      options(future.globals.maxSize = 2 * 1024^3)

      out_parts <- furrr::future_map(
        seq_len(n_r),
        function(i_r) run_one(i_r),
        .options = furrr::furrr_options(seed = TRUE)
      )
    } else {
      out_parts <- purrr::map(seq_len(n_r), function(i_r) {
        message(sprintf(
          "Parameter set %d of %d (%d of %d total)",
          i_r, n_r, i_paramset + i_r - 1, n_combos
        ))
        run_one(i_r)
      })
    }
    out <- dplyr::bind_rows(out_parts)

    outpt <- list(
      results = out,
      parameterselections = list(
        trialdesigns = trialdesigns,
        respparamsets = respparamsets,
        blparamsets = blparamsets,
        censorparams = censorparams,
        modelparams = modelparams,
        analysisparams = analysisparams,
        simparam = simparam
      ),
      cells = cells,
      infeasible = infeasible
    )

    if (simparam$progressiveSave) {
      saveRDS(
        outpt,
        file.path(
          simparam$savedir,
          paste(simparam$basesavename, i_ll, sep = "_save")
        )
      )
    }
  }

  outpt
}


#' Smallest eigenvalue of the positive-variance block of a covariance
#'
#' @param built Output of [build_sigma_matrix()].
#' @return The smallest eigenvalue of the correlation-scaled
#'   covariance restricted to components with positive variance.
#' @keywords internal
.min_eigenvalue <- function(built) {
  pos <- built$standard_deviations > 0
  sub <- built$sigma[pos, pos, drop = FALSE]
  if (anyNA(sub) || any(!is.finite(sub))) return(NA_real_)
  min(eigen(stats::cov2cor(sub), symmetric = TRUE,
            only.values = TRUE)$values)
}

#' Realized estimand of a cell: mean drug response on drug
#'
#' The mean, over every on-drug follow-up visit of every path (each
#' visit counts once, whatever the path), of the mean drug-response
#' component, with the sign of the symptom change (negative when the
#' drug reduces symptoms). It is a property of the design and the
#' response curve, not of an analysis model, so it is common to all
#' designs.
#'
#' @param vpg_row One row of the cell grid.
#' @param td_paths List of per-path trial-design tibbles.
#' @param sigma_cache Cache of [build_sigma_matrix()] results.
#' @return A single number.
#' @keywords internal
.realized_delta_on <- function(vpg_row, td_paths, sigma_cache) {
  num <- 0
  den <- 0
  for (i_p in seq_along(td_paths)) {
    ck <- paste(vpg_row[["trialdesign"]], vpg_row[["respparamset"]],
                vpg_row[["blparamset"]], vpg_row[["modelparamset"]],
                i_p, sep = "_")
    sg <- sigma_cache[[ck]]
    on <- prepare_trial_data(
      ensure_timepoint_name(td_paths[[i_p]])
    )$on_drug
    br <- utils::tail(sg$means, sg$nP)
    num <- num + sum(br[on])
    den <- den + sum(on)
  }
  if (den == 0) NA_real_ else -num / den
}

#' Resolve a response-parameter wrapper for one trial design
#'
#' Returns `respset$param`, with the drug-response ceiling replaced
#' when the wrapper carries `target_delta` (see
#' [generate_simulated_results()]). Because the on-drug drug-response
#' mean is proportional to the ceiling, one evaluation at ceiling 1
#' fixes it.
#'
#' @keywords internal
.resolve_resp <- function(respset, td_paths, bpp, mpp, lambda_cor,
                          dgp_architecture) {
  rp <- respset$param
  target <- respset$target_delta
  if (is.null(target)) return(rp)
  if (target == 0) {
    rp$max[rp$cat == "br"] <- 0
    return(rp)
  }
  rp1 <- rp
  rp1$max[rp1$cat == "br"] <- 1
  num <- 0
  den <- 0
  for (i_p in seq_along(td_paths)) {
    td_path <- ensure_timepoint_name(td_paths[[i_p]])
    sg <- build_sigma_matrix(
      mpp, rp1, bpp, td_path, lambda_cor = lambda_cor,
      dgp_architecture = dgp_architecture, stop_on_nonpd = FALSE
    )
    on <- prepare_trial_data(td_path)$on_drug
    num <- num + sum(utils::tail(sg$means, sg$nP)[on])
    den <- den + sum(on)
  }
  rp$max[rp$cat == "br"] <- abs(target) / (num / den)
  rp
}
