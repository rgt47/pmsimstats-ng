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
#'   [generate_data()].
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
#' @return A list with elements `results` (tibble of per-replicate
#'   beta/SE/p) and `parameterselections` (record of inputs).
#' @export
generate_simulated_results <- function(
    trialdesigns, respparamsets, blparamsets,
    censorparams, modelparams, simparam,
    analysisparams, rawdataout = FALSE,
    lambda_cor = NA, n_cores = 1, dgp_architecture = "mvn") {

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
  for (i_r in seq_len(n_combos)) {
    td_paths <- extract_paths(
      trialdesigns[[vpg_master[i_r, "trialdesign"]]]
    )
    rp <- respparamsets[[vpg_master[i_r, "respparamset"]]]$param
    bpp <- blparamsets[[vpg_master[i_r, "blparamset"]]]$param
    mpp <- modelparams[vpg_master[i_r, "modelparamset"], ]

    for (i_p in seq_along(td_paths)) {
      cache_key <- paste(
        vpg_master[i_r, "trialdesign"],
        vpg_master[i_r, "respparamset"],
        vpg_master[i_r, "blparamset"],
        vpg_master[i_r, "modelparamset"],
        i_p, sep = "_"
      )
      if (is.null(sigma_cache[[cache_key]])) {
        td_path <- ensure_timepoint_name(td_paths[[i_p]])
        sigma_cache[[cache_key]] <- build_sigma_matrix(
          mpp, rp, bpp, td_path,
          lambda_cor = lambda_cor,
          dgp_architecture = dgp_architecture
        )
      }
    }
  }
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
      rp <- respparamsets[[vpg[i_r, "respparamset"]]]$param
      bpp <- blparamsets[[vpg[i_r, "blparamset"]]]$param
      mpp <- modelparams[vpg[i_r, "modelparamset"], ]

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

      make_result_row <- function(ap_row, et_out, a_out, i_s,
                                  censorparamset, frac_na) {
        dplyr::bind_cols(
          tibble::as_tibble(vpg[i_r, , drop = FALSE]),
          tibble::as_tibble(as.list(mpp_copy)),
          tibble::tibble(
            censorparamset = censorparamset,
            use_DE = ap_row$useDE,
            t_random_slope = ap_row$t_random_slope,
            irep = (i_r - 1) * simparam$Nreps + i_s,
            frac_NA = frac_na,
            ETbeta = et_out$beta,
            ETbetaSE = et_out$betaSE,
            beta = a_out$beta,
            betaSE = a_out$betaSE,
            p = a_out$p,
            issingular = a_out$issingular,
            warning = a_out$warning
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
            ap_row, et_out, a_out, i_s, 0, 0
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
                ap_row, et_out, a_out, i_s, i_c, frac_na
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
      )
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
