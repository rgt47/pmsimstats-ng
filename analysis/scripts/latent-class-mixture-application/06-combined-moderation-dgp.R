## Combined-moderation DGP for Study 2 (mean / covariance /
## combined triptych).
##
## Phase A infrastructure for paper 03 path (b). Implements the
## Arminger-Stein-Wittenberg-style simultaneous mean and
## covariance moderation per the ADEMP pre-registration
## (00-ademp-pre-registration.md, §"Study 2"):
##
##   - Pure mean moderation (Architecture A) when alpha = 0.
##   - Pure covariance moderation (Architecture B) when alpha = 1.
##   - Combined form for alpha in (0, 1) with simultaneous mean
##     scaling c_bm_mean = (1 - alpha) * c.bm and covariance
##     moderation c_bm_cov = alpha * c.bm.
##
## The calibration sub-study (07-calibration.R, forthcoming)
## empirically fixes the on-drug marginal Cor(B, BR) constant
## across alpha so that the moderation strength does not vary
## with alpha. The wrapper accepts both raw c.bm values (which
## yield uncalibrated combined moderation) and calibrated
## c_bm_mean / c_bm_cov pairs (which the sub-study produces).
##
## Implemented as a wrapper around generateData() rather than as
## a new dgp_architecture value in the package source, so it can
## be iterated without touching R/generateData.R until the DGP
## stabilises.

library(data.table)

##' Generate data under the Study 2 combined-moderation DGP.
##'
##' @param modelparam Same one-row data.table as generateData(),
##'   with these extra fields (NULL-safe defaults applied):
##'   \itemize{
##'     \item alpha Mixing parameter in [0, 1]. alpha = 0 yields
##'           pure mean moderation (Architecture A); alpha = 1
##'           yields pure covariance moderation (Architecture B).
##'           Default 0.5.
##'     \item c_bm_mean Calibrated mean-moderation strength.
##'           If NULL, computed as (1 - alpha) * modelparam$c.bm.
##'     \item c_bm_cov Calibrated covariance-moderation strength.
##'           If NULL, computed as alpha * modelparam$c.bm.
##'     \item calibration_mode One of 'naive' (default) or
##'           'calibrated'. In 'naive' mode the c_bm_mean and
##'           c_bm_cov values are used directly. In 'calibrated'
##'           mode the values are read from the calibration table
##'           at `analysis/data/quick-sim/03-calibration-table.rds`
##'           if present.
##'   }
##' @param respparam As for generateData().
##' @param blparam As for generateData().
##' @param trialdesign As for generateData().
##' @param empirical,makePositiveDefinite,seed,lambda_cor,verbose,
##'   cached_sigma As for generateData().
##' @return Same long-format data.table as generateData(), with
##'   two extra columns: alpha (mixing parameter used) and
##'   calibration_mode (mode used).
generate_combined_moderation_data <- function(
  modelparam, respparam, blparam, trialdesign,
  empirical = FALSE, makePositiveDefinite = TRUE,
  seed = NA, lambda_cor = NA, verbose = FALSE,
  cached_sigma = NULL
) {

  defaults <- list(
    alpha = 0.5,
    calibration_mode = 'naive'
  )
  for (nm in names(defaults)) {
    if (is.null(modelparam[[nm]])) {
      modelparam[[nm]] <- defaults[[nm]]
    }
  }

  alpha <- modelparam$alpha
  if (alpha < 0 || alpha > 1) {
    stop('alpha must be in [0, 1]; got ', alpha)
  }

  ## Determine c_bm_mean and c_bm_cov.
  if (is.null(modelparam$c_bm_mean) || is.null(modelparam$c_bm_cov)) {
    if (modelparam$calibration_mode == 'calibrated') {
      cal_path <- file.path('analysis', 'data', 'quick-sim',
                            '03-calibration-table.rds')
      if (file.exists(cal_path)) {
        cal <- readRDS(cal_path)
        ## Look up the (alpha, c.bm) pair in the calibration
        ## table; fall back to naive computation if no match.
        match_row <- cal[abs(cal$alpha - alpha) < 1e-6 &
                          abs(cal$c_bm_target - modelparam$c.bm) < 1e-6, ]
        if (nrow(match_row) == 1) {
          modelparam$c_bm_mean <- match_row$c_bm_mean
          modelparam$c_bm_cov  <- match_row$c_bm_cov
        } else {
          warning('calibration entry not found; using naive split')
          modelparam$c_bm_mean <- (1 - alpha) * modelparam$c.bm
          modelparam$c_bm_cov  <- alpha * modelparam$c.bm
        }
      } else {
        warning('calibration table not found; using naive split')
        modelparam$c_bm_mean <- (1 - alpha) * modelparam$c.bm
        modelparam$c_bm_cov  <- alpha * modelparam$c.bm
      }
    } else {
      modelparam$c_bm_mean <- (1 - alpha) * modelparam$c.bm
      modelparam$c_bm_cov  <- alpha * modelparam$c.bm
    }
  }

  ## Stage 1. Draw the MVN backbone with c.bm = c_bm_cov so the
  ## covariance-moderation channel uses the alpha-scaled
  ## correlation.
  mp <- copy(modelparam)
  mp$c.bm <- modelparam$c_bm_cov
  dat <- generateData(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, empirical = empirical,
    makePositiveDefinite = makePositiveDefinite, seed = seed,
    lambda_cor = lambda_cor, verbose = verbose,
    cached_sigma = cached_sigma, dgp_architecture = 'mvn'
  )

  ## Stage 2. Apply the mean-moderation post-draw shift on BR at
  ## on-drug timepoints, scaled by c_bm_mean and the centred
  ## biomarker. This is exactly the Architecture A mechanism
  ## but with c_bm_mean substituted for c.bm.
  if (modelparam$c_bm_mean != 0) {
    bm_mean <- blparam[cat == 'bm']$m
    bm_sd   <- blparam[cat == 'bm']$sd
    bm_z    <- (dat$bm - bm_mean) / bm_sd
    sd_BR   <- respparam[cat == 'br']$sd

    d <- data.table(trialdesign)
    d[, onDrug := (tod > 0)]
    tnames <- trialdesign$timeptnames

    shift_per_subject <- modelparam$c_bm_mean * bm_z * sd_BR

    for (tp in seq_along(tnames)) {
      if (d[tp]$onDrug) {
        br_col <- paste(tnames[tp], 'br', sep = '.')
        dat[, (br_col) := get(br_col) + shift_per_subject]
      }
    }
  }

  ## Stage 3. Recompute total drug-effect columns and symptom-
  ## score columns (BL - sum of components).
  cl <- c('br', 'tv', 'pb')
  tnames <- trialdesign$timeptnames
  for (tp in seq_along(tnames)) {
    n <- tnames[tp]
    comps <- paste(n, cl, sep = '.')
    dat[, (paste0('D_', n)) := rowSums(.SD), .SDcols = comps]
    dat[, (n) := BL - rowSums(.SD), .SDcols = comps]
  }

  dat[, alpha := alpha]
  dat[, calibration_mode := modelparam$calibration_mode]

  dat
}
