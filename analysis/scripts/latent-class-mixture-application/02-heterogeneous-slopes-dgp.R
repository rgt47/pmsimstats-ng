## Heterogeneous-random-slopes DGP wrapper.
##
## Phase 1 implementation of the Study 3 DGP: a two-component
## mixture over the random-slope distribution on the drug
## indicator. Each participant is assigned to a latent class
## via a logistic gating function on the (standardised)
## biomarker; class-conditional drug response is scaled by
## class-specific multipliers on the BR factor.
##
## Implemented as a wrapper around generateData() rather than as
## a new dgp_architecture value in the package source, so it can
## be iterated without touching R/generateData.R until the DGP
## stabilises. Promotion into the package is straightforward
## once Study 3 cells are run and the parameterisation is
## confirmed.

library(data.table)

##' Generate data under the heterogeneous-random-slopes DGP.
##'
##' @param modelparam Same one-row data.table as generateData(),
##'   with these extra fields (NULL-safe defaults applied):
##'   \itemize{
##'     \item beta_R_factor BR-factor multiplier for class 1
##'           (responder). Default 1.0 (full drug response).
##'     \item beta_NR_factor BR-factor multiplier for class 0
##'           (non-responder). Default 0.0 (no drug response).
##'     \item gating_intercept Logistic-gating intercept; the
##'           log-odds of class 1 at the population mean of bm.
##'           Default 0 (50/50 at the biomarker mean).
##'     \item gating_slope Logistic-gating slope on standardised
##'           biomarker. Default 1.0.
##'   }
##' @param respparam As for generateData().
##' @param blparam As for generateData().
##' @param trialdesign As for generateData().
##' @param empirical,makePositiveDefinite,seed,lambda_cor,verbose,
##'   cached_sigma As for generateData().
##' @return Same long-format data.table as generateData() under
##'   the MVN architecture, with two extra columns: latent_class
##'   (0/1 per participant) and gating_prob (P(class=1) per
##'   participant). Symptom scores incorporate the class-specific
##'   BR scaling at on-drug timepoints.
generate_random_slopes_data <- function(
  modelparam, respparam, blparam, trialdesign,
  empirical = FALSE, makePositiveDefinite = TRUE,
  seed = NA, lambda_cor = NA, verbose = FALSE,
  cached_sigma = NULL
) {

  ## Apply defaults for the new fields.
  defaults <- list(
    beta_R_factor = 1.0, beta_NR_factor = 0.0,
    gating_intercept = 0.0, gating_slope = 1.0
  )
  for (nm in names(defaults)) {
    if (is.null(modelparam[[nm]])) {
      modelparam[[nm]] <- defaults[[nm]]
    }
  }

  ## Stage 1. Draw the MVN backbone with c.bm forced to 0 so the
  ## biomarker-response covariance is generated entirely by the
  ## class-membership channel rather than by the MVN coupling. The
  ## downstream calibration sub-study fixes the marginal Cor(B, BR)
  ## to a target by adjusting gating_slope, beta_R_factor, and
  ## beta_NR_factor jointly; c.bm is forced to 0 here so that
  ## generateData()'s MVN cross-correlation does not double-count.
  mp <- copy(modelparam)
  mp$c.bm <- 0
  ## Suppress the package's mean_moderation path; the random-
  ## slopes DGP modulates BR through class membership, not through
  ## the mean-shift mechanism in generateData.
  dat <- generateData(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, empirical = empirical,
    makePositiveDefinite = makePositiveDefinite, seed = seed,
    lambda_cor = lambda_cor, verbose = verbose,
    cached_sigma = cached_sigma, dgp_architecture = "mvn"
  )

  ## Stage 2. Draw latent class labels.
  bm_mean <- blparam[cat == "bm"]$m
  bm_sd   <- blparam[cat == "bm"]$sd
  bm_z    <- (dat$bm - bm_mean) / bm_sd
  logit_p <- modelparam$gating_intercept +
             modelparam$gating_slope * bm_z
  p_class1 <- plogis(logit_p)
  z <- rbinom(length(p_class1), size = 1, prob = p_class1)
  dat[, gating_prob := p_class1]
  dat[, latent_class := z]

  ## Stage 3. Apply class-specific BR scaling at on-drug timepoints.
  ## Class 1 multiplies BR by beta_R_factor; class 0 by beta_NR_factor.
  ## At off-drug timepoints, BR is left as drawn from the MVN.
  ## We modify the per-timepoint .br component, then reconstruct
  ## the symptom-score columns via the same arithmetic that
  ## generateData() uses.
  d <- data.table(trialdesign)
  d[, onDrug := (tod > 0)]
  tnames <- trialdesign$timeptnames

  ## buildSigma's labels for the three response factors.
  cl <- c("br", "tv", "pb")

  for (tp in seq_along(tnames)) {
    if (d[tp]$onDrug) {
      br_col <- paste(tnames[tp], "br", sep = ".")
      ## Multiplicative class-specific scaling.
      scale_vec <- ifelse(z == 1, modelparam$beta_R_factor,
                                  modelparam$beta_NR_factor)
      dat[, (br_col) := get(br_col) * scale_vec]
    }
  }

  ## Recompute total drug-effect columns and symptom-score columns
  ## (BL - sum of components), mirroring generateData()'s tail.
  for (tp in seq_along(tnames)) {
    n <- tnames[tp]
    comps <- paste(n, cl, sep = ".")
    dat[, (paste0("D_", n)) := rowSums(.SD), .SDcols = comps]
    dat[, (n) := BL - rowSums(.SD), .SDcols = comps]
  }

  dat
}
