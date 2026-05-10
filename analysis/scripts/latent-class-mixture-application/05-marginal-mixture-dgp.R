## Marginal-mixture DGP for Study 1 (mixture-vs-MVN factorial).
##
## Phase A infrastructure for paper 03 path (b). Implements the
## two-class Bernoulli mixture described in the ADEMP pre-reg
## (00-ademp-pre-registration.md, §"Study 1") with the PB-class-
## correlation axis that the pre-reg requires:
##
##   - PB-class-independent: PB drawn from the population
##     distribution independently of class label; this is the
##     original Hendrickson specification.
##
##   - PB-class-correlated: responder class draws maxr^PB from a
##     distribution shifted upward relative to the non-responder
##     class, with shift calibrated so that Cor(BR, PB) at on-drug
##     timepoints equals a specified target (0.3 or 0.6).
##
## The PB-class-correlated cell is the natural mechanism for
## generating apparent subadditivity in a one-component additive
## analysis. Reporting both cells lets the analyst distinguish
## 'true mechanistic subadditivity' from 'apparent subadditivity
## emerging from unmodelled latent-class structure'.
##
## Implemented as a wrapper around generateData() rather than as a
## new dgp_architecture value in the package source, so it can be
## iterated without touching R/generateData.R until the DGP
## stabilises.

library(data.table)

##' Generate data under the Study 1 marginal-mixture DGP.
##'
##' @param modelparam Same one-row data.table as generateData(),
##'   with these extra fields (NULL-safe defaults applied):
##'   \itemize{
##'     \item beta_R_factor BR-factor multiplier for class 1
##'           (responder). Default 1.0.
##'     \item beta_NR_factor BR-factor multiplier for class 0
##'           (non-responder). Default 0.0.
##'     \item gating_intercept Logistic-gating intercept. Default 0
##'           (class prevalence 0.5 at the biomarker mean).
##'     \item gating_slope Logistic-gating slope on standardised
##'           biomarker. Default 1.0.
##'     \item pb_class_correlation One of 'independent' (default;
##'           PB drawn independently of class) or 'correlated' (PB
##'           shifted by class with calibrated target).
##'     \item pb_class_target Target Cor(BR, PB) at on-drug
##'           timepoints when pb_class_correlation = 'correlated'.
##'           Default 0.3. Only consulted when correlation =
##'           'correlated'.
##'   }
##' @param respparam As for generateData().
##' @param blparam As for generateData().
##' @param trialdesign As for generateData().
##' @param empirical,makePositiveDefinite,seed,lambda_cor,verbose,
##'   cached_sigma As for generateData().
##' @return Same long-format data.table as generateData() under the
##'   MVN architecture, with three extra columns: latent_class
##'   (0/1 per participant), gating_prob (P(class=1) per
##'   participant), and pb_class_shift (numeric shift applied to PB
##'   for class-correlated cells; 0 for class-independent cells).
##'   Symptom scores incorporate the class-specific BR scaling at
##'   on-drug timepoints.
generate_marginal_mixture_data <- function(
  modelparam, respparam, blparam, trialdesign,
  empirical = FALSE, makePositiveDefinite = TRUE,
  seed = NA, lambda_cor = NA, verbose = FALSE,
  cached_sigma = NULL
) {

  defaults <- list(
    beta_R_factor = 1.0, beta_NR_factor = 0.0,
    gating_intercept = 0.0, gating_slope = 1.0,
    pb_class_correlation = 'independent',
    pb_class_target = 0.3
  )
  for (nm in names(defaults)) {
    if (is.null(modelparam[[nm]])) {
      modelparam[[nm]] <- defaults[[nm]]
    }
  }

  ## Stage 1. Draw the MVN backbone with c.bm forced to 0 so the
  ## biomarker-response covariance is generated entirely by the
  ## class-membership channel rather than by the MVN coupling.
  mp <- copy(modelparam)
  mp$c.bm <- 0
  dat <- generateData(
    modelparam = mp, respparam = respparam, blparam = blparam,
    trialdesign = trialdesign, empirical = empirical,
    makePositiveDefinite = makePositiveDefinite, seed = seed,
    lambda_cor = lambda_cor, verbose = verbose,
    cached_sigma = cached_sigma, dgp_architecture = 'mvn'
  )

  ## Stage 2. Draw latent class labels via logistic gating on the
  ## standardised biomarker.
  bm_mean <- blparam[cat == 'bm']$m
  bm_sd   <- blparam[cat == 'bm']$sd
  bm_z    <- (dat$bm - bm_mean) / bm_sd
  logit_p <- modelparam$gating_intercept +
             modelparam$gating_slope * bm_z
  p_class1 <- plogis(logit_p)
  z <- rbinom(length(p_class1), size = 1, prob = p_class1)
  dat[, gating_prob := p_class1]
  dat[, latent_class := z]

  ## Stage 3. Apply class-specific BR scaling at on-drug timepoints.
  ## Class 1 (responder) multiplies BR by beta_R_factor; class 0
  ## (non-responder) by beta_NR_factor. At off-drug timepoints, BR
  ## is left as drawn from the MVN.
  d <- data.table(trialdesign)
  d[, onDrug := (tod > 0)]
  tnames <- trialdesign$timeptnames

  for (tp in seq_along(tnames)) {
    if (d[tp]$onDrug) {
      br_col <- paste(tnames[tp], 'br', sep = '.')
      scale_vec <- ifelse(z == 1, modelparam$beta_R_factor,
                                  modelparam$beta_NR_factor)
      dat[, (br_col) := get(br_col) * scale_vec]
    }
  }

  ## Stage 4. PB-class correlation extension. When
  ## pb_class_correlation = 'correlated', shift each participant's
  ## PB component upward (responders) or downward (non-responders)
  ## by an amount calibrated to deliver a target on-drug
  ## Cor(BR, PB).
  ##
  ## The shift is computed analytically rather than via grid search:
  ## given the marginal SDs of BR and PB at on-drug timepoints, the
  ## class-membership channel that already shifts BR (via
  ## beta_R_factor / beta_NR_factor) introduces a Cor(BR, class)
  ## of magnitude approximately
  ##     rho_BR_class = (beta_R_factor - beta_NR_factor) * sqrt(p (1-p)) / sd_BR
  ## where p is the marginal class prevalence. To achieve a target
  ## Cor(BR, PB) = rho_target via class membership alone, we shift
  ## PB by class with a magnitude that contributes the matching
  ## rho_PB_class. Concretely, the shift on PB for the responder
  ## class is computed so that
  ##     rho_PB_class = rho_target / rho_BR_class
  ## (assuming both correlations are positive). When the prevalence
  ## p departs from 0.5 the formula adjusts accordingly. The formula
  ## is approximate; the calibration sub-study refines it
  ## empirically for production runs.
  pb_class_shift <- 0
  if (modelparam$pb_class_correlation == 'correlated' &&
      modelparam$pb_class_target != 0) {
    pb_sd <- respparam[cat == 'pb']$sd
    if (length(pb_sd) == 0 || is.na(pb_sd)) pb_sd <- 2
    p_mean <- mean(p_class1)
    factor_diff <- modelparam$beta_R_factor - modelparam$beta_NR_factor
    sd_BR_marginal <- respparam[cat == 'br']$sd
    rho_BR_class <- abs(factor_diff) * sqrt(p_mean * (1 - p_mean)) /
                    pmax(sd_BR_marginal, 1e-6)
    rho_target <- modelparam$pb_class_target
    rho_PB_class <- rho_target / pmax(rho_BR_class, 1e-6)
    pb_class_shift <- rho_PB_class * pb_sd /
                      sqrt(p_mean * (1 - p_mean))
    pb_class_shift <- pmin(abs(pb_class_shift),
                           3 * pb_sd) * sign(pb_class_shift)

    ## Apply the shift to PB at on-drug timepoints.
    for (tp in seq_along(tnames)) {
      if (d[tp]$onDrug) {
        pb_col <- paste(tnames[tp], 'pb', sep = '.')
        shift_vec <- ifelse(z == 1, pb_class_shift, -pb_class_shift)
        dat[, (pb_col) := get(pb_col) + shift_vec]
      }
    }
  }
  dat[, pb_class_shift := pb_class_shift]

  ## Stage 5. Recompute total drug-effect columns and symptom-score
  ## columns (BL - sum of components).
  cl <- c('br', 'tv', 'pb')
  for (tp in seq_along(tnames)) {
    n <- tnames[tp]
    comps <- paste(n, cl, sep = '.')
    dat[, (paste0('D_', n)) := rowSums(.SD), .SDcols = comps]
    dat[, (n) := BL - rowSums(.SD), .SDcols = comps]
  }

  dat
}
