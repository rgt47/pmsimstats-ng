## True heterogeneous-random-slopes DGP for Study 3 (Verbeke-
## Lesaffre / Proust-Lima parameterisation).
##
## Phase A infrastructure for paper 03 path (b). Implements the
## parameter-level mixture per the ADEMP pre-registration
## (00-ademp-pre-registration.md, §"Study 3"):
##
##   beta_R   ~ N(mu_R,  sigma_R^2)    (responder class)
##   beta_NR  ~ N(0,     sigma_NR^2)   (non-responder class)
##   pi(B_i) = logit^{-1}(alpha + beta * standardised(B_i))
##
## Each participant's drug-effect random slope on D_bc is drawn
## from a class-conditional normal distribution; class membership
## is Bernoulli with logistic gating on the standardised
## biomarker. The marginal random-slope distribution is a
## two-component normal mixture: this is the formal mixture-of-
## normals at the parameter level, distinct from the
## marginal-mixture DGP at 05-marginal-mixture-dgp.R (which
## applies a deterministic class-specific BR-factor scaling).
##
## Note: the existing file at 02-heterogeneous-slopes-dgp.R
## implements the deterministic class-specific scaling (the
## marginal-mixture form) under the file name "heterogeneous-
## slopes". This 07- file is the parameter-level mixture proper,
## per the pre-reg's distinction between Study 1 (outcome-level
## mixture) and Study 3 (parameter-level mixture).

library(data.table)

##' Generate data under the Study 3 heterogeneous-random-slopes
##' DGP.
##'
##' @param modelparam Same one-row data.table as generateData(),
##'   with these extra fields (NULL-safe defaults applied):
##'   \itemize{
##'     \item mu_R Mean of the responder-class random slope.
##'           Default 1.0.
##'     \item sigma_R Standard deviation of responder-class
##'           random slope. Default 0.3.
##'     \item mu_NR Mean of the non-responder-class random slope.
##'           Default 0.0.
##'     \item sigma_NR Standard deviation of non-responder-class
##'           random slope. Default 0.3.
##'     \item gating_intercept Logistic-gating intercept. Default
##'           0 (class prevalence 0.5 at the biomarker mean).
##'     \item gating_slope Logistic-gating slope on standardised
##'           biomarker. Default 1.0.
##'   }
##' @param respparam As for generateData().
##' @param blparam As for generateData().
##' @param trialdesign As for generateData().
##' @param empirical,makePositiveDefinite,seed,lambda_cor,verbose,
##'   cached_sigma As for generateData().
##' @return Same long-format data.table as generateData() with
##'   four extra columns: latent_class (0/1), gating_prob, and
##'   beta_slope (the realised random slope on D_bc per subject).
generate_true_random_slopes_data <- function(
  modelparam, respparam, blparam, trialdesign,
  empirical = FALSE, makePositiveDefinite = TRUE,
  seed = NA, lambda_cor = NA, verbose = FALSE,
  cached_sigma = NULL
) {

  defaults <- list(
    mu_R = 1.0, sigma_R = 0.3,
    mu_NR = 0.0, sigma_NR = 0.3,
    gating_intercept = 0.0, gating_slope = 1.0
  )
  for (nm in names(defaults)) {
    if (is.null(modelparam[[nm]])) {
      modelparam[[nm]] <- defaults[[nm]]
    }
  }

  ## Stage 1. Draw the MVN backbone with c.bm = 0 so the
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

  ## Stage 3. Draw class-conditional random slopes on D_bc.
  ## Responders draw from N(mu_R, sigma_R^2); non-responders
  ## draw from N(mu_NR, sigma_NR^2). The marginal distribution
  ## of beta_slope is the two-component normal mixture that
  ## defines the Study 3 DGP.
  beta_slope <- ifelse(
    z == 1,
    rnorm(length(z), mean = modelparam$mu_R, sd = modelparam$sigma_R),
    rnorm(length(z), mean = modelparam$mu_NR, sd = modelparam$sigma_NR)
  )
  dat[, beta_slope := beta_slope]

  ## Stage 4. Apply the random-slope scaling to BR at on-drug
  ## timepoints. Each subject's BR contribution at on-drug
  ## timepoints is multiplied by their drawn beta_slope, which
  ## under the responder-class mean of mu_R = 1.0 recovers the
  ## standard MVN draw (no scaling) and under the non-responder
  ## mean of mu_NR = 0.0 zeros out the BR contribution. The
  ## within-class sigma_R / sigma_NR contributes additional
  ## within-class variation that distinguishes this DGP from
  ## the deterministic-scaling form at 02-* / 05-*.
  d <- data.table(trialdesign)
  d[, onDrug := (tod > 0)]
  tnames <- trialdesign$timeptnames

  for (tp in seq_along(tnames)) {
    if (d[tp]$onDrug) {
      br_col <- paste(tnames[tp], 'br', sep = '.')
      dat[, (br_col) := get(br_col) * beta_slope]
    }
  }

  ## Stage 5. Recompute total drug-effect columns and symptom-
  ## score columns (BL - sum of components).
  cl <- c('br', 'tv', 'pb')
  for (tp in seq_along(tnames)) {
    n <- tnames[tp]
    comps <- paste(n, cl, sep = '.')
    dat[, (paste0('D_', n)) := rowSums(.SD), .SDcols = comps]
    dat[, (n) := BL - rowSums(.SD), .SDcols = comps]
  }

  dat
}
