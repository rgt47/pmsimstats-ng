## lcmm_analysis: class-aware companion to pmsimstatsng::lme_analysis
##
## Wraps hlme to provide a two-class joint latent-class
## linear mixed model for the biomarker-treatment interaction
## inference question. Input/output contract matches
## lme_analysis() so that downstream simulation summaries can
## consume either analysis without branching.

## Assumes the host script has already loaded pmsimstats (via
## devtools::load_all('.') from the package root, or library()).
library(data.table)
library(nlme)
library(lcmm)

##' Class-aware analysis wrapper.
##'
##' @param trialdesign_set As accepted by lme_analysis().
##' @param dat As produced by generateData().
##' @param op Options list, with all lme_analysis fields plus:
##'   \itemize{
##'     \item ng        Number of latent classes (default 2).
##'     \item gating    Right-hand-side formula for class-membership
##'                     submodel (default ~bm).
##'     \item mixture   One-sided formula listing fixed-effects terms
##'                     whose coefficients differ across classes
##'                     (default ~Dbc, i.e. responder/non-responder
##'                     differ in main drug effect; bm:Dbc held
##'                     class-invariant so a single moderation
##'                     coefficient is identified for comparison
##'                     with the linear-mixed baseline).
##'     \item nproc     Number of processes for hlme grid search
##'                     (default 1).
##'     \item B         Initial-value matrix for hlme. If NULL,
##'                     hlme runs its default grid-search start.
##'   }
##' @return data.table with columns matching lme_analysis() plus
##'   class-aware diagnostics. Three moderation-test channels are
##'   reported because no single test is uniformly appropriate at
##'   trial-relevant sample sizes:
##'   \itemize{
##'     \item beta, betaSE, p extract the gating-coefficient
##'           Wald z-test on the biomarker. This is interpretable
##'           when the data support a 2-class fit, but is known to
##'           reject at inflated rates under the null at small N
##'           because the conditional SE does not absorb labelling
##'           uncertainty (Bauer-Curran-Lubke-Muthén overextraction).
##'     \item beta_within, betaSE_within, p_within extract the
##'           within-class bm:Dbc coefficient (class-invariant when
##'           mixture=~Dbc), retained for inferential symmetry with
##'           the linear-mixed baseline.
##'     \item delta_bic = BIC(ng=1) - BIC(ng=2). Positive values
##'           favour the latent-class model; the practical class-
##'           selection rule is delta_bic > 0 (or > 6 for moderate
##'           evidence per Raftery 1995).
##'     \item lrt_stat = 2 * (logLik(ng=2) - logLik(ng=1)). Asymptotic
##'           reference distribution is chi-squared with df equal to
##'           the parameter-count difference, but the boundary
##'           nature of mixture models (gating-class probabilities
##'           on the simplex boundary under the null) means the
##'           true distribution is a chi-bar-squared mixture. The
##'           statistic is reported here; calibration is left to
##'           the simulation analysis.
##'     \item lrt_df parameter-count difference between the ng=2
##'           and ng=1 fits (used as the conventional reference df).
##'     \item bic_ng1, bic_ng2 BIC values from the ng=1 anchor and
##'           the ng=2 mixture fits.
##'     \item entropy 1 - mean(per-participant posterior entropy) /
##'           log(K), with 1 = perfect class separation and 0 =
##'           uniform posterior.
##'     \item issingular TRUE if any random-intercept variance
##'           component is near-zero.
##'     \item warning conv-status string from hlme.
##'     \item conv hlme convergence flag for the ng=2 fit
##'           (1 = converged).
##'     \item conv_ng1 convergence flag for the ng=1 anchor.
##'   }
lcmm_analysis <- function(trialdesign_set, dat, op = list()) {

  ## Apply lme_analysis defaults for shared fields.
  op_defaults <- list(
    useDE = TRUE, t_random_slope = FALSE, full_model_out = FALSE,
    carryover_t1half = 0, simplecarryover = FALSE,
    carryover_scalefactor = 1,
    ng = 2, gating = ~bm, mixture = ~Dbc, nproc = 1, B = NULL
  )
  for (nm in names(op_defaults)) {
    if (is.null(op[[nm]])) op[[nm]] <- op_defaults[[nm]]
  }

  ## Stage 1: build the long-format data.table with Dbc using the
  ## same logic as lme_analysis. Reuse lme_analysis with full_
  ## model_out = TRUE to obtain the merged data; this avoids
  ## duplicating the data-prep logic and ensures the two analyses
  ## see identical inputs.
  op_prep <- op
  op_prep$full_model_out <- TRUE
  prep <- lme_analysis(trialdesign_set, dat, op_prep)
  datamerged <- prep$datamerged

  if (is.null(datamerged) || nrow(datamerged) == 0) {
    return(data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      issingular = NA, warning = "empty datamerged",
      entropy = NA_real_, bic = NA_real_, conv = NA_integer_
    ))
  }

  ## hlme requires numeric subject IDs; ptID already is.
  ## Drop rows with NA in any model variable.
  model_vars <- c("Sx", "bm", "t", "Dbc", "ptID")
  for (mv in model_vars) {
    datamerged <- datamerged[!is.na(get(mv))]
  }

  ## hlme fixed-effects formula matches the lme baseline.
  fixed_form <- Sx ~ bm + t + Dbc + bm:Dbc
  random_form <- ~ 1
  subject_var <- "ptID"

  ## Class-membership submodel via the `classmb` argument.
  classmb_form <- op$gating

  ## Fit. The ng=1 anchor is always fitted (it provides starting
  ## values for the ng=2 gridsearch and is itself reportable as
  ## the BIC reference for class-enumeration model comparison).
  df_in <- as.data.frame(datamerged)
  m1 <- tryCatch(
    hlme(
      fixed = fixed_form, random = random_form,
      subject = subject_var, ng = 1,
      data = df_in, verbose = FALSE, maxiter = 100
    ),
    error = function(e) {
      structure(list(error = conditionMessage(e), conv = -1L),
                class = "lcmm_failure")
    }
  )

  if (op$ng == 1) {
    fit <- m1
  } else {
    fit <- tryCatch({
      if (is.null(op$B)) {
        gridsearch(
          m = hlme(
            fixed = fixed_form, mixture = op$mixture,
            random = random_form,
            subject = subject_var, ng = op$ng,
            classmb = classmb_form,
            data = df_in, verbose = FALSE, maxiter = 30
          ),
          rep = 5, maxiter = 30, minit = m1
        )
      } else {
        hlme(
          fixed = fixed_form, mixture = op$mixture,
          random = random_form,
          subject = subject_var, ng = op$ng,
          classmb = classmb_form, B = op$B,
          data = df_in, verbose = FALSE, maxiter = 100
        )
      }
    }, error = function(e) {
      structure(list(error = conditionMessage(e), conv = -1L),
                class = "lcmm_failure")
    })
  }

  m1_ok <- !is.null(m1) && !inherits(m1, "lcmm_failure") &&
           isTRUE(m1$conv == 1)

  if (is.null(fit) || inherits(fit, "lcmm_failure") ||
      isTRUE(fit$conv != 1)) {
    msg <- if (inherits(fit, "lcmm_failure")) fit$error
           else if (is.null(fit)) "NULL"
           else paste0("conv=", fit$conv)
    return(data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      beta_within = NA_real_, betaSE_within = NA_real_,
      p_within = NA_real_,
      delta_bic = NA_real_, lrt_stat = NA_real_,
      lrt_df = NA_integer_,
      bic_ng1 = if (m1_ok) m1$BIC else NA_real_,
      bic_ng2 = NA_real_,
      issingular = NA, warning = paste0("hlme ", msg),
      entropy = NA_real_, bic = NA_real_,
      conv = if (is.null(fit) || inherits(fit, "lcmm_failure")) 0L
             else fit$conv,
      conv_ng1 = if (m1_ok) m1$conv else 0L
    ))
  }

  ## Build named SE vector. hlme stores fit$V as a vech-packed
  ## numeric vector of length n*(n+1)/2 over the upper triangle in
  ## column-major order. Unpack to a dense symmetric matrix and
  ## pull the diagonal.
  est <- fit$best
  n_par <- length(est)
  V_raw <- fit$V
  if (is.matrix(V_raw) &&
      nrow(V_raw) == n_par && ncol(V_raw) == n_par) {
    V_full <- V_raw
  } else {
    V_full <- matrix(0, n_par, n_par)
    V_full[upper.tri(V_full, diag = TRUE)] <- as.numeric(V_raw)
    V_full <- V_full + t(V_full) - diag(diag(V_full))
  }
  rownames(V_full) <- colnames(V_full) <- names(est)
  ses <- suppressWarnings(sqrt(diag(V_full)))

  z_p <- function(b, s) {
    if (!is.finite(b) || !is.finite(s) || s <= 0) return(c(NA_real_, NA_real_))
    z <- b / s
    c(z, 2 * pnorm(-abs(z)))
  }

  ## Primary class-aware moderation test: gating-coefficient z-test
  ## on the biomarker. With ng=2 and gating ~ bm, this is the
  ## 'bm class1' parameter (log-odds of class 1 vs the reference
  ## class per unit increase in bm).
  gating_target <- intersect(c("bm class1", "bm class2"), names(est))
  if (length(gating_target) == 0) {
    beta <- NA_real_; betaSE <- NA_real_; p <- NA_real_
  } else {
    gating_target <- gating_target[1]
    beta <- as.numeric(est[gating_target])
    betaSE <- as.numeric(ses[gating_target])
    zp <- z_p(beta, betaSE)
    p <- zp[2]
  }

  ## Secondary moderation signal: within-class bm:Dbc coefficient.
  ## Class-invariant when mixture = ~Dbc; reported alongside the
  ## gating test for symmetry with the lme_analysis baseline.
  within_target <- intersect(c("bm:Dbc", "Dbc:bm"), names(est))
  if (length(within_target) == 0) {
    beta_within <- NA_real_; betaSE_within <- NA_real_
    p_within <- NA_real_
  } else {
    within_target <- within_target[1]
    beta_within <- as.numeric(est[within_target])
    betaSE_within <- as.numeric(ses[within_target])
    zp_w <- z_p(beta_within, betaSE_within)
    p_within <- zp_w[2]
  }

  ## Posterior class-membership entropy (1 = perfectly separated).
  pp <- fit$pprob
  prob_cols <- grep("^prob", names(pp), value = TRUE)
  P <- as.matrix(pp[, prob_cols, drop = FALSE])
  H_i <- -rowSums(ifelse(P > 0, P * log(P), 0))
  H_norm <- mean(H_i) / log(op$ng)
  entropy <- 1 - H_norm

  ## Singularity: random-intercept variance near zero (lcmm uses
  ## Cholesky parameterisation; the diagonal entry is 'cholesky 1').
  chol_idx <- grep("^cholesky", names(est))
  issingular <- if (length(chol_idx) > 0) {
    isTRUE(any(abs(est[chol_idx]) < 1e-4, na.rm = TRUE))
  } else FALSE

  ## Class-enumeration model comparison statistics.
  if (m1_ok && op$ng > 1) {
    bic_ng1 <- as.numeric(m1$BIC)
    bic_ng2 <- as.numeric(fit$BIC)
    delta_bic <- bic_ng1 - bic_ng2
    lrt_stat <- 2 * (as.numeric(fit$loglik) - as.numeric(m1$loglik))
    lrt_df <- length(fit$best) - length(m1$best)
  } else {
    bic_ng1 <- NA_real_
    bic_ng2 <- as.numeric(fit$BIC)
    delta_bic <- NA_real_
    lrt_stat <- NA_real_
    lrt_df <- NA_integer_
  }

  data.table(
    beta = beta, betaSE = betaSE, p = p,
    beta_within = beta_within, betaSE_within = betaSE_within,
    p_within = p_within,
    delta_bic = delta_bic, lrt_stat = lrt_stat, lrt_df = lrt_df,
    bic_ng1 = bic_ng1, bic_ng2 = bic_ng2,
    issingular = issingular,
    warning = "",
    entropy = entropy,
    bic = fit$BIC,
    conv = fit$conv,
    conv_ng1 = if (m1_ok) m1$conv else 0L
  )
}
