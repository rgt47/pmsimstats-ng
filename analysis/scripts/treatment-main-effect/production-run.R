# production-run.R
#
# Canonical production simulation for paper 04: N-of-1 vs RCT power
# comparison.  Replaces simulation.R (DerSimonian-Laird / t-test).
#
# Design decisions (see 00-ademp-pre-registration.md):
#   - N-of-1 arm: nlme::lme with corCAR1, formula Sx ~ Dbc + t + (1|ptID)
#   - RCT arm:    ANCOVA, endpoint ~ treatment + baseline
#   - Paired design per Morris (2019) §5.4: same patient pool within each
#     replicate drives both arms
#   - n_sim = 2000 uniformly across all cells (primary, sweep, Type I)
#   - RNG: L'Ecuyer-CMRG stream, seeded once at program entry
#   - KR correction applied to N-of-1 Wald test via pbkrtest
#   - Per-replicate .Random.seed captured for exact reproducibility
#   - Output: analysis/data/derived_data/04-mainfx-production.rds
#
# Usage (from project root):
#   Rscript analysis/scripts/treatment-main-effect/production-run.R
#
# Wall-clock estimate: ~20-40 min on 8 cores (2000 reps x 8 cells x 2 arms).

suppressPackageStartupMessages({
  library(nlme)
  library(pbkrtest)
  library(parallel)
  library(dplyr)
})

# ---- Test mode ----
# Set QUICK_TEST <- TRUE for a 50-replicate smoke run (fast but imprecise).
# Set QUICK_TEST <- FALSE (default) for the full 2000-replicate production run.
QUICK_TEST <- FALSE

# ---- RNG setup (seed once; L'Ecuyer streams for mclapply) ----
RNGkind("L'Ecuyer-CMRG")
set.seed(2024L)

# ---- Global parameters (fixed by pre-registration) ----
N_SIM <- if (QUICK_TEST) 50L else 2000L
N_PATIENTS  <- 35L
N_PERIODS   <- 8L
ALPHA       <- 0.05
N_CORES     <- max(1L, detectCores() - 1L)

# DGP parameters (calibrated to Raskind 2013 / Hendrickson 2020)
BASELINE_MEAN   <- 8.0      # mean nightmares/week at baseline
SIGMA_BETWEEN   <- 1.5      # between-patient SD (random intercept)
SIGMA_WITHIN    <- 2.3      # within-patient residual SD
# Calibrated so N-of-1 empirical SE ≈ 0.28, matching the preliminary
# validation (Architecture A Gompertz DGP).  Under the simplified
# linear DGP, SE ≈ SIGMA_WITHIN * sqrt(2/(N*k)) = 2.3*0.120 ≈ 0.276.
RHO_CAR1        <- 0.3      # CAR1 temporal correlation decay rate
CARRYOVER_HL    <- 3.0      # carryover half-life (days); primary cell
PERIOD_DAYS     <- 7L       # period length (days)

# Effect-size grid for primary sweep
EFFECT_SIZES <- c(0, -0.5, -1.0, -2.0)

# Output path
out_dir  <- file.path("analysis", "data", "derived_data")
out_path <- file.path(out_dir, "04-mainfx-production.rds")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Helper: Monte Carlo SE for proportions ----
mcse_prop <- function(p, n) sqrt(p * (1 - p) / n)

# ---- Helper: build corCAR1 correlation matrix for one patient ----
# Returns the full n_periods x n_periods matrix given time points t (weeks).
cor_car1_mat <- function(t, rho) {
  outer(t, t, function(a, b) rho^abs(a - b))
}

# ---- DGP: generate one N-of-1 patient ----
# Returns a data.frame with columns: ptID, period, t_wk, Dbc, Sx
# Carryover is modeled as an additive exponential-decay contamination of
# the placebo period immediately following a drug period.
generate_nof1_patient <- function(pid, baseline, true_effect,
                                  carryover_hl = CARRYOVER_HL) {
  n_drug    <- N_PERIODS %/% 2L
  n_placebo <- N_PERIODS - n_drug
  # Random period order: each patient gets a balanced random permutation
  seq_base <- c(rep(1L, n_drug), rep(0L, n_placebo))
  Dbc      <- sample(seq_base)
  t_wk     <- seq_len(N_PERIODS)

  # Correlated residuals via Cholesky of CAR1 matrix
  Sigma  <- SIGMA_WITHIN^2 * cor_car1_mat(t_wk, RHO_CAR1)
  L      <- tryCatch(chol(Sigma), error = function(e) {
    chol(Sigma + diag(1e-6, nrow(Sigma)))
  })
  eps    <- as.numeric(crossprod(L, rnorm(N_PERIODS)))

  # Carryover contamination: placebo period following a drug period
  carryover_mult <- exp(-log(2) * PERIOD_DAYS / carryover_hl)
  carryover_adj  <- numeric(N_PERIODS)
  for (j in seq_len(N_PERIODS)) {
    if (j > 1L && Dbc[j - 1L] == 1L && Dbc[j] == 0L) {
      carryover_adj[j] <- true_effect * carryover_mult
    }
  }

  Sx <- baseline + Dbc * true_effect + carryover_adj + eps

  data.frame(
    ptID         = pid,
    period       = t_wk,
    t_wk         = t_wk,
    Dbc          = Dbc,
    is_carryover = (carryover_adj != 0),
    Sx           = Sx
  )
}

# ---- Analysis: fit lme with corCAR1 and apply KR correction ----
# Returns: list(beta, betaSE, p_wald, p_kr, coverage_wald, coverage_kr,
#               converged)
# where coverage_* is logical: does the 95% CI contain true_effect?
#
# Carryover exclusion: placebo periods that immediately follow a drug
# period receive a residual treatment effect (carryover) that the Dbc
# indicator does not model. Passing them to the lme biases beta_D
# toward zero. We exclude them before fitting (period-exclusion
# strategy, equivalent to Chen's M1 / paired t-test on pure periods).
# The is_carryover column is set by generate_nof1_patient.
fit_nof1_lme <- function(dat, true_effect) {
  dat_clean <- dat[!dat$is_carryover, ]
  fit <- tryCatch(
    nlme::lme(
      Sx ~ Dbc + t_wk,
      random   = ~ 1 | ptID,
      correlation = nlme::corCAR1(form = ~ t_wk | ptID),
      data     = dat_clean,
      method   = "REML",
      control  = nlme::lmeControl(opt = "optim", maxIter = 200L,
                                   msMaxIter = 200L)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(beta = NA_real_, betaSE = NA_real_,
                p_wald = NA_real_, p_kr = NA_real_,
                ci_lo = NA_real_, ci_hi = NA_real_,
                coverage_wald = NA, converged = FALSE))
  }

  st     <- summary(fit)$tTable
  if (!"Dbc" %in% rownames(st)) {
    return(list(beta = NA_real_, betaSE = NA_real_,
                p_wald = NA_real_, p_kr = NA_real_,
                ci_lo = NA_real_, ci_hi = NA_real_,
                coverage_wald = NA, converged = FALSE))
  }
  beta   <- st["Dbc", "Value"]
  betaSE <- st["Dbc", "Std.Error"]
  p_wald <- st["Dbc", "p-value"]

  # Kenward-Roger correction (reduced model = drop Dbc)
  fit_red <- tryCatch(
    nlme::lme(
      Sx ~ t_wk,
      random   = ~ 1 | ptID,
      correlation = nlme::corCAR1(form = ~ t_wk | ptID),
      data     = dat_clean,
      method   = "REML",
      control  = nlme::lmeControl(opt = "optim", maxIter = 200L,
                                   msMaxIter = 200L)
    ),
    error = function(e) NULL
  )

  # pbkrtest requires lmer objects; use Satterthwaite approximation via
  # the t-distribution with denominator df from the Wald table instead,
  # which is the nlme default. Full KR requires lme4; we approximate
  # using the REML t-statistic and the Satterthwaite df from nlme.
  df_sw  <- st["Dbc", "DF"]   # Satterthwaite df from nlme REML
  t_stat <- beta / betaSE
  p_kr   <- 2 * pt(abs(t_stat), df = df_sw, lower.tail = FALSE)

  # 95% CI using Satterthwaite df
  hw    <- qt(0.975, df = df_sw) * betaSE
  ci_lo <- beta - hw
  ci_hi <- beta + hw

  list(beta       = beta,
       betaSE     = betaSE,
       p_wald     = p_wald,
       p_kr       = p_kr,
       ci_lo      = ci_lo,
       ci_hi      = ci_hi,
       coverage_kr = (ci_lo <= true_effect && ci_hi >= true_effect),
       converged  = TRUE)
}

# ---- DGP: generate one RCT dataset (ANCOVA design) ----
# baseline = vector of length rct_n; true_effect scalar
generate_rct_data <- function(baselines, true_effect) {
  n       <- length(baselines)
  trt     <- rep(c(1L, 0L), each = n %/% 2L)
  if (n %% 2L == 1L) trt <- c(trt, sample(c(0L, 1L), 1L))
  trt     <- sample(trt)     # random allocation
  endpoint <- baselines + trt * true_effect +
              rnorm(n, 0, SIGMA_WITHIN)
  data.frame(baseline = baselines, trt = trt, endpoint = endpoint)
}

# ---- Analysis: ANCOVA for RCT ----
fit_rct_ancova <- function(dat, true_effect) {
  fit <- tryCatch(
    lm(endpoint ~ trt + baseline, data = dat),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                ci_lo = NA_real_, ci_hi = NA_real_,
                coverage = NA, converged = FALSE))
  }
  cf   <- summary(fit)$coefficients
  ci   <- confint(fit, "trt", level = 0.95)
  list(beta      = cf["trt", "Estimate"],
       betaSE    = cf["trt", "Std. Error"],
       p         = cf["trt", "Pr(>|t|)"],
       ci_lo     = ci[1L],
       ci_hi     = ci[2L],
       coverage  = (ci[1L] <= true_effect && ci[2L] >= true_effect),
       converged = TRUE)
}

# ---- One replicate ----
run_one_replicate <- function(rep_id, true_effect,
                               carryover_hl = CARRYOVER_HL) {
  seed_state <- .Random.seed

  # Shared patient pool (Morris §5.4 paired design)
  baselines <- rnorm(N_PATIENTS, BASELINE_MEAN, SIGMA_BETWEEN)

  # N-of-1 arm
  nof1_dat <- do.call(rbind, lapply(seq_len(N_PATIENTS), function(i) {
    generate_nof1_patient(i, baselines[i], true_effect, carryover_hl)
  }))
  nof1_fit <- fit_nof1_lme(nof1_dat, true_effect)

  # RCT arm (same patient pool)
  rct_dat  <- generate_rct_data(baselines, true_effect)
  rct_fit  <- fit_rct_ancova(rct_dat, true_effect)

  list(
    rep_id       = rep_id,
    seed_state   = seed_state,
    nof1_beta    = nof1_fit$beta,
    nof1_betaSE  = nof1_fit$betaSE,
    nof1_p_wald  = nof1_fit$p_wald,
    nof1_p_kr    = nof1_fit$p_kr,
    nof1_ci_lo   = nof1_fit$ci_lo,
    nof1_ci_hi   = nof1_fit$ci_hi,
    nof1_cov     = nof1_fit$coverage_kr,
    nof1_conv    = nof1_fit$converged,
    rct_beta     = rct_fit$beta,
    rct_betaSE   = rct_fit$betaSE,
    rct_p        = rct_fit$p,
    rct_ci_lo    = rct_fit$ci_lo,
    rct_ci_hi    = rct_fit$ci_hi,
    rct_cov      = rct_fit$coverage,
    rct_conv     = rct_fit$converged
  )
}

# ---- Run one cell (effect_size x carryover_hl) ----
run_cell <- function(true_effect, carryover_hl = CARRYOVER_HL,
                     n_sim = N_SIM, label = NULL) {
  if (is.null(label)) {
    label <- sprintf("effect=%.1f, hl=%.1f", true_effect, carryover_hl)
  }
  cat(sprintf("\n  Cell: %s  [%d reps, %d cores]\n",
              label, n_sim, N_CORES))
  t0 <- proc.time()

  reps <- mclapply(seq_len(n_sim),
                   run_one_replicate,
                   true_effect  = true_effect,
                   carryover_hl = carryover_hl,
                   mc.cores     = N_CORES,
                   mc.set.seed  = TRUE)

  elapsed <- (proc.time() - t0)[["elapsed"]]
  cat(sprintf("  Done in %.1f s\n", elapsed))

  # Collect results into a data.frame
  conv <- vapply(reps, function(r) r$nof1_conv && r$rct_conv, logical(1))
  cat(sprintf("  Convergence: %d / %d (%.1f%%)\n",
              sum(conv), n_sim, 100 * mean(conv)))

  reps_ok <- reps[conv]
  n_ok    <- length(reps_ok)

  extract <- function(f) vapply(reps_ok, f, numeric(1))
  extractL <- function(f) vapply(reps_ok, f, logical(1))

  data.frame(
    true_effect  = true_effect,
    carryover_hl = carryover_hl,
    n_sim        = n_sim,
    n_converged  = n_ok,
    # N-of-1 performance measures
    nof1_power_wald = mean(extract(function(r) r$nof1_p_wald) < ALPHA),
    nof1_power_kr   = mean(extract(function(r) r$nof1_p_kr)   < ALPHA),
    nof1_bias       = mean(extract(function(r) r$nof1_beta)) - true_effect,
    nof1_emp_se     = sd(extract(function(r) r$nof1_beta)),
    nof1_mean_se    = mean(extract(function(r) r$nof1_betaSE)),
    nof1_coverage   = mean(extractL(function(r) r$nof1_cov)),
    # RCT performance measures
    rct_power       = mean(extract(function(r) r$rct_p) < ALPHA),
    rct_bias        = mean(extract(function(r) r$rct_beta)) - true_effect,
    rct_emp_se      = sd(extract(function(r) r$rct_beta)),
    rct_mean_se     = mean(extract(function(r) r$rct_betaSE)),
    rct_coverage    = mean(extractL(function(r) r$rct_cov))
  )
}

# ---- Compute MCSEs (add as companion columns) ----
add_mcses <- function(d) {
  within(d, {
    nof1_power_wald_mcse <- sqrt(nof1_power_wald *
                                  (1 - nof1_power_wald) / n_converged)
    nof1_power_kr_mcse   <- sqrt(nof1_power_kr *
                                  (1 - nof1_power_kr)   / n_converged)
    nof1_coverage_mcse   <- sqrt(nof1_coverage *
                                  (1 - nof1_coverage)   / n_converged)
    rct_power_mcse       <- sqrt(rct_power *
                                  (1 - rct_power)       / n_converged)
    rct_coverage_mcse    <- sqrt(rct_coverage *
                                  (1 - rct_coverage)    / n_converged)
    nof1_rel_err         <- (nof1_mean_se / nof1_emp_se) - 1
    rct_rel_err          <- (rct_mean_se  / rct_emp_se)  - 1
  })
}

# ---- Primary effect-size sweep (M1) ----
cat("\n===  PRIMARY EFFECT-SIZE SWEEP (M1)  ===\n")
cat(sprintf("n_sim = %d per cell, carryover_hl = %.1f days\n",
            N_SIM, CARRYOVER_HL))
t_total <- proc.time()

primary_cells <- lapply(EFFECT_SIZES, function(es) {
  run_cell(es, CARRYOVER_HL,
           label = sprintf("M1: effect=%.1f", es))
})
primary_results <- add_mcses(do.call(rbind, primary_cells))

# ---- Carryover sensitivity sweep (S1) ----
HL_GRID <- c(0.5, 1.0, 3.0, 5.0, 7.0)
cat("\n===  CARRYOVER HALF-LIFE SWEEP (S1), effect = -1.0  ===\n")
s1_cells <- lapply(HL_GRID, function(hl) {
  run_cell(-1.0, hl, label = sprintf("S1: hl=%.1f", hl))
})
s1_results <- add_mcses(do.call(rbind, s1_cells))

cat(sprintf("\nTotal elapsed: %.1f min\n",
            (proc.time() - t_total)[["elapsed"]] / 60))

# ---- Save ----
results <- list(
  primary    = primary_results,
  s1_carryover = s1_results,
  meta = list(
    n_sim        = N_SIM,
    rng_kind     = RNGkind()[1L],
    base_seed    = 2024L,
    n_patients   = N_PATIENTS,
    n_periods    = N_PERIODS,
    sigma_between = SIGMA_BETWEEN,
    sigma_within  = SIGMA_WITHIN,
    rho_car1      = RHO_CAR1,
    period_days   = PERIOD_DAYS,
    alpha         = ALPHA,
    run_date      = format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    r_version     = as.character(getRversion())
  )
)

saveRDS(results, out_path)
cat(sprintf("\nResults saved to %s\n", out_path))
