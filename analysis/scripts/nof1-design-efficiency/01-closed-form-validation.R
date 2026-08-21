## analysis/scripts/nof1-design-efficiency/01-closed-form-validation.R
##
## Derives a closed-form (AR(1) sample-mean-variance) approximation to
## the sampling variance and power of the biomarker-treatment
## interaction coefficient bm:Dbc under the mean-moderation
## architecture, fit under analysis specification E1 (binary drug
## indicator, carryover unmodeled) reduced to a two-group
## within-subject contrast. Validates the approximation against the
## real Monte Carlo output already produced by the dual-channel driver
## (analysis/scripts/dgp-combined/01-run-combined-factorial.R), whose
## (c_bm_a = 0.45, c_bm_b = 0) boundary is exactly the pure
## mean-moderation architecture.
##
## No new simulation is run here: this script only reads
## analysis/data/dgp-combined/01-combined-summary.rds (already on
## disk) and the trial-design schedules hard-coded in
## analysis/scripts/carryover-sensitivity/simulation-core.R
## (design_preset()), reproduced here rather than sourced, since only
## the schedule constants are needed, not the simulation machinery.
##
## Outputs:
##   analysis/data/nof1-design-efficiency/closed-form-validation.rds
##   analysis/data/nof1-design-efficiency/closed-form-validation.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

repo_root <- here::here()

## -----------------------------------------------------------------
## Trial-design schedules (copied verbatim from design_preset() in
## analysis/scripts/carryover-sensitivity/simulation-core.R, to avoid
## sourcing that file's heavier dependencies for a read-only check).
## -----------------------------------------------------------------

design_preset <- function(name) {
  switch(name,
    CO = list(
      timepoints = cumsum(rep(2.5, 8)),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
                    pathB = c(0, 0, 0, 0, 1, 1, 1, 1))
    ),
    OLBDC = list(
      timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
      ondrug = list(pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
                    pathB = c(1, 1, 1, 1, 1, 0, 0, 0))
    ),
    Hybrid = list(
      timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
                    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
                    pathD = c(1, 1, 1, 0, 0, 0, 0, 1))
    )
  )
}

## Time since the most recent on-drug timepoint, for off-drug visits;
## Inf (never yet on drug) for an off-drug visit preceding any on-drug
## exposure on that path (e.g. CO's placebo-then-drug arm).
compute_tsd <- function(timepoints, ondrug) {
  tsd <- numeric(length(timepoints))
  last_on_time <- NA_real_
  for (i in seq_along(timepoints)) {
    if (ondrug[i] == 1) {
      tsd[i] <- 0
      last_on_time <- timepoints[i]
    } else {
      tsd[i] <- if (is.na(last_on_time)) Inf else timepoints[i] - last_on_time
    }
  }
  tsd
}

## Exact variance of the mean of k stationary AR(1)-correlated
## observations at arbitrary (possibly unequally spaced) times:
## Var(Xbar) = (sigma^2 / k^2) * sum_{j,l} rho^{|t_j - t_l|}.
var_mean_ar1 <- function(times, sigma2, rho) {
  k <- length(times)
  if (k == 0) return(NA_real_)
  d <- abs(outer(times, times, "-"))
  (sigma2 / k^2) * sum(rho ^ d)
}

## Per-design, per-half-life: attenuation factor A = 1 - mean(D_bc at
## off-drug visits), and the closed-form on+off contrast variance
## v_sum = v_on + v_off (sigma2 = 1; sigma2 enters later as a single
## free scale), averaged across randomization paths.
closed_form_terms <- function(design, t1half, rho) {
  spec <- design_preset(design)
  lambda <- if (t1half > 0) log(2) / t1half else Inf
  per_path <- map(spec$ondrug, function(od) {
    tp <- spec$timepoints
    on_idx  <- which(od == 1)
    off_idx <- which(od == 0)
    tsd <- compute_tsd(tp, od)
    dbc_off <- if (length(off_idx) == 0) numeric(0) else
      if (is.infinite(lambda)) rep(0, length(off_idx)) else
        exp(-lambda * tsd[off_idx])
    tibble(
      dbar_off = if (length(off_idx) == 0) NA_real_ else mean(dbc_off),
      v_on  = var_mean_ar1(tp[on_idx],  sigma2 = 1, rho = rho),
      v_off = if (length(off_idx) == 0) NA_real_ else
        var_mean_ar1(tp[off_idx], sigma2 = 1, rho = rho)
    )
  }) |> bind_rows()
  tibble(v_sum = mean(per_path$v_on + per_path$v_off),
         A = 1 - mean(per_path$dbar_off, na.rm = TRUE))
}

## -----------------------------------------------------------------
## Real Monte Carlo output: the pure mean-moderation boundary
## (c_bm_a = 0.45, c_bm_b = 0) of the dual-channel driver, analysis
## specification E1 (binary indicator).
## -----------------------------------------------------------------

s <- readRDS(file.path(repo_root,
  "analysis/data/dgp-combined/01-combined-summary.rds"))

d <- s |>
  filter(spec == "E1", c_bm_a == 0.45, c_bm_b == 0) |>
  select(design, N, t1half, power, mcse_power, sd_beta, mean_beta)

cf <- expand.grid(design = c("CO", "OLBDC", "Hybrid"),
                   t1half = c(0, 0.5, 1.0),
                   stringsAsFactors = FALSE) |>
  rowwise() |>
  mutate(terms = list(closed_form_terms(design, t1half, rho = 0.7))) |>
  unnest(terms) |>
  ungroup()

d2 <- d |> left_join(cf, by = c("design", "t1half"))

## Single free scale parameter sigma^2, fit by OLS through the origin
## on sd_beta^2 * N ~ v_sum (the closed form predicts this product is
## constant in N and design-and-half-life dependent only through
## v_sum).
fit <- lm(I(sd_beta^2 * N) ~ 0 + v_sum, data = d2)
sigma2_hat <- unname(coef(fit)[1])

d2 <- d2 |>
  mutate(
    pred_sd_beta = sqrt(sigma2_hat * v_sum / N),
    ratio_sd     = sd_beta / pred_sd_beta,
    z_pred       = abs(mean_beta) / pred_sd_beta,
    power_pred   = pnorm(z_pred - qnorm(0.975)) +
                    pnorm(-z_pred - qnorm(0.975)),
    power_err    = power_pred - power
  )

cat(sprintf("Fitted sigma^2 (variance-scale calibration): %.4f\n",
            sigma2_hat))
cat(sprintf("Variance-model R^2 (sd_beta^2 * N ~ v_sum): %.4f\n",
            summary(fit)$r.squared))
cat(sprintf("Power prediction MAE:  %.4f\n", mean(abs(d2$power_err))))
cat(sprintf("Power prediction RMSE: %.4f\n", sqrt(mean(d2$power_err^2))))
cat(sprintf("Power prediction correlation with simulated power: %.4f\n",
            cor(d2$power, d2$power_pred)))

out_dir <- file.path(repo_root, "analysis/data/nof1-design-efficiency")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(d2, file.path(out_dir, "closed-form-validation.rds"))
write.csv(d2, file.path(out_dir, "closed-form-validation.csv"),
          row.names = FALSE)
cat("Saved to", out_dir, "\n")
