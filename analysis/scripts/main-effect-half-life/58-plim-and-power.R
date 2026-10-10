## docs/58, closed-form results for the drug main effect beta_D under the
## refined process: estimand values, the GLS probability limit of the
## coefficient of D_bc(h) under three working covariances, its
## decomposition into a null bias and an exposure-weighted drug signal,
## sandwich standard errors, normal-approximation size and power, the
## dependence on h, t1/2, rho, tau and N, anchored against recursive
## carryover, the row-group decomposition behind the expectancy
## confounding, and the expected profile deviance in h.
##
## Run from the repository root:
##   Rscript analysis/scripts/main-effect-half-life/58-plim-and-power.R
## Writes analysis/data/quick-sim/main-effect-half-life/58-cf-*.csv.
source('analysis/scripts/main-effect-half-life/58-common.R')
load_driver()
t_start <- Sys.time()
DES <- c('Hybrid', 'OL+BDC', 'CO')
TH <- c(0, 0.1, 0.2, 0.5, 1, 2)
HH <- c(0, 0.1, 0.2, 0.25, 0.5, 1, 2)
TYPES <- c('S', 'F1', 'VF')

## ---- 1. estimands, anchored and recursive, f_m = 1 ----
est <- rbindlist(lapply(DES, function(des) rbindlist(lapply(TH, function(th)
  rbindlist(lapply(c('anchored', 'recursive'), function(cy) {
    pms <- design_moments(des, th, 1, carry = cy)
    data.table(design = des, t_half = th, carry = cy, full = estimand(pms, 1, 'full'),
               exp = estimand(pms, 1, 'exp'), bl = estimand(pms, 1, 'bl'))
  }))))))
fwrite(est, file.path(OUT, '58-cf-estimands.csv'))
print(est[carry == 'anchored'])

## Working-covariance limits are computed once per (design, t1/2, h,
## type, f_m, rho, tau); cache them.
rc_cache <- new.env()
plim_one <- function(des, th, h, type, f_m, rho = 0.7, tau = 0, cov = 'rc', N = 70,
                     carry = 'anchored') {
  pms <- design_moments(des, th, f_m, tau, rho, carry)
  rc <- NULL
  if (cov == 'rc') {
    key <- paste(des, th, h, type, f_m, rho, tau, carry, sep = '|')
    if (is.null(rc_cache[[key]])) {
      Xs <- lapply(pms, design_X, h = h, type = type)
      rc_cache[[key]] <- rc_plim(Xs, pms)
    }
    rc <- rc_cache[[key]]
  }
  g <- analysis_plim(pms, h, type, cov, N, rc)
  ## decomposition: null bias b0 (same V) and drug signal S(h; t1/2)
  pm0 <- lapply(pms, function(pm) { pm$mu <- pm$mu + f_m * M_BR * pm$q; pm })
  b0 <- gls_plim(g$Xs, g$V, pm0, N = N)$b
  sig <- mean(mapply(function(om, pm) sum(om * pm$q), g$omega, pms))
  ## variance pieces: tau^2 part
  data.table(design = des, t_half = th, h = h, type = type, f_m = f_m, rho = rho,
             tau = tau, cov = cov, N = N, carry = carry, plim = g$b, b0 = b0,
             signal = sig, se = sqrt(g$var), exp = estimand(pms, f_m, 'exp'),
             dbar = dbar(pms, h), rc_phi = if (!is.null(rc)) rc$phi else NA_real_,
             rc_su2 = if (!is.null(rc)) rc$su2 else NA_real_,
             rc_se2 = if (!is.null(rc)) rc$se2 else NA_real_)
}
plim_safe <- function(...) tryCatch(plim_one(...), error = function(e) NULL)
finish <- function(dt) {
  dt[, `:=`(plim_exp = plim * dbar, bias_exp = plim * dbar - exp,
            size_or_power = rej(plim, se), z = plim / se)]
  dt[, power_t := {
    df <- N - 4
    1 - pt(qt(0.975, df), df, ncp = abs(z)) + pt(-qt(0.975, df), df, ncp = abs(z))
  }]
  dt
}

## ---- 2. null limits by design, analysis, h and working covariance ----
null <- rbindlist(lapply(DES, function(des) rbindlist(lapply(TYPES, function(ty)
  rbindlist(lapply(c(0, 0.2, 0.5, 1, 2, Inf), function(h)
    rbindlist(lapply(c('indep', 'true', 'rc'), function(cv)
      plim_safe(des, 0.2, h, ty, 0, cov = cv))))))))) |> finish()
fwrite(null, file.path(OUT, '58-cf-null.csv'))
print(null[cov == 'rc', .(design, type, h, b0 = round(b0, 3), se = round(se, 3),
                          size = round(size_or_power, 4))])

## ---- 3. alternatives: plim, attenuation and power as h departs from t1/2 ----
alt <- rbindlist(lapply(DES, function(des) rbindlist(lapply(c('F1', 'VF', 'S'), function(ty)
  rbindlist(lapply(c(0, 0.2, 0.5, 1, 2), function(th) rbindlist(lapply(HH, function(h)
    plim_one(des, th, h, ty, 0.4))))))))) |> finish()
fwrite(alt, file.path(OUT, '58-cf-alt.csv'))
print(alt[design == 'Hybrid' & type == 'F1',
          .(t_half, h, plim = round(plim, 3), b0 = round(b0, 3), sig = round(signal, 3),
            plim_exp = round(plim_exp, 3), exp = round(exp, 3), se = round(se, 3),
            power = round(size_or_power, 3))])

## ---- 4. power by rho, tau and N at h = t1/2 = 0.2, f_m = 0.4 ----
grid4 <- CJ(design = DES, type = c('F1', 'VF'), rho = c(0.5, 0.7, 0.8),
            tau = c(0, 2.75, 5.49), N = c(35, 70, 150))
pw <- rbindlist(lapply(seq_len(nrow(grid4)), function(i) with(grid4[i],
  plim_one(design, 0.2, 0.2, type, 0.4, rho = rho, tau = tau, N = N)))) |> finish()
fwrite(pw, file.path(OUT, '58-cf-power-grid.csv'))
print(pw[design == 'Hybrid', .(type, rho, tau, N, se = round(se, 3),
                               power = round(size_or_power, 3))])
## heterogeneous (weak) null: f_m = 0, tau > 0
wn <- rbindlist(lapply(DES, function(des) rbindlist(lapply(c('F1', 'VF', 'S'), function(ty)
  rbindlist(lapply(c(0, 2.75, 5.49), function(tau)
    plim_one(des, 0.2, 0.2, ty, 0, tau = tau))))))) |> finish()
fwrite(wn, file.path(OUT, '58-cf-weak-null.csv'))

## Power curves over f_m for the figure (Hybrid, OL+BDC, CO; F1 and VF;
## h = t1/2 = 0.2; N = 70).
curve <- rbindlist(lapply(DES, function(des) rbindlist(lapply(c('F1', 'VF'), function(ty)
  rbindlist(lapply(seq(0, 1, 0.1), function(f) plim_one(des, 0.2, 0.2, ty, f))))))) |> finish()
fwrite(curve, file.path(OUT, '58-cf-power-curve.csv'))

## ---- 5. anchored against recursive carryover at h = t1/2 ----
ar <- rbindlist(lapply(DES, function(des) rbindlist(lapply(c('F1', 'VF'), function(ty)
  rbindlist(lapply(c(0.2, 0.5, 1, 2), function(th) rbindlist(lapply(c('anchored', 'recursive'),
    function(cy) plim_one(des, th, th, ty, 0.4, carry = cy))))))))) |> finish()
fwrite(ar, file.path(OUT, '58-cf-anchored-recursive.csv'))
print(dcast(ar, design + type + t_half ~ carry, value.var = c('plim', 'exp', 'size_or_power')))

## ---- 6. expectancy confounding: correlation and row-group decomposition ----
pms0 <- design_moments('Hybrid', 0.2, 0)
long <- rbindlist(lapply(pms0, function(pm) data.table(D = as.numeric(pm$on), e = pm$e,
                                                       X = dbc_h(pm$on, pm$tsd, 0.2))))
cors <- data.table(quantity = c('cor(D, eta), nine rows', 'cor(D_bc(0.2), eta), nine rows',
                                'cor(D, eta), post-baseline rows'),
                   value = c(cor(long$D, long$e), cor(long$X, long$e),
                             long[e > 0, cor(D, e)]))
print(cors)
fwrite(cors, file.path(OUT, '58-cf-expectancy-cor.csv'))
## Per-row weights of the drug coefficient and the contribution of each
## row group to the null bias b0 = sum_p w_p omega_p' mu0_p (the mean
## centered within path, which changes no total because the weights are
## orthogonal to the intercept) and to the
## signal; groups: baseline, open label (all on drug), blinded on,
## blinded off.
rowdec <- rbindlist(lapply(c('S', 'F1', 'VF'), function(ty) {
  g <- analysis_plim(pms0, 0.2, ty, 'rc',
                     rc = rc_plim(lapply(pms0, design_X, h = 0.2, type = ty), pms0))
  pms1 <- design_moments('Hybrid', 0.2, 1)
  rbindlist(Map(function(om, pm, pm1, k) {
    grp <- ifelse(seq_len(pm$n) == 1, 'baseline', ifelse(pm$e == 1, 'open label, on',
                  ifelse(pm$on, 'blinded, on', 'blinded, off')))
    data.table(type = ty, path = k, group = grp, w = om / length(pms0),
               bias = om * (pm$mu - mean(pm$mu)) / length(pms0), signal = -om * M_BR * pm1$q / length(pms0))
  }, g$omega, pms0, pms1, seq_along(pms0)))
}))
rowsum <- rowdec[, .(sum_abs_w = sum(abs(w)), sum_w = sum(w), bias = sum(bias),
                     signal = sum(signal)), by = .(type, group)]
rowsum[, share_abs_w := sum_abs_w / sum(sum_abs_w), by = type]
print(rowsum)
fwrite(rowsum, file.path(OUT, '58-cf-rowgroups.csv'))

## ---- 7. how much the data say about h: expected profile deviance ----
## Under the truth (t1/2, f_m), the expected ML deviance of F1 at h is
## about N * sum_p w_p [r_p(h)' V^-1 r_p(h)] minus its minimum over h,
## with V the working covariance limit at the true h (first order).
hfine <- c(0, 0.02, 0.05, seq(0.1, 1, 0.05), seq(1.1, 3, 0.1), 4, 6, 10, Inf)
dev <- rbindlist(lapply(c(0, 0.2, 0.5, 1, 2), function(th) rbindlist(lapply(c(0.4, 1), function(f) {
  pms <- design_moments('Hybrid', th, f)
  rc <- rc_plim(lapply(pms, design_X, h = max(th, 0), type = 'F1'), pms)
  Vi <- solve(rc$V)
  rss <- sapply(hfine, function(h) {
    Xs <- lapply(pms, design_X, h = h, type = 'F1')
    g <- gls_plim(Xs, rc$V, pms)
    mean(sapply(g$res, function(r) sum(r * (Vi %*% r))))
  })
  data.table(t_half = th, f_m = f, h = hfine, dev70 = 70 * (rss - min(rss)),
             h_star = hfine[which.min(rss)])
}))))
fwrite(dev, file.path(OUT, '58-cf-profile-deviance.csv'))
print(dev[h %in% c(0, 0.2, 0.5, 1, 2, Inf), .(t_half, f_m, h, dev70 = round(dev70, 2), h_star)])

## ---- 8. the information a randomization test can use ----
## For a linear statistic with label-invariant weights, beta = C + K
## with K = wbar' Ybar_all fixed under relabeling, so a randomization
## test uses only the contrast part C = sum_p w_p (omega_p - wbar)'
## Ybar_p, whose mean depends only on between-path differences of the
## drug profile. Power of the model-based test (z on the whole
## estimate) against the randomization test (|C| against the
## permutation spread, which uses the path mixture of the covariance
## and of the means), normal approximation, N = 70, h = t1/2 = 0.2.
rpow <- rbindlist(lapply(DES, function(des) rbindlist(lapply(c('F1', 'VF', 'S'), function(ty)
  rbindlist(lapply(c(0, 0.2, 0.4, 0.6, 0.8, 1), function(f) {
    pms <- design_moments(des, 0.2, f)
    Xs <- lapply(pms, design_X, h = 0.2, type = ty)
    g <- gls_plim(Xs, rc_plim(Xs, pms)$V, pms)
    P <- length(pms); wb <- Reduce('+', g$omega) / P
    mu_bar <- Reduce('+', lapply(pms, `[[`, 'mu')) / P
    Smix <- Reduce('+', lapply(pms, function(pm) pm$S + tcrossprod(pm$mu - mu_bar))) / P
    mC <- mean(mapply(function(om, pm) sum((om - wb) * pm$mu), g$omega, pms))
    vC <- mean(mapply(function(om, pm) sum((om - wb) * (pm$S %*% (om - wb))), g$omega, pms)) / 70
    vP <- mean(sapply(g$omega, function(om) sum((om - wb) * (Smix %*% (om - wb))))) / 70
    zc <- qnorm(0.975) * sqrt(vP / vC)
    data.table(design = des, type = ty, f_m = f, plim = g$b, se = sqrt(g$var),
               power_model = rej(g$b, sqrt(g$var)), mean_C = mC, se_C = sqrt(vC),
               power_rand = pnorm(mC / sqrt(vC) - zc) + pnorm(-mC / sqrt(vC) - zc),
               share_K = 1 - mC / g$b)
  }))))))
fwrite(rpow, file.path(OUT, '58-cf-randomization-power.csv'))
print(rpow[, .(design, type, f_m, plim = round(plim, 2), mean_C = round(mean_C, 2),
               power_model = round(power_model, 3), power_rand = round(power_rand, 3))])

tm <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
fwrite(data.table(script = '58-plim-and-power.R', seconds = tm), file.path(OUT, '58-cf-timing.csv'))
cat(sprintf('total %.1f s\n', tm))
