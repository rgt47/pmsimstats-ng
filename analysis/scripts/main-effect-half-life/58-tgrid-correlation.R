## docs/58, the joint null distribution of the drug-coefficient
## statistics T(h) across the analysis half-life, in closed form:
## their correlation from the D_bc(h) regressors and the working
## covariance, the max-|t| critical value from the multivariate normal,
## the size of AIC selection with the tested term in the model (which
## for the main effect is the max-|t| statistic referred to a single
## normal critical value), Davies' bound for a continuous range of h,
## and the weak-null variance ratio that governs an unstudentized
## permutation test.
##
## Run from the repository root:
##   Rscript analysis/scripts/main-effect-half-life/58-tgrid-correlation.R
## Writes analysis/data/quick-sim/main-effect-half-life/58-tg-*.csv.
source('analysis/scripts/main-effect-half-life/58-common.R')
load_driver()
t_start <- Sys.time()
## qmvnorm and pmvnorm integrate by quasi-Monte Carlo; fix the stream.
set.seed(5800L)
DES <- c('Hybrid', 'OL+BDC', 'CO')

## Covariance of beta(h) and beta(h') for linear GLS estimators with
## fixed working covariances: (1/N) sum_p w_p omega_p(h)'(S_p + r_p(h)
## r_p(h')')omega_p(h'). With the working covariance fixed at its null
## limit for each h.
joint_null <- function(des, type, hs, th = 0.2, tau = 0, N = 70) {
  pms <- design_moments(des, th, 0, tau)
  gs <- lapply(hs, function(h) {
    Xs <- lapply(pms, design_X, h = h, type = type)
    rc <- rc_plim(Xs, pms)
    gls_plim(Xs, rc$V, pms, N = N)
  })
  K <- length(hs); C <- matrix(0, K, K)
  for (a in seq_len(K)) for (b in seq_len(K)) C[a, b] <- mean(mapply(function(oa, ob, pm, ra, rb)
    sum(oa * ((pm$S + tcrossprod(ra, rb)) %*% ob)),
    gs[[a]]$omega, gs[[b]]$omega, pms, gs[[a]]$res, gs[[b]]$res)) / N
  list(b0 = sapply(gs, `[[`, 'b'), C = C, R = cov2cor(C), se = sqrt(diag(C)))
}

summ <- list(); cors <- list()
for (des in DES) for (ty in c('F1', 'VF', 'S')) {
  J <- joint_null(des, ty, G5)
  mu <- J$b0 / J$se
  cmax <- qmvnorm(0.95, corr = J$R, tail = 'both.tails')$quantile
  size_aic <- 1 - as.numeric(pmvnorm(lower = rep(-qnorm(0.975), 5), upper = rep(qnorm(0.975), 5),
                                     mean = mu, corr = J$R))
  size_aic0 <- 1 - as.numeric(pmvnorm(lower = rep(-qnorm(0.975), 5), upper = rep(qnorm(0.975), 5),
                                      corr = J$R))
  size_max <- 1 - as.numeric(pmvnorm(lower = rep(-cmax, 5), upper = rep(cmax, 5), mean = mu,
                                     corr = J$R))
  summ[[length(summ) + 1]] <- data.table(design = des, type = ty, c_max = cmax,
    size_naive_maxt = size_aic, size_naive_maxt_nobias = size_aic0,
    size_calibrated_maxt = size_max, size_fixed_h02 = rej(J$b0[2], J$se[2]),
    min_cor = min(J$R), cor_0_2 = J$R[1, 5])
  cors[[length(cors) + 1]] <- data.table(design = des, type = ty, h1 = rep(G5, 5),
                                          h2 = rep(G5, each = 5), cor = as.vector(J$R))
}
summ <- rbindlist(summ); cors <- rbindlist(cors)
print(summ)
print(dcast(cors[design == 'Hybrid' & type == 'F1'], h1 ~ h2, value.var = 'cor'))
fwrite(summ, file.path(OUT, '58-tg-maxt-summary.csv'))
fwrite(cors, file.path(OUT, '58-tg-correlation.csv'))

## Davies' bound over a continuous range of h. For a unit-variance
## Gaussian process Z(h), P(sup Z > u) <= Phi(-u) + V exp(-u^2 / 2) /
## (2 pi), with V the total variation of the process's direction,
## approximated by the arc length sum_j arccos R(h_j, h_{j+1}) on a fine
## grid; two-sided by the union of the two tails.
hfine <- c(0, 0.02, 0.05, seq(0.1, 1, 0.05), seq(1.1, 2, 0.1), 2.5, 3, 4)
dav <- rbindlist(lapply(DES, function(des) rbindlist(lapply(c('F1', 'VF'), function(ty) {
  J <- joint_null(des, ty, hfine)
  r1 <- pmin(1, J$R[cbind(seq_len(length(hfine) - 1), 2:length(hfine))])
  arc <- cumsum(c(0, acos(r1)))
  rbindlist(lapply(c(2, 4), function(hmax) {
    V <- arc[which(hfine == hmax)]
    bound <- function(u) 2 * pnorm(-u) + V * exp(-u^2 / 2) / pi
    ub <- uniroot(function(u) bound(u) - 0.05, c(1.5, 5))$root
    kk <- which(hfine <= hmax)
    cm <- qmvnorm(0.95, corr = J$R[kk, kk], tail = 'both.tails')$quantile
    data.table(design = des, type = ty, h_max = hmax, V = V, bound_at_1.96 = bound(qnorm(0.975)),
               u_davies = ub, c_max_fine = cm)
  }))
}))))
print(dav)
fwrite(dav, file.path(OUT, '58-tg-davies.csv'))

## Weak null: an unstudentized linear statistic beta = C + K, with
## K = wbar' Ybar_all invariant to relabeling and C = sum_p w_p
## (omega_p - wbar)' Ybar_p the contrast part. The true variance of C
## uses path-specific covariances S_p = S + tau^2 q_p q_p'; the
## permutation variance uses their mixture. Ratio > 1: anticonservative.
wn <- rbindlist(lapply(DES, function(des) rbindlist(lapply(c('F1', 'VF', 'S'), function(ty)
  rbindlist(lapply(c(0, 2.75, 5.49), function(tau) {
    pms <- design_moments(des, 0.2, 0, tau)
    Xs <- lapply(pms, design_X, h = 0.2, type = ty)
    g <- gls_plim(Xs, rc_plim(Xs, pms)$V, pms)
    P <- length(pms); wb <- Reduce('+', g$omega) / P
    Smix <- Reduce('+', lapply(pms, `[[`, 'S')) / P
    vt <- mean(mapply(function(om, pm) sum((om - wb) * (pm$S %*% (om - wb))), g$omega, pms))
    vp <- mean(sapply(g$omega, function(om) sum((om - wb) * (Smix %*% (om - wb)))))
    data.table(design = des, type = ty, tau = tau, var_true = vt / 70, var_perm = vp / 70,
               ratio = vt / vp, size_unstud = 2 * pnorm(-qnorm(0.975) / sqrt(vt / vp)),
               wbar_norm = sqrt(sum(wb^2)), b0 = g$b)
  }))))))
print(wn)
fwrite(wn, file.path(OUT, '58-tg-weak-null-permutation.csv'))

## Power of the max-|t| test (critical value c_max from the null
## correlation) and of the AIC choice at 1.96, against the fixed-h
## test, under alternatives (f_m = 0.4, N = 70).
joint_alt <- function(des, type, hs, th, f, N = 70) {
  pms <- design_moments(des, th, f)
  gs <- lapply(hs, function(h) {
    Xs <- lapply(pms, design_X, h = h, type = type)
    gls_plim(Xs, rc_plim(Xs, pms)$V, pms, N = N)
  })
  K <- length(hs); C <- matrix(0, K, K)
  for (a in seq_len(K)) for (b in seq_len(K)) C[a, b] <- mean(mapply(function(oa, ob, pm, ra, rb)
    sum(oa * ((pm$S + tcrossprod(ra, rb)) %*% ob)),
    gs[[a]]$omega, gs[[b]]$omega, pms, gs[[a]]$res, gs[[b]]$res)) / N
  list(b = sapply(gs, `[[`, 'b'), R = cov2cor(C), se = sqrt(diag(C)))
}
pw <- rbindlist(lapply(c('Hybrid', 'OL+BDC', 'CO'), function(des) rbindlist(lapply(c('F1', 'VF'), function(ty)
  rbindlist(lapply(c(0.2, 1), function(th) {
    J0 <- joint_null(des, ty, G5)
    cmax <- qmvnorm(0.95, corr = J0$R, tail = 'both.tails')$quantile
    J <- joint_alt(des, ty, G5, th, 0.4)
    mu <- J$b / J$se
    data.table(design = des, type = ty, t_half = th, c_max = cmax,
               power_fixed_02 = rej(J$b[2], J$se[2]),
               power_fixed_true = rej(J$b[match(th, G5)], J$se[match(th, G5)]),
               power_maxt = 1 - as.numeric(pmvnorm(lower = rep(-cmax, 5), upper = rep(cmax, 5),
                                                   mean = mu, corr = J$R)),
               power_aic_naive = 1 - as.numeric(pmvnorm(lower = rep(-qnorm(0.975), 5),
                                                        upper = rep(qnorm(0.975), 5), mean = mu, corr = J$R)))
  }))))))
print(pw)
fwrite(pw, file.path(OUT, '58-tg-maxt-power.csv'))

tm <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
fwrite(data.table(script = '58-tgrid-correlation.R', seconds = tm), file.path(OUT, '58-tg-timing.csv'))
cat(sprintf('total %.1f s\n', tm))
