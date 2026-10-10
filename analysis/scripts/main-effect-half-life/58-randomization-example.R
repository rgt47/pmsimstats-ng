## docs/58, worked example: one simulated Hybrid trial of the refined
## process (anchored carryover, randomized allocation, separable AR(1)
## at rho = 0.7, N = 70) with a drug effect, f_m = 0.4 (m_BR = 4.39
## points) and true t1/2 = 0.2 week, no heterogeneity. The seed (5801)
## was fixed before the trial was generated, and the trial was not
## selected for any property.
##
## Reports: T(h) across the half-life for F1 and VF; the profile
## likelihood of h; fixed-h, AIC-selected and AIC-averaged estimates of
## beta_D^exp; the max-|t| p-value by Bonferroni, by the multivariate
## normal and by randomization; randomization tests by relabeling with
## refits (199) and with label-invariant null-model weights (999, by
## permutation and by re-randomization); a parametric bootstrap from
## the null model (199); and confidence intervals that account for the
## half-life choice.
##
## Run from the repository root:
##   Rscript analysis/scripts/main-effect-half-life/58-randomization-example.R [--cores K] [--refit M] [--boot M] [--inv M]
## (run for docs/58 with --refit 49 --boot 49 --inv 199 under heavy machine load)
## Writes analysis/data/quick-sim/main-effect-half-life/58-ex-*.csv.
source('analysis/scripts/main-effect-half-life/58-common.R')
args <- commandArgs(trailingOnly = TRUE)
iarg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
cores <- iarg('--cores', 7L)
M_REFIT <- iarg('--refit', 199L); M_BOOT <- iarg('--boot', 199L); M_INV <- iarg('--inv', 499L)
load_driver('anchored', 'random')
t_start <- Sys.time(); tm <- list()
tick <- function(nm, t0) tm[[nm]] <<- as.numeric(difftime(Sys.time(), t0, units = 'secs'))

cl <- mk_cl(0.2); F_M <- 0.4
set.seed(5801L)
Ns <- alloc(list(N = 70L), NULL)
d <- make_trial_me(cl, F_M, 0, Ns)
cat('path sizes:', Ns, '\n')
pms <- design_moments('Hybrid', 0.2, F_M)
truth <- data.table(f_m = F_M, m_BR = F_M * M_BR, beta_exp = estimand(pms, F_M, 'exp'),
                    beta_bl = estimand(pms, F_M, 'bl'), Ns = paste(Ns, collapse = '/'))
print(truth); fwrite(truth, file.path(OUT, '58-ex-truth.csv'))

## ---- grids and procedures ----
t0 <- Sys.time()
GF <- grid_me(d, form = FMS$F1); GV <- grid_me(d, form = FMS$VF)
PF <- procedures(GF)[, type := 'F1']; PV <- procedures(GV)[, type := 'VF']
tick('grids (F1 and VF, five half-lives)', t0)
gtab <- rbind(GF$g[, type := 'F1'], GV$g[, type := 'VF'])
print(gtab); print(rbind(PF, PV))
fwrite(gtab, file.path(OUT, '58-ex-grid.csv'))
## correlation of the five T(h) estimated from the influence vectors
RF <- cov2cor(crossprod(GF$U)); RV <- cov2cor(crossprod(GV$U))
cmaxF <- qmvnorm(0.95, corr = RF, tail = 'both.tails')$quantile
cmaxV <- qmvnorm(0.95, corr = RV, tail = 'both.tails')$quantile
fwrite(data.table(type = rep(c('F1', 'VF'), each = 25), h1 = rep(G5, 10), h2 = rep(rep(G5, each = 5), 2),
                  cor = c(as.vector(RF), as.vector(RV))), file.path(OUT, '58-ex-tcor.csv'))

## ---- fine profile of h (ML, F1), and the t-statistic over h ----
t0 <- Sys.time()
hfine <- c(0, 0.02, 0.05, seq(0.1, 1, 0.05), seq(1.1, 3, 0.1), 4, 6, 10, 20, Inf)
prof <- rbindlist(lapply(hfine, function(h) {
  dd <- with_xh(d, h)
  fF <- fit_lme(FMS$F1, dd); fV <- fit_lme(FMS$VF, dd)
  data.table(h = h, r = if (h == 0) 0 else if (!is.finite(h)) 1 else 0.5^(1 / h),
             llF = ll_or_na(fit_lme(FMS$F1, dd, 'ML')), llV = ll_or_na(fit_lme(FMS$VF, dd, 'ML')),
             tF = cr2x(fF)['t'], tV = cr2x(fV)['t'])
}))
ll0 <- ll_or_na(fit_lme(F1_0, d, 'ML'))
opt <- optimize(function(lh) ll_or_na(fit_lme(FMS$F1, with_xh(d, exp(lh)), 'ML')),
                log(c(0.01, 50)), maximum = TRUE)
lmax <- max(prof$llF, opt$objective, na.rm = TRUE)
prof[, `:=`(devF = 2 * (lmax - llF), LR = 2 * (llF - ll0))]
pint <- data.table(h_ml_continuous = exp(opt$maximum), h_ml_fine = prof$h[which.max(prof$llF)],
                   lower_chi1 = min(prof[devF <= qchisq(0.95, 1), h]),
                   upper_chi1 = max(prof[devF <= qchisq(0.95, 1), h]),
                   upper_mix = max(prof[devF <= qchisq(0.90, 1), h]),
                   dev_at_0 = prof[h == 0, devF], dev_at_1 = prof[h == 1, devF],
                   dev_at_inf = prof[!is.finite(h), devF], LR_at_0 = prof[h == 0, LR],
                   LR_max = max(prof$LR, na.rm = TRUE))
## Davies' bound for sup |T(h)| over h in [0, 4] weeks with the
## empirical variation V = sum |T(h_j) - T(h_{j-1})| (the form of
## segmented::davies.test, two tails by the union bound).
dav <- rbindlist(lapply(c('tF', 'tV'), function(k) {
  z <- prof[is.finite(h) & h <= 4][[k]]; M <- max(abs(z)); V <- sum(abs(diff(z)))
  data.table(type = if (k == 'tF') 'F1' else 'VF', M = M, V = V,
             p_davies = min(1, 2 * pnorm(-M) + V * exp(-M^2 / 2) / sqrt(2 * pi)))
}))
print(dav); fwrite(dav, file.path(OUT, '58-ex-davies.csv'))
print(pint)
tick('fine profile (54 half-lives, F1 and VF)', t0)
fwrite(prof, file.path(OUT, '58-ex-profile.csv')); fwrite(pint, file.path(OUT, '58-ex-profile-interval.csv'))

## ---- randomization tests with label-invariant null-model weights ----
t0 <- Sys.time()
nv <- null_V(d)
S <- perm_setup(d, nv$V)
RNGkind("L'Ecuyer-CMRG"); set.seed(5802L)
fast <- mclapply(c('permute', 'rerandomize'), function(sc) rand_test(S, 999, sc), mc.cores = 2)
tick('fast randomization tests (999 permutations and 999 re-randomizations)', t0)
fastp <- data.table(statistic = names(fast[[1]]$p), observed = fast[[1]]$obs[names(fast[[1]]$p)],
                    p_permute = fast[[1]]$p, p_rerandomize = fast[[2]]$p)
print(fastp); fwrite(fastp, file.path(OUT, '58-ex-randomization-fast.csv'))
fwrite(data.table(t(fast[[1]]$sims)), file.path(OUT, '58-ex-randomization-fast-sims.csv'))

## ---- randomization test with refits (199 permutations of the labels) ----
## Each relabeling refits F1 and VF at every grid half-life (REML, CR2)
## and repeats the AIC selection (ML, with the tested term), so the
## reference distribution includes the selection.
stat_refit <- function(dd) {
  a <- grid_me(dd, form = FMS$F1, infl = FALSE)$g; b <- grid_me(dd, form = FMS$VF, infl = FALSE)$g
  c(F1_fix_t = abs(a$t[a$h == 0.2]), F1_aic_t = abs(a$t[which.max(a$ll1)]), F1_max_t = max(abs(a$t)),
    VF_fix_t = abs(b$t[b$h == 0.2]), VF_aic_t = abs(b$t[which.max(b$ll1)]), VF_max_t = max(abs(b$t)))
}
t0 <- Sys.time()
obs_r <- stat_refit(d)
path0 <- unique(d[, .(ptID, path)])$path
set.seed(5803L)
sims_r <- do.call(rbind, mclapply(seq_len(M_REFIT), function(m) stat_refit(relabel(d, sample(path0))),
                                  mc.cores = cores, mc.set.seed = TRUE))
p_refit <- sapply(names(obs_r), function(k) (1 + sum(sims_r[, k] >= obs_r[k] - 1e-12)) / (M_REFIT + 1))
tick(sprintf('randomization test with refits (%d permutations)', M_REFIT), t0)
refit <- data.table(statistic = names(obs_r), observed = obs_r, p_refit = p_refit,
                    crit95 = apply(sims_r, 2, quantile, 0.95))
print(refit); fwrite(refit, file.path(OUT, '58-ex-randomization-refit.csv'))

## ---- parametric bootstrap from the null model (199) ----
t0 <- Sys.time()
f0 <- fit_lme(F1_0, d)
sim0 <- function() {
  mu <- drop(model.matrix(F1_0[-2], d) %*% fixef(f0)); y <- mu
  for (ii in split(seq_len(nrow(d)), d$ptID)) {
    tt <- d$t[ii]; Sg <- nv$su2 + nv$s2 * nv$phi^abs(outer(tt, tt, '-'))
    y[ii] <- mu[ii] + drop(crossprod(chol(Sg), rnorm(length(ii))))
  }
  dd <- copy(d); dd[, Sx := y]; dd
}
set.seed(5804L)
sims_b <- do.call(rbind, mclapply(seq_len(M_BOOT), function(m) stat_refit(sim0()), mc.cores = cores,
                                  mc.set.seed = TRUE))
p_boot <- sapply(names(obs_r), function(k) (1 + sum(sims_b[, k] >= obs_r[k] - 1e-12)) / (M_BOOT + 1))
tick(sprintf('parametric bootstrap from the null model (%d)', M_BOOT), t0)
boot <- data.table(statistic = names(obs_r), observed = obs_r, p_boot = p_boot,
                   crit95 = apply(sims_b, 2, quantile, 0.95))
print(boot); fwrite(boot, file.path(OUT, '58-ex-bootstrap.csv'))
fwrite(as.data.table(sims_r), file.path(OUT, '58-ex-randomization-refit-sims.csv'))

## ---- confidence intervals for beta_D^exp that account for the choice ----
ci_rows <- function(G, P, cmax, type) {
  g <- G$g; ee <- g$est * g$dbar; se <- g$se * g$dbar
  k2 <- which(g$h == 0.2); ka <- which.max(g$ll1)
  w <- exp(g$ll1 - max(g$ll1)); w <- w / sum(w)
  em <- sum(w * ee); sm <- sum(w * sqrt(se^2 + (ee - em)^2)); dfm <- sum(w * g$df)
  rbind(
    data.table(type = type, interval = 'fixed h = 0.2, CR2', h = 0.2, est = ee[k2],
               lo = ee[k2] - qt(0.975, g$df[k2]) * se[k2], hi = ee[k2] + qt(0.975, g$df[k2]) * se[k2]),
    data.table(type = type, interval = 'AIC choice, naive', h = g$h[ka], est = ee[ka],
               lo = ee[ka] - qt(0.975, g$df[ka]) * se[ka], hi = ee[ka] + qt(0.975, g$df[ka]) * se[ka]),
    data.table(type = type, interval = 'AIC choice, max-|t| constant', h = g$h[ka], est = ee[ka],
               lo = ee[ka] - cmax * se[ka], hi = ee[ka] + cmax * se[ka]),
    data.table(type = type, interval = 'AIC-weight average (Buckland)', h = sum(w * g$h), est = em,
               lo = em - qt(0.975, dfm) * sm, hi = em + qt(0.975, dfm) * sm),
    data.table(type = type, interval = 'union over grid, max-|t| constant', h = NA_real_, est = NA_real_,
               lo = min(ee - cmax * se), hi = max(ee + cmax * se)))
}
cis <- rbind(ci_rows(GF, PF, cmaxF, 'F1'), ci_rows(GV, PV, cmaxV, 'VF'))

## Randomization inversion under constant-shape sharp hypotheses
## H(b, h): the drug effect is b * D_bc(h) at every visit for everyone.
## Under H the adjusted outcome Y - b D_bc(h) is invariant to the path
## labels; test with the VF statistic (499 permutations), invert over
## b, convert to the beta_D^exp scale with dbar(h), and take the union
## over the grid.
t0 <- Sys.time()
bgrid <- seq(-9, 4, by = 0.25)
inv <- rbindlist(mclapply(seq_len(length(G5) * length(bgrid)), function(j) {
  h <- G5[(j - 1) %/% length(bgrid) + 1]; b <- bgrid[(j - 1) %% length(bgrid) + 1]
  dd <- with_xh(d, h); dd[, Sx := Sx - b * X]
  Sj <- perm_setup(dd, null_V(dd)$V, hs = h, types = 'VF')
  set.seed(5805L + j)
  rt <- rand_test(Sj, M_INV)
  data.table(h = h, b = b, p = rt$p[['VF_fix_t']])
}, mc.cores = cores, mc.set.seed = TRUE))
tick(sprintf('randomization inversion (5 half-lives x 53 values, %d permutations)', M_INV), t0)
fwrite(inv, file.path(OUT, '58-ex-inversion.csv'))
dbh <- GV$g[, .(h, dbar)]
inv_ci <- inv[p > 0.05, .(lo = min(b), hi = max(b)), by = h][dbh, on = 'h']
inv_ci[, `:=`(lo_exp = lo * dbar, hi_exp = hi * dbar)]
print(inv_ci)
cis <- rbind(cis, data.table(type = 'VF', interval = 'randomization inversion, union over grid',
                             h = NA_real_, est = NA_real_, lo = min(inv_ci$lo_exp), hi = max(inv_ci$hi_exp)),
             data.table(type = 'VF', interval = 'randomization inversion, h = 0.2', h = 0.2, est = NA_real_,
                        lo = inv_ci[h == 0.2, lo_exp], hi = inv_ci[h == 0.2, hi_exp]))
cis[, covers := lo <= truth$beta_exp & truth$beta_exp <= hi]
print(cis); fwrite(cis, file.path(OUT, '58-ex-intervals.csv')); fwrite(inv_ci, file.path(OUT, '58-ex-inversion-ci.csv'))
fwrite(data.table(type = c('F1', 'VF'), c_max = c(cmaxF, cmaxV)), file.path(OUT, '58-ex-cmax.csv'))

tm[['total']] <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
tmd <- data.table(step = names(tm), seconds = unlist(tm), cores = cores)
print(tmd); fwrite(tmd, file.path(OUT, '58-ex-timing.csv'))
