## docs/56, Worked example 1: one simulated Hybrid trial (N = 70) from
## the stress-test construct with a prognostic biomarker (constant
## natural-course pathway, c_bm,TV = 0.3), analyzed by S and by F1
## exactly as the stress-test driver fits them; then the same seed with
## an interaction (c_bm = 0.30). Also: the D_bc coding by path, the
## pre-exposure guard on the crossover, the drug-expectancy correlation,
## the visit-group weights, and a large trial (N = 2000) that
## approximates the probability limits of the two estimates.
##
## Run from the repository root:
##   Rscript analysis/scripts/recommended-analysis/56-example-trial.R [out_dir]
## Tables go to out_dir (default analysis/data/recommended-analysis/);
## the figure to docs/figures/56-fig1-fitted-slopes.png.
suppressMessages({
  library(data.table); library(nlme); library(clubSandwich)
})
source('analysis/scripts/recommended-analysis/56-load-driver.R')
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else 'analysis/data/recommended-analysis'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()

env <- load_driver('recursive', 'fixed')
SEED <- 20261009L
TH_A <- env$TH_A
stopifnot(TH_A == 0.5)

## ---------------- D_bc by path, and the guard ----------------
dbc_rows <- rbindlist(lapply(seq_along(env$PIS), function(k) {
  pi <- env$PIS[[k]]
  data.table(design = 'Hybrid', path = k, week = pi$w, on = as.integer(pi$on),
             eta = pi$e, tsd = pi$tsd, Dbc = env$dbc(pi$on, pi$tsd, TH_A))
}))
co <- env$dgp$designs$CO
pi_co <- env$dgp$path_info(co$t_wk, co$e, co$paths[[2]])
noguard <- function(Db, tsd, h) ifelse(Db, 1, 0.5^(tsd / h))
dbc_co <- data.table(design = 'CO', path = 2, week = pi_co$w,
                     on = as.integer(pi_co$on), eta = pi_co$e, tsd = pi_co$tsd,
                     Dbc = env$dbc(pi_co$on, pi_co$tsd, 1),
                     Dbc_noguard = noguard(pi_co$on, pi_co$tsd, 1))
fwrite(dbc_rows, file.path(out_dir, '56-dbc-hybrid.csv'))
fwrite(dbc_co, file.path(out_dir, '56-dbc-co-guard.csv'))
cat('\nD_bc at t_half_a = 0.5 week, Hybrid paths (rows = visits):\n')
print(dcast(dbc_rows, week + eta ~ path, value.var = 'Dbc'), digits = 3)
cat('\nCO placebo-first path at t_half_a = 1 week, with and without guard:\n')
print(dbc_co)

## ---------------- design facts: weights and Cor(D, eta) ----------------
Ns <- env$cell_setup(env$cell('curve-C07', 0))$Ns
long_design <- rbindlist(lapply(seq_along(env$PIS), function(k) {
  pi <- env$PIS[[k]]
  data.table(n = Ns[k], t = c(0, pi$w), D = c(0, as.integer(pi$on)),
             eta = c(0, pi$e), X = c(0, env$dbc(pi$on, pi$tsd, TH_A)))
}))
long_design[, group := fifelse(t == 0, 'baseline', fifelse(eta == 1, 'OL on',
                               fifelse(D == 1, 'blinded on', 'blinded off')))]
wcor <- function(x, y, w) cov.wt(cbind(x, y), wt = w / sum(w), cor = TRUE)$cor[1, 2]
facts <- data.table(
  item = c('allocation', 'w_baseline', 'w_OL_on', 'w_blinded_on', 'w_blinded_off',
           'cor_D_eta_9rows', 'cor_Dbc_eta_9rows', 'cor_D_eta_post'),
  value = c(paste(Ns, collapse = '/'),
            long_design[, .(w = sum(n)), by = group][
              match(c('baseline', 'OL on', 'blinded on', 'blinded off'), group), w],
            round(long_design[, wcor(D, eta, n)], 4),
            round(long_design[, wcor(X, eta, n)], 4),
            round(long_design[t > 0, wcor(D, eta, n)], 4)))
fwrite(facts, file.path(out_dir, '56-design-facts.csv'))
cat('\nDesign facts:\n'); print(facts)

## ---------------- fit helpers ----------------
fit_pair <- function(d) {
  dd <- env$with_x(d, TH_A)
  list(S = env$fit_lme(env$F_S, dd), F1 = env$fit_lme(env$F_F1, dd), data = dd)
}
inference <- function(fit, name) {
  tt <- summary(fit)$tTable
  ct <- coef_test(fit, vcov = 'CR2', test = 'Satterthwaite')
  k <- which(ct$Coef == 'bmc:X')
  data.table(analysis = name, est = tt['bmc:X', 'Value'],
             se_mb = tt['bmc:X', 'Std.Error'], df_mb = tt['bmc:X', 'DF'],
             p_mb = tt['bmc:X', 'p-value'],
             se_cr2 = ct$SE[k], df_cr2 = ct$df_Satt[k], p_cr2 = ct$p_Satt[k])
}
coefs <- function(fit, name) {
  tt <- summary(fit)$tTable
  ct <- coef_test(fit, vcov = 'CR2', test = 'Satterthwaite')
  data.table(analysis = name, term = rownames(tt), est = tt[, 'Value'],
             se_mb = tt[, 'Std.Error'], se_cr2 = ct$SE, df_cr2 = ct$df_Satt,
             p_cr2 = ct$p_Satt)
}
## Biomarker slope implied by a fit at each row, and averaged by group
slope_rows <- function(fit, dd) {
  b <- fixef(fit)
  g <- function(nm) if (nm %in% names(b)) b[[nm]] else 0
  g('bmc') + g('bmc:X') * dd$X + g('bmc:De') * dd$De
}
add_group <- function(dd) {
  dd[, path := as.integer(as.character(ptID)) %/% 1000 + 1]
  dd[, group := fifelse(t == 0, 'baseline', fifelse(De == 1, 'OL on',
                        fifelse(Db, 'blinded on', 'blinded off')))]
  dd
}
## True slope of Y on B (per raw unit) at each path and visit, from the
## construct's covariance: Cov(B, BL - TV - PB - BR) / Var(B)
true_slopes <- function(cl) {
  st <- env$cell_setup(cl)
  rbindlist(lapply(seq_along(env$PIS), function(k) {
    S <- st$bs[[k]]$Sigma; n <- st$bs[[k]]$n
    i_tv <- 3:(2 + n); i_pb <- (3 + n):(2 + 2 * n); i_br <- (3 + 2 * n):(2 + 3 * n)
    s_post <- (S[1, 2] - S[1, i_tv] - S[1, i_pb] - S[1, i_br]) / S[1, 1]
    data.table(path = k, t = c(0, env$PIS[[k]]$w), s_true = c(S[1, 2] / S[1, 1], s_post),
               n = st$Ns[k])
  }))
}

run_trial <- function(cl, label, seed, n_total = NULL) {
  if (!is.null(n_total)) {
    cl <- copy(cl)
    cl[, N := as.integer(n_total)]
  }
  st <- env$cell_setup(cl)
  set.seed(seed)
  d <- env$make_trial(cl, st$bs, env$alloc(cl, st), 1)
  t0 <- Sys.time()
  fp <- fit_pair(d)
  el <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
  dd <- add_group(fp$data)
  dd[, `:=`(slope_S = slope_rows(fp$S, dd), slope_F1 = slope_rows(fp$F1, dd))]
  ts <- true_slopes(cl)
  dd <- merge(dd, ts[, .(path, t, s_true)], by = c('path', 't'))
  emp <- dd[, .(emp = coef(lm(Sx ~ factor(paste(path, t)) + bmc))[['bmc']]), by = group]
  grp <- dd[, .(rows = .N, true = mean(s_true), S = mean(slope_S), F1 = mean(slope_F1)),
            by = group][emp, on = 'group']
  grp[, case := label]
  list(fits = fp, inf = rbind(inference(fp$S, 'S'), inference(fp$F1, 'F1'))[, case := label],
       coefs = rbind(coefs(fp$S, 'S'), coefs(fp$F1, 'F1'))[, case := label],
       groups = grp, secs = el, cl = cl)
}

cl_null <- env$cell('prog-TV-const', 0, c_tv = 0.3)
cl_alt <- env$cell('prog-TV-const', 0.30, c_tv = 0.3)
## The check in the driver itself: cr2() on the same fit gives the same
## estimate, SE and p-value as the coefficient test used here.
ex_null <- run_trial(cl_null, 'null', SEED)
ex_alt <- run_trial(cl_alt, 'alt', SEED)
chk <- env$cr2(ex_null$fits$F1)
stopifnot(abs(chk[['se']] - ex_null$inf[analysis == 'F1', se_cr2]) < 1e-12)

inf <- rbind(ex_null$inf, ex_alt$inf)
cf <- rbind(ex_null$coefs, ex_alt$coefs)
grp <- rbind(ex_null$groups, ex_alt$groups)
fwrite(inf, file.path(out_dir, '56-example-inference.csv'))
fwrite(cf, file.path(out_dir, '56-example-coefficients.csv'))
fwrite(grp, file.path(out_dir, '56-example-slopes-by-group.csv'))
cat('\nInference for bmc:X (seed', SEED, '):\n'); print(inf, digits = 4)
cat('\nCoefficients:\n'); print(cf, digits = 4)
cat('\nSlopes by visit group:\n'); print(grp, digits = 4)
cat(sprintf('\nfit time for S and F1 on one N = 70 trial: %.2f s (null), %.2f s (alt)\n',
            ex_null$secs, ex_alt$secs))

## True interaction and prognostic slope, per raw unit of the biomarker
s_prog <- -0.3 * 10 / env$dgp$bm_sd
beta_true <- -0.30 * 8 / env$dgp$bm_sd
cat(sprintf('s_prog = %.4f per unit (%.3f per SD); beta_bm:D at c_bm = 0.30 = %.4f\n',
            s_prog, s_prog * env$dgp$bm_sd, beta_true))

## ---------------- large trial ----------------
t0 <- Sys.time()
big_null <- run_trial(cl_null, 'null, N = 2000', SEED + 1L, n_total = 2000)
big_alt <- run_trial(cl_alt, 'alt, N = 2000', SEED + 1L, n_total = 2000)
stopifnot(sum(env$cell_setup(big_null$cl)$Ns) == 2000)
big <- rbind(big_null$inf, big_alt$inf)
big[, `:=`(s_prog = s_prog, ratio_to_sprog = est / s_prog)]
fwrite(big, file.path(out_dir, '56-example-large-trial.csv'))
fwrite(rbind(big_null$groups, big_alt$groups),
       file.path(out_dir, '56-example-large-slopes.csv'))
cat('\nLarge trials (N = 2000):\n'); print(big, digits = 4)
cat(sprintf('large-trial time: %.1f s\n', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

## ---------------- figure 1 ----------------
lev <- c('baseline', 'OL on', 'blinded on', 'blinded off')
labs <- c('Baseline', 'Open label\n(on drug)', 'Blinded\non drug', 'Blinded\noff drug')
grp_all <- rbind(grp, big_null$groups, big_alt$groups)
png('docs/figures/56-fig1-fitted-slopes.png', width = 2400, height = 2200, res = 220)
par(mfrow = c(2, 2), mar = c(5.2, 4.6, 3, 0.8), mgp = c(2.6, 0.7, 0), las = 1)
cols <- c(true = 'grey35', emp = 'grey60', S = '#B2182B', F1 = '#2166AC')
yl <- range(unlist(grp_all[, .(true, emp, S, F1)])) + c(-0.04, 0.04)
ttl <- c(null = 'N = 70, no interaction (c_bm = 0)', alt = 'N = 70, interaction (c_bm = 0.30)',
         `null, N = 2000` = 'N = 2000, no interaction (c_bm = 0)',
         `alt, N = 2000` = 'N = 2000, interaction (c_bm = 0.30)')
for (cs in names(ttl)) {
  g <- grp_all[case == cs][match(lev, group)]
  x <- seq_along(lev)
  plot(NA, xlim = c(0.6, 4.4), ylim = yl, xaxt = 'n',
       xlab = '', ylab = 'Slope of Y on the biomarker (points per unit)',
       main = ttl[[cs]], cex.main = 0.95)
  axis(1, at = x, labels = labs, padj = 0.5, cex.axis = 0.85)
  abline(h = 0, col = 'grey85')
  points(x, g$emp, pch = 4, col = cols['emp'], cex = 1.3, lwd = 2)
  lines(x, g$true, col = cols['true'], lwd = 2, lty = 3); points(x, g$true, pch = 15, col = cols['true'])
  lines(x - 0.06, g$S, col = cols['S'], lwd = 2); points(x - 0.06, g$S, pch = 17, col = cols['S'], cex = 1.2)
  lines(x + 0.06, g$F1, col = cols['F1'], lwd = 2, lty = 2); points(x + 0.06, g$F1, pch = 16, col = cols['F1'], cex = 1.2)
  if (cs == 'null') legend('bottomleft', bty = 'n', cex = 0.8,
    legend = c('true slope (construct)', 'observed slope in this trial', 'fitted by S', 'fitted by F1'),
    col = cols[c('true', 'emp', 'S', 'F1')], pch = c(15, 4, 17, 16), lty = c(3, NA, 1, 2), lwd = 2)
}
dev.off()
cat(sprintf('\ntotal elapsed %.1f s\n', as.numeric(difftime(Sys.time(), t_start, units = 'secs'))))
