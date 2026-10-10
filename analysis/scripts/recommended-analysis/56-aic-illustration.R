## docs/56, Section 6: a small illustration of AIC selection of the
## analysis half-life on null trials (c_bm = 0) of the stress-test
## construct, with the driver's code. For each trial:
##   - the ML AIC profile of F1 (with the biomarker-by-D_bc term, as F1A
##     selects) and of F1 without that term (the proposed F1A0), and the
##     CR2 t-statistic of the interaction, over a fine half-life grid;
##   - the half-life selected on the driver's grid {0, 0.25, 0.5, 1, 2}
##     by F1A and by F1A0, and the CR2 test after each selection.
## This is an illustration on a few trials, not a size estimate.
##
## Run from the repository root:
##   Rscript analysis/scripts/recommended-analysis/56-aic-illustration.R [out_dir] [n_sel]
## n_sel: trials per cell for the selection comparison (default 50).
suppressMessages({
  library(data.table); library(nlme); library(clubSandwich)
})
source('analysis/scripts/recommended-analysis/56-load-driver.R')
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else 'analysis/data/recommended-analysis'
n_sel <- if (length(args) >= 2) as.integer(args[2]) else 50L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
env <- load_driver('recursive', 'fixed')
F_F1_noint <- Sx ~ bmc + t + X + De + bmc:De
FINE <- c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3)

aic_ml <- function(f, dd) {
  fit <- env$fit_lme(f, dd, 'ML')
  if (is.null(fit)) NA_real_ else AIC(fit)
}
t_cr2 <- function(dd) {
  fit <- env$fit_lme(env$F_F1, dd)
  if (is.null(fit)) return(c(t = NA_real_, p = NA_real_, df = NA_real_))
  ct <- coef_test(fit, vcov = 'CR2', test = 'Satterthwaite')
  k <- which(ct$Coef == 'bmc:X')
  c(t = ct$tstat[k], p = ct$p_Satt[k], df = ct$df_Satt[k])
}
new_trial <- function(cl, seed) {
  st <- env$cell_setup(cl)
  set.seed(seed)
  env$make_trial(cl, st$bs, env$alloc(cl, st), 1)
}
cells <- list(base = env$cell('curve-C07', 0), thalf1 = env$cell('thalf1', 0, t_half = 1))

## ---------------- profiles ----------------
t0 <- Sys.time()
prof <- rbindlist(lapply(names(cells), function(cn) rbindlist(lapply(1:3, function(i) {
  d <- new_trial(cells[[cn]], 56000L + 100L * match(cn, names(cells)) + i)
  rbindlist(lapply(FINE, function(h) {
    dd <- env$with_x(d, h)
    tt <- t_cr2(dd)
    data.table(cell = cn, trial = i, h = h, aic_F1 = aic_ml(env$F_F1, dd),
               aic_noint = aic_ml(F_F1_noint, dd), t = tt[['t']], p = tt[['p']], df = tt[['df']])
  }))
}))))
prof[, `:=`(daic_F1 = aic_F1 - min(aic_F1, na.rm = TRUE),
            daic_noint = aic_noint - min(aic_noint, na.rm = TRUE)), by = .(cell, trial)]
t_prof <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
fwrite(prof, file.path(out_dir, '56-aic-profiles.csv'))
cat(sprintf('profiles: 6 trials x %d half-lives in %.1f s\n', length(FINE), t_prof))
print(prof[, .(cell, trial, h, daic_F1 = round(daic_F1, 2), daic_noint = round(daic_noint, 2),
               t = round(t, 3))], nrows = 200)

## ---------------- selection on the driver grid ----------------
t0 <- Sys.time()
sel <- rbindlist(lapply(names(cells), function(cn) rbindlist(lapply(seq_len(n_sel), function(i) {
  d <- new_trial(cells[[cn]], 57000L + 1000L * match(cn, names(cells)) + i)
  ds <- lapply(env$GRID, function(h) env$with_x(d, h))
  a1 <- sapply(ds, function(dd) aic_ml(env$F_F1, dd))
  a0 <- sapply(ds, function(dd) aic_ml(F_F1_noint, dd))
  h1 <- env$GRID[which.min(a1)]; h0 <- env$GRID[which.min(a0)]
  tf <- t_cr2(ds[[which(env$GRID == env$TH_A)]])
  ta <- t_cr2(ds[[which.min(a1)]]); tz <- t_cr2(ds[[which.min(a0)]])
  data.table(cell = cn, trial = i, h_F1A = h1, h_F1A0 = h0,
             t_F1 = tf[['t']], p_F1 = tf[['p']], t_F1A = ta[['t']], p_F1A = ta[['p']],
             t_F1A0 = tz[['t']], p_F1A0 = tz[['p']])
}))))
t_sel <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
fwrite(sel, file.path(out_dir, '56-aic-selection.csv'))
cat(sprintf('\nselection: %d trials in %.1f s\n', nrow(sel), t_sel))
cat('\nSelected half-life (counts):\n')
print(rbind(sel[, .(rule = 'F1A', .N), by = .(cell, h = h_F1A)],
            sel[, .(rule = 'F1A0', .N), by = .(cell, h = h_F1A0)])[order(cell, rule, h)])
summ <- sel[, .(trials = .N, same_h = mean(h_F1A == h_F1A0),
                mean_abs_t_F1 = mean(abs(t_F1)), mean_abs_t_F1A = mean(abs(t_F1A)),
                mean_abs_t_F1A0 = mean(abs(t_F1A0)),
                rej_F1 = sum(p_F1 < 0.05), rej_F1A = sum(p_F1A < 0.05),
                rej_F1A0 = sum(p_F1A0 < 0.05)), by = cell]
fwrite(summ, file.path(out_dir, '56-aic-selection-summary.csv'))
print(summ, digits = 3)

## ---------------- figure 3 ----------------
png('docs/figures/56-fig3-aic-profiles.png', width = 2400, height = 2000, res = 220)
par(mfrow = c(2, 2), mar = c(4.4, 4.6, 2.8, 0.8), mgp = c(2.5, 0.7, 0), las = 1)
tc <- c('#1B9E77', '#D95F02', '#7570B3')
for (cn in names(cells)) {
  z <- prof[cell == cn]
  yl <- c(0, max(c(z$daic_F1, z$daic_noint), na.rm = TRUE))
  plot(NA, xlim = range(FINE), ylim = yl, xlab = 'Analysis half-life (weeks)',
       ylab = 'AIC minus its minimum (ML)',
       main = sprintf('%s: AIC profiles', if (cn == 'base') 'True half-life 0.5 week' else 'True half-life 1 week'),
       cex.main = 0.95)
  for (i in 1:3) {
    w <- z[trial == i]
    lines(w$h, w$daic_F1, col = tc[i], lwd = 2)
    lines(w$h, w$daic_noint, col = tc[i], lwd = 2, lty = 2)
  }
  abline(h = 2, col = 'grey70', lty = 3)
  legend('topright', bty = 'n', cex = 0.75, lty = c(1, 2), lwd = 2,
         legend = c('F1 (selection of F1A)', 'F1 without bmc:X (F1A0)'))
}
for (cn in names(cells)) {
  z <- prof[cell == cn]
  plot(NA, xlim = range(FINE), ylim = range(c(z$t, -2.1, 2.1), na.rm = TRUE),
       xlab = 'Analysis half-life (weeks)', ylab = 'CR2 t-statistic of the interaction',
       main = 'Same trials: interaction t-statistic', cex.main = 0.95)
  abline(h = 0, col = 'grey80'); abline(h = c(-2, 2), col = 'grey70', lty = 3)
  for (i in 1:3) {
    w <- z[trial == i]
    lines(w$h, w$t, col = tc[i], lwd = 2)
    points(w$h[which.min(w$aic_F1)], w$t[which.min(w$aic_F1)], pch = 16, col = tc[i], cex = 1.4)
    points(w$h[which.min(w$aic_noint)], w$t[which.min(w$aic_noint)], pch = 1, col = tc[i], cex = 2, lwd = 2)
  }
  legend('bottomright', bty = 'n', cex = 0.75, pch = c(16, 1), pt.cex = c(1.4, 2),
         legend = c('half-life chosen by F1 AIC', 'half-life chosen by F1A0 AIC'))
}
dev.off()
cat('figure written\n')
