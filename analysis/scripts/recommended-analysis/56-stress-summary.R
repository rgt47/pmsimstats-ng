## docs/56, Worked example 2: Monte Carlo summaries read from the stored
## stress-test outputs (11-strawman-stress-test.R, run 2026-10-07 under
## the recursive carryover and fixed allocation; 12-strawman-stress-
## decision.R --strawman DS --no-symmetry --with-extras). No new
## simulation.
##
## Run from the repository root:
##   Rscript analysis/scripts/recommended-analysis/56-stress-summary.R [out_dir]
## Tables go to out_dir (default analysis/data/recommended-analysis/);
## the figure to docs/figures/56-fig2-size-power.png.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else 'analysis/data/recommended-analysis'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ss <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'

summ <- rbind(fread(file.path(ss, 'summary-full.csv')),
              fread(file.path(ss, 'summary-full-extra.csv'))[analysis != 'S'],
              fread(file.path(ss, 'summary-full-extra2.csv'))[analysis != 'S'])
stopifnot(!anyDuplicated(summ[, .(cell, analysis)]))
summ[, construct_label := fifelse(label == 'curve-CS', 'curve-CS', label)]
std_cells <- c('curve-C07', 'curve-CS', 'meas-error', 'N35', 'N150', 'thalf0', 'thalf1',
               'weibull0.5', 'weibull2', 'dropout-MCAR', 'dropout-MAR')
prog_cells <- c('prog-TV-const', 'prog-TV-grow', 'prog-PB', 'prog-BL')
show <- c('S', 'F1', 'ES', 'F1A', 'DS', 'G1', 'G2', 'R')

## ---------------- size by null cell ----------------
size <- summ[hyp == 'null' & analysis %in% show,
             .(label, analysis, reject, mcse, fits = fits - fails, kappa)]
size_w <- dcast(size, label ~ analysis, value.var = 'reject')
size_w <- size_w[match(c(std_cells, prog_cells), label), c('label', show), with = FALSE]
pooled <- size[label %in% std_cells, .(reject = sum(reject * fits) / sum(fits)), by = analysis]
pool_row <- as.data.table(c(list(label = 'pooled, 11 standard cells'),
                            as.list(setNames(pooled$reject, pooled$analysis)[show])))
size_w <- rbind(size_w, pool_row)
fwrite(size_w, file.path(out_dir, '56-size-by-cell.csv'))
cat('Null rejection rates (1,000 trials per cell; MCSE 0.0069 at 0.05):\n')
print(size_w, digits = 3)

## Check against the decision rule's pooled rates
scr <- fread(file.path(ss, 'decision-screen-ds-nosym-all14.csv'))
## The decision script excludes failed fits per trial; pooling the cell
## rates weighted by successful fits agrees to within one failed fit.
pool_diff <- max(abs(pooled[match(scr$analysis, analysis), reject] - scr$pooled), na.rm = TRUE)
cat(sprintf('largest difference from decision-screen pooled rates: %.1e\n', pool_diff))
stopifnot(pool_diff < 2e-4)

## ---------------- model-based versus CR2 on the same fit ----------------
pairs <- rbind(
  summ[hyp == 'null' & analysis %in% c('C1', 'C2'), .(pair = 'C1 (CR2) vs C2 (model-based)', label, analysis, reject, kappa, mean_se)],
  summ[hyp == 'null' & analysis %in% c('G1', 'G2'), .(pair = 'G1 (CR2) vs G2 (model-based)', label, analysis, reject, kappa, mean_se)])
mbcr <- dcast(pairs[label %in% c(std_cells, prog_cells)], pair + label ~ analysis,
              value.var = c('reject', 'kappa'))
fwrite(mbcr, file.path(out_dir, '56-cr2-vs-model-based.csv'))
cat('\nCR2 versus model-based inference on the same fit:\n'); print(mbcr, digits = 3)

kap <- dcast(summ[hyp == 'null' & label %in% c(std_cells, prog_cells) & analysis %in% c('S', 'F1', 'R', 'F1A')],
             label ~ analysis, value.var = 'kappa')
fwrite(kap, file.path(out_dir, '56-kappa.csv'))
cat('\nkappa = mean SE / empirical SD of the estimate, null cells:\n'); print(kap, digits = 3)

## ---------------- size-adjusted power ----------------
pw <- fread(file.path(ss, 'decision-power-ds-nosym-all14.csv'))
pw_show <- pw[analysis %in% show]
pw_w <- dcast(pw_show[cbm %in% c(0.30, 0.25)], label + cbm ~ analysis, value.var = 'power_adj')
pw_w <- pw_w[, c('label', 'cbm', show), with = FALSE]
fwrite(pw_w, file.path(out_dir, '56-power-adjusted.csv'))
cat('\nSize-adjusted power (500 trials per cell; MCSE at most 0.022):\n'); print(pw_w, digits = 3)
fwrite(pw_show, file.path(out_dir, '56-power-long.csv'))

## Mean and range of the cost of F1 relative to S over the standard
## AR(1) alternative cells at c_bm = 0.30
ar1_alt <- c('curve-C07', setdiff(std_cells, c('curve-C07', 'curve-CS')))
cost <- dcast(pw[label %in% ar1_alt & abs(cbm - 0.30) < 1e-9 & analysis %in% c('S', 'F1')],
              label ~ analysis, value.var = 'power_adj')
cost[, rel_loss := 1 - F1 / S]
fwrite(cost, file.path(out_dir, '56-cost-F1-vs-S.csv'))
cat('\nRelative loss of F1 against S, AR(1) cells at 0.30:\n'); print(cost, digits = 3)

## ---------------- F1A against F1 on the same trials ----------------
r1 <- readRDS(file.path(ss, 'replicates-full.rds'))
r3 <- readRDS(file.path(ss, 'replicates-full-extra2.rds'))
pair_chk <- merge(r1[analysis == 'S', .(cell, rep, est1 = est)],
                  r3[analysis == 'S', .(cell, rep, est3 = est)], by = c('cell', 'rep'))
cat(sprintf('\npairing check, S in main and extra2 runs: max |diff| = %.2e over %d trials\n',
            max(abs(pair_chk$est1 - pair_chk$est3), na.rm = TRUE), nrow(pair_chk)))
f <- merge(r1[analysis == 'F1', .(cell, rep, label, cbm, hyp, t_half, est_f = est, se_f = se, p_f = p)],
           r3[analysis == 'F1A', .(cell, rep, est_a = est, se_a = se, p_a = p, h = h)],
           by = c('cell', 'rep'))
f <- f[!is.na(p_f) & !is.na(p_a)]
f[, `:=`(t_f = est_f / se_f, t_a = est_a / se_a)]
nulls <- f[hyp == 'null' & label %in% std_cells]
choice <- nulls[, .(share = .N / nrow(nulls)), by = h][order(h)]
by_h <- nulls[, .(n = .N, size_F1 = mean(p_f < 0.05), size_F1A = mean(p_a < 0.05),
                  mean_abs_t_F1 = mean(abs(t_f)), mean_abs_t_F1A = mean(abs(t_a))), by = h][order(h)]
disc <- nulls[, .(trials = .N, both = sum(p_f < 0.05 & p_a < 0.05),
                  F1_only = sum(p_f < 0.05 & p_a >= 0.05), F1A_only = sum(p_f >= 0.05 & p_a < 0.05),
                  mean_abs_t_F1 = mean(abs(t_f)), mean_abs_t_F1A = mean(abs(t_a)),
                  sd_t_F1 = sd(t_f), sd_t_F1A = sd(t_a))]
th1 <- f[hyp == 'null' & label == 'thalf1']
th1_by_h <- th1[, .(n = .N, size_F1 = mean(p_f < 0.05), size_F1A = mean(p_a < 0.05)), by = h][order(h)]
choice_all <- f[, .(share = .N / nrow(f)), by = h][order(h)]
choice_cell <- dcast(f[hyp == 'null'], label ~ h, fun.aggregate = length, value.var = 'rep')
fwrite(by_h, file.path(out_dir, '56-f1a-by-selected-h.csv'))
fwrite(disc, file.path(out_dir, '56-f1a-discordance.csv'))
fwrite(th1_by_h, file.path(out_dir, '56-f1a-thalf1-by-h.csv'))
fwrite(choice_all, file.path(out_dir, '56-f1a-choice-all-cells.csv'))
cat('\nF1A selected half-life, all cells (cf. docs/47, Section 8):\n'); print(choice_all, digits = 3)
cat('\nF1A selections by null cell (counts):\n'); print(choice_cell)
cat('\nStandard null cells, by selected half-life:\n'); print(by_h, digits = 3)
cat('\nDiscordance, standard null cells (11,000 trials):\n'); print(disc, digits = 4)
cat('\nNull cell with true t_half = 1, by selected half-life:\n'); print(th1_by_h, digits = 3)

## ---------------- figure 2 ----------------
cells_fig <- c(std_cells, prog_cells)
fig_an <- c('S', 'F1', 'ES', 'F1A', 'R')
cols <- c(S = '#B2182B', F1 = '#2166AC', ES = '#4D9221', F1A = '#E08214', R = 'grey40')
pchs <- c(S = 17, F1 = 16, ES = 15, F1A = 18, R = 4)
png('docs/figures/56-fig2-size-power.png', width = 2600, height = 1350, res = 220)
layout(matrix(1:2, 1), widths = c(1.35, 1))
par(mar = c(7.6, 4.4, 2.6, 0.6), mgp = c(2.5, 0.7, 0), las = 1)
x <- seq_along(cells_fig)
plot(NA, xlim = c(0.5, length(x) + 0.5), ylim = c(0, 0.13), xaxt = 'n', xlab = '',
     ylab = 'Null rejection rate', main = 'Size, 1,000 null trials per cell')
rect(0.5, 0.05 - 2 * 0.0069, length(x) + 0.5, 0.05 + 2 * 0.0069, col = 'grey92', border = NA)
abline(h = 0.05, col = 'grey60'); abline(h = 0.075, col = 'grey40', lty = 2)
abline(v = length(std_cells) + 0.5, col = 'grey70', lty = 3)
off <- seq(-0.28, 0.28, length.out = length(fig_an))
for (j in seq_along(fig_an)) {
  y <- size_w[match(cells_fig, label)][[fig_an[j]]]
  points(x + off[j], pmin(y, 0.13), pch = pchs[fig_an[j]], col = cols[fig_an[j]], cex = 1.1, lwd = 2)
}
axis(1, at = x, labels = cells_fig, las = 2, cex.axis = 0.75)
legend('topleft', legend = fig_an, pch = pchs[fig_an], col = cols[fig_an], ncol = 5,
       bty = 'n', cex = 0.8, pt.lwd = 2)
text(length(x) + 0.4, 0.078, 'cap 0.075', adj = c(1, 0), cex = 0.7, col = 'grey30')
crv <- pw[label == 'curve-C07' & analysis %in% fig_an]
par(mar = c(7.6, 4.4, 2.6, 0.8))
plot(NA, xlim = c(-0.3, 0.45), ylim = c(0, 1), xlab = 'c_bm (AR(1) construct, N = 70)',
     ylab = 'Size-adjusted power', main = 'Power curve, 500 trials per point')
abline(h = 0.05, col = 'grey60')
for (a in fig_an) {
  z <- crv[analysis == a][order(cbm)]
  ## no null point on the stored curve, so the line joins c_bm > 0 only
  lines(z[cbm > 0]$cbm, z[cbm > 0]$power_adj, col = cols[a], lwd = 2)
  points(z$cbm, z$power_adj, col = cols[a], pch = pchs[a], cex = 0.9, lwd = 2)
}
dev.off()
cat('\nfigure written\n')
