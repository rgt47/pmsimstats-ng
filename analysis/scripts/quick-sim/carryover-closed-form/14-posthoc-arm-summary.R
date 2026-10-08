## Summarizes a post hoc arm of 11-strawman-stress-test.R (--extra or
## --extra2) against the main run: pairing check on S, size in the
## standard and prognostic null cells, and size-adjusted power alongside
## S, F1 and C1 from the main run.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/14-posthoc-arm-summary.R [extra|extra2]
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args)) args[1] else 'extra'
d <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'
main <- readRDS(file.path(d, 'replicates-full.rds'))
ex <- readRDS(file.path(d, sprintf('replicates-full-%s.rds', tag)))
## Pairing check: S must agree replicate for replicate.
chk <- merge(main[analysis == 'S', .(cell, rep, est_main = est)],
             ex[analysis == 'S', .(cell, rep, est_ex = est)], by = c('cell', 'rep'))
cat(sprintf('pairing check on S: %d replicates, max |difference| %.2e\n', nrow(chk),
            max(abs(chk$est_main - chk$est_ex), na.rm = TRUE)))
new <- setdiff(unique(ex$analysis), 'S')
all <- rbind(main[analysis %in% c('S', 'F1', 'C1', 'ES', 'DS')],
             ex[analysis %in% new], fill = TRUE)
all <- merge(all[, !c('label', 'cbm', 'hyp'), with = FALSE],
             unique(main[, .(cell, label, cbm, hyp)]), by = 'cell')
all[, rej := !is.na(p) & p < 0.05]
all[, p1 := fifelse(is.na(p), 1, p)]
prog <- c('prog-TV-const', 'prog-TV-grow', 'prog-PB', 'prog-BL')
nul <- all[hyp == 'null', .(reject = mean(rej), n = .N), by = .(label, analysis)]
std <- nul[!label %in% prog, .(pooled = sum(reject * n) / sum(n), max_cell = max(reject),
                              worst = label[which.max(reject)]), by = analysis]
cat('\nSize, standard null cells (pooled band 0.046 to 0.054; cell limit 0.075)\n')
print(std, digits = 3)
cat('\nSize, prognostic null cells\n')
print(dcast(nul[label %in% prog], analysis ~ label, value.var = 'reject'), digits = 3)
crit <- all[hyp == 'null', .(crit = quantile(p1, 0.05)), by = .(label, analysis)]
alt <- merge(all[hyp == 'alt'], crit, by = c('label', 'analysis'))
pw <- alt[, .(power_adj = mean(p1 < crit)), by = .(label, cbm, analysis)]
cat('\nSize-adjusted power, alternative cells\n')
print(dcast(pw, label + cbm ~ analysis, value.var = 'power_adj'), digits = 3, nrows = 100)
fwrite(merge(std, dcast(nul[label %in% prog], analysis ~ label, value.var = 'reject'),
             by = 'analysis', all = TRUE), file.path(d, sprintf('posthoc-%s-size.csv', tag)))
fwrite(pw, file.path(d, sprintf('posthoc-%s-power.csv', tag)))
