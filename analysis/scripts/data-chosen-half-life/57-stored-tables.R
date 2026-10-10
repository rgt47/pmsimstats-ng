## docs/57, Tables 2 and 3: F1A's half-life choices by null cell, and
## F1A against F1 on the same stored stress-test trials by the choice.
## Reads the stored replicates; run from the repository root.
suppressMessages(library(data.table))
sd <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'
x2 <- readRDS(file.path(sd, 'replicates-full-extra2.rds'))
x1 <- readRDS(file.path(sd, 'replicates-full.rds'))
a <- x2[analysis == 'F1A' & hyp == 'null']
cat('F1A choice shares by null cell:\n')
tab <- a[, .N, by = .(cell, label, t_half, construct, h)][, share := round(N / sum(N), 3), by = cell]
print(dcast(tab, cell + label + t_half + construct ~ h, value.var = 'share'), nrows = 50)
cat('pooled over all cells F1A:\n')
print(x2[analysis == 'F1A', .N, by = h][order(h)][, share := round(N / sum(N), 3)][])
m <- merge(a[, .(cell, rep, h, pA = p)], x1[analysis == 'F1', .(cell, rep, pF = p)], by = c('cell', 'rep'))
nonprog <- unique(x2[hyp == 'null' & c_tv == 0 & b_tv == 0 & c_pb == 0 & c_bl == 0, cell])
cat('non-prognostic null cells:', nonprog, '\n')
print(m[cell %in% nonprog, .(N = .N, rejF1A = mean(pA < .05, na.rm = TRUE), rejF1 = mean(pF < .05, na.rm = TRUE)), by = h][order(h)])
print(m[cell %in% nonprog, .(rejF1A = mean(pA < .05, na.rm = TRUE), rejF1 = mean(pF < .05, na.rm = TRUE), .N)])
print(m[cell %in% c(2, 20), .(rejF1A = mean(pA < .05, na.rm = TRUE), rejF1 = mean(pF < .05, na.rm = TRUE), .N), by = cell])
print(m[cell %in% c(2, 20) & rep <= 500, .(rejF1A = mean(pA < .05, na.rm = TRUE), rejF1 = mean(pF < .05, na.rm = TRUE), .N), by = cell])
print(x1[analysis %in% c('F1') & cell %in% c(6, 25), .(rej = mean(p < .05, na.rm = TRUE), .N), by = .(cell, analysis)])
print(x2[analysis %in% c('F1A') & cell %in% c(6, 25), .(rej = mean(p < .05, na.rm = TRUE), .N), by = .(cell, analysis)])
print(x1[analysis %in% c('F1') & cell %in% c(6, 25) & rep <= 500, .(rej = mean(p < .05, na.rm = TRUE), .N), by = .(cell, analysis)])
print(x2[analysis %in% c('F1A') & cell %in% c(6, 25) & rep <= 500, .(rej = mean(p < .05, na.rm = TRUE), .N), by = .(cell, analysis)])
