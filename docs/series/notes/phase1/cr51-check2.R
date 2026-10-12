suppressMessages(library(data.table))
d <- '/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/data/quick-sim/carryover-closed-form/strawman-stress'
m <- readRDS(file.path(d, 'replicates-full.rds'))
m[, failed := is.na(est) | is.na(p) | is.na(se)]
cat('cells 40, 41 labels:', unique(m[cell %in% 40:41]$label), '\n')
a <- m[cell <= 39]
cat('cells 1-39: trials', uniqueN(a[, .(cell, rep)]), ' fits', nrow(a),
  ' failed', sum(a$failed), '\n')
cat('cells 40-41 failed:', sum(m[cell >= 40]$failed), '\n')
print(m[failed == TRUE, .N, by = .(cell, label, analysis)][order(cell)])
