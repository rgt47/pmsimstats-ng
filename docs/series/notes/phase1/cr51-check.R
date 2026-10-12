suppressMessages(library(data.table))
d <- '/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/data/quick-sim/carryover-closed-form/strawman-stress'
m <- readRDS(file.path(d, 'replicates-full.rds'))
e1 <- readRDS(file.path(d, 'replicates-full-extra.rds'))
e2 <- readRDS(file.path(d, 'replicates-full-extra2.rds'))

cat('main: cells', uniqueN(m$cell), ' trials', uniqueN(m[, .(cell, rep)]),
  ' analyses', paste(sort(unique(m$analysis)), collapse = ' '), '\n')
cat('extra analyses:', paste(sort(unique(e1$analysis)), collapse = ' '),
  '; extra2 analyses:', paste(sort(unique(e2$analysis)), collapse = ' '),
  '\n')
stopifnot(nrow(m) == uniqueN(m[, .(cell, rep, analysis)]))

fail <- function(x) is.na(x$est) | is.na(x$p) | is.na(x$se)
m[, failed := fail(.SD)]
e1[, failed := fail(.SD)]
e2[, failed := fail(.SD)]
cat('\nFailures, main run, by analysis (denominator per analysis):\n')
print(m[, .(fits = .N, failed = sum(failed)), by = analysis][order(-failed)])
cat('main run total:', sum(m$failed), 'of', nrow(m), '\n')
cat('\nFailures, extra (G1, G2) and extra2 (F1A), non-S rows:\n')
print(e1[analysis != 'S', .(fits = .N, failed = sum(failed)), by = analysis])
print(e2[analysis != 'S', .(fits = .N, failed = sum(failed)), by = analysis])
cat('S in extra files failed:', sum(e1[analysis == 'S']$failed),
  sum(e2[analysis == 'S']$failed), '\n')
s_same <- merge(m[analysis == 'S', .(cell, rep, est)],
  e1[analysis == 'S', .(cell, rep, est2 = est)], by = c('cell', 'rep'))
cat('S identical across main and extra:', isTRUE(all.equal(s_same$est,
  s_same$est2)), '\n')
all14 <- rbind(m[, .(cell, rep, analysis, failed)],
  e1[analysis != 'S', .(cell, rep, analysis, failed)],
  e2[analysis != 'S', .(cell, rep, analysis, failed)])
cat('\nAll 14 analyses: total failed', sum(all14$failed), 'of', nrow(all14),
  '\n')
pc <- all14[, .(rate = mean(failed)), by = .(cell, analysis)]
cat('max per-cell failure rate:', max(pc$rate), '\n')
print(pc[rate == max(rate)])

cells <- unique(m[, .(cell, label, hyp, t_half, construct, rho, N,
  decay_k, dropout, c_tv, c_pb, c_bl, meas, reps)])
cat('\ncells by hyp:\n'); print(cells[, .N, by = hyp])
cat('cells with reps:\n'); print(cells[, .N, by = .(hyp, reps)])

sel_share <- function(x, an) {
  y <- x[analysis == an & !is.na(h)]
  y[, .(n = .N, h0 = mean(h == 0), h1or2 = mean(h %in% c(1, 2)),
    h05 = mean(h == 0.5)), by = .(cell, label, hyp, t_half, construct)]
}
for (an in c('DS', 'F1A')) {
  src <- if (an == 'DS') m else e2
  cat('\n== selection shares,', an, '==\n')
  cat('grid used:', sort(unique(src[analysis == an]$h)), '\n')
  s <- sel_share(src, an)
  print(s[t_half %in% c(0, 1) | label %like% 'base|CS'][order(hyp, cell)])
  cat('null cells with t_half 0 pooled / 1 pooled:\n')
  y <- src[analysis == an & !is.na(h)]
  print(y[, .(n = .N, h0 = mean(h == 0), h1or2 = mean(h %in% c(1, 2))),
    by = .(hyp, t_half)][order(hyp, t_half)])
}
