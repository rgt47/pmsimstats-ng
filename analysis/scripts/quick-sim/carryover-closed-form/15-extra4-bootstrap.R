## Bootstrap (trials resampled within cell, paired across analyses) of
## the maximin shortfall over the non-prognostic alternative cells, as
## 12-strawman-stress-decision.R computes it, for job 37 alone ('new')
## and pooled with jobs 15 and 21 ('pooled'). Reads the long files that
## 15-extra4-pooled-inputs.R writes to the work directory (argument 1, or
## EXTRA4_WORK). The reference set is all eight analyses, C2 included,
## so a shortfall here can differ from the decision script's when C2 is
## dropped by the size screen (ES in the pooled set). Post hoc.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args) >= 1) args[1] else Sys.getenv('EXTRA4_WORK')
stopifnot(nzchar(work))
prog <- c('prog-TV-const', 'prog-TV-grow', 'prog-PB', 'prog-BL')
an <- c('S', 'C1', 'C2', 'ES', 'F1', 'F2', 'G1', 'G2')
B <- 300L
for (set in c('new', 'pooled')) {
  x <- fread(file.path(work, sprintf('%s-long.csv', set)))[analysis %in% an]
  x[, p1 := fifelse(is.na(p), 1, p)]
  w <- dcast(x, cell + label + cbm + hyp + rep ~ analysis, value.var = 'p1')
  stat <- function(w) {
    nul <- w[hyp == 'null']; alt <- w[hyp == 'alt']
    crit <- nul[, lapply(.SD, quantile, 0.05), by = label, .SDcols = an]
    alt <- alt[!label %in% prog]
    a <- merge(alt[, c('label', 'cbm', an), with = FALSE], crit, by = 'label', suffixes = c('', '.c'))
    pw <- a[, lapply(an, function(k) mean(get(k) < get(paste0(k, '.c')))), by = .(label, cbm)]
    setnames(pw, c('label', 'cbm', an))
    M <- as.matrix(pw[, ..an]); sf <- apply(M, 1, max) - M
    c(apply(sf, 2, max), mean = colMeans(M))
  }
  obs <- stat(w)
  set.seed(20261015)
  bt <- t(replicate(B, stat(w[, .SD[sample.int(.N, .N, TRUE)], by = cell])))
  mx <- bt[, an]
  cat(sprintf('\n== %s: maximin shortfall, non-prognostic cells (bootstrap SE, B = %d)\n', set, B))
  print(data.table(analysis = an, maximin = obs[an], se = apply(mx, 2, sd))[order(maximin)], digits = 3)
  for (pr in list(c('G2', 'F1'), c('G2', 'G1'), c('G1', 'F1'))) {
    dd <- mx[, pr[2]] - mx[, pr[1]]
    cat(sprintf('%s lead over %s: %.3f (SE %.3f); share of resamples with %s ahead %.2f\n',
                pr[1], pr[2], obs[pr[2]] - obs[pr[1]], sd(dd), pr[1], mean(dd > 0)))
  }
}
