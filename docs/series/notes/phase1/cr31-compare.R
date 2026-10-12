sp <- '/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng-analysis-report-01-dgp-mean-moderation-vs-mvn/10bd8af7-1e83-4060-8bd0-85de612ff8a1/scratchpad'
old_dir <- file.path(sp, 'pre-anchored-outputs/analysis/data/quick-sim/carryover-closed-form')
new_dir <- '/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/data/quick-sim/carryover-closed-form'

cmp <- function(f, keys, val) {
  o <- read.csv(file.path(old_dir, f))
  n <- read.csv(file.path(new_dir, f))
  m <- merge(o, n, by = keys, suffixes = c('.old', '.new'))
  stopifnot(nrow(m) == nrow(o), nrow(m) == nrow(n))
  m$d <- abs(m[[paste0(val, '.new')]] - m[[paste0(val, '.old')]])
  m
}

summ <- function(m, label) {
  cat('\n==', label, '== rows', nrow(m), '\n')
  a <- aggregate(d ~ design, data = m, FUN = max)
  print(a)
  cat('t_half range:', range(m$t_half), '\n')
  sub <- m[m$t_half <= 2, ]
  cat('restricted t_half <= 2:\n')
  print(aggregate(d ~ design, data = sub, FUN = max))
  top <- m[order(-m$d), ][1:6, ]
  print(top[, intersect(c('design', 'construct', 'rule', 't_half',
    'variant', 'c.bm', 'power.old', 'power.new', 'd'), names(top))])
}

k <- c('design', 'construct', 'rule', 't_half', 'variant', 'c.bm')
m1 <- cmp('e9-closed-form-curves.csv', k, 'power')
summ(m1, 'curves power')
m2 <- cmp('e9-closed-form-curves-designs.csv', k, 'power')
summ(m2, 'curves-designs power')
for (v in c('E9', 'E9b')) {
  cat('\n-- curves-designs, variant', v, '\n')
  print(aggregate(d ~ design, data = m2[m2$variant == v, ], FUN = max))
}
cat('\n-- curves-designs, CO by construct x rule x variant\n')
print(aggregate(d ~ construct + rule + variant,
  data = m2[m2$design == 'CO', ], FUN = max))
m3 <- cmp('e9-closed-form-vs-mc.csv', k, 'power')
summ(m3, 'vs-mc closed power')
m3$dmc <- abs(m3$mc_power.new - m3$mc_power.old)
cat('\nMC power change (stored, includes reseeding) by design:\n')
print(aggregate(dmc ~ design, data = m3, FUN = max))
cat('MC power change excluding Hybrid CS c.bm 0.25 block:\n')
nb <- !(m3$design == 'Hybrid' & m3$construct == 'CS' & m3$c.bm == 0.25)
print(aggregate(dmc ~ design, data = m3[nb, ], FUN = max))
m4 <- cmp('design-power.csv', intersect(names(read.csv(file.path(new_dir,
  'design-power.csv'))), c('design', 'construct', 'rule', 't_half',
  'variant', 'c.bm', 'rho')), 'power')
cat('\ndesign-power.csv power change by design:\n')
print(aggregate(d ~ design, data = m4, FUN = max))
cat('design-power t_half values:', sort(unique(m4$t_half)), '\n')
pm <- cmp('e9-per-path.csv', intersect(names(read.csv(file.path(new_dir,
  'e9-per-path.csv'))), c('design', 'construct', 'rule', 't_half',
  'variant', 'c.bm', 'path')), 'mu')
cat('\nper-path mu change max by design:\n')
print(aggregate(d ~ design, data = pm, FUN = max))
