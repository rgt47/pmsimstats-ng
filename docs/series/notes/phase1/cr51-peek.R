d <- '/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/data/quick-sim/carryover-closed-form/strawman-stress'
for (f in c('replicates-full.rds', 'replicates-full-extra.rds',
  'replicates-full-extra2.rds')) {
  x <- readRDS(file.path(d, f))
  cat('\n==', f, 'class', class(x), '\n')
  if (is.data.frame(x)) {
    cat('dim', dim(x), '\n')
    str(x, list.len = 60, vec.len = 2)
  } else {
    str(x, max.level = 2, list.len = 20)
  }
}
