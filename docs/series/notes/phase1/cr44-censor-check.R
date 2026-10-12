suppressMessages(library(data.table))
src <- '/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/data/quick-sim/hendrickson-58b32a9/R'
source(file.path(src, 'utilities.R'))
source(file.path(src, 'censordata.R'))
set.seed(20261011)
labs <- paste0('V', 1:8)
td <- data.table(timeptname = labs, t_wk = rep(2, 8))
n <- 20000
grp <- rep(c(0, 30), each = n / 2)
dat <- data.table(ptID = 1:n)
for (k in seq_along(labs)) {
  dat[, (paste0('D_', labs[k])) := grp * k / 8]
  dat[, (labs[k]) := 50 - grp * k / 8]
}
cp <- data.table(beta0 = 0.05, beta1 = 1 - 1e-9, eb1 = 2)
out <- censordata(copy(dat), td, cp)
cat('biased only (beta1 ~ 1): share censored at last visit\n')
cat('  non-improvers (D = 0):  ', mean(is.na(out[grp == 0, V8])), '\n')
cat('  improvers (D up to 30):', mean(is.na(out[grp == 30, V8])), '\n')
cat('  share censored at first visit, all:', mean(is.na(out$V1)), '\n')
cp2 <- data.table(beta0 = 0.05, beta1 = 0.5, eb1 = 2)
out2 <- censordata(copy(dat), td, cp2)
cat('balanced (0.05, 0.5, 2): censored at V1', mean(is.na(out2$V1)),
  ' at V8', mean(is.na(out2$V8)), '\n')
mono <- out2[, all(diff(is.na(unlist(.SD))) >= 0), by = ptID,
  .SDcols = labs]
cat('monotone censoring for every participant:', all(mono$V1), '\n')
