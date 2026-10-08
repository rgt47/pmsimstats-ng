## Summarize the power simulation of 04-power-simulation.R.
##
## Reports nominal power, Type I error, size-adjusted power (critical
## value from the matching null cell: same arm, design and half-life),
## and the mean interaction estimate; and sets the published arm beside
## Hendrickson's own saved results at matching cells (N = 70).
##
## Size adjustment matters because the published analysis model (random
## intercept, no serial correlation) is correctly specified for CS data
## but not for AR(1) data, so nominal power is not comparable across
## arms.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/05-summarize-power.R [--reps R]
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
reps <- if ('--reps' %in% args) as.integer(args[which(args == '--reps') + 1]) else 500L
analysis <- if ('--analysis' %in% args) args[which(args == '--analysis') + 1] else 'published'
tag <- if (analysis == 'published') '' else paste0('-', analysis)
dir <- 'analysis/data/quick-sim/hendrickson-problems/power-sim'
r <- readRDS(file.path(dir, sprintf('replicates-reps%d%s.rds', reps, tag)))
summ <- fread(file.path(dir, sprintf('summary-reps%d%s.csv', reps, tag)))
## Merge any arm-restricted runs (e.g. arm F) written by 04 --arms.
extra <- list.files(dir, pattern = sprintf('^summary-reps%d%s-[A-F]+\\.csv$', reps, tag))
for (f in extra) {
  sfx <- sub(sprintf('^summary-reps%d%s(-[A-F]+)\\.csv$', reps, tag), '\\1', f)
  summ <- rbind(summ, fread(file.path(dir, f)), fill = TRUE)
  r <- rbind(r, readRDS(file.path(dir, sprintf('replicates-reps%d%s%s.rds', reps, tag, sfx))),
             fill = TRUE)
  cat(sprintf('merged arm-restricted run: %s\n', f))
}
cat(sprintf('analysis model: %s\n', analysis))

crit <- r[c.bm == 0, .(p_crit = quantile(p, 0.05, na.rm = TRUE),
                       type1 = mean(p < 0.05, na.rm = TRUE)),
          by = .(arm, design, t_half)]
alt <- merge(r[c.bm > 0], crit, by = c('arm', 'design', 't_half'))
pw <- alt[, .(power = mean(p < 0.05, na.rm = TRUE),
              power_adj = mean(p < p_crit[1], na.rm = TRUE),
              mean_beta = mean(beta, na.rm = TRUE),
              emp_sd = sd(beta, na.rm = TRUE),
              mean_se = mean(se, na.rm = TRUE)),
          by = .(arm, design, c.bm, t_half)]
fwrite(pw, file.path(dir, sprintf('power-summary-reps%d%s.csv', reps, tag)))
fwrite(crit, file.path(dir, sprintf('type1-reps%d%s.csv', reps, tag)))

arm_lab <- c(A = 'A published', B = 'B adj1', C = 'C adj1+2',
             D = 'D adj3', E = 'E adj1+2+3', F = 'F adj1+2 covar')
show <- function(d, v, title) {
  cat('\n', title, '\n', sep = '')
  w <- dcast(d, arm + design + c.bm ~ t_half, value.var = v)
  w[, arm := arm_lab[arm]]
  print(w, digits = 2)
}
cat(sprintf('replicates per cell: %d; infeasible cells skipped: %d\n',
            reps, sum(!summ$feasible)))
cat('\nInfeasible cells (matrix not positive definite):\n')
print(summ[feasible == FALSE, .(arm, design, c.bm, t_half)][
  , .(t_half = paste(t_half, collapse = ' ')), by = .(arm, design, c.bm)])
show(copy(crit)[, c.bm := 0], 'type1', 'Type I error (columns: half-life)')
show(pw, 'power', 'Nominal power')
show(pw, 'power_adj', 'Size-adjusted power')
show(pw, 'mean_beta', 'Mean interaction estimate')

## Published arm against Hendrickson's saved results (N = 70).
pub_file <- 'analysis/data/quick-sim/hendrickson-problems/published-power.csv'
if (file.exists(pub_file)) {
  pub <- fread(pub_file)[N == 70, .(design = sub('N-of-1', 'Hybrid', design),
                                    c.bm, t_half = carryover_t1half,
                                    published = power)]
  cmp <- merge(pub, rbind(
    pw[arm == 'A', .(design, c.bm, t_half, sim_A = power)],
    crit[arm == 'A', .(design, c.bm = 0, t_half, sim_A = type1)]),
    by = c('design', 'c.bm', 't_half'), all.x = TRUE)
  cat('\nPublished saved results (with repair) vs simulated arm A (no repair; NA = infeasible)\n')
  print(cmp[design != 'OL'], digits = 2)
}
