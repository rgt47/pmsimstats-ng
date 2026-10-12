## Combine the main, follow-up and Dbc runs of 08-covariance-study.R
## (same cells, same chunk seeds, so the arms are paired) into one table
## of analyses, and print the size and power tables of docs/45.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/09-covariance-study-tables.R \
##     [--tag T] [--reps ALT,NULL] [--dir PATH]
## --tag T (e.g. -anchored) reads the 08 runs whose tag ends in T and
## appends T to the two output files; without it the stored untagged
## runs are read as before. --reps (default 250,500) and --dir allow a
## smoke test on a small run.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
sarg <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
tag <- sarg('--tag', '')
rr <- as.integer(strsplit(sarg('--reps', '250,500'), ',')[[1]])
rt <- sprintf('alt%d-null%d', rr[1], rr[2])
d <- sarg('--dir', 'analysis/data/quick-sim/carryover-closed-form/covariance-study')
runs <- c(main = '', followup = '-followup', dbc = '-dbc')
## A pass that was not run under this tag (e.g. the anchored follow-up
## and Dbc passes, cancelled; DQ-8) is skipped, so the tables then
## cover the remaining passes only.
have <- file.exists(file.path(d, sprintf('summary-%s%s%s.csv', rt, runs, tag)))
if (!have[1]) stop('main pass missing: ', sprintf('summary-%s%s.csv', rt, tag))
if (!all(have)) cat('Passes not found for tag', sQuote(tag), 'and skipped:',
                    paste(names(runs)[!have], collapse = ', '), '\n')
runs <- runs[have]
all <- rbindlist(lapply(names(runs), function(r)
  fread(file.path(d, sprintf('summary-%s%s%s.csv', rt, runs[[r]], tag)))[, run := r]),
  fill = TRUE)
all[is.na(extra), extra := 'none']
all[is.na(expo), `:=`(expo = 'Db', engine = 'mmrm')]
all[, analysis := paste0(
  fifelse(coding == 'covariate', 'BL+', ''),
  fifelse(engine == 'lme', 'RI+CAR1', cov), '-',
  c(Satterthwaite = 'Satt', `Kenward-Roger` = 'KR', CR2 = 'CR2')[test],
  fifelse(extra == 'sepbm', '-sepbm', ''),
  fifelse(expo == 'Dbc', paste0('-Dbc', th_a), ''))]
stopifnot(!anyDuplicated(all[, .(cell, analysis, mean)]))
fwrite(all, file.path(d, sprintf('summary-%s-combined%s.csv', rt, tag)))

## Size: maximum and range of the null rejection rate over the eight null
## cells (2 constructs x 2 designs x 2 half-lives); 500 replicates per cell,
## so 0.05 has Monte Carlo SE 0.0097 and a 95% band of 0.031 to 0.069.
sz <- all[hyp == 'null', .(min_size = min(reject), max_size = max(reject),
                           mean_size = mean(reject), fails = sum(fails),
                           fits = sum(fits)), by = .(analysis, mean)]
sz[, sized := max_size <= 0.069 & min_size >= 0.031]
setorder(sz, -sized, max_size)
cat('\nSize over the 8 null cells\n')
print(sz, digits = 3, nrows = 200)

cat('\nNull rejection by cell\n')
print(dcast(all[hyp == 'null'], analysis + mean ~ construct + design + t_half,
            value.var = 'reject'), digits = 3, nrows = 200)
cat('\nPower by cell\n')
print(dcast(all[hyp == 'alt'], analysis + mean ~ construct + design + t_half,
            value.var = 'reject'), digits = 3, nrows = 200)
cat('\nMean estimate by cell (alternative)\n')
print(dcast(all[hyp == 'alt'], analysis + mean ~ construct + design + t_half,
            value.var = 'mean_est'), digits = 3, nrows = 200)
cat('\nkappa (mean SE / empirical SD), null cells\n')
print(dcast(all[hyp == 'null'], analysis + mean ~ construct + design + t_half,
            value.var = 'kappa'), digits = 3, nrows = 200)
cat('\nfailed fits:', sum(all$fails), 'of', sum(all$fits), '\n')

## Size-adjusted power: reject when p is below the empirical 5% quantile
## of the same analysis's p-values in the matching null cell (same
## construct, design and t1/2). Compares power at a common realized size;
## the null quantile rests on 500 replicates, so it is itself noisy.
reps <- rbindlist(lapply(names(runs), function(r)
  readRDS(file.path(d, sprintf('replicates-%s%s%s.rds', rt, runs[[r]], tag)))), fill = TRUE)
reps[is.na(extra), extra := 'none']
reps[is.na(expo), `:=`(expo = 'Db', engine = 'mmrm')]
reps <- merge(reps, unique(all[, .(cell, construct, design, hyp, t_half, coding, mean, cov,
                                   test, extra, expo, th_a, engine, analysis)]),
              by = c('cell', 'coding', 'mean', 'cov', 'test', 'extra', 'expo', 'th_a', 'engine'))
stopifnot(nrow(reps) == sum(all$fits))
crit <- reps[hyp == 'null', .(p_crit = quantile(p, 0.05, na.rm = TRUE)),
             by = .(construct, design, t_half, analysis, mean)]
adj <- merge(reps[hyp == 'alt'], crit, by = c('construct', 'design', 't_half', 'analysis', 'mean'))
adj <- adj[, .(power = mean(p < 0.05, na.rm = TRUE), power_adj = mean(p < p_crit, na.rm = TRUE)),
           by = .(construct, design, t_half, analysis, mean)]
fwrite(adj, file.path(d, sprintf('size-adjusted-power-%s%s.csv', rt, tag)))
cat('\nSize-adjusted power by cell\n')
print(dcast(adj, analysis + mean ~ construct + design + t_half, value.var = 'power_adj'),
      digits = 3, nrows = 200)
