sp <- '/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng-analysis-report-01-dgp-mean-moderation-vs-mvn/10bd8af7-1e83-4060-8bd0-85de612ff8a1/scratchpad'
old_dir <- file.path(sp, 'pre-anchored-outputs/analysis/data/quick-sim/carryover-closed-form')
new_dir <- '/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/data/quick-sim/carryover-closed-form'
k <- c('design', 'construct', 'rule', 't_half', 'variant', 'c.bm')
o <- read.csv(file.path(old_dir, 'e9-closed-form-vs-mc.csv'))
n <- read.csv(file.path(new_dir, 'e9-closed-form-vs-mc.csv'))
m <- merge(o, n, by = k, suffixes = c('.old', '.new'))
m$dmc <- abs(m$mc_power.new - m$mc_power.old)
m$tz <- m$t_half == 0
top <- m[order(-m$dmc), c(k, 'mc_power.old', 'mc_power.new', 'dmc')]
print(head(top, 10))
cat('max stored MC change at t_half > 0, by design:\n')
print(aggregate(dmc ~ design, data = m[!m$tz, ], FUN = max))
cat('max stored MC change at t_half = 0, by design:\n')
print(aggregate(dmc ~ design, data = m[m$tz, ], FUN = max))
cat('Hybrid E9 only, t_half > 0:\n')
h <- m[m$design == 'Hybrid' & m$variant == 'E9' & !m$tz, ]
print(head(h[order(-h$dmc), c(k, 'mc_power.old', 'mc_power.new', 'dmc')], 4))
f <- file.path(sp, 'analyze-01c.out')
if (file.exists(f)) cat(readLines(f), sep = '\n')
