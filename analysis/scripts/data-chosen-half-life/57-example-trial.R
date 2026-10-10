## docs/57, Example 1: one simulated null Hybrid trial (N = 70), refined
## process of the stress-test driver (anchored carryover, randomized
## allocation, base cell: rho = 0.7, t1/2 = 0.5 week, c_bm = 0). Reports
## the ML profile of the analysis half-life with and without the
## interaction term, the interaction t statistic across the half-life,
## and what each data-chosen half-life procedure reports for this trial,
## including a parametric bootstrap with B = 199.
##
## Run from the repository root:
##   Rscript analysis/scripts/data-chosen-half-life/57-example-trial.R [--cores K] [--boot B]
## Writes analysis/data/quick-sim/data-chosen-half-life/57-example-*.csv
## and docs/figures/57-fig1-example-trial.png.
source('analysis/scripts/data-chosen-half-life/57-common.R')
suppressMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
iarg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
cores <- iarg('--cores', 7L)
B <- iarg('--boot', 199L)
out <- 'analysis/data/quick-sim/data-chosen-half-life'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()

load_driver('anchored', 'random')
cl <- cells[id == 2]
st <- cell_setup(cl)
## The seed was fixed before the trial was looked at; the trial is not
## selected for any property.
set.seed(5701L)
Ns <- alloc(cl, st)
d <- make_trial(cl, st$bs, Ns, 1)
cat('path sizes:', Ns, '\n')

## Fine profile. h = 0 is the binary indicator; 1e6 weeks stands for
## the limit t1/2 -> infinity, in which D_bc is 1 at every visit after
## first exposure.
fine <- c(0, 0.02, 0.05, seq(0.1, 1, by = 0.05), seq(1.1, 3, by = 0.1), 3.5, 4, 5, 6, 8,
          10, 15, 20, 50, 1e6)
t0 <- Sys.time()
prof <- grid_fits(d, fine)
t_prof <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
prof[, `:=`(dev1 = 2 * (max(ll1, na.rm = TRUE) - ll1), dev0 = 2 * (max(ll0, na.rm = TRUE) - ll0))]
fwrite(prof, file.path(out, '57-example-profile.csv'))

## Continuous ML estimate of the half-life under the null model (F1A0's
## criterion) and under F1, on log(h) over [0.01, 100] weeks.
llh <- function(lh, form) ll_or_na(fit_lme(form, with_x(d, exp(lh)), 'ML'))
opt0 <- optimize(llh, c(log(0.01), log(100)), form = FM0, maximum = TRUE)
opt1 <- optimize(llh, c(log(0.01), log(100)), form = FM1, maximum = TRUE)
ll0_bin <- prof[h == 0, ll0]; ll1_bin <- prof[h == 0, ll1]
## Profile-likelihood interval: all h with deviance <= qchisq(0.95, 1),
## read from the fine grid against the larger of the grid and
## continuous maxima.
lmax0 <- max(opt0$objective, prof$ll0, na.rm = TRUE)
lmax1 <- max(opt1$objective, prof$ll1, na.rm = TRUE)
in0 <- prof[2 * (lmax0 - ll0) <= qchisq(0.95, 1), h]
in1 <- prof[2 * (lmax1 - ll1) <= qchisq(0.95, 1), h]
pl <- data.table(model = c('null (F1A0 criterion)', 'F1 (F1A criterion)'),
                 h_ml_continuous = exp(c(opt0$maximum, opt1$maximum)),
                 ll_max = c(lmax0, lmax1),
                 h_ml_fine = c(prof$h[which.max(prof$ll0)], prof$h[which.max(prof$ll1)]),
                 lower = c(min(in0), min(in1)), upper = c(max(in0), max(in1)),
                 dev_at_0 = 2 * (c(lmax0, lmax1) - c(ll0_bin, ll1_bin)),
                 dev_at_inf = 2 * (c(lmax0, lmax1) - c(prof[h == 1e6, ll0], prof[h == 1e6, ll1])))
print(pl)
fwrite(pl, file.path(out, '57-example-profile-interval.csv'))

## All methods on the five-point grid.
t0 <- Sys.time()
set.seed(5702L)
g <- grid_fits(d)
m <- rbind(grid_methods(g), split_methods(d))
t_methods <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
print(g); print(m)
fwrite(g, file.path(out, '57-example-grid.csv'))

## Parametric bootstrap from the null model at F1A0's half-life, each
## resample rerunning every selection (F1AB for F1A; likewise F1A0 and
## the max-|t| statistic).
h0 <- m[method == 'F1A0', h]
nf <- null_fit(d, h0)
t0 <- Sys.time()
RNGkind('L\'Ecuyer-CMRG'); set.seed(5703L)
bs <- rbindlist(mclapply(seq_len(B), function(b) as.list(boot_stats(sim_null(nf, d, h0))),
                         mc.cores = cores, mc.set.seed = TRUE))
t_boot <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
fwrite(bs, file.path(out, '57-example-boot.csv'))
obs <- setNames(abs(m[match(c('F1', 'F1A', 'F1A0', 'MaxT'), method), t]),
                c('F1', 'F1A', 'F1A0', 'MaxT'))
bp <- data.table(statistic = names(obs), observed = obs,
                 crit95 = sapply(names(obs), function(k) quantile(bs[[k]], 0.95, na.rm = TRUE)),
                 p_boot = sapply(names(obs), function(k)
                   (1 + sum(bs[[k]] >= obs[k], na.rm = TRUE)) / (1 + sum(!is.na(bs[[k]])))))
print(bp)
fwrite(bp, file.path(out, '57-example-boot-summary.csv'))
fwrite(m, file.path(out, '57-example-methods.csv'))

tm <- data.table(step = c('fine profile (', 'five-point grid and split', 'bootstrap'),
                 seconds = c(t_prof, t_methods, t_boot),
                 detail = c(sprintf('%d half-lives, 3 fits and one CR2 each', length(fine)),
                            'one worker', sprintf('B = %d on %d workers', B, cores)))
tm[1, step := 'fine profile']
fwrite(tm, file.path(out, '57-example-timing.csv'))
print(tm)

## Figure 1.
source('analysis/scripts/data-chosen-half-life/57-example-figure.R')
cat(sprintf('total %.1f s\n', as.numeric(difftime(Sys.time(), t_start, units = 'secs'))))
