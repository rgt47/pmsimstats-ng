## docs/58, a modest Monte Carlo check of the closed-form results for the
## drug main effect in the Hybrid design (refined process: anchored
## carryover, randomized allocation, separable AR(1) at rho = 0.7,
## N = 70). Not a stress test. Four cells: the sharp null; a weak null
## with participant-by-treatment SD tau = 5.49; and two alternatives at
## f_m = 0.4 with true t1/2 = 0.2 and 1 week. Per trial: F1 and VF at
## every grid half-life (REML, CR2), the AIC choice with the tested term,
## AIC-weight averaging on the beta_D^exp scale, max-|t| with Bonferroni
## and multivariate-normal calibration, and exact randomization tests
## (199 relabelings, fixed path counts) with label-invariant
## null-model weights.
##
## Run from the repository root:
##   Rscript analysis/scripts/main-effect-half-life/58-mc-check.R [--cores K] [--reps-null R0] [--reps-alt R1] [--summarize-only]
## Writes analysis/data/quick-sim/main-effect-half-life/58-mc-*.csv and
## checkpointed chunks in 58-mc-chunks/.
source('analysis/scripts/main-effect-half-life/58-common.R')
args <- commandArgs(trailingOnly = TRUE)
iarg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
cores <- iarg('--cores', 7L); R0 <- iarg('--reps-null', 1000L); R1 <- iarg('--reps-alt', 500L)
CH <- 25L
load_driver('anchored', 'random')
chunk_dir <- file.path(OUT, '58-mc-chunks'); dir.create(chunk_dir, showWarnings = FALSE)
t_start <- Sys.time()

cells58 <- data.table(cell = 1:4, label = c('sharp null', 'weak null', 'alt, t1/2 = 0.2', 'alt, t1/2 = 1'),
                      f_m = c(0, 0, 0.4, 0.4), tau = c(0, 5.49, 0, 0), t_half = c(0.2, 0.2, 0.2, 1),
                      reps = c(R0, R0, R1, R1))
jobs <- cells58[, .(chunk = seq_len(ceiling(reps / CH))), by = cell]

one_trial <- function(cl, f_m, tau) {
  Ns <- alloc(list(N = 70L), NULL)
  d <- make_trial_me(cl, f_m, tau, Ns)
  rows <- list()
  for (ty in c('F1', 'VF')) {
    G <- grid_me(d, form = FMS[[ty]])
    pr <- procedures(G)
    rows[[ty]] <- list(grid = G$g[, .(type = ty, h, est, se, t, df, p, dbar)],
                       proc = pr[, type := ty])
  }
  nv <- null_V(d)
  S <- perm_setup(d, nv$V)
  rt <- rand_test(S, 199)
  list(grid = rbindlist(lapply(rows, `[[`, 'grid')), proc = rbindlist(lapply(rows, `[[`, 'proc')),
       rand = as.data.table(as.list(c(rt$p, aic_h_F1 = rt$obs[['F1_aic_h']]))),
       Ns = paste(Ns, collapse = '/'))
}

run_job <- function(j) {
  cid <- jobs$cell[j]; ch <- jobs$chunk[j]
  f <- file.path(chunk_dir, sprintf('cell%d-chunk%03d.rds', cid, ch))
  if (file.exists(f)) return(readRDS(f))
  c0 <- cells58[cell == cid]
  cl <- mk_cl(c0$t_half)
  r_ids <- ((ch - 1) * CH + 1):min(ch * CH, c0$reps)
  t0 <- Sys.time()
  out <- lapply(r_ids, function(r) {
    set.seed(58000000L + 100000L * cid + r)
    x <- one_trial(cl, c0$f_m, c0$tau)
    list(grid = x$grid[, `:=`(cell = cid, rep = r)], proc = x$proc[, `:=`(cell = cid, rep = r)],
         rand = x$rand[, `:=`(cell = cid, rep = r, Ns = x$Ns)])
  })
  res <- list(grid = rbindlist(lapply(out, `[[`, 'grid')), proc = rbindlist(lapply(out, `[[`, 'proc')),
              rand = rbindlist(lapply(out, `[[`, 'rand')),
              secs = as.numeric(difftime(Sys.time(), t0, units = 'secs')), n = length(r_ids))
  saveRDS(res, f)
  res
}
## --summarize-only: summarize the chunks already on disk for this grid
## of jobs, without running any new chunk (used when a run is capped).
if ('--summarize-only' %in% args) {
  all <- Filter(Negate(is.null), lapply(seq_len(nrow(jobs)), function(j) {
    f <- file.path(chunk_dir, sprintf('cell%d-chunk%03d.rds', jobs$cell[j], jobs$chunk[j]))
    if (file.exists(f)) readRDS(f) else NULL
  }))
} else {
  all <- mclapply(seq_len(nrow(jobs)), run_job, mc.cores = cores, mc.preschedule = FALSE)
}
grid <- rbindlist(lapply(all, `[[`, 'grid')); proc <- rbindlist(lapply(all, `[[`, 'proc'))
rand <- rbindlist(lapply(all, `[[`, 'rand'))
fwrite(grid, file.path(OUT, '58-mc-grid.csv')); fwrite(proc, file.path(OUT, '58-mc-proc.csv'))
fwrite(rand, file.path(OUT, '58-mc-rand.csv'))

mcse <- function(p, n) sqrt(p * (1 - p) / n)
## Rejection rates at 0.05.
rr_grid <- grid[, .(R = .N, rej = mean(p < 0.05, na.rm = TRUE), mean_est = mean(est, na.rm = TRUE),
                    sd_est = sd(est, na.rm = TRUE), mean_se = mean(se, na.rm = TRUE),
                    mean_df = mean(df, na.rm = TRUE), mean_est_exp = mean(est * dbar, na.rm = TRUE)),
                by = .(cell, type, h)]
rr_grid[, mcse := mcse(rej, R)]
rr_proc <- proc[, .(R = .N, rej = mean(p < 0.05, na.rm = TRUE), mean_est_exp = mean(est_exp, na.rm = TRUE),
                    sd_est_exp = sd(est_exp, na.rm = TRUE)), by = .(cell, type, method)]
rr_proc[, mcse := mcse(rej, R)]
rr_rand <- melt(rand[, !c('Ns', 'aic_h_F1')], id.vars = c('cell', 'rep'))[
  , .(R = .N, rej = mean(value <= 0.05)), by = .(cell, statistic = variable)]
rr_rand[, mcse := mcse(rej, R)]
aic_choice <- proc[method == 'aic', .N, by = .(cell, type, h)][, share := N / sum(N), by = .(cell, type)]
print(merge(rr_grid, cells58[, .(cell, label)], by = 'cell'))
print(merge(rr_proc, cells58[, .(cell, label)], by = 'cell'))
print(merge(rr_rand, cells58[, .(cell, label)], by = 'cell'))
print(dcast(aic_choice, cell + type ~ h, value.var = 'share'))
fwrite(rr_grid, file.path(OUT, '58-mc-summary-grid.csv'))
fwrite(rr_proc, file.path(OUT, '58-mc-summary-procedures.csv'))
fwrite(rr_rand, file.path(OUT, '58-mc-summary-randomization.csv'))
fwrite(aic_choice, file.path(OUT, '58-mc-aic-choice.csv'))
secs <- sum(sapply(all, `[[`, 'secs')); n <- sum(sapply(all, `[[`, 'n'))
tm <- data.table(script = '58-mc-check.R', trials = n, worker_seconds = secs,
                 per_trial = secs / n, wall_seconds = as.numeric(difftime(Sys.time(), t_start, units = 'secs')),
                 cores = cores)
print(tm)
fwrite(tm, file.path(OUT, '58-mc-timing.csv'))
