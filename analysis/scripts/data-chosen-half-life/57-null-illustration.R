## docs/57, Example 2: a small Monte Carlo illustration of the
## data-chosen half-life procedures on stress-test cells. Not a
## replacement for the planned stress-test extension.
##
## Trials are those of the stress-test driver 11-strawman-stress-test.R,
## seed for seed (chunks of 25, seed 20261007 + 1000 * cell + chunk), so
## replicate r of a cell here is replicate r of the stored stress test
## when run with '--carryover recursive --allocation fixed' (the stored
## setting, and the default here). The split and the one bootstrap
## resample per replicate (warp-speed method) use a separate stream.
##
## Run from the repository root:
##   Rscript analysis/scripts/data-chosen-half-life/57-null-illustration.R \
##     [--cells 2,20,6,25] [--reps 500,500,250,250] [--cores 7] \
##     [--carryover recursive|anchored] [--allocation fixed|random]
## Writes analysis/data/quick-sim/data-chosen-half-life/57-mc-*.{rds,csv}
## and docs/figures/57-fig2-null-illustration.png. Chunks are
## checkpointed; a relaunch with the same arguments resumes.
source('analysis/scripts/data-chosen-half-life/57-common.R')
args <- commandArgs(trailingOnly = TRUE)
sarg <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
cell_ids <- as.integer(strsplit(sarg('--cells', '2,20,6,25'), ',')[[1]])
reps <- as.integer(strsplit(sarg('--reps', '500,500,250,250'), ',')[[1]])
stopifnot(length(reps) == length(cell_ids), all(reps %% 25 == 0))
cores <- as.integer(sarg('--cores', '7'))
carry <- sarg('--carryover', 'recursive'); allo <- sarg('--allocation', 'fixed')
load_driver(carry, allo)
tag <- if (carry == 'recursive' && allo == 'fixed') '' else paste0('-', carry, '-', allo)
out <- 'analysis/data/quick-sim/data-chosen-half-life'
chunk_dir <- file.path(out, paste0('57-mc-chunks', tag))
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

setups <- lapply(cell_ids, function(i) cell_setup(cells[id == i]))
names(setups) <- cell_ids
jobs <- rbindlist(lapply(seq_along(cell_ids), function(j)
  data.table(id = cell_ids[j], chunk = seq_len(reps[j] / 25))))

run_chunk <- function(job) {
  f <- file.path(chunk_dir, sprintf('cell-%03d-chunk-%02d.rds', job$id, job$chunk))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  cl <- cells[id == job$id]
  st <- setups[[as.character(job$id)]]
  set.seed(20261007L + 1000L * job$id + job$chunk)
  res <- lapply(seq_len(25), function(r) {
    rr <- (job$chunk - 1L) * 25L + r
    d <- make_trial(cl, st$bs, alloc(cl, st), r)
    a <- analyze_trial(d, 570000000L + 10000L * job$id + rr, n_boot = 1L)
    list(methods = a$methods[, `:=`(cell = job$id, rep = rr)],
         boot = if (is.null(a$boot)) NULL else a$boot[, `:=`(cell = job$id, rep = rr)],
         grid = a$grid[, `:=`(cell = job$id, rep = rr)])
  })
  o <- list(methods = rbindlist(lapply(res, `[[`, 'methods')),
            boot = rbindlist(lapply(res, `[[`, 'boot')),
            grid = rbindlist(lapply(res, `[[`, 'grid')),
            secs = as.numeric(difftime(Sys.time(), t0, units = 'secs')))
  saveRDS(o, f)
  o
}
t0 <- Sys.time()
res <- mclapply(split(jobs, seq_len(nrow(jobs))), function(j)
  tryCatch(run_chunk(j), error = function(e) conditionMessage(e)),
  mc.cores = cores, mc.preschedule = FALSE)
bad <- which(!vapply(res, is.list, logical(1)))
if (length(bad)) stop('chunks failed: ', paste(unlist(res[bad]), collapse = ' | '))
wall <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))
M <- rbindlist(lapply(res, `[[`, 'methods'))
BS <- rbindlist(lapply(res, `[[`, 'boot'))
GR <- rbindlist(lapply(res, `[[`, 'grid'))
secs <- sum(vapply(res, `[[`, numeric(1), 'secs'))
saveRDS(list(methods = M, boot = BS, grid = GR), file.path(out, sprintf('57-mc-replicates%s.rds', tag)))
cinfo <- cells[id %in% cell_ids, .(cell = id, label, cbm, hyp, t_half)]

## Nominal rejection rates.
nom <- M[, .(reject = mean(p < 0.05, na.rm = TRUE), fails = sum(is.na(p)), R = .N),
         by = .(cell, method)]
nom[, mcse := sqrt(reject * (1 - reject) / (R - fails))]
nom[, calibration := 'nominal']
## Warp-speed bootstrap (Giacomini, Politis and White 2013): one
## null-model resample per replicate; the critical value is the 95%
## quantile of the pooled bootstrap statistics of the cell.
crit <- melt(BS, id.vars = c('cell', 'rep'), variable.name = 'method', value.name = 'tstar')
crit <- crit[, .(crit = quantile(tstar, 0.95, na.rm = TRUE)), by = .(cell, method)]
crit[, method := as.character(method)]
ws <- merge(M[method %in% crit$method], crit, by = c('cell', 'method'))
ws <- ws[, .(reject = mean(abs(t) > crit, na.rm = TRUE), fails = sum(is.na(t)), R = .N,
             crit = crit[1]), by = .(cell, method)]
ws[, mcse := sqrt(reject * (1 - reject) / (R - fails))]
ws[, `:=`(calibration = 'warp-speed bootstrap', method = paste0(method, 'B'))]
summ <- merge(rbind(nom, ws, fill = TRUE), cinfo, by = 'cell')
fwrite(summ, file.path(out, sprintf('57-mc-summary%s.csv', tag)))
print(dcast(summ, method ~ label + cbm, value.var = 'reject'), digits = 3)

## Distribution of the chosen half-life.
hs <- M[method %in% c('F1A', 'F1A0', 'Split'), .N, by = .(cell, method, h)]
hs[, share := N / sum(N), by = .(cell, method)]
fwrite(hs, file.path(out, sprintf('57-mc-choices%s.csv', tag)))
## Rejection by chosen half-life (selection dependence).
byh <- M[method %in% c('F1A', 'F1A0'), .(reject = mean(p < 0.05, na.rm = TRUE), .N),
         by = .(cell, method, h)]
fwrite(byh, file.path(out, sprintf('57-mc-reject-by-choice%s.csv', tag)))
## Agreement of the two selections, and of each with the fixed value.
ag <- dcast(M[method %in% c('F1', 'F1A', 'F1A0')], cell + rep ~ method, value.var = 'p')
agr <- ag[, .(F1A_eq_F1A0 = mean(abs(F1A - F1A0) < 1e-12, na.rm = TRUE),
              F1A_eq_F1 = mean(abs(F1A - F1) < 1e-12, na.rm = TRUE),
              F1A0_eq_F1 = mean(abs(F1A0 - F1) < 1e-12, na.rm = TRUE)), by = cell]
fwrite(agr, file.path(out, sprintf('57-mc-agreement%s.csv', tag)))

## Pairing check against the stored stress test (stored setting only).
if (tag == '') {
  sd <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'
  s1 <- readRDS(file.path(sd, 'replicates-full.rds'))[analysis == 'F1' & cell %in% cell_ids,
                                                      .(cell, rep, est_s = est, p_s = p)]
  s2 <- readRDS(file.path(sd, 'replicates-full-extra2.rds'))[analysis == 'F1A' & cell %in% cell_ids,
                                                             .(cell, rep, est_s = est, p_s = p, h_s = h)]
  c1 <- merge(M[method == 'F1', .(cell, rep, est, p)], s1, by = c('cell', 'rep'))
  c2 <- merge(M[method == 'F1A', .(cell, rep, est, p, h)], s2, by = c('cell', 'rep'))
  pc <- rbind(c1[, .(what = 'F1', n = .N, max_abs_diff_est = max(abs(est - est_s), na.rm = TRUE),
                     same_decision = mean((p < .05) == (p_s < .05), na.rm = TRUE)), by = cell],
              c2[, .(what = 'F1A', n = .N, max_abs_diff_est = max(abs(est - est_s), na.rm = TRUE),
                     same_decision = mean((p < .05) == (p_s < .05), na.rm = TRUE)), by = cell])
  print(pc)
  fwrite(pc, file.path(out, '57-mc-pairing-check.csv'))
  stored <- rbind(c1[, .(method = 'F1 (stored)', reject = mean(p_s < .05, na.rm = TRUE)), by = cell],
                  c2[, .(method = 'F1A (stored)', reject = mean(p_s < .05, na.rm = TRUE)), by = cell])
  fwrite(stored, file.path(out, '57-mc-stored-rates.csv'))
}
tm <- data.table(replicates = nrow(M[method == 'F1']), worker_seconds = secs,
                 worker_seconds_per_rep = secs / nrow(M[method == 'F1']),
                 wall_minutes = wall, cores = cores)
fwrite(tm, file.path(out, sprintf('57-mc-timing%s.csv', tag)))
print(tm)

## Figure 2: rejection rates with 2-MCSE bars.
suppressMessages(library(ggplot2))
show <- c('F1', 'F1A', 'F1A0', 'MA', 'MaxT', 'Split', 'CrossAvg', 'CrossCauchy',
          'F1B', 'F1AB', 'F1A0B', 'MaxTB')
lab <- c(F1 = 'F1 fixed 0.5', F1A = 'F1A', F1A0 = 'F1A0', MA = 'AIC averaging',
         MaxT = 'max-|t| Bonferroni', Split = 'single split', CrossAvg = 'cross-fit, average',
         CrossCauchy = 'cross-fit, Cauchy', F1B = 'F1, bootstrap', F1AB = 'F1AB',
         F1A0B = 'F1A0, bootstrap', MaxTB = 'max-|t|, bootstrap')
fd <- summ[method %in% show]
fd[, mlab := factor(lab[method], levels = rev(lab[show]))]
fd[, panel := sprintf('%s, c_bm = %s, t1/2 = %s (%s, R = %d)', label, cbm, t_half, hyp, R)]
p <- ggplot(fd, aes(reject, mlab)) +
  geom_vline(data = fd[hyp == 'null'], aes(xintercept = 0.05), colour = 'grey50', linewidth = 0.3) +
  geom_errorbarh(aes(xmin = reject - 2 * mcse, xmax = reject + 2 * mcse), height = 0.25) +
  geom_point(aes(shape = calibration), size = 1.8) +
  facet_wrap(~ panel, scales = 'free_x', ncol = 2) +
  labs(x = 'rejection rate at nominal 0.05 (bars: 2 MCSE)', y = NULL, shape = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = 'bottom')
png(sprintf('docs/figures/57-fig2-null-illustration%s.png', tag), width = 9, height = 6.5,
    units = 'in', res = 200)
print(p)
invisible(dev.off())
