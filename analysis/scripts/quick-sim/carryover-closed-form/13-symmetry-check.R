## Check of the power-curve symmetry assumed by docs/48: reruns the
## +/- coupling cells of the stress test (11-strawman-stress-test.R)
## with fresh seeds and 1,000 replicates, fitting the strawman S and
## its freed-baseline version ES. Theory says power at -c_bm equals
## power at +c_bm exactly (docs/48, Section 5).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/13-symmetry-check.R [--reps R] [--cores K]
suppressMessages({
  library(data.table); library(MASS); library(mmrm); library(nlme); library(parallel)
})
drv <- 'analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R'
exprs <- parse(drv)
stop_at <- which(vapply(exprs, function(e) startsWith(paste(deparse(e), collapse = ''),
                                                      'setups <-'), logical(1)))[1]
for (e in exprs[seq_len(stop_at - 1)]) eval(e, globalenv())
reps <- arg('--reps', 1000L)
cores <- arg('--cores', 8L)
sym_cells <- rbind(cell('curve-C07', -0.30), cell('curve-C07', 0.30),
                   cell('curve-CS', -0.20, construct = 'CS', rho = 0.8),
                   cell('curve-CS', 0.20, construct = 'CS', rho = 0.8))
sym_cells[, id := .I]
CH <- 50L
jobs <- sym_cells[, .(chunk = seq_len(reps / CH)), by = id][sym_cells, on = 'id']
jobs[, seed := 20261107L + 1000L * id + chunk]
fit_two <- function(d) {
  d05 <- with_x(d, TH_A)
  rbind(as.data.table(as.list(cr2(fit_lme(F_S, d05))))[, analysis := 'S'],
        as.data.table(as.list(cr2(fit_lme(F_E, d05))))[, analysis := 'ES'])
}
out <- rbindlist(mclapply(split(jobs, seq_len(nrow(jobs))), function(j) {
  set.seed(j$seed)
  st <- cell_setup(j)
  rbindlist(lapply(seq_len(CH), function(r)
    fit_two(make_trial(j, st$bs, st$Ns, r))))[, `:=`(label = j$label, cbm = j$cbm)]
}, mc.cores = cores, mc.preschedule = FALSE))
summ <- out[, .(power = mean(!is.na(p) & p < 0.05), mean_est = mean(est, na.rm = TRUE),
                emp_sd = sd(est, na.rm = TRUE), n = .N), by = .(label, cbm, analysis)]
summ[, mcse := sqrt(power * (1 - power) / n)]
dir.create(file.path(out_dir, 'symmetry'), showWarnings = FALSE)
fwrite(summ, file.path(out_dir, 'symmetry', sprintf('summary-reps%d.csv', reps)))
print(summ[order(analysis, label, cbm)], digits = 3)
