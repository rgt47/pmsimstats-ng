## Monte Carlo check: does the form of the time adjustment in the
## published analysis model matter for the interaction test? The
## published model adjusts the natural-course (TV) and expectancy (PB)
## trends with a single linear t; both are sigmoid in the DGP, and PB
## also steps down when blinding begins. Three forms, fitted to the same
## replicates:
##   linear   Sx ~ bm + Db + t + bm*Db + (1|ptID)            (published)
##   spline   Sx ~ bm + Db + ns(t, 3) + bm*Db + (1|ptID)
##   visit    Sx ~ bm + Db + factor(t) + bm*Db + (1|ptID)
## The biomarker is centered at its sample mean so that the Db
## coefficient is the drug effect at the average biomarker; the bm:Db
## estimate and test are unchanged by centering.
##
## Data: 04-power-simulation.R DGP, CS construct, published allocation
## (Hybrid 18/18/17/17, CO 35/35), N = 70, strength 0.25; rules step
## (arm A), graded (CS matrix, graded coupling), mean (arm D).
##
## With --null the strength is 0 (no interaction; the rules coincide),
## at t1/2 = 0 and 1, so the rejection rates are Type I error rates.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/07-time-adjustment.R \
##     [--reps R] [--cores K] [--null]
suppressMessages({
  library(data.table); library(MASS); library(lmerTest); library(parallel)
  library(splines)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
reps <- arg('--reps', 500L)
cores <- arg('--cores', 8L)
null_run <- '--null' %in% args
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/time-adjustment'
tag <- if (null_run) '-null' else ''
cell_dir <- file.path(out_dir, sprintf('cells-reps%d%s', reps, tag))
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build', 'draw_path')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))
CBM <- if (null_run) 0 else 0.25

## As in 06-single-path-lme.R: the graded rule under CS takes the arm A
## matrix with the graded coupling pattern in the biomarker row.
build_rule <- function(pi, rule, th) {
  if (rule == 'step') return(dgp$build(pi, 'A', CBM, th))
  if (rule == 'mean') return(dgp$build(pi, 'D', CBM, th))
  b <- dgp$build(pi, 'A', 0, th)
  n <- b$n
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, n)
  phi[pi$on | pi$tsd <= 0] <- 0
  idx <- (3 + 2 * n):(2 + 3 * n)
  b$R[1, idx] <- b$R[idx, 1] <- CBM * ifelse(pi$on, 1, phi)
  sds <- sqrt(diag(b$Sigma))
  b$Sigma <- outer(sds, sds) * b$R
  b
}

FORMS <- list(
  linear = Sx ~ bmc + Db + t + bmc * Db + (1 | ptID),
  spline = Sx ~ bmc + Db + ns(t, df = 3) + bmc * Db + (1 | ptID),
  visit = Sx ~ bmc + Db + factor(t) + bmc * Db + (1 | ptID))

fit_one <- function(long, form) {
  tryCatch({
    fit <- suppressMessages(suppressWarnings(lmerTest::lmer(form, data = long)))
    cf <- summary(fit)$coefficients
    c(b_int = cf['bmc:DbTRUE', 'Estimate'], p_int = cf['bmc:DbTRUE', 'Pr(>|t|)'],
      b_db = cf['DbTRUE', 'Estimate'], p_db = cf['DbTRUE', 'Pr(>|t|)'])
  }, error = function(e) c(b_int = NA, p_int = NA, b_db = NA, p_db = NA))
}

run_cell <- function(cell) {
  set.seed(cell$seed)
  d <- dgp$designs[[cell$design]]
  pis <- lapply(d$paths, function(od) dgp$path_info(d$t_wk, d$e, od))
  nP <- length(pis)
  Ns <- N_TOTAL %/% nP + c(rep(1L, N_TOTAL %% nP), rep(0L, nP - N_TOTAL %% nP))
  bs <- lapply(pis, build_rule, rule = cell$rule, th = cell$t_half)
  stopifnot(min(sapply(bs, function(b) min(eigen(b$R, TRUE, TRUE)$values))) > 0)
  draws <- lapply(seq_len(nP), function(k)
    dgp$draw_path(bs[[k]], pis[[k]], CBM, Ns[k] * reps))
  res <- rbindlist(lapply(seq_len(reps), function(r) {
    long <- rbindlist(lapply(seq_len(nP), function(k) {
      rows <- ((r - 1) * Ns[k] + 1):(r * Ns[k])
      dr <- draws[[k]]; pi <- pis[[k]]
      id <- (k - 1) * 1000 + seq_len(Ns[k])
      rbind(data.table(ptID = id, bm = dr$bm[rows], Sx = dr$BL[rows], t = 0, Db = FALSE),
            data.table(ptID = rep(id, times = length(pi$w)),
                       bm = rep(dr$bm[rows], times = length(pi$w)),
                       Sx = as.vector(dr$Y[rows, , drop = FALSE]),
                       t = rep(pi$w, each = Ns[k]), Db = rep(pi$on, each = Ns[k])))
    }))
    long[, bmc := bm - mean(bm[t == 0])]
    rbindlist(lapply(names(FORMS), function(f)
      as.data.table(as.list(c(rep = r, fit_one(long, FORMS[[f]]))))[, form := f]))
  }))
  ## The drug effect the Db coefficient targets at the average biomarker:
  ## the mean BR on-drug minus off-drug difference implied by the DGP,
  ## averaged over visits and paths with the published allocation.
  truth_db <- sum(Ns * vapply(seq_len(nP), function(k) {
    br <- bs[[k]]$mu[(3 + 2 * bs[[k]]$n):(2 + 3 * bs[[k]]$n)]
    on <- pis[[k]]$on
    mean(br[on]) - mean(br[!on])
  }, numeric(1))) / sum(Ns)
  cbind(cell, res[, .(power_int = mean(p_int < 0.05, na.rm = TRUE),
                      mean_b_int = mean(b_int, na.rm = TRUE),
                      sd_b_int = sd(b_int, na.rm = TRUE),
                      power_db = mean(p_db < 0.05, na.rm = TRUE),
                      mean_b_db = mean(b_db, na.rm = TRUE),
                      sd_b_db = sd(b_db, na.rm = TRUE),
                      fails = sum(is.na(p_int)), fits = .N), by = form],
        truth_db_contrast = truth_db)
}

run_cell_checkpointed <- function(cell) {
  f <- file.path(cell_dir, sprintf('cell-%03d.rds', cell$id))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  out <- run_cell(cell)
  saveRDS(out, f)
  cat(sprintf('%s cell %3d %-6s %-6s t1/2=%.1f  %.1f min\n',
              format(Sys.time(), '%H:%M:%S'), cell$id, cell$design, cell$rule,
              cell$t_half, as.numeric(difftime(Sys.time(), t0, units = 'mins'))),
      file = file.path(cell_dir, 'progress.log'), append = TRUE)
  out
}

N_TOTAL <- dgp$N_TOTAL
cells <- if (null_run) {
  CJ(design = c('Hybrid', 'CO'), rule = 'step', t_half = c(0, 1), sorted = FALSE)
} else {
  CJ(design = c('Hybrid', 'CO'), rule = c('step', 'graded', 'mean'),
     t_half = c(0, 1, 2), sorted = FALSE)
}
cells[, id := .I]
cells[, seed := (if (null_run) 20261005L else 20261004L) + id]
cat(sprintf('cells: %d, reps: %d, forms: %d, cores: %d\n',
            nrow(cells), reps, length(FORMS), cores))
t0 <- Sys.time()
out <- mclapply(split(cells, seq_len(nrow(cells))), run_cell_checkpointed,
                mc.cores = cores, mc.preschedule = FALSE)
bad <- which(sapply(out, inherits, 'try-error'))
if (length(bad)) stop('cells failed: ', paste(bad, collapse = ', '))
summ <- rbindlist(out)
summ[, mcse_int := sqrt(power_int * (1 - power_int) / fits)]
fwrite(summ, file.path(out_dir, sprintf('summary-reps%d%s.csv', reps, tag)))
cat(sprintf('elapsed: %.1f min; failed fits: %d of %d\n',
            as.numeric(difftime(Sys.time(), t0, units = 'mins')),
            sum(summ$fails), sum(summ$fits)))
cat('\nInteraction (bm:Db) power by time adjustment\n')
print(dcast(summ, design + rule + t_half ~ form, value.var = 'power_int')[
  , .(design, rule, t_half, linear, spline, visit)], digits = 3)
cat('\nInteraction estimate: mean (sd)\n')
print(dcast(summ[, .(design, rule, t_half, form,
                     v = sprintf('%.4f (%.4f)', mean_b_int, sd_b_int))],
            design + rule + t_half ~ form, value.var = 'v')[
  , .(design, rule, t_half, linear, spline, visit)])
cat('\nDrug effect Db at the average biomarker: mean (sd), and the DGP on-off BR contrast\n')
print(dcast(summ[, .(design, rule, t_half, form, truth = round(truth_db_contrast, 2),
                     v = sprintf('%.2f (%.2f)', mean_b_db, sd_b_db))],
            design + rule + t_half + truth ~ form, value.var = 'v')[
  , .(design, rule, t_half, truth, linear, spline, visit)])
