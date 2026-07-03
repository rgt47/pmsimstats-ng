## Study B (recovery, balanced-placebo design): does the belief-covariate
## decomposition recover the unbiased biomarker-by-BR slope under a
## contaminated biomarker, when the design decouples drug exposure from
## belief?
##
## Paper 06 round-2 referee point R1. Earlier work (run-study-b-recovery.R)
## showed that on the standard hybrid OL+BDC design NO tractable
## decomposition recovers beta_bm^BR, because the drug state (Dbc) and the
## expectancy (De) are collinear there. This driver tests the
## balanced-placebo / open-hidden design, which decouples them:
##
##   OL  (3 timepoints): on drug,  full belief (e = 1)
##   COV (6 timepoints): on drug,  no belief   (e = 0)  -- covert/hidden
##   OLP (6 timepoints): off drug, full belief (e = 1)  -- open placebo
##
## In the MVN architecture the biomarker-BR coupling lives on Dbc and the
## biomarker-PB coupling scales with De; decoupling Dbc from De lets the
## bm:Dbc + bm:De decomposition isolate beta_bm^BR.
##
## Grid restricted to the PD-feasible contamination range: at c.bm.pb >=
## ~0.45 the covariance is non-positive-definite (min eigenvalue < -0.3)
## and the DGP is no longer faithful; see /tmp pd-check.
##
## USAGE (from repo root):
##   nohup Rscript \
##     analysis/scripts/component-decomposition/run-study-b-balanced.R \
##     [n_reps] > /tmp/study-b-balanced.log 2>&1 &

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table); library(nlme); library(parallel)
})

args        <- commandArgs(trailingOnly = TRUE)
N_REPS      <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
STUDY_SEED  <- 916L
CHECKPOINT_DIR <- file.path(
  'analysis', 'data', 'derived_data', 'component-decomposition',
  'study-b-balanced')
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
TRUE_BETA <- -0.45 * 5 / 1            # -c.bm * sigma_br / sigma_bm

## ------------------------------------------------------------------ ##
## Balanced-placebo design + parameter constructors
## ------------------------------------------------------------------ ##

td_balanced <- function()
  buildtrialdesign(
    name_longform  = 'balanced-placebo (OL/covert/open-placebo)',
    name_shortform = 'BAL',
    timepoints   = seq(1.5, 22.5, by = 1.5),
    timeptnames  = c(paste0('OL', 1:3), paste0('COV', 1:6),
                     paste0('OLP', 1:6)),
    expectancies = c(rep(1, 3), rep(0, 6), rep(1, 6)),
    ondrug       = list(p = c(rep(1, 3), rep(1, 6), rep(0, 6))))

make_resp_params <- function(m_PB = 6, m_TV = 0, m_BR = 10.98604)
  data.table(cat = c('tv', 'pb', 'br'), max = c(m_TV, m_PB, m_BR),
             disp = c(5, 5, 5), rate = c(0.42, 0.35, 0.42),
             sd = c(5, 2, 5))
make_bl_params <- function()
  data.table(cat = c('bm', 'BL'), m = c(0, 70), sd = c(1, 10))
make_model_params <- function(N = 150, c_bm = 0.45, c_bm_pb = 0)
  data.table(N = N, c.bm = c_bm, c.bm.pb = c_bm_pb, carryover_t1half = 1,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.1, c.cfct = 0.05)

## ------------------------------------------------------------------ ##
## Analyses: lumped one-component vs belief-covariate decomposition
## ------------------------------------------------------------------ ##

.prep <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE, full_model_out = TRUE,
             simplecarryover = FALSE, carryover_t1half = 1,
             carryover_scalefactor = 1)
  p <- tryCatch(lme_analysis(td$trialpaths, dat, op), error = function(e) NULL)
  if (is.null(p) || is.null(p$datamerged)) return(NULL)
  d <- copy(p$datamerged)
  d <- d[!is.na(Sx) & !is.na(bm) & !is.na(t) & !is.na(Dbc) & !is.na(De)]
  d[t > 0]
}
.na <- function(tag) data.table(beta = NA_real_, betaSE = NA_real_,
  p = NA_real_, converged = FALSE, formula_dropped = tag)
.get <- function(fit) {
  if (is.null(fit)) return(.na('fit_failed'))
  ct <- summary(fit)$tTable
  tg <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(ct))
  if (!length(tg)) return(.na('no_bmDbc'))
  data.table(beta = unname(ct[tg[1], 'Value']),
             betaSE = unname(ct[tg[1], 'Std.Error']),
             p = unname(ct[tg[1], 'p-value']),
             converged = TRUE, formula_dropped = 'none')
}
CTL <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)

lumped_fit <- function(d)
  .get(tryCatch(lme(Sx ~ bm + t + Dbc + bm:Dbc, random = ~1 | ptID,
    correlation = corCAR1(form = ~t | ptID), data = d, control = CTL),
    error = function(e) NULL))
decomp_fit <- function(d) {
  if (length(unique(d$De)) < 2 || length(unique(round(d$Dbc, 3))) < 2)
    return(.na('no_DbcDe_var'))
  .get(tryCatch(lme(Sx ~ bm + t + Dbc + De + bm:Dbc + bm:De,
    random = ~1 | ptID, correlation = corCAR1(form = ~t | ptID),
    data = d, control = CTL), error = function(e) NULL))
}

ANALYSES <- c('lumped', 'decomposition')

one_rep <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- td_balanced()
  mp <- make_model_params(N = cell$N, c_bm = 0.45, c_bm_pb = cell$c_bm_pb)
  dat <- tryCatch(generateData(mp, make_resp_params(), make_bl_params(),
    td$trialpaths[[1]], empirical = FALSE, makePositiveDefinite = TRUE,
    dgp_architecture = 'mvn'), error = function(e) NULL)
  if (is.null(dat))
    return(rbindlist(lapply(ANALYSES, function(a)
      cbind(analysis = a, rep_idx = rep_idx, .na('gen_failed')))))
  dat[, path := 1]
  d <- .prep(td, dat)
  if (is.null(d))
    return(rbindlist(lapply(ANALYSES, function(a)
      cbind(analysis = a, rep_idx = rep_idx, .na('prep_failed')))))
  rbindlist(list(
    cbind(analysis = 'lumped',        rep_idx = rep_idx, lumped_fit(d)),
    cbind(analysis = 'decomposition', rep_idx = rep_idx, decomp_fit(d))),
    fill = TRUE)
}

run_cell <- function(cell, n_reps, study_seed, cell_id, n_cores) {
  chunks <- mclapply(seq_len(n_reps),
    function(r) one_rep(r, cell, study_seed, cell_id), mc.cores = n_cores)
  good <- !vapply(chunks, function(x)
    inherits(x, 'try-error') || is.null(x), logical(1))
  out <- rbindlist(chunks[good], fill = TRUE)
  out[, `:=`(c_bm_pb = cell$c_bm_pb, N = cell$N, cell_id = cell_id)]
  out
}

summarise_cell <- function(d, true_beta) {
  conv <- d[converged == TRUE]; n <- nrow(conv)
  if (n == 0) return(data.table(n_reps = nrow(d), n_converged = 0L,
    conv_rate = 0, mean_beta = NA_real_, bias = NA_real_,
    mcse_bias = NA_real_, power = NA_real_, coverage = NA_real_))
  mb <- mean(conv$beta); sb <- sd(conv$beta)
  lo <- conv$beta - 1.96 * conv$betaSE; hi <- conv$beta + 1.96 * conv$betaSE
  data.table(n_reps = nrow(d), n_converged = n,
    conv_rate = mean(d$converged), mean_beta = mb, bias = mb - true_beta,
    mcse_bias = sb / sqrt(n), power = mean(conv$p < 0.05),
    coverage = mean(true_beta >= lo & true_beta <= hi))
}

## ------------------------------------------------------------------ ##
## Grid (PD-feasible contamination range) + checkpointed main loop
## ------------------------------------------------------------------ ##

cells <- CJ(c_bm_pb = c(0, 0.10, 0.20, 0.30), N = c(70L, 150L),
            sorted = FALSE)
cells[, cell_id := .I]
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf('Study B balanced: %d cells x %d reps x %d analyses | cores=%d\nOutput: %s\n',
            nrow(cells), N_REPS, length(ANALYSES), n_cores, CHECKPOINT_DIR))

t0 <- proc.time()[['elapsed']]
for (i in seq_len(nrow(cells))) {
  cid <- cells$cell_id[i]
  ck  <- file.path(CHECKPOINT_DIR, sprintf('cell-%03d.rds', cid))
  if (file.exists(ck)) { cat(sprintf('  [skip] cell %d/%d\n', i, nrow(cells))); next }
  cell <- as.list(cells[i])
  cat(sprintf('  cell %d/%d (c_bm_pb=%.2f N=%d) ...\n', i, nrow(cells),
              cell$c_bm_pb, cell$N))
  saveRDS(run_cell(cell, N_REPS, STUDY_SEED, cid, n_cores), ck)
}

all_reps <- rbindlist(lapply(cells$cell_id, function(cid) {
  ck <- file.path(CHECKPOINT_DIR, sprintf('cell-%03d.rds', cid))
  if (file.exists(ck)) readRDS(ck) else NULL }), fill = TRUE)
saveRDS(all_reps, file.path(CHECKPOINT_DIR, 'study-b-balanced-all-reps.rds'))

summary_dt <- all_reps[, summarise_cell(.SD, TRUE_BETA),
                       by = .(analysis, c_bm_pb, N)]
saveRDS(summary_dt, file.path(CHECKPOINT_DIR, 'study-b-balanced-summary.rds'))
write.table(summary_dt,
  file.path(CHECKPOINT_DIR, 'study-b-balanced-summary.txt'),
  sep = '\t', quote = FALSE, row.names = FALSE)

cat(sprintf('\nDone in %.1f min. True beta = %.3f\n',
            (proc.time()[['elapsed']] - t0) / 60, TRUE_BETA))
cat('\n=== bias by analysis x c_bm_pb x N ===\n')
print(summary_dt[order(N, c_bm_pb, analysis),
  .(analysis, c_bm_pb, N, bias = round(bias, 3),
    mcse = round(mcse_bias, 3), cov = round(coverage, 3),
    conv = round(conv_rate, 3))])
