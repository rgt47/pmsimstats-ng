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
## Parameters: the base response and baseline tables are the active
## set's (analysis/scripts/carryover-sensitivity/simulation-core.R;
## docs/74, D6), with m_PB = 6 and m_TV = 0 overriding the PB and TV
## maxima as before. The default set is the raw tinytest fixture, the
## set the stored run used. The analysis biomarker is standardized by
## the set's mean and SD (docs/74, D1), so the true slope is
## -c_bm * sigma_BR per biomarker SD (-2.25 fixture, -3.6 published).
##
## USAGE (from repo root):
##   nohup Rscript \
##     analysis/scripts/component-decomposition/run-study-b-balanced.R \
##     [n_reps] [out_dir] [--params NAME] > /tmp/study-b-balanced.log 2>&1 &
##   defaults: 1000, study-b-balanced, --params fixture (or PMSIM_PARAMS
##   when set). Under any other set the default out_dir is
##   study-b-balanced-<set>. out_dir is a subdirectory of
##   analysis/data/derived_data/component-decomposition/ unless it is
##   an absolute path; it carries a PARAM_SET marker and refuses a
##   resume under another set. meta.rds records the set and the true
##   value. SIM_CORES (env) sets the worker count; otherwise
##   detectCores() - 1.

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table); library(nlme); library(parallel)
})

## The stored paper 06 runs used the raw tinytest fixture, so this
## driver defaults to 'fixture' rather than the core's 'simcore2026':
## a run without --params then reproduces the stored one.
if (!'--params' %in% commandArgs(trailingOnly = TRUE) &&
    !nzchar(Sys.getenv('PMSIM_PARAMS')))
  Sys.setenv(PMSIM_PARAMS = 'fixture')
source(file.path('analysis', 'scripts', 'carryover-sensitivity',
                 'simulation-core.R'))

## SIM_CORES lets a cgroup-limited job size its pool; detectCores()
## reports the host, not the job.
n_workers <- function() {
  v <- Sys.getenv('SIM_CORES')
  if (nzchar(v)) as.integer(v) else max(1L, detectCores() - 1L)
}

## Any set other than the fixture writes to a sibling folder, so it
## can neither resume from nor overwrite the stored run.
set_dir_name <- function(base)
  if (param_set() == 'fixture') base else paste0(base, '-', param_set())

## guard_out_dir() reads an unmarked folder that holds files as
## 'simcore2026', but the stored paper 06 folders hold raw-fixture
## runs, so the legacy folder is marked 'fixture' before it is guarded.
guard_paper06_dir <- function(dir, legacy_dir) {
  marker <- file.path(dir, 'PARAM_SET')
  if (param_set() == 'fixture' && !file.exists(marker) &&
      dir.exists(dir) && length(list.files(dir)) > 0 &&
      normalizePath(dir) == normalizePath(legacy_dir, mustWork = FALSE))
    writeLines('fixture', marker)
  guard_out_dir(dir)
}

args        <- commandArgs(trailingOnly = TRUE)
i_params    <- which(args == '--params')
if (length(i_params)) args <- args[-c(i_params[1], i_params[1] + 1L)]
N_REPS      <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
OUT_NAME    <- if (length(args) >= 2) args[[2]] else
  set_dir_name('study-b-balanced')
STUDY_SEED  <- 916L
CD_ROOT <- file.path('analysis', 'data', 'derived_data',
                     'component-decomposition')
CHECKPOINT_DIR <- if (startsWith(OUT_NAME, '/')) OUT_NAME else
  file.path(CD_ROOT, OUT_NAME)
guard_paper06_dir(CHECKPOINT_DIR, file.path(CD_ROOT, 'study-b-balanced'))

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

## Base table from the active set (docs/74, D6); m_PB and m_TV override
## only the PB and TV maxima.
make_resp_params <- function(m_PB = 6, m_TV = 0) {
  rp <- as.data.table(default_resp_param())
  rp[cat == 'pb', max := m_PB]
  rp[cat == 'tv', max := m_TV]
  rp
}
make_bl_params <- function() as.data.table(default_baseline_param())
make_model_params <- function(N = 150, c_bm = 0.45, c_bm_pb = 0)
  data.table(N = N, c.bm = c_bm, c.bm.pb = c_bm_pb, carryover_t1half = 1,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.1, c.cfct = 0.05)

SIGMA_BR <- make_resp_params()[cat == 'br', sd]
BM_MEAN  <- make_bl_params()[cat == 'bm', m]
BM_SD    <- make_bl_params()[cat == 'bm', sd]
## Per biomarker SD (c_bm = 0.45 in every cell), because the analysis
## biomarker is standardized.
TRUE_BETA <- -0.45 * SIGMA_BR

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
  ## The generator is done with the raw biomarker; both analyses take
  ## bm from dat through .prep() (docs/74, D1).
  dat[, bm := (bm - BM_MEAN) / BM_SD]
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
n_cores <- n_workers()
cat(sprintf('Study B balanced: %d cells x %d reps x %d analyses | cores=%d | params=%s | true beta=%.3f\nOutput: %s\n',
            nrow(cells), N_REPS, length(ANALYSES), n_cores, param_set(),
            TRUE_BETA, CHECKPOINT_DIR))

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
saveRDS(list(params = param_meta(), true_beta = TRUE_BETA,
             bm_scale = c(mean = BM_MEAN, sd = BM_SD),
             seed = STUDY_SEED, n_reps = N_REPS, workers = n_cores,
             date = Sys.time()),
        file.path(CHECKPOINT_DIR, 'meta.rds'))

cat(sprintf('\nDone in %.1f min. True beta = %.3f\n',
            (proc.time()[['elapsed']] - t0) / 60, TRUE_BETA))
cat('\n=== bias by analysis x c_bm_pb x N ===\n')
print(summary_dt[order(N, c_bm_pb, analysis),
  .(analysis, c_bm_pb, N, bias = round(bias, 3),
    mcse = round(mcse_bias, 3), cov = round(coverage, 3),
    conv = round(conv_rate, 3))])
