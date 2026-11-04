## run-study-2-feasible.R  (paper 06, Study 2 corrected, 2026-10-28)
##
## Contaminated-biomarker study on feasible cells only. Replaces the
## June 2026 run (run-study-contam.R, study-contam-alt100-null100/) and
## the hybrid recovery run (run-study-b-recovery.R, study-b-recovery-
## pilot/), both of which drew most of their contaminated cells from
## covariance matrices that were not positive definite and were
## repaired without warning (referee report 2026-10-28, M1; the same
## failure mode as docs/65 for paper 02).
##
## Every cell here lies inside the feasibility bound computed by
## 05-feasibility-bounds.R (Schur complement on the biomarker row,
## c_bm = 0.45): c_bm,PB <= 0.30 on the eight-visit and balanced
## designs (bounds 0.321 and 0.315) and c_bm,PB <= 0.10 on the
## sixteen-visit design (bound 0.126). The driver checks the smallest
## eigenvalue of every cell's matrix before any replicate runs and
## stops if one is not positive. Data are generated under the package
## default stopOnNonPD = TRUE; the repair is never enabled.
##
## Designs
##   simple    eight-visit single-path OL+BDC (5 OL at eta = 1, then
##             BDC on/off/off at eta = 0.5), 2.5-week spacing
##   full      sixteen-visit two-path OL+BDC+crossover, 1.25-week spacing
##   balanced  3 OL (on, eta = 1) / 6 covert on-drug (eta = 0) /
##             6 concealed placebo (off, eta = 1), 1.5-week spacing
##
## Parameters: the base response and baseline tables are the active
## set's (analysis/scripts/carryover-sensitivity/simulation-core.R;
## docs/74, D6), with m_PB = 6 and m_TV = 0 overriding the PB and TV
## maxima as before. The default set is the raw tinytest fixture, the
## set the stored run used; the feasibility bounds above depend only on
## correlations and hold under every set.
##
## Analyses (all random intercept + corCAR1(~t|ptID), REML; bm is the
## biomarker standardized by the active set's mean and SD, docs/74 D1)
##   one_component    Sx ~ bm + t + Dbc + bm:Dbc
##   phase_augmented  Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase
##                    (the formula the June runs fitted; the Dbc:phase
##                    term is aliased on every design whose open-label
##                    visits are all on drug, so it is not attempted)
##   f1               Sx ~ bm + t + Dbc + De + bm:Dbc + bm:De
##                    (the series' F1 model, docs/56; the belief-
##                    covariate decomposition of the June Study B)
##
## Estimand: the bm:Dbc coefficient against the pharmacological slope
## per biomarker SD, -c_bm sigma_BR (-2.25 under the fixture, -3.6
## under published) for the alternative, 0 for the null.
##
## Seeds: STUDY_SEED = 2061, set.seed(STUDY_SEED + 1000L * cell_id +
## rep_idx), independent of the June seeds 642/644/916/815.
##
## USAGE (from the repository root, host R):
##   Rscript analysis/scripts/component-decomposition/run-study-2-feasible.R \
##     [n_reps_alt] [n_reps_null] [out_dir] [--params NAME]
##   defaults: 1000 2000 study-2-feasible, --params fixture (or
##   PMSIM_PARAMS when set). Under any set other than the fixture the
##   default out_dir is study-2-feasible-<set> (e.g.
##   study-2-feasible-published).
##   out_dir is a subdirectory of
##   analysis/data/derived_data/component-decomposition/ unless it is
##   an absolute path (used for the smoke test). It carries a PARAM_SET
##   marker and refuses a resume under another set.
##   SIM_CORES (env) sets the worker count; otherwise detectCores() - 1.
## Smoke test (seconds):
##   Rscript analysis/scripts/component-decomposition/run-study-2-feasible.R \
##     3 3 /tmp/study-2-smoke --params published
## Per-cell checkpoints hold every replicate; rerunning skips cells
## whose checkpoint exists. meta.rds records the parameter set and the
## true value.

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  library(nlme)
  library(parallel)
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
## reports the host, not the job (job 49 ran 11 workers on 8 cores).
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
N_REPS_ALT  <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
N_REPS_NULL <- if (length(args) >= 2) as.integer(args[[2]]) else 2000L
OUT_NAME    <- if (length(args) >= 3) args[[3]] else
  set_dir_name('study-2-feasible')
STUDY_SEED  <- 2061L
CD_ROOT <- file.path('analysis', 'data', 'derived_data',
                     'component-decomposition')
OUT_DIR <- if (startsWith(OUT_NAME, '/')) OUT_NAME else
  file.path(CD_ROOT, OUT_NAME)
guard_paper06_dir(OUT_DIR, file.path(CD_ROOT, 'study-2-feasible'))

## ------------------------------------------------------------------ ##
## Designs and parameters (identical to 05-feasibility-bounds.R)
## ------------------------------------------------------------------ ##

design_simple <- function() buildtrialdesign(
  name_longform = 'OL + BDC (simple)', name_shortform = 'OLBDC',
  timepoints = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20),
  timeptnames = c(paste0('OL', 1:5), paste0('BDC', 1:3)),
  expectancies = c(rep(1, 5), rep(0.5, 3)),
  ondrug = list(pathA = c(rep(1, 5), 1, 0, 0)))

design_full <- function() buildtrialdesign(
  name_longform = 'OL + BDC + crossover (full)', name_shortform = 'OLBDCxo',
  timepoints = seq(1, 20, by = 1.25),
  timeptnames = c(paste0('OL', 1:8), paste0('BDC', 1:4),
                  paste0('XO1_', 1:2), paste0('XO2_', 1:2)),
  expectancies = c(rep(1, 8), rep(0.5, 4), rep(0.5, 4)),
  ondrug = list(AB = c(rep(1, 8), rep(0, 4), rep(1, 2), rep(0, 2)),
                BA = c(rep(1, 8), rep(0, 4), rep(0, 2), rep(1, 2))))

design_balanced <- function() buildtrialdesign(
  name_longform = 'balanced-placebo (OL/covert/concealed placebo)',
  name_shortform = 'BAL',
  timepoints = seq(1.5, 22.5, by = 1.5),
  timeptnames = c(paste0('OL', 1:3), paste0('COV', 1:6), paste0('OLP', 1:6)),
  expectancies = c(rep(1, 3), rep(0, 6), rep(1, 6)),
  ondrug = list(p = c(rep(1, 3), rep(1, 6), rep(0, 6))))

DESIGNS <- list(simple = design_simple, full = design_full,
                balanced = design_balanced)

## Base table from the active set (docs/74, D6); m_PB and m_TV override
## only the PB and TV maxima.
make_resp_params <- function(m_PB = 6, m_TV = 0) {
  rp <- as.data.table(default_resp_param())
  rp[cat == 'pb', max := m_PB]
  rp[cat == 'tv', max := m_TV]
  rp
}
make_bl_params <- function() as.data.table(default_baseline_param())
make_model_params <- function(N, c_bm, c_bm_pb)
  data.table(N = N, c.bm = c_bm, c.bm.pb = c_bm_pb, carryover_t1half = 1,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.1, c.cfct = 0.05)
SIGMA_BR <- make_resp_params()[cat == 'br', sd]
BM_MEAN  <- make_bl_params()[cat == 'bm', m]
BM_SD    <- make_bl_params()[cat == 'bm', sd]
## Per biomarker SD, because the analysis biomarker is standardized.
TRUE_BETA <- function(c_bm) -c_bm * SIGMA_BR

## N is the total across paths (NOTATION.md rule 4), allocated with the
## remainder one per path.
allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) + c(rep(1L, n_total %% n_paths),
                         rep(0L, n_paths - n_total %% n_paths))
}

## ------------------------------------------------------------------ ##
## Cell grid
## ------------------------------------------------------------------ ##

alt_cells <- rbindlist(list(
  CJ(design = 'simple',   N = c(35L, 70L),   c_bm_pb = c(0, 0.10, 0.20, 0.30)),
  CJ(design = 'full',     N = c(140L, 300L), c_bm_pb = c(0, 0.05, 0.10)),
  CJ(design = 'balanced', N = c(70L, 150L),  c_bm_pb = c(0, 0.10, 0.20, 0.30))))
alt_cells[, `:=`(c_bm = 0.45, regime = 'alt')]

## Null cells: no interaction, with and without the prognostic
## biomarker-PB correlation (docs/65 Section 3.4 flags the
## expectancy-scaled PB channel as a possible size inflation).
null_cells <- rbindlist(list(
  data.table(design = 'simple',   N = 70L,  c_bm_pb = c(0, 0.30)),
  data.table(design = 'full',     N = 300L, c_bm_pb = c(0, 0.10)),
  data.table(design = 'balanced', N = 150L, c_bm_pb = c(0, 0.30))))
null_cells[, `:=`(c_bm = 0, regime = 'null')]

cells <- rbindlist(list(alt_cells, null_cells), use.names = TRUE)
cells[, cell_id := .I]
cells[, m_PB := 6]
ANALYSES <- c('one_component', 'phase_augmented', 'f1')

## ------------------------------------------------------------------ ##
## Pre-flight feasibility check: every cell, every path, raw matrix
## ------------------------------------------------------------------ ##

min_eig_cell <- function(design, c_bm, c_bm_pb) {
  td <- DESIGNS[[design]]()
  min(vapply(td$trialpaths, function(p) {
    s <- buildSigma(make_model_params(50, c_bm, c_bm_pb), make_resp_params(),
                    make_bl_params(), p, makePositiveDefinite = FALSE,
                    stopOnNonPD = FALSE)$sigma
    pos <- sqrt(diag(s)) > 0
    min(eigen(s[pos, pos], symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1)))
}
cells[, min_eig := mapply(min_eig_cell, design, c_bm, c_bm_pb)]
if (any(cells$min_eig <= 0)) {
  print(cells[min_eig <= 0])
  stop('Infeasible cell(s) in the grid; see the table above.')
}
cat('Pre-flight: all', nrow(cells), 'cells positive definite (smallest',
    'eigenvalue', signif(min(cells$min_eig), 3), ').\n')

## ------------------------------------------------------------------ ##
## Analyses
## ------------------------------------------------------------------ ##

CTL <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)
.na_row <- function(tag)
  data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
             converged = FALSE, formula_dropped = tag)
.extract <- function(fit, tag = 'none') {
  if (is.null(fit)) return(.na_row('fit_failed'))
  ct <- summary(fit)$tTable
  target <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(ct))
  if (length(target) == 0) return(.na_row('no_bmDbc'))
  data.table(beta = unname(ct[target[1], 'Value']),
             betaSE = unname(ct[target[1], 'Std.Error']),
             p = unname(ct[target[1], 'p-value']),
             converged = TRUE, formula_dropped = tag)
}

## Merged analysis frame (Sx, bm, t, Dbc, De) from the package's own
## lme_analysis preparation step, with D_bc built at t1/2 = 1 week.
.prep <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE, full_model_out = TRUE,
             simplecarryover = FALSE, carryover_t1half = 1,
             carryover_scalefactor = 1)
  prep <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                   error = function(e) NULL)
  if (is.null(prep) || is.null(prep$datamerged)) return(NULL)
  d <- copy(prep$datamerged)
  d <- d[!is.na(Sx) & !is.na(bm) & !is.na(t) & !is.na(Dbc) & !is.na(De)]
  d[t > 0]
}

one_component_fit <- function(d)
  .extract(tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = d, control = CTL), error = function(e) NULL))

phase_augmented_fit <- function(d) {
  d <- copy(d)
  d[, phase := factor(ifelse(De >= 0.99, 'OL', 'other'),
                      levels = c('OL', 'other'))]
  if (length(unique(d$phase)) < 2) return(.na_row('single_phase'))
  .extract(tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase,
              random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = d, control = CTL), error = function(e) NULL),
    tag = 'Dbc:phase_not_attempted')
}

f1_fit <- function(d) {
  if (length(unique(d$De)) < 2 || length(unique(round(d$Dbc, 3))) < 2)
    return(.na_row('no_DbcDe_var'))
  .extract(tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + De + bm:Dbc + bm:De, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = d, control = CTL), error = function(e) NULL))
}

## ------------------------------------------------------------------ ##
## One replicate
## ------------------------------------------------------------------ ##

one_rep <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- DESIGNS[[cell$design]]()
  rp <- make_resp_params(m_PB = cell$m_PB, m_TV = 0)
  bp <- make_bl_params()
  paths <- td$trialpaths
  n_per_path <- allocate_across_paths(cell$N, length(paths))
  dat_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
    mp <- make_model_params(n_per_path[[g]], cell$c_bm, cell$c_bm_pb)
    ## stopOnNonPD is left at its default TRUE: a non-PD matrix stops
    ## the replicate (recorded as gen_refused) and is never repaired.
    di <- tryCatch(
      generateData(mp, rp, bp, paths[[g]], empirical = FALSE,
                   makePositiveDefinite = FALSE, dgp_architecture = 'mvn'),
      error = function(e) NULL)
    if (is.null(di)) next
    di[, path := g]
    dat_list[[g]] <- di
  }
  good <- !vapply(dat_list, is.null, logical(1))
  if (!all(good))
    return(rbindlist(lapply(ANALYSES, function(an)
      cbind(analysis = an, rep_idx = rep_idx, .na_row('gen_refused')))))
  dat <- rbindlist(dat_list, fill = TRUE)
  ## The generator is done with the raw biomarker; every analysis takes
  ## bm from dat through .prep() (docs/74, D1).
  dat[, bm := (bm - BM_MEAN) / BM_SD]
  d <- .prep(td, dat)
  if (is.null(d))
    return(rbindlist(lapply(ANALYSES, function(an)
      cbind(analysis = an, rep_idx = rep_idx, .na_row('prep_failed')))))
  rbindlist(list(
    cbind(analysis = 'one_component',   rep_idx = rep_idx, one_component_fit(d)),
    cbind(analysis = 'phase_augmented', rep_idx = rep_idx, phase_augmented_fit(d)),
    cbind(analysis = 'f1',              rep_idx = rep_idx, f1_fit(d))),
    fill = TRUE)
}

run_cell <- function(cell, n_reps, study_seed, cell_id, n_cores) {
  chunks <- mclapply(seq_len(n_reps),
                     function(r) one_rep(r, cell, study_seed, cell_id),
                     mc.cores = n_cores)
  good <- !vapply(chunks, function(x)
    inherits(x, 'try-error') || is.null(x), logical(1))
  out <- rbindlist(chunks[good], fill = TRUE)
  out[, `:=`(design = cell$design, N = cell$N, c_bm_pb = cell$c_bm_pb,
             c_bm = cell$c_bm, m_PB = cell$m_PB, regime = cell$regime,
             cell_id = cell_id)]
  out
}

## ------------------------------------------------------------------ ##
## Cell summary with MCSE on every performance measure
## ------------------------------------------------------------------ ##

summarise_cell <- function(d, true_beta) {
  conv <- d[converged == TRUE]
  n <- nrow(conv)
  if (n == 0)
    return(data.table(n_reps = nrow(d), n_converged = 0L, conv_rate = 0,
                      mean_beta = NA_real_, sd_beta = NA_real_,
                      bias = NA_real_, mcse_bias = NA_real_,
                      power = NA_real_, mcse_power = NA_real_,
                      coverage = NA_real_, mcse_coverage = NA_real_,
                      n_gen_refused = sum(d$formula_dropped == 'gen_refused')))
  mb <- mean(conv$beta); sb <- sd(conv$beta)
  pw <- mean(conv$p < 0.05)
  lo <- conv$beta - 1.96 * conv$betaSE; hi <- conv$beta + 1.96 * conv$betaSE
  cv <- mean(true_beta >= lo & true_beta <= hi)
  data.table(n_reps = nrow(d), n_converged = n, conv_rate = n / nrow(d),
             mean_beta = mb, sd_beta = sb,
             bias = mb - true_beta, mcse_bias = sb / sqrt(n),
             power = pw, mcse_power = sqrt(pw * (1 - pw) / n),
             coverage = cv, mcse_coverage = sqrt(cv * (1 - cv) / n),
             n_gen_refused = sum(d$formula_dropped == 'gen_refused'))
}

## ------------------------------------------------------------------ ##
## Main loop with per-cell checkpoints and timing
## ------------------------------------------------------------------ ##

n_cores <- n_workers()
cat(sprintf('Study 2 (feasible cells): %d cells | alt %d reps, null %d reps | seed %d | cores %d | params %s | true beta %.3f\nOutput: %s\n',
            nrow(cells), N_REPS_ALT, N_REPS_NULL, STUDY_SEED, n_cores,
            param_set(), TRUE_BETA(0.45), OUT_DIR))
timing <- list()
t_wall <- proc.time()[['elapsed']]
for (i in seq_len(nrow(cells))) {
  cid <- cells$cell_id[i]
  ckpt <- file.path(OUT_DIR, sprintf('cell-%03d.rds', cid))
  if (file.exists(ckpt)) {
    cat(sprintf('  [skip] cell %02d, checkpoint exists\n', cid)); next
  }
  cell <- as.list(cells[i])
  n_reps <- if (cell$regime == 'null') N_REPS_NULL else N_REPS_ALT
  t0 <- proc.time()[['elapsed']]
  cat(sprintf('  cell %02d/%d (%s N=%d c_bm_pb=%.2f %s) ...',
              cid, nrow(cells), cell$design, cell$N, cell$c_bm_pb, cell$regime))
  res <- run_cell(cell, n_reps, STUDY_SEED, cid, n_cores)
  saveRDS(res, ckpt)
  el <- proc.time()[['elapsed']] - t0
  timing[[length(timing) + 1L]] <- data.table(cell_id = cid, n_reps = n_reps,
                                              seconds = el, cores = n_cores)
  cat(sprintf(' %.1f s (%.2f s per replicate per core)\n', el, el * n_cores / n_reps))
}
if (length(timing))
  saveRDS(rbindlist(timing), file.path(OUT_DIR, 'timing.rds'))

## ------------------------------------------------------------------ ##
## Collate and summarise
## ------------------------------------------------------------------ ##

all_reps <- rbindlist(lapply(cells$cell_id, function(cid) {
  ckpt <- file.path(OUT_DIR, sprintf('cell-%03d.rds', cid))
  if (file.exists(ckpt)) readRDS(ckpt) else NULL
}), fill = TRUE)
saveRDS(all_reps, file.path(OUT_DIR, 'study-2-all-reps.rds'))

summary_dt <- all_reps[, summarise_cell(.SD, TRUE_BETA(c_bm[1])),
                       by = .(analysis, regime, cell_id, design, N, c_bm_pb,
                              c_bm, m_PB)]
summary_dt <- merge(summary_dt, cells[, .(cell_id, min_eig)], by = 'cell_id')
saveRDS(summary_dt, file.path(OUT_DIR, 'study-2-summary.rds'))
write.table(summary_dt, file.path(OUT_DIR, 'study-2-summary.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE)
saveRDS(list(params = param_meta(), true_beta_alt = TRUE_BETA(0.45),
             bm_scale = c(mean = BM_MEAN, sd = BM_SD),
             seed = STUDY_SEED, n_reps_alt = N_REPS_ALT,
             n_reps_null = N_REPS_NULL, workers = n_cores,
             date = Sys.time()),
        file.path(OUT_DIR, 'meta.rds'))

cat(sprintf('\nDone: %d cells in %.1f min. Refused draws: %d.\n',
            nrow(cells), (proc.time()[['elapsed']] - t_wall) / 60,
            sum(summary_dt$n_gen_refused)))
options(width = 160)
print(summary_dt[order(regime, design, N, c_bm_pb, analysis),
                 .(regime, design, N, c_bm_pb, analysis, n_converged,
                   bias = round(bias, 3), mcse = round(mcse_bias, 3),
                   power = round(power, 3), coverage = round(coverage, 3))])
