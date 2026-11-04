## Study A production runner -- paper 06 component decomposition.
##
## USAGE
##   cd <repo-root>
##   Rscript analysis/scripts/component-decomposition/run-study-a-prod.R \
##     [n_reps_alt] [n_reps_null] [out_dir] [--params NAME]
##   defaults: 1000 5000, out_dir study-a-alt<alt>-null<null>, and
##   --params fixture (or PMSIM_PARAMS when set). The parameter sets
##   (published, simcore2026, fixture) are defined in
##   analysis/scripts/carryover-sensitivity/simulation-core.R; the
##   stored run used the raw fixture, which is why that is the default.
##   SIM_CORES (env) sets the worker count; otherwise detectCores() - 1.
##
## OUTPUT
##   Results land in analysis/data/derived_data/component-decomposition/
##   <out_dir>/ as one .rds per cell plus a collated summary and
##   meta.rds (parameter set and true value). Under any set other than
##   the fixture the default out_dir gains a '-<set>' suffix (e.g.
##   study-a-alt1000-null5000-published), and the folder carries a
##   PARAM_SET marker that refuses a resume under another set. out_dir
##   is taken as given when it is an absolute path.
##   Re-running skips cells whose checkpoint file already exists, so
##   the script is safe to interrupt and resume.
##
## BIOMARKER SCALE (docs/74, D1)
##   The biomarker enters every analysis standardized by the active
##   set's mean and SD, so the true bm:Dbc slope is -c_bm * sigma_BR per
##   biomarker SD (-2.25 under the fixture, -3.6 under published).
##
## PREREQUISITES (host R, not Docker):
##   install.packages(c('pkgload', 'data.table', 'nlme', 'parallel'))
##   # parallel is a base package; others install from CRAN.
##
## NOTE: uses pkgload::load_all() rather than devtools::load_all() so
##   it runs on hosts where devtools is absent or uninstalled.

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

args <- commandArgs(trailingOnly = TRUE)
i_params <- which(args == '--params')
if (length(i_params)) args <- args[-c(i_params[1], i_params[1] + 1L)]

## ---- parameters to change ---------------------------------------- ##
## replicates for alternative cells (c_bm = 0.45) and null cells (c_bm = 0)
N_REPS_ALT  <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
N_REPS_NULL <- if (length(args) >= 2) as.integer(args[[2]]) else 5000L
## ------------------------------------------------------------------ ##

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

STUDY_SEED  <- 642L
CD_ROOT <- file.path('analysis', 'data', 'derived_data',
                     'component-decomposition')
LEGACY_NAME <- sprintf('study-a-alt%d-null%d', N_REPS_ALT, N_REPS_NULL)
OUT_NAME <- if (length(args) >= 3) args[[3]] else set_dir_name(LEGACY_NAME)
CHECKPOINT_DIR <- if (startsWith(OUT_NAME, '/')) OUT_NAME else
  file.path(CD_ROOT, OUT_NAME)
guard_paper06_dir(CHECKPOINT_DIR, file.path(CD_ROOT, LEGACY_NAME))

## ------------------------------------------------------------------ ##
## Trial designs
## ------------------------------------------------------------------ ##

design_hybrid_simple <- function() {
  buildtrialdesign(
    name_longform  = 'OL + BDC (simple)',
    name_shortform = 'OLBDC',
    timepoints   = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20),
    timeptnames  = c(paste0('OL', 1:5), paste0('BDC', 1:3)),
    expectancies = c(rep(1, 5), rep(0.5, 3)),
    ondrug = list(pathA = c(rep(1, 5), 1, 0, 0))
  )
}

design_hybrid_full <- function() {
  buildtrialdesign(
    name_longform  = 'OL + BDC + crossover (full)',
    name_shortform = 'OLBDCxo',
    timepoints   = seq(1, 20, by = 1.25),
    timeptnames  = c(paste0('OL', 1:8), paste0('BDC', 1:4),
                     paste0('XO1_', 1:2), paste0('XO2_', 1:2)),
    expectancies = c(rep(1, 8), rep(0.5, 4), rep(0.5, 4)),
    ondrug = list(
      AB = c(rep(1, 8), rep(0, 4), rep(1, 2), rep(0, 2)),
      BA = c(rep(1, 8), rep(0, 4), rep(0, 2), rep(1, 2))
    )
  )
}

## ------------------------------------------------------------------ ##
## Sample-size allocation
##
## N is the TOTAL across randomization paths (NOTATION.md rule 4) and
## is allocated across them, remainder one per path. The simple design
## has one path and the full design two, so the cell grid's 200 and
## 300 put 100 and 150 on each path of the full design, which is what
## earlier versions of this driver passed to every path while calling
## it N. The draws are unchanged; only the label is corrected.
## ------------------------------------------------------------------ ##

allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) +
    c(rep(1L, n_total %% n_paths),
      rep(0L, n_paths - n_total %% n_paths))
}

## ------------------------------------------------------------------ ##
## Parameter constructors
## ------------------------------------------------------------------ ##

## The base table is the active set's (docs/74, D6); m_PB and m_TV are
## the design axes and override only the PB and TV maxima per cell.
make_resp_params <- function(m_PB = 6.5, m_TV = 0) {
  rp <- as.data.table(default_resp_param())
  rp[cat == 'pb', max := m_PB]
  rp[cat == 'tv', max := m_TV]
  rp
}

make_bl_params <- function() as.data.table(default_baseline_param())

make_model_params <- function(N = 70, c_bm = 0.45, t1half = 1) {
  data.table(N = N, c.bm = c_bm, carryover_t1half = t1half,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.1, c.cfct = 0.05)
}

SIGMA_BR <- make_resp_params()[cat == 'br', sd]
BM_MEAN  <- make_bl_params()[cat == 'bm', m]
BM_SD    <- make_bl_params()[cat == 'bm', sd]

## Per biomarker SD, because the analysis biomarker is standardized.
TRUE_BETA <- function(c_bm = 0.45) -c_bm * SIGMA_BR

## ------------------------------------------------------------------ ##
## Analysis wrappers
## ------------------------------------------------------------------ ##

one_component_fit <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE,
             full_model_out = FALSE, simplecarryover = FALSE,
             carryover_t1half = 1, carryover_scalefactor = 1)
  res <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                  error = function(e) NULL)
  if (is.null(res))
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'fit_failed'))
  data.table(beta = res$beta, betaSE = res$betaSE, p = res$p,
             converged = !is.na(res$beta), formula_dropped = 'none')
}

phase_augmented_fit <- function(td, dat, N) {
  op_prep <- list(useDE = FALSE, t_random_slope = FALSE,
                  full_model_out = TRUE, simplecarryover = FALSE,
                  carryover_t1half = 1, carryover_scalefactor = 1)
  prep <- tryCatch(lme_analysis(td$trialpaths, dat, op_prep),
                   error = function(e) NULL)
  if (is.null(prep) || is.null(prep$datamerged))
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'prep_failed'))
  d <- copy(prep$datamerged)
  d <- d[!is.na(Sx) & !is.na(bm) & !is.na(t) & !is.na(Dbc) & !is.na(De)]
  d <- d[t > 0]
  d[, phase := factor(ifelse(De >= 0.99, 'OL', 'BDC'),
                      levels = c('OL', 'BDC'))]
  if (length(unique(d$phase)) < 2)
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'single_phase'))
  ctl <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)
  formula_dropped <- 'none'
  fit <- NULL
  if (N >= 70) {
    fit <- tryCatch(
      nlme::lme(Sx ~ bm + t + Dbc + phase + Dbc:phase
                     + bm:Dbc + bm:Dbc:phase,
                random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = d, control = ctl),
      error = function(e) NULL)
    if (is.null(fit)) formula_dropped <- 'Dbc:phase'
  } else {
    formula_dropped <- 'Dbc:phase'
  }
  if (is.null(fit)) {
    fit <- tryCatch(
      nlme::lme(Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase,
                random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = d, control = ctl),
      error = function(e) NULL)
  }
  if (is.null(fit))
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = paste0(formula_dropped, '+all')))
  ct     <- summary(fit)$tTable
  target <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(ct))
  if (length(target) == 0)
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = paste0(formula_dropped, '+bmDbc')))
  data.table(beta            = unname(ct[target[1], 'Value']),
             betaSE          = unname(ct[target[1], 'Std.Error']),
             p               = unname(ct[target[1], 'p-value']),
             converged       = TRUE,
             formula_dropped = formula_dropped)
}

## ------------------------------------------------------------------ ##
## One replicate
## ------------------------------------------------------------------ ##

one_rep_06 <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- if (!is.null(cell$design) && cell$design == 'full')
          design_hybrid_full() else design_hybrid_simple()
  rp    <- make_resp_params(m_PB = cell$m_PB, m_TV = cell$m_TV)
  bp    <- make_bl_params()
  paths <- td$trialpaths
  ## cell$N is the TOTAL across paths (NOTATION.md rule 4); allocate
  ## it, remainder one per path.
  n_per_path <- allocate_across_paths(cell$N, length(paths))
  mp    <- make_model_params(N = NA_integer_, c_bm = cell$c_bm,
                             t1half = 1)
  dat_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
    mp$N <- n_per_path[[g]]
    di <- tryCatch(
      generateData(mp, rp, bp, paths[[g]],
                   empirical = FALSE, makePositiveDefinite = TRUE,
                   dgp_architecture = 'mvn'),
      error = function(e) NULL)
    if (is.null(di)) next
    di[, path := g]
    dat_list[[g]] <- di
  }
  good <- !vapply(dat_list, is.null, logical(1))
  if (!any(good)) {
    na_row <- data.table(beta = NA_real_, betaSE = NA_real_,
                         p = NA_real_, converged = FALSE,
                         formula_dropped = 'gen_failed')
    return(rbindlist(lapply(cell$analyses, function(an)
      data.table(analysis = an, rep_idx = rep_idx, na_row))))
  }
  dat <- rbindlist(dat_list[good], fill = TRUE)
  ## The generator is done with the raw biomarker; both analyses take
  ## bm from dat (docs/74, D1).
  dat[, bm := (bm - BM_MEAN) / BM_SD]
  out <- vector('list', 0L)
  for (an in cell$analyses) {
    fit <- switch(an,
      one_component   = one_component_fit(td, dat),
      ## The phase-augmented formula ladder is keyed on the per-path
      ## count, which is what this argument carried before N became
      ## the total. Passing the total here would change which
      ## formula is fitted.
      phase_augmented = phase_augmented_fit(td, dat, n_per_path[[1]]),
      data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                 converged = FALSE, formula_dropped = 'unknown'))
    out[[length(out) + 1L]] <- cbind(analysis = an,
                                     rep_idx = rep_idx, fit)
  }
  rbindlist(out, fill = TRUE)
}

## ------------------------------------------------------------------ ##
## Cell runner
## ------------------------------------------------------------------ ##

run_cell_06 <- function(cell, n_reps, study_seed, cell_id,
                         n_cores = n_workers()) {
  t0 <- proc.time()[['elapsed']]
  el <- function() proc.time()[['elapsed']] - t0
  lbl <- sprintf('cell %02d (m_PB=%2g m_TV=%2g N=%3d)',
                 cell_id, cell$m_PB, cell$m_TV, cell$N)
  cat(sprintf('[%6.1fs] %s start\n', el(), lbl))
  chunks <- mclapply(seq_len(n_reps),
                     function(r) one_rep_06(r, cell, study_seed, cell_id),
                     mc.cores = n_cores)
  good <- !vapply(chunks, function(x)
    inherits(x, 'try-error') || is.null(x), logical(1))
  out <- rbindlist(chunks[good], fill = TRUE)
  out[, `:=`(m_PB = cell$m_PB, m_TV = cell$m_TV,
             N = cell$N, c_bm = cell$c_bm, cell_id = cell_id)]
  cat(sprintf('[%6.1fs] %s done (%d/%d ok)\n',
              el(), lbl, sum(good), n_reps))
  out
}

## ------------------------------------------------------------------ ##
## Morris ADEMP cell summary
## ------------------------------------------------------------------ ##

summarise_cell_06 <- function(d, true_beta) {
  conv <- d[converged == TRUE]
  n    <- nrow(conv)
  if (n == 0)
    return(data.table(n_reps = nrow(d), n_converged = 0L,
                      conv_rate = 0, power = NA_real_,
                      mean_beta = NA_real_, bias = NA_real_,
                      mcse_power = NA_real_, mcse_bias = NA_real_,
                      coverage = NA_real_))
  pw  <- mean(conv$p < 0.05, na.rm = TRUE)
  mb  <- mean(conv$beta, na.rm = TRUE)
  sb  <- sd(conv$beta, na.rm = TRUE)
  lo  <- conv$beta - 1.96 * conv$betaSE
  hi  <- conv$beta + 1.96 * conv$betaSE
  cov <- mean(true_beta >= lo & true_beta <= hi, na.rm = TRUE)
  data.table(n_reps = nrow(d), n_converged = n,
             conv_rate = mean(d$converged, na.rm = TRUE),
             power = pw, mcse_power = sqrt(pw * (1 - pw) / n),
             mean_beta = mb, bias = mb - true_beta,
             mcse_bias = sb / sqrt(n),
             coverage = cov)
}

## ------------------------------------------------------------------ ##
## Cell grid
## ------------------------------------------------------------------ ##

alt_cells  <- CJ(m_PB = c(0, 1, 3, 6, 10),
                 m_TV = c(-1, 0, 1, 2),
                 N    = c(35L, 70L, 200L, 300L),
                 sorted = FALSE)
alt_cells[, `:=`(c_bm = 0.45, regime = 'alt')]

null_cells <- data.table(m_PB = 0, m_TV = 0,
                          N = c(35L, 70L, 200L, 300L),
                          c_bm = 0, regime = 'null')

cells <- rbindlist(list(alt_cells, null_cells),
                   use.names = TRUE, fill = TRUE)
cells[, cell_id := .I]
cells[, analyses := list(list(c('one_component', 'phase_augmented')))]

n_cores <- n_workers()
cat(sprintf(
  'Study A: %d cells | alt=%d reps, null=%d reps | seed=%d | cores=%d | params=%s | true beta=%.3f\nOutput: %s\n',
  nrow(cells), N_REPS_ALT, N_REPS_NULL, STUDY_SEED, n_cores,
  param_set(), TRUE_BETA(), CHECKPOINT_DIR))

## ------------------------------------------------------------------ ##
## Main loop with per-cell checkpointing
## ------------------------------------------------------------------ ##

t_wall <- proc.time()[['elapsed']]

for (i in seq_len(nrow(cells))) {
  cid   <- cells$cell_id[i]
  ckpt  <- file.path(CHECKPOINT_DIR,
                     sprintf('cell-%03d.rds', cid))
  if (file.exists(ckpt)) {
    cat(sprintf('[skip] cell %02d -- checkpoint exists\n', cid))
    next
  }
  cell <- list(m_PB     = cells$m_PB[i],
               m_TV     = cells$m_TV[i],
               N        = cells$N[i],
               c_bm     = cells$c_bm[i],
               analyses = cells$analyses[[i]],
               ## N is now the total across paths, so the threshold
               ## that selects the two-path design moves with it:
               ## the old per-path 100 and 150 cells are the totals
               ## 200 and 300.
               design   = if (cells$N[i] >= 200) 'full' else 'simple')
  n_reps          <- if (cells$regime[i] == 'null') N_REPS_NULL
                     else N_REPS_ALT
  reps_dt         <- run_cell_06(cell, n_reps = n_reps,
                                  study_seed = STUDY_SEED,
                                  cell_id    = cid,
                                  n_cores    = n_cores)
  reps_dt[, regime := cells$regime[i]]
  saveRDS(reps_dt, ckpt)
}

## ------------------------------------------------------------------ ##
## Collate and summarise
## ------------------------------------------------------------------ ##

all_reps <- rbindlist(
  lapply(seq_len(nrow(cells)), function(i) {
    ckpt <- file.path(CHECKPOINT_DIR,
                      sprintf('cell-%03d.rds', cells$cell_id[i]))
    if (!file.exists(ckpt)) return(NULL)
    dt <- readRDS(ckpt)
    if (!'regime' %in% names(dt))
      dt[, regime := cells$regime[i]]
    dt
  }),
  fill = TRUE)

saveRDS(all_reps,
        file.path(CHECKPOINT_DIR, 'study-a-all-reps.rds'))

summary_dt <- all_reps[, {
  tb <- if (unique(c_bm) == 0) 0 else TRUE_BETA()
  summarise_cell_06(.SD, tb)
}, by = .(analysis, regime, cell_id, m_PB, m_TV, N, c_bm)]

saveRDS(summary_dt,
        file.path(CHECKPOINT_DIR, 'study-a-summary.rds'))
write.table(summary_dt,
            file.path(CHECKPOINT_DIR, 'study-a-summary.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE)
saveRDS(list(params = param_meta(), true_beta = TRUE_BETA(),
             bm_scale = c(mean = BM_MEAN, sd = BM_SD),
             seed = STUDY_SEED, n_reps_alt = N_REPS_ALT,
             n_reps_null = N_REPS_NULL, workers = n_cores,
             date = Sys.time()),
        file.path(CHECKPOINT_DIR, 'meta.rds'))

elapsed <- proc.time()[['elapsed']] - t_wall
cat(sprintf(
  '\nDone. %d cells (alt=%d, null=%d reps) in %.1f min.\nSummary: %s\n',
  nrow(cells), N_REPS_ALT, N_REPS_NULL, elapsed / 60,
  file.path(CHECKPOINT_DIR, 'study-a-summary.rds')))
