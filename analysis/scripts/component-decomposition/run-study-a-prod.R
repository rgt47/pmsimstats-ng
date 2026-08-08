## Study A production runner -- paper 06 component decomposition.
##
## USAGE
##   cd <repo-root>
##   Rscript analysis/scripts/component-decomposition/run-study-a-prod.R
##
## CONFIGURATION -- change only N_REPS below.
##   Results land in analysis/data/derived_data/component-decomposition/
##   study-a-N<nreps>/  as one .rds per cell plus a collated summary.
##   Re-running skips cells whose checkpoint file already exists, so
##   the script is safe to interrupt and resume.
##
## PREREQUISITES (host R, not Docker):
##   install.packages(c('pkgload', 'data.table', 'nlme', 'parallel'))
##   # parallel is a base package; others install from CRAN.
##
## NOTE: uses pkgload::load_all() rather than devtools::load_all() so
##   it runs on hosts where devtools is absent or uninstalled.

## ---- parameters to change ---------------------------------------- ##
N_REPS_ALT  <- 1000L   # replicates for alternative cells (c_bm = 0.45)
N_REPS_NULL <- 5000L   # replicates for null cells (c_bm = 0)
## ------------------------------------------------------------------ ##

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  library(nlme)
  library(parallel)
})

STUDY_SEED  <- 642L
CHECKPOINT_DIR <- file.path(
  'analysis', 'data', 'derived_data', 'component-decomposition',
  sprintf('study-a-alt%d-null%d', N_REPS_ALT, N_REPS_NULL))
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

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

make_resp_params <- function(m_PB = 6.5, m_TV = 0, m_BR = 10.98604) {
  data.table(cat  = c('tv', 'pb', 'br'),
             max  = c(m_TV, m_PB, m_BR),
             disp = c(5, 5, 5),
             rate = c(0.42, 0.35, 0.42),
             sd   = c(5, 2, 5))
}

make_bl_params <- function() {
  data.table(cat = c('bm', 'BL'), m = c(0, 70), sd = c(1, 10))
}

make_model_params <- function(N = 70, c_bm = 0.45, t1half = 1) {
  data.table(N = N, c.bm = c_bm, carryover_t1half = t1half,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.1, c.cfct = 0.05)
}

TRUE_BETA <- function(c_bm = 0.45, sigma_br = 5, sigma_bm = 1) {
  -c_bm * sigma_br / sigma_bm
}

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
                         n_cores = max(1L, detectCores() - 1L)) {
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

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf(
  'Study A: %d cells | alt=%d reps, null=%d reps | seed=%d | cores=%d\nOutput: %s\n',
  nrow(cells), N_REPS_ALT, N_REPS_NULL, STUDY_SEED, n_cores,
  CHECKPOINT_DIR))

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

elapsed <- proc.time()[['elapsed']] - t_wall
cat(sprintf(
  '\nDone. %d cells (alt=%d, null=%d reps) in %.1f min.\nSummary: %s\n',
  nrow(cells), N_REPS_ALT, N_REPS_NULL, elapsed / 60,
  file.path(CHECKPOINT_DIR, 'study-a-summary.rds')))
