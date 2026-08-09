## Study B (recovery): does the component decomposition recover the
## unbiased biomarker-by-BR slope under a contaminated biomarker?
##
## Paper 06 round-2 referee point R1. The contaminated-biomarker pilot
## showed that BOTH the one-component and the phase-augmented analyses
## are biased when the biomarker is coupled to PB (c.bm.pb > 0). This
## driver adds a third analysis -- a belief-covariate decomposition --
## and tests whether it recovers the pharmacological slope (-2.25).
##
## Generative structure (see R/generateData.R buildSigma):
##   * BM-BR coupling (c.bm) is applied ONLY at on-drug timepoints, so
##     the pharmacological signal is carried by the drug-state Dbc.
##   * BM-PB coupling (c.bm.pb) is applied at EVERY timepoint; the
##     observed PB contribution is scaled by expectancy De (1 in OL,
##     0.5 in BDC).
## Because Dbc (drug state) and De (expectancy) vary independently
## across phases, interacting the biomarker with BOTH separates the
## channels: bm:Dbc isolates beta_bm^BR, bm:De absorbs the PB channel.
##
## USAGE (from repo root):
##   Rscript analysis/scripts/component-decomposition/run-study-b-recovery.R \
##     [n_reps]            # default 100 (pilot)
##
## Writes summary to
##   analysis/data/derived_data/component-decomposition/study-b-recovery-pilot/

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  library(nlme)
  library(parallel)
})

args        <- commandArgs(trailingOnly = TRUE)
N_REPS      <- if (length(args) >= 1) as.integer(args[[1]]) else 100L
STUDY_SEED  <- 815L
CHECKPOINT_DIR <- file.path(
  'analysis', 'data', 'derived_data', 'component-decomposition',
  'study-b-recovery-pilot')
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------ ##
## Trial designs / parameter constructors (shared with run-study-contam)
## ------------------------------------------------------------------ ##

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
      BA = c(rep(1, 8), rep(0, 4), rep(0, 2), rep(1, 2))))
}

## N is the TOTAL across randomization paths (NOTATION.md rule 4) and
## is allocated across them, remainder one per path. Unlike the study-A
## and contamination drivers, this study uses the two-path full design
## for every cell, so there is no design threshold keyed on N: the
## grid's 140 and 300 put 70 and 150 on each path, which are the values
## earlier versions passed to every path while calling them N. The
## draws are unchanged; only the label is corrected.

allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) +
    c(rep(1L, n_total %% n_paths),
      rep(0L, n_paths - n_total %% n_paths))
}

make_resp_params <- function(m_PB = 6, m_TV = 0, m_BR = 10.98604) {
  data.table(cat  = c('tv', 'pb', 'br'),
             max  = c(m_TV, m_PB, m_BR),
             disp = c(5, 5, 5),
             rate = c(0.42, 0.35, 0.42),
             sd   = c(5, 2, 5))
}
make_bl_params <- function()
  data.table(cat = c('bm', 'BL'), m = c(0, 70), sd = c(1, 10))
make_model_params <- function(N = 150, c_bm = 0.45, c_bm_pb = 0,
                              t1half = 1)
  data.table(N = N, c.bm = c_bm, c.bm.pb = c_bm_pb,
             carryover_t1half = t1half,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.1, c.cfct = 0.05)
TRUE_BETA <- function(c_bm = 0.45, sigma_br = 5, sigma_bm = 1)
  -c_bm * sigma_br / sigma_bm

## ------------------------------------------------------------------ ##
## Analyses: one-component, phase-augmented, decomposition
## ------------------------------------------------------------------ ##

.prep <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE,
             full_model_out = TRUE, simplecarryover = FALSE,
             carryover_t1half = 1, carryover_scalefactor = 1)
  prep <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                   error = function(e) NULL)
  if (is.null(prep) || is.null(prep$datamerged)) return(NULL)
  d <- copy(prep$datamerged)
  d <- d[!is.na(Sx) & !is.na(bm) & !is.na(t) & !is.na(Dbc) & !is.na(De)]
  d[t > 0]
}

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

CTL <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)

one_component_fit <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE,
             full_model_out = FALSE, simplecarryover = FALSE,
             carryover_t1half = 1, carryover_scalefactor = 1)
  res <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                  error = function(e) NULL)
  if (is.null(res)) return(.na_row('fit_failed'))
  data.table(beta = res$beta, betaSE = res$betaSE, p = res$p,
             converged = !is.na(res$beta), formula_dropped = 'none')
}

phase_augmented_fit <- function(td, dat, N) {
  d <- .prep(td, dat); if (is.null(d)) return(.na_row('prep_failed'))
  d[, phase := factor(ifelse(De >= 0.99, 'OL', 'BDC'),
                      levels = c('OL', 'BDC'))]
  if (length(unique(d$phase)) < 2) return(.na_row('single_phase'))
  fit <- tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + phase + Dbc:phase
                   + bm:Dbc + bm:Dbc:phase,
              random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = d, control = CTL), error = function(e) NULL)
  if (is.null(fit))
    fit <- tryCatch(
      nlme::lme(Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase,
                random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = d, control = CTL), error = function(e) NULL)
  .extract(fit)
}

## Belief-covariate decomposition: interact bm with BOTH the drug state
## (Dbc) and the expectancy (De). bm:Dbc -> beta_bm^BR; bm:De -> PB
## channel. Identified because (Dbc, De) take distinct combinations
## across OL (1,1), BDC-on (1,0.5), BDC-off (0,0.5).
decomposition_fit <- function(td, dat, N) {
  d <- .prep(td, dat); if (is.null(d)) return(.na_row('prep_failed'))
  if (length(unique(d$De)) < 2 ||
      length(unique(round(d$Dbc, 3))) < 2)
    return(.na_row('no_DbcDe_var'))
  fit <- tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + De + bm:Dbc + bm:De,
              random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = d, control = CTL), error = function(e) NULL)
  .extract(fit)
}

## Blinded-stratum decomposition: the manuscript's own BR-identification
## logic. beta_bm^BR is the biomarker's coupling to the drug-on/off
## contrast WITHIN a constant-belief stratum, where PB is held fixed
## (eta = 0.5) and therefore cannot contaminate the contrast. Restrict
## to the blinded / crossover timepoints (De < 1) and read bm:Dbc.
decomposition_blinded_fit <- function(td, dat, N) {
  d <- .prep(td, dat); if (is.null(d)) return(.na_row('prep_failed'))
  d <- d[De < 0.99]                       # constant-belief stratum
  if (nrow(d) < 10 || length(unique(round(d$Dbc, 3))) < 2 ||
      length(unique(d$ptID)) < 5)
    return(.na_row('no_blinded_Dbc_var'))
  fit <- tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc,
              random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = d, control = CTL), error = function(e) NULL)
  .extract(fit)
}

## ------------------------------------------------------------------ ##
## One replicate / cell runner / summary
## ------------------------------------------------------------------ ##

ANALYSES <- c('one_component', 'phase_augmented', 'decomposition',
              'decomposition_blinded')

one_rep <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- design_hybrid_full()
  rp <- make_resp_params(m_PB = cell$m_PB, m_TV = 0)
  bp <- make_bl_params()
  paths <- td$trialpaths
  ## cell$N is the TOTAL across paths (NOTATION.md rule 4).
  n_per_path <- allocate_across_paths(cell$N, length(paths))
  mp <- make_model_params(N = NA_integer_, c_bm = cell$c_bm,
                          c_bm_pb = cell$c_bm_pb)
  dat_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
    mp$N <- n_per_path[[g]]
    di <- tryCatch(
      generateData(mp, rp, bp, paths[[g]], empirical = FALSE,
                   makePositiveDefinite = TRUE, dgp_architecture = 'mvn'),
      error = function(e) NULL)
    if (is.null(di)) next
    di[, path := g]; dat_list[[g]] <- di
  }
  good <- !vapply(dat_list, is.null, logical(1))
  if (!any(good))
    return(rbindlist(lapply(ANALYSES, function(an)
      cbind(analysis = an, rep_idx = rep_idx, .na_row('gen_failed')))))
  dat <- rbindlist(dat_list[good], fill = TRUE)
  rbindlist(lapply(ANALYSES, function(an) {
    fit <- switch(an,
      one_component         = one_component_fit(td, dat),
      ## These three take the per-path count, which is what this
      ## argument carried before N became the total. None of them
      ## currently reads it, but keeping it per-path preserves the
      ## meaning if a formula ladder keyed on it is reintroduced, as
      ## the study-A and contamination drivers have.
      phase_augmented       = phase_augmented_fit(td, dat,
                                                  n_per_path[[1]]),
      decomposition         = decomposition_fit(td, dat,
                                                n_per_path[[1]]),
      decomposition_blinded = decomposition_blinded_fit(
                                td, dat, n_per_path[[1]]))
    cbind(analysis = an, rep_idx = rep_idx, fit)
  }), fill = TRUE)
}

run_cell <- function(cell, n_reps, study_seed, cell_id, n_cores) {
  chunks <- mclapply(seq_len(n_reps),
                     function(r) one_rep(r, cell, study_seed, cell_id),
                     mc.cores = n_cores)
  good <- !vapply(chunks, function(x)
    inherits(x, 'try-error') || is.null(x), logical(1))
  out <- rbindlist(chunks[good], fill = TRUE)
  out[, `:=`(c_bm_pb = cell$c_bm_pb, m_PB = cell$m_PB, N = cell$N,
             c_bm = cell$c_bm, cell_id = cell_id)]
  out
}

summarise_cell <- function(d, true_beta) {
  conv <- d[converged == TRUE]; n <- nrow(conv)
  if (n == 0)
    return(data.table(n_reps = nrow(d), n_converged = 0L, conv_rate = 0,
                      mean_beta = NA_real_, bias = NA_real_,
                      mcse_bias = NA_real_, power = NA_real_))
  mb <- mean(conv$beta, na.rm = TRUE); sb <- sd(conv$beta, na.rm = TRUE)
  data.table(n_reps = nrow(d), n_converged = n,
             conv_rate = mean(d$converged, na.rm = TRUE),
             mean_beta = mb, bias = mb - true_beta,
             mcse_bias = sb / sqrt(n),
             power = mean(conv$p < 0.05, na.rm = TRUE))
}

## ------------------------------------------------------------------ ##
## Pilot grid: contamination axis at N=300 total (full design)
##
## N is the total across the full design's two paths, so 140 and 300
## are 70 and 150 per path, the values the old grid recorded as N.
## ------------------------------------------------------------------ ##

cells <- CJ(c_bm_pb = c(0, 0.15, 0.30, 0.45), m_PB = 0,
            N = c(140L, 300L), sorted = FALSE)
cells[, `:=`(c_bm = 0.45, cell_id = .I)]

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf('Study B recovery: %d cells x %d reps x %d analyses | cores=%d\nOutput: %s\n',
            nrow(cells), N_REPS, length(ANALYSES), n_cores, CHECKPOINT_DIR))

t0 <- proc.time()[['elapsed']]
for (i in seq_len(nrow(cells))) {
  cid  <- cells$cell_id[i]
  ckpt <- file.path(CHECKPOINT_DIR, sprintf('cell-%03d.rds', cid))
  if (file.exists(ckpt)) {
    cat(sprintf('  [skip] cell %d/%d -- checkpoint exists\n', i, nrow(cells)))
    next
  }
  cell <- as.list(cells[i])
  cat(sprintf('  cell %d/%d (c_bm_pb=%.2f N=%d) ...\n', i, nrow(cells),
              cell$c_bm_pb, cell$N))
  saveRDS(run_cell(cell, N_REPS, STUDY_SEED, cid, n_cores), ckpt)
}

all_reps <- rbindlist(lapply(cells$cell_id, function(cid) {
  ckpt <- file.path(CHECKPOINT_DIR, sprintf('cell-%03d.rds', cid))
  if (file.exists(ckpt)) readRDS(ckpt) else NULL
}), fill = TRUE)

saveRDS(all_reps, file.path(CHECKPOINT_DIR, 'study-b-recovery-all-reps.rds'))

tb <- TRUE_BETA()
summary_dt <- all_reps[, summarise_cell(.SD, tb),
                       by = .(analysis, c_bm_pb, m_PB, N, c_bm)]
saveRDS(summary_dt, file.path(CHECKPOINT_DIR, 'study-b-recovery-summary.rds'))
write.table(summary_dt,
            file.path(CHECKPOINT_DIR, 'study-b-recovery-summary.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE)

cat(sprintf('\nDone in %.1f min. True beta = %.3f\n',
            (proc.time()[['elapsed']] - t0) / 60, tb))
cat('\n=== Recovery comparison: bias by analysis x c_bm_pb x N ===\n')
print(dcast(summary_dt, c_bm_pb + N ~ analysis,
            value.var = 'bias')[order(N, c_bm_pb)])
cat('\n--- full summary ---\n')
print(summary_dt[order(analysis, c_bm_pb),
                 .(analysis, c_bm_pb, mean_beta = round(mean_beta, 3),
                   bias = round(bias, 3), mcse_bias = round(mcse_bias, 3),
                   conv_rate)])
