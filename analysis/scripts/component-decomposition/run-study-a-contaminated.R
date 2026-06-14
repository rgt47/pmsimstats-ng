## Study A (contaminated biomarker) runner -- paper 06.
##
## Decisive cell for the biomarker-validation failure. Study A as
## pre-registered coupled the biomarker to BR only, so varying the PB and
## TV *magnitudes* left the lumped biomarker-by-treatment estimate
## unbiased (the section 6.1 identity: bias = w_PB*beta_bm^PB +
## w_TV*beta_bm^TV = 0 when beta_bm^PB = beta_bm^TV = 0). This run varies
## the biomarker-PB *coupling* c_bm_pb at fixed m_PB, which is the axis the
## identity predicts to drive the bias.
##
## Mechanism. On-drug observations are predominantly open-label (high
## belief); off-drug observations are blinded (lower belief). A biomarker
## coupled to PB therefore loads onto the bm:Dbc contrast in the
## one-component analysis (which pools open-label-on-drug against
## blinded-off-drug), but is separable by the bm:Dbc:phase term, because
## within the blinded phase the on/off-drug contrast holds belief fixed.
## The clean BR estimate from the phase-augmented model is therefore the
## within-blinded contrast beta_bdc = (bm:Dbc) + (bm:Dbc:phaseBDC).
##
## PREDICTIONS (pre-registered here):
##   1. alt cells (true beta_bm:BR = -2.25): one-component bias grows
##      (more negative) with c_bm_pb; phase-augmented beta_bdc stays near
##      -2.25.
##   2. contam_null cells (true beta = 0, BR uncoupled, PB contaminated):
##      one-component manufactures a spurious non-zero bm:Dbc (rejection
##      rate above nominal 0.05, rising with c_bm_pb); phase-augmented
##      beta_bdc stays near 0 with nominal rejection.
##
## USAGE
##   cd <repo-root>
##   Rscript analysis/scripts/component-decomposition/run-study-a-contaminated.R
##   # change only N_REPS below. Resume-safe: existing cell-NNN.rds skipped.

## ---- single parameter to change ---------------------------------- ##
N_REPS <- 250L
## ------------------------------------------------------------------ ##

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  library(nlme)
  library(parallel)
})

STUDY_SEED <- 991L
CHECKPOINT_DIR <- file.path(
  'analysis', 'data', 'derived_data', 'component-decomposition',
  paste0('study-a-contaminated-N', N_REPS))
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------ ##
## Trial designs (identical to run-study-a-prod.R)
## ------------------------------------------------------------------ ##

design_hybrid_simple <- function() {
  buildtrialdesign(
    name_longform  = 'OL + BDC (simple)',
    name_shortform = 'OLBDC',
    timepoints   = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20),
    timeptnames  = c(paste0('OL', 1:5), paste0('BDC', 1:3)),
    expectancies = c(rep(1, 5), rep(0.5, 3)),
    ondrug = list(pathA = c(rep(1, 5), 1, 0, 0)))
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
      BA = c(rep(1, 8), rep(0, 4), rep(0, 2), rep(1, 2))))
}

## ------------------------------------------------------------------ ##
## Parameter constructors
## ------------------------------------------------------------------ ##

make_resp_params <- function(m_PB = 6, m_TV = 0, m_BR = 10.98604) {
  data.table(cat  = c('tv', 'pb', 'br'),
             max  = c(m_TV, m_PB, m_BR),
             disp = c(5, 5, 5),
             rate = c(0.42, 0.35, 0.42),
             sd   = c(5, 2, 5))
}

make_bl_params <- function() {
  data.table(cat = c('bm', 'BL'), m = c(0, 70), sd = c(1, 10))
}

## c_bm_pb is the new biomarker-PB coupling; 0 reproduces Study A.
make_model_params <- function(N = 70, c_bm = 0.45, c_bm_pb = 0,
                              t1half = 1) {
  data.table(N = N, c.bm = c_bm, c.bm.pb = c_bm_pb,
             carryover_t1half = t1half,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.1, c.cfct = 0.05)
}

## True biomarker-by-BR slope: the pharmacological estimand. Unchanged by
## PB contamination, which is the bias we are measuring against it.
TRUE_BETA <- function(c_bm = 0.45, sigma_br = 5, sigma_bm = 1) {
  -c_bm * sigma_br / sigma_bm
}

## ------------------------------------------------------------------ ##
## Analysis wrappers
## ------------------------------------------------------------------ ##

## One-component: pooled bm:Dbc. This is the (predicted) contaminated
## estimate.
one_component_fit <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE,
             full_model_out = FALSE, simplecarryover = FALSE,
             carryover_t1half = 1, carryover_scalefactor = 1)
  res <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                  error = function(e) NULL)
  if (is.null(res))
    return(data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                      beta_contam = NA_real_, p_contam = NA_real_,
                      converged = FALSE, formula_dropped = 'fit_failed'))
  data.table(beta = res$beta, betaSE = res$betaSE, p = res$p,
             beta_contam = NA_real_, p_contam = NA_real_,
             converged = !is.na(res$beta), formula_dropped = 'none')
}

## Phase-augmented: report the clean within-blinded contrast
## beta_bdc = (bm:Dbc) + (bm:Dbc:phaseBDC) as the BR estimate, with its
## SE from the fixed-effects covariance, plus the bm:Dbc:phase term as the
## contamination diagnostic.
phase_augmented_fit <- function(td, dat, N) {
  op_prep <- list(useDE = FALSE, t_random_slope = FALSE,
                  full_model_out = TRUE, simplecarryover = FALSE,
                  carryover_t1half = 1, carryover_scalefactor = 1)
  prep <- tryCatch(lme_analysis(td$trialpaths, dat, op_prep),
                   error = function(e) NULL)
  fail <- function(tag)
    data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
               beta_contam = NA_real_, p_contam = NA_real_,
               converged = FALSE, formula_dropped = tag)
  if (is.null(prep) || is.null(prep$datamerged)) return(fail('prep_failed'))
  d <- copy(prep$datamerged)
  d <- d[!is.na(Sx) & !is.na(bm) & !is.na(t) & !is.na(Dbc) & !is.na(De)]
  d <- d[t > 0]
  d[, phase := factor(ifelse(De >= 0.99, 'OL', 'BDC'),
                      levels = c('OL', 'BDC'))]
  if (length(unique(d$phase)) < 2) return(fail('single_phase'))
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
  if (is.null(fit)) return(fail(paste0(formula_dropped, '+all')))

  fe <- nlme::fixef(fit)
  V  <- as.matrix(vcov(fit))
  nm <- names(fe)
  i_pair   <- which(grepl('bm', nm) & grepl('Dbc', nm) & !grepl('phase', nm))
  i_triple <- which(grepl('bm', nm) & grepl('Dbc', nm) & grepl('phase', nm))
  if (length(i_pair) != 1) return(fail(paste0(formula_dropped, '+no_bmDbc')))

  ## Clean BR estimate = within-blinded contrast.
  L <- numeric(length(fe))
  L[i_pair] <- 1
  if (length(i_triple) >= 1) L[i_triple] <- 1
  beta_bdc <- sum(fe[c(i_pair, i_triple)])
  se_bdc   <- sqrt(as.numeric(t(L) %*% V %*% L))
  p_bdc    <- 2 * pnorm(-abs(beta_bdc / se_bdc))

  ## Contamination diagnostic = the bm:Dbc:phase term.
  if (length(i_triple) >= 1) {
    bc  <- unname(fe[i_triple[1]])
    sec <- sqrt(V[i_triple[1], i_triple[1]])
    pc  <- 2 * pnorm(-abs(bc / sec))
  } else {
    bc <- NA_real_; pc <- NA_real_
  }
  data.table(beta = beta_bdc, betaSE = se_bdc, p = p_bdc,
             beta_contam = bc, p_contam = pc,
             converged = TRUE, formula_dropped = formula_dropped)
}

## ------------------------------------------------------------------ ##
## One replicate
## ------------------------------------------------------------------ ##

one_rep_06 <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- if (!is.null(cell$design) && cell$design == 'full')
          design_hybrid_full() else design_hybrid_simple()
  mp <- make_model_params(N = cell$N, c_bm = cell$c_bm,
                          c_bm_pb = cell$c_bm_pb, t1half = 1)
  rp <- make_resp_params(m_PB = cell$m_PB, m_TV = cell$m_TV)
  bp <- make_bl_params()
  paths <- td$trialpaths
  dat_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
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
    na_row <- data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                         beta_contam = NA_real_, p_contam = NA_real_,
                         converged = FALSE, formula_dropped = 'gen_failed')
    return(rbindlist(lapply(cell$analyses, function(an)
      data.table(analysis = an, rep_idx = rep_idx, na_row))))
  }
  dat <- rbindlist(dat_list[good], fill = TRUE)
  out <- vector('list', 0L)
  for (an in cell$analyses) {
    fit <- switch(an,
      one_component   = one_component_fit(td, dat),
      phase_augmented = phase_augmented_fit(td, dat, cell$N),
      data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                 beta_contam = NA_real_, p_contam = NA_real_,
                 converged = FALSE, formula_dropped = 'unknown'))
    out[[length(out) + 1L]] <- cbind(analysis = an, rep_idx = rep_idx, fit)
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
  lbl <- sprintf('cell %02d (c_bm=%g c_bm_pb=%g N=%3d)',
                 cell_id, cell$c_bm, cell$c_bm_pb, cell$N)
  cat(sprintf('[%6.1fs] %s start\n', el(), lbl))
  chunks <- mclapply(seq_len(n_reps),
                     function(r) one_rep_06(r, cell, study_seed, cell_id),
                     mc.cores = n_cores)
  good <- !vapply(chunks, function(x)
    inherits(x, 'try-error') || is.null(x), logical(1))
  out <- rbindlist(chunks[good], fill = TRUE)
  out[, `:=`(m_PB = cell$m_PB, m_TV = cell$m_TV, N = cell$N,
             c_bm = cell$c_bm, c_bm_pb = cell$c_bm_pb, cell_id = cell_id)]
  cat(sprintf('[%6.1fs] %s done (%d/%d ok)\n',
              el(), lbl, sum(good), n_reps))
  out
}

## ------------------------------------------------------------------ ##
## Morris ADEMP cell summary
## ------------------------------------------------------------------ ##

summarise_cell_06 <- function(d, true_beta) {
  conv <- d[converged == TRUE]
  n <- nrow(conv)
  if (n == 0)
    return(data.table(n_reps = nrow(d), n_converged = 0L, conv_rate = 0,
                      power = NA_real_, mcse_power = NA_real_,
                      mean_beta = NA_real_, bias = NA_real_,
                      mcse_bias = NA_real_, coverage = NA_real_,
                      contam_detect = NA_real_))
  pw  <- mean(conv$p < 0.05, na.rm = TRUE)
  mb  <- mean(conv$beta, na.rm = TRUE)
  sb  <- sd(conv$beta, na.rm = TRUE)
  lo  <- conv$beta - 1.96 * conv$betaSE
  hi  <- conv$beta + 1.96 * conv$betaSE
  cov <- mean(true_beta >= lo & true_beta <= hi, na.rm = TRUE)
  cd  <- if (all(is.na(conv$p_contam))) NA_real_
         else mean(conv$p_contam < 0.05, na.rm = TRUE)
  data.table(n_reps = nrow(d), n_converged = n,
             conv_rate = mean(d$converged, na.rm = TRUE),
             power = pw, mcse_power = sqrt(pw * (1 - pw) / n),
             mean_beta = mb, bias = mb - true_beta,
             mcse_bias = sb / sqrt(n), coverage = cov,
             contam_detect = cd)
}

## ------------------------------------------------------------------ ##
## Cell grid
## ------------------------------------------------------------------ ##

cbmpb <- c(0, 0.15, 0.30, 0.45, 0.60)
Ns    <- c(35L, 70L, 100L, 150L)

## alt: true BR coupling present; vary PB contamination.
alt_cells <- CJ(c_bm_pb = cbmpb, N = Ns, sorted = FALSE)
alt_cells[, `:=`(m_PB = 6, m_TV = 0, c_bm = 0.45, regime = 'alt')]

## contam_null: NO true BR coupling, PB contamination present. Falsifier:
## does contamination alone manufacture a spurious pharmacological signal?
contam_null <- CJ(c_bm_pb = c(0.30, 0.60), N = Ns, sorted = FALSE)
contam_null[, `:=`(m_PB = 6, m_TV = 0, c_bm = 0, regime = 'contam_null')]

## null: no coupling of any kind. Type I baseline.
null_cells <- data.table(c_bm_pb = 0, m_PB = 6, m_TV = 0, N = Ns,
                          c_bm = 0, regime = 'null')

cells <- rbindlist(list(alt_cells, contam_null, null_cells),
                   use.names = TRUE, fill = TRUE)
cells[, cell_id := .I]
cells[, analyses := list(list(c('one_component', 'phase_augmented')))]

n_cores <- max(1L, detectCores() - 1L)
cat(sprintf(
  'Study A (contaminated): %d cells x %d reps | seed=%d | cores=%d\nOutput: %s\n',
  nrow(cells), N_REPS, STUDY_SEED, n_cores, CHECKPOINT_DIR))

## ------------------------------------------------------------------ ##
## Main loop with per-cell checkpointing
## ------------------------------------------------------------------ ##

t_wall <- proc.time()[['elapsed']]

for (i in seq_len(nrow(cells))) {
  cid  <- cells$cell_id[i]
  ckpt <- file.path(CHECKPOINT_DIR, sprintf('cell-%03d.rds', cid))
  if (file.exists(ckpt)) {
    cat(sprintf('[skip] cell %02d -- checkpoint exists\n', cid))
    next
  }
  cell <- list(m_PB = cells$m_PB[i], m_TV = cells$m_TV[i],
               N = cells$N[i], c_bm = cells$c_bm[i],
               c_bm_pb = cells$c_bm_pb[i],
               analyses = cells$analyses[[i]],
               design = if (cells$N[i] >= 100) 'full' else 'simple')
  reps_dt <- run_cell_06(cell, n_reps = N_REPS, study_seed = STUDY_SEED,
                         cell_id = cid, n_cores = n_cores)
  reps_dt[, regime := cells$regime[i]]
  saveRDS(reps_dt, ckpt)
}

## ------------------------------------------------------------------ ##
## Collate and summarise
## ------------------------------------------------------------------ ##

all_reps <- rbindlist(
  lapply(seq_len(nrow(cells)), function(i) {
    ckpt <- file.path(CHECKPOINT_DIR, sprintf('cell-%03d.rds', cells$cell_id[i]))
    if (!file.exists(ckpt)) return(NULL)
    dt <- readRDS(ckpt)
    if (!'regime' %in% names(dt)) dt[, regime := cells$regime[i]]
    dt
  }), fill = TRUE)

saveRDS(all_reps, file.path(CHECKPOINT_DIR, 'study-a-contaminated-all-reps.rds'))

summary_dt <- all_reps[, {
  tb <- if (unique(c_bm) == 0) 0 else TRUE_BETA()
  summarise_cell_06(.SD, tb)
}, by = .(analysis, regime, cell_id, m_PB, m_TV, N, c_bm, c_bm_pb)]

saveRDS(summary_dt, file.path(CHECKPOINT_DIR, 'study-a-contaminated-summary.rds'))
write.table(summary_dt,
            file.path(CHECKPOINT_DIR, 'study-a-contaminated-summary.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE)

elapsed <- proc.time()[['elapsed']] - t_wall
cat(sprintf('\nDone. %d cells x %d reps in %.1f min.\nSummary: %s\n',
            nrow(cells), N_REPS, elapsed / 60,
            file.path(CHECKPOINT_DIR, 'study-a-contaminated-summary.rds')))
