## Study 5 production driver — small-N identifiability and Type I
## error pre-flight for paper 03 path (b).
##
## This driver is the Phase B kill-switch: it establishes whether
## lcmm::hlme is usable at trial-relevant N at all. Per the ADEMP
## pre-registration (00-ademp-pre-registration.md, §"Study 5"),
## three DGP families are exercised:
##
##   Cell 1: Single-component Gaussian (c.bm = 0).
##           Type I error sub-cell: target rejection 0.05.
##   Cell 2: Single-component, mildly non-Gaussian residuals on BR.
##           Two sub-families: t5-distributed and skew-normal
##           (shape = 4). Class-count selection target K = 1.
##   Cell 3: True two-class mixture (Study 1 DGP) varied across
##           class separation Delta in {0, 0.5, 1.0, 1.5, 2.0}.
##           Class-count selection target K = 2 at large Delta.
##
## Sample sizes: N in {35, 70, 100}.
##
## Pre-flight reps schedule (kill-switch budget; full-production
## reps are committed only after this run validates feasibility):
##   - Cell 1:  N x 2000 reps  (Type I error, MCSE 0.005 at 0.05)
##   - Cell 2:  N x 500 reps   (class-count selection, MCSE 0.022)
##   - Cell 3:  N x 500 reps   (class-count selection, MCSE 0.022)
##
## Total cells: 3 (cell 1) + 6 (cell 2) + 15 (cell 3) = 24.
## Total replicates: 6000 + 3000 + 7500 = 16500.
## Per-replicate cost: lme_analysis (~0.05 s) + lcmm_analysis
## (K = 1, K = 2 fits, ~3-7 s).
## Wall-clock estimate at 8 cores: 4-6 hours.
##
## Output: analysis/data/study5/study5-results.rds (per-replicate)
##         analysis/data/study5/study5-summary.rds (Morris cells)
##         analysis/data/study5/study5-runlog.txt
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/10-study5-production.R
##
## Optionally NCORES, NREPS_TYPE1, NREPS_KSEL via environment.

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
})

source('analysis/scripts/latent-class-mixture-application/01-lcmm-wrapper.R')
source('analysis/scripts/latent-class-mixture-application/05-marginal-mixture-dgp.R')

DATA_DIR <- file.path('analysis', 'data', 'study5')
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

NCORES <- as.integer(Sys.getenv('NCORES', '8'))
NREPS_TYPE1 <- as.integer(Sys.getenv('NREPS_TYPE1', '2000'))
NREPS_KSEL <- as.integer(Sys.getenv('NREPS_KSEL', '500'))
SEED_BASE <- 20260509L

## ---- Common trial design and parameters ----
td <- buildtrialdesign(
  name_longform  = 'open label',
  name_shortform = 'OL',
  timepoints     = cumulative(rep(2.5, 8)),
  timeptnames    = paste0('OL', 1:8),
  expectancies   = rep(1, 8),
  ondrug         = list(pathA = rep(1, 8))
)
TRIAL_PATHS <- td$trialpaths

RP <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)
BP <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)
OP <- list(
  useDE = TRUE, t_random_slope = FALSE, full_model_out = FALSE,
  carryover_t1half = 1.0, simplecarryover = FALSE,
  carryover_scalefactor = 1
)

## ---- Helper: apply heavy-tailed or skewed transformation to BR ----
##
## Operates on the column-wise BR draws. The marginal SD is
## preserved by re-standardising after the transformation. This
## is a post-hoc reshape of the MVN draw rather than a fully
## non-Gaussian DGP, which is acceptable here because Study 5's
## non-Gaussian cells are diagnostic (does lcmm break under
## mild departure from normality?), not the primary substantive
## DGP.
apply_residual_transform <- function(dat, transform, trialdesign) {
  if (transform == 'gaussian') return(dat)
  tnames <- trialdesign$timeptnames
  for (tp in tnames) {
    br_col <- paste(tp, 'br', sep = '.')
    x <- dat[[br_col]]
    z <- (x - mean(x)) / sd(x)
    if (transform == 't5') {
      ## Standardise to a t_5 cdf percentile, then re-shape.
      u <- pnorm(z)
      u <- pmin(pmax(u, 1e-9), 1 - 1e-9)
      x_new <- qt(u, df = 5)
      ## t_5 has variance 5/(5-2) = 5/3 ~= 1.667; rescale to
      ## SD of the original column to preserve marginal SD.
      x_new <- x_new / sd(x_new) * sd(x) + mean(x)
    } else if (transform == 'skewnorm4') {
      ## Skew-normal with shape alpha = 4 via inversion of a
      ## standardised target. Use direct draw approach: mix
      ## abs(N) tilt onto the existing residual.
      u <- pnorm(z)
      u <- pmin(pmax(u, 1e-9), 1 - 1e-9)
      ## qsn from sn package would be cleaner; use a simple
      ## stochastic representation: skew = delta*|U| + sqrt(1
      ## - delta^2)*V, with U, V independent standard normals.
      ## Here we approximate by adding a |z|-proportional
      ## positive shift then re-standardising.
      delta <- 4 / sqrt(1 + 4^2)  ## delta from shape 4
      u_norm <- abs(qnorm(u))
      v_norm <- z
      x_new <- delta * u_norm + sqrt(1 - delta^2) * v_norm
      x_new <- x_new / sd(x_new) * sd(x) + mean(x)
    } else {
      stop('unknown transform: ', transform)
    }
    dat[, (br_col) := x_new]
  }
  ## Recompute total drug-effect columns and symptom-score
  ## columns (BL - sum of components).
  cl <- c('br', 'tv', 'pb')
  for (tp in tnames) {
    comps <- paste(tp, cl, sep = '.')
    dat[, (paste0('D_', tp)) := rowSums(.SD), .SDcols = comps]
    dat[, (tp) := BL - rowSums(.SD), .SDcols = comps]
  }
  dat
}

## ---- DGP dispatcher ----
##
## cell$dgp_family in {'gaussian_null', 't5_null', 'skewnorm4_null',
##                     'two_class'}
## cell$delta numeric: class separation parameter for two_class.
generate_data_for_cell <- function(cell, rep_idx) {
  set.seed(SEED_BASE + rep_idx + 1e6 * cell$cell_id)
  mp <- data.table(
    N = cell$N, c.bm = 0,
    carryover_t1half = 1.0,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )
  if (cell$dgp_family == 'gaussian_null') {
    dat <- generateData(
      modelparam = mp, respparam = RP, blparam = BP,
      trialdesign = TRIAL_PATHS[[1]],
      empirical = FALSE, makePositiveDefinite = TRUE,
      dgp_architecture = 'mvn'
    )
  } else if (cell$dgp_family == 't5_null') {
    dat <- generateData(
      modelparam = mp, respparam = RP, blparam = BP,
      trialdesign = TRIAL_PATHS[[1]],
      empirical = FALSE, makePositiveDefinite = TRUE,
      dgp_architecture = 'mvn'
    )
    dat <- apply_residual_transform(dat, 't5', TRIAL_PATHS[[1]])
  } else if (cell$dgp_family == 'skewnorm4_null') {
    dat <- generateData(
      modelparam = mp, respparam = RP, blparam = BP,
      trialdesign = TRIAL_PATHS[[1]],
      empirical = FALSE, makePositiveDefinite = TRUE,
      dgp_architecture = 'mvn'
    )
    dat <- apply_residual_transform(dat, 'skewnorm4', TRIAL_PATHS[[1]])
  } else if (cell$dgp_family == 'two_class') {
    mp_mm <- copy(mp)
    mp_mm$beta_R_factor  <- cell$delta  ## responder scaling
    mp_mm$beta_NR_factor <- 0.0
    mp_mm$gating_intercept <- 0.0
    mp_mm$gating_slope <- 1.0
    dat <- generate_marginal_mixture_data(
      modelparam = mp_mm, respparam = RP, blparam = BP,
      trialdesign = TRIAL_PATHS[[1]]
    )
  } else {
    stop('unknown dgp_family: ', cell$dgp_family)
  }
  dat[, path := 1]
  dat
}

## ---- Per-replicate worker ----
run_one_rep <- function(cell, rep_idx) {
  res <- tryCatch({
    dat <- generate_data_for_cell(cell, rep_idx)
    t_lme0 <- Sys.time()
    out_lme <- tryCatch(
      lme_analysis(TRIAL_PATHS, dat, OP),
      error = function(e) data.table(
        beta = NA_real_, betaSE = NA_real_, p = NA_real_,
        issingular = NA, warning = paste0('err: ',
                                          conditionMessage(e))
      )
    )
    t_lme <- as.numeric(difftime(Sys.time(), t_lme0,
                                  units = 'secs'))
    t_lcmm0 <- Sys.time()
    out_lcmm <- tryCatch(
      lcmm_analysis(TRIAL_PATHS, dat, OP),
      error = function(e) data.table(
        beta = NA_real_, betaSE = NA_real_, p = NA_real_,
        beta_within = NA_real_, betaSE_within = NA_real_,
        p_within = NA_real_, delta_bic = NA_real_,
        lrt_stat = NA_real_, lrt_df = NA_integer_,
        bic_ng1 = NA_real_, bic_ng2 = NA_real_,
        issingular = NA, warning = paste0('err: ',
                                          conditionMessage(e)),
        entropy = NA_real_, bic = NA_real_, conv = 0L,
        conv_ng1 = 0L
      )
    )
    t_lcmm <- as.numeric(difftime(Sys.time(), t_lcmm0,
                                   units = 'secs'))
    data.table(
      cell_id = cell$cell_id,
      dgp_family = cell$dgp_family,
      delta = cell$delta, N = cell$N,
      rep_idx = rep_idx,
      lme_beta = out_lme$beta, lme_p = out_lme$p,
      lme_t = t_lme,
      lcmm_gating_beta = out_lcmm$beta,
      lcmm_gating_p = out_lcmm$p,
      lcmm_within_beta = out_lcmm$beta_within,
      lcmm_within_p = out_lcmm$p_within,
      lcmm_delta_bic = out_lcmm$delta_bic,
      lcmm_bic_ng1 = out_lcmm$bic_ng1,
      lcmm_bic_ng2 = out_lcmm$bic_ng2,
      lcmm_entropy = out_lcmm$entropy,
      lcmm_conv = out_lcmm$conv,
      lcmm_t = t_lcmm
    )
  }, error = function(e) {
    data.table(
      cell_id = cell$cell_id,
      dgp_family = cell$dgp_family,
      delta = cell$delta, N = cell$N,
      rep_idx = rep_idx,
      lme_beta = NA_real_, lme_p = NA_real_, lme_t = NA_real_,
      lcmm_gating_beta = NA_real_, lcmm_gating_p = NA_real_,
      lcmm_within_beta = NA_real_, lcmm_within_p = NA_real_,
      lcmm_delta_bic = NA_real_,
      lcmm_bic_ng1 = NA_real_, lcmm_bic_ng2 = NA_real_,
      lcmm_entropy = NA_real_, lcmm_conv = 0L,
      lcmm_t = NA_real_
    )
  })
  res
}

## ---- Build cell list ----
N_VALUES <- c(35, 70, 100)
DELTAS <- c(0, 0.5, 1.0, 1.5, 2.0)

cells <- rbindlist(list(
  ## Cell 1: Gaussian null
  CJ(dgp_family = 'gaussian_null', delta = 0,
     N = N_VALUES, n_reps = NREPS_TYPE1),
  ## Cell 2a: t5
  CJ(dgp_family = 't5_null', delta = 0,
     N = N_VALUES, n_reps = NREPS_KSEL),
  ## Cell 2b: skew-normal
  CJ(dgp_family = 'skewnorm4_null', delta = 0,
     N = N_VALUES, n_reps = NREPS_KSEL),
  ## Cell 3: two-class
  CJ(dgp_family = 'two_class', delta = DELTAS,
     N = N_VALUES, n_reps = NREPS_KSEL)
))
cells[, cell_id := seq_len(.N)]
setcolorder(cells, c('cell_id', 'dgp_family', 'delta',
                     'N', 'n_reps'))

cat(sprintf('Study 5 production driver started at %s\n',
            format(Sys.time())))
cat(sprintf('Cores: %d\n', NCORES))
cat(sprintf('Total cells: %d\n', nrow(cells)))
cat(sprintf('Total replicates: %d\n', sum(cells$n_reps)))
cat('Cell breakdown:\n')
print(cells)

## ---- Build job list ----
jobs <- vector('list', sum(cells$n_reps))
k <- 1L
for (i in seq_len(nrow(cells))) {
  cell <- as.list(cells[i])
  for (r in seq_len(cell$n_reps)) {
    jobs[[k]] <- list(cell = cell, rep_idx = r)
    k <- k + 1L
  }
}
cat(sprintf('Total jobs queued: %d\n', length(jobs)))

## ---- Run via mclapply with chunked checkpointing ----
worker <- function(job) run_one_rep(job$cell, job$rep_idx)

CHUNK_SIZE <- 500L
n_chunks <- ceiling(length(jobs) / CHUNK_SIZE)
results_list <- vector('list', n_chunks)

t_run0 <- Sys.time()
for (ch in seq_len(n_chunks)) {
  i0 <- (ch - 1L) * CHUNK_SIZE + 1L
  i1 <- min(ch * CHUNK_SIZE, length(jobs))
  chunk_jobs <- jobs[i0:i1]
  t_ch0 <- Sys.time()
  chunk_res <- mclapply(chunk_jobs, worker,
                         mc.cores = NCORES,
                         mc.preschedule = FALSE,
                         mc.set.seed = FALSE)
  results_list[[ch]] <- rbindlist(chunk_res, fill = TRUE)
  t_ch <- as.numeric(difftime(Sys.time(), t_ch0, units = 'secs'))
  t_total <- as.numeric(difftime(Sys.time(), t_run0,
                                   units = 'secs'))
  cat(sprintf('  chunk %d/%d: %d jobs in %.1f s; cumulative %.1f min\n',
              ch, n_chunks, length(chunk_jobs), t_ch, t_total / 60))
  ## Progressive save after each chunk.
  saveRDS(rbindlist(results_list[seq_len(ch)], fill = TRUE),
          file.path(DATA_DIR,
                     sprintf('study5-progress.rds')))
}

results <- rbindlist(results_list, fill = TRUE)
t_run <- as.numeric(difftime(Sys.time(), t_run0, units = 'secs'))

## ---- Final outputs ----
saveRDS(results, file.path(DATA_DIR, 'study5-results.rds'))
cat(sprintf('\nPer-replicate output -> %s (rows=%d)\n',
            file.path(DATA_DIR, 'study5-results.rds'),
            nrow(results)))

## Morris-format cell summary.
binom_mcse <- function(rate, n) {
  if (is.na(rate) || is.na(n) || n == 0) return(NA_real_)
  sqrt(rate * (1 - rate) / n)
}

summary_tab <- results[, {
  n_total <- .N
  n_lme <- sum(!is.na(lme_p))
  n_lcmm <- sum(!is.na(lcmm_gating_p))
  conv_lme <- mean(!is.na(lme_p))
  conv_lcmm <- mean(lcmm_conv == 1, na.rm = TRUE)
  rej_lme <- mean(lme_p < 0.05, na.rm = TRUE)
  rej_lcmm_gating <- mean(lcmm_gating_p < 0.05, na.rm = TRUE)
  k2_pref <- mean(lcmm_delta_bic > 0, na.rm = TRUE)
  k2_strong <- mean(lcmm_delta_bic > 6, na.rm = TRUE)
  list(
    n_reps = n_total,
    conv_lme = conv_lme, conv_lcmm = conv_lcmm,
    rej_lme = rej_lme,
    rej_lcmm_gating = rej_lcmm_gating,
    p_select_k2 = k2_pref,
    p_select_k2_strong = k2_strong,
    mean_lme_t = mean(lme_t, na.rm = TRUE),
    mean_lcmm_t = mean(lcmm_t, na.rm = TRUE),
    median_lcmm_t = median(lcmm_t, na.rm = TRUE),
    mcse_rej_lme = binom_mcse(rej_lme, n_lme),
    mcse_rej_lcmm = binom_mcse(rej_lcmm_gating, n_lcmm),
    mcse_k2 = binom_mcse(k2_pref, n_lcmm)
  )
}, by = .(cell_id, dgp_family, delta, N)]

saveRDS(summary_tab, file.path(DATA_DIR, 'study5-summary.rds'))
cat(sprintf('Cell summary -> %s\n',
            file.path(DATA_DIR, 'study5-summary.rds')))
cat('\nMorris cell summary:\n')
print(summary_tab)

## Run log.
log_path <- file.path(DATA_DIR, 'study5-runlog.txt')
writeLines(c(
  sprintf('run_started: %s', format(Sys.time())),
  sprintf('NCORES: %d', NCORES),
  sprintf('NREPS_TYPE1: %d', NREPS_TYPE1),
  sprintf('NREPS_KSEL: %d', NREPS_KSEL),
  sprintf('total_cells: %d', nrow(cells)),
  sprintf('total_reps: %d', sum(cells$n_reps)),
  sprintf('total_wall_sec: %.1f', t_run),
  sprintf('total_wall_hr: %.2f', t_run / 3600),
  sprintf('mean_lcmm_fit_s: %.2f',
          mean(results$lcmm_t, na.rm = TRUE))
), log_path)
cat(sprintf('Run log -> %s\n', log_path))
cat(sprintf('\nTotal wall: %.1f min (%.2f hr)\n',
            t_run / 60, t_run / 3600))
cat('Done.\n')
