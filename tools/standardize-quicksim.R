## tools/standardize-quicksim.R
##
## One-pass standardisation of the nine quick-sim summary tables.
## Reads each `<slug>-replicates.rds` from
## `analysis/data/quick-sim/`, recomputes the per-cell summary
## with canonical Morris columns, and overwrites
## `<slug>-summary.txt`.
##
## Canonical column structure:
##   <axes>... | n_reps | rejection_rate | mcse_rate
##           | mean_beta | sd_beta | mcse_mean_beta
##           | coverage | mcse_coverage
##           | conv_rate
##
## Coverage convention: nominal-95% CI coverage of the cell's
## empirical mean (Morris's model-SE calibration check). For
## each replicate, does [beta - 1.96*SE, beta + 1.96*SE] contain
## the cell's mean(beta)? Coverage = mean of indicator.
##
## For paper 08 the per-replicate output preserves no betaSE so
## coverage and mean_beta-MCSE are reported as NA.
##
## Run: Rscript tools/standardize-quicksim.R
##      from the project root.

suppressPackageStartupMessages({
  library(data.table)
})

QSDIR <- 'analysis/data/quick-sim'
stopifnot(dir.exists(QSDIR))

mcse_prop <- function(p, n) sqrt(p * (1 - p) / n)
mcse_mean <- function(s, n) s / sqrt(n)
mcse_cov <- function(c, n) sqrt(c * (1 - c) / n)

## Per-cell summariser. Takes one cell's worth of replicates
## (a data.table with at least: beta, betaSE, p, converged) and
## returns a one-row data.table with the canonical Morris columns.
summarise_cell <- function(d, alpha = 0.05) {
  d <- copy(as.data.table(d))
  if (!('converged' %in% names(d))) {
    d[, converged := TRUE]
  }
  if (!('betaSE' %in% names(d))) {
    d[, betaSE := NA_real_]
  }
  n_reps <- nrow(d)
  conv <- d$converged
  ## Use only converged replicates for performance measures.
  use <- d[as.logical(conv)]
  n_conv <- nrow(use)
  rej <- if (n_conv == 0) NA_real_ else mean(use$p < alpha, na.rm = TRUE)

  if (all(is.na(use$beta))) {
    mb <- NA_real_; sb <- NA_real_; mcse_b <- NA_real_
  } else {
    mb <- mean(use$beta, na.rm = TRUE)
    sb <- sd(use$beta, na.rm = TRUE)
    mcse_b <- mcse_mean(sb, n_conv)
  }

  if (all(is.na(use$betaSE)) || is.na(mb)) {
    cov_emp <- NA_real_; mcse_cov_v <- NA_real_
  } else {
    use[, ci_lo := beta - 1.96 * betaSE]
    use[, ci_hi := beta + 1.96 * betaSE]
    use[, covers_mean := !is.na(ci_lo) & !is.na(ci_hi) &
                          ci_lo <= mb & mb <= ci_hi]
    cov_emp <- mean(use$covers_mean, na.rm = TRUE)
    mcse_cov_v <- mcse_cov(cov_emp, n_conv)
  }

  data.table(
    n_reps = n_reps,
    rejection_rate = rej,
    mcse_rate = if (is.na(rej)) NA_real_ else mcse_prop(rej, n_conv),
    mean_beta = mb,
    sd_beta = sb,
    mcse_mean_beta = mcse_b,
    coverage = cov_emp,
    mcse_coverage = mcse_cov_v,
    conv_rate = mean(as.logical(conv), na.rm = TRUE)
  )
}

## Wrapper that handles the per-paper column-name differences.
standardise_paper <- function(slug, axes, beta_col, se_col, p_col,
                               conv_col = NULL) {
  reps_path <- file.path(QSDIR, paste0(slug, '-replicates.rds'))
  d <- as.data.table(readRDS(reps_path))
  ## Rename to canonical names for `summarise_cell`.
  setnames(d, beta_col, 'beta')
  if (!is.null(se_col) && se_col %in% names(d)) {
    setnames(d, se_col, 'betaSE')
  }
  setnames(d, p_col, 'p')
  if (!is.null(conv_col) && conv_col %in% names(d)) {
    setnames(d, conv_col, 'converged')
  }
  out <- d[, summarise_cell(.SD), by = axes,
           .SDcols = intersect(names(d),
                                c('beta', 'betaSE', 'p', 'converged'))]
  out_path <- file.path(QSDIR, paste0(slug, '-summary.txt'))
  fwrite(out, out_path, sep = '\t', quote = FALSE,
         na = 'NA')
  message(sprintf('  %s: %d cells -> %s', slug, nrow(out), out_path))
  invisible(out)
}

cat('Standardising quick-sim summaries:\n')

## Paper 01: dgp-mean-moderation-vs-mvn
standardise_paper('01-dgp',
                   axes = c('architecture', 'c.bm', 't1half'),
                   beta_col = 'beta', se_col = 'betaSE', p_col = 'p',
                   conv_col = NULL)

## Paper 02: carryover-sensitivity
standardise_paper('02-carryover',
                   axes = c('analysis_spec', 'c.bm'),
                   beta_col = 'beta', se_col = 'betaSE', p_col = 'p',
                   conv_col = 'converged')

## Paper 03: latent-class-mixture (3 tests per cell, but the per-
## replicate file has one row per analysis-strategy, with `beta`
## meaning bm:Dbc for lme rows and gating-coefficient for lcmm rows;
## summary therefore reports per (gating_slope, analysis) cell,
## with the rejection-rate column reflecting whichever analysis
## the row records. Three-test summary will be produced separately.)
slug03 <- '03-latent'
d03 <- as.data.table(readRDS(file.path(QSDIR,
                              paste0(slug03, '-replicates.rds'))))
setnames(d03, 'conv', 'converged')
## Build three rows per gating_slope: lme bm:Dbc, lcmm gating
## unconditional, lcmm BIC-conditional gating. Per-group summarise
## then bind, prepending the test name.
do_test <- function(filter_expr, p_col, test_name) {
  sub <- d03[eval(filter_expr)]
  sub_dt <- data.table(
    beta = sub$beta, betaSE = sub$betaSE,
    p = sub[[p_col]], converged = sub$converged,
    gating_slope = sub$gating_slope
  )
  res <- sub_dt[, summarise_cell(.SD), by = .(gating_slope)]
  res[, test := test_name]
  res
}
out03_lme <- d03[analysis == 'lme',
  summarise_cell(data.table(beta = beta, betaSE = betaSE,
                            p = p, converged = converged)),
  by = .(gating_slope)]
out03_lme[, test := 'lme_bmDbc']
out03_gate <- d03[analysis == 'lcmm',
  summarise_cell(data.table(beta = beta, betaSE = betaSE,
                            p = gating_p, converged = converged)),
  by = .(gating_slope)]
out03_gate[, test := 'lcmm_gating_unconditional']
out03_bic <- d03[analysis == 'lcmm',
  summarise_cell(data.table(beta = beta, betaSE = betaSE,
                            p = ifelse(delta_bic > 0, gating_p, 1),
                            converged = converged)),
  by = .(gating_slope)]
out03_bic[, test := 'lcmm_gating_BIC_conditional']
out03 <- rbindlist(list(out03_lme, out03_gate, out03_bic))
setcolorder(out03, c('gating_slope', 'test'))
out03 <- out03[order(gating_slope, test)]
fwrite(out03, file.path(QSDIR, paste0(slug03, '-summary.txt')),
       sep = '\t', quote = FALSE, na = 'NA')
message(sprintf('  %s: %d cells -> %s-summary.txt', slug03,
                nrow(out03), slug03))

## Paper 04: treatment-main-effect
standardise_paper('04-mainfx',
                   axes = c('design', 'true_effect'),
                   beta_col = 'beta_main', se_col = 'betaSE_main',
                   p_col = 'p_main', conv_col = 'converged')

## Paper 05: nof1-design-sensitivity
standardise_paper('05-design',
                   axes = c('c_br'),
                   beta_col = 'beta_Dbc', se_col = 'betaSE_Dbc',
                   p_col = 'p_Dbc', conv_col = NULL)

## Paper 06: component-decomposition
standardise_paper('06-decomp',
                   axes = c('pb_strength', 'analysis'),
                   beta_col = 'beta_bmDbc', se_col = 'betaSE_bmDbc',
                   p_col = 'p_bmDbc', conv_col = 'converged')

## Paper 07: gompertz-evaluation
standardise_paper('07-gompertz',
                   axes = c('family'),
                   beta_col = 'beta_bmDbc', se_col = 'betaSE_bmDbc',
                   p_col = 'p_bmDbc', conv_col = 'converged')

## Paper 08: test-procedure (no betaSE; coverage will be NA)
standardise_paper('08-test',
                   axes = c('procedure', 'c.bm'),
                   beta_col = 'beta', se_col = NULL,
                   p_col = 'p_interaction', conv_col = 'converged')

## Paper 09: informative-dropout
standardise_paper('09-dropout',
                   axes = c('design', 'dropout_pattern'),
                   beta_col = 'beta_bmDbc', se_col = 'betaSE_bmDbc',
                   p_col = 'p_bmDbc', conv_col = NULL)

cat('Done.\n')
