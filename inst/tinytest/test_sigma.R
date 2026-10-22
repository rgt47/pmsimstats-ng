make_ol_design <- function() {
  buildtrialdesign(
    name_longform = 'open label',
    name_shortform = 'OL',
    timepoints = cumulative(rep(2.5, 8)),
    timeptnames = paste0('OL', 1:8),
    expectancies = rep(1, 8),
    ondrug = list(pathA = rep(1, 8))
  )
}

make_default_params <- function(c.bm = 0.3, carryover_t1half = 0) {
  mp <- data.table(
    N = 35, c.bm = c.bm,
    carryover_t1half = carryover_t1half,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )
  rp <- data.table(
    cat = c('tv', 'pb', 'br'),
    max = c(10.98604, 6.50647, 10.98604),
    disp = c(5, 5, 5),
    rate = c(0.42, 0.35, 0.42),
    sd = c(5, 2, 5)
  )
  bp <- data.table(
    cat = c('bm', 'BL'),
    m = c(0, 70),
    sd = c(1, 10)
  )
  list(mp = mp, rp = rp, bp = bp)
}

{
  td <- make_ol_design()
  p <- make_default_params()
  result <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  expect_true(is.positive.definite(result$sigma),
              info = 'buildSigma returns PD matrix for OL design')
}

{
  td <- make_ol_design()
  p <- make_default_params()
  result <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  expect_equal(nrow(result$sigma), 26L,
               info = 'buildSigma: 26 rows for 8-timepoint OL')
  expect_equal(ncol(result$sigma), 26L,
               info = 'buildSigma: 26 cols for 8-timepoint OL')
  expect_equal(length(result$labels), 26L,
               info = 'buildSigma: 26 labels')
  expect_equal(length(result$means), 26L,
               info = 'buildSigma: 26 means')
}

{
  td <- make_ol_design()
  p <- make_default_params()
  result <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  expect_equal(result$labels[1], 'bm',
               info = 'buildSigma label[1] is bm')
  expect_equal(result$labels[2], 'BL',
               info = 'buildSigma label[2] is BL')
  expect_equal(result$labels[3], 'OL1.tv',
               info = 'buildSigma label[3] is first tv')
  expect_equal(result$labels[10], 'OL8.tv',
               info = 'buildSigma label[10] is last tv')
  expect_equal(result$labels[11], 'OL1.pb',
               info = 'buildSigma label[11] is first pb')
  expect_equal(result$labels[19], 'OL1.br',
               info = 'buildSigma label[19] is first br')
}

{
  td <- make_ol_design()
  p <- make_default_params(c.bm = 0)
  result <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  bm_idx <- 1
  br_idx <- which(grepl('\\.br$', result$labels))
  for (i in br_idx) {
    corr_bm_br <- result$sigma[bm_idx, i] /
      (sqrt(result$sigma[bm_idx, bm_idx]) * sqrt(result$sigma[i, i]))
    expect_equal(corr_bm_br, 0,
                 info = paste0('c.bm=0 yields zero BM-BR correlation at ',
                                result$labels[i]))
  }
}

{
  td <- make_ol_design()
  p <- make_default_params()
  result <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  corr <- cov2cor(result$sigma)
  tv1 <- which(result$labels == 'OL1.tv')
  tv2 <- which(result$labels == 'OL2.tv')
  tv8 <- which(result$labels == 'OL8.tv')
  expect_true(abs(corr[tv1, tv2]) > abs(corr[tv1, tv8]),
              info = 'AR(1) correlation decays with time gap')
}

{
  td <- buildtrialdesign(
    name_longform = 'crossover',
    name_shortform = 'CO',
    timepoints = cumulative(rep(2.5, 8)),
    timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
    expectancies = rep(0.5, 8),
    ondrug = list(
      pathA = c(1, 1, 1, 1, 0, 0, 0, 0)
    )
  )
  p <- make_default_params(c.bm = 0.3, carryover_t1half = 1.0)
  result <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  corr <- cov2cor(result$sigma)
  bm_idx <- 1
  br_on <- which(result$labels == 'COa1.br')
  br_off1 <- which(result$labels == 'COb1.br')
  br_off4 <- which(result$labels == 'COb4.br')
  expect_true(abs(corr[bm_idx, br_on]) > abs(corr[bm_idx, br_off1]),
              info = 'BM-BR correlation attenuated at first off-drug tp')
  expect_true(abs(corr[bm_idx, br_off1]) > abs(corr[bm_idx, br_off4]),
              info = 'BM-BR correlation decays across off-drug tps')
}

{
  td <- make_ol_design()
  p <- make_default_params()
  result <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  expect_true(!is.null(result$chol_sigma),
              info = 'buildSigma returns Cholesky factor')
  expect_equal(nrow(result$chol_sigma), 26L,
               info = 'Cholesky factor has 26 rows')
  reconstructed <- t(result$chol_sigma) %*% result$chol_sigma
  expect_equal(reconstructed, result$sigma, tolerance = 1e-8,
               info = 'Cholesky factor reconstructs sigma')
}

{
  # stopOnNonPD: the refined-process recommendation (paper 01, Table
  # tab:spec and Appendix B.2) is to raise an error instead of
  # silently repairing a non-positive-definite covariance matrix.
  td <- make_ol_design()
  p_nonpd <- make_default_params(c.bm = 0.6, carryover_t1half = 0)
  expect_true(!is.positive.definite(
    pmsimstats:::buildSigma(p_nonpd$mp, p_nonpd$rp, p_nonpd$bp,
                            td$trialpaths[[1]],
                            makePositiveDefinite = FALSE)$sigma),
    info = 'c.bm = 0.6 produces a non-PD matrix at this design/rho (precondition for the tests below)')

  expect_error(
    pmsimstats:::buildSigma(p_nonpd$mp, p_nonpd$rp, p_nonpd$bp,
                            td$trialpaths[[1]], stopOnNonPD = TRUE),
    info = 'stopOnNonPD = TRUE raises an error on a non-PD matrix instead of repairing it')

  repaired <- pmsimstats:::buildSigma(p_nonpd$mp, p_nonpd$rp, p_nonpd$bp,
                                      td$trialpaths[[1]],
                                      stopOnNonPD = FALSE,
                                      makePositiveDefinite = TRUE)
  expect_true(is.positive.definite(repaired$sigma),
              info = 'stopOnNonPD = FALSE (default) preserves the prior silent-repair behavior')

  p_pd <- make_default_params()
  expect_silent(
    pmsimstats:::buildSigma(p_pd$mp, p_pd$rp, p_pd$bp,
                            td$trialpaths[[1]], stopOnNonPD = TRUE))
}

{
  # c.bm.bl: optional biomarker-baseline-severity correlation
  # (docs/57 Section 7.3, stress-test cells 40-41). Absent or zero
  # must leave bm and BL uncorrelated, as before; a nonzero value
  # sets exactly that correlation and nothing else, and the matrix
  # stays positive definite at a moderate coupling.
  td <- make_ol_design()
  p <- make_default_params()

  result_default <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp,
                                             td$trialpaths[[1]])
  expect_equal(result_default$sigma['BL', 'bm'], 0,
               info = 'bm-BL correlation is 0 when c.bm.bl is absent')

  p$mp$c.bm.bl <- 0.3
  result_coupled <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp,
                                            td$trialpaths[[1]],
                                            stopOnNonPD = TRUE)
  corr_bl <- result_coupled$sigma['BL', 'bm'] /
    sqrt(result_coupled$sigma['BL', 'BL'] * result_coupled$sigma['bm', 'bm'])
  expect_equal(corr_bl, 0.3, tolerance = 1e-10,
               info = 'c.bm.bl = 0.3 sets exactly that bm-BL correlation')

  br_cols <- grepl('br$', result_default$labels)
  expect_equal(result_coupled$sigma['bm', br_cols],
               result_default$sigma['bm', br_cols],
               info = 'c.bm.bl leaves bm-BR couplings untouched')

  p$mp$c.bm.bl <- 0
  result_zero <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp,
                                         td$trialpaths[[1]])
  expect_equal(result_zero$sigma['BL', 'bm'], 0,
               info = 'c.bm.bl = 0 reproduces the uncoupled default exactly')
}

{
  # c.bm.tv: optional biomarker-TV correlation (prog-TV-const cell,
  # docs/57 Section 7.3). Absent or zero must leave bm and TV
  # uncorrelated, as before; a nonzero value sets exactly that
  # correlation at every timepoint and nothing else.
  td <- make_ol_design()
  p <- make_default_params()

  result_default <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp,
                                             td$trialpaths[[1]])
  tv_cols <- grepl('tv$', result_default$labels)
  expect_true(all(result_default$sigma['bm', tv_cols] == 0),
              info = 'bm-TV correlation is 0 when c.bm.tv is absent')

  p$mp$c.bm.tv <- 0.3
  result_coupled <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp,
                                            td$trialpaths[[1]],
                                            stopOnNonPD = TRUE)
  corr_tv <- result_coupled$sigma['bm', tv_cols] /
    sqrt(result_coupled$sigma['bm', 'bm'] *
           diag(result_coupled$sigma)[tv_cols])
  expect_equal(unname(corr_tv), rep(0.3, sum(tv_cols)), tolerance = 1e-10,
               info = 'c.bm.tv = 0.3 sets exactly that bm-TV correlation at every timepoint')

  br_cols <- grepl('br$', result_default$labels)
  expect_equal(result_coupled$sigma['bm', br_cols],
               result_default$sigma['bm', br_cols],
               info = 'c.bm.tv leaves bm-BR couplings untouched')

  p$mp$c.bm.tv <- 0
  result_zero <- pmsimstats:::buildSigma(p$mp, p$rp, p$bp,
                                         td$trialpaths[[1]])
  expect_true(all(result_zero$sigma['bm', tv_cols] == 0),
              info = 'c.bm.tv = 0 reproduces the uncoupled default exactly')
}
