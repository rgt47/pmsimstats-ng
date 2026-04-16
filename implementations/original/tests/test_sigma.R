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
  result <- buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  expect_true(is.positive.definite(result$sigma),
              info = 'buildSigma returns PD matrix for OL design')
}

{
  td <- make_ol_design()
  p <- make_default_params()
  result <- buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
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
  result <- buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
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
  result <- buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
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
  result <- buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
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
  result <- buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
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
  result <- buildSigma(p$mp, p$rp, p$bp, td$trialpaths[[1]])
  expect_true(!is.null(result$chol_sigma),
              info = 'buildSigma returns Cholesky factor')
  expect_equal(nrow(result$chol_sigma), 26L,
               info = 'Cholesky factor has 26 rows')
  reconstructed <- t(result$chol_sigma) %*% result$chol_sigma
  expect_equal(reconstructed, result$sigma, tolerance = 1e-8,
               info = 'Cholesky factor reconstructs sigma')
}
