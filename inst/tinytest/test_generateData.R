make_ol_setup <- function(N = 35) {
  td <- buildtrialdesign(
    name_longform = 'open label',
    name_shortform = 'OL',
    timepoints = cumulative(rep(2.5, 8)),
    timeptnames = paste0('OL', 1:8),
    expectancies = rep(1, 8),
    ondrug = list(pathA = rep(1, 8))
  )
  mp <- data.table(
    N = N, c.bm = 0.3,
    carryover_t1half = 0,
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
  list(td = td, mp = mp, rp = rp, bp = bp)
}

{
  s <- make_ol_setup(N = 50)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  expect_true(inherits(dat, 'data.table'),
              info = 'generateData returns data.table')
  expect_equal(nrow(dat), 50L,
               info = 'generateData row count matches N')
}

{
  s <- make_ol_setup(N = 20)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  expect_true('bm' %in% names(dat),
              info = 'generateData produces bm column')
  expect_true('BL' %in% names(dat),
              info = 'generateData produces BL column')
  expect_true('ptID' %in% names(dat),
              info = 'generateData produces ptID column')
  expect_true('OL1' %in% names(dat),
              info = 'generateData produces first outcome column')
  expect_true('OL8' %in% names(dat),
              info = 'generateData produces last outcome column')
}

{
  s <- make_ol_setup(N = 30)
  sigma_obj <- pmsimstats:::buildSigma(s$mp, s$rp, s$bp, s$td$trialpaths[[1]])
  set.seed(99)
  dat1 <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                       empirical = FALSE, makePositiveDefinite = TRUE,
                       seed = NA)
  set.seed(99)
  dat2 <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                       empirical = FALSE, makePositiveDefinite = TRUE,
                       cached_sigma = sigma_obj)
  expect_equal(dat1$bm, dat2$bm,
               info = 'cached_sigma path matches direct call: bm column')
  expect_equal(dat1$OL1, dat2$OL1,
               info = 'cached_sigma path matches direct call: OL1 outcome')
}

{
  s <- make_ol_setup(N = 10)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  factor_sum <- dat$OL1.tv + dat$OL1.pb + dat$OL1.br
  expect_equal(dat$OL1, dat$BL - factor_sum, tolerance = 1e-10,
               info = 'outcome column equals BL minus factor sum')
}

{
  # b.bm.tv: optional biomarker-TV mean shift (prog-TV-grow cell,
  # docs/57 Section 7.3). Absent or zero must leave TV unshifted;
  # a nonzero value shifts every TV column by exactly
  # b.bm.tv * bm_z * sigma_TV and leaves BR untouched.
  s <- make_ol_setup(N = 200)
  set.seed(7)
  dat_default <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                               empirical = FALSE, makePositiveDefinite = TRUE)

  s$mp$b.bm.tv <- 2
  set.seed(7)
  dat_shifted <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                               empirical = FALSE, makePositiveDefinite = TRUE)

  bm_z <- (dat_default$bm - s$bp[cat == 'bm']$m) / s$bp[cat == 'bm']$sd
  expected_shift <- 2 * bm_z * s$rp[cat == 'tv']$sd
  expect_equal(dat_shifted$OL1.tv - dat_default$OL1.tv, expected_shift,
               tolerance = 1e-10,
               info = 'b.bm.tv = 2 shifts TV by exactly b.bm.tv * bm_z * sigma_TV')
  expect_equal(dat_shifted$OL1.br, dat_default$OL1.br,
               info = 'b.bm.tv leaves BR untouched')

  s$mp$b.bm.tv <- 0
  set.seed(7)
  dat_zero <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                           empirical = FALSE, makePositiveDefinite = TRUE)
  expect_equal(dat_zero$OL1.tv, dat_default$OL1.tv,
               info = 'b.bm.tv = 0 reproduces the unshifted default exactly')
}
