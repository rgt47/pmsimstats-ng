make_ol_setup <- function(N = 35) {
  td <- buildtrialdesign(
    name_longform = 'open label',
    name_shortform = 'OL',
    timepoints = cumulative(rep(2.5, 8)),
    timeptname = paste0('OL', 1:8),
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

test_that('generateData returns data.table with correct rows', {
  s <- make_ol_setup(N = 50)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  expect_s3_class(dat, 'data.table')
  expect_equal(nrow(dat), 50)
})

test_that('generateData has expected columns', {
  s <- make_ol_setup(N = 20)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  expect_true('bm' %in% names(dat))
  expect_true('BL' %in% names(dat))
  expect_true('ptID' %in% names(dat))
  expect_true('OL1' %in% names(dat))
  expect_true('OL8' %in% names(dat))
})

test_that('generateData with cached_sigma matches direct call', {
  s <- make_ol_setup(N = 30)
  sigma_obj <- buildSigma(s$mp, s$rp, s$bp, s$td$trialpaths[[1]])
  set.seed(99)
  dat1 <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                       empirical = FALSE, makePositiveDefinite = TRUE,
                       seed = NA)
  set.seed(99)
  dat2 <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                       empirical = FALSE, makePositiveDefinite = TRUE,
                       cached_sigma = sigma_obj)
  expect_equal(dat1$bm, dat2$bm)
  expect_equal(dat1$OL1, dat2$OL1)
})

test_that('generateData outcome columns equal BL minus factor sums', {
  s <- make_ol_setup(N = 10)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  factor_sum <- dat$OL1.tv + dat$OL1.pb + dat$OL1.br
  expect_equal(dat$OL1, dat$BL - factor_sum, tolerance = 1e-10)
})
