make_analysis_setup <- function(N = 35, c.bm = 0.3) {
  td <- buildtrialdesign(
    name_longform = 'open label',
    name_shortform = 'OL',
    timepoints = cumulative(rep(2.5, 8)),
    timeptname = paste0('OL', 1:8),
    expectancies = rep(1, 8),
    ondrug = list(pathA = rep(1, 8))
  )
  mp <- data.table(
    N = N, c.bm = c.bm,
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

test_that('lme_analysis returns expected columns', {
  s <- make_analysis_setup(N = 35)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  dat[, path := 1]
  result <- lme_analysis(s$td$trialpaths, dat)
  expect_true('beta' %in% names(result))
  expect_true('betaSE' %in% names(result))
  expect_true('p' %in% names(result))
})

test_that('lme_analysis p-value is between 0 and 1', {
  s <- make_analysis_setup(N = 35)
  set.seed(42)
  dat <- generateData(s$mp, s$rp, s$bp, s$td$trialpaths[[1]],
                      empirical = FALSE, makePositiveDefinite = TRUE)
  dat[, path := 1]
  result <- lme_analysis(s$td$trialpaths, dat)
  if (!is.na(result$p)) {
    expect_true(result$p >= 0 && result$p <= 1)
  }
})

test_that('lme_analysis works with CO design (varInDb)', {
  td <- buildtrialdesign(
    name_longform = 'crossover',
    name_shortform = 'CO',
    timepoints = cumulative(rep(2.5, 8)),
    timeptname = c(paste0('COa', 1:4), paste0('COb', 1:4)),
    expectancies = rep(0.5, 8),
    ondrug = list(
      pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
      pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
    )
  )
  mp <- data.table(
    N = 35, c.bm = 0.3,
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

  set.seed(42)
  dat1 <- generateData(mp, rp, bp, td$trialpaths[[1]],
                       empirical = FALSE, makePositiveDefinite = TRUE)
  dat1[, path := 1]
  dat2 <- generateData(mp, rp, bp, td$trialpaths[[2]],
                       empirical = FALSE, makePositiveDefinite = TRUE)
  dat2[, path := 2]
  dat <- rbind(dat1, dat2, fill = TRUE)

  result <- lme_analysis(td$trialpaths, dat)
  expect_true('beta' %in% names(result))
  if (!is.na(result$p)) {
    expect_true(result$p >= 0 && result$p <= 1)
  }
})
