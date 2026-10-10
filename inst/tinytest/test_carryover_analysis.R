## analyze_trial_extended(): the predicted drug effect and the clinical
## threshold must be on the scale of the fitted model. lme_analysis()
## fits the uncentered biomarker, so a biomarker with a nonzero mean
## distinguishes a correct implementation from one that centers it.

make_co_data <- function(bm_mean = 100, bm_sd = 15, seed = 42) {
  td <- buildtrialdesign(
    name_longform = 'crossover',
    name_shortform = 'CO',
    timepoints = cumulative(rep(2.5, 8)),
    timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
    expectancies = rep(0.5, 8),
    ondrug = list(
      pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
      pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
    )
  )
  mp <- data.table(
    N = 40, c.bm = 0.3,
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
  bp <- data.table(cat = c('bm', 'BL'), m = c(bm_mean, 70),
                   sd = c(bm_sd, 10))
  set.seed(seed)
  dat1 <- generateData(mp, rp, bp, td$trialpaths[[1]],
                       empirical = FALSE, makePositiveDefinite = TRUE)
  dat1[, path := 1]
  dat2 <- generateData(mp, rp, bp, td$trialpaths[[2]],
                       empirical = FALSE, makePositiveDefinite = TRUE)
  dat2[, path := 2]
  list(td = td, dat = rbind(dat1, dat2, fill = TRUE))
}

{
  s <- make_co_data()
  res <- analyze_trial_extended(s$td$trialpaths, s$dat, threshold = -3)
  expect_false(is.null(res$model), info = 'extended analysis fitted a model')

  fit <- res$model
  d <- res$data
  b_values <- c(res$biomarker_mean, res$biomarker_mean + res$biomarker_sd)
  for (b in b_values) {
    nd <- d[c(1, 1)]
    nd[, `:=`(bm = b, Dbc = c(1, 0))]
    contrast <- diff(rev(as.numeric(predict(fit, newdata = nd, level = 0))))
    expect_equal(res$predicted_effect(b), contrast, tolerance = 1e-8,
                 info = sprintf('predicted_effect(%.1f) equals the model contrast Dbc 1 vs 0', b))
  }

  expect_equal(res$biomarker_sd, sd(d$bm), tolerance = 1e-10,
               info = 'biomarker_sd is the SD of the analysis data')
  expect_equal(res$predicted_effect(res$threshold_bm), -3, tolerance = 1e-8,
               info = 'the predicted effect at threshold_bm equals the threshold')
}
