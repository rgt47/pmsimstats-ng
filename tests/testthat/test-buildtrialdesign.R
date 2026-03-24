test_that('OL design has correct structure', {
  td <- buildtrialdesign(
    name_longform = 'open label',
    name_shortform = 'OL',
    timepoints = cumulative(rep(2.5, 8)),
    timeptname = paste0('OL', 1:8),
    expectancies = rep(1, 8),
    ondrug = list(pathA = rep(1, 8))
  )
  expect_equal(td$metadata$name_shortform, 'OL')
  expect_length(td$trialpaths, 1)
  path <- td$trialpaths[[1]]
  expect_equal(nrow(path), 8)
  expect_true(all(path$tod > 0))
  expect_true(all(path$tsd == 0))
})

test_that('CO design has two paths', {
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
  expect_length(td$trialpaths, 2)
  pa <- td$trialpaths[[1]]
  pb <- td$trialpaths[[2]]
  expect_true(all(pa$tod[1:4] > 0))
  expect_true(all(pa$tod[5:8] == 0))
  expect_true(all(pb$tod[1:4] == 0))
  expect_true(all(pb$tod[5:8] > 0))
})

test_that('CO path A accumulates tsd in off-drug phase', {
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
  pa <- td$trialpaths[[1]]
  expect_equal(pa$tsd[5:8], c(2.5, 5, 7.5, 10))
})

test_that('Hybrid design has four paths', {
  td <- buildtrialdesign(
    name_longform = 'hybrid',
    name_shortform = 'HY',
    timepoints = cumulative(c(4, 4, 1, 1, 1, 1, 4, 4)),
    timeptname = c(paste0('HYa', 1:4), paste0('HYb', 1:4)),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    ondrug = list(
      pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
      pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
      pathC = c(1, 1, 0, 1, 1, 0, 1, 0),
      pathD = c(1, 1, 0, 1, 1, 0, 0, 1)
    )
  )
  expect_length(td$trialpaths, 4)
})

test_that('tpb accumulates only during positive expectancy', {
  td <- buildtrialdesign(
    name_longform = 'test',
    name_shortform = 'T',
    timepoints = cumulative(rep(2.5, 4)),
    timeptname = paste0('T', 1:4),
    expectancies = c(1, 1, 0, 0),
    ondrug = list(pathA = c(1, 1, 0, 0))
  )
  path <- td$trialpaths[[1]]
  expect_equal(path$tpb[1], 2.5)
  expect_equal(path$tpb[2], 5)
  expect_equal(path$tpb[3], 0)
  expect_equal(path$tpb[4], 0)
})
