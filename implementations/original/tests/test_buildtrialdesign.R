{
  td <- buildtrialdesign(
    name_longform = 'open label',
    name_shortform = 'OL',
    timepoints = cumulative(rep(2.5, 8)),
    timeptnames = paste0('OL', 1:8),
    expectancies = rep(1, 8),
    ondrug = list(pathA = rep(1, 8))
  )
  expect_equal(td$metadata$name_shortform, 'OL',
               info = 'OL design metadata shortform')
  expect_equal(length(td$trialpaths), 1L,
               info = 'OL design has one path')
  path <- td$trialpaths[[1]]
  expect_equal(nrow(path), 8L,
               info = 'OL path has 8 timepoints')
  expect_true(all(path$tod > 0),
              info = 'OL path is always on drug')
  expect_true(all(path$tsd == 0),
              info = 'OL path has no time since discontinuation')
}

{
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
  expect_equal(length(td$trialpaths), 2L,
               info = 'CO design has two paths')
  pa <- td$trialpaths[[1]]
  pb <- td$trialpaths[[2]]
  expect_true(all(pa$tod[1:4] > 0),
              info = 'CO path A: first four timepoints on drug')
  expect_true(all(pa$tod[5:8] == 0),
              info = 'CO path A: last four timepoints off drug')
  expect_true(all(pb$tod[1:4] == 0),
              info = 'CO path B: first four timepoints off drug')
  expect_true(all(pb$tod[5:8] > 0),
              info = 'CO path B: last four timepoints on drug')

  expect_equal(pa$tsd[5:8], c(2.5, 5, 7.5, 10),
               info = 'CO path A: tsd accumulates during washout')
}

{
  td <- buildtrialdesign(
    name_longform = 'test',
    name_shortform = 'T',
    timepoints = cumulative(rep(2.5, 4)),
    timeptnames = paste0('T', 1:4),
    expectancies = c(1, 1, 0, 0),
    ondrug = list(pathA = c(1, 1, 0, 0))
  )
  path <- td$trialpaths[[1]]
  expect_equal(path$tpb[1], 2.5,
               info = 'tpb accumulates at t1 with positive expectancy')
  expect_equal(path$tpb[2], 5,
               info = 'tpb accumulates at t2 with positive expectancy')
  expect_equal(path$tpb[3], 0,
               info = 'tpb is zero at t3 with zero expectancy')
  expect_equal(path$tpb[4], 0,
               info = 'tpb is zero at t4 with zero expectancy')
}
