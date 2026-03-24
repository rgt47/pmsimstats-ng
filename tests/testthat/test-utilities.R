test_that('cumulative converts intervals to running time', {
  expect_equal(cumulative(c(2.5, 2.5, 2.5, 2.5)), c(2.5, 5, 7.5, 10))
  expect_equal(cumulative(c(4, 4, 1, 1, 1, 1, 4, 4)), c(4, 8, 9, 10, 11, 12, 16, 20))
})

test_that('cumulative handles single element', {
  expect_equal(cumulative(c(3)), c(3))
})

test_that('modgompertz returns 0 at t=0', {
  expect_equal(modgompertz(0, 10, 5, 0.42), 0)
  expect_equal(modgompertz(0, 6.50647, 5, 0.35), 0)
})

test_that('modgompertz asymptotes at maxr', {
  expect_equal(modgompertz(200, 10, 5, 0.42), 10, tolerance = 0.01)
})

test_that('modgompertz is monotonically increasing for positive maxr', {
  t_seq <- seq(0, 30, by = 0.5)
  vals <- modgompertz(t_seq, 10, 5, 0.42)
  expect_true(all(diff(vals) >= 0))
})

test_that('modgompertz is vectorized over t', {
  vals <- modgompertz(c(0, 5, 10), 10, 5, 0.42)
  expect_length(vals, 3)
  expect_equal(vals[1], 0)
  expect_true(vals[2] < vals[3])
})
