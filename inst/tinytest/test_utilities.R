expect_equal(cumulative(c(2.5, 2.5, 2.5, 2.5)), c(2.5, 5, 7.5, 10),
             info = 'cumulative converts equal intervals to running time')
expect_equal(cumulative(c(4, 4, 1, 1, 1, 1, 4, 4)),
             c(4, 8, 9, 10, 11, 12, 16, 20),
             info = 'cumulative handles variable intervals')

expect_equal(cumulative(c(3)), c(3),
             info = 'cumulative handles single element')

expect_equal(modgompertz(0, 10, 5, 0.42), 0,
             info = 'modgompertz returns 0 at t=0 (case 1)')
expect_equal(modgompertz(0, 6.50647, 5, 0.35), 0,
             info = 'modgompertz returns 0 at t=0 (case 2)')

expect_equal(modgompertz(200, 10, 5, 0.42), 10, tolerance = 0.01,
             info = 'modgompertz asymptotes at maxr')

{
  t_seq <- seq(0, 30, by = 0.5)
  vals <- modgompertz(t_seq, 10, 5, 0.42)
  expect_true(all(diff(vals) >= 0),
              info = 'modgompertz is monotonically non-decreasing')
}

{
  vals <- modgompertz(c(0, 5, 10), 10, 5, 0.42)
  expect_equal(length(vals), 3L,
               info = 'modgompertz vectorised: length matches input')
  expect_equal(vals[1], 0,
               info = 'modgompertz vectorised: first element is 0')
  expect_true(vals[2] < vals[3],
              info = 'modgompertz vectorised: increasing within curve')
}
