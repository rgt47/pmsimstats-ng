library(tinytest)

# cumulative forwards to cumsum
expect_equal(
  nof1power::cumulative(c(1, 2, 3, 4)),
  c(1, 3, 6, 10)
)

expect_equal(
  nof1power::cumulative(numeric(0)),
  numeric(0)
)

expect_equal(
  nof1power::cumulative(c(-1, 2, -3)),
  c(-1, 1, -2)
)

#---------------------------------------------------------------------
# mod_gompertz: origin-passing Gompertz
#---------------------------------------------------------------------

# y(0) = 0 exactly (origin-passing Hendrickson form).
expect_equal(
  nof1power::mod_gompertz(0, maxr = 10, disp = 1, rate = 0.5),
  0
)

# y -> maxr as t -> Inf.
expect_equal(
  nof1power::mod_gompertz(Inf, maxr = 10, disp = 1, rate = 0.5),
  10
)

# Vectorised call returns vector of same length.
expect_equal(
  length(nof1power::mod_gompertz(seq(0, 5, 1), maxr = 1,
                                 disp = 1, rate = 0.5)),
  6
)

# modgompertz wrapper is equivalent to mod_gompertz.
expect_equal(
  nof1power::modgompertz(2, max = 10, disp = 1, rate = 0.5),
  nof1power::mod_gompertz(2, maxr = 10, disp = 1, rate = 0.5)
)
