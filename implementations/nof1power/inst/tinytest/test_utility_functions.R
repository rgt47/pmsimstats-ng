library(tinytest)

defaults <- list(a = 1, b = 2, c = 3)
override <- list(b = 20, d = 40)

result <- nof1power::modify_list(defaults, override)
expect_equal(result$a, 1)
expect_equal(result$b, 20)
expect_equal(result$c, 3)
expect_equal(result$d, 40)

expect_equal(
  nof1power::modify_list(NULL, list(a = 1)),
  list(a = 1)
)

expect_equal(
  nof1power::modify_list(list(a = 1), NULL),
  list(a = 1)
)

result <- nof1power::modify_list(list(x = "old"), list(x = "new"))
expect_equal(result$x, "new")
expect_equal(length(result), 1)
