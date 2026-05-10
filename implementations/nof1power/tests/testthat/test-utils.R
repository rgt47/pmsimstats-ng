library(tinytest)
library(nof1power)

# create_ar1_corr generates valid AR(1) correlation matrix
rho <- 0.5
n <- 4
corr_matrix <- create_ar1_corr(n, rho)

expect_equal(dim(corr_matrix), c(n, n))
expect_true(all.equal(corr_matrix, t(corr_matrix)))
expect_true(all(diag(corr_matrix) == 1))
expect_equal(corr_matrix[1, 2], rho)
expect_equal(corr_matrix[1, 3], rho^2)
