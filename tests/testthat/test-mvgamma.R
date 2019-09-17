context("Multivariate Gamma function")

test_that("Comparison with gamma", {
  expect_equal(mvgamma(-1.5,1), gamma(-1.5))
  expect_equal(mvgamma(-2.5,1), gamma(-2.5))
})
