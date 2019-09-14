context("2F1")

test_that("2F1 for scalar x", {
  obtained <- hypergeomPFQ(m = 25, a = c(1,2), b = c(3), x = 0.5)
  expected <- gsl::hyperg_2F1(1, 2, 3, 0.5)
  expect_equal(obtained, expected)
})

test_that("A value for 2F1", {
  obtained <- hypergeomPFQ(m = 10, a = c(1,2), b = c(3), x = c(0.2,0.5))
  expect_equal(obtained, 1.79412894456143)
})
