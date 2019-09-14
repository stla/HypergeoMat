context("1F0")

test_that("1F0 is det(I-X)^(-a)", {
  X <- toeplitz(c(3,2,1))/100
  obtained <- hypergeomPFQ(m = 15, a = 3, b = NULL, x = X)
  expected <- det(diag(3)-X)^(-3)
  expect_equal(obtained, expected)
})


