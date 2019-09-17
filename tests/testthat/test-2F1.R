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

test_that("Gauss formula", {
  a <- 1
  b <- 2
  c <- 9
  o1 <- mvgamma(c,3)*mvgamma(c-a-b,3)/mvgamma(c-a,3)/mvgamma(c-b,3)
  o2 <- hypergeomPFQ(300, c(a,b), c, c(1,1,1))
  expect_equal(o1, o2)
})

test_that("Reflection formula", {
  a <- 1.1
  b <- 2.2
  c <- 8.1
  X <- diag(3)*0.99
  o1 <- hypergeomPFQ(100, c(a,b), c, X)
  o2 <- mvgamma(c,3)*mvgamma(c-a-b,3)/mvgamma(c-a,3)/mvgamma(c-b,3) *
    hypergeomPFQ(100, c(a,b), a+b-c+(3+1)/2, diag(3)-X) +
    mvgamma(c,3)*mvgamma(a+b-c,3)/mvgamma(a,3)/mvgamma(b,3) *
    hypergeomPFQ(100, c(c-a,c-b), (3+1)/2+c-a-b, diag(3)-X) *
    det(diag(3)-X)^(c-a-b)
  expect_equal(o1, o2, tolerance = 1e-4)
})
