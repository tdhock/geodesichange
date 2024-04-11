library(testthat)
test_that("param is 1.5", {
  (result <- geodesichange::geodesicFPOP_vec(c(1,1.5,3), Inf))
  expect_equal(result$segments$param, 1.5)
})
test_that("param is 0.1", {
  (result <- geodesichange::geodesicFPOP_vec(c(0.1,1,6), Inf))
  expect_equal(result$segments$param, 0.1)
})
test_that("param is 6.1", {
  (result <- geodesichange::geodesicFPOP_vec(c(6.1,1,6), Inf))
  expect_equal(result$segments$param, 6.1)
})
test_that("param is same as data", {
  angle.vec <- c(0.1, 6, 6.1, 6.2, 3, 3.1, 2.9)
  (result <- geodesichange::geodesicFPOP_vec(angle.vec, 0))
  expect_equal(result$segments$param, angle.vec)
})
test_that("two params for reasonable penalty", {
  angle.vec <- c(0.1, 6, 6.1, 3, 3.1, 2.9)
  angle.vec <- c(1.1,1.2,1.3, 3, 3.1, 2.9)
  angle.vec <- c(1.1,1.2,3)
  (result <- geodesichange::geodesicFPOP_vec(angle.vec, 0.1))
  result$segments
  expect_equal(rev(result$segments$param), angle.vec)
})

