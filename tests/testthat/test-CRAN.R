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
test_that("params for reasonable penalties", {
  angle.vec <- c(0.1, 6, 6.1, 3, 3.1, 2.9)
  ##angle.vec <- c(1.1,1.2,1.3, 3, 3.1, 2.9)
  ##angle.vec <- c(1.1,1.3,3)
  (result <- geodesichange::geodesicFPOP_vec(angle.vec, 1))
  expect_equal(result$segments$param, c(6.1,3))
  (result <- geodesichange::geodesicFPOP_vec(angle.vec, 0.001))
  expect_equal(result$segments$param, angle.vec)
  (result <- geodesichange::geodesicFPOP_vec(angle.vec, 1000))
  expect_equal(nrow(result$segments),1)
})

