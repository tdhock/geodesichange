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

iloop.vec <- c(
  4.9823357159891, 2.23772298397808, 1.41985124082998, 2.00970131165917, 
  3.46980108993953, 3.32069086176044, 0.255812812453875, 2.08157188638142, 
  3.37763278103609, 1.95651413585301, 6.100145143647, 2.18852682570223, 
  0.43141096152997, 1.75116363420413, 2.79256556789621, 5.45589383326984, 
  6.11380144189501, 2.09255995819222, 2.97204096127947, 4.42737237250248, 
  2.6466996877735, 1.59789500382623, 0.984757668637457, 2.3626926741305, 
  6.05497231303033, 2.43138205917179, 2.05383986573757, 0.369984111202642, 
  4.07747369163124)
test_that("params for reasonable penalties", {
  iloop.vec <- c(
    4.9823357159891, 2.23772298397808, 1.41985124082998, 2.00970131165917, 
    3.46980108993953, 3.32069086176044, 0.255812812453875, 2.08157188638142, 
    3.37763278103609, 1.95651413585301, 6.100145143647, 2.18852682570223, 
    0.43141096152997, 1.75116363420413, 2.79256556789621, 5.45589383326984, 
    6.11380144189501, 2.09255995819222, 2.97204096127947, 4.42737237250248, 
    2.6466996877735, 1.59789500382623, 0.984757668637457, 2.3626926741305, 
    6.05497231303033, 2.43138205917179, 2.05383986573757, 0.369984111202642, 
    4.07747369163124)
  (result <- geodesichange::geodesicFPOP_vec(iloop.vec, 1))
  expect_is(result$segments$param, "numeric")
  (result <- geodesichange::geodesicFPOP_vec(iloop.vec, 0))
  expect_lt(sum(abs(result$segments$param-iloop.vec)), 1e-3)
  (result <- geodesichange::geodesicFPOP_vec(iloop.vec, Inf))
  expect_equal(length(result$segments$param), 1)
})
