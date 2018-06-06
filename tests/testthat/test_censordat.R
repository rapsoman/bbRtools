context("Config")

library(bbRtools)

test_that("Symmetric censoring works correctly",{
  dat = 1:1000
  quant = 0.99
  out_dat = bbRtools::censor_dat(dat, quant=quant, symmetric=T)
  expect_equal(max(out_dat), as.numeric(quantile(dat, 0.995)))
  expect_equal(min(out_dat),as.numeric(quantile(dat, 0.005)))
})