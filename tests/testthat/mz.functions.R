context("mz.* functions")

test_that("mz.ppm output", {
  result10 <- mz.ppm(135.5, 10)
  result100 <- mz.ppm(135.5, 100)
  expect_length(result10, 2)
  expect_that(result10, is_a("numeric"))
})
