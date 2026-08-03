test_that("dominant_eigen returns largest eigenvalue", {
  m <- matrix(c(2, 0,
                0, 1), nrow = 2)

  expect_equal(dominant_eigen(m), 2)
})
