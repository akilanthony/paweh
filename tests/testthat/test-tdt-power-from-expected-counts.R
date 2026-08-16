test_that("TDT power function runs without error", {
  res <- tdt_power_from_expected_counts(ET = 140, ENT = 100, alpha = 0.05)
  expect_type(res, "list")
  expect_true("Power" %in% names(res))
})

