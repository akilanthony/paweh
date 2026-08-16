test_that("tdt_mssn_from_model computes reasonable output and structure", {
  res <- tdt_mssn_from_model(
    power = 0.8, alpha = 0.05, df = 1,
    pd = 0.25, prev = 0.005, R1 = 2, R2 = 2,
    delta_prime = 1, pi = 1
  )

  # Ensure correct object structure
  expect_type(res, "list")
  expect_true(all(c(
    "Non-Centrality Parameter (lambda_star)",
    "Expected Transmission (gT_star)",
    "Expected Non-Transmission (gNT_star)",
    "Required Number of Trios (N_star)"
  ) %in% names(res)))

  # Check logical and numeric constraints
  expect_true(res$`Non-Centrality Parameter (lambda_star)` > 0)
  expect_true(res$`Expected Transmission (gT_star)` > res$`Expected Non-Transmission (gNT_star)`)
  expect_true(res$`Required Number of Trios (N_star)` > 0)

  # Check reproducibility (rounded values should be stable)
  expect_equal(round(res$`Expected Transmission (gT_star)`, 5), 0.26087, tolerance = 1e-3)
})

