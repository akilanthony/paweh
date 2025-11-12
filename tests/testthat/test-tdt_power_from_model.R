test_that("tdt_power_from_model runs correctly and returns expected structure", {
  res <- tdt_power_from_model(
    pd = 0.25, N = 10000, delta_prime = 1,
    prev = 0.05, R1 = 1, R2 = 1.1, alpha = 0.05
  )

  # Basic structure
  expect_type(res, "list")
  expect_true(all(c("Power", "Non-Centrality Parameter (lambda)") %in% names(res)))
  expect_true(res$`Power` > 0 && res$`Power` <= 1)

  # Check numeric stability and sensible output
  expect_true(res$`Non-Centrality Parameter (lambda)` > 0)
  expect_true(res$`Expected Transmissions (ET)` > res$`Expected Non-Transmissions (ENT)`)
})
