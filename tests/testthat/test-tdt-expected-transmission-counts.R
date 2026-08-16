test_that("tdt_expected_transmission_counts computes valid ET* and ENT*", {
  res <- tdt_expected_transmission_counts(
    N_star = 1000, pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1, pi = 1,
    verbose = FALSE
  )

  # Structure and expected keys
  expect_type(res, "list")
  expect_true(all(c("ET_star", "ENT_star", "C", "D") %in% names(res)))

  # Numeric sanity
  expect_true(res$ET_star > 0)
  expect_true(res$ENT_star > 0)
  expect_true(res$ET_star > res$ENT_star)
  expect_true(abs(res$C) < 1)
  expect_true(res$D > 0)

  # Heterogeneity sensitivity
  res_het <- tdt_expected_transmission_counts(
    N_star = 1000, pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1, pi = 0.7,
    verbose = FALSE
  )

  # Check that ET* changes under heterogeneity
  expect_true(abs(res_het$ET_star - res$ET_star) > 0)
  expect_true(abs(res_het$ENT_star - res$ENT_star) > 0)
})

