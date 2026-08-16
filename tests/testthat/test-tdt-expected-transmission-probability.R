test_that("tdt_expected_transmission_probability produces valid gT* output", {
  res <- tdt_expected_transmission_probability(
    pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1,
    pi01 = 0,
    verbose = FALSE
  )

  # Structure and key names
  expect_type(res, "list")
  expect_true(all(c("gT_star", "C", "D") %in% names(res)))

  # Reasonable numeric ranges
  expect_true(res$gT_star > 0)
  expect_true(res$gT_star < 1)
  expect_true(res$D > 0)
  expect_true(abs(res$C) < 1)

  # Compare misclassification effect (pi01 = 0.1)
  res_mis <- tdt_expected_transmission_probability(
    pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1,
    pi01 = 0.1,
    verbose = FALSE
  )

  # Ensure misclassification changes gT*
  expect_true(abs(res_mis$gT_star - res$gT_star) > 0)
})
