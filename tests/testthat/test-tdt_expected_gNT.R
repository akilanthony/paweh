test_that("tdt_expected_gNT produces valid gNT* output", {
  res <- tdt_expected_gNT(
    pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1,
    pi01 = 0,
    verbose = FALSE
  )

  # Structure and naming
  expect_type(res, "list")
  expect_true(all(c("gNT_star", "C", "D") %in% names(res)))

  # Basic numerical sanity checks
  expect_true(res$gNT_star > 0)
  expect_true(res$gNT_star < 1)
  expect_true(res$D > 0)
  expect_true(abs(res$C) < 1)

  # Misclassification sensitivity
  res_mis <- tdt_expected_gNT(
    pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1,
    pi01 = 0.1,
    verbose = FALSE
  )

  # Expect gNT* to increase slightly when misclassification > 0
  expect_true(abs(res_mis$gNT_star - res$gNT_star) > 0)
})

