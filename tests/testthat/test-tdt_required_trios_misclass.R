test_that("tdt_required_trios_misclass runs correctly and produces valid outputs", {
  res <- tdt_required_trios_misclass(
    lambda_star = 7.8488,
    pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1,
    pi01 = 0.1
  )

  # Structure and names
  expect_type(res, "list")
  expect_true(all(c("gT_star", "gNT_star", "N_required", "lambda_star") %in% names(res)))

  # Value sanity checks
  expect_true(res$gT_star > res$gNT_star)
  expect_true(res$gT_star > 0 && res$gT_star < 1)
  expect_true(res$gNT_star > 0 && res$gNT_star < 1)
  expect_true(res$N_required > 0)

  # Misclassification adjustment: gT_star should drop slightly when pi01 > 0
  res_nomc <- tdt_required_trios_misclass(
    lambda_star = 7.8488,
    pd = 0.25, prev = 0.005,
    R1 = 2, R2 = 2,
    delta_prime = 1,
    pi01 = 0
  )

  expect_true(res$gT_star < res_nomc$gT_star)
})

