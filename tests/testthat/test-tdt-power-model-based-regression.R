test_that("tdt_power returns expected structure and class", {
  out <- tdt_power(
    N             = 500,
    pd            = 0.30,
    prev          = 0.05,
    R1            = 1.5,
    R2            = 2.25,
    alpha         = 0.05,
    delta_prime   = 1,
    misclass_rate = 0.05,
    heter_rate    = 0.10,
    verbose       = FALSE
  )

  expect_s3_class(out, "tdt_power")
  expect_type(out$lambda, "list")
  expect_type(out$power, "list")
  expect_true(all(c("no_error", "misclassification", "heterogeneity") %in% names(out$power)))
})


test_that("tdt_power enforces misclass_rate and heter_rate bounds", {
  expect_error(
    tdt_power(
      N             = 500,
      pd            = 0.3,
      prev          = 0.05,
      R1            = 1.5,
      R2            = 2.25,
      misclass_rate = -0.01,
      heter_rate    = 0.01,
      verbose       = FALSE
    ),
    regexp = "misclass_rate must be in \\[0, 1\\)"
  )

  expect_error(
    tdt_power(
      N             = 500,
      pd            = 0.3,
      prev          = 0.05,
      R1            = 1.5,
      R2            = 2.25,
      misclass_rate = 0.01,
      heter_rate    = 1,
      verbose       = FALSE
    ),
    regexp = "heter_rate must be in \\[0, 1\\)"
  )
})


test_that("misclassification and heterogeneity reduce power for fixed N", {
  base <- tdt_power(
    N             = 500,
    pd            = 0.30,
    prev          = 0.05,
    R1            = 1.5,
    R2            = 2.25,
    alpha         = 0.05,
    delta_prime   = 1,
    misclass_rate = 0,
    heter_rate    = 0,
    verbose       = FALSE
  )

  err  <- tdt_power(
    N             = 500,
    pd            = 0.30,
    prev          = 0.05,
    R1            = 1.5,
    R2            = 2.25,
    alpha         = 0.05,
    delta_prime   = 1,
    misclass_rate = 0.10,
    heter_rate    = 0.20,
    verbose       = FALSE
  )

  expect_lt(err$power$misclassification, base$power$no_error)
  expect_lt(err$power$heterogeneity,    base$power$no_error)
})
