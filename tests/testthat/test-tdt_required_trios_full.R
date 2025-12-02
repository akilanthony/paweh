test_that("tdt_required_trios_full returns expected structure and class", {
  out <- tdt_required_trios_full(
    target_power  = 0.80,
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

  expect_s3_class(out, "tdt_required_trios_full")
  expect_type(out$N, "list")
  expect_true(all(c("no_error", "misclassification", "heterogeneity") %in% names(out$N)))
  expect_type(out$percent_increase, "list")
})

test_that("tdt_required_trios_full returns expected structure and class", {
  out <- tdt_required_trios_full(
    target_power  = 0.80,
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

  expect_s3_class(out, "tdt_required_trios_full")
  expect_type(out$N, "list")
  expect_true(all(c("no_error", "misclassification", "heterogeneity") %in% names(out$N)))
  expect_type(out$percent_increase, "list")
})


test_that("required trios increase under misclassification and heterogeneity", {
  out <- tdt_required_trios_full(
    target_power  = 0.80,
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

  expect_lt(out$N$no_error,         out$N$misclassification)
  expect_lt(out$N$no_error,         out$N$heterogeneity)
  expect_gte(out$percent_increase$misclassification, 0)
  expect_gte(out$percent_increase$heterogeneity,    0)
})

