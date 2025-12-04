test_that("tdt_plot_power returns a ggplot/patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  p <- tdt_plot_power(
    N    = 200,
    pd   = 0.30,
    prev = 0.05,
    R1   = 1.5,
    R2   = 2.25,
    alpha = 0.05,
    delta_prime   = 1,
    misclass_seq  = c(0, 0.05, 0.10),
    heter_seq     = c(0, 0.10, 0.20),
    heter_fixed   = 0,
    misclass_fixed = 0
  )

  # patchwork objects usually inherit from "gg"/"ggplot"
  expect_s3_class(p, c("gg", "ggplot"), exact = FALSE)
})

test_that("tdt_plot_sample_size returns a ggplot/patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  p <- tdt_plot_sample_size(
    target_power = 0.80,
    pd   = 0.30,
    prev = 0.05,
    R1   = 1.5,
    R2   = 2.25,
    alpha = 0.05,
    delta_prime   = 1,
    misclass_seq  = c(0, 0.05, 0.10),
    heter_seq     = c(0, 0.10, 0.20),
    heter_fixed   = 0,
    misclass_fixed = 0
  )

  expect_s3_class(p, c("gg", "ggplot"), exact = FALSE)
})

test_that("tdt_plot_* handle obviously bad input", {
  expect_error(
    tdt_plot_power(
      N    = 200,
      pd   = 0.30,
      prev = 0.05,
      R1   = 1.5,
      R2   = 2.25,
      alpha = 0.05,
      delta_prime   = 1,
      misclass_seq  = c(-0.01, 0.05),  # invalid misclassification
      heter_seq     = c(0, 0.10),
      heter_fixed   = 0,
      misclass_fixed = 0
    ),
    regexp = "misclass|misclassification",
    ignore.case = TRUE
  )
})
