threshold_example <- list(
  qtl_var = 0.025,
  tau = 0.5,
  pd = 0.15,
  x_upper = 5,
  x_lower = 5
)

test_that("threshold parameters reproduce the Chapter 6 selection example", {
  out <- do.call(
    qtl_falconer_threshold_parameters,
    c(threshold_example, list(verbose = FALSE))
  )

  expect_s3_class(out, "qtl_falconer_threshold_parameters")
  expect_equal(out$upper_threshold, 1.644854, tolerance = 1e-6)
  expect_equal(out$lower_threshold, -1.644854, tolerance = 1e-6)
  expect_equal(unname(out$f_affected), c(0.0388, 0.0781, 0.0963), tolerance = 1e-3)
  expect_equal(unname(out$f_unaffected), c(0.0585, 0.0278, 0.0212), tolerance = 3e-3)
  expect_equal(out$prev_affected, 0.0501, tolerance = 1e-3)
  expect_equal(out$prev_unaffected, 0.0499, tolerance = 1e-3)
  expect_equal(unname(out$g_affected), c(0.5596, 0.3972, 0.0432), tolerance = 1e-3)
  expect_equal(unname(out$g_unaffected), c(0.8481, 0.1424, 0.0096), tolerance = 1e-3)
  expect_equal(sum(out$g_affected), 1, tolerance = 1e-12)
  expect_equal(sum(out$g_unaffected), 1, tolerance = 1e-12)
})

test_that("direct thresholds agree with their percentile equivalents", {
  percentile <- do.call(
    qtl_falconer_threshold_parameters,
    c(threshold_example, list(verbose = FALSE))
  )
  direct <- qtl_falconer_threshold_parameters(
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    threshold_mode = "direct",
    upper_threshold = percentile$upper_threshold,
    lower_threshold = percentile$lower_threshold,
    verbose = FALSE
  )

  expect_equal(direct$f_affected, percentile$f_affected, tolerance = 1e-12)
  expect_equal(direct$g_unaffected, percentile$g_unaffected, tolerance = 1e-12)
})

test_that("threshold-selected chi-square MSSN reproduces the published example", {
  out <- do.call(qtl_threshold_chisq_mssn, c(
    list(power = 0.80, alpha = 0.0001), threshold_example,
    list(k = 1, verbose = FALSE)
  ))

  expect_s3_class(out, "qtl_threshold_chisq_mssn")
  expect_equal(out$N_case, 126)
  expect_equal(out$N_control, 126)
  expect_equal(out$N_total, 252)
  expect_equal(out$lambda_star, 25.2520, tolerance = 1e-3)
  expect_gt(out$expected_population_screened_cases, out$N_case)
  expect_gt(out$expected_population_screened_controls, out$N_control)
})

test_that("threshold chi-square power increases with selected case sample size", {
  common <- c(list(alpha = 0.01), threshold_example, list(k = 1, verbose = FALSE))
  small <- do.call(qtl_threshold_chisq_power, c(list(N_case = 100), common))
  large <- do.call(qtl_threshold_chisq_power, c(list(N_case = 300), common))

  expect_gt(large$power, small$power)
  expect_gt(large$lambda, small$lambda)
})

test_that("threshold chi-square MSSN responds to alpha and QTL variance", {
  common <- c(list(power = 0.80), threshold_example, list(k = 1, verbose = FALSE))
  loose_alpha <- do.call(qtl_threshold_chisq_mssn, c(common, list(alpha = 0.05)))
  strict_alpha <- do.call(qtl_threshold_chisq_mssn, c(common, list(alpha = 0.0001)))
  small_qtl <- do.call(qtl_threshold_chisq_mssn, c(
    list(power = 0.80, alpha = 0.01, qtl_var = 0.01, tau = 0.5, pd = 0.15,
         x_upper = 5, x_lower = 5, k = 1, verbose = FALSE)
  ))
  large_qtl <- do.call(qtl_threshold_chisq_mssn, c(
    list(power = 0.80, alpha = 0.01, qtl_var = 0.05, tau = 0.5, pd = 0.15,
         x_upper = 5, x_lower = 5, k = 1, verbose = FALSE)
  ))

  expect_gt(strict_alpha$N_case, loose_alpha$N_case)
  expect_gt(small_qtl$N_case, large_qtl$N_case)
})

test_that("threshold chi-square has near-null power for a near-zero QTL variance", {
  out <- qtl_threshold_chisq_power(
    N_case = 1000, alpha = 0.05, qtl_var = 1e-8, tau = 0.5, pd = 0.15,
    x_upper = 5, x_lower = 5, k = 1, verbose = FALSE
  )

  expect_equal(out$power, 0.05, tolerance = 1e-3)
  expect_lt(out$lambda, 1e-3)
})

test_that("Falconer threshold functions validate thresholds and test inputs", {
  expect_error(
    qtl_falconer_threshold_parameters(0.025, 0.5, 0.15, x_upper = 0, x_lower = 5),
    "x_upper"
  )
  expect_error(
    qtl_falconer_threshold_parameters(0.025, 0.5, 0.15, x_upper = 5, x_lower = 100),
    "x_lower"
  )
  expect_error(
    qtl_falconer_threshold_parameters(
      0.025, 0.5, 0.15, threshold_mode = "direct",
      upper_threshold = 0, lower_threshold = 0
    ),
    "lower_threshold"
  )
  expect_error(
    qtl_threshold_chisq_power(0, 0.05, 0.025, 0.5, 0.15, 5, 5, verbose = FALSE),
    "N_case"
  )
  expect_error(
    qtl_threshold_chisq_power(100, 0, 0.025, 0.5, 0.15, 5, 5, verbose = FALSE),
    "alpha"
  )
  expect_error(
    qtl_threshold_chisq_mssn(0, 0.05, 0.025, 0.5, 0.15, 5, 5, verbose = FALSE),
    "power"
  )
  expect_error(
    qtl_threshold_chisq_mssn(0.8, 0.05, 0.025, 0.5, 0.15, 5, 5, k = 0, verbose = FALSE),
    "k"
  )
})
