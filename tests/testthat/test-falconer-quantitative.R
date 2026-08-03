falconer_example <- list(
  qtl_var = 0.025,
  tau = 0.5,
  pd = 0.15,
  alpha = 0.0001
)

test_that("Falconer parameters reproduce the Chapter 6 univariate example", {
  out <- do.call(
    falconer_parameters,
    c(falconer_example[c("qtl_var", "tau", "pd")], list(verbose = FALSE))
  )

  expect_s3_class(out, "falconer_parameters")
  expect_equal(out$a, 0.2279819, tolerance = 1e-6)
  expect_equal(out$delta, 0.1139909, tolerance = 1e-6)
  expect_equal(out$m, 0.1305196, tolerance = 1e-6)
  # The displayed equations give -0.09746; the published -0.0985 appears to
  # be a printed rounding/typographical discrepancy.
  expect_equal(unname(out$mu), c(-0.0974622, 0.2445105, 0.3585015), tolerance = 1e-6)
  expect_equal(unname(out$pi), c(0.7225, 0.255, 0.0225), tolerance = 1e-12)
  expect_equal(out$residual_variance, 0.975, tolerance = 1e-12)
  expect_equal(out$residual_sd, 0.9874209, tolerance = 1e-6)
})

test_that("Falconer parameters centre and standardize the trait", {
  out <- falconer_parameters(
    qtl_var = 0.10, tau = -0.35, pd = 0.27, verbose = FALSE
  )

  expect_equal(sum(out$pi), 1, tolerance = 1e-12)
  expect_equal(sum(out$pi * out$mu), 0, tolerance = 1e-12)
  expect_equal(out$weighted_mean, 0, tolerance = 1e-12)
  expect_equal(out$total_variance, 1, tolerance = 1e-12)
})

test_that("ANOVA power and MSSN reproduce the published example", {
  power_out <- do.call(qtl_anova_power, c(
    list(N = 996), falconer_example,
    list(count_method = "rounded", verbose = FALSE)
  ))
  mssn_out <- do.call(qtl_anova_mssn, c(
    list(power = 0.80), falconer_example,
    list(count_method = "rounded", multiple_of_three = TRUE, verbose = FALSE)
  ))

  expect_s3_class(power_out, "qtl_anova_power")
  expect_equal(power_out$df1, 2)
  expect_equal(power_out$df2, 993)
  expect_equal(unname(power_out$genotype_counts), c(720, 254, 22))
  expect_equal(power_out$lambda, 25.4894, tolerance = 1e-3)
  expect_equal(power_out$power, 0.80005, tolerance = 1e-4)

  expect_s3_class(mssn_out, "qtl_anova_mssn")
  expect_equal(mssn_out$N, 996)
  expect_equal(mssn_out$achieved_power, 0.80005, tolerance = 1e-4)
  expect_true(mssn_out$achieved_power >= mssn_out$target_power)
})

test_that("ANOVA power increases with sample size", {
  small <- qtl_anova_power(
    N = 600, alpha = 0.01, qtl_var = 0.025, tau = 0.5, pd = 0.15,
    verbose = FALSE
  )
  large <- qtl_anova_power(
    N = 1200, alpha = 0.01, qtl_var = 0.025, tau = 0.5, pd = 0.15,
    verbose = FALSE
  )

  expect_gt(large$power, small$power)
})

test_that("ANOVA MSSN responds to alpha and QTL variance", {
  common <- list(power = 0.80, tau = 0.5, pd = 0.15, verbose = FALSE)
  loose_alpha <- do.call(qtl_anova_mssn, c(common, list(alpha = 0.05, qtl_var = 0.025)))
  strict_alpha <- do.call(qtl_anova_mssn, c(common, list(alpha = 0.0001, qtl_var = 0.025)))
  small_qtl <- do.call(qtl_anova_mssn, c(common, list(alpha = 0.01, qtl_var = 0.01)))
  large_qtl <- do.call(qtl_anova_mssn, c(common, list(alpha = 0.01, qtl_var = 0.05)))

  expect_gt(strict_alpha$N, loose_alpha$N)
  expect_gt(small_qtl$N, large_qtl$N)
})

test_that("ANOVA has near-null power for a near-zero QTL variance", {
  out <- qtl_anova_power(
    N = 1000, alpha = 0.05, qtl_var = 1e-8, tau = 0.5, pd = 0.15,
    count_method = "expected", verbose = FALSE
  )

  expect_equal(out$power, 0.05, tolerance = 1e-4)
  expect_lt(out$lambda, 1e-3)
})

test_that("Falconer quantitative functions validate their inputs", {
  expect_error(falconer_parameters(qtl_var = 0, tau = 0.5, pd = 0.15), "qtl_var")
  expect_error(falconer_parameters(qtl_var = 1, tau = 0.5, pd = 0.15), "qtl_var")
  expect_error(falconer_parameters(qtl_var = 0.025, tau = 0.5, pd = 0), "pd")
  expect_error(qtl_anova_power(3, 0.05, 0.025, 0.5, 0.15, verbose = FALSE), "N")
  expect_error(qtl_anova_power(10.5, 0.05, 0.025, 0.5, 0.15, verbose = FALSE), "integer")
  expect_error(qtl_anova_power(100, 1, 0.025, 0.5, 0.15, verbose = FALSE), "alpha")
  expect_error(qtl_anova_mssn(1, 0.05, 0.025, 0.5, 0.15, verbose = FALSE), "power")
})
