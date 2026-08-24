mv_example <- list(
  qtl_var = c(0.01, 0.005),
  tau = c(0, 0.5),
  pd = 0.25,
  cor_matrix = matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
)

mv_threshold_args <- c(
  mv_example,
  list(x_upper = c(10, 10), x_lower = c(15, 15))
)

test_that("multivariate model validates vector and correlation inputs", {
  base <- c(list(N = 100, alpha = 0.05, test = "pillai", verbose = FALSE), mv_example)
  expect_error(do.call(qtl_multivariate_power_full, within(base, tau <- 0)), "tau")
  expect_error(do.call(qtl_multivariate_power_full, within(base, qtl_var <- c(0.01, 1))), "qtl_var")
  expect_error(do.call(qtl_multivariate_power_full, within(base, cor_matrix <- diag(3))), "p x p")

  nonsymmetric <- matrix(c(1, 0.1, 0.2, 1), 2)
  expect_error(do.call(qtl_multivariate_power_full, within(base, cor_matrix <- nonsymmetric)), "symmetric")
  bad_diag <- matrix(c(0.9, 0.1, 0.1, 1), 2)
  expect_error(do.call(qtl_multivariate_power_full, within(base, cor_matrix <- bad_diag)), "unit diagonal")
  invalid_cor <- matrix(c(1, 1.1, 1.1, 1), 2)
  expect_error(do.call(qtl_multivariate_power_full, within(base, cor_matrix <- invalid_cor)), "correlations")
  singular_cor <- matrix(1, 2, 2)
  expect_error(do.call(qtl_multivariate_power_full, within(base, cor_matrix <- singular_cor)), "positive definite")
  expect_error(do.call(qtl_multivariate_power_full, within(base, N <- 100.5)), "integer")
})

test_that("Chapter 6.2 Falconer parameters and matrix orientation are reproduced", {
  out <- do.call(
    qtl_multivariate_power_full,
    c(list(N = 4514, alpha = 5e-8, test = "pillai", verbose = FALSE), mv_example)
  )
  model <- out$falconer
  expect_identical(dim(model$mean_matrix), c(2L, 3L))
  expect_identical(model$mean_matrix_orientation, "rows = phenotypes; columns = genotypes")
  expect_equal(model$parameters$a, c(0.16329932, 0.08972354), tolerance = 1e-7)
  expect_equal(model$parameters$delta, c(0, 0.04486177), tolerance = 1e-7)
  expect_equal(model$parameters$m, c(0.08164966, 0.02803861), tolerance = 1e-7)
  expect_equal(
    unname(model$mean_matrix),
    matrix(c(
      -0.08164966, 0.08164966, 0.24494897,
      -0.06168494, 0.07290038, 0.11776215
    ), 2, byrow = TRUE),
    tolerance = 1e-7
  )
  expect_equal(unname(model$genotype_frequencies), c(0.5625, 0.375, 0.0625))
})

test_that("variance decomposition and residual covariance are correct", {
  model <- do.call(paweh:::.falconer_mv_parameters, mv_example)
  expect_equal(model$parameters$additive_variance, c(0.01, 0.004716981), tolerance = 1e-8)
  expect_equal(model$parameters$dominance_variance, c(0, 0.000283018867924528), tolerance = 1e-12)
  expect_equal(
    model$parameters$additive_variance + model$parameters$dominance_variance,
    mv_example$qtl_var,
    tolerance = 1e-12
  )
  expect_equal(unname(diag(model$residual_covariance_matrix)), c(0.99, 0.995))
  expect_equal(model$residual_covariance_matrix[1, 2], 0.15 * sqrt(0.99 * 0.995), tolerance = 1e-12)
})

test_that("joint thresholds and deterministic MVN penetrances reproduce the example", {
  out <- do.call(
    qtl_multivariate_mssn_full,
    c(list(power = 0.95, alpha = 5e-8, test = "threshold_chisq", k = 1, verbose = FALSE), mv_threshold_args)
  )
  th <- out$thresholds
  expect_equal(th$upper_threshold, rep(stats::qnorm(0.9), 2))
  expect_equal(th$lower_threshold, rep(stats::qnorm(0.15), 2))
  expect_equal(
    unname(th$penetrances$affected),
    c(0.011904559517903, 0.0190035838467114, 0.0257671363375846),
    tolerance = 1e-10
  )
  expect_equal(
    unname(th$penetrances$unaffected),
    c(0.0377034348674891, 0.0248592441194994, 0.0181630964441482),
    tolerance = 1e-10
  )
  expect_equal(unname(th$prevalences), c(0.0154331046924363, 0.0316655921855342), tolerance = 1e-10)
  expect_equal(sum(th$frequencies$case), 1, tolerance = 1e-12)
  expect_equal(sum(th$frequencies$control), 1, tolerance = 1e-12)
  expect_true(all(vapply(th$integration$affected, `[[`, character(1), "algorithm") == "Miwa"))

  # The independently validated modern MVN values are the regression targets.
  # The historical C values are descriptive only because its Riemann-sum
  # integration gives the coarser printed approximations below.
  expect_lt(max(abs(unname(th$penetrances$affected) - c(0.011836, 0.018901, 0.025635))), 5e-4)
  expect_lt(max(abs(unname(th$penetrances$unaffected) - c(0.037275, 0.024555, 0.017930))), 5e-4)
})

test_that("independent direct mvtnorm calculation confirms bivariate penetrances", {
  model <- do.call(paweh:::.falconer_mv_parameters, mv_example)
  direct_affected <- vapply(1:3, function(j) {
    as.numeric(mvtnorm::pmvnorm(
      lower = rep(stats::qnorm(0.9), 2), upper = rep(Inf, 2),
      mean = model$mean_matrix[, j], sigma = model$residual_covariance_matrix,
      algorithm = mvtnorm::Miwa(steps = 128L)
    ))
  }, numeric(1))
  out <- do.call(
    qtl_multivariate_power_full,
    c(list(N_case = 100, alpha = 0.05, test = "threshold_chisq", verbose = FALSE), mv_threshold_args)
  )
  expect_equal(unname(out$thresholds$penetrances$affected), direct_affected, tolerance = 1e-12)
})

test_that("threshold chi-square MSSN, NCP, counts, and sparse diagnostics are auditable", {
  out <- do.call(
    qtl_multivariate_mssn_full,
    c(list(power = 0.95, alpha = 5e-8, test = "threshold_chisq", k = 1, verbose = FALSE), mv_threshold_args)
  )
  expect_s3_class(out, "qtl_multivariate_mssn_full")
  expect_equal(out$target_noncentrality_parameter, 54.26847, tolerance = 1e-5)
  expect_lt(abs(out$historical_fractional_cases - 448.8045), 1e-4)
  expect_identical(out$N_case, 449)
  expect_identical(out$N_control, 449)
  expect_identical(out$N_total, 898)
  expect_gte(out$achieved_power, 0.95)
  expect_equal(
    unname(out$expected_counts),
    matrix(c(194.8179, 207.3289, 46.8532, 300.7199, 132.1837, 16.0964), 2, byrow = TRUE),
    tolerance = 1e-4
  )
  expect_identical(out$cells_below_one, 0L)
  expect_equal(out$percent_cells_below_one, 0)
  expect_gt(out$expected_population_screened_cases, out$N_case)
  expect_gt(out$expected_population_screened_controls, out$N_control)
})

test_that("Pillai derivation reproduces Chapter 6.2 NCP and integer MSSN", {
  historical <- do.call(
    paweh:::.falconer_mv_pillai_components,
    list(N = 4513.6, model = do.call(paweh:::.falconer_mv_parameters, mv_example))
  )
  expect_equal(historical$lambda, 59.88577, tolerance = 1e-5)
  expect_identical(dim(historical$phi_star), c(2L, 2L))
  expect_equal(historical$eigenvalues, c(0.0131266207741616, 0.000223917176545547), tolerance = 1e-11)

  out <- do.call(
    qtl_multivariate_mssn_full,
    c(list(power = 0.95, alpha = 5e-8, test = "pillai", verbose = FALSE), mv_example)
  )
  expect_identical(out$N, 4514L)
  expect_equal(out$historical_fractional_mssn, 4513.5796, tolerance = 1e-3)
  expect_equal(out$noncentrality_parameter, 59.89108, tolerance = 1e-5)
  expect_identical(out$numerator_df, 4L)
  expect_equal(out$denominator_df, 9022)
  expect_equal(out$achieved_power, 0.9500367, tolerance = 1e-7)
})

test_that("power and MSSN move in expected directions", {
  p_small <- do.call(qtl_multivariate_power_full, c(
    list(N = 1000, alpha = 5e-8, test = "pillai", verbose = FALSE), mv_example
  ))
  p_large <- do.call(qtl_multivariate_power_full, c(
    list(N = 5000, alpha = 5e-8, test = "pillai", verbose = FALSE), mv_example
  ))
  expect_gt(p_large$power, p_small$power)

  m_loose <- do.call(qtl_multivariate_mssn_full, c(
    list(power = 0.8, alpha = 1e-4, test = "pillai", verbose = FALSE), mv_example
  ))
  m_strict <- do.call(qtl_multivariate_mssn_full, c(
    list(power = 0.8, alpha = 5e-8, test = "pillai", verbose = FALSE), mv_example
  ))
  expect_gt(m_strict$N, m_loose$N)

  larger_effect <- mv_example
  larger_effect$qtl_var <- 2 * larger_effect$qtl_var
  m_effect <- do.call(qtl_multivariate_mssn_full, c(
    list(power = 0.8, alpha = 5e-8, test = "pillai", verbose = FALSE), larger_effect
  ))
  expect_lt(m_effect$N, m_strict$N)
})

test_that("phenotype correlation changes threshold and Pillai results", {
  independent <- mv_example
  independent$cor_matrix <- diag(2)
  threshold_cor <- do.call(qtl_multivariate_power_full, c(
    list(N_case = 300, alpha = 0.05, test = "threshold_chisq", verbose = FALSE), mv_threshold_args
  ))
  threshold_ind <- do.call(qtl_multivariate_power_full, c(
    list(N_case = 300, alpha = 0.05, test = "threshold_chisq", x_upper = c(10, 10), x_lower = c(15, 15), verbose = FALSE), independent
  ))
  expect_false(isTRUE(all.equal(threshold_cor$power, threshold_ind$power)))

  pillai_cor <- do.call(qtl_multivariate_power_full, c(
    list(N = 2000, alpha = 0.05, test = "pillai", verbose = FALSE), mv_example
  ))
  pillai_ind <- do.call(qtl_multivariate_power_full, c(
    list(N = 2000, alpha = 0.05, test = "pillai", verbose = FALSE), independent
  ))
  expect_false(isTRUE(all.equal(pillai_cor$power, pillai_ind$power)))
})

test_that("p = 1 threshold behavior agrees with the existing single-trait backend", {
  old <- qtl_falconer_threshold_parameters(
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    x_upper = 5, x_lower = 5, verbose = FALSE
  )
  new <- qtl_multivariate_power_full(
    N_case = 126, alpha = 0.0001,
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    cor_matrix = matrix(1), test = "threshold_chisq",
    x_upper = 5, x_lower = 5, verbose = FALSE
  )
  old_power <- qtl_threshold_chisq_power(
    N_case = 126, alpha = 0.0001,
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    x_upper = 5, x_lower = 5, verbose = FALSE
  )
  expect_equal(new$thresholds$penetrances$affected, old$f_affected, tolerance = 1e-14)
  expect_equal(new$thresholds$penetrances$unaffected, old$f_unaffected, tolerance = 1e-14)
  expect_equal(new$power, old_power$power, tolerance = 1e-14)
})

test_that("threshold inputs are validated", {
  common <- c(
    list(N_case = 100, alpha = 0.05, test = "threshold_chisq", verbose = FALSE),
    mv_example
  )
  expect_error(do.call(qtl_multivariate_power_full, c(common, list(x_upper = 10, x_lower = c(15, 15)))), "x_upper")
  expect_error(do.call(qtl_multivariate_power_full, c(common, list(x_upper = c(0, 10), x_lower = c(15, 15)))), "x_upper")
  expect_error(do.call(qtl_multivariate_power_full, c(common, list(x_upper = c(60, 60), x_lower = c(60, 60)))), "lower threshold")
})

test_that("verbose is silent when false and printing does not change results", {
  args <- c(
    list(N = 4514, alpha = 5e-8, test = "pillai"),
    mv_example
  )
  expect_silent(quiet <- do.call(qtl_multivariate_power_full, c(args, list(verbose = FALSE))))
  printed <- suppressMessages(do.call(qtl_multivariate_power_full, c(args, list(verbose = TRUE))))
  expect_equal(unclass(quiet), unclass(printed))

  threshold_args <- c(
    list(power = 0.95, alpha = 5e-8, test = "threshold_chisq", k = 1),
    mv_threshold_args
  )
  expect_silent(quiet_t <- do.call(qtl_multivariate_mssn_full, c(threshold_args, list(verbose = FALSE))))
  printed_t <- suppressMessages(do.call(qtl_multivariate_mssn_full, c(threshold_args, list(verbose = TRUE))))
  expect_equal(unclass(quiet_t), unclass(printed_t))
})
