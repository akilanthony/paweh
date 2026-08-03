plot_falconer_fixed <- list(
  alpha = 0.0001,
  qtl_var = 0.025,
  tau = 0.5,
  pd = 0.15
)

test_that("Falconer power wrappers return ggplot objects", {
  anova <- do.call(plot_qtl_anova_power, c(
    list(x_var = "N", x_values = c(900, 996, 1100)),
    plot_falconer_fixed
  ))
  threshold <- do.call(plot_qtl_threshold_chisq_power, c(
    list(x_var = "N_case", x_values = c(100, 126, 150)),
    plot_falconer_fixed,
    list(x_upper = 5, x_lower = 5, k = 1)
  ))

  expect_s3_class(anova, "ggplot")
  expect_s3_class(threshold, "ggplot")
  expect_equal(anova$labels$x, "Total sample size")
  expect_equal(threshold$labels$x, "Number of selected cases")
})

test_that("Falconer MSSN wrappers return ggplot objects", {
  anova <- plot_qtl_anova_mssn(
    x_var = "qtl_var", x_values = c(0.02, 0.025, 0.03),
    power = 0.8, alpha = 0.0001, tau = 0.5, pd = 0.15
  )
  threshold <- plot_qtl_threshold_chisq_mssn(
    x_var = "qtl_var", x_values = c(0.02, 0.025, 0.03),
    sample_size = "case", power = 0.8, alpha = 0.0001,
    tau = 0.5, pd = 0.15, x_upper = 5, x_lower = 5, k = 1
  )

  expect_s3_class(anova, "ggplot")
  expect_s3_class(threshold, "ggplot")
  expect_equal(anova$labels$y, "Required total sample size")
  expect_equal(threshold$labels$y, "Required selected cases")
})

test_that("Falconer power plotting data match direct backend calls", {
  N_values <- c(900, 996, 1100)
  anova_dat <- do.call(plot_qtl_anova_power, c(
    list(x_var = "N", x_values = N_values, return_data = TRUE),
    plot_falconer_fixed
  ))
  expected_anova <- vapply(N_values, function(N) {
    do.call(qtl_anova_power, c(
      list(N = N, verbose = FALSE), plot_falconer_fixed
    ))$power
  }, numeric(1))

  case_values <- c(100, 126, 150)
  threshold_dat <- do.call(plot_qtl_threshold_chisq_power, c(
    list(x_var = "N_case", x_values = case_values, return_data = TRUE),
    plot_falconer_fixed,
    list(x_upper = 5, x_lower = 5, k = 1)
  ))
  expected_threshold <- vapply(case_values, function(N_case) {
    do.call(qtl_threshold_chisq_power, c(
      list(N_case = N_case, verbose = FALSE), plot_falconer_fixed,
      list(x_upper = 5, x_lower = 5, k = 1)
    ))$power
  }, numeric(1))

  expect_s3_class(anova_dat, "data.frame")
  expect_equal(names(anova_dat), c("N", "power"))
  expect_equal(anova_dat$power, expected_anova)
  expect_s3_class(threshold_dat, "data.frame")
  expect_equal(names(threshold_dat), c("N_case", "power"))
  expect_equal(threshold_dat$power, expected_threshold)
})

test_that("Falconer MSSN plotting data match direct backend calls", {
  qtl_values <- c(0.02, 0.025, 0.03)
  anova_dat <- plot_qtl_anova_mssn(
    x_var = "qtl_var", x_values = qtl_values,
    power = 0.8, alpha = 0.0001, tau = 0.5, pd = 0.15,
    return_data = TRUE
  )
  threshold_dat <- plot_qtl_threshold_chisq_mssn(
    x_var = "qtl_var", x_values = qtl_values,
    sample_size = "total", power = 0.8, alpha = 0.0001,
    tau = 0.5, pd = 0.15, x_upper = 5, x_lower = 5, k = 1,
    return_data = TRUE
  )

  expected_anova <- vapply(qtl_values, function(qtl_var) {
    qtl_anova_mssn(
      power = 0.8, alpha = 0.0001, qtl_var = qtl_var,
      tau = 0.5, pd = 0.15, verbose = FALSE
    )$N
  }, numeric(1))
  expected_threshold <- vapply(qtl_values, function(qtl_var) {
    qtl_threshold_chisq_mssn(
      power = 0.8, alpha = 0.0001, qtl_var = qtl_var,
      tau = 0.5, pd = 0.15, x_upper = 5, x_lower = 5,
      k = 1, verbose = FALSE
    )$N_total
  }, numeric(1))

  expect_equal(names(anova_dat), c("qtl_var", "required_N"))
  expect_equal(anova_dat$required_N, expected_anova)
  expect_equal(names(threshold_dat), c("qtl_var", "MSSN_total"))
  expect_equal(threshold_dat$MSSN_total, expected_threshold)
})

test_that("Falconer plot wrappers validate sweep inputs", {
  expect_error(
    plot_qtl_anova_power(
      "bad", c(1, 2), alpha = 0.05, qtl_var = 0.025, tau = 0.5, pd = 0.15
    ),
    "Unsupported x_var"
  )
  expect_error(
    plot_qtl_anova_power(
      "N", 100, alpha = 0.05, qtl_var = 0.025, tau = 0.5, pd = 0.15
    ),
    "length at least 2"
  )
  expect_error(
    plot_qtl_anova_mssn(
      "N", c(100, 200), power = 0.8, alpha = 0.05,
      qtl_var = 0.025, tau = 0.5, pd = 0.15
    ),
    "Unsupported x_var"
  )
  expect_error(
    plot_qtl_threshold_chisq_mssn(
      "N_case", c(100, 200), power = 0.8, alpha = 0.05,
      qtl_var = 0.025, tau = 0.5, pd = 0.15,
      x_upper = 5, x_lower = 5
    ),
    "Unsupported x_var"
  )
})
