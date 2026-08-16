test_that("renamed specialized TDT 2D plots use the canonical backends", {
  skip_if_not_installed("ggplot2")

  common <- list(
    pd = 0.3, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1
  )

  power_misclassification <- suppressMessages(do.call(
    plot_tdt_power_phenotype_misclassification,
    c(common, list(N = 300, misclass_seq = c(0, 0.02)))
  ))
  power_heterogeneity <- suppressMessages(do.call(
    plot_tdt_power_locus_heterogeneity,
    c(common, list(N = 300, heter_seq = c(0, 0.1)))
  ))
  mssn_misclassification <- suppressMessages(do.call(
    plot_tdt_mssn_phenotype_misclassification,
    c(common, list(target_power = 0.8, misclass_seq = c(0, 0.02)))
  ))
  mssn_heterogeneity <- suppressMessages(do.call(
    plot_tdt_mssn_locus_heterogeneity,
    c(common, list(target_power = 0.8, heter_seq = c(0, 0.1)))
  ))

  expect_s3_class(power_misclassification, "ggplot")
  expect_s3_class(power_heterogeneity, "ggplot")
  expect_s3_class(mssn_misclassification, "ggplot")
  expect_s3_class(mssn_heterogeneity, "ggplot")
})


test_that("renamed transitional TDT 3D plots preserve Plotly output", {
  skip_if_not_installed("plotly")

  common <- list(
    pd_seq = c(0.2, 0.3), prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1
  )

  power_misclassification <- suppressMessages(do.call(
    plot_tdt_power_phenotype_misclassification_3d,
    c(common, list(N = 300, misclass_seq = c(0, 0.02)))
  ))
  power_heterogeneity <- suppressMessages(do.call(
    plot_tdt_power_locus_heterogeneity_3d,
    c(common, list(N = 300, heter_seq = c(0, 0.1)))
  ))
  mssn_misclassification <- suppressMessages(do.call(
    plot_tdt_mssn_phenotype_misclassification_3d,
    c(common, list(target_power = 0.8, misclass_seq = c(0, 0.02)))
  ))
  mssn_heterogeneity <- suppressMessages(do.call(
    plot_tdt_mssn_locus_heterogeneity_3d,
    c(common, list(target_power = 0.8, heter_seq = c(0, 0.1)))
  ))

  expect_s3_class(power_misclassification, "plotly")
  expect_s3_class(power_heterogeneity, "plotly")
  expect_s3_class(mssn_misclassification, "plotly")
  expect_s3_class(mssn_heterogeneity, "plotly")
})
