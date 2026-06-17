plot_g_aff <- c((1 - 0.05)^2, 2 * 0.05 * (1 - 0.05), 0.05^2)
plot_g_unaff <- c((1 - 0.15)^2, 2 * 0.15 * (1 - 0.15), 0.15^2)

test_that("cc_plot_power returns ggplot for phenotype misclassification sweep", {
  p <- cc_plot_power(
    x_var = "phi",
    x_values = c(0, 0.01, 0.02),
    test = "genotypes",
    input_mode = "model_free",
    N_case = 250,
    alpha = 0.01,
    g1 = plot_g_aff,
    g0 = plot_g_unaff,
    prev = 0.05,
    pheno_misclass = TRUE,
    theta = 0,
    k = 1
  )

  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "Phenotype misclassification")
})

test_that("cc_plot_mssn returns ggplot for model-based genotype error sweep", {
  p <- cc_plot_mssn(
    x_var = "geno_error_multiplier",
    x_values = c(0, 0.5, 1),
    test = "trend",
    input_mode = "model_based",
    power = 0.80,
    alpha = 0.05,
    prev = 0.10,
    pd = 0.25,
    R2 = 2,
    MOI = "M",
    geno_misclass = "3p",
    e01_base = 0.02,
    e02_base = 0.01,
    e03_base = 0.005,
    k = 1
  )

  expect_s3_class(p, "ggplot")
})

test_that("cc_plot_power compare_tests returns plot and readable data labels", {
  x_values <- c(0, 0.01, 0.02)

  p <- cc_plot_power(
    x_var = "phi",
    x_values = x_values,
    input_mode = "model_free",
    compare_tests = TRUE,
    N_case = 250,
    alpha = 0.01,
    g1 = plot_g_aff,
    g0 = plot_g_unaff,
    prev = 0.05,
    pheno_misclass = TRUE,
    theta = 0,
    k = 1
  )

  dat <- cc_plot_power(
    x_var = "phi",
    x_values = x_values,
    input_mode = "model_free",
    compare_tests = TRUE,
    N_case = 250,
    alpha = 0.01,
    g1 = plot_g_aff,
    g0 = plot_g_unaff,
    prev = 0.05,
    pheno_misclass = TRUE,
    theta = 0,
    k = 1,
    return_data = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), length(x_values) * 3)
  expect_setequal(
    dat$test,
    c("Genotype chi-square", "Allelic chi-square", "Trend test")
  )
})

test_that("return_data TRUE returns simple data frames", {
  x_values <- c(0, 0.01, 0.02)

  dat <- cc_plot_power(
    x_var = "phi",
    x_values = x_values,
    test = "genotypes",
    input_mode = "model_free",
    N_case = 250,
    alpha = 0.01,
    g1 = plot_g_aff,
    g0 = plot_g_unaff,
    prev = 0.05,
    pheno_misclass = TRUE,
    theta = 0,
    k = 1,
    return_data = TRUE
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), length(x_values))
  expect_equal(names(dat), c("phi", "test", "power"))
  expect_equal(unique(dat$test), "Genotype chi-square")
})

test_that("plot wrappers use clear default labels", {
  cc_p <- cc_plot_power(
    x_var = "phi",
    x_values = c(0, 0.01),
    test = "genotypes",
    input_mode = "model_free",
    N_case = 250,
    alpha = 0.01,
    g1 = plot_g_aff,
    g0 = plot_g_unaff,
    prev = 0.05,
    pheno_misclass = TRUE,
    theta = 0,
    k = 1
  )

  expect_match(cc_p$labels$x, "Phenotype misclassification")
})

test_that("plot wrappers report unsupported arguments clearly", {
  expect_error(
    cc_plot_power(
      x_var = "not_supported",
      x_values = c(0, 1),
      input_mode = "model_free"
    ),
    "Unsupported x_var"
  )

  expect_error(
    cc_plot_power(
      x_var = "phi",
      x_values = c(0, 0.01),
      test = "bad_test",
      input_mode = "model_free"
    ),
    "'arg' should be one of"
  )

})

test_that("plot wrappers support genotype and phenotype error multipliers", {
  expect_no_error(
    cc_plot_mssn(
      x_var = "geno_error_multiplier",
      x_values = c(0, 0.5, 1),
      test = "trend",
      input_mode = "model_based",
      power = 0.80,
      alpha = 0.05,
      prev = 0.10,
      pd = 0.25,
      R2 = 2,
      MOI = "M",
      geno_misclass = "3p",
      e01_base = 0.02,
      e02_base = 0.01,
      e03_base = 0.005,
      k = 1
    )
  )

  expect_no_error(
    cc_plot_power(
      x_var = "pheno_error_multiplier",
      x_values = c(0, 0.5, 1),
      test = "genotypes",
      input_mode = "model_free",
      N_case = 250,
      alpha = 0.01,
      g1 = plot_g_aff,
      g0 = plot_g_unaff,
      prev = 0.05,
      theta_base = 0,
      phi_base = 0.01,
      k = 1
    )
  )
})
