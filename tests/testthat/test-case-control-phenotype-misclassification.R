cc_pheno_g_aff <- c((1 - 0.05)^2, 2 * 0.05 * (1 - 0.05), 0.05^2)
cc_pheno_g_unaff <- c((1 - 0.15)^2, 2 * 0.15 * (1 - 0.15), 0.15^2)
cc_pheno_g_case <- c(0.25, 0.50, 0.25)
cc_pheno_g_ctrl <- c(0.36, 0.48, 0.16)

test_that("phenotype misclassification defaults preserve full-function outputs", {
  old_power <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_pheno_g_case,
    g0 = cc_pheno_g_ctrl,
    verbose = FALSE
  )

  explicit_power <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_pheno_g_case,
    g0 = cc_pheno_g_ctrl,
    pheno_misclass = FALSE,
    theta = 0,
    phi = 0,
    verbose = FALSE
  )

  expect_equal(explicit_power$tests$genotypes, old_power$tests$genotypes)
  expect_equal(explicit_power$tests$trend, old_power$tests$trend)

  old_mssn <- cc_mssn_conditional_full(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_pheno_g_case,
    g0 = cc_pheno_g_ctrl,
    verbose = FALSE
  )

  explicit_mssn <- cc_mssn_conditional_full(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_pheno_g_case,
    g0 = cc_pheno_g_ctrl,
    pheno_misclass = FALSE,
    theta = 0,
    phi = 0,
    verbose = FALSE
  )

  expect_equal(explicit_mssn$tests$genotypes, old_mssn$tests$genotypes)
  expect_equal(explicit_mssn$tests$trend, old_mssn$tests$trend)
})

test_that("phenotype misclassification matches Edwards et al Figure 2 genotype power", {
  cases <- list(
    list(prev = 0.05, theta = 0, phi = 0.01, expected = 0.9134913),
    list(prev = 0.05, theta = 0, phi = 0.02, expected = 0.7633919),
    list(prev = 0.01, theta = 0, phi = 0.01, expected = 0.3320457),
    list(prev = 0.01, theta = 0, phi = 0.02, expected = 0.1098263),
    list(prev = 0.05, theta = 0.15, phi = 0, expected = 0.988711),
    list(prev = 0.01, theta = 0.15, phi = 0, expected = 0.9894593)
  )

  for (case in cases) {
    out <- cc_power_conditional_full(
      N_case = 250,
      alpha = 0.01,
      k = 1,
      input_mode = "model_free",
      g1 = cc_pheno_g_aff,
      g0 = cc_pheno_g_unaff,
      prev = case$prev,
      pheno_misclass = TRUE,
      theta = case$theta,
      phi = case$phi,
      verbose = FALSE
    )

    expect_equal(out$tests$genotypes$power, case$expected, tolerance = 1e-6)
  }
})

test_that("phenotype misclassification MSSN round trip reaches target power", {
  target_power <- 0.80

  mssn <- cc_mssn_conditional_full(
    power = target_power,
    alpha = 0.01,
    input_mode = "model_free",
    g1 = cc_pheno_g_aff,
    g0 = cc_pheno_g_unaff,
    prev = 0.05,
    pheno_misclass = TRUE,
    theta = 0,
    phi = 0.01,
    verbose = FALSE
  )

  power <- cc_power_conditional_full(
    N_case = mssn$tests$genotypes$MSSN_case,
    alpha = 0.01,
    input_mode = "model_free",
    g1 = cc_pheno_g_aff,
    g0 = cc_pheno_g_unaff,
    prev = 0.05,
    pheno_misclass = TRUE,
    theta = 0,
    phi = 0.01,
    verbose = FALSE
  )

  expect_gte(power$tests$genotypes$power, target_power)
})

test_that("phenotype misclassification validates full-function inputs", {
  expect_error(
    cc_power_conditional_full(
      N_case = 250,
      alpha = 0.01,
      input_mode = "model_free",
      g1 = cc_pheno_g_aff,
      g0 = cc_pheno_g_unaff,
      prev = 0.05,
      pheno_misclass = "yes",
      verbose = FALSE
    ),
    "pheno_misclass must be TRUE or FALSE"
  )

  expect_error(
    cc_power_conditional_full(
      N_case = 250,
      alpha = 0.01,
      input_mode = "model_free",
      g1 = cc_pheno_g_aff,
      g0 = cc_pheno_g_unaff,
      prev = 0.05,
      pheno_misclass = TRUE,
      theta = 1,
      verbose = FALSE
    ),
    "theta must be a single number in \\[0,1\\)"
  )

  expect_error(
    cc_power_conditional_full(
      N_case = 250,
      alpha = 0.01,
      input_mode = "model_free",
      g1 = cc_pheno_g_aff,
      g0 = cc_pheno_g_unaff,
      prev = 0.05,
      pheno_misclass = TRUE,
      phi = 1,
      verbose = FALSE
    ),
    "phi must be a single number in \\[0,1\\)"
  )

  expect_error(
    cc_power_conditional_full(
      N_case = 250,
      alpha = 0.01,
      input_mode = "model_free",
      g1 = cc_pheno_g_aff,
      g0 = cc_pheno_g_unaff,
      pheno_misclass = TRUE,
      theta = 0,
      phi = 0.01,
      verbose = FALSE
    ),
    "prev must be a single number in \\(0,1\\)"
  )
})

test_that("specialized phenotype misclassification functions run and return key fields", {
  power <- cc_pheno_power_test(
    N_case = 250,
    alpha = 0.01,
    g_aff = cc_pheno_g_aff,
    g_unaff = cc_pheno_g_unaff,
    prev = 0.05,
    theta = 0,
    phi = 0.01
  )

  expect_equal(power$power, 0.9134913, tolerance = 1e-7)
  expect_gt(power$lambda, 0)
  expect_equal(sum(power$g_case_obs), 1, tolerance = 1e-12)
  expect_equal(sum(power$g_ctrl_obs), 1, tolerance = 1e-12)

  mssn <- cc_pheno_mssn_test(
    target_power = 0.80,
    alpha = 0.01,
    g_aff = cc_pheno_g_aff,
    g_unaff = cc_pheno_g_unaff,
    prev = 0.05,
    theta = 0,
    phi = 0.01
  )

  expect_gt(mssn$N_case, 0)
  expect_gt(mssn$N_ctrl, 0)
  expect_gt(mssn$N_total, 0)
  expect_gt(mssn$lambda_star, 0)
})
