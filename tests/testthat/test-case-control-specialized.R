cc_spec_prev <- 0.10
cc_spec_pd <- 0.30
cc_spec_R2 <- 1.80
cc_spec_g_case <- c(0.25, 0.50, 0.25)
cc_spec_g_ctrl <- c(0.36, 0.48, 0.16)

test_that("1p genotype misclassification specialized functions return stable structures", {
  mssn <- cc_mssn_genotypes_conditional_misclass(
    power = 0.80,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    e = 0.02,
    verbose = FALSE
  )

  expect_s3_class(mssn, "cc_mssn_genotypes_conditional_misclass")
  expect_equal(mssn$misclassification$model, "symmetric")
  expect_equal(rowSums(mssn$misclassification$M), rep(1, 3), tolerance = 1e-12)
  expect_gt(mssn$N_case, 0)
  expect_gt(mssn$S, 0)

  power <- cc_power_genotypes_conditional_misclass(
    N_case = mssn$N_case,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    e = 0.02,
    verbose = FALSE
  )

  expect_s3_class(power, "cc_power_genotypes_conditional_misclass")
  expect_gt(power$lambda, 0)
  expect_gt(power$power, 0)
  expect_lte(power$power, 1)
})

test_that("2p genotype misclassification specialized functions return stable structures", {
  mssn <- cc_mssn_genotypes_conditional_misclass_2p(
    power = 0.80,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    e1 = 0.01,
    e2 = 0.02,
    verbose = FALSE
  )

  expect_s3_class(mssn, "cc_mssn_genotypes_conditional_misclass_2p")
  expect_equal(mssn$misclassification$model, "two_param_hom_het")
  expect_equal(rowSums(mssn$misclassification$M), rep(1, 3), tolerance = 1e-12)
  expect_gt(mssn$N_case, 0)

  power <- cc_power_genotypes_conditional_misclass_2p(
    N_case = mssn$N_case,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    e1 = 0.01,
    e2 = 0.02,
    verbose = FALSE
  )

  expect_s3_class(power, "cc_power_genotypes_conditional_misclass_2p")
  expect_gt(power$lambda, 0)
  expect_gt(power$power, 0)
  expect_lte(power$power, 1)
})

test_that("3p genotype misclassification specialized functions return stable structures", {
  mssn <- cc_mssn_genotypes_conditional_misclass_3p(
    power = 0.80,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    e01 = 0.02,
    e02 = 0.01,
    e03 = 0.005,
    verbose = FALSE
  )

  expect_s3_class(mssn, "cc_mssn_genotypes_conditional_misclass_3p")
  expect_equal(mssn$misclassification$model, "3p")
  expect_equal(rowSums(mssn$misclassification$M), rep(1, 3), tolerance = 1e-12)
  expect_gt(mssn$MSSN_case, 0)

  power <- cc_power_genotypes_conditional_misclass_3p(
    N_case = mssn$MSSN_case,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    e01 = 0.02,
    e02 = 0.01,
    e03 = 0.005,
    verbose = FALSE
  )

  expect_s3_class(power, "cc_power_genotypes_conditional_misclass_3p")
  expect_gt(power$lambda, 0)
  expect_gt(power$power, 0)
  expect_lte(power$power, 1)
})

test_that("differential 3p genotype misclassification specialized functions preserve matrices", {
  mssn <- cc_mssn_genotypes_conditional_diffmisclass(
    power = 0.80,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    case_e01 = 0.02,
    case_e02 = 0.01,
    case_e03 = 0.005,
    ctrl_e01 = 0.01,
    ctrl_e02 = 0.005,
    ctrl_e03 = 0.002,
    verbose = FALSE
  )

  expect_s3_class(mssn, "cc_mssn_genotypes_conditional_diffmisclass")
  expect_equal(mssn$misclassification$model, "differential_3p")
  expect_equal(rowSums(mssn$misclassification$case_matrix), rep(1, 3), tolerance = 1e-12)
  expect_equal(rowSums(mssn$misclassification$ctrl_matrix), rep(1, 3), tolerance = 1e-12)
  expect_gt(mssn$MSSN_case, 0)
  expect_equal(
    mssn$misclassification$case_params,
    c(e01 = 0.02, e02 = 0.01, e03 = 0.005),
    tolerance = 1e-12
  )

  power <- cc_power_genotypes_conditional_diffmisclass(
    N_case = mssn$MSSN_case,
    alpha = 0.05,
    prev = cc_spec_prev,
    pd = cc_spec_pd,
    R2 = cc_spec_R2,
    MOI = "D",
    case_e01 = 0.02,
    case_e02 = 0.01,
    case_e03 = 0.005,
    ctrl_e01 = 0.01,
    ctrl_e02 = 0.005,
    ctrl_e03 = 0.002,
    verbose = FALSE
  )

  expect_s3_class(power, "cc_power_genotypes_conditional_diffmisclass")
  expect_gt(power$lambda, 0)
  expect_gt(power$power, 0)
  expect_lte(power$power, 1)
})

test_that("locus heterogeneity specialized genotype functions adjust case frequencies", {
  mssn <- cc_mssn_locus_het_genotypes(
    power = 0.80,
    alpha = 0.05,
    g_case_assoc = cc_spec_g_case,
    g_ctrl = cc_spec_g_ctrl,
    pi = 0.75,
    verbose = FALSE
  )

  expected <- 0.75 * cc_spec_g_case + 0.25 * cc_spec_g_ctrl
  expect_s3_class(mssn, "cc_mssn_locus_het_genotypes")
  expect_equal(mssn$freqs$g_case_het, expected, tolerance = 1e-12)
  expect_gt(mssn$N_case, 0)

  power <- cc_power_locus_het_genotypes(
    N_case = mssn$N_case,
    alpha = 0.05,
    g_case_assoc = cc_spec_g_case,
    g_ctrl = cc_spec_g_ctrl,
    pi = 0.75,
    verbose = FALSE
  )

  expect_s3_class(power, "cc_power_locus_het_genotypes")
  expect_gt(power$lambda, 0)
  expect_gt(power$power, 0)
})

test_that("locus heterogeneity specialized allele functions return allele frequencies", {
  mssn <- cc_mssn_locus_het_alleles(
    power = 0.80,
    alpha = 0.05,
    g_case_assoc = cc_spec_g_case,
    g_ctrl = cc_spec_g_ctrl,
    pi = 0.75,
    verbose = FALSE
  )

  expect_s3_class(mssn, "cc_mssn_locus_het_alleles")
  expect_named(mssn$freqs$p_case_het, c("q", "p"))
  expect_named(mssn$freqs$p_ctrl_het, c("q", "p"))
  expect_gt(mssn$N_case, 0)

  power <- cc_power_locus_het_alleles(
    N_case = mssn$N_case,
    alpha = 0.05,
    g_case_assoc = cc_spec_g_case,
    g_ctrl = cc_spec_g_ctrl,
    pi = 0.75,
    verbose = FALSE
  )

  expect_s3_class(power, "cc_power_locus_het_alleles")
  expect_gt(power$lambda, 0)
  expect_gt(power$power, 0)
})

test_that("locus heterogeneity specialized trend functions return numerator and denominator", {
  mssn <- cc_mssn_locus_het_trend(
    power = 0.80,
    alpha = 0.05,
    g_case_assoc = cc_spec_g_case,
    g_ctrl = cc_spec_g_ctrl,
    pi = 0.75,
    w = c(0, 1, 2),
    verbose = FALSE
  )

  expect_s3_class(mssn, "cc_mssn_locus_het_trend")
  expect_gt(mssn$numerator, 0)
  expect_gt(mssn$denominator, 0)
  expect_gt(mssn$N_case, 0)

  power <- cc_power_locus_het_trend(
    N_case = mssn$N_case,
    alpha = 0.05,
    g_case_assoc = cc_spec_g_case,
    g_ctrl = cc_spec_g_ctrl,
    pi = 0.75,
    w = c(0, 1, 2),
    verbose = FALSE
  )

  expect_s3_class(power, "cc_power_locus_het_trend")
  expect_gt(power$lambda, 0)
  expect_gt(power$power, 0)
})

test_that("specialized helper-facing validation errors remain active", {
  expect_error(
    cc_mssn_locus_het_genotypes(
      power = 0.80,
      alpha = 0.05,
      g_case_assoc = c(0.2, 0.2, 0.2),
      g_ctrl = cc_spec_g_ctrl,
      pi = 0.75,
      verbose = FALSE
    ),
    "must sum to 1"
  )

  expect_error(
    cc_mssn_genotypes_conditional_misclass_3p(
      power = 0.80,
      alpha = 0.05,
      prev = cc_spec_prev,
      pd = cc_spec_pd,
      R2 = cc_spec_R2,
      MOI = "D",
      e01 = 0.8,
      e02 = 0.01,
      e03 = 0.3,
      verbose = FALSE
    ),
    "Need e01 \\+ e03 <= 1"
  )
})
