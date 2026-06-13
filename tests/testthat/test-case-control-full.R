cc_g_case <- c(0.25, 0.50, 0.25)
cc_g_ctrl <- c(0.36, 0.48, 0.16)

test_that("cc_power_conditional_full returns expected structure and class", {
  out <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    verbose = FALSE
  )

  expect_s3_class(out, "cc_power_conditional_full")
  expect_type(out$tests, "list")
  expect_type(out$freqs, "list")
  expect_true(all(c("genotypes", "alleles", "trend") %in% names(out$tests)))
})

test_that("cc_mssn_conditional_full returns expected structure and class", {
  out <- cc_mssn_conditional_full(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    verbose = FALSE
  )

  expect_s3_class(out, "cc_mssn_conditional_full")
  expect_type(out$tests, "list")
  expect_type(out$freqs, "list")
  expect_true(all(c("genotypes", "alleles", "trend") %in% names(out$tests)))
})

test_that("model_free and model_based inputs run without error", {
  expect_no_error(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      geno_misclass = "none",
      verbose = FALSE
    )
  )

  expect_no_error(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_based",
      prev = 0.05,
      pd = 0.30,
      R2 = 1.8,
      MOI = "M",
      geno_misclass = "none",
      verbose = FALSE
    )
  )
})

test_that("locus heterogeneity transforms true case frequencies", {
  out <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    locus_het = TRUE,
    pi = 0.75,
    verbose = FALSE
  )

  expected <- 0.75 * out$freqs$g_base_case + 0.25 * out$freqs$g_base_ctrl
  expect_equal(out$freqs$g_true_case, expected, tolerance = 1e-12)
})

test_that("pi boundary values transform true case frequencies", {
  pi_one <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    locus_het = TRUE,
    pi = 1,
    verbose = FALSE
  )

  expect_equal(pi_one$freqs$g_true_case, pi_one$freqs$g_base_case, tolerance = 1e-12)

  pi_zero <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    locus_het = TRUE,
    pi = 0,
    geno_misclass = "diff3p",
    diff_source = "explicit",
    case_e01 = 0.02,
    case_e02 = 0.01,
    case_e03 = 0.005,
    ctrl_e01 = 0,
    ctrl_e02 = 0,
    ctrl_e03 = 0,
    verbose = FALSE
  )

  expect_equal(pi_zero$freqs$g_true_case, pi_zero$freqs$g_base_ctrl, tolerance = 1e-12)
})

test_that("3p genotype misclassification stores model and prints status", {
  out <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    geno_misclass = "3p",
    e01 = 0.02,
    e02 = 0.01,
    e03 = 0.005,
    verbose = FALSE
  )

  expect_equal(out$errors$genotype_misclass$model, "3p_homhet_homhom")

  printed <- capture.output(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      geno_misclass = "3p",
      e01 = 0.02,
      e02 = 0.01,
      e03 = 0.005,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_true(any(grepl("3-parameter", printed, fixed = TRUE)))
})

test_that("diff3p multiplier shortcut scales parameters from case or control", {
  from_case <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    geno_misclass = "diff3p",
    diff_source = "case",
    diff_multiplier = 0.5,
    case_e01 = 0.02,
    case_e02 = 0.01,
    case_e03 = 0.004,
    verbose = FALSE
  )

  expect_equal(
    from_case$errors$genotype_misclass$ctrl_params,
    c(e01 = 0.01, e02 = 0.005, e03 = 0.002),
    tolerance = 1e-12
  )

  from_ctrl <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    geno_misclass = "diff3p",
    diff_source = "ctrl",
    diff_multiplier = 2,
    ctrl_e01 = 0.01,
    ctrl_e02 = 0.005,
    ctrl_e03 = 0.002,
    verbose = FALSE
  )

  expect_equal(
    from_ctrl$errors$genotype_misclass$case_params,
    c(e01 = 0.02, e02 = 0.01, e03 = 0.004),
    tolerance = 1e-12
  )
})

test_that("include_allelic FALSE omits allelic result and printed output", {
  out <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    include_allelic = FALSE,
    verbose = FALSE
  )

  expect_null(out$tests$alleles)

  printed_power <- capture.output(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      include_allelic = FALSE,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_false(any(grepl("Alleles:", printed_power, fixed = TRUE)))
  expect_false(any(grepl("allelic", printed_power, ignore.case = TRUE)))

  printed_mssn <- capture.output(
    cc_mssn_conditional_full(
      power = 0.80,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      include_allelic = FALSE,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_false(any(grepl("Alleles:", printed_mssn, fixed = TRUE)))
  expect_false(any(grepl("allelic", printed_mssn, ignore.case = TRUE)))
})

test_that("invalid inputs throw errors", {
  expect_error(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = c(0.2, 0.2, 0.2),
      g0 = cc_g_ctrl,
      verbose = FALSE
    ),
    "must sum to 1"
  )

  expect_error(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      locus_het = TRUE,
      pi = 1.1,
      verbose = FALSE
    ),
    "pi must be a single number in \\[0,1\\]"
  )

  expect_error(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      w = c(1, 1, 1),
      verbose = FALSE
    ),
    "Trend weights w cannot all be equal"
  )

  expect_error(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      geno_misclass = "3p",
      e01 = 0.8,
      e02 = 0.01,
      e03 = 0.3,
      verbose = FALSE
    ),
    "Need e01 \\+ e03 <= 1"
  )
})

test_that("returned numeric fields are finite and positive", {
  power_out <- cc_power_conditional_full(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    verbose = FALSE
  )

  expect_true(is.finite(power_out$N_case))
  expect_true(is.finite(power_out$N_ctrl))
  expect_true(is.finite(power_out$N_total))
  expect_gt(power_out$N_case, 0)
  expect_gt(power_out$N_ctrl, 0)
  expect_gt(power_out$N_total, 0)
  expect_gt(power_out$tests$genotypes$lambda, 0)
  expect_gt(power_out$tests$genotypes$power, 0)
  expect_gt(power_out$tests$trend$lambda, 0)
  expect_gt(power_out$tests$trend$power, 0)

  mssn_out <- cc_mssn_conditional_full(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    verbose = FALSE
  )

  expect_gt(mssn_out$tests$genotypes$MSSN_case, 0)
  expect_gt(mssn_out$tests$genotypes$MSSN_ctrl, 0)
  expect_gt(mssn_out$tests$genotypes$MSSN_total, 0)
  expect_gt(mssn_out$tests$trend$MSSN_case, 0)
  expect_gt(mssn_out$tests$trend$MSSN_ctrl, 0)
  expect_gt(mssn_out$tests$trend$MSSN_total, 0)
})

test_that("verbose output is clean and reports modifiers", {
  printed <- capture.output(
    cc_power_conditional_full(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      locus_het = TRUE,
      pi = 0.75,
      geno_misclass = "3p",
      e01 = 0.02,
      e02 = 0.01,
      e03 = 0.005,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_false(any(grepl("S component", printed, fixed = TRUE)))
  expect_true(any(grepl("Locus heterogeneity:", printed, fixed = TRUE)))
  expect_true(any(grepl("Genotype misclassification:", printed, fixed = TRUE)))
})
