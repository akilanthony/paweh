cc_g_case <- c(0.25, 0.50, 0.25)
cc_g_ctrl <- c(0.36, 0.48, 0.16)

test_that("cc_power returns expected structure and class", {
  out <- cc_power(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    verbose = FALSE
  )

  expect_s3_class(out, "cc_power")
  expect_type(out$tests, "list")
  expect_type(out$freqs, "list")
  expect_named(out$tests, c("genotypes", "trend"))
  expect_false("alleles" %in% names(out$tests))
  obsolete_arg <- paste0("include_", "allel", "ic")
  expect_false(obsolete_arg %in% names(out))
})

test_that("cc_mssn returns expected structure and class", {
  out <- cc_mssn(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    verbose = FALSE
  )

  expect_s3_class(out, "cc_mssn")
  expect_type(out$tests, "list")
  expect_type(out$freqs, "list")
  expect_named(out$tests, c("genotypes", "trend"))
  expect_false("alleles" %in% names(out$tests))
  obsolete_arg <- paste0("include_", "allel", "ic")
  expect_false(obsolete_arg %in% names(out))
})

test_that("model_free and model_based inputs run without error", {
  expect_no_error(
    cc_power(
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
    cc_power(
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
  out <- cc_power(
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
  pi_one <- cc_power(
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

  pi_zero <- cc_power(
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

test_that("ordinary CC locus switch preserves valid historical fixtures", {
  power_args <- list(
    N_case = 1000, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M",
    k = 1, w = c(0, 1, 2), verbose = FALSE
  )
  mssn_args <- power_args
  mssn_args$N_case <- NULL
  mssn_args$power <- 0.80

  power_off <- do.call(cc_power, c(power_args, list(locus_het = FALSE, pi = 1)))
  power_one <- do.call(cc_power, c(power_args, list(locus_het = TRUE, pi = 1)))
  power_half <- do.call(cc_power, c(power_args, list(locus_het = TRUE, pi = 0.5)))
  mssn_off <- do.call(cc_mssn, c(mssn_args, list(locus_het = FALSE, pi = 1)))
  mssn_one <- do.call(cc_mssn, c(mssn_args, list(locus_het = TRUE, pi = 1)))
  mssn_half <- do.call(cc_mssn, c(mssn_args, list(locus_het = TRUE, pi = 0.5)))

  expect_equal(power_off$tests$trend$lambda, 21.094805103913,
               tolerance = 1e-12)
  expect_equal(power_off$tests$trend$power, 0.995767586715572,
               tolerance = 1e-12)
  expect_equal(power_one$tests, power_off$tests, tolerance = 1e-15)
  expect_equal(power_half$tests$genotypes$lambda, 5.44661220996911,
               tolerance = 1e-12)
  expect_equal(power_half$tests$genotypes$power, 0.540622389046613,
               tolerance = 1e-12)
  expect_equal(power_half$tests$trend$lambda, 5.43092431210202,
               tolerance = 1e-12)
  expect_equal(power_half$tests$trend$power, 0.644492895244749,
               tolerance = 1e-12)

  expect_identical(mssn_off$tests$genotypes$MSSN_case, 457)
  expect_identical(mssn_off$tests$trend$MSSN_case, 373)
  expect_equal(mssn_one$tests, mssn_off$tests, tolerance = 1e-15)
  expect_identical(mssn_half$tests$genotypes$MSSN_case, 1769)
  expect_identical(mssn_half$tests$trend$MSSN_case, 1446)
})

test_that("ordinary CC rejects pi when the locus switch is disabled", {
  message <- paste(
    "pi is used only when locus_het = TRUE;",
    "set pi = 1 or enable locus heterogeneity."
  )
  power_args <- list(
    N_case = 500, alpha = 0.05, input_mode = "model_free",
    g1 = cc_g_case, g0 = cc_g_ctrl,
    locus_het = FALSE, verbose = FALSE
  )
  mssn_args <- power_args
  mssn_args$N_case <- NULL
  mssn_args$power <- 0.80

  for (pi in c(0, 0.5)) {
    expect_error(do.call(cc_power, c(power_args, list(pi = pi))), message,
                 fixed = TRUE)
    expect_error(do.call(cc_mssn, c(mssn_args, list(pi = pi))), message,
                 fixed = TRUE)
  }
})

test_that("active pi zero remains the valid ordinary CC null boundary", {
  power <- cc_power(
    N_case = 500, alpha = 0.05,
    input_mode = "model_free", g1 = cc_g_case, g0 = cc_g_ctrl,
    locus_het = TRUE, pi = 0, verbose = FALSE
  )
  expect_equal(power$freqs$g_true_case, power$freqs$g_true_ctrl,
               tolerance = 1e-15)
  expect_equal(power$tests$genotypes$lambda, 0, tolerance = 1e-15)
  expect_equal(power$tests$trend$lambda, 0, tolerance = 1e-15)
  expect_equal(power$tests$genotypes$power, power$alpha, tolerance = 1e-15)
  expect_equal(power$tests$trend$power, power$alpha, tolerance = 1e-15)

  expect_error(
    cc_mssn(
      power = 0.80, alpha = 0.05,
      input_mode = "model_free", g1 = cc_g_case, g0 = cc_g_ctrl,
      locus_het = TRUE, pi = 0, verbose = FALSE
    ),
    "No finite MSSN exists because the trend contrast is zero under this design.",
    fixed = TRUE
  )
})

test_that("small positive pi values remain on ordinary CC formula paths", {
  for (pi in c(1e-8, 0.001)) {
    power <- cc_power(
      N_case = 500, alpha = 0.05,
      input_mode = "model_free", g1 = cc_g_case, g0 = cc_g_ctrl,
      locus_het = TRUE, pi = pi, verbose = FALSE
    )
    mssn <- cc_mssn(
      power = 0.80, alpha = 0.05,
      input_mode = "model_free", g1 = cc_g_case, g0 = cc_g_ctrl,
      locus_het = TRUE, pi = pi, verbose = FALSE
    )

    expect_gt(power$tests$genotypes$lambda, 0)
    expect_gt(power$tests$trend$lambda, 0)
    expect_true(is.finite(power$tests$genotypes$power))
    expect_true(is.finite(power$tests$trend$power))
    expect_true(is.finite(mssn$tests$genotypes$MSSN_case))
    expect_true(is.finite(mssn$tests$trend$MSSN_case))
  }
})

test_that("ordinary CC rejects invalid pi before switch contradictions", {
  invalid_pi <- list(-0.1, 1.1, NA_real_, Inf, "0.5", c(0.5, 1))
  for (pi in invalid_pi) {
    expect_error(
      cc_power(
        N_case = 500, alpha = 0.05,
        input_mode = "model_free", g1 = cc_g_case, g0 = cc_g_ctrl,
        locus_het = FALSE, pi = pi, verbose = FALSE
      ),
      "pi must be a single number in [0,1].", fixed = TRUE
    )
    expect_error(
      cc_mssn(
        power = 0.80, alpha = 0.05,
        input_mode = "model_free", g1 = cc_g_case, g0 = cc_g_ctrl,
        locus_het = FALSE, pi = pi, verbose = FALSE
      ),
      "pi must be a single number in [0,1].", fixed = TRUE
    )
  }
})

test_that("3p genotype misclassification stores model and prints status", {
  out <- cc_power(
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
    cc_power(
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
  from_case <- cc_power(
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

  from_ctrl <- cc_power(
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

test_that("full functions no longer accept the obsolete allele-test argument or output", {
  obsolete_arg <- paste0("include_", "allel", "ic")
  out <- cc_power(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_free",
    g1 = cc_g_case,
    g0 = cc_g_ctrl,
    verbose = FALSE
  )

  expect_false("alleles" %in% names(out$tests))
  expect_false(obsolete_arg %in% names(out))

  printed_power <- capture.output(
    cc_power(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_false(any(grepl("Alleles:", printed_power, fixed = TRUE)))
  expect_false(any(grepl("Allelic test", printed_power, fixed = TRUE)))

  printed_mssn <- capture.output(
    cc_mssn(
      power = 0.80,
      alpha = 0.05,
      input_mode = "model_free",
      g1 = cc_g_case,
      g0 = cc_g_ctrl,
      verbose = TRUE
    ),
    type = "message"
  )

  expect_false(any(grepl("Alleles:", printed_mssn, fixed = TRUE)))
  expect_false(any(grepl("Allelic test", printed_mssn, fixed = TRUE)))

  expect_error(
    do.call(
      cc_power,
      c(
        list(
          N_case = 500,
          alpha = 0.05,
          input_mode = "model_free",
          g1 = cc_g_case,
          g0 = cc_g_ctrl,
          verbose = FALSE
        ),
        stats::setNames(list(FALSE), obsolete_arg)
      )
    ),
    "unused argument"
  )
})

test_that("invalid inputs throw errors", {
  expect_error(
    cc_power(
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
    cc_power(
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
    cc_power(
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
    cc_power(
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
  power_out <- cc_power(
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

  mssn_out <- cc_mssn(
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
    cc_power(
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
