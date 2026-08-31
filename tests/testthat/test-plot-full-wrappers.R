plot_g_aff <- c((1 - 0.05)^2, 2 * 0.05 * (1 - 0.05), 0.05^2)
plot_g_unaff <- c((1 - 0.15)^2, 2 * 0.15 * (1 - 0.15), 0.15^2)

test_that("plot_cc_power returns ggplot for phenotype misclassification sweep", {
  p <- plot_cc_power(
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

test_that("plot_cc_mssn returns ggplot for model-based genotype error sweep", {
  p <- plot_cc_mssn(
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

test_that("TDT full wrappers return ggplot objects", {
  power_plot <- plot_tdt_power(
    x_var = "heter_rate",
    x_values = c(0, 0.05, 0.10),
    scenario = "heterogeneity",
    N = 600,
    pd = 0.30,
    prev = 0.05,
    R1 = 1.5,
    R2 = 2.25,
    alpha = 0.05,
    delta_prime = 1
  )

  mssn_plot <- plot_tdt_mssn(
    x_var = "misclass_rate",
    x_values = c(0, 0.01, 0.02),
    scenario = "misclassification",
    target_power = 0.80,
    pd = 0.30,
    prev = 0.05,
    R1 = 1.5,
    R2 = 2.25,
    alpha = 0.05,
    delta_prime = 1
  )

  expect_s3_class(power_plot, "ggplot")
  expect_s3_class(mssn_plot, "ggplot")
})

test_that("plot_cc_power compare_tests returns plot and readable data labels", {
  x_values <- c(0, 0.01, 0.02)

  p <- plot_cc_power(
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

  dat <- plot_cc_power(
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
  expect_equal(nrow(dat), length(x_values) * 2)
  expect_setequal(
    dat$test,
    c("Genotype chi-square", "Trend test")
  )
})

test_that("return_data TRUE returns simple data frames", {
  x_values <- c(0, 0.01, 0.02)

  dat <- plot_cc_power(
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
  tdt_p <- plot_tdt_power(
    x_var = "heter_rate",
    x_values = c(0, 0.05),
    scenario = "heterogeneity",
    N = 600,
    pd = 0.30,
    prev = 0.05,
    R1 = 1.5,
    R2 = 2.25,
    alpha = 0.05,
    delta_prime = 1
  )

  cc_p <- plot_cc_power(
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

  expect_equal(tdt_p$labels$x, "Heterogeneity rate")
  expect_match(cc_p$labels$x, "Phenotype misclassification")
})

test_that("plot wrappers report unsupported arguments clearly", {
  expect_error(
    plot_cc_power(
      x_var = "not_supported",
      x_values = c(0, 1),
      input_mode = "model_free"
    ),
    "Unsupported x_var"
  )

  expect_error(
    plot_cc_power(
      x_var = "phi",
      x_values = c(0, 0.01),
      test = "bad_test",
      input_mode = "model_free"
    ),
    "'arg' should be one of"
  )

  expect_error(
    plot_tdt_power(
      x_var = "heter_rate",
      x_values = c(0, 0.1),
      scenario = "bad_scenario"
    ),
    "'arg' should be one of"
  )
})

test_that("plot_tdt_power model_based data matches tdt_power() directly", {
  x_values <- c(0, 0.05, 0.10, 0.20)

  dat <- plot_tdt_power(
    x_var = "heter_rate",
    x_values = x_values,
    scenario = "heterogeneity",
    input_mode = "model_based",
    N = 600, pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1,
    return_data = TRUE
  )

  expected <- vapply(x_values, function(hr) {
    suppressMessages(tdt_power(
      N = 600, pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
      alpha = 0.05, delta_prime = 1, heter_rate = hr,
      verbose = FALSE
    ))$power$heterogeneity
  }, numeric(1))

  expect_equal(dat$power, expected)
})

test_that("plot_tdt_mssn model_based data matches tdt_mssn() directly", {
  x_values <- c(0, 0.01, 0.02, 0.05)

  dat <- plot_tdt_mssn(
    x_var = "misclass_rate",
    x_values = x_values,
    scenario = "misclassification",
    input_mode = "model_based",
    target_power = 0.80, pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1,
    return_data = TRUE
  )

  expected <- vapply(x_values, function(mr) {
    suppressMessages(tdt_mssn(
      target_power = 0.80, pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
      alpha = 0.05, delta_prime = 1, misclass_rate = mr,
      verbose = FALSE
    ))$N$misclassification
  }, numeric(1))

  expect_equal(dat$required_trios, expected)
})

test_that("plot_tdt_power works in model_free mode sweeping ET and ENT", {
  dat_et <- plot_tdt_power(
    x_var = "ET",
    x_values = c(120, 140, 160),
    scenario = "no_error",
    input_mode = "model_free",
    N = 500, ENT = 100,
    return_data = TRUE
  )
  expect_equal(nrow(dat_et), 3)
  expect_true(all(diff(dat_et$power) > 0))

  dat_ent <- plot_tdt_power(
    x_var = "ENT",
    x_values = c(60, 80, 100),
    scenario = "no_error",
    input_mode = "model_free",
    N = 500, ET = 160,
    return_data = TRUE
  )
  expect_equal(nrow(dat_ent), 3)
  expect_true(all(diff(dat_ent$power) < 0))
})

test_that("plot_tdt_mssn works in model_free mode sweeping ET, ENT, and n_trios", {
  expect_no_error(
    plot_tdt_mssn(
      x_var = "ET",
      x_values = c(120, 140, 160),
      scenario = "no_error",
      input_mode = "model_free",
      target_power = 0.80, ENT = 100, n_trios = 120
    )
  )

  expect_no_error(
    plot_tdt_mssn(
      x_var = "n_trios",
      x_values = c(100, 120, 150),
      scenario = "no_error",
      input_mode = "model_free",
      target_power = 0.80, ET = 140, ENT = 100
    )
  )
})

test_that("TDT plot wrappers validate x_var against input_mode", {
  expect_error(
    plot_tdt_power(
      x_var = "pd",
      x_values = c(0.1, 0.2),
      input_mode = "model_free",
      N = 500, ET = 140, ENT = 100
    ),
    "model_based"
  )

  expect_error(
    plot_tdt_power(
      x_var = "ET",
      x_values = c(100, 150),
      input_mode = "model_based",
      N = 500, pd = 0.3, prev = 0.05, R1 = 1.5, R2 = 2.25
    ),
    "model_free"
  )

  expect_error(
    plot_tdt_mssn(
      x_var = "prev",
      x_values = c(0.02, 0.05),
      input_mode = "model_free",
      target_power = 0.8, ET = 140, ENT = 100, n_trios = 120
    ),
    "model_based"
  )

  expect_error(
    plot_tdt_mssn(
      x_var = "n_trios",
      x_values = c(100, 150),
      input_mode = "model_based",
      target_power = 0.8, pd = 0.3, prev = 0.05, R1 = 1.5, R2 = 2.25
    ),
    "model_free"
  )
})

test_that("TDT plot wrappers require pd/prev for model_free heterogeneity/misclassification scenarios", {
  expect_error(
    plot_tdt_power(
      x_var = "heter_rate",
      x_values = c(0, 0.1),
      scenario = "heterogeneity",
      input_mode = "model_free",
      N = 500, ET = 140, ENT = 100
    ),
    "pd"
  )

  expect_error(
    plot_tdt_power(
      x_var = "misclass_rate",
      x_values = c(0, 0.1),
      scenario = "misclassification",
      input_mode = "model_free",
      N = 500, ET = 140, ENT = 100, pd = 0.3
    ),
    "prev"
  )

  expect_no_error(
    plot_tdt_power(
      x_var = "heter_rate",
      x_values = c(0, 0.1),
      scenario = "heterogeneity",
      input_mode = "model_free",
      N = 500, ET = 140, ENT = 100, pd = 0.3
    )
  )
})

test_that("model_free round-trip reproduces a model_based power sweep exactly", {
  x_values <- c(0, 0.05, 0.10, 0.20)

  mb <- plot_tdt_power(
    x_var = "heter_rate",
    x_values = x_values,
    scenario = "heterogeneity",
    input_mode = "model_based",
    N = 600, pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1,
    return_data = TRUE
  )

  baseline <- tdt_power(
    N = 600, input_mode = "model_based",
    pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1,
    verbose = FALSE
  )

  mf <- plot_tdt_power(
    x_var = "heter_rate",
    x_values = x_values,
    scenario = "heterogeneity",
    input_mode = "model_free",
    N = 600, ET = baseline$ET$no_error, ENT = baseline$ENT$no_error,
    pd = 0.30,
    return_data = TRUE
  )

  expect_equal(mf$power, mb$power, tolerance = 1e-6)
})

test_that("model_free round-trip reproduces a model_based MSSN sweep exactly", {
  x_values <- c(0, 0.02, 0.05, 0.10)

  mb <- plot_tdt_mssn(
    x_var = "misclass_rate",
    x_values = x_values,
    scenario = "misclassification",
    input_mode = "model_based",
    target_power = 0.80, pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1,
    return_data = TRUE
  )

  baseline <- tdt_mssn(
    target_power = 0.80, input_mode = "model_based",
    pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1,
    verbose = FALSE
  )

  n_trios_used <- baseline$N$no_error
  ET_val  <- baseline$gT_star$no_error  * 2 * n_trios_used
  ENT_val <- baseline$gNT_star$no_error * 2 * n_trios_used

  mf <- plot_tdt_mssn(
    x_var = "misclass_rate",
    x_values = x_values,
    scenario = "misclassification",
    input_mode = "model_free",
    target_power = 0.80, ET = ET_val, ENT = ENT_val, n_trios = n_trios_used,
    pd = 0.30, prev = 0.05,
    return_data = TRUE
  )

  expect_equal(mf$required_trios, mb$required_trios, tolerance = 1e-4)
})

test_that("plot wrappers support genotype and phenotype error multipliers", {
  expect_no_error(
    plot_cc_mssn(
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
    plot_cc_power(
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

test_that("ordinary CC pi power sweeps delegate exactly and retain the null", {
  pis <- c(1, 0.75, 0.5, 0.25, 0)
  fixed <- list(
    N_case = 1000, alpha = 0.05, prev = 0.05, pd = 0.30,
    R2 = 1.8, MOI = "M", k = 1, locus_het = TRUE, verbose = FALSE
  )
  dat <- do.call(plot_cc_power, c(list(
    x_var = "pi", x_values = pis, test = "trend", return_data = TRUE
  ), fixed))
  direct <- lapply(pis, function(pi_i) {
    do.call(cc_power, c(fixed, list(pi = pi_i)))$tests$trend
  })

  expect_equal(dat$power, vapply(direct, `[[`, numeric(1), "power"))
  expect_equal(dat$lambda, vapply(direct, `[[`, numeric(1), "lambda"))
  expect_equal(dat$power[dat$pi == 0], fixed$alpha, tolerance = 1e-14)
  expect_equal(dat$lambda[dat$pi == 0], 0)

  p <- do.call(plot_cc_power, c(list(
    x_var = "pi", x_values = pis, test = "trend"
  ), fixed))
  expect_s3_class(p, "ggplot")
  expect_true(any(p$data$x == 0))
  expect_equal(p$data$lambda[p$data$x == 0], 0)
})

test_that("ordinary CC pi MSSN sweeps retain structural non-finiteness", {
  pis <- c(1, 0.75, 0.5, 0.25, 0)
  fixed <- list(
    power = 0.80, alpha = 0.05, prev = 0.05, pd = 0.30,
    R2 = 1.8, MOI = "M", k = 1, locus_het = TRUE, verbose = FALSE
  )
  dat <- do.call(plot_cc_mssn, c(list(
    x_var = "pi", x_values = pis, test = "trend",
    sample_size = "total", return_data = TRUE
  ), fixed))
  direct <- vapply(pis[pis > 0], function(pi_i) {
    do.call(cc_mssn, c(fixed, list(pi = pi_i)))$tests$trend$MSSN_total
  }, numeric(1))

  expect_equal(dat$MSSN_total[dat$pi > 0], direct)
  null <- dat[dat$pi == 0, ]
  expect_true(is.na(null$MSSN_total))
  expect_false(null$finite_mssn)
  expect_identical(null$status, "no finite MSSN")

  p <- do.call(plot_cc_mssn, c(list(
    x_var = "pi", x_values = pis, test = "trend", sample_size = "case"
  ), fixed))
  expect_s3_class(p, "ggplot")
  expect_true(any(!p$data$finite_mssn))
  expect_true(is.na(p$data$y[p$data$x == 0]))
})

test_that("ordinary CC plot master switch and unexpected errors propagate", {
  common <- list(
    alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
    MOI = "M", k = 1, locus_het = FALSE
  )
  expect_error(
    do.call(plot_cc_power, c(list(
      x_var = "pi", x_values = c(1, 0.5), test = "trend", N_case = 1000
    ), common)),
    "pi is used only when locus_het = TRUE; set pi = 1 or enable locus heterogeneity\\."
  )
  expect_error(
    do.call(plot_cc_mssn, c(list(
      x_var = "pi", x_values = c(1, 1.1), test = "trend", power = 0.8
    ), common)),
    "pi must be"
  )
})

test_that("ordinary CC historical non-locus plot data structure is unchanged", {
  dat <- plot_cc_power(
    x_var = "N_case", x_values = c(500, 1000), test = "trend",
    alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
    MOI = "M", k = 1, locus_het = FALSE, pi = 1,
    return_data = TRUE
  )
  expect_identical(names(dat), c("N_case", "test", "power"))
  direct <- vapply(c(500, 1000), function(n) {
    cc_power(
      N_case = n, alpha = 0.05, prev = 0.05, pd = 0.30,
      R2 = 1.8, MOI = "M", k = 1, locus_het = FALSE, pi = 1,
      verbose = FALSE
    )$tests$trend$power
  }, numeric(1))
  expect_equal(dat$power, direct)
})
