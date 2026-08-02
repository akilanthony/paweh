# Tests for the generalized conditional TDT framework
#
# Anchor parameter set: pd = 0.25, prev = 0.005, R1 = R2 = 2, delta_prime = 1.
# The existing test for tdt_required_trios() fixes gT* = 0.26087 for these
# inputs, so the new framework must reproduce that value.

tdt_pars <- list(
  prev = 0.005,
  pd = 0.25,
  R1 = 2,
  R2 = 2,
  delta_prime = 1
)


test_that("tdt_power_conditional_full returns expected structure and class", {
  out <- tdt_power_conditional_full(
    n_trios = 500,
    alpha = 0.05,
    input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    verbose = FALSE
  )

  expect_s3_class(out, "tdt_power_conditional_full")
  expect_type(out$tests, "list")
  expect_type(out$transmissions, "list")
  expect_named(out$tests, "tdt")
  expect_equal(out$tests$tdt$df, 1)
  expect_true(all(c("locus_het", "errors", "model_info") %in% names(out)))
  expect_true(out$tests$tdt$power > 0 && out$tests$tdt$power <= 1)
})


test_that("tdt_mssn_conditional_full returns expected structure and class", {
  out <- tdt_mssn_conditional_full(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    verbose = FALSE
  )

  expect_s3_class(out, "tdt_mssn_conditional_full")
  expect_type(out$tests, "list")
  expect_type(out$transmissions, "list")
  expect_named(out$tests, "tdt")
  expect_equal(out$tests$tdt$df, 1)
  expect_true(all(c("locus_het", "errors", "model_info") %in% names(out)))
  expect_true(out$tests$tdt$MSSN_trios > 0)
  expect_equal(out$tests$tdt$MSSN_trios, ceiling(out$tests$tdt$MSSN_trios))
})


test_that("model_based and model_free inputs both run without error", {
  expect_no_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05,
      input_mode = "model_based",
      prev = 0.05, pd = 0.25, R2 = 1.5, MOI = "M",
      verbose = FALSE
    )
  )

  expect_no_error(
    tdt_power_conditional_full(
      n_trios = 120, alpha = 0.05,
      input_mode = "model_free",
      ET = 140, ENT = 100,
      verbose = FALSE
    )
  )

  expect_no_error(
    tdt_mssn_conditional_full(
      power = 0.8, alpha = 0.05,
      input_mode = "model_free",
      ET = 140, ENT = 100, n_trios = 120,
      verbose = FALSE
    )
  )
})


test_that("model-based transmission probabilities match the documented anchor", {
  out <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05,
    input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    verbose = FALSE
  )

  expect_equal(round(out$transmissions$gT_obs, 5), 0.26087, tolerance = 1e-3)
  expect_true(out$transmissions$gT_obs > out$transmissions$gNT_obs)
  expect_equal(out$transmissions$ET, 2 * 500 * out$transmissions$gT_obs)
  expect_equal(out$transmissions$ENT, 2 * 500 * out$transmissions$gNT_obs)
})


test_that("prevalence recomputed from penetrances equals the supplied prevalence", {
  out <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05, pd = 0.10, R1 = 1.5, R2 = 2.25,
    delta_prime = 0.8,
    verbose = FALSE
  )

  expect_equal(out$model_info$phi1_check, 0.05, tolerance = 1e-12)
})


test_that("power agrees with tdt_power_from_model", {
  new <- tdt_power_conditional_full(
    n_trios = 10000, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05, pd = 0.25, R1 = 1, R2 = 1.1,
    delta_prime = 1,
    verbose = FALSE
  )

  old <- suppressMessages(
    tdt_power_from_model(
      pd = 0.25, N = 10000, delta_prime = 1,
      prev = 0.05, R1 = 1, R2 = 1.1, alpha = 0.05
    )
  )

  expect_equal(new$tests$tdt$power, old$`Power`, tolerance = 1e-10)
  expect_equal(
    new$tests$tdt$lambda,
    old$`Non-Centrality Parameter (lambda)`,
    tolerance = 1e-10
  )
  expect_equal(
    new$transmissions$ET,
    old$`Expected Transmissions (ET)`,
    tolerance = 1e-8
  )
  expect_equal(
    new$transmissions$ENT,
    old$`Expected Non-Transmissions (ENT)`,
    tolerance = 1e-8
  )
})


test_that("MSSN agrees with tdt_required_trios at pi = 1", {
  new <- tdt_mssn_conditional_full(
    power = 0.8, alpha = 0.05,
    input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    verbose = FALSE
  )

  old <- suppressMessages(
    tdt_required_trios(
      power = 0.8, alpha = 0.05, df = 1,
      pd = tdt_pars$pd, prev = tdt_pars$prev,
      R1 = tdt_pars$R1, R2 = tdt_pars$R2,
      delta_prime = tdt_pars$delta_prime, pi = 1
    )
  )

  expect_equal(
    new$tests$tdt$MSSN_trios,
    ceiling(old$`Required Number of Trios (N_star)`)
  )
  expect_equal(
    new$transmissions$gT_obs,
    old$`Expected Transmission (gT_star)`,
    tolerance = 1e-10
  )
  expect_equal(
    new$transmissions$gNT_obs,
    old$`Expected Non-Transmission (gNT_star)`,
    tolerance = 1e-10
  )
})


test_that("model_free reproduces model_based when fed its own ET and ENT", {
  mb <- tdt_power_conditional_full(
    n_trios = 800, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.20, pd = 0.40, R1 = 1.2, R2 = 1.44,
    delta_prime = 0.5,
    verbose = FALSE
  )

  mf <- tdt_power_conditional_full(
    n_trios = 800, alpha = 0.05,
    input_mode = "model_free",
    ET = mb$transmissions$ET,
    ENT = mb$transmissions$ENT,
    verbose = FALSE
  )

  expect_equal(mf$tests$tdt$lambda, mb$tests$tdt$lambda, tolerance = 1e-10)
  expect_equal(mf$tests$tdt$power, mb$tests$tdt$power, tolerance = 1e-10)

  mssn_mb <- tdt_mssn_conditional_full(
    power = 0.8, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.20, pd = 0.40, R1 = 1.2, R2 = 1.44,
    delta_prime = 0.5,
    verbose = FALSE
  )

  mssn_mf <- tdt_mssn_conditional_full(
    power = 0.8, alpha = 0.05,
    input_mode = "model_free",
    ET = mb$transmissions$ET,
    ENT = mb$transmissions$ENT,
    n_trios = 800,
    verbose = FALSE
  )

  expect_equal(mssn_mf$tests$tdt$MSSN_trios, mssn_mb$tests$tdt$MSSN_trios)
})


test_that("MSSN and power are mutually consistent", {
  target <- 0.8

  mssn <- tdt_mssn_conditional_full(
    power = target, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05, pd = 0.10, R1 = 1.5, R2 = 2.25,
    delta_prime = 0.8,
    verbose = FALSE
  )

  n_star <- mssn$tests$tdt$MSSN_trios

  at_mssn <- tdt_power_conditional_full(
    n_trios = n_star, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05, pd = 0.10, R1 = 1.5, R2 = 2.25,
    delta_prime = 0.8,
    verbose = FALSE
  )

  below_mssn <- tdt_power_conditional_full(
    n_trios = n_star - 1, alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05, pd = 0.10, R1 = 1.5, R2 = 2.25,
    delta_prime = 0.8,
    verbose = FALSE
  )

  expect_gte(at_mssn$tests$tdt$power, target)
  expect_lt(below_mssn$tests$tdt$power, target)
})


test_that("mode of inheritance derives R1 as in the case-control framework", {
  mult <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05, input_mode = "model_based",
    prev = 0.05, pd = 0.25, R2 = 1.5, MOI = "M", verbose = FALSE
  )
  dom <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05, input_mode = "model_based",
    prev = 0.05, pd = 0.25, R2 = 1.5, MOI = "D", verbose = FALSE
  )
  rec <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05, input_mode = "model_based",
    prev = 0.05, pd = 0.25, R2 = 1.5, MOI = "Rec", verbose = FALSE
  )

  expect_equal(mult$model_info$R1, sqrt(1.5))
  expect_equal(dom$model_info$R1, 1.5)
  expect_equal(rec$model_info$R1, 1)

  explicit <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05, input_mode = "model_based",
    prev = 0.05, pd = 0.25, R1 = 1.3, R2 = 1.5, MOI = "M", verbose = FALSE
  )
  expect_equal(explicit$model_info$R1, 1.3)
})


test_that("power increases with trios and MSSN increases with target power", {
  pow <- vapply(
    c(100, 200, 400, 800),
    function(n) {
      tdt_power_conditional_full(
        n_trios = n, alpha = 0.05, input_mode = "model_based",
        prev = 0.05, pd = 0.25, R1 = 1.5, R2 = 2.25,
        delta_prime = 1, verbose = FALSE
      )$tests$tdt$power
    },
    numeric(1)
  )
  expect_true(all(diff(pow) > 0))

  mssn <- vapply(
    c(0.5, 0.7, 0.8, 0.9, 0.95),
    function(p) {
      tdt_mssn_conditional_full(
        power = p, alpha = 0.05, input_mode = "model_based",
        prev = 0.05, pd = 0.25, R1 = 1.5, R2 = 2.25,
        delta_prime = 1, verbose = FALSE
      )$tests$tdt$MSSN_trios
    },
    numeric(1)
  )
  expect_true(all(diff(mssn) > 0))
})


test_that("locus heterogeneity is disabled by default", {
  out <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05, input_mode = "model_based",
    prev = 0.05, pd = 0.25, R1 = 1.5, R2 = 2.25, verbose = FALSE
  )

  expect_false(out$locus_het$enabled)
  expect_false(out$errors$phenotype_misclass$enabled)
  expect_false(out$errors$genotype_misclass$enabled)
  expect_equal(out$transmissions$gT_base, out$transmissions$gT_obs)
  expect_equal(out$transmissions$gNT_base, out$transmissions$gNT_obs)
})


test_that("locus_het = TRUE, pi = 1 reproduces the no-heterogeneity baseline exactly", {
  baseline <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05, input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    verbose = FALSE
  )
  het_homogeneous <- tdt_power_conditional_full(
    n_trios = 500, alpha = 0.05, input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    locus_het = TRUE, pi = 1,
    verbose = FALSE
  )

  expect_true(het_homogeneous$locus_het$enabled)
  expect_equal(het_homogeneous$transmissions$gT_obs, baseline$transmissions$gT_obs)
  expect_equal(het_homogeneous$transmissions$gNT_obs, baseline$transmissions$gNT_obs)
  expect_equal(het_homogeneous$tests$tdt$power, baseline$tests$tdt$power)

  mssn_baseline <- tdt_mssn_conditional_full(
    power = 0.8, alpha = 0.05, input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    verbose = FALSE
  )
  mssn_het_homogeneous <- tdt_mssn_conditional_full(
    power = 0.8, alpha = 0.05, input_mode = "model_based",
    prev = tdt_pars$prev, pd = tdt_pars$pd,
    R1 = tdt_pars$R1, R2 = tdt_pars$R2,
    delta_prime = tdt_pars$delta_prime,
    locus_het = TRUE, pi = 1,
    verbose = FALSE
  )

  expect_equal(
    mssn_het_homogeneous$tests$tdt$MSSN_trios,
    mssn_baseline$tests$tdt$MSSN_trios
  )
})


test_that("locus heterogeneity power matches tdt_power_full's heterogeneity branch", {
  for (pi_val in c(1, 0.99, 0.95, 0.9, 0.7)) {
    new <- tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = tdt_pars$prev, pd = tdt_pars$pd,
      R1 = tdt_pars$R1, R2 = tdt_pars$R2,
      delta_prime = tdt_pars$delta_prime,
      locus_het = TRUE, pi = pi_val,
      verbose = FALSE
    )

    old <- suppressMessages(
      tdt_power_full(
        N = 500, pd = tdt_pars$pd, prev = tdt_pars$prev,
        R1 = tdt_pars$R1, R2 = tdt_pars$R2,
        alpha = 0.05, delta_prime = tdt_pars$delta_prime,
        misclass_rate = 0, heter_rate = 1 - pi_val,
        verbose = FALSE
      )
    )

    expect_equal(
      new$tests$tdt$power, old$power$heterogeneity,
      tolerance = 1e-10, info = paste("pi =", pi_val)
    )
    expect_equal(
      new$tests$tdt$lambda, old$lambda$heterogeneity,
      tolerance = 1e-10, info = paste("pi =", pi_val)
    )
  }
})


test_that("locus heterogeneity MSSN matches ceiling(tdt_required_trios(pi = ...))", {
  for (pi_val in c(1, 0.99, 0.95, 0.9, 0.7)) {
    new <- tdt_mssn_conditional_full(
      power = 0.8, alpha = 0.05, input_mode = "model_based",
      prev = tdt_pars$prev, pd = tdt_pars$pd,
      R1 = tdt_pars$R1, R2 = tdt_pars$R2,
      delta_prime = tdt_pars$delta_prime,
      locus_het = TRUE, pi = pi_val,
      verbose = FALSE
    )

    old <- suppressMessages(
      tdt_required_trios(
        power = 0.8, alpha = 0.05, df = 1,
        pd = tdt_pars$pd, prev = tdt_pars$prev,
        R1 = tdt_pars$R1, R2 = tdt_pars$R2,
        delta_prime = tdt_pars$delta_prime, pi = pi_val
      )
    )

    expect_equal(
      new$tests$tdt$MSSN_trios,
      ceiling(old$`Required Number of Trios (N_star)`),
      info = paste("pi =", pi_val)
    )
  }
})


test_that("locus_het = TRUE requires input_mode = 'model_based'", {
  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_free",
      ET = 140, ENT = 100,
      locus_het = TRUE, pi = 0.9,
      verbose = FALSE
    ),
    "model_based"
  )

  expect_error(
    tdt_mssn_conditional_full(
      power = 0.8, alpha = 0.05, input_mode = "model_free",
      ET = 140, ENT = 100, n_trios = 120,
      locus_het = TRUE, pi = 0.9,
      verbose = FALSE
    ),
    "model_based"
  )
})


test_that("locus heterogeneity argument validation rejects malformed arguments", {
  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R2 = 1.5,
      locus_het = "yes", pi = 0.9,
      verbose = FALSE
    ),
    "locus_het"
  )

  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R2 = 1.5,
      locus_het = TRUE, pi = 1.5,
      verbose = FALSE
    ),
    "pi"
  )

  expect_error(
    tdt_mssn_conditional_full(
      power = 0.8, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R2 = 1.5,
      locus_het = TRUE, pi = -0.1,
      verbose = FALSE
    ),
    "pi"
  )
})


test_that("verbose output reports locus heterogeneity when enabled", {
  expect_message(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = tdt_pars$prev, pd = tdt_pars$pd,
      R1 = tdt_pars$R1, R2 = tdt_pars$R2,
      delta_prime = tdt_pars$delta_prime,
      locus_het = TRUE, pi = 0.9,
      verbose = TRUE
    ),
    "enabled, pi=0.900"
  )
})


test_that("a null effect gives lambda 0 and power equal to alpha", {
  expect_warning(
    out <- tdt_power_conditional_full(
      n_trios = 5000, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R1 = 1.5, R2 = 2.25,
      delta_prime = 0, verbose = FALSE
    )
  )
  expect_equal(out$tests$tdt$lambda, 0)
  expect_equal(out$tests$tdt$power, 0.05, tolerance = 1e-8)

  expect_error(
    tdt_mssn_conditional_full(
      power = 0.8, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R1 = 1.5, R2 = 2.25,
      delta_prime = 0, verbose = FALSE
    ),
    "non-centrality parameter is 0"
  )
})


test_that("input validation rejects malformed arguments", {
  expect_error(
    tdt_power_conditional_full(
      n_trios = -1, alpha = 0.05, input_mode = "model_free",
      ET = 140, ENT = 100, verbose = FALSE
    ),
    "n_trios"
  )

  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 1.2, input_mode = "model_free",
      ET = 140, ENT = 100, verbose = FALSE
    ),
    "alpha"
  )

  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_free",
      ET = 140, verbose = FALSE
    ),
    "ET and ENT"
  )

  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, verbose = FALSE
    ),
    "R2 must be supplied"
  )

  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 1.5, R2 = 1.5, verbose = FALSE
    ),
    "pd"
  )

  expect_error(
    tdt_mssn_conditional_full(
      power = 0.8, alpha = 0.05, input_mode = "model_free",
      ET = 140, ENT = 100, verbose = FALSE
    ),
    "n_trios must be supplied"
  )

  expect_error(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R2 = 1.5, delta_prime = 2, verbose = FALSE
    ),
    "delta_prime"
  )
})


test_that("verbose output is emitted as messages and returns invisibly", {
  expect_message(
    tdt_power_conditional_full(
      n_trios = 500, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R1 = 1.5, R2 = 2.25, verbose = TRUE
    ),
    "Family-Based \\(TDT\\): Power"
  )

  expect_message(
    tdt_mssn_conditional_full(
      power = 0.8, alpha = 0.05, input_mode = "model_based",
      prev = 0.05, pd = 0.25, R1 = 1.5, R2 = 2.25, verbose = TRUE
    ),
    "Minimum Sample Size Necessary"
  )
})
