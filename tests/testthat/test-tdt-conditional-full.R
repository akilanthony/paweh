# Tests for the conditional TDT framework (R/tdt_conditional.R)
#
# tdt_power_conditional_full() is a copy of tdt_power_full(), and
# tdt_mssn_conditional_full() is a copy of tdt_required_trios_full(), both from
# R/tdt_functions.R, with input_mode = "model_based" (default, numerically
# identical to the copied function) or "model_free" (ET/ENT supplied
# directly) added on top. These tests pin model_based against the copied
# functions, and confirm model_free round-trips exactly against model_based
# when fed that same run's own ET/ENT/pd/prev.

anchor <- list(
  pd = 0.30,
  prev = 0.05,
  R1 = 1.5,
  R2 = 2.25,
  delta_prime = 1
)


# ---- model_based numerically identical to tdt_power_full() -----------------

test_that("model_based power matches tdt_power_full() at default rates", {
  new <- tdt_power_conditional_full(
    N = 600, input_mode = "model_based",
    pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
    verbose = FALSE
  )
  old <- suppressMessages(
    tdt_power_full(
      N = 600, pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
      verbose = FALSE
    )
  )

  expect_equal(new$power, old$power)
  expect_equal(new$lambda, old$lambda)
  expect_equal(new$gT_star, old$gT_star)
  expect_equal(new$gNT_star, old$gNT_star)
  expect_equal(new$ET, old$ET)
  expect_equal(new$ENT, old$ENT)
  expect_equal(new$power_loss, old$power_loss)
})


test_that("model_based power matches tdt_power_full() with custom heter_rate and misclass_rate", {
  for (heter_rate in c(0, 0.05, 0.2)) {
    for (misclass_rate in c(0, 0.02, 0.15)) {
      new <- tdt_power_conditional_full(
        N = 800, input_mode = "model_based",
        pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
        delta_prime = anchor$delta_prime,
        misclass_rate = misclass_rate, heter_rate = heter_rate,
        verbose = FALSE
      )
      old <- suppressMessages(
        tdt_power_full(
          N = 800, pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
          delta_prime = anchor$delta_prime,
          misclass_rate = misclass_rate, heter_rate = heter_rate,
          verbose = FALSE
        )
      )
      info <- paste("heter_rate =", heter_rate, "misclass_rate =", misclass_rate)
      expect_equal(new$power, old$power, info = info)
      expect_equal(new$lambda, old$lambda, info = info)
      expect_equal(new$gT_star, old$gT_star, info = info)
      expect_equal(new$gNT_star, old$gNT_star, info = info)
    }
  }
})


test_that("model_based class and structure match the copied function's layout", {
  out <- tdt_power_conditional_full(
    N = 600, input_mode = "model_based",
    pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
    verbose = FALSE
  )
  expect_s3_class(out, "tdt_power_conditional_full")
  expect_equal(out$input_mode, "model_based")
  expect_true(all(c("alpha", "N", "lambda", "power", "power_loss", "gT_star",
                     "gNT_star", "ET", "ENT", "model_parameters") %in% names(out)))
  expect_named(out$power, c("no_error", "misclassification", "heterogeneity"))
})


# ---- model_based numerically identical to tdt_required_trios_full() --------

test_that("model_based MSSN matches tdt_required_trios_full() at default rates", {
  new <- tdt_mssn_conditional_full(
    target_power = 0.80, input_mode = "model_based",
    pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
    verbose = FALSE
  )
  old <- suppressMessages(
    tdt_required_trios_full(
      target_power = 0.80, pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
      verbose = FALSE
    )
  )

  expect_equal(new$N, old$N)
  expect_equal(new$lambda_star, old$lambda_star)
  expect_equal(new$percent_increase, old$percent_increase)
  expect_equal(new$power_at_N_no_error, old$power_at_N_no_error)
  expect_equal(new$power_loss_at_N_no_error, old$power_loss_at_N_no_error)
  expect_equal(new$gT_star, old$gT_star)
  expect_equal(new$gNT_star, old$gNT_star)
})


test_that("model_based MSSN matches tdt_required_trios_full() with custom heter_rate and misclass_rate", {
  for (heter_rate in c(0, 0.05, 0.2)) {
    for (misclass_rate in c(0, 0.02, 0.15)) {
      new <- tdt_mssn_conditional_full(
        target_power = 0.80, input_mode = "model_based",
        pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
        delta_prime = anchor$delta_prime,
        misclass_rate = misclass_rate, heter_rate = heter_rate,
        verbose = FALSE
      )
      old <- suppressMessages(
        tdt_required_trios_full(
          target_power = 0.80, pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
          delta_prime = anchor$delta_prime,
          misclass_rate = misclass_rate, heter_rate = heter_rate,
          verbose = FALSE
        )
      )
      info <- paste("heter_rate =", heter_rate, "misclass_rate =", misclass_rate)
      expect_equal(new$N, old$N, info = info)
      expect_equal(new$lambda_star, old$lambda_star, info = info)
      expect_equal(new$gT_star, old$gT_star, info = info)
      expect_equal(new$gNT_star, old$gNT_star, info = info)
    }
  }
})


test_that("model_based MSSN class and structure match the copied function's layout", {
  out <- tdt_mssn_conditional_full(
    target_power = 0.80, input_mode = "model_based",
    pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
    verbose = FALSE
  )
  expect_s3_class(out, "tdt_mssn_conditional_full")
  expect_equal(out$input_mode, "model_based")
  expect_true(all(c("alpha", "target_power", "lambda_star", "N", "percent_increase",
                     "power_at_N_no_error", "power_loss_at_N_no_error", "gT_star",
                     "gNT_star", "model_parameters") %in% names(out)))
  expect_named(out$N, c("no_error", "misclassification", "heterogeneity"))
})


# ---- model_free round-trips against model_based -----------------------------

test_that("model_free power round-trips exactly against model_based", {
  for (heter_rate in c(0, 0.03, 0.1)) {
    for (misclass_rate in c(0, 0.01, 0.1)) {
      mb <- tdt_power_conditional_full(
        N = 700, input_mode = "model_based",
        pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
        delta_prime = anchor$delta_prime,
        misclass_rate = misclass_rate, heter_rate = heter_rate,
        verbose = FALSE
      )

      mf <- tdt_power_conditional_full(
        N = 700, input_mode = "model_free",
        ET = mb$ET$no_error, ENT = mb$ENT$no_error,
        pd = anchor$pd, prev = anchor$prev,
        misclass_rate = misclass_rate, heter_rate = heter_rate,
        verbose = FALSE
      )

      info <- paste("heter_rate =", heter_rate, "misclass_rate =", misclass_rate)
      expect_equal(mf$power$no_error, mb$power$no_error, tolerance = 1e-8, info = info)
      expect_equal(mf$power$misclassification, mb$power$misclassification, tolerance = 1e-8, info = info)
      expect_equal(mf$power$heterogeneity, mb$power$heterogeneity, tolerance = 1e-8, info = info)
      expect_equal(mf$lambda$no_error, mb$lambda$no_error, tolerance = 1e-6, info = info)
      expect_equal(mf$lambda$misclassification, mb$lambda$misclassification, tolerance = 1e-6, info = info)
      expect_equal(mf$lambda$heterogeneity, mb$lambda$heterogeneity, tolerance = 1e-6, info = info)
    }
  }
})


test_that("model_free MSSN round-trips exactly against model_based", {
  for (heter_rate in c(0, 0.03, 0.1)) {
    for (misclass_rate in c(0, 0.01, 0.1)) {
      mb <- tdt_mssn_conditional_full(
        target_power = 0.80, input_mode = "model_based",
        pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
        delta_prime = anchor$delta_prime,
        misclass_rate = misclass_rate, heter_rate = heter_rate,
        verbose = FALSE
      )

      n_trios_used <- mb$N$no_error
      ET_val  <- mb$gT_star$no_error  * 2 * n_trios_used
      ENT_val <- mb$gNT_star$no_error * 2 * n_trios_used

      mf <- tdt_mssn_conditional_full(
        target_power = 0.80, input_mode = "model_free",
        ET = ET_val, ENT = ENT_val, n_trios = n_trios_used,
        pd = anchor$pd, prev = anchor$prev,
        misclass_rate = misclass_rate, heter_rate = heter_rate,
        verbose = FALSE
      )

      info <- paste("heter_rate =", heter_rate, "misclass_rate =", misclass_rate)
      expect_equal(mf$N$no_error, mb$N$no_error, tolerance = 1e-6, info = info)
      expect_equal(mf$N$misclassification, mb$N$misclassification, tolerance = 1e-6, info = info)
      expect_equal(mf$N$heterogeneity, mb$N$heterogeneity, tolerance = 1e-6, info = info)
      expect_equal(mf$lambda_star, mb$lambda_star, tolerance = 1e-8, info = info)
    }
  }
})


# ---- model_free: rate = 0 needs neither pd nor prev -------------------------

test_that("model_free with heter_rate = 0 and misclass_rate = 0 needs no pd or prev", {
  expect_no_error(
    out <- tdt_power_conditional_full(
      N = 500, input_mode = "model_free",
      ET = 140, ENT = 100,
      misclass_rate = 0, heter_rate = 0,
      verbose = FALSE
    )
  )
  expect_equal(out$power$no_error, out$power$misclassification)
  expect_equal(out$power$no_error, out$power$heterogeneity)
  expect_null(out$model_parameters$pd)

  expect_no_error(
    tdt_mssn_conditional_full(
      target_power = 0.8, input_mode = "model_free",
      ET = 140, ENT = 100, n_trios = 120,
      misclass_rate = 0, heter_rate = 0,
      verbose = FALSE
    )
  )
})


# ---- model_free: pd solved from ET/ENT when omitted -------------------------

test_that("model_free solves pd from ET/ENT when pd is omitted and a rate is non-zero", {
  mb <- tdt_power_conditional_full(
    N = 700, input_mode = "model_based",
    pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
    heter_rate = 0.1, misclass_rate = 0,
    verbose = FALSE
  )

  expect_message(
    mf <- tdt_power_conditional_full(
      N = 700, input_mode = "model_free",
      ET = mb$ET$no_error, ENT = mb$ENT$no_error,
      heter_rate = 0.1, misclass_rate = 0,
      verbose = FALSE
    ),
    "pd not supplied"
  )

  expect_equal(mf$model_parameters$pd, anchor$pd, tolerance = 1e-6)
  expect_equal(mf$power$heterogeneity, mb$power$heterogeneity, tolerance = 1e-6)
})


# ---- input validation -------------------------------------------------------

test_that("input validation rejects malformed arguments", {
  expect_error(
    tdt_power_conditional_full(
      N = 500, input_mode = "model_based",
      prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
      verbose = FALSE
    ),
    "pd, prev, R1, and R2"
  )

  expect_error(
    tdt_power_conditional_full(
      N = 500, input_mode = "model_free",
      ET = 140,
      verbose = FALSE
    ),
    "ET and ENT"
  )

  expect_error(
    tdt_power_conditional_full(
      N = 500, input_mode = "model_free",
      ET = 140, ENT = 100, pd = 0.3,
      misclass_rate = 0.05,
      verbose = FALSE
    ),
    "prev must be supplied"
  )

  expect_error(
    tdt_mssn_conditional_full(
      target_power = 0.8, input_mode = "model_free",
      ET = 140, ENT = 100,
      verbose = FALSE
    ),
    "n_trios must be supplied"
  )

  expect_error(
    tdt_power_conditional_full(
      N = 500, input_mode = "model_based",
      pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
      misclass_rate = 1,
      verbose = FALSE
    ),
    "misclass_rate"
  )

  expect_error(
    tdt_power_conditional_full(
      N = 500, input_mode = "model_based",
      pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
      heter_rate = -0.1,
      verbose = FALSE
    ),
    "heter_rate"
  )

  # gT == gNT == 0.5 (A = 0, zero effect size) gives a negative discriminant:
  # there is no pd consistent with a null-effect ET/ENT pair.
  expect_error(
    tdt_power_conditional_full(
      N = 500, input_mode = "model_free",
      ET = 500, ENT = 500,
      heter_rate = 0.05, misclass_rate = 0,
      verbose = FALSE
    ),
    "Cannot solve for pd"
  )
})


# ---- verbose output ----------------------------------------------------------

test_that("verbose output is emitted as messages and returns invisibly", {
  expect_message(
    tdt_power_conditional_full(
      N = 600, input_mode = "model_based",
      pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
      verbose = TRUE
    ),
    "POWER FOR A FIXED SAMPLE SIZE"
  )

  expect_message(
    tdt_mssn_conditional_full(
      target_power = 0.8, input_mode = "model_based",
      pd = anchor$pd, prev = anchor$prev, R1 = anchor$R1, R2 = anchor$R2,
      verbose = TRUE
    ),
    "MINIMUM SAMPLE SIZE NECESSARY FOR A FIXED POWER"
  )
})
