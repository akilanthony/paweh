test_that("tdt_power_full matches tdt_power_from_model in no-error case", {
  pd   <- 0.30
  prev <- 0.05
  R1   <- 1.5
  R2   <- 2.25
  N    <- 600
  a    <- 0.05
  dp   <- 1

  full <- tdt_power_full(
    N             = N,
    pd            = pd,
    prev          = prev,
    R1            = R1,
    R2            = R2,
    alpha         = a,
    delta_prime   = dp,
    misclass_rate = 0,
    heter_rate    = 0,
    verbose       = FALSE
  )

  base <- tdt_power_from_model(
    pd          = pd,
    N           = N,
    delta_prime = dp,
    prev        = prev,
    R1          = R1,
    R2          = R2,
    alpha       = a
  )

  expect_equal(
    full$lambda$no_error,
    base[["Non-Centrality Parameter (lambda)"]],
    tolerance = 1e-8
  )
  expect_equal(
    full$power$no_error,
    base[["Power"]],
    tolerance = 1e-8
  )
  expect_equal(
    full$ET$no_error,
    base[["Expected Transmissions (ET)"]],
    tolerance = 1e-8
  )
  expect_equal(
    full$ENT$no_error,
    base[["Expected Non-Transmissions (ENT)"]],
    tolerance = 1e-8
  )
})

test_that("misclassification gT* / gNT* agree with tdt_required_trios_misclass", {
  pd   <- 0.30
  prev <- 0.05
  R1   <- 1.5
  R2   <- 2.25
  N    <- 600
  a    <- 0.05
  dp   <- 1
  pi01 <- 0.10

  full <- tdt_power_full(
    N             = N,
    pd            = pd,
    prev          = prev,
    R1            = R1,
    R2            = R2,
    alpha         = a,
    delta_prime   = dp,
    misclass_rate = pi01,
    heter_rate    = 0,
    verbose       = FALSE
  )

  # gT* and gNT* do not depend on lambda_star, so any positive value is fine.
  misc <- tdt_required_trios_misclass(
    lambda_star  = 1,
    pd           = pd,
    prev         = prev,
    R1           = R1,
    R2           = R2,
    delta_prime  = dp,
    pi01         = pi01
  )

  expect_equal(full$gT_star$misclassification,  misc$gT_star,  tolerance = 1e-8)
  expect_equal(full$gNT_star$misclassification, misc$gNT_star, tolerance = 1e-8)

  # ET/ENT should be 2N * gT* / gNT*
  expect_equal(full$ET$misclassification,  2 * N * misc$gT_star,  tolerance = 1e-8)
  expect_equal(full$ENT$misclassification, 2 * N * misc$gNT_star, tolerance = 1e-8)
})

test_that("heterogeneity gT* / gNT* agree with tdt_required_trios", {
  pd   <- 0.30
  prev <- 0.05
  R1   <- 1.5
  R2   <- 2.25
  N    <- 600
  a    <- 0.05
  dp   <- 1
  het  <- 0.20              # 20% of trios not due to the locus
  pi   <- 1 - het

  full <- tdt_power_full(
    N             = N,
    pd            = pd,
    prev          = prev,
    R1            = R1,
    R2            = R2,
    alpha         = a,
    delta_prime   = dp,
    misclass_rate = 0,
    heter_rate    = het,
    verbose       = FALSE
  )

  baseN <- tdt_required_trios(
    power       = 0.80,     # any valid power is fine for gT*/gNT*
    alpha       = a,
    df          = 1,
    pd          = pd,
    prev        = prev,
    R1          = R1,
    R2          = R2,
    delta_prime = dp,
    pi          = pi
  )

  expect_equal(
    full$gT_star$heterogeneity,
    baseN[["Expected Transmission (gT_star)"]],
    tolerance = 1e-8
  )
  expect_equal(
    full$gNT_star$heterogeneity,
    baseN[["Expected Non-Transmission (gNT_star)"]],
    tolerance = 1e-8
  )
})

test_that("tdt_required_trios_full matches tdt_required_trios in no-error case", {
  pd   <- 0.30
  prev <- 0.05
  R1   <- 1.5
  R2   <- 2.25
  a    <- 0.05
  dp   <- 1
  pow  <- 0.80

  fullN <- tdt_required_trios_full(
    target_power  = pow,
    pd            = pd,
    prev          = prev,
    R1            = R1,
    R2            = R2,
    alpha         = a,
    delta_prime   = dp,
    misclass_rate = 0,
    heter_rate    = 0,
    verbose       = FALSE
  )

  baseN <- tdt_required_trios(
    power       = pow,
    alpha       = a,
    df          = 1,
    pd          = pd,
    prev        = prev,
    R1          = R1,
    R2          = R2,
    delta_prime = dp,
    pi          = 1
  )

  expect_equal(
    fullN$N$no_error,
    baseN[["Required Number of Trios (N_star)"]],
    tolerance = 1e-6
  )
})

test_that("penetrance checks / warnings behave as expected", {
  # A parameter set that forces f1 or f2 > 1
  pd   <- 0.50
  prev <- 0.30
  R1   <- 10
  R2   <- 20

  expect_warning(
    tdt_power_full(
      N       = 100,
      pd      = pd,
      prev    = prev,
      R1      = R1,
      R2      = R2,
      verbose = FALSE
    ),
    regexp = "penetrances f1 or f2 exceed 1",
    fixed  = TRUE
  )

  expect_error(
    tdt_required_trios_full(
      target_power = 0.80,
      pd           = pd,
      prev         = prev,
      R1           = R1,
      R2           = R2,
      verbose      = FALSE
    ),
    regexp = "Invalid penetrances",
    fixed  = TRUE
  )
})

