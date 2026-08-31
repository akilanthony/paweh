cc_ngs_mssn_test_args <- function(...) {
  modifyList(
    list(
      power = 0.80, alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
      coverage = 20, seq_error = 0.01, MOI = "M", k = 1.2,
      verbose = FALSE
    ),
    list(...)
  )
}

cc_ngs_mssn_reference_ncp <- function(power, alpha) {
  critical <- qchisq(1 - alpha, df = 1)
  uniroot(
    function(lambda) {
      pchisq(
        critical, df = 1, ncp = lambda, lower.tail = FALSE
      ) - power
    },
    lower = 0,
    upper = 1e6
  )$root
}

cc_ngs_mssn_design_power <- function(result, N_case) {
  N_control <- ceiling(result$k * N_case)
  lambda <- .cc_ahn_trend_ncp(
    result$freqs$case_called, result$freqs$control_called,
    N_case, N_control, result$scores
  )
  list(
    N_control = N_control,
    lambda = lambda,
    power = pchisq(
      qchisq(1 - result$alpha, df = 1),
      df = 1, ncp = lambda, lower.tail = FALSE
    )
  )
}

test_that("target NCP inversion matches an independent uniroot", {
  settings <- list(
    c(power = 0.80, alpha = 0.05),
    c(power = 0.90, alpha = 0.001),
    c(power = 0.90, alpha = 5e-8),
    c(power = 0.95, alpha = 0.01)
  )

  for (setting in settings) {
    expected <- cc_ngs_mssn_reference_ncp(
      setting[["power"]], setting[["alpha"]]
    )
    observed <- .cc_ngs_target_ncp(
      setting[["power"]], setting[["alpha"]]
    )
    expect_equal(observed, expected, tolerance = 1e-10)
  }
})

test_that("continuous case MSSN follows the direct Ahn identity", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  scores <- c(0, 1, 2)
  k <- 1.4
  lambda_target <- cc_ngs_mssn_reference_ncp(0.80, 0.05)
  D <- sum(scores * (g_case - g_control))
  pooled <- g_case + k * g_control
  Q <- sum(scores^2 * pooled) -
    sum(scores * pooled)^2 / (1 + k)
  expected <- lambda_target * Q / (k * D^2)
  observed <- .cc_ngs_mssn_components(
    g_case, g_control, k, scores, lambda_target
  )

  expect_equal(observed$D, D, tolerance = 1e-15)
  expect_equal(observed$Q, Q, tolerance = 1e-15)
  expect_equal(observed$N_case_continuous, expected, tolerance = 1e-12)
})

test_that("MSSN round trips to target power with actual integer allocation", {
  settings <- list(
    cc_ngs_mssn_test_args(MOI = "M", k = 1, coverage = 4),
    cc_ngs_mssn_test_args(
      MOI = "D", pd = 0.20, R2 = 2.0, prev = 0.02,
      k = 1.5, coverage = 20, seq_error = 0.04, power = 0.90
    ),
    cc_ngs_mssn_test_args(
      MOI = "Rec", pd = 0.40, R2 = 2.5, prev = 0.10,
      k = 0.75, coverage = 50, seq_error = 0.001, alpha = 0.001
    )
  )

  for (args in settings) {
    result <- do.call(cc_ngs_mssn, args)
    achieved <- cc_ngs_mssn_design_power(result, result$MSSN_case)

    expect_equal(result$MSSN_ctrl, achieved$N_control)
    expect_equal(result$achieved_lambda, achieved$lambda, tolerance = 1e-12)
    expect_equal(result$achieved_power, achieved$power, tolerance = 1e-14)
    expect_gte(result$achieved_power + 1e-12, result$power_target)
  }
})

test_that("planned integer MSSN is minimal under the ceiling convention", {
  settings <- list(
    cc_ngs_mssn_test_args(k = 0.5, coverage = 4),
    cc_ngs_mssn_test_args(k = 1, coverage = 20),
    cc_ngs_mssn_test_args(k = 1.7, coverage = 50, MOI = "D"),
    cc_ngs_mssn_test_args(k = 2, coverage = 10, MOI = "Rec", pd = 0.4)
  )

  for (args in settings) {
    result <- do.call(cc_ngs_mssn, args)
    achieved <- cc_ngs_mssn_design_power(result, result$MSSN_case)
    expect_gte(achieved$power + 1e-12, result$power_target)

    if (result$MSSN_case > 1) {
      previous <- cc_ngs_mssn_design_power(result, result$MSSN_case - 1)
      expect_lt(previous$power, result$power_target)
    }
  }
})

test_that("high-depth MSSN converges to canonical trend MSSN", {
  for (MOI in c("M", "D", "Rec")) {
    args <- cc_ngs_mssn_test_args(
      MOI = MOI, coverage = 500, seq_error = 0.01, k = 1
    )
    sequencing <- do.call(cc_ngs_mssn, args)
    canonical <- cc_mssn(
      power = args$power, alpha = args$alpha,
      input_mode = "model_based",
      prev = args$prev, pd = args$pd, R2 = args$R2, MOI = MOI,
      k = args$k, w = .cc_ngs_scores_from_moi(MOI), verbose = FALSE
    )

    expect_equal(
      sequencing$MSSN_case,
      canonical$tests$trend$MSSN_case,
      tolerance = 1
    )
    expect_equal(
      sequencing$N_case_continuous,
      canonical$tests$trend$lambda_star /
        (args$k * canonical$tests$trend$S),
      tolerance = 1e-9
    )
  }
})

test_that("finite-depth zero-error sequencing can increase MSSN", {
  args <- cc_ngs_mssn_test_args(
    MOI = "M", coverage = 4, seq_error = 0, k = 1
  )
  sequencing <- do.call(cc_ngs_mssn, args)
  canonical <- cc_mssn(
    power = args$power, alpha = args$alpha,
    input_mode = "model_based",
    prev = args$prev, pd = args$pd, R2 = args$R2, MOI = args$MOI,
    k = args$k, w = c(0, 1, 2), verbose = FALSE
  )

  expect_gt(sequencing$MSSN_case, canonical$tests$trend$MSSN_case)
  expect_gt(
    sequencing$N_case_continuous,
    canonical$tests$trend$lambda_star /
      (args$k * canonical$tests$trend$S)
  )
})

test_that("higher coverage approaches ordinary MSSN", {
  coverages <- c(2, 4, 20, 100)
  results <- lapply(
    coverages,
    function(coverage) {
      do.call(
        cc_ngs_mssn,
        cc_ngs_mssn_test_args(coverage = coverage, seq_error = 0.01, k = 1)
      )
    }
  )
  canonical <- cc_mssn(
    power = 0.80, alpha = 0.05, input_mode = "model_based",
    prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M", k = 1,
    w = c(0, 1, 2), verbose = FALSE
  )$tests$trend
  continuous <- vapply(
    results,
    function(x) x[["N_case_continuous"]],
    numeric(1)
  )
  canonical_continuous <- canonical$lambda_star / canonical$S

  expect_lt(
    abs(continuous[[4]] - canonical_continuous),
    abs(continuous[[1]] - canonical_continuous)
  )
  expect_equal(continuous[[4]], canonical_continuous, tolerance = 1e-9)
})

test_that("control allocations follow the PAWEH ceiling convention", {
  for (k in c(0.5, 1, 2)) {
    result <- do.call(cc_ngs_mssn, cc_ngs_mssn_test_args(k = k))

    expect_gt(result$MSSN_case, 0)
    expect_equal(result$MSSN_case, ceiling(result$MSSN_case))
    expect_equal(result$MSSN_ctrl, ceiling(k * result$MSSN_case))
    expect_equal(result$MSSN_total, result$MSSN_case + result$MSSN_ctrl)
    expect_true(is.finite(result$achieved_power))
    expect_true(result$achieved_power >= 0 && result$achieved_power <= 1)
  }
})

test_that("MSSN reuses the S5 MOI score mapping", {
  mappings <- list(
    M = c(0, 1, 2),
    D = c(0, 1, 1),
    Rec = c(0, 0, 1)
  )
  for (MOI in names(mappings)) {
    result <- do.call(cc_ngs_mssn, cc_ngs_mssn_test_args(MOI = MOI))
    expect_identical(result$scores, mappings[[MOI]])
  }
})

test_that("cc_ngs_mssn returns a transparent invisible result", {
  visible <- withVisible(do.call(cc_ngs_mssn, cc_ngs_mssn_test_args()))
  result <- visible$value
  expected_fields <- c(
    "power_target", "alpha", "MSSN_case", "MSSN_ctrl", "MSSN_total",
    "N_case_continuous", "achieved_power", "achieved_lambda",
    "lambda_target", "initial_MSSN_case", "rounding_adjustment", "k",
    "MOI", "scores", "coverage", "seq_error", "model_info", "freqs",
    "transition_matrix"
  )

  expect_false(visible$visible)
  expect_identical(class(result), "cc_ngs_mssn")
  expect_true(all(expected_fields %in% names(result)))
  expect_true(is.finite(result$lambda_target))
  expect_true(is.finite(result$achieved_lambda))
  expect_gte(result$lambda_target, 0)
  expect_gte(result$achieved_lambda, 0)
  expect_identical(dim(result$transition_matrix), c(3L, 3L))
  expect_equal(
    unname(rowSums(result$transition_matrix)),
    rep(1, 3),
    tolerance = 1e-14
  )
  expect_true(is.list(result$model_info))
  expect_true(all(c("penetrances", "R1", "R2", "MOI") %in%
                    names(result$model_info)))

  for (g in result$freqs) {
    expect_true(all(is.finite(g)))
    expect_true(all(g >= 0 & g <= 1))
    expect_equal(sum(g), 1, tolerance = 1e-14)
  }
})

test_that("MSSN verbose output is detailed and optional", {
  expect_output(
    invisible(do.call(cc_ngs_mssn, cc_ngs_mssn_test_args(verbose = FALSE))),
    NA
  )
  output <- capture.output(
    invisible(do.call(cc_ngs_mssn, cc_ngs_mssn_test_args(verbose = TRUE))),
    type = "message"
  )

  expect_true(any(grepl("Minimum Sample Size Necessary", output, fixed = TRUE)))
  expect_true(any(grepl("Required Sample Size", output, fixed = TRUE)))
  expect_true(any(grepl("Total MSSN:", output, fixed = TRUE)))
  expect_true(any(grepl("Achieved power:", output, fixed = TRUE)))
  expect_false(any(grepl("transition_matrix", output, fixed = TRUE)))
})

test_that("cc_ngs_mssn validates public inputs", {
  invalid_values <- list(
    power = list(0, 1, -0.1, NA_real_, NaN, Inf, c(0.8, 0.9)),
    alpha = list(0, 1, -0.1, NA_real_, NaN, Inf, c(0.01, 0.05)),
    pd = list(0, 1, -0.1, NA_real_, NaN, Inf, c(0.2, 0.3)),
    R2 = list(0, -1, NA_real_, NaN, Inf, c(1, 2)),
    prev = list(0, 1, -0.1, NA_real_, NaN, Inf, c(0.1, 0.2)),
    k = list(0, -1, NA_real_, NaN, Inf, c(1, 2)),
    coverage = list(0, -1, 1.5, NA_real_, NaN, Inf, c(10, 20)),
    seq_error = list(-0.1, 0.5, 1, NA_real_, NaN, Inf, c(0.01, 0.02))
  )

  for (argument in names(invalid_values)) {
    for (value in invalid_values[[argument]]) {
      args <- cc_ngs_mssn_test_args()
      args[[argument]] <- value
      expect_error(do.call(cc_ngs_mssn, args), argument)
    }
  }

  expect_error(
    do.call(cc_ngs_mssn, cc_ngs_mssn_test_args(MOI = "additive")),
    "arg"
  )
  expect_error(
    do.call(cc_ngs_mssn, cc_ngs_mssn_test_args(verbose = NA)),
    "verbose"
  )
})

test_that("zero trend contrast gives an informative finite-MSSN error", {
  equal_frequencies <- c(0.25, 0.50, 0.25)
  expect_error(
    .cc_ngs_mssn_components(
      equal_frequencies, equal_frequencies, k = 1,
      scores = c(0, 1, 2), lambda_target = 7.8488
    ),
    "No finite MSSN exists because the trend contrast is zero"
  )

  expect_error(
    do.call(cc_ngs_mssn, cc_ngs_mssn_test_args(R2 = 1)),
    "No finite MSSN exists because the trend contrast is zero"
  )
})

test_that("CC-NGS MSSN API excludes deferred sequencing features", {
  arguments <- names(formals(cc_ngs_mssn))
  excluded <- c(
    "N_case", "pheno_misclass", "theta", "phi", "case_seq_error",
    "ctrl_seq_error", "coverage_dist", "g1", "g0", "raw_reads", "trios",
    "plot", "shiny"
  )

  expect_false(any(excluded %in% arguments))
})
