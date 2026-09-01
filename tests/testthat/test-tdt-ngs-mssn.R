tdt_ngs_mssn_args <- function(...) {
  utils::modifyList(
    list(
      power = 0.80, pd = 0.325, R1 = 1.2,
      coverage = 12, seq_error = 0.005,
      alpha = 0.05, verbose = FALSE
    ),
    list(...)
  )
}

tdt_ngs_mssn_quiet <- function(...) {
  do.call(tdt_ngs_mssn, tdt_ngs_mssn_args(...))
}

tdt_ngs_reference_target_ncp <- function(power, alpha) {
  if (power <= alpha) {
    return(0)
  }
  critical <- stats::qchisq(1 - alpha, df = 1)
  stats::uniroot(
    function(lambda) {
      stats::pchisq(
        critical, df = 1, ncp = lambda, lower.tail = FALSE
      ) - power
    },
    lower = 0,
    upper = 1e6
  )$root
}

test_that("target NCP inversion matches an independent uniroot", {
  settings <- list(
    c(power = 0.80, alpha = 0.05),
    c(power = 0.90, alpha = 0.001),
    c(power = 0.80, alpha = 5e-8)
  )
  for (setting in settings) {
    expected <- tdt_ngs_reference_target_ncp(
      setting[["power"]], setting[["alpha"]]
    )
    observed <- .tdt_ngs_target_ncp(
      setting[["power"]], setting[["alpha"]]
    )
    expect_equal(observed, expected, tolerance = 1e-11)
  }
  expect_identical(.tdt_ngs_target_ncp(0.01, 0.05), 0)
})

test_that("continuous trio MSSN follows the exact analytic identity", {
  args <- tdt_ngs_mssn_args(alpha = 5e-8)
  information <- .tdt_ngs_information(
    args$pd, args$coverage, args$seq_error
  )
  lambda_target <- tdt_ngs_reference_target_ncp(args$power, args$alpha)
  expected <- lambda_target /
    (log(args$R1)^2 * information$efficient_information)
  observed <- do.call(tdt_ngs_mssn, args)

  expect_equal(observed$lambda_target, lambda_target, tolerance = 1e-11)
  expect_equal(observed$delta, log(args$R1), tolerance = 2e-15)
  expect_equal(
    observed$efficient_information,
    information$efficient_information,
    tolerance = 1e-14
  )
  expect_equal(observed$N_trios_continuous, expected, tolerance = 1e-10)
  expect_equal(
    observed$lambda_target / observed$N_trios_continuous,
    observed$delta^2 * observed$efficient_information,
    tolerance = 1e-13
  )
})

test_that("MSSN performs one expensive information evaluation", {
  original <- .tdt_ngs_information
  evaluations <- 0L
  testthat::local_mocked_bindings(
    .tdt_ngs_information = function(...) {
      evaluations <<- evaluations + 1L
      original(...)
    },
    .package = "paweh"
  )

  result <- tdt_ngs_mssn_quiet(coverage = 4)
  expect_s3_class(result, "tdt_ngs_mssn")
  expect_identical(evaluations, 1L)
})

test_that("integer MSSN round trips through public TDT1-NGS power", {
  designs <- list(
    tdt_ngs_mssn_args(pd = 0.15, R1 = 1.1, coverage = 4),
    tdt_ngs_mssn_args(
      power = 0.90, pd = 0.325, R1 = 1.2, coverage = 12,
      seq_error = 0.01, alpha = 0.001
    ),
    tdt_ngs_mssn_args(
      pd = 0.5, R1 = 1.3, coverage = 44,
      seq_error = 0.005, alpha = 5e-8
    ),
    tdt_ngs_mssn_args(
      power = 0.90, pd = 0.675, R1 = 0.9, coverage = 20,
      seq_error = 0
    )
  )

  for (args in designs) {
    mssn <- do.call(tdt_ngs_mssn, args)
    achieved <- tdt_ngs_power(
      N = mssn$MSSN_trios,
      pd = args$pd,
      R1 = args$R1,
      coverage = args$coverage,
      seq_error = args$seq_error,
      alpha = args$alpha,
      verbose = FALSE
    )
    expect_equal(mssn$achieved_lambda, achieved$lambda, tolerance = 1e-12)
    expect_equal(mssn$achieved_power, achieved$power, tolerance = 1e-14)
    expect_gte(mssn$achieved_power, args$power)
  }
})

test_that("the immediately smaller integer trio design misses target", {
  designs <- list(
    tdt_ngs_mssn_args(pd = 0.15, R1 = 1.1, coverage = 4),
    tdt_ngs_mssn_args(pd = 0.325, R1 = 1.2, coverage = 12),
    tdt_ngs_mssn_args(
      power = 0.90, pd = 0.5, R1 = 1.3,
      coverage = 20, alpha = 5e-8
    )
  )
  for (args in designs) {
    result <- do.call(tdt_ngs_mssn, args)
    expect_gte(result$achieved_power, args$power)
    expect_equal(result$MSSN_trios, ceiling(result$N_trios_continuous))
    expect_identical(result$rounding_adjustment, 0)

    if (result$MSSN_trios > 1) {
      previous <- tdt_ngs_power(
        N = result$MSSN_trios - 1,
        pd = args$pd,
        R1 = args$R1,
        coverage = args$coverage,
        seq_error = args$seq_error,
        alpha = args$alpha,
        verbose = FALSE
      )
      expect_lt(previous$power, args$power)
    }
  }
})

test_that("Chapter 5-style MSSN settings are finite and minimal", {
  designs <- list(
    tdt_ngs_mssn_args(
      pd = 0.15, R1 = 1.1, coverage = 4,
      seq_error = 0.005, alpha = 5e-8
    ),
    tdt_ngs_mssn_args(
      power = 0.90, pd = 0.325, R1 = 1.2, coverage = 12,
      seq_error = 0.01, alpha = 5e-8
    ),
    tdt_ngs_mssn_args(
      pd = 0.5, R1 = 1.3, coverage = 20,
      seq_error = 0.005, alpha = 5e-8
    ),
    tdt_ngs_mssn_args(
      power = 0.90, pd = 0.325, R1 = 1.2, coverage = 44,
      seq_error = 0.005, alpha = 5e-8
    )
  )
  expected_MSSN <- c(46616, 6361, 2304, 6217)
  expected_continuous <- c(
    46615.5909640973,
    6360.22638377445,
    2303.34409038314,
    6216.36032773523
  )
  for (i in seq_along(designs)) {
    args <- designs[[i]]
    result <- do.call(tdt_ngs_mssn, args)
    expect_true(is.finite(result$N_trios_continuous))
    expect_gt(result$MSSN_trios, 0)
    expect_gte(result$achieved_power, args$power)
    expect_identical(result$rounding_adjustment, 0)
    expect_equal(result$MSSN_trios, expected_MSSN[[i]])
    expect_equal(
      result$N_trios_continuous,
      expected_continuous[[i]],
      tolerance = 2e-11
    )
  }
})

test_that("zero effect has no finite MSSN above alpha", {
  expect_error(
    tdt_ngs_mssn_quiet(R1 = 1, power = 0.80, alpha = 0.05),
    "No finite MSSN exists because R1 = 1 implies zero transmission effect",
    fixed = TRUE
  )

  supported_null <- tdt_ngs_mssn_quiet(R1 = 1, power = 0.01, alpha = 0.05)
  expect_identical(supported_null$lambda_target, 0)
  expect_identical(supported_null$MSSN_trios, 1)
  expect_identical(supported_null$achieved_lambda, 0)
  expect_lt(abs(supported_null$achieved_power - 0.05), 1e-14)
})

test_that("near-null effects require increasingly many trios", {
  risks <- c(1.01, 1.02, 1.05)
  results <- lapply(risks, function(R1) {
    tdt_ngs_mssn_quiet(R1 = R1, coverage = 12)
  })
  requirements <- vapply(results, `[[`, numeric(1), "MSSN_trios")
  expect_true(all(is.finite(requirements)))
  expect_true(requirements[[1]] > requirements[[2]])
  expect_true(requirements[[2]] > requirements[[3]])
})

test_that("coverage approaches a stable high-depth MSSN plateau", {
  coverage <- c(2, 4, 12, 20, 44)
  results <- lapply(coverage, function(value) {
    tdt_ngs_mssn_quiet(coverage = value, alpha = 5e-8)
  })
  continuous <- vapply(
    results, `[[`, numeric(1), "N_trios_continuous"
  )
  expect_true(all(is.finite(continuous)))
  expect_gt(continuous[[1]], continuous[[5]])
  expect_lt(
    abs(continuous[[5]] - continuous[[4]]),
    abs(continuous[[2]] - continuous[[1]])
  )
})

test_that("supported sequencing errors yield finite MSSN designs", {
  results <- lapply(c(0, 0.005, 0.01), function(error) {
    tdt_ngs_mssn_quiet(seq_error = error, alpha = 5e-8)
  })
  for (result in results) {
    expect_true(is.finite(result$N_trios_continuous))
    expect_gt(result$MSSN_trios, 0)
    expect_gte(result$achieved_power, result$power_target)
  }
})

test_that("MSSN respects disease-allele reflection symmetry", {
  for (pd in c(0.15, 0.325)) {
    left <- tdt_ngs_mssn_quiet(pd = pd, coverage = 4)
    right <- tdt_ngs_mssn_quiet(pd = 1 - pd, coverage = 4)
    expect_equal(
      left$N_trios_continuous,
      right$N_trios_continuous,
      tolerance = 2e-11
    )
    expect_equal(left$MSSN_trios, right$MSSN_trios)
    expect_identical(left$pd, pd)
    expect_identical(right$pd, 1 - pd)
  }
})

test_that("more stringent alpha requires no fewer trios", {
  requirements <- vapply(c(0.05, 0.001, 5e-8), function(alpha) {
    tdt_ngs_mssn_quiet(alpha = alpha)$MSSN_trios
  }, numeric(1))
  expect_true(requirements[[1]] <= requirements[[2]])
  expect_true(requirements[[2]] <= requirements[[3]])
})

test_that("higher target power requires no fewer trios", {
  requirements <- vapply(c(0.80, 0.90, 0.95), function(power) {
    tdt_ngs_mssn_quiet(power = power)$MSSN_trios
  }, numeric(1))
  expect_true(requirements[[1]] <= requirements[[2]])
  expect_true(requirements[[2]] <= requirements[[3]])
})

test_that("MSSN result is transparent, valid, and returned invisibly", {
  visible <- withVisible(tdt_ngs_mssn_quiet(alpha = 5e-8))
  result <- visible$value

  expect_false(visible$visible)
  expect_s3_class(result, "tdt_ngs_mssn", exact = TRUE)
  expect_named(
    result,
    c(
      "power_target", "alpha", "MSSN_trios", "total_individuals",
      "N_trios_continuous", "achieved_power", "achieved_lambda",
      "lambda_target", "ncp_per_trio", "initial_MSSN_trios",
      "rounding_adjustment", "pd", "R1", "R2", "t", "delta",
      "coverage", "seq_error", "efficient_information",
      "information_matrix", "nuisance_rcond", "score_mean", "model_info"
    )
  )
  expect_gt(result$MSSN_trios, 0)
  expect_equal(result$MSSN_trios, ceiling(result$MSSN_trios))
  expect_identical(result$total_individuals, 3 * result$MSSN_trios)
  expect_true(is.finite(result$N_trios_continuous))
  expect_true(result$achieved_power >= 0 && result$achieved_power <= 1)
  expect_gte(result$achieved_power, result$power_target)
  expect_true(all(is.finite(c(
    result$achieved_lambda, result$lambda_target, result$ncp_per_trio
  ))))
  expect_gt(result$efficient_information, 0)
  expect_equal(dim(result$information_matrix), c(11, 11))
  expect_true(all(is.finite(result$information_matrix)))
  expect_equal(result$information_matrix, t(result$information_matrix),
               tolerance = 1e-14)
  expect_equal(result$R2, result$R1^2, tolerance = 1e-15)
  expect_equal(result$t, result$R1 / (1 + result$R1), tolerance = 1e-15)
  expect_equal(result$delta, log(result$R1), tolerance = 2e-15)
  expect_identical(result$model_info$test, "TDT1-NGS")
  expect_identical(result$model_info$objective, "MSSN")
  expect_identical(result$model_info$sampling_unit, "complete_trios")
  expect_identical(result$model_info$sample_size_solution, "analytic")
  expect_identical(
    result$model_info$likelihood,
    "raw_read_counts_with_latent_trio_states"
  )
})

test_that("MSSN verbose report is detailed while print remains concise", {
  expect_output(tdt_ngs_mssn_quiet(verbose = FALSE), NA)
  verbose_text <- capture.output(
    tdt_ngs_mssn_quiet(verbose = TRUE),
    type = "message"
  )
  expect_true(any(grepl("Minimum Sample Size Necessary", verbose_text,
                        fixed = TRUE)))
  expect_true(any(grepl("Required Sample Size", verbose_text, fixed = TRUE)))
  expect_true(any(grepl("Efficient information", verbose_text, fixed = TRUE)))
  result <- tdt_ngs_mssn_quiet()
  printed <- capture.output(print(result))
  expect_lte(length(printed), 7)
  expect_true(any(grepl("Required trios", printed, fixed = TRUE)))
  expect_false(any(grepl("information_matrix", printed, fixed = TRUE)))
})

test_that("coverage one uses the same public identifiability error", {
  expect_error(
    tdt_ngs_mssn_quiet(coverage = 1),
    "coverage = 1 is unsupported.*not identifiable"
  )
})

test_that("MSSN validation rejects unsupported inputs", {
  invalid <- list(
    list(power = 0), list(power = -0.1), list(power = 1),
    list(power = NA_real_), list(power = Inf),
    list(pd = 0), list(pd = -0.1), list(pd = 1), list(pd = 1.1),
    list(pd = NA_real_), list(pd = Inf),
    list(R1 = 0), list(R1 = -1), list(R1 = NA_real_), list(R1 = Inf),
    list(coverage = 0), list(coverage = -1), list(coverage = 1.5),
    list(coverage = NA_real_), list(coverage = Inf),
    list(seq_error = -0.01), list(seq_error = 0.5), list(seq_error = 1),
    list(seq_error = NA_real_), list(seq_error = Inf),
    list(alpha = 0), list(alpha = -0.01), list(alpha = 1),
    list(alpha = NA_real_), list(alpha = Inf),
    list(verbose = 1), list(verbose = NA), list(verbose = c(TRUE, FALSE))
  )
  for (override in invalid) {
    args <- utils::modifyList(tdt_ngs_mssn_args(), override)
    expect_error(do.call(tdt_ngs_mssn, args))
  }
})

test_that("MSSN source retains the intended analytic raw-read scope", {
  source_text <- paste(deparse(body(tdt_ngs_mssn)), collapse = "\n")
  expect_equal(
    sum(gregexpr(".tdt_ngs_information(", source_text, fixed = TRUE)[[1]] > 0),
    1
  )
  forbidden <- c(
    "ngs_genotype_error_matrix", ".ngs_call_genotype_ml",
    "tdt_power(", "tdt_mssn(", "lambda_from_gTgNT",
    ".tdt_ngs_ncp("
  )
  for (name in forbidden) {
    expect_false(grepl(name, source_text, fixed = TRUE))
  }
  expect_identical(
    names(formals(tdt_ngs_mssn)),
    c("power", "pd", "R1", "coverage", "seq_error", "alpha", "verbose")
  )
})
