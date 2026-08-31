tdt_ngs_power_args <- function(...) {
  utils::modifyList(
    list(
      N = 5000, pd = 0.325, R1 = 1.2,
      coverage = 12, seq_error = 0.005,
      alpha = 0.05, verbose = FALSE
    ),
    list(...)
  )
}

tdt_ngs_power_quiet <- function(...) {
  do.call(tdt_ngs_power, tdt_ngs_power_args(...))
}

test_that("public power is a thin wrapper around the frozen S8 NCP", {
  designs <- list(
    list(N = 800, pd = 0.15, R1 = 1.1, coverage = 2, seq_error = 0),
    list(N = 1200, pd = 0.325, R1 = 1.2, coverage = 4,
         seq_error = 0.005),
    list(N = 500, pd = 0.5, R1 = 1.3, coverage = 12, seq_error = 0.01),
    list(N = 250, pd = 0.675, R1 = 0.9, coverage = 20,
         seq_error = 0.005)
  )

  for (design in designs) {
    public <- do.call(
      tdt_ngs_power,
      c(design, list(alpha = 0.001, verbose = FALSE))
    )
    kernel <- do.call(.tdt_ngs_ncp, design)
    expect_equal(public$lambda, kernel$lambda, tolerance = 1e-14)
    expect_equal(
      public$efficient_information,
      kernel$efficient_information,
      tolerance = 1e-14
    )
    expect_equal(public$information_matrix, kernel$information_matrix,
                 tolerance = 1e-14)
    expect_equal(public$nuisance_rcond, kernel$nuisance_rcond,
                 tolerance = 1e-14)
  }

  source_text <- paste(deparse(body(tdt_ngs_power)), collapse = "\n")
  calls <- gregexpr(".tdt_ngs_ncp(", source_text, fixed = TRUE)[[1]]
  expect_equal(sum(calls > 0), 1)
})

test_that("power is the independent one-df noncentral chi-square tail", {
  for (alpha in c(0.05, 0.001, 5e-8)) {
    result <- tdt_ngs_power_quiet(alpha = alpha)
    critical <- stats::qchisq(1 - alpha, df = 1)
    expected <- stats::pchisq(
      critical, df = 1, ncp = result$lambda, lower.tail = FALSE
    )
    expect_equal(result$power, expected, tolerance = 1e-15)
  }
})

test_that("the multiplicative null gives lambda zero and power alpha", {
  designs <- list(
    list(pd = 0.15, coverage = 2, seq_error = 0, alpha = 0.05),
    list(pd = 0.325, coverage = 4, seq_error = 0.005, alpha = 0.001),
    list(pd = 0.5, coverage = 12, seq_error = 0.01, alpha = 5e-8)
  )
  for (design in designs) {
    result <- do.call(
      tdt_ngs_power,
      c(list(N = 1000, R1 = 1, verbose = FALSE), design)
    )
    expect_identical(result$lambda, 0)
    expect_lt(abs(result$power - design$alpha), 5e-16)
  }
})

test_that("Chapter 5 coverage fixtures remain fixed through the public API", {
  fixtures <- data.frame(
    coverage = c(4, 12, 20, 44),
    information = c(
      0.166296146079271,
      0.216362516007223,
      0.219166852960543,
      0.219374904964103
    ),
    lambda = c(
      27.6393757408911,
      35.9606943225104,
      36.426791250096,
      36.4613706894614
    )
  )

  for (i in seq_len(nrow(fixtures))) {
    result <- tdt_ngs_power_quiet(
      coverage = fixtures$coverage[[i]], alpha = 5e-8
    )
    expect_equal(result$efficient_information, fixtures$information[[i]],
                 tolerance = 2e-13)
    expect_equal(result$lambda, fixtures$lambda[[i]], tolerance = 2e-12)
    expected_power <- stats::pchisq(
      stats::qchisq(1 - 5e-8, df = 1),
      df = 1, ncp = fixtures$lambda[[i]], lower.tail = FALSE
    )
    expect_equal(result$power, expected_power, tolerance = 1e-14)
  }
})

test_that("Chapter 5 sequencing-error fixtures remain fixed", {
  fixtures <- data.frame(
    seq_error = c(0, 0.005, 0.01),
    information = c(
      0.219190883720878,
      0.216362516007223,
      0.214412722729257
    ),
    lambda = c(
      36.4307853006492,
      35.9606943225104,
      35.6366274677014
    )
  )

  for (i in seq_len(nrow(fixtures))) {
    result <- tdt_ngs_power_quiet(seq_error = fixtures$seq_error[[i]])
    expect_equal(result$efficient_information, fixtures$information[[i]],
                 tolerance = 2e-13)
    expect_equal(result$lambda, fixtures$lambda[[i]], tolerance = 2e-12)
  }
})

test_that("lambda scales linearly with trio count through the public API", {
  base <- tdt_ngs_power_quiet(N = 700)
  twice <- tdt_ngs_power_quiet(N = 1400)
  fivefold <- tdt_ngs_power_quiet(N = 3500)

  expect_equal(twice$lambda, 2 * base$lambda, tolerance = 1e-13)
  expect_equal(fivefold$lambda, 5 * base$lambda, tolerance = 1e-13)
  expect_gte(twice$power, base$power)
  expect_gte(fivefold$power, twice$power)
})

test_that("R1 metadata follows the multiplicative transmission model", {
  for (R1 in c(1, 1.1, 1.2, 1.3)) {
    result <- tdt_ngs_power_quiet(R1 = R1, coverage = 4)
    expect_equal(result$t, R1 / (1 + R1), tolerance = 1e-15)
    expect_equal(result$delta, log(R1), tolerance = 2e-15)
    expect_equal(result$R2, R1^2, tolerance = 1e-15)
  }
})

test_that("allele-label reflection is preserved without recoding pd", {
  for (pd in c(0.15, 0.325)) {
    left <- tdt_ngs_power_quiet(pd = pd, coverage = 4)
    right <- tdt_ngs_power_quiet(pd = 1 - pd, coverage = 4)
    expect_identical(left$pd, pd)
    expect_identical(right$pd, 1 - pd)
    expect_equal(left$lambda, right$lambda, tolerance = 2e-12)
    expect_equal(left$power, right$power, tolerance = 2e-13)
  }
})

test_that("coverage one has a public identifiability error", {
  expect_error(
    tdt_ngs_power_quiet(coverage = 1),
    "coverage = 1 is unsupported.*not identifiable"
  )
})

test_that("coverage two is identifiable over the representative public grid", {
  grid <- expand.grid(
    pd = c(0.15, 0.325, 0.5),
    R1 = c(1.1, 1.2, 1.3),
    seq_error = c(0, 0.005, 0.01),
    KEEP.OUT.ATTRS = FALSE
  )
  for (i in seq_len(nrow(grid))) {
    result <- tdt_ngs_power_quiet(
      N = 200,
      pd = grid$pd[[i]],
      R1 = grid$R1[[i]],
      coverage = 2,
      seq_error = grid$seq_error[[i]]
    )
    expect_true(is.finite(result$lambda))
    expect_gte(result$lambda, 0)
    expect_true(is.finite(result$efficient_information))
    expect_gt(result$efficient_information, 0)
    expect_true(is.finite(result$nuisance_rcond))
    expect_gt(result$nuisance_rcond, 1e-12)
  }
})

test_that("result structure is transparent, finite, and returned invisibly", {
  visible <- withVisible(tdt_ngs_power_quiet())
  result <- visible$value

  expect_false(visible$visible)
  expect_s3_class(result, "tdt_ngs_power", exact = TRUE)
  expect_named(
    result,
    c(
      "N", "alpha", "power", "lambda", "pd", "R1", "R2", "t",
      "delta", "coverage", "seq_error", "efficient_information",
      "information_matrix", "nuisance_rcond", "score_mean", "model_info"
    )
  )
  expect_true(is.finite(result$power) && result$power >= 0 && result$power <= 1)
  expect_true(is.finite(result$lambda) && result$lambda >= 0)
  expect_true(is.finite(result$efficient_information))
  expect_gte(result$efficient_information, 0)
  expect_equal(dim(result$information_matrix), c(11, 11))
  expect_true(all(is.finite(result$information_matrix)))
  expect_equal(result$information_matrix, t(result$information_matrix),
               tolerance = 1e-14)
  expect_identical(result$N, 5000)
  expect_identical(result$model_info$test, "TDT1-NGS")
  expect_identical(result$model_info$inheritance, "multiplicative")
  expect_identical(result$model_info$coverage_model, "equal_fixed")
  expect_identical(result$model_info$sequencing_error, "symmetric")
  expect_identical(
    result$model_info$trio_type,
    "father-mother-affected-child"
  )
  expect_identical(
    result$model_info$likelihood,
    "raw_read_counts_with_latent_trio_states"
  )
  expect_identical(result$model_info$information_evaluation, "null")
})

test_that("verbose and print methods are concise", {
  expect_output(
    tdt_ngs_power_quiet(verbose = FALSE),
    NA
  )
  expect_output(
    tdt_ngs_power_quiet(verbose = TRUE),
    "TDT1-NGS analytic power"
  )
  result <- tdt_ngs_power_quiet()
  printed <- capture.output(print(result))
  expect_lte(length(printed), 6)
  expect_true(any(grepl("Affected-child trios", printed, fixed = TRUE)))
  expect_false(any(grepl("information_matrix", printed, fixed = TRUE)))
})

test_that("less stringent alpha gives no less power", {
  powers <- vapply(c(0.05, 0.001, 5e-8), function(alpha) {
    tdt_ngs_power_quiet(alpha = alpha)$power
  }, numeric(1))
  expect_true(powers[[1]] >= powers[[2]])
  expect_true(powers[[2]] >= powers[[3]])
})

test_that("public input validation rejects invalid designs", {
  invalid <- list(
    list(N = 0), list(N = -1), list(N = 1.5), list(N = NA_real_),
    list(N = Inf),
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
    args <- utils::modifyList(tdt_ngs_power_args(), override)
    expect_error(do.call(tdt_ngs_power, args))
  }
})

test_that("public scope contains no hard-call or ordinary-TDT substitution", {
  source_text <- paste(deparse(body(tdt_ngs_power)), collapse = "\n")
  forbidden <- c(
    "ngs_genotype_error_matrix", ".ngs_call_genotype_ml",
    "tdt_power(", "tdt_mssn(", "lambda_from_gTgNT"
  )
  for (name in forbidden) {
    expect_false(grepl(name, source_text, fixed = TRUE))
  }
  expect_identical(
    names(formals(tdt_ngs_power)),
    c("N", "pd", "R1", "coverage", "seq_error", "alpha", "verbose")
  )
})
