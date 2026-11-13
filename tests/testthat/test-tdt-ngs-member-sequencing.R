tdt_ngs_member_args <- function(...) {
  utils::modifyList(list(
    N = 1200, pd = 0.325, R1 = 1.2,
    coverage = 10, seq_error = 0.01,
    alpha = 0.001, verbose = FALSE
  ), list(...))
}

test_that("explicit equal trio depth preserves legacy power and MSSN", {
  legacy_power <- do.call(tdt_ngs_power, tdt_ngs_member_args())
  explicit_power <- tdt_ngs_power(
    N = 1200, pd = 0.325, R1 = 1.2,
    father_coverage = 10, mother_coverage = 10, child_coverage = 10,
    epsilon0 = 0.01, epsilon1 = 0.01,
    alpha = 0.001, verbose = FALSE
  )
  legacy_mssn <- tdt_ngs_mssn(
    power = 0.8, pd = 0.325, R1 = 1.2,
    coverage = 10, seq_error = 0.01, alpha = 0.001, verbose = FALSE
  )
  explicit_mssn <- tdt_ngs_mssn(
    power = 0.8, pd = 0.325, R1 = 1.2,
    father_coverage = 10, mother_coverage = 10, child_coverage = 10,
    epsilon0 = 0.01, epsilon1 = 0.01,
    alpha = 0.001, verbose = FALSE
  )

  for (field in c("lambda", "power", "efficient_information",
                  "information_matrix")) {
    expect_identical(explicit_power[[field]], legacy_power[[field]])
  }
  for (field in c("MSSN_trios", "N_trios_continuous", "achieved_power",
                  "efficient_information", "information_matrix")) {
    expect_identical(explicit_mssn[[field]], legacy_mssn[[field]])
  }
  expect_null(explicit_power$coverage)
  expect_null(explicit_power$seq_error)
  expect_true(explicit_power$symmetric_error)
})

test_that("member overrides fall back independently to common coverage", {
  child <- do.call(tdt_ngs_power, tdt_ngs_member_args(child_coverage = 20))
  father <- do.call(tdt_ngs_power, tdt_ngs_member_args(father_coverage = 6))

  expect_identical(
    c(child$father_coverage, child$mother_coverage, child$child_coverage),
    c(10, 10, 20)
  )
  expect_identical(
    c(father$father_coverage, father$mother_coverage, father$child_coverage),
    c(6, 10, 10)
  )
  expect_identical(child$coverage, 10)
  expect_identical(father$coverage, 10)
  expect_false(child$sequencing$equal_coverage)
  expect_false(father$sequencing$equal_coverage)
})

test_that("unequal member coverage reaches the lower information kernel", {
  public <- tdt_ngs_power(
    N = 900, pd = 0.30, R1 = 1.25,
    father_coverage = 6, mother_coverage = 10, child_coverage = 20,
    seq_error = 0.01, alpha = 0.01, verbose = FALSE
  )
  kernel <- .tdt_ngs_ncp(
    N = 900, pd = 0.30, R1 = 1.25,
    coverage = c(father = 6, mother = 10, child = 20),
    seq_error = 0.01
  )

  expect_identical(public$lambda, kernel$lambda)
  expect_identical(public$efficient_information,
                   kernel$efficient_information)
  expect_identical(public$information_matrix, kernel$information_matrix)
  expect_identical(kernel$effective_coverage,
                   c(father = 6, mother = 10, child = 20))
})

test_that("symmetric and partial directional error fallbacks resolve exactly", {
  legacy <- do.call(tdt_ngs_power, tdt_ngs_member_args())
  explicit <- do.call(tdt_ngs_power, tdt_ngs_member_args(
    epsilon0 = 0.01, epsilon1 = 0.01
  ))
  left <- do.call(tdt_ngs_power, tdt_ngs_member_args(epsilon0 = 0.005))
  right <- do.call(tdt_ngs_power, tdt_ngs_member_args(epsilon1 = 0.015))

  expect_identical(explicit$lambda, legacy$lambda)
  expect_identical(explicit$information_matrix, legacy$information_matrix)
  expect_identical(c(left$epsilon0, left$epsilon1), c(0.005, 0.01))
  expect_identical(c(right$epsilon0, right$epsilon1), c(0.01, 0.015))
  expect_false(left$symmetric_error)
  expect_false(right$symmetric_error)
})

test_that("directional error follows the existing read-probability model", {
  epsilon0 <- 0.005
  epsilon1 <- 0.015
  expected <- epsilon0 + ((1 - epsilon0 - epsilon1) / 2) * 0:2
  observed <- vapply(
    0:2, .tdt_ngs_read_probability, numeric(1),
    epsilon0 = epsilon0, epsilon1 = epsilon1
  )
  expect_identical(observed, expected)

  public <- tdt_ngs_power(
    N = 1000, pd = 0.325, R1 = 1.2, coverage = 12,
    epsilon0 = epsilon0, epsilon1 = epsilon1,
    alpha = 0.05, verbose = FALSE
  )
  kernel <- .tdt_ngs_ncp(
    1000, 0.325, 1.2, 12, epsilon0 = epsilon0, epsilon1 = epsilon1
  )
  expect_identical(public$lambda, kernel$lambda)
  expect_identical(public$information_matrix, kernel$information_matrix)
  expect_identical(public$model_info$sequencing_error, "directional")
})

test_that("unequal coverage and directional error work simultaneously", {
  public <- tdt_ngs_power(
    N = 1000, pd = 0.325, R1 = 1.2,
    father_coverage = 6, mother_coverage = 12, child_coverage = 20,
    epsilon0 = 0.005, epsilon1 = 0.015,
    alpha = 5e-8, verbose = FALSE
  )
  kernel <- .tdt_ngs_ncp(
    N = 1000, pd = 0.325, R1 = 1.2,
    coverage = c(father = 6, mother = 12, child = 20),
    epsilon0 = 0.005, epsilon1 = 0.015
  )
  expected_power <- stats::pchisq(
    stats::qchisq(1 - 5e-8, df = 1),
    df = 1, ncp = kernel$lambda, lower.tail = FALSE
  )

  expect_identical(public$lambda, kernel$lambda)
  expect_identical(public$efficient_information,
                   kernel$efficient_information)
  expect_equal(public$power, expected_power, tolerance = 1e-15)
  expect_identical(public$model_info$coverage_model,
                   "member_specific_fixed")
  expect_identical(public$model_info$sequencing_error, "directional")
})

test_that("directional errors and member depths reject invalid inputs", {
  base <- list(N = 100, pd = 0.3, R1 = 1.2, coverage = 6,
               seq_error = 0.01, verbose = FALSE)
  errors <- list(
    list(epsilon0 = -0.01), list(epsilon1 = -0.01),
    list(epsilon0 = 0.6, epsilon1 = 0.4),
    list(epsilon0 = NA_real_), list(epsilon1 = Inf),
    list(epsilon0 = "0.01")
  )
  for (values in errors) {
    expect_error(do.call(tdt_ngs_power, utils::modifyList(base, values)),
                 "epsilon")
  }

  coverage_errors <- list(
    list(father_coverage = 1), list(mother_coverage = 1),
    list(child_coverage = 1), list(mother_coverage = 0),
    list(child_coverage = 1.5), list(father_coverage = NA_real_),
    list(mother_coverage = Inf), list(child_coverage = "6")
  )
  for (values in coverage_errors) {
    member <- names(values)
    expect_error(do.call(tdt_ngs_power, utils::modifyList(base, values)),
                 member, fixed = TRUE)
  }
  expect_error(do.call(tdt_ngs_power, utils::modifyList(
    base, list(seq_error = 0.5)
  )), "seq_error", fixed = TRUE)

  mssn <- c(list(power = 0.8), base[setdiff(names(base), "N")])
  expect_error(do.call(tdt_ngs_mssn, utils::modifyList(
    mssn, list(child_coverage = 1)
  )), "child_coverage", fixed = TRUE)
  expect_error(do.call(tdt_ngs_mssn, utils::modifyList(
    mssn, list(epsilon0 = 0.7, epsilon1 = 0.3)
  )), "epsilon0 + epsilon1", fixed = TRUE)
})

test_that("MSSN attains target under the three new sequencing designs", {
  designs <- list(
    list(
      father_coverage = 6, mother_coverage = 10, child_coverage = 20,
      seq_error = 0.01
    ),
    list(coverage = 10, epsilon0 = 0.005, epsilon1 = 0.015),
    list(
      father_coverage = 6, mother_coverage = 12, child_coverage = 20,
      epsilon0 = 0.005, epsilon1 = 0.015
    )
  )
  for (design in designs) {
    common <- c(list(
      pd = 0.325, R1 = 1.2, alpha = 0.05, verbose = FALSE
    ), design)
    result <- do.call(tdt_ngs_mssn, c(list(power = 0.8), common))
    achieved <- do.call(tdt_ngs_power, c(
      list(N = result$MSSN_trios), common
    ))
    expect_gte(achieved$power, 0.8)
    expect_equal(achieved$lambda, result$achieved_lambda, tolerance = 2e-15)
    if (result$MSSN_trios > 1) {
      previous <- do.call(tdt_ngs_power, c(
        list(N = result$MSSN_trios - 1), common
      ))
      expect_lt(previous$power, 0.8)
    }
  }
})

test_that("directional member read distributions normalize", {
  coverage <- c(father = 6, mother = 12, child = 20)
  epsilon0 <- 0.005
  epsilon1 <- 0.015
  for (depth in coverage) {
    for (genotype in 0:2) {
      q <- .tdt_ngs_read_probability(genotype, epsilon0, epsilon1)
      probabilities <- stats::dbinom(0:depth, depth, q)
      expect_equal(sum(probabilities), 1, tolerance = 1e-15)
      expect_true(all(probabilities >= 0 & probabilities <= 1))
    }
  }
  genotypes <- c(father = 0, mother = 1, child = 2)
  grid <- expand.grid(
    father = 0:coverage[["father"]],
    mother = 0:coverage[["mother"]],
    child = 0:coverage[["child"]]
  )
  joint <- apply(grid, 1, function(x) {
    .tdt_ngs_read_likelihood(
      x, coverage, genotypes, epsilon0, epsilon1
    )
  })
  expect_equal(sum(joint), 1, tolerance = 1e-14)
})

test_that("high-quality sequencing approaches ordinary no-error TDT", {
  ngs <- tdt_ngs_power(
    N = 1000, pd = 0.325, R1 = 1.01,
    father_coverage = 100, mother_coverage = 100, child_coverage = 100,
    epsilon0 = 0, epsilon1 = 0,
    alpha = 0.05, verbose = FALSE
  )
  ordinary <- tdt_power(
    N = 1000, input_mode = "model_based",
    pd = 0.325, prev = 0.05, R1 = 1.01, R2 = 1.01^2,
    alpha = 0.05, verbose = FALSE
  )

  expect_equal(ngs$power, ordinary$power$no_error, tolerance = 1e-4)
  expect_equal(ngs$lambda, ordinary$lambda$no_error, tolerance = 0.002)
})

test_that("console output reports member-specific directional settings", {
  out <- tdt_ngs_power(
    N = 100, pd = 0.325, R1 = 1.2,
    father_coverage = 6, mother_coverage = 12, child_coverage = 20,
    epsilon0 = 0.005, epsilon1 = 0.015, verbose = FALSE
  )
  printed <- capture.output(.paweh_print_tdt_ngs_power(out), type = "message")
  expect_true(any(grepl("Father coverage", printed, fixed = TRUE)))
  expect_true(any(grepl("epsilon0", printed, fixed = TRUE)))
  expect_true(any(grepl("epsilon1", printed, fixed = TRUE)))
})
