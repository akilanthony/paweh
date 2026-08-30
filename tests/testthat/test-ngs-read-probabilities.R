test_that("per-read alternate probabilities follow the symmetric error model", {
  for (seq_error in c(0, 0.001, 0.01, 0.04, 0.49)) {
    expect_equal(.ngs_alt_read_probability(0, seq_error), seq_error)
    expect_equal(.ngs_alt_read_probability(1, seq_error), 0.5)
    expect_equal(.ngs_alt_read_probability(2, seq_error), 1 - seq_error)
  }
})

test_that("heterozygote read probabilities and distributions are error-invariant", {
  errors <- c(0, 0.001, 0.01, 0.04, 0.25, 0.49)

  for (seq_error in errors) {
    expect_equal(.ngs_alt_read_probability(1, seq_error), 0.5)
  }

  reference <- .ngs_read_count_distribution(10, 1, errors[[1]])$probability
  for (seq_error in errors[-1]) {
    observed <- .ngs_read_count_distribution(10, 1, seq_error)$probability
    expect_equal(observed, reference, tolerance = 1e-15)
  }
})

test_that("complete read-count distributions are valid binomial PMFs", {
  combinations <- expand.grid(
    coverage = c(1, 4, 10, 50),
    genotype = 0:2,
    seq_error = c(0, 0.001, 0.01, 0.04),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(combinations))) {
    inputs <- combinations[i, ]
    distribution <- .ngs_read_count_distribution(
      coverage = inputs$coverage,
      genotype = inputs$genotype,
      seq_error = inputs$seq_error
    )

    expect_equal(distribution$alt_reads, 0:inputs$coverage)
    expect_true(all(is.finite(distribution$probability)))
    expect_true(all(distribution$probability >= 0))
    expect_true(all(distribution$probability <= 1))
    expect_equal(sum(distribution$probability), 1, tolerance = 1e-14)
  }
})

test_that("homozygote distributions have the required symmetry", {
  for (coverage in c(1, 4, 10, 50)) {
    for (seq_error in c(0, 0.001, 0.01, 0.04, 0.49)) {
      g0 <- .ngs_read_count_distribution(coverage, 0, seq_error)$probability
      g2 <- .ngs_read_count_distribution(coverage, 2, seq_error)$probability
      expect_equal(g0, rev(g2), tolerance = 1e-14)
    }
  }
})

test_that("zero sequencing error gives the required limiting distributions", {
  coverage <- 10
  g0 <- .ngs_read_count_distribution(coverage, 0, 0)$probability
  g1 <- .ngs_read_count_distribution(coverage, 1, 0)$probability
  g2 <- .ngs_read_count_distribution(coverage, 2, 0)$probability

  expect_equal(g0, c(1, rep(0, coverage)))
  expect_equal(g2, c(rep(0, coverage), 1))
  expect_equal(g1, stats::dbinom(0:coverage, coverage, 0.5))
})

test_that("coverage-four heterozygote probabilities match the known expansion", {
  expected <- c(1, 4, 6, 4, 1) / 16

  for (seq_error in c(0, 0.01, 0.20, 0.49)) {
    observed <- .ngs_read_count_distribution(4, 1, seq_error)$probability
    expect_equal(observed, expected, tolerance = 1e-15)
  }
})

test_that("single-read probabilities agree with the per-read model", {
  seq_error <- 0.04

  expect_equal(.ngs_read_count_prob(1, 1, 0, seq_error), seq_error)
  expect_equal(.ngs_read_count_prob(1, 1, 1, seq_error), 0.5)
  expect_equal(.ngs_read_count_prob(1, 1, 2, seq_error), 1 - seq_error)
})

test_that("single-count helper agrees with the complete distribution", {
  coverage <- 12

  for (genotype in 0:2) {
    for (seq_error in c(0, 0.01, 0.20)) {
      distribution <- .ngs_read_count_distribution(
        coverage, genotype, seq_error
      )$probability

      for (alt_reads in 0:coverage) {
        expect_equal(
          .ngs_read_count_prob(
            alt_reads, coverage, genotype, seq_error
          ),
          distribution[[alt_reads + 1L]],
          tolerance = 1e-15
        )
      }
    }
  }
})

test_that("invalid genotypes and sequencing-error rates are rejected", {
  invalid_genotypes <- list(-1, 3, 1.5, c(0, 1), NA_real_, NaN, Inf)
  for (genotype in invalid_genotypes) {
    expect_error(
      .ngs_alt_read_probability(genotype, 0.01),
      "genotype"
    )
  }

  invalid_errors <- list(-0.01, 0.5, 0.51, NA_real_, NaN, Inf, c(0.01, 0.02))
  for (seq_error in invalid_errors) {
    expect_error(
      .ngs_alt_read_probability(1, seq_error),
      "seq_error"
    )
  }
})

test_that("invalid coverage and alternate-read counts are rejected", {
  invalid_coverage <- list(0, -1, 1.5, NA_real_, Inf)
  for (coverage in invalid_coverage) {
    expect_error(
      .ngs_read_count_distribution(coverage, 1, 0.01),
      "coverage"
    )
    expect_error(
      .ngs_read_count_prob(0, coverage, 1, 0.01),
      "coverage"
    )
  }

  invalid_alt_reads <- list(-1, 11, 1.5, NA_real_, Inf)
  for (alt_reads in invalid_alt_reads) {
    expect_error(
      .ngs_read_count_prob(alt_reads, 10, 1, 0.01),
      "alt_reads"
    )
  }
})
