test_that("ML caller returns a genotype attaining the maximum likelihood", {
  combinations <- expand.grid(
    coverage = c(1, 2, 3, 4, 10, 20),
    seq_error = c(0, 0.001, 0.01, 0.04, 0.20, 0.49),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(combinations))) {
    coverage <- combinations$coverage[[i]]
    seq_error <- combinations$seq_error[[i]]
    q <- c(seq_error, 0.5, 1 - seq_error)

    for (alt_reads in 0:coverage) {
      called <- .ngs_call_genotype_ml(alt_reads, coverage, seq_error)
      log_likelihood <- stats::dbinom(
        alt_reads, coverage, q, log = TRUE
      )

      expect_true(called %in% 0:2)
      expect_gte(
        log_likelihood[[called + 1L]],
        max(log_likelihood) - 1e-12
      )
    }
  }
})

test_that("documented numerical tie rule is deterministic and symmetric", {
  coverage <- 4
  lower_tie <- stats::uniroot(
    function(seq_error) {
      stats::dbinom(1, coverage, seq_error, log = TRUE) -
        stats::dbinom(1, coverage, 0.5, log = TRUE)
    },
    interval = c(1e-12, 0.25),
    tol = 1e-14
  )$root

  lower_likelihood <- stats::dbinom(
    1, coverage, c(lower_tie, 0.5, 1 - lower_tie), log = TRUE
  )
  upper_likelihood <- stats::dbinom(
    3, coverage, c(lower_tie, 0.5, 1 - lower_tie), log = TRUE
  )

  expect_equal(lower_likelihood[[1]], lower_likelihood[[2]], tolerance = 1e-12)
  expect_equal(upper_likelihood[[2]], upper_likelihood[[3]], tolerance = 1e-12)
  expect_equal(.ngs_call_genotype_ml(1, coverage, lower_tie), 0L)
  expect_equal(.ngs_call_genotype_ml(3, coverage, lower_tie), 2L)
})

test_that("transition matrices have the documented structure and orientation", {
  E <- ngs_genotype_error_matrix(coverage = 10, seq_error = 0.01)

  expect_type(E, "double")
  expect_equal(dim(E), c(3L, 3L))
  expect_equal(rownames(E), c("true_0", "true_1", "true_2"))
  expect_equal(colnames(E), c("called_0", "called_1", "called_2"))
  expect_true(all(is.finite(E)))
  expect_true(all(E >= 0 & E <= 1))

  g_true <- c(0.70, 0.25, 0.05)
  g_called <- as.numeric(t(E) %*% g_true)
  independently_assembled <- colSums(E * g_true)

  expect_true(all(g_called >= 0 & g_called <= 1))
  expect_equal(sum(g_called), 1, tolerance = 1e-14)
  expect_equal(g_called, unname(independently_assembled), tolerance = 1e-15)
})

test_that("every tested transition matrix is row-stochastic", {
  for (coverage in c(1, 2, 3, 4, 10, 20, 50)) {
    for (seq_error in c(0, 0.001, 0.01, 0.04, 0.20, 0.499)) {
      E <- ngs_genotype_error_matrix(coverage, seq_error)

      expect_equal(unname(rowSums(E)), rep(1, 3), tolerance = 1e-14)
      expect_true(all(is.finite(E)))
      expect_true(all(E >= 0 & E <= 1))
    }
  }
})

test_that("zero-error matrices retain finite-depth heterozygote ambiguity", {
  for (coverage in c(1, 2, 3, 4, 10, 20)) {
    E <- ngs_genotype_error_matrix(coverage, 0)
    heterozygote_tail <- 2^(-coverage)
    expected <- rbind(
      c(1, 0, 0),
      c(heterozygote_tail, 1 - 2 * heterozygote_tail, heterozygote_tail),
      c(0, 0, 1)
    )
    dimnames(expected) <- dimnames(E)

    expect_equal(E, expected, tolerance = 1e-15)
  }

  coverage_one <- ngs_genotype_error_matrix(1, 0)
  expect_equal(unname(coverage_one["true_1", ]), c(0.5, 0, 0.5))
  expect_false(isTRUE(all.equal(coverage_one, diag(3))))
})

test_that("homozygote reflection symmetry holds", {
  for (coverage in c(1, 2, 3, 4, 10, 20, 50)) {
    for (seq_error in c(0, 0.001, 0.01, 0.04, 0.20, 0.499)) {
      E <- ngs_genotype_error_matrix(coverage, seq_error)

      expect_equal(E["true_0", "called_0"], E["true_2", "called_2"], tolerance = 1e-14)
      expect_equal(E["true_0", "called_1"], E["true_2", "called_1"], tolerance = 1e-14)
      expect_equal(E["true_0", "called_2"], E["true_2", "called_0"], tolerance = 1e-14)
    }
  }
})

test_that("heterozygote call distribution is symmetric", {
  for (coverage in c(1, 2, 3, 4, 10, 20, 50)) {
    for (seq_error in c(0, 0.001, 0.01, 0.04, 0.20, 0.499)) {
      E <- ngs_genotype_error_matrix(coverage, seq_error)
      expect_equal(
        E["true_1", "called_0"],
        E["true_1", "called_2"],
        tolerance = 1e-14
      )
    }
  }
})

test_that("high fixed coverage approaches perfect genotype classification", {
  for (seq_error in c(0.001, 0.01, 0.04)) {
    low <- ngs_genotype_error_matrix(4, seq_error)
    medium <- ngs_genotype_error_matrix(20, seq_error)
    high <- ngs_genotype_error_matrix(100, seq_error)

    expect_true(all(diag(high) > diag(low)))
    expect_true(all(diag(high) >= diag(medium)))
    expect_true(all(diag(high) > 0.999))
    expect_lt(max(abs(high - diag(3))), 0.001)
  }
})

test_that("matrix and caller inputs are validated", {
  invalid_coverage <- list(0, -1, 1.5, NA_real_, NaN, Inf, c(2, 3))
  for (coverage in invalid_coverage) {
    expect_error(
      ngs_genotype_error_matrix(coverage, 0.01),
      "coverage"
    )
  }

  invalid_errors <- list(-0.01, 0.5, 0.51, NA_real_, NaN, Inf, c(0.01, 0.02))
  for (seq_error in invalid_errors) {
    expect_error(
      ngs_genotype_error_matrix(10, seq_error),
      "seq_error"
    )
  }

  invalid_alt_reads <- list(-1, 11, 1.5, NA_real_, NaN, Inf, c(1, 2))
  for (alt_reads in invalid_alt_reads) {
    expect_error(
      .ngs_call_genotype_ml(alt_reads, 10, 0.01),
      "alt_reads"
    )
  }
})
