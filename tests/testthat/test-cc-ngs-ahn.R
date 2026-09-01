ahn_reference_ncp <- function(g_case, g_control, N_case, N_control, scores) {
  numerator <- N_case * N_control *
    sum(scores * (g_case - g_control))^2
  weighted_counts <- N_case * g_case + N_control * g_control
  denominator <- sum(scores^2 * weighted_counts) -
    sum(scores * weighted_counts)^2 / (N_case + N_control)
  numerator / denominator
}

test_that("Ahn Equation (1) matches a direct hand calculation", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  N_case <- 1000
  N_control <- 1200
  scores <- c(0, 1, 2)

  expected <- ahn_reference_ncp(
    g_case, g_control, N_case, N_control, scores
  )
  observed <- .cc_ahn_trend_ncp(
    g_case, g_control, N_case, N_control, scores
  )

  expect_equal(observed, expected, tolerance = 1e-14)
})

test_that("Ahn NCP equals the pooled theta parameterization", {
  settings <- list(
    list(c(0.45, 0.40, 0.15), c(0.55, 0.35, 0.10), 1000, 1200),
    list(c(0.25, 0.50, 0.25), c(0.36, 0.48, 0.16), 500, 500),
    list(c(0.70, 0.25, 0.05), c(0.60, 0.30, 0.10), 750, 1500)
  )
  scores <- c(0, 1, 2)

  for (setting in settings) {
    g_case <- setting[[1]]
    g_control <- setting[[2]]
    N_case <- setting[[3]]
    N_control <- setting[[4]]
    N <- N_case + N_control
    theta <- N_case / N
    pbar <- theta * g_case + (1 - theta) * g_control
    pooled_variance <- sum(scores^2 * pbar) - sum(scores * pbar)^2
    expected <- N * theta * (1 - theta) *
      sum(scores * (g_case - g_control))^2 / pooled_variance

    observed <- .cc_ahn_trend_ncp(
      g_case, g_control, N_case, N_control, scores
    )
    expect_equal(observed, expected, tolerance = 1e-12)
  }
})

test_that("Ahn NCP agrees with canonical PAWEH trend NCP", {
  settings <- list(
    list(c(0.45, 0.40, 0.15), c(0.55, 0.35, 0.10), 1000, 1200),
    list(c(0.25, 0.50, 0.25), c(0.36, 0.48, 0.16), 500, 750),
    list(c(0.70, 0.25, 0.05), c(0.60, 0.30, 0.10), 900, 450)
  )

  for (setting in settings) {
    g_case <- setting[[1]]
    g_control <- setting[[2]]
    N_case <- setting[[3]]
    N_control <- setting[[4]]
    scores <- c(0, 1, 2)
    existing <- .cc_power_test_results(
      g_case = g_case,
      g_ctrl = g_control,
      k = N_control / N_case,
      w = scores,
      N_case = N_case,
      alpha = 0.05
    )$trend$lambda

    observed <- .cc_ahn_trend_ncp(
      g_case, g_control, N_case, N_control, scores
    )
    expect_equal(observed, existing, tolerance = 1e-12)
  }
})

test_that("null genotype distributions have zero NCP", {
  g <- c(0.45, 0.40, 0.15)
  for (scores in list(c(0, 1, 2), c(0, 0, 1), c(-1, 0, 1))) {
    expect_equal(.cc_ahn_trend_ncp(g, g, 1000, 1200, scores), 0)
  }
})

test_that("trend NCP is invariant to affine score transformations", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  score_sets <- list(
    c(0, 1, 2),
    c(1, 2, 3),
    c(0, 2, 4),
    c(-1, 0, 1)
  )
  lambda <- vapply(
    score_sets,
    function(scores) {
      .cc_ahn_trend_ncp(g_case, g_control, 1000, 1200, scores)
    },
    numeric(1)
  )

  expect_equal(lambda, rep(lambda[[1]], length(lambda)), tolerance = 1e-12)
})

test_that("standard Ahn trend score vectors are accepted", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  score_sets <- list(
    recessive = c(0, 0, 1),
    dominant = c(0, 1, 1),
    additive = c(0, 1, 2)
  )

  for (scores in score_sets) {
    lambda <- .cc_ahn_trend_ncp(
      g_case, g_control, 1000, 1200, scores
    )
    expect_true(is.finite(lambda))
    expect_gte(lambda, 0)
  }
})

test_that("sequencing called frequencies match direct matrix multiplication", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  E <- ngs_genotype_error_matrix(coverage = 10, seq_error = 0.01)
  observed <- .cc_ngs_called_frequencies(
    g_case, g_control, coverage = 10, seq_error = 0.01
  )

  expect_equal(observed$E, E, tolerance = 1e-15)
  expect_equal(observed$case_true, g_case)
  expect_equal(observed$control_true, g_control)
  expect_equal(observed$case_called, as.numeric(t(E) %*% g_case), tolerance = 1e-15)
  expect_equal(observed$control_called, as.numeric(t(E) %*% g_control), tolerance = 1e-15)
})

test_that("sequencing called frequencies remain valid probability vectors", {
  genotype_vectors <- list(
    c(1, 0, 0),
    c(0, 1, 0),
    c(0, 0, 1),
    c(0.45, 0.40, 0.15),
    c(0.70, 0.25, 0.05)
  )

  for (g_case in genotype_vectors) {
    for (g_control in rev(genotype_vectors)) {
      for (coverage in c(1, 4, 20, 50)) {
        for (seq_error in c(0, 0.01, 0.04, 0.20, 0.499)) {
          observed <- .cc_ngs_called_frequencies(
            g_case, g_control, coverage, seq_error
          )
          for (g_called in list(observed$case_called, observed$control_called)) {
            expect_true(all(is.finite(g_called)))
            expect_true(all(g_called >= 0 & g_called <= 1))
            expect_equal(sum(g_called), 1, tolerance = 1e-14)
          }
        }
      }
    }
  }
})

test_that("sequencing NCP composition matches an independent calculation", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  N_case <- 1000
  N_control <- 1200
  coverage <- 12
  seq_error <- 0.04
  scores <- c(0, 1, 2)

  E <- ngs_genotype_error_matrix(coverage, seq_error)
  case_called <- as.numeric(t(E) %*% g_case)
  control_called <- as.numeric(t(E) %*% g_control)
  expected <- ahn_reference_ncp(
    case_called, control_called, N_case, N_control, scores
  )
  observed <- .cc_ngs_ahn_ncp(
    g_case, g_control, N_case, N_control,
    coverage, seq_error, scores
  )

  expect_equal(observed$lambda, expected, tolerance = 1e-12)
  expect_equal(observed$E, E, tolerance = 1e-15)
  expect_equal(observed$case_called, case_called, tolerance = 1e-15)
  expect_equal(observed$control_called, control_called, tolerance = 1e-15)
  expect_equal(observed$scores, scores)
  expect_equal(observed$N_case, N_case)
  expect_equal(observed$N_control, N_control)
})

test_that("high sequencing depth approaches the true-genotype NCP", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  direct <- .cc_ahn_trend_ncp(g_case, g_control, 1000, 1200)
  sequencing <- .cc_ngs_ahn_ncp(
    g_case, g_control, 1000, 1200,
    coverage = 500, seq_error = 0.01
  )$lambda

  expect_equal(sequencing, direct, tolerance = 1e-10)
})

test_that("zero-error NCP retains finite-depth heterozygote sampling", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  coverage <- 4
  tail <- 2^(-coverage)
  expected_E <- rbind(
    c(1, 0, 0),
    c(tail, 1 - 2 * tail, tail),
    c(0, 0, 1)
  )
  dimnames(expected_E) <- list(
    c("true_0", "true_1", "true_2"),
    c("called_0", "called_1", "called_2")
  )
  case_called <- as.numeric(t(expected_E) %*% g_case)
  control_called <- as.numeric(t(expected_E) %*% g_control)
  expected <- ahn_reference_ncp(
    case_called, control_called, 1000, 1200, c(0, 1, 2)
  )
  observed <- .cc_ngs_ahn_ncp(
    g_case, g_control, 1000, 1200,
    coverage = coverage, seq_error = 0
  )
  ordinary <- .cc_ahn_trend_ncp(g_case, g_control, 1000, 1200)

  expect_equal(observed$E, expected_E, tolerance = 1e-15)
  expect_equal(observed$lambda, expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(observed$lambda, ordinary, tolerance = 1e-12)))
})

test_that("Ahn NCP scales linearly with common sample-size scaling", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  baseline <- .cc_ahn_trend_ncp(g_case, g_control, 1000, 1200)

  for (scale in c(0.5, 2, 3.5)) {
    scaled <- .cc_ahn_trend_ncp(
      g_case, g_control, scale * 1000, scale * 1200
    )
    expect_equal(scaled, scale * baseline, tolerance = 1e-12)
  }
})

test_that("swapping case and control labels leaves the NCP unchanged", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  forward <- .cc_ahn_trend_ncp(g_case, g_control, 1000, 1200)
  reverse <- .cc_ahn_trend_ncp(g_control, g_case, 1200, 1000)

  expect_equal(forward, reverse, tolerance = 1e-14)
})

test_that("Ahn and sequencing helpers reject malformed inputs", {
  valid <- c(0.45, 0.40, 0.15)
  invalid_frequencies <- list(
    c(0.5, 0.5),
    c(-0.1, 0.6, 0.5),
    c(1.1, 0, -0.1),
    c(0.4, 0.4, 0.1),
    c(NA_real_, 0.5, 0.5),
    c(NaN, 0.5, 0.5),
    c(Inf, 0, 0)
  )
  for (g in invalid_frequencies) {
    expect_error(.cc_ahn_trend_ncp(g, valid, 1000, 1200), "g_case")
    expect_error(.cc_ngs_called_frequencies(valid, g, 10, 0.01), "g_control")
  }

  for (N_case in list(0, -1, NA_real_, NaN, Inf, c(100, 200))) {
    expect_error(
      .cc_ahn_trend_ncp(valid, valid, N_case, 1200),
      "N_case"
    )
  }
  for (N_control in list(0, -1, NA_real_, NaN, Inf, c(100, 200))) {
    expect_error(
      .cc_ahn_trend_ncp(valid, valid, 1000, N_control),
      "N_control"
    )
  }

  invalid_scores <- list(c(0, 1), c(0, NA, 2), c(0, Inf, 2), c(1, 1, 1))
  for (scores in invalid_scores) {
    expect_error(
      .cc_ahn_trend_ncp(valid, valid, 1000, 1200, scores),
      "scores"
    )
  }
})
