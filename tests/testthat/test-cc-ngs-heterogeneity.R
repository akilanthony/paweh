cc_ngs_het_power_args <- function(...) {
  utils::modifyList(
    list(
      N_case = 1000, alpha = 0.05,
      prev = 0.05, pd = 0.30, R2 = 1.8,
      coverage = 20, seq_error = 0.01,
      MOI = "M", k = 1, verbose = FALSE
    ),
    list(...)
  )
}

cc_ngs_het_mssn_args <- function(...) {
  args <- cc_ngs_het_power_args(...)
  args$N_case <- NULL
  args$power <- 0.80
  args
}

cc_ngs_het_scores <- function(MOI) {
  switch(MOI, M = c(0, 1, 2), D = c(0, 1, 1), Rec = c(0, 0, 1))
}

cc_ngs_het_achieved_power <- function(result, N_case) {
  N_control <- ceiling(result$k * N_case)
  lambda <- .cc_ahn_trend_ncp(
    g_case = result$freqs$case_called,
    g_control = result$freqs$control_called,
    N_case = N_case,
    N_control = N_control,
    scores = result$scores
  )
  .cc_ngs_chisq_power(lambda, result$alpha)
}

test_that("CC-NGS uses the canonical ordinary-CC locus mixture", {
  for (MOI in c("M", "D", "Rec")) {
    model <- .cc_model_genotype_frequencies(
      pd = 0.30, R2 = 1.8, MOI = MOI, prev = 0.05
    )
    for (pi in c(0, 0.25, 0.5, 0.75, 1)) {
      observed <- .cc_ngs_apply_locus_heterogeneity(
        model$case, model$control, locus_het = TRUE, pi = pi
      )
      canonical_helper <- cc_apply_locus_het(
        g_case_assoc = model$case,
        g_ctrl = model$control,
        pi = pi
      )
      expect_equal(
        observed$g_case_after_locus_het,
        canonical_helper$g_case_het,
        tolerance = 1e-15
      )
      if (pi > 0) {
        canonical_public <- cc_power(
          N_case = 1000, alpha = 0.05,
          input_mode = "model_based",
          prev = 0.05, pd = 0.30, R2 = 1.8, MOI = MOI,
          locus_het = TRUE, pi = pi,
          k = 1, w = cc_ngs_het_scores(MOI),
          verbose = FALSE
        )
        expect_equal(
          observed$g_case_after_locus_het,
          canonical_public$freqs$g_true_case,
          tolerance = 1e-15
        )
        expect_equal(
          observed$g_ctrl_after_locus_het,
          canonical_public$freqs$g_true_ctrl,
          tolerance = 1e-15
        )
      }
    }
  }
})

test_that("the homogeneous boundary preserves pre-S12 power results", {
  for (MOI in c("M", "D", "Rec")) {
    baseline <- do.call(cc_ngs_power, cc_ngs_het_power_args(MOI = MOI))
    explicit <- do.call(
      cc_ngs_power,
      cc_ngs_het_power_args(MOI = MOI, locus_het = TRUE, pi = 1)
    )

    expect_equal(explicit$lambda, baseline$lambda, tolerance = 1e-15)
    expect_equal(explicit$power, baseline$power, tolerance = 1e-15)
    expect_equal(explicit$freqs$case_true, baseline$freqs$case_true,
                 tolerance = 1e-15)
    expect_equal(explicit$freqs$control_true, baseline$freqs$control_true,
                 tolerance = 1e-15)
    expect_equal(explicit$freqs$case_called, baseline$freqs$case_called,
                 tolerance = 1e-15)
    expect_equal(explicit$freqs$control_called, baseline$freqs$control_called,
                 tolerance = 1e-15)
    expect_equal(explicit$transition_matrix, baseline$transition_matrix,
                 tolerance = 1e-15)
  }
})

test_that("the homogeneous boundary preserves pre-S12 MSSN results", {
  for (MOI in c("M", "D", "Rec")) {
    baseline <- do.call(cc_ngs_mssn, cc_ngs_het_mssn_args(MOI = MOI))
    explicit <- do.call(
      cc_ngs_mssn,
      cc_ngs_het_mssn_args(MOI = MOI, locus_het = TRUE, pi = 1)
    )

    expect_equal(explicit$N_case_continuous, baseline$N_case_continuous,
                 tolerance = 1e-13)
    expect_identical(explicit$MSSN_case, baseline$MSSN_case)
    expect_identical(explicit$MSSN_ctrl, baseline$MSSN_ctrl)
    expect_identical(explicit$MSSN_total, baseline$MSSN_total)
    expect_equal(explicit$achieved_lambda, baseline$achieved_lambda,
                 tolerance = 1e-13)
    expect_equal(explicit$achieved_power, baseline$achieved_power,
                 tolerance = 1e-14)
    expect_equal(explicit$freqs$case_called, baseline$freqs$case_called,
                 tolerance = 1e-15)
  }
})

test_that("CC-NGS locus switch preserves fixed valid historical fixtures", {
  power_off <- do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(locus_het = FALSE, pi = 1)
  )
  power_one <- do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(locus_het = TRUE, pi = 1)
  )
  power_half <- do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(locus_het = TRUE, pi = 0.5)
  )
  mssn_off <- do.call(
    cc_ngs_mssn,
    cc_ngs_het_mssn_args(locus_het = FALSE, pi = 1)
  )
  mssn_one <- do.call(
    cc_ngs_mssn,
    cc_ngs_het_mssn_args(locus_het = TRUE, pi = 1)
  )
  mssn_half <- do.call(
    cc_ngs_mssn,
    cc_ngs_het_mssn_args(locus_het = TRUE, pi = 0.5)
  )

  expect_equal(power_off$lambda, 21.0597925597, tolerance = 1e-12)
  expect_equal(power_off$power, 0.995719830844489, tolerance = 1e-12)
  expect_equal(power_one$lambda, power_off$lambda, tolerance = 1e-15)
  expect_equal(power_one$power, power_off$power, tolerance = 1e-15)
  expect_equal(power_half$lambda, 5.42155552569551, tolerance = 1e-12)
  expect_equal(power_half$power, 0.643743647669686, tolerance = 1e-12)
  expect_identical(mssn_off$MSSN_case, 373)
  expect_identical(mssn_one$MSSN_case, 373)
  expect_identical(mssn_half$MSSN_case, 1448)
})

test_that("CC-NGS rejects pi when the locus switch is disabled", {
  message <- paste(
    "pi is used only when locus_het = TRUE;",
    "set pi = 1 or enable locus heterogeneity."
  )
  for (pi in c(0, 0.5)) {
    expect_error(
      do.call(
        cc_ngs_power,
        cc_ngs_het_power_args(locus_het = FALSE, pi = pi)
      ),
      message, fixed = TRUE
    )
    expect_error(
      do.call(
        cc_ngs_mssn,
        cc_ngs_het_mssn_args(locus_het = FALSE, pi = pi)
      ),
      message, fixed = TRUE
    )
  }
})

test_that("complete locus heterogeneity removes the CC contrast", {
  for (MOI in c("M", "D", "Rec")) {
    result <- do.call(
      cc_ngs_power,
      cc_ngs_het_power_args(
        MOI = MOI, alpha = 5e-8, locus_het = TRUE, pi = 0
      )
    )

    expect_equal(result$freqs$case_true, result$freqs$control_true,
                 tolerance = 1e-15)
    expect_equal(result$freqs$case_called, result$freqs$control_called,
                 tolerance = 1e-15)
    expect_equal(result$lambda, 0, tolerance = 1e-15)
    expect_lt(abs(result$power - result$alpha), 1e-15)

    expect_error(
      do.call(
        cc_ngs_mssn,
        cc_ngs_het_mssn_args(
          MOI = MOI, alpha = 5e-8, locus_het = TRUE, pi = 0
        )
      ),
      "No finite MSSN exists because the trend contrast is zero"
    )
  }
})

test_that("intermediate locus heterogeneity is the exact convex mixture", {
  model <- .cc_model_genotype_frequencies(
    pd = 0.25, R2 = 2, MOI = "M", prev = 0.04
  )
  for (pi in c(0.25, 0.5, 0.75)) {
    observed <- .cc_ngs_apply_locus_heterogeneity(
      model$case, model$control, locus_het = TRUE, pi = pi
    )
    expected <- pi * model$case + (1 - pi) * model$control
    expect_equal(observed$g_case_after_locus_het, expected,
                 tolerance = 1e-15)
    expect_equal(observed$g_ctrl_after_locus_het, model$control,
                 tolerance = 1e-15)
    expect_equal(sum(observed$g_case_after_locus_het), 1,
                 tolerance = 1e-15)
  }
})

test_that("heterogeneity and sequencing commute under one common matrix", {
  # This identity depends on the same nondifferential transition matrix E
  # being applied to cases and controls. It must not be generalized to future
  # differential case/control sequencing-error models.
  genotype_pairs <- list(
    list(case = c(0.45, 0.40, 0.15), control = c(0.55, 0.35, 0.10)),
    list(case = c(0.20, 0.50, 0.30), control = c(0.65, 0.30, 0.05))
  )
  for (pair in genotype_pairs) {
    for (coverage in c(2L, 4L, 20L, 100L)) {
      for (seq_error in c(0, 0.005, 0.01)) {
        E <- ngs_genotype_error_matrix(coverage, seq_error)
        case_called <- as.numeric(t(E) %*% pair$case)
        control_called <- as.numeric(t(E) %*% pair$control)
        for (pi in c(0, 0.25, 0.5, 0.75, 1)) {
          before <- as.numeric(
            t(E) %*% (pi * pair$case + (1 - pi) * pair$control)
          )
          after <- pi * case_called + (1 - pi) * control_called
          expect_equal(before, after, tolerance = 2e-15)
        }
      }
    }
  }
})

test_that("power and NCP attenuate toward alpha with dilution", {
  designs <- expand.grid(
    MOI = c("M", "D", "Rec"),
    coverage = c(4L, 20L),
    seq_error = c(0, 0.01),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(designs))) {
    pi_values <- c(0, 0.25, 0.5, 0.75, 1)
    results <- lapply(pi_values, function(pi) {
      do.call(
        cc_ngs_power,
        cc_ngs_het_power_args(
          MOI = designs$MOI[[i]],
          coverage = designs$coverage[[i]],
          seq_error = designs$seq_error[[i]],
          locus_het = TRUE,
          pi = pi
        )
      )
    })
    lambda <- vapply(results, `[[`, numeric(1), "lambda")
    power <- vapply(results, `[[`, numeric(1), "power")
    expect_true(all(diff(lambda) >= -1e-13))
    expect_true(all(diff(power) >= -1e-14))
    expect_equal(lambda[[1]], 0, tolerance = 1e-15)
    expect_equal(power[[1]], 0.05, tolerance = 2e-15)
  }
})

test_that("locus dilution inflates analytic MSSN", {
  designs <- list(
    list(MOI = "M", coverage = 4L, alpha = 0.05),
    list(MOI = "D", coverage = 100L, alpha = 5e-8),
    list(MOI = "Rec", coverage = 20L, alpha = 0.05)
  )
  for (design in designs) {
    results <- lapply(c(0.25, 0.5, 0.75, 1), function(pi) {
      do.call(
        cc_ngs_mssn,
        cc_ngs_het_mssn_args(
          MOI = design$MOI,
          coverage = design$coverage,
          alpha = design$alpha,
          locus_het = TRUE,
          pi = pi
        )
      )
    })
    continuous <- vapply(
      results, `[[`, numeric(1), "N_case_continuous"
    )
    integer <- vapply(results, `[[`, numeric(1), "MSSN_case")
    expect_true(all(diff(continuous) <= 1e-10))
    expect_true(all(diff(integer) <= 0))
  }
})

test_that("heterogeneous MSSN retains the exact Ahn analytic identity", {
  result <- do.call(
    cc_ngs_mssn,
    cc_ngs_het_mssn_args(
      power = 0.90, alpha = 5e-8,
      locus_het = TRUE, pi = 0.6, k = 1.5
    )
  )
  g_case <- result$freqs$case_called
  g_control <- result$freqs$control_called
  scores <- result$scores
  D <- sum(scores * (g_case - g_control))
  pooled <- g_case + result$k * g_control
  Q <- sum(scores^2 * pooled) -
    sum(scores * pooled)^2 / (1 + result$k)
  expected <- result$lambda_target * Q / (result$k * D^2)

  expect_equal(result$N_case_continuous, expected, tolerance = 2e-10)
})

test_that("heterogeneous integer MSSNs attain target and are minimal", {
  designs <- expand.grid(
    MOI = c("M", "D", "Rec"),
    coverage = c(4L, 100L),
    alpha = c(0.05, 5e-8),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(designs))) {
    result <- do.call(
      cc_ngs_mssn,
      cc_ngs_het_mssn_args(
        MOI = designs$MOI[[i]],
        coverage = designs$coverage[[i]],
        alpha = designs$alpha[[i]],
        locus_het = TRUE,
        pi = 0.6
      )
    )
    expect_gte(
      cc_ngs_het_achieved_power(result, result$MSSN_case),
      result$power_target - 1e-12
    )
    if (result$MSSN_case > 1) {
      expect_lt(
        cc_ngs_het_achieved_power(result, result$MSSN_case - 1),
        result$power_target
      )
    }
  }
})

test_that("high-depth heterogeneous CC-NGS converges to canonical CC", {
  for (MOI in c("M", "D", "Rec")) {
    scores <- cc_ngs_het_scores(MOI)
    ngs_power <- do.call(
      cc_ngs_power,
      cc_ngs_het_power_args(
        MOI = MOI, coverage = 500, seq_error = 0.01,
        locus_het = TRUE, pi = 0.5
      )
    )
    canonical_power <- cc_power(
      N_case = 1000, alpha = 0.05,
      input_mode = "model_based",
      prev = 0.05, pd = 0.30, R2 = 1.8, MOI = MOI,
      locus_het = TRUE, pi = 0.5,
      k = 1, w = scores, verbose = FALSE
    )
    ngs_mssn <- do.call(
      cc_ngs_mssn,
      cc_ngs_het_mssn_args(
        MOI = MOI, coverage = 500, seq_error = 0.01,
        locus_het = TRUE, pi = 0.5
      )
    )
    canonical_mssn <- cc_mssn(
      power = 0.80, alpha = 0.05,
      input_mode = "model_based",
      prev = 0.05, pd = 0.30, R2 = 1.8, MOI = MOI,
      locus_het = TRUE, pi = 0.5,
      k = 1, w = scores, verbose = FALSE
    )

    expect_equal(
      ngs_power$freqs$case_called,
      canonical_power$freqs$g_true_case,
      tolerance = 2e-12
    )
    expect_equal(ngs_power$power, canonical_power$tests$trend$power,
                 tolerance = 2e-11)
    expect_equal(
      ngs_mssn$MSSN_case,
      canonical_mssn$tests$trend$MSSN_case
    )
  }
})

test_that("finite-depth sequencing and locus dilution coexist", {
  settings <- list(
    homogeneous_high = list(coverage = 100, pi = 1),
    heterogeneous_high = list(coverage = 100, pi = 0.5),
    homogeneous_low = list(coverage = 2, pi = 1),
    heterogeneous_low = list(coverage = 2, pi = 0.5)
  )
  results <- lapply(settings, function(setting) {
    do.call(
      cc_ngs_power,
      cc_ngs_het_power_args(
        coverage = setting$coverage,
        seq_error = 0.01,
        locus_het = TRUE,
        pi = setting$pi
      )
    )
  })
  powers <- vapply(results, `[[`, numeric(1), "power")
  expect_true(all(is.finite(powers)))
  expect_gte(powers[["homogeneous_high"]], powers[["heterogeneous_high"]])
  expect_gte(powers[["homogeneous_high"]], powers[["homogeneous_low"]])
  expect_lte(powers[["heterogeneous_low"]], powers[["heterogeneous_high"]])
  expect_lte(powers[["heterogeneous_low"]], powers[["homogeneous_low"]])
})

test_that("zero read error retains finite-depth uncertainty with heterogeneity", {
  low <- do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(
      coverage = 2, seq_error = 0, locus_het = TRUE, pi = 0.6
    )
  )
  high <- do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(
      coverage = 100, seq_error = 0, locus_het = TRUE, pi = 0.6
    )
  )
  expected_low_case <- as.numeric(
    t(low$transition_matrix) %*% low$freqs$case_true
  )

  expect_equal(low$freqs$case_called, expected_low_case, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(
    low$freqs$case_called, low$freqs$case_true, tolerance = 1e-12
  )))
  expect_lt(
    max(abs(high$freqs$case_called - high$freqs$case_true)),
    max(abs(low$freqs$case_called - low$freqs$case_true))
  )
})

test_that("heterogeneous MSSN respects allocation ratios and minimality", {
  for (k in c(0.5, 1, 2)) {
    result <- do.call(
      cc_ngs_mssn,
      cc_ngs_het_mssn_args(k = k, locus_het = TRUE, pi = 0.5)
    )
    expect_equal(result$MSSN_ctrl, ceiling(k * result$MSSN_case))
    expect_equal(result$MSSN_total, result$MSSN_case + result$MSSN_ctrl)
    expect_gte(
      cc_ngs_het_achieved_power(result, result$MSSN_case),
      result$power_target - 1e-12
    )
    if (result$MSSN_case > 1) {
      expect_lt(
        cc_ngs_het_achieved_power(result, result$MSSN_case - 1),
        result$power_target
      )
    }
  }
})

test_that("heterogeneity metadata preserves old CC-NGS result fields", {
  power <- do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(locus_het = TRUE, pi = 0.6)
  )
  mssn <- do.call(
    cc_ngs_mssn,
    cc_ngs_het_mssn_args(locus_het = TRUE, pi = 0.6)
  )
  for (result in list(power, mssn)) {
    expect_true(all(c(
      "locus_het", "model_info", "freqs", "transition_matrix"
    ) %in% names(result)))
    expect_identical(result$locus_het$enabled, TRUE)
    expect_identical(result$locus_het$pi, 0.6)
    expect_equal(
      result$freqs$case_true_pre_heterogeneity,
      result$locus_het$g_case_before_locus_het,
      tolerance = 1e-15
    )
    expect_equal(
      result$freqs$case_true,
      result$locus_het$g_case_after_locus_het,
      tolerance = 1e-15
    )
    for (g in result$freqs) {
      expect_true(all(is.finite(g)))
      expect_true(all(g >= 0 & g <= 1))
      expect_equal(sum(g), 1, tolerance = 1e-14)
    }
    expect_equal(unname(rowSums(result$transition_matrix)), rep(1, 3),
                 tolerance = 1e-14)
  }
  expect_s3_class(power, "cc_ngs_power", exact = TRUE)
  expect_s3_class(mssn, "cc_ngs_mssn", exact = TRUE)
})

test_that("heterogeneity reporting is detailed and optional", {
  expect_output(
    do.call(
      cc_ngs_power,
      cc_ngs_het_power_args(locus_het = TRUE, pi = 0.6, verbose = FALSE)
    ),
    NA
  )
  expect_output(
    do.call(
      cc_ngs_mssn,
      cc_ngs_het_mssn_args(locus_het = TRUE, pi = 0.6, verbose = FALSE)
    ),
    NA
  )
  power_output <- capture.output(do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(locus_het = TRUE, pi = 0.6, verbose = TRUE)
  ), type = "message")
  mssn_output <- capture.output(do.call(
    cc_ngs_mssn,
    cc_ngs_het_mssn_args(locus_het = TRUE, pi = 0.6, verbose = TRUE)
  ), type = "message")
  default_output <- capture.output(do.call(
    cc_ngs_power,
    cc_ngs_het_power_args(verbose = TRUE)
  ), type = "message")

  expect_true(any(grepl("Locus heterogeneity:", power_output,
                         fixed = TRUE)))
  expect_true(any(grepl("Locus-homogeneity fraction (pi):", power_output,
                         fixed = TRUE)))
  expect_true(any(grepl("Locus heterogeneity:", mssn_output,
                         fixed = TRUE)))
  expect_false(any(grepl("Locus heterogeneity", default_output,
                          fixed = TRUE)))
  expect_false(any(grepl("transition_matrix", power_output, fixed = TRUE)))
})

test_that("public locus-heterogeneity inputs are validated", {
  invalid_pi <- list(-0.1, 1.1, NA_real_, NaN, Inf, "0.5", c(0.5, 0.6))
  invalid_flag <- list(NA, 0, 1, "TRUE", c(TRUE, FALSE))
  for (value in invalid_pi) {
    expect_error(
      do.call(cc_ngs_power, cc_ngs_het_power_args(pi = value)),
      "pi"
    )
    expect_error(
      do.call(cc_ngs_mssn, cc_ngs_het_mssn_args(pi = value)),
      "pi"
    )
  }
  for (value in invalid_flag) {
    expect_error(
      do.call(cc_ngs_power, cc_ngs_het_power_args(locus_het = value)),
      "locus_het"
    )
    expect_error(
      do.call(cc_ngs_mssn, cc_ngs_het_mssn_args(locus_het = value)),
      "locus_het"
    )
  }
})

test_that("a deterministic heterogeneity grid remains valid", {
  grid <- expand.grid(
    MOI = c("M", "D", "Rec"),
    coverage = c(2L, 20L, 100L),
    seq_error = c(0, 0.01),
    pi = c(0.25, 0.75, 1),
    alpha = c(0.05, 5e-8),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(grid))) {
    args <- cc_ngs_het_power_args(
      MOI = grid$MOI[[i]],
      coverage = grid$coverage[[i]],
      seq_error = grid$seq_error[[i]],
      pi = grid$pi[[i]],
      alpha = grid$alpha[[i]],
      locus_het = TRUE
    )
    first <- do.call(cc_ngs_power, args)
    second <- do.call(cc_ngs_power, args)
    expect_true(is.finite(first$lambda) && first$lambda >= 0)
    expect_true(is.finite(first$power) && first$power >= 0 && first$power <= 1)
    expect_equal(first$lambda, second$lambda, tolerance = 1e-14)
    expect_equal(first$power, second$power, tolerance = 1e-15)
    expect_equal(sum(first$freqs$case_true), 1, tolerance = 1e-14)
    expect_equal(sum(first$freqs$case_called), 1, tolerance = 1e-14)
  }
})
