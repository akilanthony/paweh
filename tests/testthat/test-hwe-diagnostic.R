test_that("Pearson HWE fits the counted-allele frequency from the sample", {
  counts <- c("0" = 40, "1" = 30, "2" = 30)
  result <- hwe_test(counts, method = "chisq")
  expected <- c("0" = 30.25, "1" = 49.5, "2" = 20.25)
  statistic <- sum((counts - expected)^2 / expected)

  expect_equal(result$N, 100)
  expect_equal(result$observed_probs, counts / 100)
  expect_equal(result$allele_frequency, .45)
  expect_equal(result$expected_counts, expected)
  expect_equal(result$statistic, statistic)
  expect_identical(result$df, 1L)
  expect_equal(result$p_value, pchisq(statistic, df = 1, lower.tail = FALSE))
  expect_equal(result$F, 1 - 30 / 49.5)
  expect_false(result$sparse_expected)

  reversed <- hwe_test(unname(rev(counts)), method = "chisq")
  expect_equal(reversed$allele_frequency, .55)
  expect_equal(reversed$p_value, result$p_value)
  expect_identical(hwe_test(c(36, 48, 16), method = "chisq")$p_value, 1)
})

test_that("exact HWE probabilities agree with direct conditional enumeration", {
  conditional_p <- function(counts) {
    N <- sum(counts)
    counted <- counts[2] + 2 * counts[3]
    h <- seq(counted %% 2, min(counted, 2 * N - counted), by = 2)
    n2 <- (counted - h) / 2
    n0 <- N - h - n2
    # Conditional weights from the multinomial likelihood at p = 1/2.
    weights <- vapply(seq_along(h), function(i) {
      dmultinom(c(n0[i], h[i], n2[i]), prob = c(.25, .50, .25))
    }, numeric(1))
    observed <- weights[h == counts[2]]
    sum(weights[weights <= observed * (1 + 1e-12)]) / sum(weights)
  }

  cases <- expand.grid(n0 = 0:5, h = 0:5, n2 = 0:5)
  segregating <- with(cases, h + n0 > 0 & h + n2 > 0)
  cases <- as.matrix(cases[segregating, ])
  actual <- apply(cases, 1, function(x) hwe_test(unname(x))$p_value)
  expected <- apply(cases, 1, conditional_p)
  expect_equal(actual, expected, tolerance = 1e-13)

  # Rohlfs and Weir (2008), Table 2: N = 100, 14 counted alleles.
  example <- hwe_test(c(88, 10, 2))
  expect_equal(round(example$p_value, 4), .0654)
  expect_equal(round(example$statistic, 2), 5.38)
  expect_equal(hwe_test(c(2, 4, 0))$p_value, 1) # conditional mode tie
  expect_equal(hwe_test(c(1e9, 1, 0))$p_value, 1)
})

test_that("sparse and monomorphic tests report interpretable results", {
  expect_warning(result <- hwe_test(c(98, 1, 1), method = "chisq"), "below 5")
  expect_true(result$sparse_expected)
  expect_true(is.finite(result$p_value))
  expect_no_warning(hwe_test(c(98, 1, 1)))
  expect_warning(exact <- hwe_test(c(0, 0, 25)), "monomorphic")
  expect_warning(asymptotic <- hwe_test(c(25, 0, 0), method = "chisq"), "monomorphic")
  expect_equal(exact$allele_frequency, 1)
  expect_equal(asymptotic$allele_frequency, 0)
  expect_equal(exact$p_value, 1)
  expect_true(is.na(exact$statistic) && is.na(exact$df) && is.na(exact$F))
  expect_true(is.na(asymptotic$p_value))
  expect_error(hwe_test(c(1e6, 0, 1e6)), "one million conditional tables")
})

test_that("notebook matrix maps true rows to observed columns", {
  M <- rbind(c(.985, .010, .005), c(.001, .998, .001), c(.005, .010, .985))
  dimnames(M) <- list(c("0", "1", "2"), c("0", "1", "2"))
  g <- c("0" = .35588, "1" = .48424, "2" = .15988)
  power <- hwe_power(1000, allele_frequency = .4, misclassification_matrix = M)
  direct <- hwe_power(1000, genotype_probs = g)

  expect_equal(power$true_genotype_probs, c("0" = .36, "1" = .48, "2" = .16))
  expect_equal(power$observed_probs, g)
  expect_equal(power$allele_frequency, .402)
  expect_equal(power$expected_probs, c("0" = .357604, "1" = .480792, "2" = .161604))
  expect_equal(power$expected_counts, power$expected_probs * 1000)
  contrast <- sum((g - power$expected_probs)^2 / power$expected_probs)
  expect_equal(power$lambda_per_person, contrast, tolerance = 1e-13)
  expect_equal(power$lambda, 1000 * contrast, tolerance = 1e-13)
  expect_equal(power$F, 1 - .48424 / .480792)
  expect_equal(power$power, pchisq(qchisq(.95, 1), 1, ncp = 1000 * contrast,
                                 lower.tail = FALSE), tolerance = 1e-11)
  expect_equal(power$power, .0559122397772, tolerance = 1e-11)
  expect_equal(power$power, direct$power)
  expect_error(hwe_power(1000, allele_frequency = .4, misclassification_matrix = t(M)),
               "rows must sum to 1")
})

test_that("MSSN is the first integer attaining the requested asymptotic power", {
  g <- c(.35588, .48424, .15988)
  result <- hwe_mssn(genotype_probs = g)
  expect_equal(result$N, 152612)
  expect_gte(result$power, .80)
  expect_lt(result$previous_power, .80)
  expect_equal(result$power, hwe_power(result$N, genotype_probs = g)$power)
  expect_equal(result$previous_power, hwe_power(result$N - 1, genotype_probs = g)$power)

  for (g in list(c(.40, .20, .40), c(.20, .60, .20), c(.68, .27, .05))) {
    target <- .90
    result <- suppressWarnings(hwe_mssn(target, alpha = .01, genotype_probs = g))
    expect_equal(result$N, ceiling(result$N))
    expect_gte(result$power, target)
    expect_lt(result$previous_power, target)
    expect_equal(result$power,
                 suppressWarnings(hwe_power(result$N, .01, g))$power)
  }
})

test_that("zero departure and allele-frequency boundaries have defined behavior", {
  null <- c(.36, .48, .16)
  expect_equal(hwe_power(2000, alpha = .15, genotype_probs = null)$power, .15)
  expect_equal(hwe_power(2000, genotype_probs = null)$lambda, 0)
  infinite <- hwe_mssn(genotype_probs = null)
  expect_equal(infinite$N, Inf)
  expect_true(is.na(infinite$power) && is.na(infinite$lambda))
  expect_warning(trivial <- hwe_mssn(.05, genotype_probs = null), "below 5")
  expect_equal(trivial$N, 1)
  expect_equal(trivial$power, .05)
  expect_true(is.na(trivial$previous_power))

  identity <- diag(3)
  expect_equal(hwe_power(2000, allele_frequency = .9, misclassification_matrix = identity)$lambda, 0)
  expect_error(hwe_power(100, allele_frequency = 1, misclassification_matrix = identity),
               "monomorphic")
  expect_error(hwe_mssn(genotype_probs = c(1, 0, 0)), "monomorphic")
  M <- rbind(c(.7, .2, .1), c(0, 1, 0), c(.1, .3, .6))
  edge <- hwe_power(1000, allele_frequency = 0, misclassification_matrix = M)
  expect_equal(unname(edge$observed_probs), M[1, ])
  expect_warning(rare <- hwe_power(100, genotype_probs = c(.999, .001, 0)), "below 5")
  expect_true(rare$sparse_expected)
  expect_lt(rare$F, 0)
  expect_error(hwe_mssn(genotype_probs = c(.25 + 1e-10, .50 - 2e-10, .25 + 1e-10)),
               "exceeds the supported")
  expect_equal(hwe_power(1e10, genotype_probs = c(.5, 0, .5))$power, 1)
})

test_that("invalid HWE counts, probabilities, scales and matrix labels are rejected", {
  for (x in list(c(1, 2), c(-1, 2, 3), c(1, NA, 3), c(Inf, 1, 1), matrix(1:3))) {
    expect_error(hwe_test(x), "three finite, nonnegative")
  }
  expect_error(hwe_test(c(0, 0, 0)), "whole numbers")
  expect_error(hwe_test(c(.5, 2, 3)), "whole numbers")
  expect_error(hwe_test(c("2" = 1, "1" = 2, "0" = 3)), "ordered")
  expect_error(hwe_power(100), "Supply genotype_probs")
  expect_error(hwe_power(100, genotype_probs = c(.4, .4, .4)), "sum to 1")
  expect_error(hwe_power(100, genotype_probs = c(1.1, -.1, 0)), "nonnegative")
  expect_error(hwe_power(100, genotype_probs = c(.25, .5, .25), allele_frequency = .5),
               "Supply genotype_probs")
  expect_error(hwe_power(1.5, genotype_probs = c(.2, .6, .2)), "whole number")
  expect_error(hwe_power(0, genotype_probs = c(.2, .6, .2)), "N must")
  expect_error(hwe_power(2^52 + 1, genotype_probs = c(.2, .6, .2)), "N must")
  expect_error(hwe_power(100, alpha = 1, genotype_probs = c(.2, .6, .2)), "alpha")
  expect_error(hwe_mssn(target_power = 0, genotype_probs = c(.2, .6, .2)), "target_power")

  invalid <- list(diag(2), matrix(NA_real_, 3, 3), matrix(Inf, 3, 3),
                  matrix(-.1, 3, 3), matrix(2, 3, 3))
  for (M in invalid) {
    expect_error(hwe_power(100, allele_frequency = .5, misclassification_matrix = M),
                 "3 by 3 matrix")
  }
  M <- diag(c(1, 0, 1))
  expect_error(hwe_power(100, allele_frequency = .5, misclassification_matrix = M),
               "rows must sum")
  expect_error(hwe_power(100, allele_frequency = -1, misclassification_matrix = diag(3)),
               "allele_frequency")
  dimnames(M) <- list(c("0", "1", "2"), c("2", "1", "0"))
  expect_error(hwe_power(100, allele_frequency = .5, misclassification_matrix = M),
               "columns.*ordered")
})
