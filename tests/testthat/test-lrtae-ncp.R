lrtae_reference_arguments <- function() {
  list(
    true_control_genotype_frequencies = c(0.81, 0.18, 0.01),
    true_case_genotype_frequencies = c(0.64, 0.32, 0.04),
    true_sampling_proportions = c(0.5, 0.5),
    phenotype_misclassification_matrix = matrix(
      c(0.95, 0.05,
        0.05, 0.95),
      nrow = 2,
      byrow = TRUE
    ),
    genotype_misclassification_matrix = matrix(
      c(0.950, 0.050, 0.000,
        0.025, 0.950, 0.025,
        0.000, 0.050, 0.950),
      nrow = 3,
      byrow = TRUE
    ),
    n_fallible = 700,
    n_pheno_validated = 100,
    n_geno_validated = 100,
    n_both_validated = 100
  )
}

test_that("professor-supplied corrected LRTae example is reproduced", {
  reference_ncp <- 29.173050672171374
  observed_ncp <- do.call(lrtae_ncp, lrtae_reference_arguments())

  expect_equal(observed_ncp, reference_ncp, tolerance = 1e-12)
  expect_equal(abs(observed_ncp - reference_ncp), 0, tolerance = 1e-12)
})

test_that("the LRTae NCP is zero under equal genotype frequencies", {
  arguments <- lrtae_reference_arguments()
  arguments$true_case_genotype_frequencies <-
    arguments$true_control_genotype_frequencies

  expect_equal(do.call(lrtae_ncp, arguments), 0, tolerance = 1e-12)
})

test_that("the LRTae NCP scales linearly with every sampling-tier count", {
  arguments <- lrtae_reference_arguments()
  baseline <- do.call(lrtae_ncp, arguments)
  count_names <- c(
    "n_fallible", "n_pheno_validated", "n_geno_validated",
    "n_both_validated"
  )
  arguments[count_names] <- lapply(arguments[count_names], `*`, 2)

  expect_equal(do.call(lrtae_ncp, arguments), 2 * baseline, tolerance = 1e-11)
})

test_that("corrected information agrees with an independent score calculation", {
  arguments <- lrtae_reference_arguments()
  p <- rbind(
    arguments$true_control_genotype_frequencies,
    arguments$true_case_genotype_frequencies
  )
  q <- arguments$true_sampling_proportions
  pi_matrix <- arguments$phenotype_misclassification_matrix
  xi_matrix <- arguments$genotype_misclassification_matrix
  counts <- unlist(arguments[c(
    "n_fallible", "n_pheno_validated", "n_geno_validated",
    "n_both_validated"
  )])
  names(counts) <- c(
    "n_fallible", "n_pheno_validated", "n_geno_validated",
    "n_both_validated"
  )
  k <- ncol(p)

  probabilities <- .lrtae_sampling_probabilities(p, q, pi_matrix, xi_matrix)
  full_information <- .lrtae_expected_information(
    p, q, pi_matrix, xi_matrix, probabilities, counts
  )
  transformed_information <- crossprod(
    .lrtae_transformation_matrix(k),
    full_information %*% .lrtae_transformation_matrix(k)
  )

  free_parameters <- c(
    p[1, seq_len(k - 1)] - p[2, seq_len(k - 1)],
    p[2, seq_len(k - 1)], q[1]
  )
  tier_probabilities <- function(theta) {
    psi <- theta[seq_len(k - 1)]
    case_free <- theta[k:(2 * k - 2)]
    q_control <- theta[2 * k - 1]
    p_case <- c(case_free, 1 - sum(case_free))
    p_control <- c(psi + case_free, 1 - sum(psi + case_free))
    p_theta <- rbind(p_control, p_case)
    q_theta <- c(q_control, 1 - q_control)
    sampled <- .lrtae_sampling_probabilities(
      p_theta, q_theta, pi_matrix, xi_matrix
    )
    p3 <- array(0, dim = c(2, 2, k, k))
    for (true_pheno in 1:2) {
      for (observed_pheno in 1:2) {
        for (true_geno in seq_len(k)) {
          for (observed_geno in seq_len(k)) {
            p3[true_pheno, observed_pheno, true_geno, observed_geno] <-
              pi_matrix[true_pheno, observed_pheno] *
              xi_matrix[true_geno, observed_geno] *
              q_theta[true_pheno] * p_theta[true_pheno, true_geno]
          }
        }
      }
    }
    list(sampled$p0, sampled$p1, sampled$p2, p3)
  }

  step <- 1e-6
  base_probabilities <- tier_probabilities(free_parameters)
  numerical_information <- matrix(0, nrow = length(free_parameters),
                                  ncol = length(free_parameters))
  for (tier in seq_along(base_probabilities)) {
    jacobian <- vapply(seq_along(free_parameters), function(column) {
      upper <- lower <- free_parameters
      upper[column] <- upper[column] + step
      lower[column] <- lower[column] - step
      (tier_probabilities(upper)[[tier]] -
         tier_probabilities(lower)[[tier]]) / (2 * step)
    }, numeric(length(base_probabilities[[tier]])))
    positive <- as.numeric(base_probabilities[[tier]]) > 0
    numerical_information <- numerical_information + counts[[tier]] *
      crossprod(
        jacobian[positive, , drop = FALSE] /
          sqrt(as.numeric(base_probabilities[[tier]])[positive])
      )
  }

  expect_equal(transformed_information, numerical_information,
               tolerance = 2e-7)
})

test_that("malformed LRTae probability inputs are rejected", {
  arguments <- lrtae_reference_arguments()

  bad <- arguments
  bad$true_control_genotype_frequencies <- c(0.8, NA_real_, 0.2)
  expect_error(do.call(lrtae_ncp, bad), "finite numeric probability vector")

  bad <- arguments
  bad$true_sampling_proportions <- c(0.5, 0.4, 0.1)
  expect_error(do.call(lrtae_ncp, bad), "length 2")

  bad <- arguments
  bad$phenotype_misclassification_matrix[1, 1] <- 1.1
  expect_error(do.call(lrtae_ncp, bad), "probabilities in \\[0, 1\\]")

  bad <- arguments
  bad$genotype_misclassification_matrix[1, ] <- c(0.8, 0.1, 0.05)
  expect_error(do.call(lrtae_ncp, bad), "rows must sum to 1")
})

test_that("LRTae genotype frequencies must sum to one", {
  arguments <- lrtae_reference_arguments()
  arguments$true_case_genotype_frequencies <- c(0.64, 0.31, 0.04)

  expect_error(do.call(lrtae_ncp, arguments), "must sum to 1")
})

test_that("LRTae sampling-tier counts must be nonnegative", {
  count_names <- c(
    "n_fallible", "n_pheno_validated", "n_geno_validated",
    "n_both_validated"
  )
  for (count_name in count_names) {
    arguments <- lrtae_reference_arguments()
    arguments[[count_name]] <- -1
    expect_error(do.call(lrtae_ncp, arguments), "nonnegative")
  }
})

test_that("LRTae matrix dimensions must match the probability vectors", {
  arguments <- lrtae_reference_arguments()

  bad <- arguments
  bad$phenotype_misclassification_matrix <- diag(3)
  expect_error(do.call(lrtae_ncp, bad), "2 x 2 matrix")

  bad <- arguments
  bad$genotype_misclassification_matrix <- diag(2)
  expect_error(do.call(lrtae_ncp, bad), "3 x 3 matrix")

  bad <- arguments
  bad$true_case_genotype_frequencies <- c(0.8, 0.2)
  expect_error(do.call(lrtae_ncp, bad), "must have the same length")
})

test_that("professor example gives independently calculated LRTae power", {
  arguments <- lrtae_reference_arguments()
  alpha <- 0.05
  reference_ncp <- 29.173050672171374
  df <- length(arguments$true_control_genotype_frequencies) - 1L
  critical <- stats::qchisq(1 - alpha, df = df)
  reference_power <- stats::pchisq(
    critical, df = df, ncp = reference_ncp, lower.tail = FALSE
  )
  observed <- do.call(lrtae_power, c(arguments, list(alpha = alpha)))

  expect_s3_class(observed, "lrtae_power", exact = TRUE)
  expect_equal(observed$ncp, reference_ncp, tolerance = 1e-12)
  expect_identical(observed$df, df)
  expect_equal(observed$critical_value, critical, tolerance = 1e-14)
  expect_equal(observed$power, reference_power, tolerance = 1e-14)
  expect_equal(abs(observed$power - reference_power), 0, tolerance = 1e-14)
})

test_that("LRTae null power equals alpha", {
  arguments <- lrtae_reference_arguments()
  arguments$true_case_genotype_frequencies <-
    arguments$true_control_genotype_frequencies
  observed <- do.call(lrtae_power, c(arguments, list(alpha = 0.05)))

  expect_equal(observed$ncp, 0, tolerance = 1e-12)
  expect_equal(observed$power, observed$alpha, tolerance = 1e-14)
})

test_that("increasing all LRTae tier counts cannot reduce power", {
  arguments <- lrtae_reference_arguments()
  baseline <- do.call(lrtae_power, c(arguments, list(alpha = 0.05)))
  count_names <- c(
    "n_fallible", "n_pheno_validated", "n_geno_validated",
    "n_both_validated"
  )
  arguments[count_names] <- lapply(arguments[count_names], `*`, 2)
  increased <- do.call(lrtae_power, c(arguments, list(alpha = 0.05)))

  expect_gte(increased$power, baseline$power)
  expect_equal(increased$ncp, 2 * baseline$ncp, tolerance = 1e-11)
})

test_that("LRTae power validates alpha and delegates model validation", {
  arguments <- lrtae_reference_arguments()
  for (bad_alpha in c(0, 1, -0.1, NA_real_)) {
    expect_error(
      do.call(lrtae_power, c(arguments, list(alpha = bad_alpha))),
      "alpha.*strictly between 0 and 1"
    )
  }

  arguments$true_control_genotype_frequencies <- c(0.8, 0.1, 0.05)
  expect_error(
    do.call(lrtae_power, c(arguments, list(alpha = 0.05))),
    "true_control_genotype_frequencies must sum to 1"
  )
})

test_that("Hamilton tier allocation is exact and resolves ties by tier order", {
  proportions <- c(
    n_fallible = 0.25,
    n_pheno_validated = 0.25,
    n_geno_validated = 0.25,
    n_both_validated = 0.25
  )

  expect_identical(
    .lrtae_allocate_tiers(7, proportions),
    c(
      n_fallible = 2L,
      n_pheno_validated = 2L,
      n_geno_validated = 2L,
      n_both_validated = 1L
    )
  )
  expect_equal(sum(.lrtae_allocate_tiers(101, proportions)), 101)
})

test_that("LRTae MSSN attains target and N minus one does not", {
  model_arguments <- lrtae_reference_arguments()[1:5]
  allocation <- c(
    prop_fallible = 0.7,
    prop_pheno_validated = 0.1,
    prop_geno_validated = 0.1,
    prop_both_validated = 0.1
  )
  observed <- do.call(
    lrtae_mssn,
    c(model_arguments, list(target_power = 0.80, alpha = 0.05),
      as.list(allocation))
  )

  expect_s3_class(observed, "lrtae_mssn", exact = TRUE)
  expect_gte(observed$achieved_power, observed$target_power)
  returned_counts <- unlist(observed[c(
    "n_fallible", "n_pheno_validated", "n_geno_validated",
    "n_both_validated"
  )])
  expect_identical(sum(returned_counts), observed$MSSN_total)

  direct <- do.call(
    lrtae_power,
    c(model_arguments, as.list(returned_counts), list(alpha = observed$alpha))
  )
  expect_equal(observed$achieved_power, direct$power, tolerance = 1e-14)
  expect_equal(observed$achieved_ncp, direct$ncp, tolerance = 1e-14)

  previous_counts <- .lrtae_allocate_tiers(
    observed$MSSN_total - 1L,
    c(
      n_fallible = allocation[["prop_fallible"]],
      n_pheno_validated = allocation[["prop_pheno_validated"]],
      n_geno_validated = allocation[["prop_geno_validated"]],
      n_both_validated = allocation[["prop_both_validated"]]
    )
  )
  previous <- do.call(
    lrtae_power,
    c(model_arguments, as.list(previous_counts), list(alpha = observed$alpha))
  )
  expect_lt(previous$power, observed$target_power)
})

test_that("LRTae MSSN is monotone in target power and significance stringency", {
  model_arguments <- lrtae_reference_arguments()[1:5]
  design <- list(
    prop_fallible = 0.7,
    prop_pheno_validated = 0.1,
    prop_geno_validated = 0.1,
    prop_both_validated = 0.1
  )
  mssn_at <- function(target_power, alpha) {
    do.call(
      lrtae_mssn,
      c(model_arguments, list(target_power = target_power, alpha = alpha),
        design)
    )$MSSN_total
  }

  expect_gte(mssn_at(0.90, 0.05), mssn_at(0.80, 0.05))
  expect_gte(mssn_at(0.80, 0.01), mssn_at(0.80, 0.05))
})

test_that("LRTae MSSN validates target power and tier proportions", {
  model_arguments <- lrtae_reference_arguments()[1:5]
  valid_design <- list(
    prop_fallible = 0.7,
    prop_pheno_validated = 0.1,
    prop_geno_validated = 0.1,
    prop_both_validated = 0.1
  )
  call_mssn <- function(target_power = 0.80, alpha = 0.05,
                        design = valid_design) {
    do.call(
      lrtae_mssn,
      c(model_arguments, list(target_power = target_power, alpha = alpha),
        design)
    )
  }

  expect_error(call_mssn(target_power = 0.05),
               "target_power must be strictly greater than alpha")
  expect_error(call_mssn(target_power = 1),
               "target_power.*strictly between 0 and 1")

  negative <- valid_design
  negative$prop_fallible <- -0.1
  negative$prop_both_validated <- 0.9
  expect_error(call_mssn(design = negative), "four finite numbers")

  wrong_sum <- valid_design
  wrong_sum$prop_fallible <- 0.6
  expect_error(call_mssn(design = wrong_sum), "must sum to 1")
})

test_that("zero-allocation LRTae tiers are allowed when identifiable", {
  model_arguments <- lrtae_reference_arguments()[1:5]
  observed <- do.call(
    lrtae_mssn,
    c(
      model_arguments,
      list(
        target_power = 0.80,
        alpha = 0.05,
        prop_fallible = 0,
        prop_pheno_validated = 0,
        prop_geno_validated = 0,
        prop_both_validated = 1
      )
    )
  )

  expect_gt(observed$MSSN_total, 0)
  expect_identical(observed$n_fallible, 0L)
  expect_identical(observed$n_pheno_validated, 0L)
  expect_identical(observed$n_geno_validated, 0L)
  expect_identical(observed$n_both_validated, observed$MSSN_total)
  expect_identical(
    observed$n_fallible + observed$n_pheno_validated +
      observed$n_geno_validated + observed$n_both_validated,
    observed$MSSN_total
  )
})
