#' Noncentrality Parameter for Prospective LRTae Designs
#'
#' Computes the analytic noncentrality parameter (NCP) for the likelihood-ratio
#' test allowing for error (LRTae) in a prospective case-control design. The
#' calculation incorporates fallible measurements and validation subsamples for
#' phenotype, genotype, or both. It does not analyze observed study data.
#'
#' @param true_control_genotype_frequencies Numeric vector of true genotype
#'   frequencies among controls.
#' @param true_case_genotype_frequencies Numeric vector of true genotype
#'   frequencies among cases, in the same category order as
#'   `true_control_genotype_frequencies`.
#' @param true_sampling_proportions Numeric vector of length two containing the
#'   true control and case sampling proportions, in that order.
#' @param phenotype_misclassification_matrix Numeric two-by-two matrix whose
#'   rows are true phenotype categories (control, case) and whose columns are
#'   observed phenotype categories (control, case).
#' @param genotype_misclassification_matrix Numeric square matrix whose rows
#'   are true genotype categories and whose columns are observed genotype
#'   categories, in the genotype-frequency order.
#' @param n_fallible Number sampled with fallible phenotype and genotype
#'   measurements only.
#' @param n_pheno_validated Number whose phenotype is validated while genotype
#'   is measured fallibly.
#' @param n_geno_validated Number whose genotype is validated while phenotype
#'   is measured fallibly.
#' @param n_both_validated Number whose phenotype and genotype are both
#'   validated.
#'
#' @return A single nonnegative numeric value: the LRTae NCP.
#' @export
#'
#' @examples
#' lrtae_ncp(
#'   true_control_genotype_frequencies = c(0.81, 0.18, 0.01),
#'   true_case_genotype_frequencies = c(0.64, 0.32, 0.04),
#'   true_sampling_proportions = c(0.5, 0.5),
#'   phenotype_misclassification_matrix = matrix(
#'     c(0.95, 0.05, 0.05, 0.95), nrow = 2, byrow = TRUE
#'   ),
#'   genotype_misclassification_matrix = matrix(
#'     c(0.950, 0.050, 0.000,
#'       0.025, 0.950, 0.025,
#'       0.000, 0.050, 0.950),
#'     nrow = 3, byrow = TRUE
#'   ),
#'   n_fallible = 700,
#'   n_pheno_validated = 100,
#'   n_geno_validated = 100,
#'   n_both_validated = 100
#' )
lrtae_ncp <- function(
    true_control_genotype_frequencies,
    true_case_genotype_frequencies,
    true_sampling_proportions,
    phenotype_misclassification_matrix,
    genotype_misclassification_matrix,
    n_fallible,
    n_pheno_validated,
    n_geno_validated,
    n_both_validated
) {
  p_control <- .lrtae_validate_probability_vector(
    true_control_genotype_frequencies,
    "true_control_genotype_frequencies",
    minimum_length = 2L,
    strictly_positive = TRUE
  )
  p_case <- .lrtae_validate_probability_vector(
    true_case_genotype_frequencies,
    "true_case_genotype_frequencies",
    minimum_length = 2L,
    strictly_positive = TRUE
  )
  if (length(p_control) != length(p_case)) {
    stop("True control and case genotype-frequency vectors must have the same length.")
  }

  q <- .lrtae_validate_probability_vector(
    true_sampling_proportions,
    "true_sampling_proportions",
    required_length = 2L,
    strictly_positive = TRUE
  )
  pi_matrix <- .lrtae_validate_probability_matrix(
    phenotype_misclassification_matrix,
    "phenotype_misclassification_matrix",
    required_dimensions = c(2L, 2L)
  )
  k <- length(p_control)
  xi_matrix <- .lrtae_validate_probability_matrix(
    genotype_misclassification_matrix,
    "genotype_misclassification_matrix",
    required_dimensions = c(k, k)
  )

  counts <- c(
    n_fallible = n_fallible,
    n_pheno_validated = n_pheno_validated,
    n_geno_validated = n_geno_validated,
    n_both_validated = n_both_validated
  )
  for (count_name in names(counts)) {
    .lrtae_validate_count(counts[[count_name]], count_name)
  }
  if (sum(counts) == 0) {
    return(0)
  }

  # Rows are true phenotype categories (control, case); columns are true
  # genotype categories. This matches the corrected i-prime/j-prime notation.
  p <- rbind(p_control, p_case)
  probabilities <- .lrtae_sampling_probabilities(p, q, pi_matrix, xi_matrix)
  information_full <- .lrtae_expected_information(
    p = p,
    q = q,
    pi_matrix = pi_matrix,
    xi_matrix = xi_matrix,
    probabilities = probabilities,
    counts = counts
  )

  transformation <- .lrtae_transformation_matrix(k)
  information <- crossprod(transformation, information_full %*% transformation)
  information <- (information + t(information)) / 2

  psi_indices <- seq_len(k - 1L)
  nuisance_indices <- k:(2L * k - 1L)
  i_psi_psi <- information[psi_indices, psi_indices, drop = FALSE]
  i_psi_nuisance <- information[psi_indices, nuisance_indices, drop = FALSE]
  i_nuisance_nuisance <- information[
    nuisance_indices, nuisance_indices, drop = FALSE
  ]

  nuisance_solution <- tryCatch(
    solve(i_nuisance_nuisance, t(i_psi_nuisance)),
    error = function(e) {
      stop(
        "The nuisance-parameter information matrix is singular for this design.",
        call. = FALSE
      )
    }
  )
  efficient_information <- i_psi_psi -
    i_psi_nuisance %*% nuisance_solution
  efficient_information <- (efficient_information + t(efficient_information)) / 2

  delta <- p_control[psi_indices] - p_case[psi_indices]
  ncp <- as.numeric(crossprod(delta, efficient_information %*% delta))
  scale <- max(1, max(abs(efficient_information)))
  if (ncp < 0 && abs(ncp) <= 1e-10 * scale) {
    ncp <- 0
  }
  if (!is.finite(ncp) || ncp < 0) {
    stop("The calculated LRTae NCP is not finite and nonnegative.")
  }

  ncp
}

.lrtae_validate_probability_vector <- function(
    x, name, minimum_length = NULL, required_length = NULL,
    strictly_positive = FALSE, tolerance = 1e-8
) {
  valid <- is.numeric(x) && is.null(dim(x)) && all(is.finite(x))
  if (!is.null(required_length)) {
    valid <- valid && length(x) == required_length
  }
  if (!is.null(minimum_length)) {
    valid <- valid && length(x) >= minimum_length
  }
  if (!valid) {
    length_text <- if (!is.null(required_length)) {
      paste0(" of length ", required_length)
    } else {
      paste0(" with length at least ", minimum_length)
    }
    stop(name, " must be a finite numeric probability vector", length_text, ".")
  }
  if (any(x < 0) || any(x > 1)) {
    stop(name, " must contain probabilities in [0, 1].")
  }
  if (strictly_positive && any(x <= 0)) {
    stop(name, " must contain strictly positive probabilities.")
  }
  if (abs(sum(x) - 1) > tolerance) {
    stop(name, " must sum to 1.")
  }
  as.numeric(x)
}

.lrtae_validate_probability_matrix <- function(
    x, name, required_dimensions, tolerance = 1e-8
) {
  if (!is.matrix(x) || !is.numeric(x) ||
      !identical(dim(x), as.integer(required_dimensions))) {
    stop(
      name, " must be a numeric ", required_dimensions[[1L]], " x ",
      required_dimensions[[2L]], " matrix."
    )
  }
  if (any(!is.finite(x)) || any(x < 0) || any(x > 1)) {
    stop(name, " must contain finite probabilities in [0, 1].")
  }
  if (any(abs(rowSums(x) - 1) > tolerance)) {
    stop(name, " rows must sum to 1.")
  }
  x
}

.lrtae_validate_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    stop(name, " must be a single finite nonnegative number.")
  }
  invisible(TRUE)
}

.lrtae_sampling_probabilities <- function(p, q, pi_matrix, xi_matrix) {
  k <- ncol(p)
  p0 <- matrix(0, nrow = 2L, ncol = k)
  p1 <- array(0, dim = c(2L, 2L, k))
  p2 <- array(0, dim = c(2L, k, k))

  for (true_pheno in seq_len(2L)) {
    for (observed_pheno in seq_len(2L)) {
      for (observed_geno in seq_len(k)) {
        p1[true_pheno, observed_pheno, observed_geno] <- sum(
          pi_matrix[true_pheno, observed_pheno] *
            xi_matrix[, observed_geno] * q[true_pheno] *
            p[true_pheno, ]
        )
      }
    }
  }
  for (observed_pheno in seq_len(2L)) {
    for (true_geno in seq_len(k)) {
      for (observed_geno in seq_len(k)) {
        p2[observed_pheno, true_geno, observed_geno] <- sum(
          pi_matrix[, observed_pheno] *
            xi_matrix[true_geno, observed_geno] * q * p[, true_geno]
        )
      }
    }
    for (observed_geno in seq_len(k)) {
      p0[observed_pheno, observed_geno] <- sum(
        p1[, observed_pheno, observed_geno]
      )
    }
  }

  if (any(!is.finite(c(p0, p1, p2))) || any(c(p0, p1, p2) < 0)) {
    stop("Corrected LRTae sampling probabilities are numerically invalid.")
  }
  list(p0 = p0, p1 = p1, p2 = p2)
}

.lrtae_expected_information <- function(
    p, q, pi_matrix, xi_matrix, probabilities, counts
) {
  k <- ncol(p)
  information <- matrix(0, nrow = 2L * k + 2L, ncol = 2L * k + 2L)
  p0 <- probabilities$p0
  p1 <- probabilities$p1
  p2 <- probabilities$p2

  # Corrected block (a): genotype-frequency by genotype-frequency.
  for (true_pheno in seq_len(2L)) {
    for (true_geno in seq_len(k)) {
      row_index <- (true_pheno - 1L) * k + true_geno
      for (other_pheno in seq_len(2L)) {
        for (other_geno in seq_len(k)) {
          column_index <- (other_pheno - 1L) * k + other_geno
          term0 <- term1 <- term2 <- term3 <- 0
          for (observed_pheno in seq_len(2L)) {
            for (observed_geno in seq_len(k)) {
              if (p0[observed_pheno, observed_geno] > 0) {
                term0 <- term0 +
                  (pi_matrix[true_pheno, observed_pheno] *
                     xi_matrix[true_geno, observed_geno] * q[true_pheno]) *
                  (pi_matrix[other_pheno, observed_pheno] *
                     xi_matrix[other_geno, observed_geno] * q[other_pheno]) /
                  p0[observed_pheno, observed_geno]
              }
              if (true_pheno == other_pheno &&
                  p1[true_pheno, observed_pheno, observed_geno] > 0) {
                term1 <- term1 +
                  (pi_matrix[true_pheno, observed_pheno] *
                     xi_matrix[true_geno, observed_geno] * q[true_pheno]) *
                  (pi_matrix[true_pheno, observed_pheno] *
                     xi_matrix[other_geno, observed_geno] * q[true_pheno]) /
                  p1[true_pheno, observed_pheno, observed_geno]
              }
              if (true_geno == other_geno &&
                  p2[observed_pheno, true_geno, observed_geno] > 0) {
                term2 <- term2 +
                  (pi_matrix[true_pheno, observed_pheno] *
                     xi_matrix[true_geno, observed_geno] * q[true_pheno]) *
                  (pi_matrix[other_pheno, observed_pheno] *
                     xi_matrix[true_geno, observed_geno] * q[other_pheno]) /
                  p2[observed_pheno, true_geno, observed_geno]
              }
            }
          }
          if (true_pheno == other_pheno && true_geno == other_geno) {
            term3 <- q[true_pheno] / p[true_pheno, true_geno]
          }
          information[row_index, column_index] <-
            counts[["n_fallible"]] * term0 +
            counts[["n_pheno_validated"]] * term1 +
            counts[["n_geno_validated"]] * term2 +
            counts[["n_both_validated"]] * term3
        }
      }
    }
  }

  # Corrected block (b): genotype-frequency by sampling proportion. The
  # phenotype-validated and both-validated likelihoods have zero cross block.
  for (true_pheno in seq_len(2L)) {
    for (true_geno in seq_len(k)) {
      row_index <- (true_pheno - 1L) * k + true_geno
      for (other_pheno in seq_len(2L)) {
        column_index <- 2L * k + other_pheno
        term0 <- term2 <- 0
        for (observed_pheno in seq_len(2L)) {
          for (observed_geno in seq_len(k)) {
            if (p0[observed_pheno, observed_geno] > 0) {
              other_derivative <- sum(
                pi_matrix[other_pheno, observed_pheno] *
                  xi_matrix[, observed_geno] * p[other_pheno, ]
              )
              term0 <- term0 +
                pi_matrix[true_pheno, observed_pheno] *
                xi_matrix[true_geno, observed_geno] * q[true_pheno] *
                other_derivative / p0[observed_pheno, observed_geno]
            }
            if (p2[observed_pheno, true_geno, observed_geno] > 0) {
              term2 <- term2 +
                pi_matrix[true_pheno, observed_pheno] *
                xi_matrix[true_geno, observed_geno] * q[true_pheno] *
                pi_matrix[other_pheno, observed_pheno] *
                xi_matrix[true_geno, observed_geno] *
                p[other_pheno, true_geno] /
                p2[observed_pheno, true_geno, observed_geno]
            }
          }
        }
        if (true_pheno == other_pheno) {
          term0 <- term0 - 1
          term2 <- term2 - 1
        }
        value <- counts[["n_fallible"]] * term0 +
          counts[["n_geno_validated"]] * term2
        information[row_index, column_index] <- value
        information[column_index, row_index] <- value
      }
    }
  }

  # Corrected block (c): sampling proportion by sampling proportion.
  for (true_pheno in seq_len(2L)) {
    for (other_pheno in seq_len(2L)) {
      row_index <- 2L * k + true_pheno
      column_index <- 2L * k + other_pheno
      term0 <- term2 <- 0
      for (observed_pheno in seq_len(2L)) {
        for (observed_geno in seq_len(k)) {
          if (p0[observed_pheno, observed_geno] > 0) {
            first_derivative <- sum(
              pi_matrix[true_pheno, observed_pheno] *
                xi_matrix[, observed_geno] * p[true_pheno, ]
            )
            second_derivative <- sum(
              pi_matrix[other_pheno, observed_pheno] *
                xi_matrix[, observed_geno] * p[other_pheno, ]
            )
            term0 <- term0 + first_derivative * second_derivative /
              p0[observed_pheno, observed_geno]
          }
          for (true_geno in seq_len(k)) {
            if (p2[observed_pheno, true_geno, observed_geno] > 0) {
              term2 <- term2 +
                (pi_matrix[true_pheno, observed_pheno] *
                   xi_matrix[true_geno, observed_geno] *
                   p[true_pheno, true_geno]) *
                (pi_matrix[other_pheno, observed_pheno] *
                   xi_matrix[true_geno, observed_geno] *
                   p[other_pheno, true_geno]) /
                p2[observed_pheno, true_geno, observed_geno]
            }
          }
        }
      }
      diagonal_term <- if (true_pheno == other_pheno) {
        (counts[["n_pheno_validated"]] +
           counts[["n_both_validated"]]) / q[true_pheno]
      } else {
        0
      }
      information[row_index, column_index] <-
        counts[["n_fallible"]] * term0 +
        counts[["n_geno_validated"]] * term2 + diagonal_term
    }
  }

  information
}

.lrtae_transformation_matrix <- function(k) {
  k_minus_one <- k - 1L
  transformation <- matrix(
    0, nrow = 2L * k + 2L, ncol = 2L * k - 1L
  )
  psi_columns <- seq_len(k_minus_one)
  genotype_nuisance_columns <- k:(2L * k - 2L)
  sampling_column <- 2L * k - 1L

  transformation[psi_columns, psi_columns] <- diag(k_minus_one)
  transformation[psi_columns, genotype_nuisance_columns] <- diag(k_minus_one)
  transformation[k, psi_columns] <- -1
  transformation[k, genotype_nuisance_columns] <- -1
  transformation[k + psi_columns, genotype_nuisance_columns] <- diag(k_minus_one)
  transformation[2L * k, genotype_nuisance_columns] <- -1
  transformation[2L * k + 1L, sampling_column] <- 1
  transformation[2L * k + 2L, sampling_column] <- -1

  transformation
}

#' Prospective Power for an LRTae Design
#'
#' Computes prospective power for the likelihood-ratio test allowing for error
#' (LRTae). The function uses [lrtae_ncp()] and the corrected LRTae information
#' matrix; it does not perform an observed-data hypothesis test.
#'
#' @inheritParams lrtae_ncp
#' @param alpha Significance level, strictly between zero and one.
#'
#' @details For `k` genotype categories, the asymptotic reference distribution
#' has `k - 1` degrees of freedom. The rejection threshold is the upper
#' `alpha` quantile of the central chi-square distribution, and power is the
#' corresponding upper-tail probability from a noncentral chi-square
#' distribution with the NCP returned by [lrtae_ncp()].
#'
#' @return An object of class `lrtae_power`, containing `power`, `ncp`, `df`,
#'   `alpha`, `critical_value`, the four sampling-tier counts, and `N_total`.
#' @importFrom stats pchisq qchisq
#' @export
#'
#' @examples
#' lrtae_power(
#'   true_control_genotype_frequencies = c(0.81, 0.18, 0.01),
#'   true_case_genotype_frequencies = c(0.64, 0.32, 0.04),
#'   true_sampling_proportions = c(0.5, 0.5),
#'   phenotype_misclassification_matrix = matrix(
#'     c(0.95, 0.05, 0.05, 0.95), 2, byrow = TRUE
#'   ),
#'   genotype_misclassification_matrix = matrix(
#'     c(0.950, 0.050, 0.000,
#'       0.025, 0.950, 0.025,
#'       0.000, 0.050, 0.950),
#'     3, byrow = TRUE
#'   ),
#'   n_fallible = 700,
#'   n_pheno_validated = 100,
#'   n_geno_validated = 100,
#'   n_both_validated = 100,
#'   alpha = 0.05
#' )
lrtae_power <- function(
    true_control_genotype_frequencies,
    true_case_genotype_frequencies,
    true_sampling_proportions,
    phenotype_misclassification_matrix,
    genotype_misclassification_matrix,
    n_fallible,
    n_pheno_validated,
    n_geno_validated,
    n_both_validated,
    alpha
) {
  .lrtae_validate_open_probability(alpha, "alpha")

  ncp <- lrtae_ncp(
    true_control_genotype_frequencies =
      true_control_genotype_frequencies,
    true_case_genotype_frequencies = true_case_genotype_frequencies,
    true_sampling_proportions = true_sampling_proportions,
    phenotype_misclassification_matrix =
      phenotype_misclassification_matrix,
    genotype_misclassification_matrix = genotype_misclassification_matrix,
    n_fallible = n_fallible,
    n_pheno_validated = n_pheno_validated,
    n_geno_validated = n_geno_validated,
    n_both_validated = n_both_validated
  )
  df <- length(true_control_genotype_frequencies) - 1L
  critical_value <- stats::qchisq(1 - alpha, df = df)
  power <- stats::pchisq(
    critical_value, df = df, ncp = ncp, lower.tail = FALSE
  )

  structure(
    list(
      power = as.numeric(power),
      ncp = ncp,
      df = df,
      alpha = alpha,
      critical_value = as.numeric(critical_value),
      n_fallible = n_fallible,
      n_pheno_validated = n_pheno_validated,
      n_geno_validated = n_geno_validated,
      n_both_validated = n_both_validated,
      N_total = n_fallible + n_pheno_validated + n_geno_validated +
        n_both_validated
    ),
    class = "lrtae_power"
  )
}

#' Minimum Total Sample Size for an LRTae Design
#'
#' Finds the minimum integer total sample size required to attain a requested
#' prospective LRTae power while preserving a requested allocation among the
#' four measurement and validation tiers as closely as integer counts permit.
#' This is a prospective-design calculation and is not an observed-data test.
#'
#' @inheritParams lrtae_ncp
#' @param target_power Requested power, strictly greater than `alpha` and less
#'   than one.
#' @param alpha Significance level, strictly between zero and one.
#' @param prop_fallible Proportion allocated to fallible phenotype and genotype
#'   measurements only.
#' @param prop_pheno_validated Proportion allocated to phenotype validation
#'   with fallible genotype measurement.
#' @param prop_geno_validated Proportion allocated to genotype validation with
#'   fallible phenotype measurement.
#' @param prop_both_validated Proportion allocated to validation of both
#'   phenotype and genotype.
#'
#' @details The four allocation proportions must be nonnegative and sum to one.
#' For every candidate total `N`, integer counts are assigned by the
#' largest-remainder (Hamilton) rule: floor `N * proportion` for each tier,
#' then distribute the remaining individuals in descending order of fractional
#' remainder. Exact ties are resolved in the fixed order: fallible, phenotype
#' validated, genotype validated, both validated. Thus allocated counts always
#' sum exactly to `N`.
#'
#' The search doubles an upper bound until it attains `target_power`, applies
#' integer binary search within the bracket, and explicitly verifies the
#' neighboring design at `N - 1`. Power at each candidate is evaluated through
#' [lrtae_power()], which uses the validated corrected [lrtae_ncp()] engine,
#' `k - 1` degrees of freedom, and a noncentral chi-square distribution.
#'
#' @return An object of class `lrtae_mssn`, containing `MSSN_total`,
#'   `target_power`, `achieved_power`, `achieved_ncp`, `df`, `alpha`,
#'   `critical_value`, requested allocation proportions, and the four allocated
#'   integer tier counts.
#' @export
#'
#' @examples
#' lrtae_mssn(
#'   true_control_genotype_frequencies = c(0.81, 0.18, 0.01),
#'   true_case_genotype_frequencies = c(0.64, 0.32, 0.04),
#'   true_sampling_proportions = c(0.5, 0.5),
#'   phenotype_misclassification_matrix = matrix(
#'     c(0.95, 0.05, 0.05, 0.95), 2, byrow = TRUE
#'   ),
#'   genotype_misclassification_matrix = matrix(
#'     c(0.950, 0.050, 0.000,
#'       0.025, 0.950, 0.025,
#'       0.000, 0.050, 0.950),
#'     3, byrow = TRUE
#'   ),
#'   target_power = 0.80,
#'   alpha = 0.05,
#'   prop_fallible = 0.7,
#'   prop_pheno_validated = 0.1,
#'   prop_geno_validated = 0.1,
#'   prop_both_validated = 0.1
#' )
lrtae_mssn <- function(
    true_control_genotype_frequencies,
    true_case_genotype_frequencies,
    true_sampling_proportions,
    phenotype_misclassification_matrix,
    genotype_misclassification_matrix,
    target_power,
    alpha,
    prop_fallible,
    prop_pheno_validated,
    prop_geno_validated,
    prop_both_validated
) {
  .lrtae_validate_open_probability(alpha, "alpha")
  .lrtae_validate_open_probability(target_power, "target_power")
  if (target_power <= alpha) {
    stop("target_power must be strictly greater than alpha.")
  }

  proportions <- c(
    n_fallible = prop_fallible,
    n_pheno_validated = prop_pheno_validated,
    n_geno_validated = prop_geno_validated,
    n_both_validated = prop_both_validated
  )
  .lrtae_validate_tier_proportions(proportions)
  proportions <- proportions / sum(proportions)

  # A fractional one-person design validates all shared model inputs and tests
  # whether this allocation has positive information under the alternative.
  ncp_per_person <- lrtae_ncp(
    true_control_genotype_frequencies =
      true_control_genotype_frequencies,
    true_case_genotype_frequencies = true_case_genotype_frequencies,
    true_sampling_proportions = true_sampling_proportions,
    phenotype_misclassification_matrix =
      phenotype_misclassification_matrix,
    genotype_misclassification_matrix = genotype_misclassification_matrix,
    n_fallible = proportions[["n_fallible"]],
    n_pheno_validated = proportions[["n_pheno_validated"]],
    n_geno_validated = proportions[["n_geno_validated"]],
    n_both_validated = proportions[["n_both_validated"]]
  )
  if (ncp_per_person <= 0) {
    stop(
      "target_power cannot be reached because this design has zero LRTae NCP.",
      call. = FALSE
    )
  }

  evaluate_total <- function(N) {
    allocated <- .lrtae_allocate_tiers(N, proportions)
    result <- tryCatch(
      lrtae_power(
        true_control_genotype_frequencies =
          true_control_genotype_frequencies,
        true_case_genotype_frequencies = true_case_genotype_frequencies,
        true_sampling_proportions = true_sampling_proportions,
        phenotype_misclassification_matrix =
          phenotype_misclassification_matrix,
        genotype_misclassification_matrix =
          genotype_misclassification_matrix,
        n_fallible = allocated[["n_fallible"]],
        n_pheno_validated = allocated[["n_pheno_validated"]],
        n_geno_validated = allocated[["n_geno_validated"]],
        n_both_validated = allocated[["n_both_validated"]],
        alpha = alpha
      ),
      error = function(e) {
        if (grepl("information matrix is singular", conditionMessage(e),
                  fixed = TRUE)) {
          return(NULL)
        }
        stop(e)
      }
    )
    list(
      N = N,
      counts = allocated,
      result = result,
      power = if (is.null(result)) -Inf else result$power
    )
  }

  lower <- 0L
  upper <- 1L
  upper_result <- evaluate_total(upper)
  while (upper_result$power < target_power) {
    lower <- upper
    if (upper > .Machine$integer.max %/% 2L) {
      stop("Unable to bracket the requested LRTae sample size.", call. = FALSE)
    }
    upper <- as.integer(upper * 2L)
    upper_result <- evaluate_total(upper)
  }

  while (lower + 1L < upper) {
    midpoint <- as.integer(floor((lower + upper) / 2))
    midpoint_result <- evaluate_total(midpoint)
    if (midpoint_result$power >= target_power) {
      upper <- midpoint
      upper_result <- midpoint_result
    } else {
      lower <- midpoint
    }
  }

  selected <- upper_result
  if (selected$N != upper) {
    selected <- evaluate_total(upper)
  }
  while (selected$N > 1L) {
    previous <- evaluate_total(selected$N - 1L)
    if (previous$power < target_power) {
      break
    }
    selected <- previous
  }
  while (selected$power < target_power) {
    selected <- evaluate_total(selected$N + 1L)
  }

  structure(
    list(
      MSSN_total = selected$N,
      target_power = target_power,
      achieved_power = selected$result$power,
      achieved_ncp = selected$result$ncp,
      df = selected$result$df,
      alpha = alpha,
      critical_value = selected$result$critical_value,
      prop_fallible = proportions[["n_fallible"]],
      prop_pheno_validated = proportions[["n_pheno_validated"]],
      prop_geno_validated = proportions[["n_geno_validated"]],
      prop_both_validated = proportions[["n_both_validated"]],
      n_fallible = selected$counts[["n_fallible"]],
      n_pheno_validated = selected$counts[["n_pheno_validated"]],
      n_geno_validated = selected$counts[["n_geno_validated"]],
      n_both_validated = selected$counts[["n_both_validated"]]
    ),
    class = "lrtae_mssn"
  )
}

.lrtae_validate_open_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop(name, " must be a single finite number strictly between 0 and 1.")
  }
  invisible(TRUE)
}

.lrtae_validate_tier_proportions <- function(proportions, tolerance = 1e-8) {
  if (!is.numeric(proportions) || length(proportions) != 4L ||
      any(!is.finite(proportions)) || any(proportions < 0) ||
      any(proportions > 1)) {
    stop("LRTae tier proportions must be four finite numbers in [0, 1].")
  }
  if (!any(proportions > 0)) {
    stop("At least one LRTae tier proportion must be positive.")
  }
  if (abs(sum(proportions) - 1) > tolerance) {
    stop("LRTae tier proportions must sum to 1.")
  }
  invisible(TRUE)
}

.lrtae_allocate_tiers <- function(N, proportions) {
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
      N < 0 || N != floor(N) || N > .Machine$integer.max) {
    stop("N must be a single nonnegative integer.")
  }
  .lrtae_validate_tier_proportions(proportions)
  proportions <- proportions / sum(proportions)

  raw_counts <- N * proportions
  allocated <- floor(raw_counts)
  remaining <- as.integer(N - sum(allocated))
  if (remaining > 0L) {
    fractional_remainders <- raw_counts - allocated
    fixed_tier_order <- seq_along(proportions)
    recipients <- order(-fractional_remainders, fixed_tier_order)[
      seq_len(remaining)
    ]
    allocated[recipients] <- allocated[recipients] + 1
  }
  storage.mode(allocated) <- "integer"
  allocated
}
