# Private analytic bridge from sequencing-derived genotype calls to the
# Ahn/Chapman-Nam Cochran-Armitage trend-test NCP. This is not the raw-read
# LTTae,NGS likelihood/EM method.

.cc_ngs_validate_genotype_frequencies <- function(g, name) {
  if (!is.numeric(g) || length(g) != 3L) {
    stop(name, " must be a numeric vector of length 3.")
  }
  if (any(!is.finite(g))) {
    stop(name, " must contain only finite values.")
  }
  if (any(g < 0) || any(g > 1)) {
    stop(name, " must contain probabilities in [0, 1].")
  }
  if (abs(sum(g) - 1) > 1e-6) {
    stop(name, " must sum to 1.")
  }

  invisible(TRUE)
}

.cc_ahn_trend_ncp <- function(g_case, g_control, N_case, N_control,
                              scores = c(0, 1, 2)) {
  .cc_ngs_validate_genotype_frequencies(g_case, "g_case")
  .cc_ngs_validate_genotype_frequencies(g_control, "g_control")

  if (!is.numeric(N_case) || length(N_case) != 1L ||
      !is.finite(N_case) || N_case <= 0) {
    stop("N_case must be a single finite positive number.")
  }
  if (!is.numeric(N_control) || length(N_control) != 1L ||
      !is.finite(N_control) || N_control <= 0) {
    stop("N_control must be a single finite positive number.")
  }
  if (!is.numeric(scores) || length(scores) != 3L ||
      any(!is.finite(scores))) {
    stop("scores must be a finite numeric vector of length 3.")
  }
  if (length(unique(scores)) == 1L) {
    stop("scores cannot all be equal.")
  }

  contrast <- sum(scores * (g_case - g_control))
  weighted_counts <- N_case * g_case + N_control * g_control
  denominator <- sum(scores^2 * weighted_counts) -
    sum(scores * weighted_counts)^2 / (N_case + N_control)

  if (!is.finite(denominator) || denominator <= 0) {
    stop("Ahn trend-test variance denominator must be finite and positive.")
  }
  if (contrast == 0) {
    return(0)
  }

  lambda <- N_case * N_control * contrast^2 / denominator
  if (!is.finite(lambda) || lambda < 0) {
    stop("Ahn trend-test NCP must be finite and nonnegative.")
  }

  as.numeric(lambda)
}

.cc_ngs_called_frequencies <- function(g_case, g_control, coverage,
                                       seq_error) {
  .cc_ngs_validate_genotype_frequencies(g_case, "g_case")
  .cc_ngs_validate_genotype_frequencies(g_control, "g_control")

  E <- ngs_genotype_error_matrix(
    coverage = coverage,
    seq_error = seq_error
  )
  E <- .validate_genotype_misclassification_matrix(E, tolerance = 1e-12)
  case_called <- as.numeric(t(E) %*% g_case)
  control_called <- as.numeric(t(E) %*% g_control)

  .cc_ngs_validate_genotype_frequencies(case_called, "case_called")
  .cc_ngs_validate_genotype_frequencies(control_called, "control_called")

  list(
    E = E,
    case_true = as.numeric(g_case),
    control_true = as.numeric(g_control),
    case_called = case_called,
    control_called = control_called
  )
}

.cc_ngs_ahn_ncp <- function(g_case, g_control, N_case, N_control, coverage,
                            seq_error, scores = c(0, 1, 2)) {
  called <- .cc_ngs_called_frequencies(
    g_case = g_case,
    g_control = g_control,
    coverage = coverage,
    seq_error = seq_error
  )
  lambda <- .cc_ahn_trend_ncp(
    g_case = called$case_called,
    g_control = called$control_called,
    N_case = N_case,
    N_control = N_control,
    scores = scores
  )

  list(
    lambda = lambda,
    E = called$E,
    case_true = called$case_true,
    control_true = called$control_true,
    case_called = called$case_called,
    control_called = called$control_called,
    scores = as.numeric(scores),
    N_case = N_case,
    N_control = N_control
  )
}
