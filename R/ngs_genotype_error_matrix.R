.ngs_call_genotype_ml <- function(alt_reads, coverage, seq_error) {
  if (!is.numeric(coverage) || length(coverage) != 1L ||
      !is.finite(coverage) || coverage < 1 || coverage != floor(coverage)) {
    stop("coverage must be a single finite integer greater than or equal to 1.")
  }
  if (!is.numeric(alt_reads) || length(alt_reads) != 1L ||
      !is.finite(alt_reads) || alt_reads < 0 ||
      alt_reads > coverage || alt_reads != floor(alt_reads)) {
    stop("alt_reads must be a single finite integer between 0 and coverage.")
  }

  q <- vapply(
    0:2,
    .ngs_alt_read_probability,
    numeric(1),
    seq_error = seq_error
  )
  log_likelihood <- stats::dbinom(
    alt_reads,
    size = coverage,
    prob = q,
    log = TRUE
  )

  tie_tolerance <- 1e-12
  candidates <- which(log_likelihood >= max(log_likelihood) - tie_tolerance) - 1L

  if (length(candidates) == 1L) {
    return(candidates)
  }

  read_fraction <- alt_reads / coverage
  if (read_fraction < 0.5) {
    return(min(candidates))
  }
  if (read_fraction > 0.5) {
    return(max(candidates))
  }

  candidates[[which.min(abs(candidates - 1L))]]
}

#' Sequencing-Derived Genotype Error Matrix
#'
#' Constructs the analytic genotype transition matrix induced by fixed
#' sequencing coverage, symmetric per-read sequencing error, and deterministic
#' maximum-likelihood genotype calling.
#'
#' @param coverage A single finite integer greater than or equal to 1 giving
#'   sequencing coverage (read depth).
#' @param seq_error A single finite numeric sequencing-error probability in
#'   \code{[0, 0.5)}.
#'
#' @details
#' For true genotype \eqn{G \in \{0,1,2\}}, the alternate-read count follows
#' \deqn{X \mid G,V,\epsilon \sim \mathrm{Binomial}(V,q_G),}
#' where \eqn{q_0=\epsilon}, \eqn{q_1=0.5}, and
#' \eqn{q_2=1-\epsilon}. For every possible alternate-read count, the called
#' genotype maximizes its binomial likelihood among the three genotype models.
#'
#' Numerical likelihood ties within \code{1e-12} are resolved symmetrically:
#' below an alternate-read fraction of one half, the lower tied genotype is
#' selected; above one half, the higher tied genotype is selected; at exactly
#' one half, the tied genotype closest to 1 is selected.
#'
#' Rows correspond to true genotypes and columns to called genotypes. Thus, if
#' \code{E} is the returned matrix and \code{g_true} is a length-three vector of
#' true genotype probabilities, called-genotype probabilities are obtained as
#' \code{as.numeric(t(E) \%*\% g_true)}. Every row sums to one within numerical
#' precision.
#'
#' At finite coverage, even with zero sequencing error, a true heterozygote can
#' produce all-reference or all-alternate reads. Consequently the finite-depth
#' matrix need not be exactly the identity, especially at very low coverage.
#'
#' @return A numeric \code{3 x 3} matrix. Rows \code{true_0}, \code{true_1}, and
#'   \code{true_2} identify the true genotype; columns \code{called_0},
#'   \code{called_1}, and \code{called_2} identify the maximum-likelihood call.
#'
#' @examples
#' E <- ngs_genotype_error_matrix(coverage = 10, seq_error = 0.01)
#' g_true <- c(0.70, 0.25, 0.05)
#' as.numeric(t(E) %*% g_true)
#'
#' @export
ngs_genotype_error_matrix <- function(coverage, seq_error) {
  distributions <- lapply(0:2, function(genotype) {
    .ngs_read_count_distribution(
      coverage = coverage,
      genotype = genotype,
      seq_error = seq_error
    )
  })

  calls <- vapply(
    0:coverage,
    .ngs_call_genotype_ml,
    integer(1),
    coverage = coverage,
    seq_error = seq_error
  )

  E <- matrix(0, nrow = 3L, ncol = 3L)
  for (genotype in 0:2) {
    probabilities <- distributions[[genotype + 1L]]$probability
    for (called_genotype in 0:2) {
      E[genotype + 1L, called_genotype + 1L] <-
        sum(probabilities[calls == called_genotype])
    }
  }

  dimnames(E) <- list(
    c("true_0", "true_1", "true_2"),
    c("called_0", "called_1", "called_2")
  )

  boundary_tolerance <- 1e-12
  if (!identical(dim(E), c(3L, 3L))) {
    stop("Internal error: genotype transition matrix must be 3 x 3.")
  }
  if (any(!is.finite(E)) || any(E < -boundary_tolerance) ||
      any(E > 1 + boundary_tolerance)) {
    stop("Internal error: genotype transition probabilities must be finite and in [0, 1].")
  }
  # Snap only values that crossed a boundary through floating-point noise;
  # retain valid small tail probabilities and do not renormalize matrix rows.
  E[E < 0] <- 0
  E[E > 1] <- 1
  .validate_genotype_misclassification_matrix(E, tolerance = 1e-12)
}
