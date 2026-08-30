# Analytic read-probability helpers for the symmetric sequencing-error model.
#
# If G is the true alternate-allele count, V is coverage, X is the observed
# alternate-read count, and epsilon is the per-read sequencing-error rate, then
# X | G, V, epsilon follows Binomial(V, q_G), where q_0 = epsilon,
# q_1 = 0.5, and q_2 = 1 - epsilon.

.ngs_alt_read_probability <- function(genotype, seq_error) {
  if (!is.numeric(genotype) || length(genotype) != 1L ||
      !is.finite(genotype) || !(genotype %in% 0:2)) {
    stop("genotype must be a single finite value in {0, 1, 2}.")
  }
  if (!is.numeric(seq_error) || length(seq_error) != 1L ||
      !is.finite(seq_error) || seq_error < 0 || seq_error >= 0.5) {
    stop("seq_error must be a single finite number in [0, 0.5).")
  }

  ((2 - genotype) / 2) * seq_error +
    (genotype / 2) * (1 - seq_error)
}

.ngs_read_count_prob <- function(alt_reads, coverage, genotype, seq_error) {
  if (!is.numeric(coverage) || length(coverage) != 1L ||
      !is.finite(coverage) || coverage < 1 || coverage != floor(coverage)) {
    stop("coverage must be a single finite integer greater than or equal to 1.")
  }
  if (!is.numeric(alt_reads) || length(alt_reads) != 1L ||
      !is.finite(alt_reads) || alt_reads < 0 ||
      alt_reads > coverage || alt_reads != floor(alt_reads)) {
    stop("alt_reads must be a single finite integer between 0 and coverage.")
  }

  q <- .ngs_alt_read_probability(genotype, seq_error)
  probability <- stats::dbinom(alt_reads, size = coverage, prob = q)

  if (length(probability) != 1L || !is.finite(probability) ||
      probability < 0 || probability > 1) {
    stop("Computed read-count probability is not a finite value in [0, 1].")
  }

  probability
}

.ngs_read_count_distribution <- function(coverage, genotype, seq_error) {
  if (!is.numeric(coverage) || length(coverage) != 1L ||
      !is.finite(coverage) || coverage < 1 || coverage != floor(coverage)) {
    stop("coverage must be a single finite integer greater than or equal to 1.")
  }

  q <- .ngs_alt_read_probability(genotype, seq_error)
  alt_reads <- 0:coverage

  data.frame(
    alt_reads = alt_reads,
    probability = stats::dbinom(alt_reads, size = coverage, prob = q)
  )
}
