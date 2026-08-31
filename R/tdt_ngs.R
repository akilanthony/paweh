# Public prospective-power interface for the TDT1-NGS kernel implemented from
# Kim (2015), Appendix B. All statistical computation is delegated to the
# frozen raw-read likelihood kernel in tdt_ngs_ncp.R.

#' Analytic Power for a TDT1-NGS Sequencing Study
#'
#' Computes prospective analytic power for a single-variant TDT1-NGS study of
#' complete father-mother-affected-child trios. The implementation uses the
#' published latent-state, sequencing-read-count likelihood rather than hard
#' genotype calls.
#'
#' @param N A single finite integer greater than or equal to 1. Number of
#'   complete father-mother-affected-child trios.
#' @param pd A single finite disease/risk-allele frequency strictly between 0
#'   and 1. The null information is reflection-symmetric in \code{pd} and
#'   \code{1 - pd}; the supplied allele labeling is retained.
#' @param R1 A single finite positive heterozygote genotype relative risk under
#'   the multiplicative model. The homozygote relative risk is
#'   \eqn{R_2 = R_1^2}.
#' @param coverage A single finite integer greater than or equal to 2. This is
#'   the equal fixed read depth for the father, mother, and affected child.
#' @param seq_error A single finite symmetric per-read sequencing-error
#'   probability in \eqn{[0,0.5)}. Internally, the directional error parameters
#'   are evaluated at \eqn{\epsilon_0 = \epsilon_1 =} \code{seq_error}.
#' @param alpha A single finite significance level in \eqn{(0,1)}. Defaults to
#'   0.05.
#' @param verbose Logical scalar. If \code{TRUE}, print a concise result
#'   summary.
#'
#' @details
#' TDT1-NGS is evaluated for one biallelic variant under Hardy-Weinberg
#' parental genotype frequencies, random mating, a multiplicative disease
#' model, equal fixed coverage, and symmetric public sequencing error. With
#' \eqn{t = R_1/(1+R_1)}, the transmission parameter is
#' \eqn{\delta = \log\{t/(1-t)\} = \log(R_1)} and \eqn{R_2 = R_1^2}.
#'
#' Kim's Appendix B noncentrality parameter is
#' \deqn{\lambda = N \delta^2 I_{eff},}
#' where \eqn{I_{eff}} is the nuisance-adjusted per-trio information evaluated
#' under the null from raw sequencing read-count probabilities and the 15
#' latent Mendelian trio states. Power is the upper-tail probability beyond
#' the central one-degree-of-freedom chi-square critical value under a
#' noncentral chi-square distribution with NCP \eqn{\lambda}.
#'
#' Coverage 1 is unsupported under the full published nuisance model. Its
#' eight observable read-count triples provide at most seven independent
#' probability dimensions for an 11-parameter information model, so efficient
#' information is not identifiable without changing the model or using a
#' generalized inverse.
#'
#' This prospective calculation performs no simulation or EM fitting. It does
#' not implement TDT2-NGS, unequal member-specific coverage, locus
#' heterogeneity, phenotype misclassification, or multi-locus testing.
#'
#' @return Invisibly, an object of class \code{"tdt_ngs_power"} containing the
#'   design inputs, power and NCP, multiplicative-model parameters, efficient
#'   information, the 11 by 11 information matrix, compact numerical
#'   diagnostics, and model metadata.
#'
#' @references
#' Kim, W. (2015). Transmission disequilibrium tests based on read counts for
#' low-coverage next-generation sequence data. \emph{Human Heredity}, 80(1),
#' 36--49. \doi{10.1159/000434645}.
#'
#' @examples
#' tdt_ngs_power(
#'   N = 5000, pd = 0.325, R1 = 1.2,
#'   coverage = 12, seq_error = 0.005,
#'   alpha = 5e-8, verbose = FALSE
#' )
#'
#' @importFrom stats pchisq qchisq
#' @export
tdt_ngs_power <- function(
    N,
    pd,
    R1,
    coverage,
    seq_error,
    alpha = 0.05,
    verbose = TRUE
) {
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
      N < 1 || N != floor(N)) {
    stop("N must be a single finite integer greater than or equal to 1.")
  }
  if (!is.numeric(pd) || length(pd) != 1L || !is.finite(pd) ||
      pd <= 0 || pd >= 1) {
    stop("pd must be a single finite number strictly between 0 and 1.")
  }
  if (!is.numeric(R1) || length(R1) != 1L || !is.finite(R1) || R1 <= 0) {
    stop("R1 must be a single finite positive number.")
  }
  if (!is.numeric(coverage) || length(coverage) != 1L ||
      !is.finite(coverage) || coverage != floor(coverage) || coverage < 1) {
    stop("coverage must be a single finite integer greater than or equal to 2.")
  }
  if (coverage == 1) {
    stop(
      "coverage = 1 is unsupported: TDT1-NGS efficient information is not ",
      "identifiable under the implemented 11-parameter nuisance model."
    )
  }
  if (!is.numeric(seq_error) || length(seq_error) != 1L ||
      !is.finite(seq_error) || seq_error < 0 || seq_error >= 0.5) {
    stop("seq_error must be a single finite number in [0, 0.5).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number in (0, 1).")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  fit <- .tdt_ngs_ncp(
    N = N,
    pd = pd,
    R1 = R1,
    coverage = coverage,
    seq_error = seq_error
  )
  critical <- stats::qchisq(1 - alpha, df = 1)
  power <- stats::pchisq(
    critical,
    df = 1,
    ncp = fit$lambda,
    lower.tail = FALSE
  )

  if (!is.finite(power) || power < 0 || power > 1) {
    stop("Computed TDT1-NGS power must be a finite probability in [0, 1].")
  }

  out <- list(
    N = N,
    alpha = alpha,
    power = as.numeric(power),
    lambda = fit$lambda,
    pd = pd,
    R1 = fit$R1,
    R2 = fit$R2,
    t = fit$t,
    delta = fit$delta,
    coverage = coverage,
    seq_error = seq_error,
    efficient_information = fit$efficient_information,
    information_matrix = fit$information_matrix,
    nuisance_rcond = fit$nuisance_rcond,
    score_mean = fit$score_mean,
    model_info = list(
      test = "TDT1-NGS",
      inheritance = "multiplicative",
      coverage_model = "equal_fixed",
      sequencing_error = "symmetric",
      trio_type = "father-mother-affected-child",
      likelihood = "raw_read_counts_with_latent_trio_states",
      information_evaluation = "null"
    )
  )
  class(out) <- "tdt_ngs_power"

  if (isTRUE(verbose)) {
    print(out)
  }
  invisible(out)
}

#' @export
print.tdt_ngs_power <- function(x, ...) {
  cat("TDT1-NGS analytic power\n")
  cat(sprintf("Affected-child trios: %s\n",
              formatC(x$N, format = "d", big.mark = ",")))
  cat(sprintf("Disease allele frequency: %.4g; R1: %.4g\n", x$pd, x$R1))
  cat(sprintf("Coverage: %s; sequencing error: %.4g\n",
              formatC(x$coverage, format = "d"), x$seq_error))
  cat(sprintf("Alpha: %.4g; NCP: %.4f; power: %.1f%%\n",
              x$alpha, x$lambda, 100 * x$power))
  invisible(x)
}
