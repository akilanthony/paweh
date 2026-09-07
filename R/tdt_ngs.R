# Public prospective-power interface for the TDT1-NGS kernel implemented from
# Kim (2015), Appendix B. All statistical computation is delegated to the
# frozen raw-read likelihood kernel in tdt_ngs_ncp.R.

.tdt_ngs_resolve_sequencing_design <- function(
    coverage, seq_error,
    father_coverage, mother_coverage, child_coverage,
    epsilon0, epsilon1,
    father_supplied, mother_supplied, child_supplied,
    epsilon0_supplied, epsilon1_supplied
) {
  coverage_values <- list(
    father = if (father_supplied) father_coverage else coverage,
    mother = if (mother_supplied) mother_coverage else coverage,
    child = if (child_supplied) child_coverage else coverage
  )
  coverage_sources <- c(
    father = if (father_supplied) "father_coverage" else "coverage",
    mother = if (mother_supplied) "mother_coverage" else "coverage",
    child = if (child_supplied) "child_coverage" else "coverage"
  )
  for (member in names(coverage_values)) {
    value <- coverage_values[[member]]
    source <- coverage_sources[[member]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < 1) {
      stop(source, " must be a single finite integer greater than or equal to 2.")
    }
    if (value == 1) {
      stop(
        source, " = 1 is unsupported: TDT1-NGS efficient information is not ",
        "identifiable under the implemented 11-parameter nuisance model."
      )
    }
  }
  effective_coverage <- unlist(coverage_values, use.names = TRUE)

  if (!epsilon0_supplied || !epsilon1_supplied) {
    if (!is.numeric(seq_error) || length(seq_error) != 1L ||
        !is.finite(seq_error) || seq_error < 0 || seq_error >= 0.5) {
      stop("seq_error must be a single finite number in [0, 0.5).")
    }
  }
  effective_epsilon0 <- if (epsilon0_supplied) epsilon0 else seq_error
  effective_epsilon1 <- if (epsilon1_supplied) epsilon1 else seq_error
  if (!is.numeric(effective_epsilon0) || length(effective_epsilon0) != 1L ||
      !is.finite(effective_epsilon0) || effective_epsilon0 < 0 ||
      effective_epsilon0 >= 1) {
    stop("epsilon0 must be a single finite number in [0, 1).")
  }
  if (!is.numeric(effective_epsilon1) || length(effective_epsilon1) != 1L ||
      !is.finite(effective_epsilon1) || effective_epsilon1 < 0 ||
      effective_epsilon1 >= 1) {
    stop("epsilon1 must be a single finite number in [0, 1).")
  }
  if (effective_epsilon0 + effective_epsilon1 >= 1) {
    stop("epsilon0 + epsilon1 must be less than 1.")
  }

  list(
    coverage = effective_coverage,
    father_coverage = unname(effective_coverage[["father"]]),
    mother_coverage = unname(effective_coverage[["mother"]]),
    child_coverage = unname(effective_coverage[["child"]]),
    epsilon0 = effective_epsilon0,
    epsilon1 = effective_epsilon1,
    equal_coverage = length(unique(effective_coverage)) == 1L,
    symmetric_error = identical(effective_epsilon0, effective_epsilon1)
  )
}

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
#' @param coverage A single finite integer greater than or equal to 2. Common
#'   fixed-depth shorthand for the father, mother, and affected child. It may
#'   be omitted when all three member-specific depths are supplied.
#' @param seq_error A single finite symmetric per-read sequencing-error
#'   probability in \eqn{[0,0.5)}. Internally, the directional error parameters
#'   default to \eqn{\epsilon_0 = \epsilon_1 =} \code{seq_error}. It may be
#'   omitted when both directional parameters are supplied.
#' @param father_coverage,mother_coverage,child_coverage Optional fixed depths
#'   for each trio member. An explicitly supplied member depth overrides
#'   \code{coverage}; otherwise that member inherits \code{coverage}.
#' @param epsilon0 Directional probability that a reference-allele read is
#'   observed as the alternative allele. Defaults to \code{seq_error}.
#' @param epsilon1 Directional probability that an alternative-allele read is
#'   observed as the reference allele. Defaults to \code{seq_error}.
#' @param alpha A single finite significance level in \eqn{(0,1)}. Defaults to
#'   0.05.
#' @param verbose Logical scalar. If \code{TRUE}, print a concise result
#'   summary.
#'
#' @details
#' TDT1-NGS is evaluated for one biallelic variant under Hardy-Weinberg
#' parental genotype frequencies, random mating, a multiplicative disease
#' model, fixed member-specific coverage, and common directional sequencing
#' error. With
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
#' The existing directional read model is
#' \deqn{q_G=\epsilon_0 + (1-\epsilon_0-\epsilon_1)G/2,}
#' for genotype \eqn{G\in\{0,1,2\}}. Thus \eqn{q_0=\epsilon_0},
#' \eqn{q_2=1-\epsilon_1}, and
#' \eqn{q_1=(1+\epsilon_0-\epsilon_1)/2}; the heterozygote probability is
#' one-half in the symmetric case. The effective errors must be nonnegative
#' and sum to less than one.
#'
#' Coverage 1 is unsupported under the full published nuisance model. Its
#' eight observable read-count triples provide at most seven independent
#' probability dimensions for an 11-parameter information model, so efficient
#' information is not identifiable without changing the model or using a
#' generalized inverse.
#'
#' This prospective calculation performs no simulation or EM fitting. It does
#' not implement TDT2-NGS, locus heterogeneity, phenotype misclassification,
#' conventional genotype misclassification/TDTae, or multi-locus testing.
#' Coverage is fixed for each member rather than random or sample-specific.
#'
#' @return Invisibly, an object of class \code{"tdt_ngs_power"} containing the
#'   design inputs, power and NCP, multiplicative-model parameters, efficient
#'   information, the 11 by 11 information matrix, compact numerical
#'   diagnostics, effective member-specific coverage and directional-error
#'   metadata, and model metadata. Legacy \code{coverage} and \code{seq_error}
#'   fields retain the supplied common shorthand values.
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
    coverage = NULL,
    seq_error = NULL,
    alpha = 0.05,
    verbose = TRUE,
    father_coverage = coverage,
    mother_coverage = coverage,
    child_coverage = coverage,
    epsilon0 = seq_error,
    epsilon1 = seq_error
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
  sequencing <- .tdt_ngs_resolve_sequencing_design(
    coverage = coverage, seq_error = seq_error,
    father_coverage = father_coverage,
    mother_coverage = mother_coverage,
    child_coverage = child_coverage,
    epsilon0 = epsilon0, epsilon1 = epsilon1,
    father_supplied = !missing(father_coverage),
    mother_supplied = !missing(mother_coverage),
    child_supplied = !missing(child_coverage),
    epsilon0_supplied = !missing(epsilon0),
    epsilon1_supplied = !missing(epsilon1)
  )
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
    coverage = sequencing$coverage,
    seq_error = seq_error,
    epsilon0 = sequencing$epsilon0,
    epsilon1 = sequencing$epsilon1
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
    father_coverage = sequencing$father_coverage,
    mother_coverage = sequencing$mother_coverage,
    child_coverage = sequencing$child_coverage,
    epsilon0 = sequencing$epsilon0,
    epsilon1 = sequencing$epsilon1,
    symmetric_error = sequencing$symmetric_error,
    sequencing = sequencing,
    efficient_information = fit$efficient_information,
    information_matrix = fit$information_matrix,
    nuisance_rcond = fit$nuisance_rcond,
    score_mean = fit$score_mean,
    model_info = list(
      test = "TDT1-NGS",
      inheritance = "multiplicative",
      coverage_model = if (sequencing$equal_coverage) {
        "equal_fixed"
      } else {
        "member_specific_fixed"
      },
      sequencing_error = if (sequencing$symmetric_error) {
        "symmetric"
      } else {
        "directional"
      },
      trio_type = "father-mother-affected-child",
      likelihood = "raw_read_counts_with_latent_trio_states",
      information_evaluation = "null"
    )
  )
  class(out) <- "tdt_ngs_power"

  if (isTRUE(verbose)) {
    .paweh_print_tdt_ngs_power(out)
  }
  invisible(out)
}

#' @export
print.tdt_ngs_power <- function(x, ...) {
  cat("TDT1-NGS analytic power\n")
  cat(sprintf("Affected-child trios: %s\n",
              formatC(x$N, format = "d", big.mark = ",")))
  cat(sprintf("Disease allele frequency: %.4g; R1: %.4g\n", x$pd, x$R1))
  if (is.null(x$sequencing) ||
      (isTRUE(x$sequencing$equal_coverage) &&
       isTRUE(x$sequencing$symmetric_error) &&
       !is.null(x$coverage) && !is.null(x$seq_error))) {
    cat(sprintf("Coverage: %s; sequencing error: %.4g\n",
                formatC(x$coverage, format = "d"), x$seq_error))
  } else {
    cat(sprintf("Coverage (father/mother/child): %s/%s/%s\n",
                x$father_coverage, x$mother_coverage, x$child_coverage))
    if (isTRUE(x$symmetric_error)) {
      cat(sprintf("Symmetric sequencing error: %.4g\n", x$epsilon0))
    } else {
      cat(sprintf(
        "Directional sequencing error: epsilon0 = %.4g; epsilon1 = %.4g\n",
        x$epsilon0, x$epsilon1
      ))
    }
  }
  cat(sprintf("Alpha: %.4g; NCP: %.4f; power: %.1f%%\n",
              x$alpha, x$lambda, 100 * x$power))
  invisible(x)
}

.tdt_ngs_target_ncp <- function(power, alpha) {
  if (power <= alpha) {
    return(0)
  }
  cc_chisq_ncp_target(power = power, alpha = alpha, df = 1)
}

#' Analytic MSSN for a TDT1-NGS Sequencing Study
#'
#' Computes the prospective minimum sample size necessary (MSSN), measured in
#' complete father-mother-affected-child trios, for the same single-variant
#' TDT1-NGS design as \code{tdt_ngs_power()}.
#'
#' @param power A single finite target power strictly between 0 and 1.
#' @param pd A single finite disease/risk-allele frequency strictly between 0
#'   and 1. The supplied allele labeling is retained.
#' @param R1 A single finite positive heterozygote genotype relative risk under
#'   the multiplicative model. The homozygote relative risk is
#'   \eqn{R_2 = R_1^2}.
#' @param coverage A single finite integer greater than or equal to 2. Common
#'   fixed-depth shorthand for the father, mother, and affected child. It may
#'   be omitted when all three member-specific depths are supplied.
#' @param seq_error A single finite symmetric per-read sequencing-error
#'   probability in \eqn{[0,0.5)}. Internally, the directional error parameters
#'   default to \eqn{\epsilon_0 = \epsilon_1 =} \code{seq_error}. It may be
#'   omitted when both directional parameters are supplied.
#' @param father_coverage,mother_coverage,child_coverage Optional fixed depths
#'   for each trio member. An explicitly supplied member depth overrides
#'   \code{coverage}; otherwise that member inherits \code{coverage}.
#' @param epsilon0 Directional probability that a reference-allele read is
#'   observed as the alternative allele. Defaults to \code{seq_error}.
#' @param epsilon1 Directional probability that an alternative-allele read is
#'   observed as the reference allele. Defaults to \code{seq_error}.
#' @param alpha A single finite significance level in \eqn{(0,1)}. Defaults to
#'   0.05.
#' @param verbose Logical scalar. If \code{TRUE}, print a concise result
#'   summary.
#'
#' @details
#' The calculation uses Kim's Appendix B TDT1-NGS noncentrality parameter
#' \deqn{\lambda = N \delta^2 I_{eff},}
#' where \eqn{N} is the number of complete trios, \eqn{I_{eff}} is the
#' nuisance-adjusted per-trio information under the null, and under the
#' multiplicative model \eqn{\delta = \log(R_1)}. The target one-degree-of-
#' freedom NCP, \eqn{\lambda_*}, is obtained by numerical inversion of the
#' noncentral chi-square power function. Trio MSSN is then solved analytically:
#' \deqn{N_{continuous} = \frac{\lambda_*}{\log(R_1)^2 I_{eff}}.}
#' The planned MSSN is the ceiling of this quantity, with achieved power and
#' the immediately smaller design checked at the integer boundary.
#'
#' If target \code{power} is no greater than \code{alpha}, the target NCP is
#' zero and the minimum supported design is one trio. If \code{R1 = 1} and
#' target power exceeds alpha, no finite MSSN exists because the transmission
#' effect is zero.
#'
#' Each member uses its resolved fixed coverage, and the directional read model
#' is \eqn{q_G=\epsilon_0+(1-\epsilon_0-\epsilon_1)G/2}. See
#' \code{\link{tdt_ngs_power}} for parameter precedence and interpretation.
#' Coverage 1 is unsupported for every member because efficient information is
#' not identifiable under the implemented 11-parameter nuisance model. This prospective
#' calculation uses raw sequencing read-count information and latent trio
#' genotype states; it performs no simulation, genotype calling, or EM fitting.
#'
#' @return Invisibly, an object of class \code{"tdt_ngs_mssn"} containing the
#'   target, continuous and integer trio requirements, achieved power and NCP,
#'   per-trio NCP coefficient, multiplicative-model parameters, efficient
#'   information, the 11 by 11 information matrix, numerical diagnostics,
#'   effective sequencing-design metadata, and model metadata. Legacy
#'   \code{coverage} and \code{seq_error} fields retain supplied shorthand.
#'
#' @references
#' Kim, W. (2015). Transmission disequilibrium tests based on read counts for
#' low-coverage next-generation sequence data. \emph{Human Heredity}, 80(1),
#' 36--49. \doi{10.1159/000434645}.
#'
#' @seealso \code{\link{tdt_ngs_power}}
#'
#' @examples
#' tdt_ngs_mssn(
#'   power = 0.80, pd = 0.325, R1 = 1.2,
#'   coverage = 12, seq_error = 0.005,
#'   alpha = 5e-8, verbose = FALSE
#' )
#'
#' @importFrom stats pchisq qchisq
#' @export
tdt_ngs_mssn <- function(
    power,
    pd,
    R1,
    coverage = NULL,
    seq_error = NULL,
    alpha = 0.05,
    verbose = TRUE,
    father_coverage = coverage,
    mother_coverage = coverage,
    child_coverage = coverage,
    epsilon0 = seq_error,
    epsilon1 = seq_error
) {
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0 || power >= 1) {
    stop("power must be a single finite number in (0, 1).")
  }
  if (!is.numeric(pd) || length(pd) != 1L || !is.finite(pd) ||
      pd <= 0 || pd >= 1) {
    stop("pd must be a single finite number strictly between 0 and 1.")
  }
  if (!is.numeric(R1) || length(R1) != 1L || !is.finite(R1) || R1 <= 0) {
    stop("R1 must be a single finite positive number.")
  }
  sequencing <- .tdt_ngs_resolve_sequencing_design(
    coverage = coverage, seq_error = seq_error,
    father_coverage = father_coverage,
    mother_coverage = mother_coverage,
    child_coverage = child_coverage,
    epsilon0 = epsilon0, epsilon1 = epsilon1,
    father_supplied = !missing(father_coverage),
    mother_supplied = !missing(mother_coverage),
    child_supplied = !missing(child_coverage),
    epsilon0_supplied = !missing(epsilon0),
    epsilon1_supplied = !missing(epsilon1)
  )
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number in (0, 1).")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  information <- .tdt_ngs_information(
    pd = pd,
    coverage = sequencing$coverage,
    seq_error = seq_error,
    epsilon0 = sequencing$epsilon0,
    epsilon1 = sequencing$epsilon1
  )
  efficient_information <- information$efficient_information
  if (!is.finite(efficient_information) || efficient_information <= 0) {
    stop(
      "No finite TDT1-NGS MSSN can be computed because efficient ",
      "information is not positive for this design."
    )
  }

  t <- R1 / (1 + R1)
  delta <- log(t / (1 - t))
  coefficient <- delta^2 * efficient_information
  lambda_target <- .tdt_ngs_target_ncp(power, alpha)

  if (coefficient == 0 && lambda_target > 0) {
    stop(
      "No finite MSSN exists because R1 = 1 implies zero transmission effect."
    )
  }
  N_continuous <- if (lambda_target == 0) {
    0
  } else {
    lambda_target / coefficient
  }
  if (!is.finite(N_continuous) || N_continuous < 0) {
    stop("The analytic TDT1-NGS continuous MSSN must be finite and nonnegative.")
  }

  critical <- stats::qchisq(1 - alpha, df = 1)
  evaluate_integer <- function(N) {
    lambda <- N * coefficient
    achieved <- stats::pchisq(
      critical,
      df = 1,
      ncp = lambda,
      lower.tail = FALSE
    )
    list(N = N, lambda = as.numeric(lambda), power = as.numeric(achieved))
  }

  initial_MSSN <- max(1, ceiling(N_continuous))
  planned <- evaluate_integer(initial_MSSN)

  # Only an adjacent floating-point boundary correction is possible because
  # lambda is exactly linear in the single integer sampling dimension N.
  if (planned$power < power) {
    planned <- evaluate_integer(planned$N + 1)
  }
  while (planned$N > 1) {
    previous <- evaluate_integer(planned$N - 1)
    if (previous$power < power) {
      break
    }
    planned <- previous
  }

  out <- list(
    power_target = power,
    alpha = alpha,
    MSSN_trios = planned$N,
    total_individuals = 3 * planned$N,
    N_trios_continuous = N_continuous,
    achieved_power = planned$power,
    achieved_lambda = planned$lambda,
    lambda_target = lambda_target,
    ncp_per_trio = coefficient,
    initial_MSSN_trios = initial_MSSN,
    rounding_adjustment = planned$N - initial_MSSN,
    pd = pd,
    R1 = R1,
    R2 = R1^2,
    t = t,
    delta = delta,
    coverage = coverage,
    seq_error = seq_error,
    father_coverage = sequencing$father_coverage,
    mother_coverage = sequencing$mother_coverage,
    child_coverage = sequencing$child_coverage,
    epsilon0 = sequencing$epsilon0,
    epsilon1 = sequencing$epsilon1,
    symmetric_error = sequencing$symmetric_error,
    sequencing = sequencing,
    efficient_information = efficient_information,
    information_matrix = information$information_matrix,
    nuisance_rcond = information$nuisance_rcond,
    score_mean = information$score_mean,
    model_info = list(
      test = "TDT1-NGS",
      objective = "MSSN",
      sampling_unit = "complete_trios",
      inheritance = "multiplicative",
      coverage_model = if (sequencing$equal_coverage) {
        "equal_fixed"
      } else {
        "member_specific_fixed"
      },
      sequencing_error = if (sequencing$symmetric_error) {
        "symmetric"
      } else {
        "directional"
      },
      trio_type = "father-mother-affected-child",
      likelihood = "raw_read_counts_with_latent_trio_states",
      information_evaluation = "null",
      sample_size_solution = "analytic"
    )
  )
  class(out) <- "tdt_ngs_mssn"

  if (isTRUE(verbose)) {
    .paweh_print_tdt_ngs_mssn(out)
  }
  invisible(out)
}

#' @export
print.tdt_ngs_mssn <- function(x, ...) {
  cat("TDT1-NGS analytic MSSN\n")
  cat(sprintf("Target power: %.1f%%; alpha: %.4g\n",
              100 * x$power_target, x$alpha))
  cat(sprintf("Required trios: %s; total individuals: %s\n",
              formatC(x$MSSN_trios, format = "d", big.mark = ","),
              formatC(x$total_individuals, format = "d", big.mark = ",")))
  cat(sprintf("Disease allele frequency: %.4g; R1: %.4g\n", x$pd, x$R1))
  if (is.null(x$sequencing) ||
      (isTRUE(x$sequencing$equal_coverage) &&
       isTRUE(x$sequencing$symmetric_error) &&
       !is.null(x$coverage) && !is.null(x$seq_error))) {
    cat(sprintf("Coverage: %s; sequencing error: %.4g\n",
                formatC(x$coverage, format = "d"), x$seq_error))
  } else {
    cat(sprintf("Coverage (father/mother/child): %s/%s/%s\n",
                x$father_coverage, x$mother_coverage, x$child_coverage))
    if (isTRUE(x$symmetric_error)) {
      cat(sprintf("Symmetric sequencing error: %.4g\n", x$epsilon0))
    } else {
      cat(sprintf(
        "Directional sequencing error: epsilon0 = %.4g; epsilon1 = %.4g\n",
        x$epsilon0, x$epsilon1
      ))
    }
  }
  cat(sprintf("Achieved power: %.1f%%\n", 100 * x$achieved_power))
  invisible(x)
}
