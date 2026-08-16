#' Transmission Disequilibrium Test (TDT) Power from Expected Transmissions and Non-Transmissions
#'
#' Computes the statistical power of the Transmission Disequilibrium Test (TDT)
#' given the expected number of transmissions (ET) and non-transmissions (ENT)
#' under a specified significance level. Implements Eq. 1.25 (Chapter 1,
#' Section 1.6.1.3, p. 27) of Gordon, Finch, and Kim (2020).

#' @param ET Numeric. Expected number of transmissions.
#' @param ENT Numeric. Expected number of non-transmissions.
#' @param alpha Numeric. Significance level (default = 0.05).
#'
#'
#' @details
#' The function calculates the non-centrality parameter and statistical power for the TDT
#' using the chi-square distribution with 1 degree of freedom.
#' \code{ET} and \code{ENT} are expected counts, not probabilities, and must
#' refer to the same affected-trio design.
#'
#' The non-centrality parameter is computed as:
#' \deqn{\lambda = \frac{(ET - ENT)^2}{ET + ENT}}
#'
#' Power is then obtained as:
#' \deqn{1 - P(\chi^2_{1,\lambda} < q_{\chi^2_{1,1-\alpha}})}
#'
#' @return A list containing:
#' \item{lambda}{Non-centrality parameter.}
#' \item{power}{Computed power at the given alpha level.}
#' \item{ET}{Expected transmissions.}
#' \item{ENT}{Expected non-transmissions.}
#'
#'
#' @importFrom stats pchisq qchisq uniroot
#'
#' @examples
#' # Example: compute power for ET = 140 and ENT = 100
#' tdt_power_from_expected_counts(ET = 140, ENT = 100, alpha = 0.05)
#'
#' @references
#' Spielman, R. S., McGinnis, R. E., & Ewens, W. J. (1993). Transmission test
#' for linkage disequilibrium: the insulin gene region and insulin-dependent
#' diabetes mellitus. \emph{American Journal of Human Genetics}, 52(3),
#' 506--516. PMID: 8447318; PMCID: PMC1682161.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' @seealso \code{\link{tdt_power}} and \code{\link{tdt_power_from_model}}.
#'
#'@export

tdt_power_from_expected_counts <- function(ET, ENT, alpha = 0.05) {
  lambda <- (ET - ENT)^2 / (ET + ENT)
  thr <- qchisq(1 - alpha, df = 1)
  power <- 1 - pchisq(thr, df = 1, ncp = lambda)

  message("\n--- Transmission Disequilibrium Test (TDT) ---")
  message("Equation: 1.25  |  Input: Expected Transmissions and Non-Transmissions")
  message("-----------------------------------------------------------")

  # Use sprintf() for alignment (value column starts at same width)
  message(sprintf("%-38s %10.4f", "Expected Transmissions (ET):", ET))
  message(sprintf("%-38s %10.4f", "Expected Non-Transmissions (ENT):", ENT))
  message(sprintf("%-38s %10.4f", "Non-Centrality Parameter (lambda):", lambda))
  message(sprintf("%-38s %10.4f", paste0("Power at alpha = ", alpha, ":"), power))

  message("-----------------------------------------------------------")

  invisible(list(
    `Expected Transmissions (ET)` = ET,
    `Expected Non-Transmissions (ENT)` = ENT,
    `Non-Centrality Parameter (lambda)` = lambda,
    `Power` = power
  ))
}


#' Transmission Disequilibrium Test (TDT) Power from Genetic Model Parameters
#'
#' Computes the statistical power of the Transmission Disequilibrium Test (TDT)
#' using genetic model parameters such as allele frequency, relative risks,
#' disease prevalence, and the number of affected trios. Implements the
#' penetrance construction in Eqs. 1.6--1.7 (p. 13) and TDT Eq. 1.25
#' (p. 27) of Gordon, Finch, and Kim (2020).
#'
#' @param pd Numeric. Frequency of the disease-associated allele.
#' @param N Numeric. Number of affected trios.
#' @param delta_prime Numeric. Linkage disequilibrium (LD) scale factor (default = 1).
#' @param f0,f1,f2 Optional. Penetrances for genotypes with 0, 1, and 2 risk alleles.
#'   If not provided, they are computed internally using \code{prev}, \code{R1}, and \code{R2}.
#' @param prev Numeric. Disease prevalence.
#' @param R1 Numeric. Relative risk for heterozygotes.
#' @param R2 Numeric. Relative risk for homozygotes.
#' @param alpha Numeric. Significance level (default = 0.05).
#'
#' @details
#' When penetrances (\code{f0}, \code{f1}, \code{f2}) are not provided, they are
#' derived from the model parameters using:
#' \deqn{f0 = prev / Z, \quad f1 = R1 * f0, \quad f2 = R2 * f0}
#' where
#' \deqn{Z = (1 - p_d)^2 + 2 * p_d * (1 - p_d) * R1 + p_d^2 * R2.}
#'
#' Expected transmission (\eqn{ET}) and non-transmission (\eqn{ENT}) counts are
#' computed based on the allele frequency and penetrance model, and the power
#' is derived from the non-central chi-square distribution with 1 degree of freedom.
#' Here \code{N} is the number of affected-child trios. \code{delta_prime}
#' scales \eqn{D=p_d(1-p_d)D'} in the implemented model; the disease-locus and
#' marker-locus assumptions should be considered when interpreting \code{pd}.
#'
#' @return A list containing:
#' \item{lambda}{Non-centrality parameter.}
#' \item{power}{Computed power at the given significance level.}
#' \item{ET}{Expected transmissions.}
#' \item{ENT}{Expected non-transmissions.}
#' \item{Penetrances}{Vector of computed penetrances (f0, f1, f2).}
#'
#' @examples
#' # Example: model-based power computation
#' tdt_power_from_model(
#'   pd = 0.25, N = 10000, delta_prime = 1,
#'   prev = 0.05, R1 = 1, R2 = 1.1, alpha = 0.05
#' )
#'
#' @references
#' Spielman, R. S., McGinnis, R. E., & Ewens, W. J. (1993). Transmission test
#' for linkage disequilibrium: the insulin gene region and insulin-dependent
#' diabetes mellitus. \emph{American Journal of Human Genetics}, 52(3),
#' 506--516. PMID: 8447318; PMCID: PMC1682161.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' @seealso \code{\link{tdt_power}},
#' \code{\link{tdt_power_from_expected_counts}}, and
#' \code{\link{tdt_mssn_from_model}}.
#'
#' @importFrom stats pchisq qchisq
#' @export

tdt_power_from_model <- function(pd, N, delta_prime,
                                 f0 = NULL, f1 = NULL, f2 = NULL,
                                 prev = NULL, R1 = NULL, R2 = NULL,
                                 alpha = 0.05) {
  p_plus <- 1 - pd

  if (is.null(f0) || is.null(f1) || is.null(f2)) {
    if (is.null(prev) || is.null(R1) || is.null(R2))
      stop("Provide either f0,f1,f2 OR prev,R1,R2.")
    Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
    f0 <- prev / Z
    f1 <- R1 * f0
    f2 <- R2 * f0
  }

  C     <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  prev_ <- pd^2 * f2 + 2 * pd * p_plus * f1 + p_plus^2 * f0
  delta <- delta_prime * pd * p_plus

  ET  <- 2 * N * (pd * p_plus + delta * p_plus * C / prev_)
  ENT <- 2 * N * (pd * p_plus - delta * pd    * C / prev_)

  lambda <- (ET - ENT)^2 / (ET + ENT)
  power  <- 1 - pchisq(qchisq(1 - alpha, df = 1), df = 1, ncp = lambda)

  message("\n--- Transmission Disequilibrium Test (Model-Based) ---")
  message("Equation: 1.25  |  Inputs: Allele Frequency and Genetic Model Parameters")
  message("-----------------------------------------------------------")

  message(sprintf("%-38s %10.4f", "Allele Frequency (p_d):", pd))
  message(sprintf("%-38s %10.4f", "Prevalence (phi1):", prev))
  message(sprintf("%-38s %10s",   "Relative Risks (R1,R2):",
                  paste0(R1, ", ", R2)))
  message(sprintf("%-38s %10.0f", "Number of Affected Trios (N):", N))

  message("-----------------------------------------------------------")
  message(sprintf("%-38s %10.4f", "Non-Centrality Parameter (lambda):", lambda))
  message(sprintf("%-38s %10.4f", paste0("Power at alpha = ", alpha, ":"), power))
  message(sprintf("%-38s %10.2f  |  Expected ENT: %10.2f", "Expected ET:", ET, ENT))
  message("-----------------------------------------------------------")

  invisible(list(
    `Non-Centrality Parameter (lambda)` = lambda,
    `Power` = power,
    `Expected Transmissions (ET)` = ET,
    `Expected Non-Transmissions (ENT)` = ENT,
    `Penetrances` = c(f0 = f0, f1 = f1, f2 = f2),
    `C` = C,
    `Delta` = delta
  ))
}


#' TDT Minimum Sample Size Necessary from Genetic Model Parameters
#'
#' Computes the required number of affected trios (\eqn{N^*}) needed to achieve
#' a specified statistical power in the Transmission Disequilibrium Test (TDT),
#' given model parameters for allele frequency, relative risks, disease prevalence,
#' and heterogeneity. Implements the probability and MSSN chain in
#' Eqs. 5.33--5.34b (Chapter 5, Section 5.3.3, pp. 293--294) of Gordon,
#' Finch, and Kim (2020).
#'
#' @param power Numeric. Desired power (e.g., 0.8).
#' @param alpha Numeric. Significance level (e.g., 0.05).
#' @param df Integer. Degrees of freedom (typically 1 for TDT).
#' @param pd Numeric. Frequency of the disease-associated allele.
#' @param prev Numeric. Disease prevalence.
#' @param R1 Numeric. Relative risk for heterozygotes.
#' @param R2 Numeric. Relative risk for homozygotes.
#' @param delta_prime Numeric. Linkage disequilibrium (LD) scale factor (default = 1).
#' @param pi Numeric in \eqn{[0,1]}. Linked/homogeneous trio fraction; \code{1}
#'   is complete homogeneity and \eqn{1-\pi} is the heterogeneous fraction.
#'
#' @details
#' This function determines the non-centrality parameter (\eqn{\lambda^*}) via
#' root-finding (\code{uniroot}) such that the test power equals the desired level.
#' It then computes the expected transmission (\eqn{gT^*}) and non-transmission
#' (\eqn{gNT^*}) probabilities, followed by the required number of trios using:
#' \deqn{N^* = \frac{\lambda^*}{2} \frac{(gT^* + gNT^*)}{(gT^* - gNT^*)^2}}
#'
#' The expected transmission and non-transmission components are calculated
#' under allele frequency and penetrance model assumptions in Eq. 5.33. Their
#' NCP and MSSN are Eqs. 5.34a--b.
#'
#' @return A list containing:
#' \item{lambda_star}{Non-centrality parameter (\eqn{\lambda^*}).}
#' \item{gT_star}{Expected transmission probability.}
#' \item{gNT_star}{Expected non-transmission probability.}
#' \item{N_star}{Required number of trios.}
#'
#' @examples
#' # Example: compute required trios for 80% power at alpha = 0.05
#' tdt_mssn_from_model(
#'   power = 0.8, alpha = 0.05, df = 1,
#'   pd = 0.25, prev = 0.005, R1 = 2, R2 = 2,
#'   delta_prime = 1, pi = 1
#' )
#'
#' @references
#' Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
#' (2009). Transmission disequilibrium test power and sample size in the
#' presence of locus heterogeneity. \emph{Statistical Applications in Genetics
#' and Molecular Biology}, 8, Article 44. \doi{10.2202/1544-6115.1501}.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' @seealso \code{\link{tdt_mssn}}, \code{\link{tdt_power_from_model}}, and
#' \code{\link{tdt_expected_transmission_counts}}.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_mssn_from_model <- function(power, alpha, df,
                               pd, prev, R1, R2,
                               delta_prime = 1,
                               pi = 1) {
  crit <- qchisq(1 - alpha, df)
  f <- function(ncp) pchisq(crit, df, ncp = ncp, lower.tail = FALSE) - power
  ncp_star <- uniroot(f, c(0, 1e4))$root

  p_plus <- 1 - pd
  Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
  f0 <- prev / Z; f1 <- R1 * f0; f2 <- R2 * f0
  C <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  delta <- delta_prime * pd * p_plus

  gT_star  <- pi * ( p_plus * ( delta * C / prev + pd ) ) +
    (1 - pi) * ( delta * (p_plus - 0.5) * C / prev + pd * p_plus )
  gNT_star <- pi * ( pd * ( p_plus - delta * C / prev ) ) +
    (1 - pi) * ( delta * (0.5 - pd) * C / prev + pd * p_plus )

  N_star <- (ncp_star / 2) * ((gT_star + gNT_star) / (gT_star - gNT_star)^2)

  message("\n--- Transmission Disequilibrium Test (Trios) ---")
  message("Equation: 5.34b  |  Computes Required Number of Trios (N_star)")
  message("-----------------------------------------------------------")
  message(sprintf("%-38s %10.3f  |  %s %6.3f",
                  "Desired Power:", power, "Significance Level (alpha):", alpha))
  message(sprintf("%-38s %10.3f  |  %s %6.3f",
                  "Allele Frequency (p_d):", pd, "Prevalence (phi1):", prev))
  message(sprintf("%-38s %10s", "Relative Risks (R1,R2):", paste0(R1, ", ", R2)))
  message(sprintf("%-38s %10.3f", "Heterogeneity Parameter (pi):", pi))
  message("-----------------------------------------------------------")
  message(sprintf("%-38s %10.4f", "Non-Centrality Parameter (lambda_star):", ncp_star))
  message(sprintf("%-38s %10.5f", "Expected Transmission (gT_star):", gT_star))
  message(sprintf("%-38s %10.5f", "Expected Non-Transmission (gNT_star):", gNT_star))
  message(sprintf("%-38s %10.0f", "Required Number of Trios (N_star):", ceiling(N_star)))
  message("-----------------------------------------------------------")

  invisible(list(
    `Non-Centrality Parameter (lambda_star)` = ncp_star,
    `Expected Transmission (gT_star)` = gT_star,
    `Expected Non-Transmission (gNT_star)` = gNT_star,
    `Required Number of Trios (N_star)` = N_star
  ))
}


#' TDT Minimum Sample Size Necessary from an NCP under Phenotype Misclassification
#'
#' Computes the expected transmissions (\eqn{gT^*}) and non-transmissions (\eqn{gNT^*})
#' as well as the required number of trios (\eqn{N^*}) for a specified non-centrality parameter
#' (\eqn{\lambda^*}) under a misclassification model. Implements Eq. 5.26 and
#' Eqs. 5.27a--b (Chapter 5, Section 5.2.6, pp. 284--285) of Gordon, Finch,
#' and Kim (2020). Eq. 5.28a--b is the book's numerical worked example.
#'
#' @param lambda_star Numeric. Non-centrality parameter (\eqn{\lambda^*}) derived from
#' desired power (e.g., from \code{tdt_mssn_from_model()}).
#' @param pd Numeric. Frequency of the disease-associated allele.
#' @param prev Numeric. Disease prevalence (\eqn{\phi_1}).
#' @param R1 Numeric. Relative risk for heterozygotes.
#' @param R2 Numeric. Relative risk for homozygotes.
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter (\eqn{D'}) (default = 1).
#' @param pi01 Numeric. Misclassification rate for controls (default = 0).
#' @param digits Integer. Number of digits to round intermediate printed results (default = 5).
#'
#' @details
#' This function adjusts the Transmission Disequilibrium Test (TDT) for possible
#' phenotype misclassification. When \eqn{\pi_{01} > 0}, a fraction of unaffected
#' individuals are incorrectly classified as affected, inflating the apparent
#' transmission probability. The expected transmission (\eqn{gT^*}) and non-transmission
#' (\eqn{gNT^*}) are calculated as:
#'
#' \deqn{gT^* = p_d p_+ + \frac{D(p_T - \theta_1)C(1 - \pi_{01})}{\phi_1 + \pi_{01}\phi_0}}
#' \deqn{gNT^* = p_d p_+ + \frac{D(p_A - \theta_1)C(\pi_{01} - 1)}{\phi_1 + \pi_{01}\phi_0}}
#'
#' The required number of trios is then:
#' \deqn{N^* = \frac{\lambda^* (gT^* + gNT^*)}{2 (gT^* - gNT^*)^2}}
#'
#' Setting \eqn{\pi_{01} = 0} reproduces the standard (non-misclassified) TDT.
#'
#' @return A list containing:
#' \item{gT_star}{Expected transmission probability.}
#' \item{gNT_star}{Expected non-transmission probability.}
#' \item{N_required}{Required number of trios (\eqn{N^*}).}
#' \item{lambda_star}{Non-centrality parameter.}
#' \item{C, f0, f1, f2}{Intermediate derived values used in Eq. 5.26.}
#'
#' @examples
#' # Example: Compute N* with misclassification adjustment (pi01 = 0.1)
#' tdt_mssn_phenotype_misclassification_from_ncp(
#'   lambda_star = 7.8488,
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1,
#'   pi01 = 0.1
#' )
#'
#' @references
#' Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is not
#' a case: Effects of phenotype misclassification on power and sample size
#' requirements for the transmission disequilibrium test with affected child
#' trios. \emph{Human Heredity}, 67(4), 287--292.
#' \doi{10.1159/000194981}. PMID: 19172087.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' @seealso \code{\link{tdt_mssn}},
#' \code{\link{tdt_expected_transmission_probability}}, and
#' \code{\link{tdt_expected_nontransmission_probability}}.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_mssn_phenotype_misclassification_from_ncp <- function(
    lambda_star,
    pd,
    prev,
    R1, R2,
    delta_prime = 1,
    pi01 = 0,
    digits = 5
) {
  p_plus <- 1 - pd
  phi0   <- 1 - prev

  Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
  f0 <- prev / Z
  f1 <- R1 * f0
  f2 <- R2 * f0

  C <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0

  D   <- delta_prime * pd * p_plus
  DpT <- D * p_plus
  DpA <- D * pd

  # General adjusted probabilities in Eq. 5.26; Eq. 5.28 is a worked example.
  denom    <- prev + pi01 * phi0
  gT_star  <- (pd * p_plus) + (DpT * C * (1 - pi01)) / denom
  gNT_star <- (pd * p_plus) + (DpA * C * (pi01 - 1)) / denom

  N_star <- (lambda_star * (gT_star + gNT_star)) / (2 * (gT_star - gNT_star)^2)

  message("\n--- Transmission Disequilibrium Test (Trios) with Misclassification ---")
  message("Implements Eq. 5.26 (gT_star, gNT_star) and Eq. 5.27b for N_star")
  message("-----------------------------------------------------------------------")
  message("Parameters")
  message(sprintf("%-38s %10.3f", "Misclassification Rate (pi01):", pi01))
  message(sprintf("%-38s %10.3f  |  %s %8.3f",
                  "Allele Frequency (p_d):", pd, "Prevalence (phi1):", prev))
  message(sprintf("%-38s %10s", "Relative Risks (R1,R2):", paste0(R1, ", ", R2)))
  message(sprintf("%-38s %10.3f", "LD scale c (delta_prime):", delta_prime))
  message("-----------------------------------------------------------------------")
  message("Derived")
  message(sprintf("%-38s %10.5f  |  %s %8.5f",
                  "p_plus (transmitting allele freq):", p_plus, "phi0:", phi0))
  message(sprintf("%-38s %10.5f  |  f1: %8.5f  |  f2: %8.5f",
                  "f0:", f0, f1, f2))
  message(sprintf("%-38s %10.5f  |  DpT: %8.5f  |  DpA: %8.5f",
                  "C:", C, DpT, DpA))
  message("-----------------------------------------------------------------------")
  message("Expectations")
  message(sprintf("%-38s %10.5f", "gT_star:", gT_star))
  message(sprintf("%-38s %10.5f", "gNT_star:", gNT_star))
  message("-----------------------------------------------------------------------")
  message(sprintf("%-38s %10.5f", "lambda_star (TDT):", lambda_star))
  message(sprintf("%-38s %10.0f", "Required Number of Trios (N_star):", ceiling(N_star)))
  message("-----------------------------------------------------------------------")

  invisible(list(
    gT_star = gT_star,
    gNT_star = gNT_star,
    N_required = N_star,
    lambda_star = lambda_star,
    C = C, f0 = f0, f1 = f1, f2 = f2,
    pd = pd, p_plus = p_plus, prev = prev, phi0 = phi0,
    R1 = R1, R2 = R2, delta_prime = delta_prime, pi01 = pi01
  ))
}


#' Expected Transmission Probability (gT*) in the Transmission Disequilibrium Test (TDT)
#'
#' Computes the expected transmission probability (\eqn{g_T^*}) for the
#' Transmission Disequilibrium Test (TDT), incorporating effects of linkage
#' disequilibrium (LD), relative risks, and phenotype misclassification.
#' Implements Equation 5.24 from *Gordon et al.* (2020),
#' *Heterogeneity in Statistical Genetics*.
#'
#' @param pd Numeric. Frequency of the disease-associated allele.
#' @param prev Numeric. Disease prevalence (\eqn{\phi_1}).
#' @param R1 Numeric. Relative risk for heterozygotes.
#' @param R2 Numeric. Relative risk for homozygotes.
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter (\eqn{D'}) (default = 1).
#' @param pi01 Numeric. Misclassification rate (default = 0).
#' @param theta1 Numeric. Population allele frequency (defaults to `pd` if `NULL`).
#' @param digits Integer. Number of digits for printed results (default = 6).
#' @param verbose Logical. If `TRUE`, prints intermediate quantities (default = TRUE).
#'
#' @details
#' The function computes the penetrance-weighted genotype frequencies
#' \eqn{f_0,f_1,f_2} and the contrast term \eqn{C}, representing genotype-specific
#' effects. Linkage disequilibrium is modeled as \eqn{D = D' p_d p_+}.
#'
#' The expected transmission probability is given by:
#'
#' \deqn{
#' g_T^* =
#' \frac{
#'   p_d p_+ \phi_1 +
#'   D(p_+ - \theta_1)C +
#'   \pi_{01}(p_d p_+ \phi_0 + D(\theta_1 - p_+)C)
#' }{\phi_1 + \pi_{01}\phi_0}
#' }
#'
#' When \eqn{\pi_{01}=0}, this reduces to the standard non-misclassified case.
#' This is Eq. 5.24 in Chapter 5, Section 5.2.6 (p. 284).
#' \code{gT_star} is a probability, not an expected count.
#'
#' @return A list containing:
#' \item{gT_star}{Expected transmission probability (\eqn{g_T^*}).}
#' \item{C, D}{Contrast and LD terms.}
#' \item{f0, f1, f2}{Derived penetrance frequencies.}
#' \item{pd, phi1, phi0, theta1}{Input and derived parameters.}
#'
#' @examples
#' # Example 1 – no misclassification
#' tdt_expected_transmission_probability(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0
#' )
#'
#' # Example 2 – 10% misclassification
#' tdt_expected_transmission_probability(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0.1
#' )
#'
#' @references
#' Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is not
#' a case: Effects of phenotype misclassification on power and sample size
#' requirements for the transmission disequilibrium test with affected child
#' trios. \emph{Human Heredity}, 67(4), 287--292.
#' \doi{10.1159/000194981}. PMID: 19172087.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' @seealso \code{\link{tdt_expected_nontransmission_probability}},
#' \code{\link{tdt_expected_transmission_counts}}, and \code{\link{tdt_power}}.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_expected_transmission_probability <- function(pd, prev, R1, R2,
                            delta_prime = 1,
                            pi01 = 0,        # misclassification rate (default = 0)
                            theta1 = NULL,   # population allele freq (defaults to pd)
                            digits = 6,
                            verbose = TRUE) {
  # --- Derived quantities ---
  p_plus <- 1 - pd
  phi1   <- prev
  phi0   <- 1 - phi1
  if (is.null(theta1)) theta1 <- pd

  Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
  f0 <- prev / Z
  f1 <- R1 * f0
  f2 <- R2 * f0
  C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  D  <- delta_prime * pd * p_plus   # LD scale term

  # --- Equation 5.24 ---
  num <- (pd * p_plus * phi1) +
    D * (p_plus - theta1) * C +
    pi01 * (pd * p_plus * phi0 + D * (theta1 - p_plus) * C)
  denom <- phi1 + pi01 * phi0
  gT_star <- num / denom

  if (isTRUE(verbose)) {
    message("\n--- Transmission Disequilibrium Test: Expected Transmission (gT*) ---")
    message("Implements Equation 5.24 from Gordon et al. (2020)")
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.6f", "Allele Frequency (p_d):", pd))
    message(sprintf("%-38s %10.6f", "Prevalence (phi1):", phi1))
    message(sprintf("%-38s %10.6f", "LD Scale (delta_prime):", delta_prime))
    message(sprintf("%-38s %10.6f", "Misclassification Rate (pi01):", pi01))
    message(sprintf("%-38s %10.6f", "Population Allele Freq (theta1):", theta1))
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.6f", "Contrast Term (C):", C))
    message(sprintf("%-38s %10.6f", "LD Term (D):", D))
    message(sprintf("%-38s %10.6f", "Expected Transmission (gT*):", gT_star))
    message("-----------------------------------------------------------")
  }

  invisible(list(
    gT_star = gT_star,
    pd = pd,
    p_plus = p_plus,
    phi1 = phi1,
    phi0 = phi0,
    C = C,
    D = D,
    f0 = f0, f1 = f1, f2 = f2,
    delta_prime = delta_prime,
    pi01 = pi01,
    theta1 = theta1
  ))
}


#' Expected Non-Transmission Probability (gNT*) in the Transmission Disequilibrium Test (TDT)
#'
#' Computes the expected non-transmission probability (\eqn{g_{NT}^*}) for the
#' Transmission Disequilibrium Test (TDT) under a general model including linkage
#' disequilibrium (LD), relative risks, and phenotype misclassification.
#' Implements Equation 5.25 from *Gordon et al.* (2020),
#' *Heterogeneity in Statistical Genetics*.
#'
#' @param pd Numeric. Frequency of the disease-associated allele.
#' @param prev Numeric. Disease prevalence (\eqn{\phi_1}).
#' @param R1 Numeric. Relative risk for heterozygotes.
#' @param R2 Numeric. Relative risk for homozygotes.
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter (\eqn{D'}) (default = 1).
#' @param pi01 Numeric. Misclassification rate (default = 0).
#' @param theta1 Numeric. Population allele frequency (defaults to `pd` if `NULL`).
#' @param digits Integer. Number of digits for printing (default = 6).
#' @param verbose Logical. If `TRUE`, prints intermediate quantities (default = TRUE).
#'
#' @details
#' The function derives intermediate penetrance-weighted frequencies \eqn{f_0,f_1,f_2}
#' and the contrast term \eqn{C}, which captures differences in genotype contributions.
#' Linkage disequilibrium is represented by \eqn{D = D' p_d p_+}.
#' The expected non-transmission probability is then:
#'
#' \deqn{
#' g_{NT}^* =
#' \frac{
#'   p_d p_+ \phi_1 +
#'   D(\theta_1 - p_d)C +
#'   \pi_{01}(p_d p_+ \phi_0 + D(p_d - \theta_1)C)
#' }{\phi_1 + \pi_{01}\phi_0}
#' }
#'
#' When \eqn{\pi_{01}=0}, this reduces to the standard non-misclassified case.
#' This is Eq. 5.25 in Chapter 5, Section 5.2.6 (pp. 284--285).
#' \code{gNT_star} is a probability, not an expected count.
#'
#' @return A list containing:
#' \item{gNT_star}{Expected non-transmission probability (\eqn{g_{NT}^*}).}
#' \item{C, D}{Contrast and LD terms.}
#' \item{f0, f1, f2}{Derived penetrance frequencies.}
#' \item{pd, phi1, phi0, theta1}{Input and derived parameters.}
#'
#' @examples
#' # Example: compute gNT* under no misclassification
#' tdt_expected_nontransmission_probability(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0
#' )
#'
#' # Example: with 10% misclassification
#' tdt_expected_nontransmission_probability(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0.1
#' )
#'
#' @references
#' Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is not
#' a case: Effects of phenotype misclassification on power and sample size
#' requirements for the transmission disequilibrium test with affected child
#' trios. \emph{Human Heredity}, 67(4), 287--292.
#' \doi{10.1159/000194981}. PMID: 19172087.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' @seealso \code{\link{tdt_expected_transmission_probability}},
#' \code{\link{tdt_expected_transmission_counts}}, and \code{\link{tdt_power}}.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_expected_nontransmission_probability <- function(pd, prev, R1, R2,
                             delta_prime = 1,
                             pi01 = 0,       # misclassification rate
                             theta1 = NULL,   # population allele freq (defaults to pd)
                             digits = 6,
                             verbose = TRUE) {
  # --- Derived quantities ---
  p_plus <- 1 - pd
  phi1   <- prev
  phi0   <- 1 - phi1
  if (is.null(theta1)) theta1 <- pd

  Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
  f0 <- prev / Z
  f1 <- R1 * f0
  f2 <- R2 * f0
  C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  D  <- delta_prime * pd * p_plus   # LD scale term

  # --- Equation 5.25 ---
  num <- (pd * p_plus * phi1) +
    D * (theta1 - pd) * C +
    pi01 * (pd * p_plus * phi0 + D * (pd - theta1) * C)
  denom <- phi1 + pi01 * phi0
  gNT_star <- num / denom

  if (isTRUE(verbose)) {
    message("\n--- Transmission Disequilibrium Test: Expected Non-Transmission (gNT*) ---")
    message("Implements Equation 5.25 from Gordon et al. (2020)")
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.6f", "Allele Frequency (p_d):", pd))
    message(sprintf("%-38s %10.6f", "Prevalence (phi1):", phi1))
    message(sprintf("%-38s %10.6f", "LD Scale (delta_prime):", delta_prime))
    message(sprintf("%-38s %10.6f", "Misclassification Rate (pi01):", pi01))
    message(sprintf("%-38s %10.6f", "Population Allele Freq (theta1):", theta1))
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.6f", "Contrast Term (C):", C))
    message(sprintf("%-38s %10.6f", "LD Term (D):", D))
    message(sprintf("%-38s %10.6f", "Expected Non-Transmission (gNT*):", gNT_star))
    message("-----------------------------------------------------------")
  }

  invisible(list(
    gNT_star = gNT_star,
    pd = pd,
    p_plus = p_plus,
    phi1 = phi1,
    phi0 = phi0,
    C = C,
    D = D,
    f0 = f0, f1 = f1, f2 = f2,
    delta_prime = delta_prime,
    pi01 = pi01,
    theta1 = theta1
  ))
}


#' Expected Transmissions (ET*) and Non-Transmissions (ENT*) in the Transmission Disequilibrium Test (TDT)
#'
#' Computes the expected number of transmissions (\eqn{ET^*}) and non-transmissions (\eqn{ENT^*})
#' for a given number of affected trios under heterogeneity and linkage disequilibrium models.
#' Implements Equations 5.31-5.32 from *Gordon et al.* (2020),
#' *Heterogeneity in Statistical Genetics*.
#'
#' @param N_star Numeric. Number of affected trios.
#' @param pd Numeric. Frequency of the disease-associated allele.
#' @param prev Numeric. Disease prevalence (\eqn{\phi_1}).
#' @param R1 Numeric. Relative risk for heterozygotes.
#' @param R2 Numeric. Relative risk for homozygotes.
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter (\eqn{D'}) (default = 1).
#' @param pi Numeric in \eqn{[0,1]}. Linked/homogeneous trio fraction
#'   (default \code{1}); \eqn{1-\pi} is the heterogeneous fraction.
#' @param theta1 Numeric. Population allele frequency (defaults to `pd` if `NULL`).
#' @param digits Integer. Number of digits for printing (default = 6).
#' @param verbose Logical. If `TRUE`, prints intermediate quantities (default = TRUE).
#'
#' @details
#' The expected transmissions (\eqn{ET^*}) and non-transmissions (\eqn{ENT^*}) are computed as:
#'
#' \deqn{
#' ET^* = 2N^* \left[
#'   \pi \left( \frac{D(p_+ - \theta_1)C}{\phi_1} + p_d p_+ \right)
#'   + (1 - \pi)\left( \frac{D(p_+ - 0.5)C}{\phi_1} + p_d p_+ \right)
#' \right]
#' }
#'
#' \deqn{
#' ENT^* = 2N^* \left[
#'   \pi \left( \frac{D(\theta_1 - p_d)C}{\phi_1} + p_d p_+ \right)
#'   + (1 - \pi)\left( \frac{D(0.5 - p_d)C}{\phi_1} + p_d p_+ \right)
#' \right]
#' }
#'
#' These are Eqs. 5.31--5.32 in Chapter 5, Section 5.3.3 (pp. 293--294).
#' They are expected counts over \code{N_star} affected trios. Dividing by
#' \eqn{2N^*} gives the corresponding \eqn{g_T^*} and \eqn{g_{NT}^*}
#' probabilities. The canonical \code{heter_rate} used by \code{tdt_power()}
#' and \code{tdt_mssn()} equals \eqn{1-\pi}.
#'
#' @return A list containing:
#' \item{ET_star}{Expected transmissions (\eqn{ET^*}).}
#' \item{ENT_star}{Expected non-transmissions (\eqn{ENT^*}).}
#' \item{C, D}{Contrast and LD terms.}
#' \item{f0, f1, f2}{Derived penetrance frequencies.}
#' \item{pd, prev, R1, R2, pi, theta1}{Input parameters.}
#'
#' @examples
#' # Example 1: Homogeneous model (pi = 1)
#' tdt_expected_transmission_counts(
#'   N_star = 1000, pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi = 1
#' )
#'
#' # Example 2: Heterogeneous model (pi = 0.7)
#' tdt_expected_transmission_counts(
#'   N_star = 1000, pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi = 0.7
#' )
#'
#' @references
#' Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
#' (2009). Transmission disequilibrium test power and sample size in the
#' presence of locus heterogeneity. \emph{Statistical Applications in Genetics
#' and Molecular Biology}, 8, Article 44. \doi{10.2202/1544-6115.1501}.
#' PMID: 19883370.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' @seealso \code{\link{tdt_power}}, \code{\link{tdt_mssn}}, and
#' \code{\link{tdt_mssn_from_model}}.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_expected_transmission_counts <- function(N_star, pd, prev, R1, R2,
                                       delta_prime = 1, pi = 1,
                                       theta1 = NULL, digits = 6,
                                       verbose = TRUE) {
  # Derived quantities
  p_plus <- 1 - pd
  phi1   <- prev
  if (is.null(theta1)) theta1 <- pd

  # Penetrances and contrast
  Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
  f0 <- prev / Z
  f1 <- R1 * f0
  f2 <- R2 * f0
  C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  D  <- delta_prime * pd * p_plus

  # --- Equation 5.31: ET* ---
  ET_star <- 2 * N_star * (
    pi * ( (D * (p_plus - theta1) * C / phi1) + pd * p_plus ) +
      (1 - pi) * ( (D * (p_plus - 0.5) * C / phi1) + pd * p_plus )
  )

  # --- Equation 5.32: ENT* ---
  ENT_star <- 2 * N_star * (
    pi * ( (D * (theta1 - pd) * C / phi1) + pd * p_plus ) +
      (1 - pi) * ( (D * (0.5 - pd) * C / phi1) + pd * p_plus )
  )

  if (isTRUE(verbose)) {
    message("\n--- Transmission Disequilibrium Test: Expected ET* and ENT* ---")
    message("Implements Equations 5.31-5.32 (Heterogeneity model)")
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.0f", "Number of Trios (N*):", N_star))
    message(sprintf("%-38s %10.6f", "Allele Frequency (p_d):", pd))
    message(sprintf("%-38s %10.6f", "Prevalence (phi1):", phi1))
    message(sprintf("%-38s %10.6f", "LD Scale (delta_prime):", delta_prime))
    message(sprintf("%-38s %10.6f", "Heterogeneity Parameter (pi):", pi))
    message(sprintf("%-38s %10.6f", "Population Allele Freq (theta1):", theta1))
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.6f", "Contrast Term (C):", C))
    message(sprintf("%-38s %10.6f", "LD Term (D):", D))
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.6f", "Expected Transmissions (ET*):", ET_star))
    message(sprintf("%-38s %10.6f", "Expected Non-Transmissions (ENT*):", ENT_star))
    message("-----------------------------------------------------------")
  }

  invisible(list(
    ET_star = ET_star,
    ENT_star = ENT_star,
    C = C,
    D = D,
    f0 = f0, f1 = f1, f2 = f2,
    pd = pd, p_plus = p_plus,
    prev = prev,
    R1 = R1, R2 = R2,
    delta_prime = delta_prime,
    pi = pi,
    theta1 = theta1
  ))
}


#' Plot TDT Power Sensitivity to Phenotype Misclassification
#'
#' Plot TDT power vs phenotype misclassification rate (pi01),
#' holding heterogeneity fixed.
#'
#' @param N Integer. Number of affected trios.
#' @param pd Numeric in (0,1). High-risk allele frequency at the marker.
#' @param prev Numeric in (0,1). Disease prevalence (phi1).
#' @param R1,R2 Numeric. Genotype relative risks.
#' @param alpha Numeric. Significance level (default 0.05).
#' @param delta_prime Numeric. LD scale parameter (default 1).
#' @param misclass_seq Numeric vector. Sequence of misclassification rates.
#' @param heter_fixed Numeric. Heterogeneity rate held fixed (default 0).
#' @param title Character. Plot title.
#'
#' @details The x-axis is phenotype misclassification probability
#' \eqn{\pi_{01}} and the y-axis is the corresponding power from
#' \code{tdt_power()}. \code{N} is a fixed number of affected trios;
#' \code{heter_fixed} and all genetic-model parameters remain fixed while
#' \code{misclass_seq} is swept.
#'
#' @return A \code{ggplot} object, returned invisibly after it is printed.
#' @importFrom ggplot2 ggplot geom_line geom_point labs theme_bw aes
#' @importFrom rlang .data
#' @export
#' @examples
#' # Power sensitivity to phenotype misclassification
#' plot_tdt_power_phenotype_misclassification(
#'   N    = 600,
#'   pd   = 0.30,
#'   prev = 0.05,
#'   R1   = 1.5,
#'   R2   = 2.25,
#'   misclass_seq = c(0, 0.03, 0.06),
#'   heter_fixed  = 0
#' )
#'
#' @seealso \code{\link{tdt_power}},
#' \code{\link{plot_tdt_power_locus_heterogeneity}}, and
#' \code{\link{plot_tdt_power}}.

plot_tdt_power_phenotype_misclassification <- function(
    pd, prev, R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    N,
    misclass_seq  = seq(0, 0.15, by = 0.01),
    heter_fixed   = 0,
    title = "TDT power vs misclassification (pi01)"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.", call. = FALSE)
  }

  power_mis <- vapply(misclass_seq, function(mis) {
    res <- tdt_power(
      N = N,
      pd = pd, prev = prev, R1 = R1, R2 = R2,
      alpha = alpha, delta_prime = delta_prime,
      misclass_rate = mis,
      heter_rate   = heter_fixed,
      verbose = FALSE
    )
    res$power$misclassification
  }, numeric(1))

  df_mis <- data.frame(
    misclass_rate = misclass_seq,
    power = power_mis
  )

  p <- ggplot2::ggplot(df_mis, ggplot2::aes(x = .data$misclass_rate, y = .data$power)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Misclassification rate (pi01)",
      y = "Power",
      title = title
    ) +
    ggplot2::theme_bw()

  print(p)
  invisible(p)
}


#' Plot TDT Power Sensitivity to Locus Heterogeneity
#'
#' Plot TDT power vs locus heterogeneity rate (1 - pi),
#' holding misclassification fixed.
#' @param N Number of affected trios.
#' @param pd Disease allele frequency.
#' @param prev Disease prevalence.
#' @param R1 Genotype relative risk for heterozygotes.
#' @param R2 Genotype relative risk for risk homozygotes.
#' @param alpha Significance level.
#' @param delta_prime Linkage disequilibrium scaling parameter.
#' @param title Plot title.
#' @param heter_seq Numeric vector. Sequence of heterogeneity rates.
#' @param misclass_fixed Numeric. Misclassification rate held fixed (default 0).
#' (Other params same meaning as in plot_tdt_power_phenotype_misclassification.)
#'
#' @details The x-axis is the heterogeneous trio fraction \eqn{1-\pi}; the
#' homogeneous/linked fraction is \eqn{\pi}. The y-axis is heterogeneity-scenario
#' power from \code{tdt_power()}. \code{N}, \code{misclass_fixed}, and all
#' genetic-model parameters remain fixed while \code{heter_seq} is swept.
#'
#' @return A \code{ggplot} object, returned invisibly after it is printed.
#' @importFrom ggplot2 ggplot geom_line geom_point labs theme_bw aes
#' @importFrom rlang .data
#' @export
#' @examples
#' # Power sensitivity to locus heterogeneity
#' plot_tdt_power_locus_heterogeneity(
#'   N    = 600,
#'   pd   = 0.30,
#'   prev = 0.05,
#'   R1   = 1.5,
#'   R2   = 2.25,
#'   heter_seq      = c(0, 0.1, 0.2),
#'   misclass_fixed = 0
#' )
#'
#' @seealso \code{\link{tdt_power}},
#' \code{\link{plot_tdt_power_phenotype_misclassification}}, and
#' \code{\link{plot_tdt_power}}.

plot_tdt_power_locus_heterogeneity <- function(
    pd, prev, R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    N,
    heter_seq     = seq(0, 0.50, by = 0.05),
    misclass_fixed = 0,
    title = "TDT power vs heterogeneity (1 - pi)"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.", call. = FALSE)
  }

  power_het <- vapply(heter_seq, function(h) {
    res <- tdt_power(
      N = N,
      pd = pd, prev = prev, R1 = R1, R2 = R2,
      alpha = alpha, delta_prime = delta_prime,
      misclass_rate = misclass_fixed,
      heter_rate   = h,
      verbose = FALSE
    )
    res$power$heterogeneity
  }, numeric(1))

  df_het <- data.frame(
    heter_rate = heter_seq,
    power = power_het
  )

  p <- ggplot2::ggplot(df_het, ggplot2::aes(x = .data$heter_rate, y = .data$power)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Heterogeneity rate (1 - pi)",
      y = "Power",
      title = title
    ) +
    ggplot2::theme_bw()

  print(p)
  invisible(p)
}


#' Plot TDT MSSN Sensitivity to Phenotype Misclassification
#'
#' Plot required number of trios vs misclassification rate (pi01),
#' holding heterogeneity fixed.
#' @param pd Disease allele frequency.
#' @param prev Disease prevalence.
#' @param R1 Genotype relative risk for heterozygotes.
#' @param R2 Genotype relative risk for risk homozygotes.
#' @param alpha Significance level.
#' @param delta_prime Linkage disequilibrium scaling parameter.
#' @param title Plot title.
#' @param target_power Numeric. Desired power for sample size calculation.
#' @param misclass_seq Numeric vector. Sequence of misclassification rates.
#' @param heter_fixed Numeric. Heterogeneity rate held fixed (default 0).
#' (Other params same meaning as in plot_tdt_power_phenotype_misclassification.)
#'
#' @details The x-axis is phenotype misclassification probability
#' \eqn{\pi_{01}}. The y-axis is MSSN, the required number of affected trios
#' returned by \code{tdt_mssn()} for fixed \code{target_power}.
#' \code{heter_fixed} and all genetic-model parameters remain fixed.
#'
#' @return A \code{ggplot} object, returned invisibly after it is printed.
#' @importFrom ggplot2 ggplot geom_line geom_point labs theme_bw aes
#' @importFrom rlang .data
#' @export
#' @examples
#' # Required sample size sensitivity to phenotype misclassification
#' plot_tdt_mssn_phenotype_misclassification(
#'   target_power = 0.80,
#'   pd   = 0.30,
#'   prev = 0.05,
#'   R1   = 1.5,
#'   R2   = 2.25,
#'   misclass_seq = c(0, 0.03, 0.06),
#'   heter_fixed  = 0
#' )
#'
#' @seealso \code{\link{tdt_mssn}},
#' \code{\link{plot_tdt_mssn_locus_heterogeneity}}, and
#' \code{\link{plot_tdt_mssn}}.

plot_tdt_mssn_phenotype_misclassification <- function(
    pd, prev, R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    target_power,
    misclass_seq  = seq(0, 0.15, by = 0.01),
    heter_fixed   = 0,
    title = "Required trios vs misclassification (pi01)"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.", call. = FALSE)
  }

  N_mis <- vapply(misclass_seq, function(mis) {
    res <- tdt_mssn(
      target_power = target_power,
      pd = pd, prev = prev, R1 = R1, R2 = R2,
      alpha = alpha, delta_prime = delta_prime,
      misclass_rate = mis,
      heter_rate   = heter_fixed,
      verbose = FALSE
    )
    res$N$misclassification
  }, numeric(1))

  df_mis <- data.frame(
    misclass_rate = misclass_seq,
    N_required = N_mis
  )

  p <- ggplot2::ggplot(df_mis, ggplot2::aes(x = .data$misclass_rate, y = .data$N_required)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Misclassification rate (pi01)",
      y = "Required number of trios",
      title = title
    ) +
    ggplot2::theme_bw()

  print(p)
  invisible(p)
}


#' Plot TDT MSSN Sensitivity to Locus Heterogeneity
#'
#' Plot required number of trios vs heterogeneity rate (1 - pi),
#' holding misclassification fixed.
#' @param target_power Desired statistical power.
#' @param pd Disease allele frequency.
#' @param prev Disease prevalence.
#' @param R1 Genotype relative risk for heterozygotes.
#' @param R2 Genotype relative risk for risk homozygotes.
#' @param alpha Significance level.
#' @param delta_prime Linkage disequilibrium scaling parameter.
#' @param title Plot title.
#' @param heter_seq Numeric vector. Sequence of heterogeneity rates.
#' @param misclass_fixed Numeric. Misclassification rate held fixed (default 0).
#' (Other params same meaning as in plot_tdt_power_phenotype_misclassification.)
#'
#' @details The x-axis is the heterogeneous trio fraction \eqn{1-\pi}; the
#' homogeneous/linked fraction is \eqn{\pi}. The y-axis is MSSN, the required
#' number of affected trios returned by \code{tdt_mssn()} for fixed
#' \code{target_power}. \code{misclass_fixed} and genetic-model parameters
#' remain fixed.
#'
#' @return A \code{ggplot} object, returned invisibly after it is printed.
#' @importFrom ggplot2 ggplot geom_line geom_point labs theme_bw aes
#' @importFrom rlang .data
#' @export
#' @examples
#' # Required sample size sensitivity to locus heterogeneity
#' plot_tdt_mssn_locus_heterogeneity(
#'   target_power = 0.80,
#'   pd   = 0.30,
#'   prev = 0.05,
#'   R1   = 1.5,
#'   R2   = 2.25,
#'   heter_seq      = c(0, 0.1, 0.2),
#'   misclass_fixed = 0
#' )
#'
#' @seealso \code{\link{tdt_mssn}},
#' \code{\link{plot_tdt_mssn_phenotype_misclassification}}, and
#' \code{\link{plot_tdt_mssn}}.

plot_tdt_mssn_locus_heterogeneity <- function(
    pd, prev, R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    target_power,
    heter_seq     = seq(0, 0.50, by = 0.05),
    misclass_fixed = 0,
    title = "Required trios vs heterogeneity (1 - pi)"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.", call. = FALSE)
  }

  N_het <- vapply(heter_seq, function(h) {
    res <- tdt_mssn(
      target_power = target_power,
      pd = pd, prev = prev, R1 = R1, R2 = R2,
      alpha = alpha, delta_prime = delta_prime,
      misclass_rate = misclass_fixed,
      heter_rate   = h,
      verbose = FALSE
    )
    res$N$heterogeneity
  }, numeric(1))

  df_het <- data.frame(
    heter_rate = heter_seq,
    N_required = N_het
  )

  p <- ggplot2::ggplot(df_het, ggplot2::aes(x = .data$heter_rate, y = .data$N_required)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Heterogeneity rate (1 - pi)",
      y = "Required number of trios",
      title = title
    ) +
    ggplot2::theme_bw()

  print(p)
  invisible(p)
}

#' Plot a 3D TDT Power Surface for Phenotype Misclassification
#'
#' This function requires the optional package \pkg{plotly}.
#' @param pd_seq Vector of allele frequencies.
#' @param misclass_seq Vector of misclassification rates (pi01).
#' @param N Number of affected trios.
#' @param prev Disease prevalence.
#' @param R1,R2 Genotype relative risks.
#' @param alpha Significance level.
#' @param delta_prime LD scale.
#' @param heter_rate Heterogeneity rate to hold fixed while varying misclassification.
#' @param title Plot title.
#'
#' @details The x-axis is risk-allele frequency \code{pd}, the y-axis is
#' phenotype misclassification probability \eqn{\pi_{01}}, and the surface
#' height is misclassification-scenario TDT power from \code{tdt_power()}.
#' The affected-trio count \code{N}, \code{heter_rate}, and other model
#' parameters remain fixed. This transitional function is restricted to these
#' two swept axes and requires the optional package \pkg{plotly}.
#'
#' @return An object inheriting from \code{"plotly"} and \code{"htmlwidget"}.
#' @export
#'
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   p <- plot_tdt_power_phenotype_misclassification_3d(
#'     pd_seq = c(0.2, 0.3, 0.4),
#'     misclass_seq = c(0, 0.03, 0.06),
#'     N = 600, prev = 0.05, R1 = 1.5, R2 = 2.25,
#'     heter_rate = 0
#'   )
#' }
#'
#' @seealso \code{\link{tdt_power}},
#' \code{\link{plot_tdt_power_phenotype_misclassification}}, and
#' \code{\link{plot_tdt_power_locus_heterogeneity_3d}}.
plot_tdt_power_phenotype_misclassification_3d <- function(
    pd_seq,
    misclass_seq,
    N,
    prev,
    R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    heter_rate = 0,
    title = "TDT power vs allele frequency and misclassification"
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this plot.", call. = FALSE)
  }

  z <- vapply(misclass_seq, function(mis) {
    vapply(pd_seq, function(pd) {
      res <- tdt_power(
        N = N, pd = pd, prev = prev, R1 = R1, R2 = R2,
        alpha = alpha, delta_prime = delta_prime,
        misclass_rate = mis,
        heter_rate = heter_rate,
        verbose = FALSE
      )
      res$power$misclassification
    }, numeric(1))
  }, numeric(length(pd_seq)))

  z <- matrix(z, nrow = length(misclass_seq), ncol = length(pd_seq), byrow = TRUE)

  plotly::plot_ly(
    x = pd_seq,
    y = misclass_seq,
    z = z,
    type = "surface"
  ) |>
    plotly::layout(
      title = list(text = title),
      scene = list(
        xaxis = list(title = "Allele frequency (p_d)"),
        yaxis = list(title = "Misclassification rate (pi01)"),
        zaxis = list(title = "Power")
      )
    )
}


#' Plot a 3D TDT Power Surface for Locus Heterogeneity
#'
#' This function requires the optional package \pkg{plotly}.
#' @param pd_seq Vector of allele frequencies.
#' @param heter_seq Vector of heterogeneity rates (1 - pi).
#' @param N Number of affected trios.
#' @param prev Disease prevalence.
#' @param R1,R2 Genotype relative risks.
#' @param alpha Significance level.
#' @param delta_prime LD scale.
#' @param misclass_rate Misclassification rate to hold fixed while varying heterogeneity.
#' @param title Plot title.
#'
#' @details The x-axis is risk-allele frequency \code{pd}, the y-axis is the
#' heterogeneous trio fraction \eqn{1-\pi}, and the surface height is
#' heterogeneity-scenario TDT power from \code{tdt_power()}. The affected-trio
#' count \code{N}, \code{misclass_rate}, and other model parameters remain
#' fixed. This transitional function is restricted to these swept axes.
#'
#' @return An object inheriting from \code{"plotly"} and \code{"htmlwidget"}.
#' @export
#'
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   p <- plot_tdt_power_locus_heterogeneity_3d(
#'     pd_seq = c(0.2, 0.3, 0.4),
#'     heter_seq = c(0, 0.1, 0.2),
#'     N = 600, prev = 0.05, R1 = 1.5, R2 = 2.25,
#'     misclass_rate = 0
#'   )
#' }
#'
#' @seealso \code{\link{tdt_power}},
#' \code{\link{plot_tdt_power_locus_heterogeneity}}, and
#' \code{\link{plot_tdt_power_phenotype_misclassification_3d}}.
plot_tdt_power_locus_heterogeneity_3d <- function(
    pd_seq,
    heter_seq,
    N,
    prev,
    R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    misclass_rate = 0,
    title = "TDT power vs allele frequency and heterogeneity"
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this plot.", call. = FALSE)
  }

  z <- vapply(heter_seq, function(h) {
    vapply(pd_seq, function(pd) {
      res <- tdt_power(
        N = N, pd = pd, prev = prev, R1 = R1, R2 = R2,
        alpha = alpha, delta_prime = delta_prime,
        misclass_rate = misclass_rate,
        heter_rate = h,
        verbose = FALSE
      )
      res$power$heterogeneity
    }, numeric(1))
  }, numeric(length(pd_seq)))

  z <- matrix(z, nrow = length(heter_seq), ncol = length(pd_seq), byrow = TRUE)

  plotly::plot_ly(
    x = pd_seq,
    y = heter_seq,
    z = z,
    type = "surface"
  ) |>
    plotly::layout(
      title = list(text = title),
      scene = list(
        xaxis = list(title = "Allele frequency (p_d)"),
        yaxis = list(title = "Heterogeneity rate (1 - pi)"),
        zaxis = list(title = "Power")
      )
    )
}


#' Plot a 3D TDT MSSN Surface for Locus Heterogeneity
#'
#' This function requires the optional package \pkg{plotly}.
#' @param pd_seq Vector of allele frequencies.
#' @param heter_seq Vector of heterogeneity rates (1 - pi).
#' @param target_power Desired power.
#' @param prev Disease prevalence.
#' @param R1,R2 Genotype relative risks.
#' @param alpha Significance level.
#' @param delta_prime LD scale.
#' @param misclass_rate Misclassification rate held fixed.
#' @param ceiling_N Logical; if TRUE (default) plots ceiling(N).
#' @param title Plot title.
#'
#' @details The x-axis is risk-allele frequency \code{pd}, the y-axis is the
#' heterogeneous trio fraction \eqn{1-\pi}, and the height is MSSN in affected
#' trios from \code{tdt_mssn()} for fixed \code{target_power}. Other model
#' parameters, including \code{misclass_rate}, remain fixed. With
#' \code{ceiling_N = TRUE}, displayed heights are integer ceilings.
#'
#' @return An object inheriting from \code{"plotly"} and \code{"htmlwidget"}.
#' @export
#'
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   p <- plot_tdt_mssn_locus_heterogeneity_3d(
#'     pd_seq = c(0.2, 0.3, 0.4),
#'     heter_seq = c(0, 0.1, 0.2),
#'     target_power = 0.80, prev = 0.05, R1 = 1.5, R2 = 2.25,
#'     misclass_rate = 0
#'   )
#' }
#'
#' @seealso \code{\link{tdt_mssn}},
#' \code{\link{plot_tdt_mssn_locus_heterogeneity}}, and
#' \code{\link{plot_tdt_mssn_phenotype_misclassification_3d}}.
plot_tdt_mssn_locus_heterogeneity_3d <- function(
    pd_seq,
    heter_seq,
    target_power,
    prev,
    R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    misclass_rate = 0,
    ceiling_N = TRUE,
    title = "Required trios vs allele frequency and heterogeneity"
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this plot.", call. = FALSE)
  }

  z <- vapply(heter_seq, function(h) {
    vapply(pd_seq, function(pd) {
      res <- tdt_mssn(
        target_power = target_power,
        pd = pd, prev = prev, R1 = R1, R2 = R2,
        alpha = alpha, delta_prime = delta_prime,
        misclass_rate = misclass_rate,
        heter_rate = h,
        verbose = FALSE
      )
      n <- res$N$heterogeneity
      if (isTRUE(ceiling_N)) ceiling(n) else n
    }, numeric(1))
  }, numeric(length(pd_seq)))

  z <- matrix(z, nrow = length(heter_seq), ncol = length(pd_seq), byrow = TRUE)

  plotly::plot_ly(
    x = pd_seq,
    y = heter_seq,
    z = z,
    type = "surface"
  ) |>
    plotly::layout(
      title = list(text = title),
      scene = list(
        xaxis = list(title = "Allele frequency (p_d)"),
        yaxis = list(title = "Heterogeneity rate (1 - pi)"),
        zaxis = list(title = "Required number of trios (N*)")
      )
    )
}

#' Plot a 3D TDT MSSN Surface for Phenotype Misclassification
#'
#' This function requires the optional package \pkg{plotly}.
#' @param pd_seq Vector of allele frequencies.
#' @param misclass_seq Vector of misclassification rates (pi01).
#' @param target_power Desired power.
#' @param prev Disease prevalence.
#' @param R1,R2 Genotype relative risks.
#' @param alpha Significance level.
#' @param delta_prime LD scale.
#' @param heter_rate Heterogeneity rate held fixed.
#' @param ceiling_N Logical; if TRUE (default) plots ceiling(N).
#' @param title Plot title.
#'
#' @details The x-axis is risk-allele frequency \code{pd}, the y-axis is
#' phenotype misclassification probability \eqn{\pi_{01}}, and the height is
#' MSSN in affected trios from \code{tdt_mssn()} for fixed
#' \code{target_power}. Other model parameters, including \code{heter_rate},
#' remain fixed. With \code{ceiling_N = TRUE}, displayed heights are integer
#' ceilings.
#'
#' @return An object inheriting from \code{"plotly"} and \code{"htmlwidget"}.
#' @export
#'
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   p <- plot_tdt_mssn_phenotype_misclassification_3d(
#'     pd_seq = c(0.2, 0.3, 0.4),
#'     misclass_seq = c(0, 0.03, 0.06),
#'     target_power = 0.80, prev = 0.05, R1 = 1.5, R2 = 2.25,
#'     heter_rate = 0
#'   )
#' }
#'
#' @seealso \code{\link{tdt_mssn}},
#' \code{\link{plot_tdt_mssn_phenotype_misclassification}}, and
#' \code{\link{plot_tdt_mssn_locus_heterogeneity_3d}}.
plot_tdt_mssn_phenotype_misclassification_3d <- function(
    pd_seq,
    misclass_seq,
    target_power,
    prev,
    R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    heter_rate = 0,
    ceiling_N = TRUE,
    title = "Required trios vs allele frequency and misclassification"
) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this plot.", call. = FALSE)
  }

  z <- vapply(misclass_seq, function(mis) {
    vapply(pd_seq, function(pd) {
      res <- tdt_mssn(
        target_power = target_power,
        pd = pd, prev = prev, R1 = R1, R2 = R2,
        alpha = alpha, delta_prime = delta_prime,
        misclass_rate = mis,
        heter_rate = heter_rate,
        verbose = FALSE
      )
      n <- res$N$misclassification
      if (isTRUE(ceiling_N)) ceiling(n) else n
    }, numeric(1))
  }, numeric(length(pd_seq)))

  z <- matrix(z, nrow = length(misclass_seq), ncol = length(pd_seq), byrow = TRUE)

  plotly::plot_ly(
    x = pd_seq,
    y = misclass_seq,
    z = z,
    type = "surface"
  ) |>
    plotly::layout(
      title = list(text = title),
      scene = list(
        xaxis = list(title = "Allele frequency (p_d)"),
        yaxis = list(title = "Misclassification rate (pi01)"),
        zaxis = list(title = "Required number of trios (N*)")
      )
    )
}
