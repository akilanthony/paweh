#' Transmission Disequilibrium Test (TDT) Power from Expected Transmissions and Non-Transmissions
#'
#' Computes the statistical power of the Transmission Disequilibrium Test (TDT)
#' given the expected number of transmissions (ET) and non-transmissions (ENT)
#' under a specified significance level. Implements Equation 1.25 from
#' *Gordon et al. (2020), Heterogeneity in Statistical Genetics*.

#' @param ET Numeric. Expected number of transmissions.
#' @param ENT Numeric. Expected number of non-transmissions.
#' @param alpha Numeric. Significance level (default = 0.05).
#'
#'
#' @details
#' The function calculates the non-centrality parameter and statistical power for the TDT
#' using the chi-square distribution with 1 degree of freedom.
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
#' tdt_power_from_ET_ENT(ET = 140, ENT = 100, alpha = 0.05)
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#'@export

tdt_power_from_ET_ENT <- function(ET, ENT, alpha = 0.05) {
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
#' disease prevalence, and the number of affected trios. Implements Equation 1.25
#' from *Gordon et al. (2020), Heterogeneity in Statistical Genetics*.
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
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
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


#' Required Number of Trios for Desired Power in the Transmission Disequilibrium Test (TDT)
#'
#' Computes the required number of affected trios (\eqn{N^*}) needed to achieve
#' a specified statistical power in the Transmission Disequilibrium Test (TDT),
#' given model parameters for allele frequency, relative risks, disease prevalence,
#' and heterogeneity. Implements Equation 5.34b from
#' *Gordon et al. (2020), Heterogeneity in Statistical Genetics*.
#'
#' @param power Numeric. Desired power (e.g., 0.8).
#' @param alpha Numeric. Significance level (e.g., 0.05).
#' @param df Integer. Degrees of freedom (typically 1 for TDT).
#' @param pd Numeric. Frequency of the disease-associated allele.
#' @param prev Numeric. Disease prevalence.
#' @param R1 Numeric. Relative risk for heterozygotes.
#' @param R2 Numeric. Relative risk for homozygotes.
#' @param delta_prime Numeric. Linkage disequilibrium (LD) scale factor (default = 1).
#' @param pi Numeric. Heterogeneity parameter, where 1 represents full homogeneity
#' and values between 0–1 allow for mixed genetic effects (default = 1).
#'
#' @details
#' This function determines the non-centrality parameter (\eqn{\lambda^*}) via
#' root-finding (\code{uniroot}) such that the test power equals the desired level.
#' It then computes the expected transmission (\eqn{gT^*}) and non-transmission
#' (\eqn{gNT^*}) probabilities, followed by the required number of trios using:
#' \deqn{N^* = \frac{\lambda^*}{2} \frac{(gT^* + gNT^*)}{(gT^* - gNT^*)^2}}
#'
#' The expected transmission and non-transmission components are calculated
#' under allele frequency and penetrance model assumptions (see Eq. 5.34b).
#'
#' @return A list containing:
#' \item{lambda_star}{Non-centrality parameter (\eqn{\lambda^*}).}
#' \item{gT_star}{Expected transmission probability.}
#' \item{gNT_star}{Expected non-transmission probability.}
#' \item{N_star}{Required number of trios.}
#'
#' @examples
#' # Example: compute required trios for 80% power at alpha = 0.05
#' tdt_required_trios(
#'   power = 0.8, alpha = 0.05, df = 1,
#'   pd = 0.25, prev = 0.005, R1 = 2, R2 = 2,
#'   delta_prime = 1, pi = 1
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_required_trios <- function(power, alpha, df,
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


#' Required Number of Trios in the Transmission Disequilibrium Test (TDT) with Misclassification
#'
#' Computes the expected transmissions (\eqn{gT^*}) and non-transmissions (\eqn{gNT^*})
#' as well as the required number of trios (\eqn{N^*}) for a specified non-centrality parameter
#' (\eqn{\lambda^*}) under a misclassification model. Implements Equations 5.26-5.28
#' and Equation 5.27b from *Gordon et al. (2020), Heterogeneity in Statistical Genetics*.
#'
#' @param lambda_star Numeric. Non-centrality parameter (\eqn{\lambda^*}) derived from
#' desired power (e.g., from \code{tdt_required_trios()}).
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
#' \item{C, f0, f1, f2}{Intermediate derived values used in Eq. 5.26-5.28.}
#'
#' @examples
#' # Example: Compute N* with misclassification adjustment (pi01 = 0.1)
#' tdt_required_trios_misclass(
#'   lambda_star = 7.8488,
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1,
#'   pi01 = 0.1
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_required_trios_misclass <- function(
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

  # Eqs. (5.28a,b)
  denom    <- prev + pi01 * phi0
  gT_star  <- (pd * p_plus) + (DpT * C * (1 - pi01)) / denom
  gNT_star <- (pd * p_plus) + (DpA * C * (pi01 - 1)) / denom

  N_star <- (lambda_star * (gT_star + gNT_star)) / (2 * (gT_star - gNT_star)^2)

  message("\n--- Transmission Disequilibrium Test (Trios) with Misclassification ---")
  message("Implements Eqs. 5.26-5.28 (gT_star, gNT_star) and Eq. 5.27b for N_star")
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
#'
#' @return A list containing:
#' \item{gT_star}{Expected transmission probability (\eqn{g_T^*}).}
#' \item{C, D}{Contrast and LD terms.}
#' \item{f0, f1, f2}{Derived penetrance frequencies.}
#' \item{pd, phi1, phi0, theta1}{Input and derived parameters.}
#'
#' @examples
#' # Example 1 – no misclassification
#' tdt_expected_gT(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0
#' )
#'
#' # Example 2 – 10% misclassification
#' tdt_expected_gT(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0.1
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_expected_gT <- function(pd, prev, R1, R2,
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
#'
#' @return A list containing:
#' \item{gNT_star}{Expected non-transmission probability (\eqn{g_{NT}^*}).}
#' \item{C, D}{Contrast and LD terms.}
#' \item{f0, f1, f2}{Derived penetrance frequencies.}
#' \item{pd, phi1, phi0, theta1}{Input and derived parameters.}
#'
#' @examples
#' # Example: compute gNT* under no misclassification
#' tdt_expected_gNT(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0
#' )
#'
#' # Example: with 10% misclassification
#' tdt_expected_gNT(
#'   pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi01 = 0.1
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_expected_gNT <- function(pd, prev, R1, R2,
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
#' @param pi Numeric. Heterogeneity parameter controlling the proportion of homogeneous trios (default = 1).
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
#' These expressions generalize the basic TDT expectations (Eq. 5.24-5.25) to include
#' heterogeneity across trios through the parameter \eqn{\pi}.
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
#' tdt_expected_transmissions(
#'   N_star = 1000, pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi = 1
#' )
#'
#' # Example 2: Heterogeneous model (pi = 0.7)
#' tdt_expected_transmissions(
#'   N_star = 1000, pd = 0.25, prev = 0.005,
#'   R1 = 2, R2 = 2,
#'   delta_prime = 1, pi = 0.7
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export

tdt_expected_transmissions <- function(N_star, pd, prev, R1, R2,
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

#' Transmission Disequilibrium Test Power with Error Components
#'
#' Computes the power of the Transmission Disequilibrium Test (TDT) for a
#' fixed number of affected trios under three scenarios:
#' (i) no error, (ii) phenotype misclassification only, and
#' (iii) locus heterogeneity only.  Misclassification and heterogeneity are
#' modeled using the framework in Gordon et al. (2020), Chapter 5.
#'
#' @param N Integer. Number of affected trios.
#' @param pd Numeric in (0,1). Frequency of the disease (high-risk) allele at
#'   the marker locus.
#' @param prev Numeric in (0,1). Disease prevalence (\eqn{\phi_1}).
#' @param R1 Numeric \eqn{> 0}. Genotype relative risk for heterozygotes.
#' @param R2 Numeric \eqn{> 0}. Genotype relative risk for homozygotes.
#' @param alpha Numeric in (0,1). Significance level for the TDT
#'   (default \code{alpha = 0.05}).
#' @param delta_prime Numeric. Linkage disequilibrium (LD) scale parameter
#'   \eqn{D'} (default \code{1}).
#' @param misclass_rate Numeric in \eqn{[0,1)}. Phenotype misclassification
#'   rate for controls (\eqn{\pi_{01}}). A value of \code{0} corresponds to
#'   no misclassification.
#' @param heter_rate Numeric in \eqn{[0,1)}. Proportion of trios whose
#'   affection status is \emph{not} due to the locus of interest
#'   (\eqn{1 - \pi}). A value of \code{0} corresponds to complete
#'   homogeneity.
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted
#'   summary of non-centrality parameters, power, and expected transmission
#'   probabilities.
#'
#' @details
#' Penetrances \eqn{f_0, f_1, f_2} for genotypes with 0, 1, and 2 risk alleles
#' are derived from \code{prev}, \code{R1}, \code{R2}, and \code{pd} via
#' the standard normalization:
#' \deqn{
#' Z = (1 - p_d)^2 + 2 p_d (1 - p_d) R_1 + p_d^2 R_2, \quad
#' f_0 = \phi_1 / Z, \quad f_1 = R_1 f_0, \quad f_2 = R_2 f_0.
#' }
#' The contrast term \eqn{C} and LD term \eqn{D = D' p_d (1-p_d)} are used
#' to construct the expected transmission and non-transmission probabilities
#' under:
#' \enumerate{
#'   \item No error (baseline TDT model).
#'   \item Misclassification only, using Equations 5.26–5.28.
#'   \item Heterogeneity only, using Equation 5.34.
#' }
#' For each scenario, the non-centrality parameter
#' \eqn{\lambda = 2N (g_T^* - g_{NT}^*)^2 / (g_T^* + g_{NT}^*)} is computed
#' and power is obtained from the non-central chi-square distribution with
#' one degree of freedom.
#'
#' If the derived penetrances \code{f1} or \code{f2} exceed 1, a warning is
#' issued, indicating that the combination of \code{prev}, \code{R1},
#' \code{R2}, and \code{pd} is not coherent for a penetrance model.
#'
#' @return
#' An object of class \code{"tdt_power_full"}: a list with components
#' \describe{
#'   \item{alpha}{Significance level used.}
#'   \item{N}{Number of affected trios.}
#'   \item{lambda}{List of non-centrality parameters with elements
#'     \code{no_error}, \code{misclassification}, and \code{heterogeneity}.}
#'   \item{power}{List of powers for each scenario:
#'     \code{no_error}, \code{misclassification}, and \code{heterogeneity}.}
#'   \item{power_loss}{List with elements \code{misclassification} and
#'     \code{heterogeneity}, giving the loss in power relative to the
#'     no-error case.}
#'   \item{gT_star}{List of expected transmission probabilities
#'     (\eqn{g_T^*}) for each scenario.}
#'   \item{gNT_star}{List of expected non-transmission probabilities
#'     (\eqn{g_{NT}^*}) for each scenario.}
#'   \item{ET, ENT}{Lists of expected transmissions and non-transmissions
#'     under each scenario.}
#'   \item{model_parameters}{List of the input model parameters
#'     (\code{pd}, \code{prev}, \code{R1}, \code{R2}, \code{delta_prime},
#'     \code{misclass_rate}, \code{heter_rate}).}
#' }
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' \emph{Heterogeneity in Statistical Genetics}. Springer Nature.
#'
#' @examples
#' # Power for N = 600 trios under a modest effect
#' tdt_power_full(
#'   N            = 600,
#'   pd           = 0.30,
#'   prev         = 0.05,
#'   R1           = 1.5,
#'   R2           = 2.25,
#'   alpha        = 0.05,
#'   delta_prime  = 1,
#'   misclass_rate = 0.01,
#'   heter_rate    = 0.10,
#'   verbose      = TRUE
#' )
#'
#' # Extract power without printing
#' tdt_power_full(
#'   N = 600, pd = 0.3, prev = 0.05,
#'   R1 = 1.5, R2 = 2.25,
#'   verbose = FALSE
#' )$power$no_error
#'
#' @export
tdt_power_full <- function(
    N,                 # number of affected trios
    pd,                # disease/high-risk allele frequency at marker
    prev,              # prevalence (phi1)
    R1, R2,            # genotype relative risks
    alpha = 0.05,
    delta_prime = 1,       # LD scale
    misclass_rate = 0.01,  # phenotype misclassification (pi01)
    heter_rate   = 0.01,   # proportion of trios NOT due to the locus (1 - pi)
    verbose = TRUE
) {
  ## ---- basic argument checks ----
  if (misclass_rate < 0 || misclass_rate >= 1)
    stop("misclass_rate must be in [0, 1).")
  if (heter_rate < 0 || heter_rate >= 1)
    stop("heter_rate must be in [0, 1).")

  ## ---------- internal helper: gT* and gNT* with misclassification ----------
  calc_gTgNT_misclass <- function(pd, prev, R1, R2,
                                  delta_prime, pi01) {
    p_plus <- 1 - pd
    phi1   <- prev
    phi0   <- 1 - phi1

    Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
    f0 <- prev / Z
    f1 <- R1 * f0
    f2 <- R2 * f0

    if (f1 > 1 || f2 > 1) {
      warning("Computed penetrances f1 or f2 exceed 1; check prev/R1/R2 inputs.")
    }

    C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0

    D   <- delta_prime * pd * p_plus
    DpT <- D * p_plus
    DpA <- D * pd

    denom <- phi1 + pi01 * phi0

    # Eq. 5.28a,b
    gT_star  <- (pd * p_plus) + (DpT * C * (1 - pi01)) / denom
    gNT_star <- (pd * p_plus) + (DpA * C * (pi01 - 1)) / denom

    list(
      gT_star = gT_star,
      gNT_star = gNT_star,
      p_plus = p_plus,
      phi1 = phi1,
      phi0 = phi0,
      C = C,
      D = D,
      f0 = f0,
      f1 = f1,
      f2 = f2
    )
  }

  ## ---------- internal helper: gT* and gNT* with heterogeneity ----------
  calc_gTgNT_heter <- function(pd, prev, R1, R2,
                               delta_prime, heter_rate) {
    p_plus <- 1 - pd
    Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
    f0 <- prev / Z
    f1 <- R1 * f0
    f2 <- R2 * f0

    if (f1 > 1 || f2 > 1) {
      warning("Computed penetrances f1 or f2 exceed 1; check prev/R1/R2 inputs.")
    }

    C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
    delta <- delta_prime * pd * p_plus

    # pi = fraction of trios truly due to this locus
    pi <- 1 - heter_rate

    gT_star  <- pi * ( p_plus * ( delta * C / prev + pd ) ) +
      (1 - pi) * ( delta * (p_plus - 0.5) * C / prev + pd * p_plus )

    gNT_star <- pi * ( pd * ( p_plus - delta * C / prev ) ) +
      (1 - pi) * ( delta * (0.5 - pd) * C / prev + pd * p_plus )

    list(
      gT_star = gT_star,
      gNT_star = gNT_star,
      p_plus = p_plus,
      C = C,
      f0 = f0,
      f1 = f1,
      f2 = f2,
      delta = delta,
      pi = pi
    )
  }

  lambda_from_gTgNT <- function(N, gT_star, gNT_star) {
    2 * N * (gT_star - gNT_star)^2 / (gT_star + gNT_star)
  }

  crit <- qchisq(1 - alpha, df = 1)

  ## ===== (1) NO MISCLASSIFICATION, NO HETEROGENEITY =====
  p_plus <- 1 - pd
  Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
  f0 <- prev / Z
  f1 <- R1 * f0
  f2 <- R2 * f0

  if (f1 > 1 || f2 > 1) {
    warning("Computed penetrances f1 or f2 exceed 1; check prev/R1/R2 inputs.")
  }

  C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  prev_ <- pd^2 * f2 + 2 * pd * p_plus * f1 + p_plus^2 * f0
  delta <- delta_prime * pd * p_plus

  ET_nomisc  <- 2 * N * (pd * p_plus + delta * p_plus * C / prev_)
  ENT_nomisc <- 2 * N * (pd * p_plus - delta * pd    * C / prev_)

  lambda_nomisc <- (ET_nomisc - ENT_nomisc)^2 / (ET_nomisc + ENT_nomisc)
  power_nomisc  <- 1 - pchisq(crit, df = 1, ncp = lambda_nomisc)

  # define gT*, gNT* for reporting
  gT_nomisc  <- ET_nomisc  / (2 * N)
  gNT_nomisc <- ENT_nomisc / (2 * N)

  ## ===== (2) MISCLASSIFICATION ONLY =====
  g_misc <- calc_gTgNT_misclass(
    pd = pd, prev = prev, R1 = R1, R2 = R2,
    delta_prime = delta_prime,
    pi01 = misclass_rate
  )
  gT_misc  <- g_misc$gT_star
  gNT_misc <- g_misc$gNT_star

  lambda_misc <- lambda_from_gTgNT(N, gT_misc, gNT_misc)
  power_misc  <- 1 - pchisq(crit, df = 1, ncp = lambda_misc)
  ET_misc     <- 2 * N * gT_misc
  ENT_misc    <- 2 * N * gNT_misc

  ## ===== (3) HETEROGENEITY ONLY =====
  g_het <- calc_gTgNT_heter(
    pd = pd, prev = prev, R1 = R1, R2 = R2,
    delta_prime = delta_prime,
    heter_rate = heter_rate
  )
  gT_het  <- g_het$gT_star
  gNT_het <- g_het$gNT_star

  lambda_het <- lambda_from_gTgNT(N, gT_het, gNT_het)
  power_het  <- 1 - pchisq(crit, df = 1, ncp = lambda_het)
  ET_het     <- 2 * N * gT_het
  ENT_het    <- 2 * N * gNT_het

  ## ----- Losses -----
  power_loss_misc <- power_nomisc - power_misc
  power_loss_het  <- power_nomisc - power_het

  # ----- Printed summary -----
  if (isTRUE(verbose)) {
    message("\nPOWER FOR A FIXED SAMPLE SIZE\n")

    # nice two-column header with aligned right-hand values
    fmt_two_col <- "%-32s %10s  |  %-28s %10s"
    message(sprintf(
      fmt_two_col,
      "Number of Trios (N):",
      formatC(N, format = "d"),
      "Significance Level (alpha):",
      formatC(alpha, format = "f", digits = 3)
    ))
    message(sprintf(
      fmt_two_col,
      "Allele Frequency (p_d):",
      formatC(pd, format = "f", digits = 3),
      "Prevalence (phi1):",
      formatC(prev, format = "f", digits = 3)
    ))
    message(sprintf("%-32s %10s", "Relative Risks (R1,R2):",
                    paste0(R1, ", ", R2)))
    message(sprintf("%-32s %10.3f", "LD scale (delta_prime):", delta_prime))
    message(sprintf("%-32s %10.3f", "Misclassification Rate (pi01):", misclass_rate))
    message(sprintf("%-32s %10.3f\n", "Heterogeneity Rate (1 - pi):", heter_rate))

    message("Non-Centrality Parameters (lambda)")
    message("---------------------------------------------------------------")
    message(sprintf("  %-18s %10.4f", "No error:",        lambda_nomisc))
    message(sprintf("  %-18s %10.4f", "Misclassification:", lambda_misc))
    message(sprintf("  %-18s %10.4f\n", "Heterogeneity:",    lambda_het))

    message("Power")
    message(sprintf("  %-18s %10.3f", "No error:",        power_nomisc))
    message(sprintf("  %-18s %10.3f  (loss = %6.3f)",
                    "Misclassification:", power_misc, power_loss_misc))
    message(sprintf("  %-18s %10.3f  (loss = %6.3f)\n",
                    "Heterogeneity:",    power_het,  power_loss_het))

    message("Expected transmission/non-transmission probabilities (gT*, gNT*)")
    fmt_g <- "  %-18s %10.5f  |  %-18s %10.5f"
    message(sprintf(fmt_g,
                    "gT* (no error):", gT_nomisc,
                    "gNT* (no error):", gNT_nomisc))
    message(sprintf(fmt_g,
                    "gT* (misclass):", gT_misc,
                    "gNT* (misclass):", gNT_misc))
    message(sprintf(fmt_g,
                    "gT* (heter):", gT_het,
                    "gNT* (heter):", gNT_het))
  }

  # ----- Clean return object -----
  out <- list(
    alpha = alpha,
    N = N,
    lambda = list(
      no_error = lambda_nomisc,
      misclassification = lambda_misc,
      heterogeneity = lambda_het
    ),
    power = list(
      no_error = power_nomisc,
      misclassification = power_misc,
      heterogeneity = power_het
    ),
    power_loss = list(
      misclassification = power_loss_misc,
      heterogeneity = power_loss_het
    ),
    gT_star = list(
      no_error = gT_nomisc,
      misclassification = gT_misc,
      heterogeneity = gT_het
    ),
    gNT_star = list(
      no_error = gNT_nomisc,
      misclassification = gNT_misc,
      heterogeneity = gNT_het
    ),
    ET = list(
      no_error = ET_nomisc,
      misclassification = ET_misc,
      heterogeneity = ET_het
    ),
    ENT = list(
      no_error = ENT_nomisc,
      misclassification = ENT_misc,
      heterogeneity = ENT_het
    ),
    model_parameters = list(
      pd = pd,
      qd = 1 - pd,
      prev = prev,
      R1 = R1,
      R2 = R2,
      delta_prime = delta_prime,
      misclass_rate = misclass_rate,
      heter_rate = heter_rate
    )
  )

  class(out) <- "tdt_power_full"
  invisible(out)
}


#' Required Number of Trios for TDT with Error Components
#'
#' Computes the minimum number of affected trios required to achieve a
#' specified power for the Transmission Disequilibrium Test (TDT) under
#' three scenarios:
#' (i) no error, (ii) phenotype misclassification only, and
#' (iii) locus heterogeneity only.
#'
#' @param target_power Numeric in (0,1). Desired power for the TDT.
#' @param pd Numeric in (0,1). Frequency of the disease (high-risk) allele
#'   at the marker locus.
#' @param prev Numeric in (0,1). Disease prevalence (\eqn{\phi_1}).
#' @param R1 Numeric \eqn{> 0}. Genotype relative risk for heterozygotes.
#' @param R2 Numeric \eqn{> 0}. Genotype relative risk for homozygotes.
#' @param alpha Numeric in (0,1). Significance level for the TDT
#'   (default \code{0.05}).
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter
#'   \eqn{D'} (default \code{1}).
#' @param misclass_rate Numeric in \eqn{[0,1)}. Phenotype misclassification
#'   rate for controls (\eqn{\pi_{01}}) in the misclassification scenario.
#' @param heter_rate Numeric in \eqn{[0,1)}. Proportion of trios whose
#'   affection status is not due to the locus of interest (\eqn{1 - \pi})
#'   in the heterogeneity scenario.
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted
#'   summary of the required numbers of trios, percent inflation, and
#'   the resulting power when the no-error design is used.
#'
#' @details
#' The function first solves for the non-centrality parameter
#' \eqn{\lambda^*} that yields \code{target_power} for a non-central
#' chi-square distribution with one degree of freedom and significance
#' level \code{alpha}.  For each scenario, the expected transmission and
#' non-transmission probabilities \eqn{g_T^*} and \eqn{g_{NT}^*} are
#' computed and the required number of trios is:
#' \deqn{
#' N^* = \frac{\lambda^* (g_T^* + g_{NT}^*)}
#'            {2 (g_T^* - g_{NT}^*)^2}.
#' }
#'
#' Penetrances are derived from \code{prev}, \code{R1}, and \code{R2} as in
#' \code{\link{tdt_power_full}}.  For coherence, the function stops with an
#' error if any of \eqn{f_0, f_1, f_2} fall outside the interval \eqn{[0,1]}.
#'
#' Percent increases are reported as ratios relative to the no-error design,
#' minus 1; e.g., a value of \code{0.25} corresponds to a 25\% increase in
#' the required number of trios.
#'
#' In addition, the function evaluates the power achieved under each error
#' scenario when using the no-error sample size \code{N(no\_error)}, and
#' reports the corresponding power loss.
#'
#' @return
#' An object of class \code{"tdt_required_trios_full"}: a list with components
#' \describe{
#'   \item{alpha}{Significance level used.}
#'   \item{target_power}{Requested power.}
#'   \item{lambda_star}{Non-centrality parameter \eqn{\lambda^*} solving
#'     for \code{target_power} under no error.}
#'   \item{N}{List with elements \code{no_error}, \code{misclassification},
#'     and \code{heterogeneity} giving the required number of trios in each
#'     scenario.}
#'   \item{percent_increase}{List with elements \code{misclassification} and
#'     \code{heterogeneity}. Each value is
#'     \code{N_scenario / N_no_error - 1}.}
#'   \item{power_at_N_no_error}{List with elements \code{no_error},
#'     \code{misclassification}, and \code{heterogeneity} giving the power
#'     attained when the design uses \code{N(no\_error)} trios.}
#'   \item{power_loss_at_N_no_error}{List with elements
#'     \code{misclassification} and \code{heterogeneity} giving the loss in
#'     power relative to the no-error case at the same sample size.}
#'   \item{gT_star}{List of expected transmission probabilities
#'     (\eqn{g_T^*}) under each scenario.}
#'   \item{gNT_star}{List of expected non-transmission probabilities
#'     (\eqn{g_{NT}^*}) under each scenario.}
#'   \item{model_parameters}{List of input model parameters
#'     (\code{pd}, \code{prev}, \code{R1}, \code{R2}, \code{delta_prime},
#'     \code{misclass_rate}, \code{heter_rate}).}
#' }
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' \emph{Heterogeneity in Statistical Genetics}. Springer Nature.
#'
#' @examples
#' # Required trios for 80% power under a modest effect
#' tdt_required_trios_full(
#'   target_power  = 0.80,
#'   pd            = 0.30,
#'   prev          = 0.05,
#'   R1            = 1.5,
#'   R2            = 2.25,
#'   alpha         = 0.05,
#'   delta_prime   = 1,
#'   misclass_rate = 0.01,
#'   heter_rate    = 0.10,
#'   verbose       = TRUE
#' )
#'
#' # Required trios without printing
#' tdt_required_trios_full(
#'   target_power = 0.80,
#'   pd = 0.3, prev = 0.05,
#'   R1 = 1.5, R2 = 2.25,
#'   verbose = FALSE
#' )$N$no_error
#'
#' @export
tdt_required_trios_full <- function(
    target_power,
    pd,
    prev,
    R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    misclass_rate = 0.01,
    heter_rate   = 0.01,
    verbose = TRUE
) {
  ## ---- basic argument checks ----
  if (misclass_rate < 0 || misclass_rate >= 1)
    stop("misclass_rate must be in [0, 1).")
  if (heter_rate < 0 || heter_rate >= 1)
    stop("heter_rate must be in [0, 1).")

  ## ---------- internal helper: gT* and gNT* with misclassification (Eq. 5.28a,b) ----------
  calc_gTgNT_misclass <- function(pd, prev, R1, R2,
                                  delta_prime, pi01) {
    p_plus <- 1 - pd
    phi1   <- prev
    phi0   <- 1 - phi1

    Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
    f0 <- prev / Z
    f1 <- R1 * f0
    f2 <- R2 * f0
    if (any(c(f0, f1, f2) < 0) || any(c(f1, f2) > 1)) {
      stop("Invalid penetrances in misclassification component: f0,f1,f2 must be in [0,1]. Check prev, R1, R2, and pd.")
    }
    C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0

    D   <- delta_prime * pd * p_plus
    DpT <- D * p_plus
    DpA <- D * pd

    denom <- phi1 + pi01 * phi0

    gT_star  <- (pd * p_plus) + (DpT * C * (1 - pi01)) / denom
    gNT_star <- (pd * p_plus) + (DpA * C * (pi01 - 1)) / denom

    list(
      gT_star = gT_star,
      gNT_star = gNT_star,
      p_plus = p_plus,
      phi1 = phi1,
      phi0 = phi0,
      C = C,
      D = D,
      f0 = f0,
      f1 = f1,
      f2 = f2
    )
  }

  ## ---------- internal helper: gT* and gNT* with heterogeneity (Eq. 5.34) ----------
  calc_gTgNT_heter <- function(pd, prev, R1, R2,
                               delta_prime, heter_rate) {
    p_plus <- 1 - pd
    Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
    f0 <- prev / Z
    f1 <- R1 * f0
    f2 <- R2 * f0
    if (any(c(f0, f1, f2) < 0) || any(c(f1, f2) > 1)) {
      stop("Invalid penetrances in heterogeneity component: f0,f1,f2 must be in [0,1]. Check prev, R1, R2, and pd.")
    }
    C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
    delta <- delta_prime * pd * p_plus

    # pi = fraction of trios truly due to this locus
    pi <- 1 - heter_rate

    gT_star  <- pi * ( p_plus * ( delta * C / prev + pd ) ) +
      (1 - pi) * ( delta * (p_plus - 0.5) * C / prev + pd * p_plus )

    gNT_star <- pi * ( pd * ( p_plus - delta * C / prev ) ) +
      (1 - pi) * ( delta * (0.5 - pd) * C / prev + pd * p_plus )

    list(
      gT_star = gT_star,
      gNT_star = gNT_star,
      p_plus = p_plus,
      C = C,
      f0 = f0,
      f1 = f1,
      f2 = f2,
      delta = delta,
      pi = pi
    )
  }

  ## ---------- helper: N from lambda* and gT*, gNT* ----------
  N_from_lambda <- function(lambda_star, gT_star, gNT_star) {
    (lambda_star * (gT_star + gNT_star)) /
      (2 * (gT_star - gNT_star)^2)
  }

  ## ---------- solve for lambda_star from target power ----------
  crit <- qchisq(1 - alpha, df = 1)
  f_lambda <- function(lambda) {
    1 - pchisq(crit, df = 1, ncp = lambda) - target_power
  }
  lambda_star <- uniroot(f_lambda, c(0, 1e4))$root

  ## ===== (1) NO MISCLASSIFICATION, NO HETEROGENEITY =====
  p_plus <- 1 - pd
  Z  <- p_plus^2 + 2 * pd * p_plus * R1 + R2 * pd^2
  f0 <- prev / Z
  f1 <- R1 * f0
  f2 <- R2 * f0
  if (any(c(f0, f1, f2) < 0) || any(c(f1, f2) > 1)) {
    stop("Invalid penetrances in no-error component: f0,f1,f2 must be in [0,1]. Check prev, R1, R2, and pd.")
  }
  C  <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  prev_ <- pd^2 * f2 + 2 * pd * p_plus * f1 + p_plus^2 * f0
  delta <- delta_prime * pd * p_plus

  # gT* and gNT* implied by the no-error ET/ENT formulas in tdt_power_full
  gT_nomisc  <- (pd * p_plus + delta * p_plus * C / prev_)
  gNT_nomisc <- (pd * p_plus - delta * pd    * C / prev_)

  N_nomisc <- N_from_lambda(lambda_star, gT_nomisc, gNT_nomisc)

  ## ===== (2) MISCLASSIFICATION ONLY =====
  g_misc <- calc_gTgNT_misclass(
    pd = pd, prev = prev, R1 = R1, R2 = R2,
    delta_prime = delta_prime,
    pi01 = misclass_rate
  )
  gT_misc  <- g_misc$gT_star
  gNT_misc <- g_misc$gNT_star

  N_misc <- N_from_lambda(lambda_star, gT_misc, gNT_misc)

  ## ===== (3) HETEROGENEITY ONLY =====
  g_het <- calc_gTgNT_heter(
    pd = pd, prev = prev, R1 = R1, R2 = R2,
    delta_prime = delta_prime,
    heter_rate = heter_rate
  )
  gT_het  <- g_het$gT_star
  gNT_het <- g_het$gNT_star

  N_het <- N_from_lambda(lambda_star, gT_het, gNT_het)

  ## ===== Percent increase relative to no-error N =====
  infl_misc <- N_misc / N_nomisc
  infl_het  <- N_het  / N_nomisc
  perc_increase_misc <- infl_misc - 1
  perc_increase_het  <- infl_het  - 1

  ## ===== Power loss if you DON'T inflate N (design at N_nomisc) =====
  lambda_nomisc_fixed <- 2 * N_nomisc * (gT_nomisc - gNT_nomisc)^2 /
    (gT_nomisc + gNT_nomisc)
  power_nomisc_fixed <- 1 - pchisq(crit, df = 1, ncp = lambda_nomisc_fixed)

  lambda_misc_fixed <- 2 * N_nomisc * (gT_misc - gNT_misc)^2 /
    (gT_misc + gNT_misc)
  power_misc_fixed <- 1 - pchisq(crit, df = 1, ncp = lambda_misc_fixed)

  lambda_het_fixed <- 2 * N_nomisc * (gT_het - gNT_het)^2 /
    (gT_het + gNT_het)
  power_het_fixed <- 1 - pchisq(crit, df = 1, ncp = lambda_het_fixed)

  power_loss_misc <- power_nomisc_fixed - power_misc_fixed
  power_loss_het  <- power_nomisc_fixed - power_het_fixed

  if (isTRUE(verbose)) {
    message("\nMINIMUM SAMPLE SIZE NECESSARY FOR A FIXED POWER\n")

    message(sprintf("%-38s %10.3f  |  %s %6.3f",
                    "Desired Power:", target_power,
                    "Significance Level (alpha):", alpha))
    message(sprintf("%-38s %10.3f  |  %s %6.3f",
                    "Allele Frequency (p_d):", pd,
                    "Prevalence (phi1):", prev))
    message(sprintf("%-38s %10s",
                    "Relative Risks (R1,R2):", paste0(R1, ", ", R2)))
    message(sprintf("%-38s %10.3f",
                    "LD scale (delta_prime):", delta_prime))
    message(sprintf("%-38s %10.3f",
                    "Misclassification Rate (pi01):", misclass_rate))
    message(sprintf("%-38s %10.3f\n",
                    "Heterogeneity Rate (1 - pi):", heter_rate))

    message(sprintf("%-38s %10.4f",
                    "Non-Centrality Parameter (lambda_star):", lambda_star))
    message("-----------------------------------------------------------")
    message(sprintf("%-38s %10.0f",
                    "Required Trios (no error):", ceiling(N_nomisc)))
    message(sprintf("%-38s %10.0f  |  Percent Increase: %8.3f",
                    "Required Trios (misclassification):",
                    ceiling(N_misc), perc_increase_misc))
    message(sprintf("%-38s %10.0f  |  Percent Increase: %8.3f\n",
                    "Required Trios (heterogeneity):",
                    ceiling(N_het), perc_increase_het))

    message("Power at N(no-error design):")
    fmt_pow  <- "  %-18s %7.3f"
    fmt_loss <- "  %-18s %7.3f  (loss = %6.3f)"
    message(sprintf(fmt_pow,  "No error:",        power_nomisc_fixed))
    message(sprintf(fmt_loss, "Misclassification:", power_misc_fixed, power_loss_misc))
    message(sprintf(fmt_loss, "Heterogeneity:",     power_het_fixed,  power_loss_het))
    message("")

    message("Expected transmission/non-transmission probabilities (gT*, gNT*)")
    message(sprintf("%-18s %10.5f  |  %s %10.5f",
                    "gT* (no error):", gT_nomisc,
                    "gNT* (no error):", gNT_nomisc))
    message(sprintf("%-18s %10.5f  |  %s %10.5f",
                    "gT* (misclass):", gT_misc,
                    "gNT* (misclass):", gNT_misc))
    message(sprintf("%-18s %10.5f  |  %s %10.5f",
                    "gT* (heter):",    gT_het,
                    "gNT* (heter):",   gNT_het))
  }

  out <- list(
    alpha = alpha,
    target_power = target_power,
    lambda_star = lambda_star,
    N = list(
      no_error = N_nomisc,
      misclassification = N_misc,
      heterogeneity = N_het
    ),
    percent_increase = list(
      misclassification = perc_increase_misc,
      heterogeneity = perc_increase_het
    ),
    power_at_N_no_error = list(
      no_error = power_nomisc_fixed,
      misclassification = power_misc_fixed,
      heterogeneity = power_het_fixed
    ),
    power_loss_at_N_no_error = list(
      misclassification = power_loss_misc,
      heterogeneity = power_loss_het
    ),
    gT_star = list(
      no_error = gT_nomisc,
      misclassification = gT_misc,
      heterogeneity = gT_het
    ),
    gNT_star = list(
      no_error = gNT_nomisc,
      misclassification = gNT_misc,
      heterogeneity = gNT_het
    ),
    model_parameters = list(
      pd = pd,
      qd = 1 - pd,
      prev = prev,
      R1 = R1,
      R2 = R2,
      delta_prime = delta_prime,
      misclass_rate = misclass_rate,
      heter_rate = heter_rate
    )
  )

  class(out) <- "tdt_required_trios_full"
  invisible(out)
}

#' Power sensitivity to misclassification and heterogeneity
#'
#' Create a 1×2 panel plot showing how TDT power changes as a function of
#' phenotype misclassification and locus heterogeneity, for a fixed design.
#' The left panel varies the misclassification rate with heterogeneity held
#' fixed; the right panel varies the heterogeneity rate with misclassification
#' held fixed.
#'
#' This function is a graphical wrapper around [tdt_power_full()], looping
#' over user–supplied grids of error rates and extracting the power under
#' misclassification and heterogeneity.
#'
#' @param N Integer. Number of affected trios in the design.
#' @param pd Numeric in (0, 1). High–risk allele frequency at the marker.
#' @param prev Numeric in (0, 1). Disease prevalence \eqn{\phi_1}.
#' @param R1,R2 Numeric. Genotype relative risks for the heterozygote and
#'   homozygote, respectively.
#' @param alpha Numeric. Type I error rate for the TDT (default `0.05`).
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter
#'   \eqn{\delta'} (default `1`).
#' @param misclass_seq Numeric vector. Sequence of phenotype
#'   misclassification rates (pi01) to use on the x–axis of the
#'   misclassification panel.
#' @param heter_seq Numeric vector. Sequence of heterogeneity rates
#'   (proportion of trios not due to the locus, \eqn{1-\pi}) to use on the
#'   x–axis of the heterogeneity panel.
#' @param heter_fixed Numeric scalar. Heterogeneity rate to hold fixed while
#'   varying `misclass_seq` in the misclassification panel (typically `0`).
#' @param misclass_fixed Numeric scalar. Misclassification rate to hold fixed
#'   while varying `heter_seq` in the heterogeneity panel (typically `0`).
#' @param title Character. Title to use for the combined plot.
#'
#' @details
#' In the left panel, each value of `misclass_seq` is passed to
#' [tdt_power_full()] via the `misclass_rate` argument, with
#' `heter_rate = heter_fixed`. The right panel does the reverse, passing
#' `heter_seq` into `heter_rate` with `misclass_rate = misclass_fixed`.
#'
#' @return A `ggplot`/`patchwork` object showing two panels side–by–side:
#'   power vs. misclassification and power vs. heterogeneity. The function is
#'   mainly called for its side–effect of plotting.
#'
#' @seealso [tdt_power_full()], [tdt_required_trios_full()]
#'
#' @examples
#' \donttest{
#' # Basic example: power sensitivity for a design with 600 trios
#' tdt_plot_power(
#'   N    = 600,
#'   pd   = 0.30,
#'   prev = 0.05,
#'   R1   = 1.5,
#'   R2   = 2.25,
#'   alpha = 0.05,
#'   delta_prime = 1,
#'   misclass_seq  = seq(0, 0.15, by = 0.01),
#'   heter_seq     = seq(0, 0.50, by = 0.05),
#'   heter_fixed   = 0,
#'   misclass_fixed = 0
#' )
#' }
#' @importFrom ggplot2 ggplot geom_line geom_point facet_wrap labs theme_bw aes
#' @export
tdt_plot_power <- function(
    pd, prev, R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    N,
    misclass_seq  = seq(0, 0.15, by = 0.01),
    heter_seq     = seq(0, 0.50, by = 0.05),
    heter_fixed   = 0,      # heterogeneity when varying misclassification
    misclass_fixed = 0,     # misclassification when varying heterogeneity
    title = "TDT power under misclassification and heterogeneity"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.", call. = FALSE)
  }

  power_mis <- vapply(misclass_seq, function(mis) {
    res <- tdt_power_full(
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
    error_type = "Misclassification (pi01)",
    error_rate = misclass_seq,
    power      = power_mis
  )

  power_het <- vapply(heter_seq, function(h) {
    res <- tdt_power_full(
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
    error_type = "Heterogeneity (1 - pi)",
    error_rate = heter_seq,
    power      = power_het
  )

  df_all <- rbind(df_mis, df_het)

  p <- ggplot2::ggplot(
    df_all,
    ggplot2::aes(x = .data$error_rate, y = .data$power)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ error_type, nrow = 1, scales = "free_x") +
    ggplot2::labs(
      x = "Error rate",
      y = "Power",
      title = title
    ) +
    ggplot2::theme_bw()


  print(p)
  invisible(p)
}


#' Required sample size sensitivity to misclassification and heterogeneity
#'
#' Create a 1×2 panel plot showing how the required number of trios changes
#' as a function of phenotype misclassification and locus heterogeneity,
#' for a fixed target power.
#'
#' This function is a graphical wrapper around [tdt_required_trios_full()],
#' looping over user–supplied grids of error rates and extracting the
#' required sample size under misclassification and heterogeneity.
#'
#' @param target_power Numeric in (0, 1). Desired power for the TDT.
#' @param pd Numeric in (0, 1). High–risk allele frequency at the marker.
#' @param prev Numeric in (0, 1). Disease prevalence \eqn{\phi_1}.
#' @param R1,R2 Numeric. Genotype relative risks for the heterozygote and
#'   homozygote, respectively.
#' @param alpha Numeric. Type I error rate for the TDT (default `0.05`).
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter
#'   \eqn{\delta'} (default `1`).
#' @param misclass_seq Numeric vector. Sequence of phenotype
#'   misclassification rates (pi01) to use on the x–axis of the
#'   misclassification panel.
#' @param heter_seq Numeric vector. Sequence of heterogeneity rates
#'   (proportion of trios not due to the locus, \eqn{1-\pi}) to use on the
#'   x–axis of the heterogeneity panel.
#' @param heter_fixed Numeric scalar. Heterogeneity rate to hold fixed while
#'   varying `misclass_seq` in the misclassification panel (typically `0`).
#' @param misclass_fixed Numeric scalar. Misclassification rate to hold fixed
#'   while varying `heter_seq` in the heterogeneity panel (typically `0`).
#' @param title Character. Title to use for the combined plot.
#'
#' @details
#' In the left panel, each value of `misclass_seq` is passed to
#' [tdt_required_trios_full()] via the `misclass_rate` argument, with
#' `heter_rate = heter_fixed`. The right panel does the reverse, passing
#' `heter_seq` into `heter_rate` with `misclass_rate = misclass_fixed`.
#'
#' @return A `ggplot`/`patchwork` object showing two panels side–by–side:
#'   required trios vs. misclassification and required trios vs. heterogeneity.
#'   The function is mainly called for its side–effect of plotting.
#'
#' @seealso [tdt_required_trios_full()], [tdt_power_full()]
#'
#' @examples
#' \donttest{
#' # Sample size sensitivity for power 0.80
#' tdt_plot_sample_size(
#'   target_power = 0.80,
#'   pd   = 0.30,
#'   prev = 0.05,
#'   R1   = 1.5,
#'   R2   = 2.25,
#'   alpha = 0.05,
#'   delta_prime = 1,
#'   misclass_seq  = seq(0, 0.15, by = 0.01),
#'   heter_seq     = seq(0, 0.50, by = 0.05),
#'   heter_fixed   = 0,
#'   misclass_fixed = 0
#' )
#' }
#'
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_line geom_point facet_wrap labs theme_bw aes
#' @export
tdt_plot_sample_size <- function(
    pd, prev, R1, R2,
    alpha = 0.05,
    delta_prime = 1,
    target_power,
    misclass_seq  = seq(0, 0.15, by = 0.01),
    heter_seq     = seq(0, 0.50, by = 0.05),
    heter_fixed   = 0,
    misclass_fixed = 0,
    title = "Required trios under misclassification and heterogeneity"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plot.", call. = FALSE)
  }

  N_mis <- vapply(misclass_seq, function(mis) {
    res <- tdt_required_trios_full(
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
    error_type = "Misclassification (pi01)",
    error_rate = misclass_seq,
    N_required = N_mis
  )

  N_het <- vapply(heter_seq, function(h) {
    res <- tdt_required_trios_full(
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
    error_type = "Heterogeneity (1 - pi)",
    error_rate = heter_seq,
    N_required = N_het
  )

  df_all <- rbind(df_mis, df_het)

  p <- ggplot2::ggplot(
    df_all,
    ggplot2::aes(x = .data$error_rate, y = .data$N_required)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ error_type, nrow = 1, scales = "free_x") +
    ggplot2::labs(
      x = "Error rate",
      y = "Required number of trios",
      title = title
    ) +
    ggplot2::theme_bw()


  print(p)
  invisible(p)
}











