#' Transmission Disequilibrium Test (TDT) Power from Expected Transmissions and Non-Transmissions
#'
#' Computes the statistical power of the Transmission Disequilibrium Test (TDT)
#' given the expected number of transmissions (ET) and non-transmissions (ENT)
#' under a specified significance level. Implements Equation 1.25 from
#' "Gordon et al. (2020), Heterogeneity in Statistical Genetics."

#' @param ET Numeric. Expected number of transmissions.
#'
#' @param ENT Numeric. Expected number of non-transmissions.
#'
#' @param alpha Numeric. Significance level (default = 0.05).
#'
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
#' "Heterogeneity in Statistical Genetics". Springer Nature.
#'
#' @export

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
  message(sprintf("%-38s %10.3f", "False-Positive Rate (pi01):", pi01))
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


tdt_expected_gNT <- function(pd, prev, R1, R2,
                             delta_prime = 1,
                             pi01 = 0,        # misclassification rate
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














