# ==============================================================================
# Canonical TDT framework
#
# The validated model-based formulas and three-scenario reporting (no error /
# misclassification / heterogeneity) are preserved. The public interface also
# supports input_mode = c("model_based", "model_free"):
#
#   - "model_based" (default) uses genetic-model parameters.
#   - "model_free" lets a user who already has expected transmission and
#     non-transmission counts (ET, ENT) supply them directly instead of
#     prev/R1/R2. gT = ET / (2 * n_trios), gNT = ENT / (2 * n_trios), where
#     n_trios is N itself for the power function (ET/ENT are assumed to be
#     for the same N trios power is computed for) and a separate n_trios
#     argument for the MSSN function (ET/ENT are absolute counts that embed
#     a sample size, which need not equal the N* being solved for).
#
# heter_rate and misclass_rate apply in both modes. In model_free mode they
# use closed-form identities (verified algebraically equivalent to the
# calc_gTgNT_heter()/calc_gTgNT_misclass() identities in terms of gT, gNT,
# pd, and prev alone -- see the
# roxygen details below. R1/R2 are never used in model_free mode.
#
# ==============================================================================


#' Family-Based (TDT) Power
#'
#' Computes power for the transmission disequilibrium test (TDT) at a fixed
#' number of affected trios under three scenarios: (i) no error, (ii)
#' phenotype misclassification only, and (iii) locus heterogeneity only.
#' \code{input_mode} lets the transmission probabilities come either from a
#' genetic model (\code{"model_based"}, the default) or directly from
#' user-supplied expected transmission and non-transmission counts
#' (\code{"model_free"}).
#'
#' @param N Numeric \eqn{> 0}. Number of affected trios.
#' @param input_mode Character. One of \code{"model_based"} (default) or
#'   \code{"model_free"}. See Details.
#' @param pd Numeric in (0,1). Frequency of the disease/high-risk allele at
#'   the marker locus. Required for \code{input_mode = "model_based"}.
#'   Optional for \code{"model_free"}: required only if \code{heter_rate} or
#'   \code{misclass_rate} is non-zero, and if omitted in that case it is
#'   solved from \code{ET}/\code{ENT} (see Details).
#' @param prev Numeric in (0,1). Disease prevalence (\eqn{\phi_1}). Required
#'   for \code{input_mode = "model_based"}. Required for \code{"model_free"}
#'   only if \code{misclass_rate} is non-zero.
#' @param R1 Numeric \eqn{> 0}. Genotype relative risk for heterozygotes.
#'   Required for \code{input_mode = "model_based"}; unused for
#'   \code{"model_free"}.
#' @param R2 Numeric \eqn{> 0}. Genotype relative risk for homozygotes.
#'   Required for \code{input_mode = "model_based"}; unused for
#'   \code{"model_free"}.
#' @param alpha Numeric in (0,1). Significance level for the TDT
#'   (default \code{0.05}).
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter
#'   \eqn{D'} (default \code{1}). Used for \code{input_mode = "model_based"}
#'   only.
#' @param misclass_rate Numeric in \eqn{[0,1)}. Phenotype misclassification
#'   rate for controls (\eqn{\pi_{01}}). A value of \code{0} corresponds to
#'   no misclassification and is the default.
#' @param heter_rate Numeric in \eqn{[0,1)}. Proportion of trios whose
#'   affection status is \emph{not} due to the locus of interest
#'   (\eqn{1 - \pi}). A value of \code{0} corresponds to complete
#'   homogeneity and is the default.
#' @param ET,ENT Numeric \eqn{\ge 0}. Expected transmission and
#'   non-transmission counts for \code{input_mode = "model_free"}, assumed to
#'   be accumulated over the same \code{N} trios that power is computed for.
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted
#'   summary of no-error power and any requested modifier-specific power.
#'
#' @details
#' With \code{input_mode = "model_based"}, penetrances
#' \eqn{f_0, f_1, f_2} are derived
#' from \code{prev}, \code{R1}, \code{R2}, and \code{pd}, and the expected
#' transmission and non-transmission probabilities under each of the three
#' scenarios are computed from the standard TDT formulation in Eq. 1.25
#' (Chapter 1, Section 1.6.1.3, p. 27), the phenotype-misclassification
#' probabilities in Eqs. 5.24--5.27 (Section 5.2.6, pp. 284--285), and the
#' locus-heterogeneity construction in Eqs. 5.30--5.34a (Section 5.3.3,
#' pp. 293--294) of Gordon, Finch, and Kim (2020). Eq. 5.28 is a numerical
#' worked example, not the general symbolic formula.
#' Phenotype misclassification and locus heterogeneity are reported as
#' separate sensitivity scenarios; supplying both rates does not create a
#' combined-error scenario.
#'
#' \code{N} is the number of affected-child trios, with both parents genotyped.
#' \code{ET} and \code{ENT} are expected counts accumulated over trios, whereas
#' \eqn{g_T} and \eqn{g_{NT}} are per-parental-allele transmission and
#' non-transmission probabilities. \code{heter_rate} is the heterogeneous
#' fraction; in lower-level formulas written with the linked/homogeneous
#' fraction \eqn{\pi}; equivalently, the heterogeneous fraction is
#' \eqn{1-\pi}.
#'
#' With \code{input_mode = "model_free"}, the no-error scenario uses
#' \eqn{g_T = ET / (2N)} and \eqn{g_{NT} = ENT / (2N)} directly. The
#' misclassification and heterogeneity scenarios are then computed from the
#' following identities, which are algebraically equivalent to the
#' the corresponding model-based identities. Let
#' \eqn{A = g_T - g_{NT}} (no-error) and
#' \eqn{p_+ = 1 - p_d}:
#' \deqn{\text{heterogeneity: } g_T = p_d p_+ + A (p_+ - 0.5 \times
#'   \text{heter\_rate}), \quad
#'   g_{NT} = p_d p_+ + A (-p_d + 0.5 \times \text{heter\_rate}),}
#' \deqn{\text{misclassification: } m = \frac{\phi_1 (1 - \pi_{01})}
#'   {\phi_1 + \pi_{01} (1 - \phi_1)}, \quad
#'   g_T = p_d p_+ + p_+ A m, \quad g_{NT} = p_d p_+ - p_d A m.}
#' A scenario is only computed from these identities when its rate is
#' non-zero; at a rate of exactly \code{0} the scenario reuses the no-error
#' \eqn{g_T}/\eqn{g_{NT}} directly, so \code{pd} is not required unless at
#' least one rate is non-zero, and \code{prev} is not required unless
#' \code{misclass_rate} is non-zero.
#'
#' If \code{pd} is needed but not supplied, it is solved from \code{ET} and
#' \code{ENT} via \eqn{2 p_d^2 - 2 p_d (1 - A) + (g_T + g_{NT} - A) = 0}. This
#' quadratic has two roots summing to \eqn{1 - A}; the root in \eqn{(0, 0.5)}
#' is used, with a message reporting the derived value. An error is raised if
#' no such unique root exists -- supplying \code{pd} directly is preferred.
#'
#' @return An object of class \code{"tdt_power"}, containing:
#' \describe{
#' \item{alpha, N, input_mode}{Significance level, affected-trio count, and
#' input mode.}
#' \item{lambda, power}{Named \code{no_error}, \code{misclassification}, and
#' \code{heterogeneity} NCP and power values.}
#' \item{power_loss}{Absolute power loss for each modifier relative to the
#' no-error scenario.}
#' \item{gT_star, gNT_star}{Scenario-specific transmission and
#' non-transmission probabilities.}
#' \item{ET, ENT}{Scenario-specific expected transmission and
#' non-transmission counts for the supplied \code{N} trios.}
#' \item{model_parameters}{Supplied and derived genetic-model, LD,
#' misclassification, and heterogeneity parameters.}
#' }
#'
#' @examples
#' # Model-based input
#' tdt_power(
#'   N = 600, input_mode = "model_based",
#'   pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
#'   misclass_rate = 0.01, heter_rate = 0.10,
#'   verbose = FALSE
#' )$power$no_error
#'
#' # model_free: supply expected transmissions/non-transmissions directly
#' tdt_power(
#'   N = 600, input_mode = "model_free",
#'   ET = 140, ENT = 100,
#'   pd = 0.30, prev = 0.05,
#'   misclass_rate = 0.01, heter_rate = 0.10,
#'   verbose = FALSE
#' )$power$no_error
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
#' Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is not
#' a case: Effects of phenotype misclassification on power and sample size
#' requirements for the transmission disequilibrium test with affected child
#' trios. \emph{Human Heredity}, 67(4), 287--292.
#' \doi{10.1159/000194981}.
#'
#' Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
#' (2009). Transmission disequilibrium test power and sample size in the
#' presence of locus heterogeneity. \emph{Statistical Applications in Genetics
#' and Molecular Biology}, 8, Article 44. \doi{10.2202/1544-6115.1501}.
#'
#' @seealso \code{\link{tdt_mssn}} for the MSSN counterpart.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export
tdt_power <- function(
    N,
    input_mode = c("model_based", "model_free"),
    pd   = NULL,
    prev = NULL,
    R1   = NULL,
    R2   = NULL,
    alpha = 0.05,
    delta_prime = 1,
    misclass_rate = 0,
    heter_rate   = 0,
    ET  = NULL,
    ENT = NULL,
    verbose = TRUE
) {
  input_mode <- match.arg(input_mode)

  ## ---- basic argument checks ----
  if (misclass_rate < 0 || misclass_rate >= 1)
    stop("misclass_rate must be in [0, 1).")
  if (heter_rate < 0 || heter_rate >= 1)
    stop("heter_rate must be in [0, 1).")

  if (input_mode == "model_based") {
    if (is.null(pd) || is.null(prev) || is.null(R1) || is.null(R2))
      stop("For input_mode='model_based', pd, prev, R1, and R2 must all be supplied.")
  } else {
    if (is.null(ET) || is.null(ENT))
      stop("For input_mode='model_free', ET and ENT must both be supplied.")
    if (misclass_rate != 0 && is.null(prev))
      stop("For input_mode='model_free' with misclass_rate != 0, prev must be supplied.")
  }

  ## ---------- internal helper: gT* and gNT* with misclassification (Eq. 5.26) ----------
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

  ## ---------- internal helper: gT* and gNT* with heterogeneity (Eqs. 5.33-5.34a) ----------
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

  ## ---------- internal helper: solve pd from gT, gNT when not supplied ----------
  solve_pd_from_gTgNT <- function(gT, gNT) {
    A <- gT - gNT
    qa <- 2
    qb <- -2 * (1 - A)
    qc <- gT + gNT - A
    disc <- qb^2 - 4 * qa * qc
    if (!is.finite(disc) || disc < 0)
      stop("Cannot solve for pd from the supplied ET/ENT: no real root exists. ",
           "Supply pd directly.")
    roots <- (-qb + c(1, -1) * sqrt(disc)) / (2 * qa)
    candidates <- roots[roots > 0 & roots < 0.5]
    if (length(candidates) != 1)
      stop("Cannot uniquely determine pd in (0, 0.5) from the supplied ET/ENT; ",
           "supply pd directly.")
    message(
      "pd not supplied for input_mode='model_free'; using pd = ",
      formatC(candidates, format = "f", digits = 6),
      " solved from the supplied ET/ENT (root in (0, 0.5))."
    )
    candidates
  }

  lambda_from_gTgNT <- function(N, gT_star, gNT_star) {
    2 * N * (gT_star - gNT_star)^2 / (gT_star + gNT_star)
  }

  crit <- qchisq(1 - alpha, df = 1)

  if (input_mode == "model_based") {

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

  } else {

    ## ===== model-free: (1) NO MISCLASSIFICATION, NO HETEROGENEITY =====
    ET_nomisc  <- ET
    ENT_nomisc <- ENT
    gT_nomisc  <- ET  / (2 * N)
    gNT_nomisc <- ENT / (2 * N)

    lambda_nomisc <- (ET_nomisc - ENT_nomisc)^2 / (ET_nomisc + ENT_nomisc)
    power_nomisc  <- 1 - pchisq(crit, df = 1, ncp = lambda_nomisc)

    A <- gT_nomisc - gNT_nomisc

    if (is.null(pd) && (heter_rate != 0 || misclass_rate != 0)) {
      pd <- solve_pd_from_gTgNT(gT_nomisc, gNT_nomisc)
    }

    ## ===== (2) MISCLASSIFICATION ONLY =====
    if (misclass_rate == 0) {
      gT_misc  <- gT_nomisc
      gNT_misc <- gNT_nomisc
    } else {
      p_plus <- 1 - pd
      m <- prev * (1 - misclass_rate) / (prev + misclass_rate * (1 - prev))
      gT_misc  <- pd * p_plus + p_plus * A * m
      gNT_misc <- pd * p_plus - pd     * A * m
    }

    lambda_misc <- lambda_from_gTgNT(N, gT_misc, gNT_misc)
    power_misc  <- 1 - pchisq(crit, df = 1, ncp = lambda_misc)
    ET_misc     <- 2 * N * gT_misc
    ENT_misc    <- 2 * N * gNT_misc

    ## ===== (3) HETEROGENEITY ONLY =====
    if (heter_rate == 0) {
      gT_het  <- gT_nomisc
      gNT_het <- gNT_nomisc
    } else {
      p_plus <- 1 - pd
      gT_het  <- pd * p_plus + A * (p_plus - 0.5 * heter_rate)
      gNT_het <- pd * p_plus + A * (-pd    + 0.5 * heter_rate)
    }

    lambda_het <- lambda_from_gTgNT(N, gT_het, gNT_het)
    power_het  <- 1 - pchisq(crit, df = 1, ncp = lambda_het)
    ET_het     <- 2 * N * gT_het
    ENT_het    <- 2 * N * gNT_het
  }

  ## ----- Losses -----
  power_loss_misc <- power_nomisc - power_misc
  power_loss_het  <- power_nomisc - power_het

  # ----- Printed summary -----
  if (isTRUE(verbose)) {
    message("TDT power")
    message(sprintf("Affected-child trios: %s", formatC(N, format = "d")))
    message("")
    message("No-error design")
    message(sprintf("  Power: %.1f%%", 100 * power_nomisc))

    if (misclass_rate > 0) {
      message("")
      message("Phenotype misclassification")
      message(sprintf("  Misclassification rate: %.1f%%", 100 * misclass_rate))
      message(sprintf("  Power: %.1f%%", 100 * power_misc))
      message(sprintf("  Absolute power loss: %.1f percentage points",
                      100 * power_loss_misc))
    }

    if (heter_rate > 0) {
      message("")
      message("Locus heterogeneity")
      message(sprintf("  Heterogeneity rate: %.1f%%", 100 * heter_rate))
      message(sprintf("  Power: %.1f%%", 100 * power_het))
      message(sprintf("  Absolute power loss: %.1f percentage points",
                      100 * power_loss_het))
    }
  }

  # ----- Clean return object -----
  out <- list(
    alpha = alpha,
    N = N,
    input_mode = input_mode,
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
      qd = if (!is.null(pd)) 1 - pd else NA_real_,
      prev = prev,
      R1 = R1,
      R2 = R2,
      delta_prime = delta_prime,
      misclass_rate = misclass_rate,
      heter_rate = heter_rate
    )
  )

  class(out) <- "tdt_power"
  invisible(out)
}


#' Family-Based (TDT) Minimum Sample Size Necessary
#'
#' Computes the minimum number of affected trios required to achieve a
#' specified power for the transmission disequilibrium test (TDT) under three
#' scenarios: (i) no error, (ii) phenotype misclassification only, and (iii)
#' locus heterogeneity only. \code{input_mode} lets the transmission
#' probabilities come either from a genetic model
#' (\code{"model_based"}, the default) or directly from user-supplied
#' expected transmission and non-transmission counts (\code{"model_free"}).
#'
#' @param target_power Numeric in (0,1). Desired power for the TDT.
#' @param input_mode Character. One of \code{"model_based"} (default) or
#'   \code{"model_free"}. See Details.
#' @param pd Numeric in (0,1). Frequency of the disease (high-risk) allele at
#'   the marker locus. Required for \code{input_mode = "model_based"}.
#'   Optional for \code{"model_free"}: required only if \code{heter_rate} or
#'   \code{misclass_rate} is non-zero, and if omitted in that case it is
#'   solved from \code{ET}/\code{ENT} (see Details).
#' @param prev Numeric in (0,1). Disease prevalence (\eqn{\phi_1}). Required
#'   for \code{input_mode = "model_based"}. Required for \code{"model_free"}
#'   only if \code{misclass_rate} is non-zero.
#' @param R1 Numeric \eqn{> 0}. Genotype relative risk for heterozygotes.
#'   Required for \code{input_mode = "model_based"}; unused for
#'   \code{"model_free"}.
#' @param R2 Numeric \eqn{> 0}. Genotype relative risk for homozygotes.
#'   Required for \code{input_mode = "model_based"}; unused for
#'   \code{"model_free"}.
#' @param alpha Numeric in (0,1). Significance level for the TDT
#'   (default \code{0.05}).
#' @param delta_prime Numeric. Linkage disequilibrium scale parameter
#'   \eqn{D'} (default \code{1}). Used for \code{input_mode = "model_based"}
#'   only.
#' @param misclass_rate Numeric in \eqn{[0,1)}. Phenotype misclassification
#'   rate for controls (\eqn{\pi_{01}}) in the misclassification scenario.
#'   Zero means no misclassification and is the default.
#' @param heter_rate Numeric in \eqn{[0,1)}. Proportion of trios whose
#'   affection status is not due to the locus of interest (\eqn{1 - \pi})
#'   in the heterogeneity scenario. Zero means complete homogeneity and is the
#'   default.
#' @param ET,ENT Numeric \eqn{\ge 0}. Expected transmission and
#'   non-transmission counts for \code{input_mode = "model_free"}.
#' @param n_trios Numeric \eqn{> 0}. Number of affected trios that the
#'   supplied \code{ET} and \code{ENT} correspond to. Required for
#'   \code{input_mode = "model_free"} (unlike the power function, there is no
#'   \code{N} argument here to reuse, since this function solves for the
#'   required number of trios).
#' @param verbose Logical. If \code{TRUE} (default), prints a formatted
#'   summary of the no-error required number of trios and any requested
#'   modifier-specific required counts and percentage increases.
#'
#' @details
#' With \code{input_mode = "model_based"}, penetrances are derived from
#' \code{prev}, \code{R1}, \code{R2}, and \code{pd}, and
#' \deqn{N^* = \frac{\lambda^* (g_T^* + g_{NT}^*)}{2 (g_T^* - g_{NT}^*)^2}}
#' is computed for each of the three scenarios.
#' The no-error calculation follows Eq. 1.25 (p. 27); phenotype
#' misclassification follows Eqs. 5.24--5.27 (pp. 284--285); and locus
#' heterogeneity follows Eqs. 5.30--5.34b (pp. 293--294). Eq. 5.28 is a
#' numerical example. Here \code{heter_rate} is the heterogeneous fraction,
#' equal to \eqn{1-\pi} when lower-level formulas use homogeneous fraction
#' \eqn{\pi}.
#' Phenotype misclassification and locus heterogeneity are reported as
#' separate sensitivity scenarios; supplying both rates does not create a
#' combined-error scenario.
#'
#' With \code{input_mode = "model_free"}, the no-error scenario uses
#' \eqn{g_T = ET / (2\,n_{trios})} and \eqn{g_{NT} = ENT / (2\,n_{trios})}.
#' The misclassification and heterogeneity scenarios then use the same
#' closed-form identities as \code{\link{tdt_power}} (see its
#' Details for the formulas and the \code{pd}-solving fallback), applied to
#' this no-error \eqn{g_T}/\eqn{g_{NT}} pair.
#'
#' @return An object of class \code{"tdt_mssn"}, containing:
#' \describe{
#' \item{alpha, target_power, lambda_star, input_mode}{Design targets, target
#' one-df NCP, and input mode.}
#' \item{N}{Named \code{no_error}, \code{misclassification}, and
#' \code{heterogeneity} required affected-trio counts.}
#' \item{percent_increase}{Modifier-specific conventional percentage inflation,
#' \eqn{100(N_{adjusted}/N_{no-error}-1)}, relative to the no-error required
#' count. Values are \code{NA} when the no-error required count is not finite
#' and positive.}
#' \item{power_at_N_no_error, power_loss_at_N_no_error}{Power and power loss
#' obtained if the no-error design size is used under each modifier.}
#' \item{gT_star, gNT_star}{Scenario-specific transmission and
#' non-transmission probabilities used in the MSSN formulas.}
#' \item{model_parameters}{Supplied and derived genetic-model, LD,
#' misclassification, heterogeneity, and model-free count-scale information.}
#' }
#'
#' @examples
#' # Model-based input
#' tdt_mssn(
#'   target_power = 0.80, input_mode = "model_based",
#'   pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
#'   verbose = FALSE
#' )$N$no_error
#'
#' # model_free: supply expected transmissions/non-transmissions directly
#' tdt_mssn(
#'   target_power = 0.80, input_mode = "model_free",
#'   ET = 140, ENT = 100, n_trios = 120,
#'   pd = 0.30, prev = 0.05,
#'   verbose = FALSE
#' )$N$no_error
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
#' Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is not
#' a case: Effects of phenotype misclassification on power and sample size
#' requirements for the transmission disequilibrium test with affected child
#' trios. \emph{Human Heredity}, 67(4), 287--292.
#' \doi{10.1159/000194981}.
#'
#' Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
#' (2009). Transmission disequilibrium test power and sample size in the
#' presence of locus heterogeneity. \emph{Statistical Applications in Genetics
#' and Molecular Biology}, 8, Article 44. \doi{10.2202/1544-6115.1501}.
#'
#' @seealso \code{\link{tdt_power}} for the power counterpart.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export
tdt_mssn <- function(
    target_power,
    input_mode = c("model_based", "model_free"),
    pd   = NULL,
    prev = NULL,
    R1   = NULL,
    R2   = NULL,
    alpha = 0.05,
    delta_prime = 1,
    misclass_rate = 0,
    heter_rate   = 0,
    ET  = NULL,
    ENT = NULL,
    n_trios = NULL,
    verbose = TRUE
) {
  input_mode <- match.arg(input_mode)

  ## ---- basic argument checks ----
  if (misclass_rate < 0 || misclass_rate >= 1)
    stop("misclass_rate must be in [0, 1).")
  if (heter_rate < 0 || heter_rate >= 1)
    stop("heter_rate must be in [0, 1).")

  if (input_mode == "model_based") {
    if (is.null(pd) || is.null(prev) || is.null(R1) || is.null(R2))
      stop("For input_mode='model_based', pd, prev, R1, and R2 must all be supplied.")
  } else {
    if (is.null(ET) || is.null(ENT))
      stop("For input_mode='model_free', ET and ENT must both be supplied.")
    if (is.null(n_trios))
      stop("For input_mode='model_free', n_trios must be supplied: it is the ",
           "number of affected trios that the given ET and ENT correspond to.")
    if (misclass_rate != 0 && is.null(prev))
      stop("For input_mode='model_free' with misclass_rate != 0, prev must be supplied.")
  }

  ## ---------- internal helper: gT* and gNT* with misclassification (Eq. 5.26) ----------
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

  ## ---------- internal helper: gT* and gNT* with heterogeneity (Eqs. 5.33-5.34a) ----------
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

  ## ---------- internal helper: solve pd from gT, gNT when not supplied ----------
  solve_pd_from_gTgNT <- function(gT, gNT) {
    A <- gT - gNT
    qa <- 2
    qb <- -2 * (1 - A)
    qc <- gT + gNT - A
    disc <- qb^2 - 4 * qa * qc
    if (!is.finite(disc) || disc < 0)
      stop("Cannot solve for pd from the supplied ET/ENT: no real root exists. ",
           "Supply pd directly.")
    roots <- (-qb + c(1, -1) * sqrt(disc)) / (2 * qa)
    candidates <- roots[roots > 0 & roots < 0.5]
    if (length(candidates) != 1)
      stop("Cannot uniquely determine pd in (0, 0.5) from the supplied ET/ENT; ",
           "supply pd directly.")
    message(
      "pd not supplied for input_mode='model_free'; using pd = ",
      formatC(candidates, format = "f", digits = 6),
      " solved from the supplied ET/ENT (root in (0, 0.5))."
    )
    candidates
  }

  ## ---------- helper: N from lambda* and gT*, gNT* ----------
  N_from_lambda <- function(lambda_star, gT_star, gNT_star) {
    (lambda_star * (gT_star + gNT_star)) /
      (2 * (gT_star - gNT_star)^2)
  }

  percent_increase_from_baseline <- function(adjusted, baseline) {
    if (!is.finite(baseline) || baseline <= 0)
      return(NA_real_)
    100 * (adjusted / baseline - 1)
  }

  ## ---------- solve for lambda_star from target power ----------
  crit <- qchisq(1 - alpha, df = 1)
  f_lambda <- function(lambda) {
    1 - pchisq(crit, df = 1, ncp = lambda) - target_power
  }
  lambda_star <- uniroot(f_lambda, c(0, 1e4))$root

  if (input_mode == "model_based") {

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

    # gT* and gNT* implied by the no-error ET/ENT formulas in tdt_power
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

  } else {

    ## ===== model-free: (1) NO MISCLASSIFICATION, NO HETEROGENEITY =====
    gT_nomisc  <- ET  / (2 * n_trios)
    gNT_nomisc <- ENT / (2 * n_trios)

    N_nomisc <- N_from_lambda(lambda_star, gT_nomisc, gNT_nomisc)

    A <- gT_nomisc - gNT_nomisc

    if (is.null(pd) && (heter_rate != 0 || misclass_rate != 0)) {
      pd <- solve_pd_from_gTgNT(gT_nomisc, gNT_nomisc)
    }

    ## ===== (2) MISCLASSIFICATION ONLY =====
    if (misclass_rate == 0) {
      gT_misc  <- gT_nomisc
      gNT_misc <- gNT_nomisc
    } else {
      p_plus <- 1 - pd
      m <- prev * (1 - misclass_rate) / (prev + misclass_rate * (1 - prev))
      gT_misc  <- pd * p_plus + p_plus * A * m
      gNT_misc <- pd * p_plus - pd     * A * m
    }

    N_misc <- N_from_lambda(lambda_star, gT_misc, gNT_misc)

    ## ===== (3) HETEROGENEITY ONLY =====
    if (heter_rate == 0) {
      gT_het  <- gT_nomisc
      gNT_het <- gNT_nomisc
    } else {
      p_plus <- 1 - pd
      gT_het  <- pd * p_plus + A * (p_plus - 0.5 * heter_rate)
      gNT_het <- pd * p_plus + A * (-pd    + 0.5 * heter_rate)
    }

    N_het <- N_from_lambda(lambda_star, gT_het, gNT_het)
  }

  ## ===== Percent increase relative to no-error N =====
  perc_increase_misc <- percent_increase_from_baseline(N_misc, N_nomisc)
  perc_increase_het <- percent_increase_from_baseline(N_het, N_nomisc)

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
    format_trios <- function(x) {
      if (is.finite(x)) {
        formatC(ceiling(x), format = "d", big.mark = ",")
      } else {
        as.character(x)
      }
    }
    format_percent_increase <- function(x) {
      if (is.na(x)) "not defined" else sprintf("%.1f%%", x)
    }

    message("TDT minimum sample size")
    message(sprintf("Target power: %.1f%%", 100 * target_power))
    message("")
    message("No-error design")
    message(sprintf("  Required trios: %s", format_trios(N_nomisc)))

    if (misclass_rate > 0) {
      message("")
      message("Phenotype misclassification")
      message(sprintf("  Misclassification rate: %.1f%%", 100 * misclass_rate))
      message(sprintf("  Required trios: %s", format_trios(N_misc)))
      message(sprintf("  Percent increase: %s",
                      format_percent_increase(perc_increase_misc)))
    }

    if (heter_rate > 0) {
      message("")
      message("Locus heterogeneity")
      message(sprintf("  Heterogeneity rate: %.1f%%", 100 * heter_rate))
      message(sprintf("  Required trios: %s", format_trios(N_het)))
      message(sprintf("  Percent increase: %s",
                      format_percent_increase(perc_increase_het)))
    }
  }

  out <- list(
    alpha = alpha,
    target_power = target_power,
    input_mode = input_mode,
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
      qd = if (!is.null(pd)) 1 - pd else NA_real_,
      prev = prev,
      R1 = R1,
      R2 = R2,
      delta_prime = delta_prime,
      misclass_rate = misclass_rate,
      heter_rate = heter_rate
    )
  )

  class(out) <- "tdt_mssn"
  invisible(out)
}
