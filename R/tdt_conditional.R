# ==============================================================================
# Generalized conditional TDT framework
#
# This file is the family-based analogue of the generalized case-control
# framework in R/case_control.R. It provides a single computational backend for
# TDT power and minimum-sample-size-necessary (MSSN) calculations:
#
#   tdt_power_conditional_full()   -- power for a fixed number of affected trios
#   tdt_mssn_conditional_full()    -- MSSN (required trios) for a target power
#
# Both functions accept input_mode = "model_based" or input_mode = "model_free",
# mirroring cc_power_conditional_full() and cc_mssn_conditional_full().
#
# Scope of this version: the standard TDT only. The transmission pipeline in
# .tdt_transmission_pipeline() is staged so that locus heterogeneity, phenotype
# misclassification, and genotype misclassification can be added later without
# changing the surrounding architecture or the returned object layout.
# ==============================================================================


# ---- internal helpers --------------------------------------------------------

#' Target non-centrality parameter for a chi-square test
#'
#' @param power Numeric in (0,1). Target power.
#' @param alpha Numeric in (0,1). Significance level.
#' @param df Integer. Degrees of freedom.
#'
#' @return Numeric non-centrality parameter achieving \code{power} at
#'   \code{alpha}.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @noRd
.tdt_chisq_ncp_target <- function(power, alpha, df = 1) {
  crit <- qchisq(1 - alpha, df = df)
  f <- function(lambda) {
    pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - power
  }
  uniroot(f, lower = 0, upper = 1e6)$root
}


#' Validate a single number lying in an interval
#'
#' @param x Value to check.
#' @param name Character. Argument name used in the error message.
#' @param lower,upper Numeric bounds.
#' @param strict Logical. If \code{TRUE}, the bounds are exclusive.
#'
#' @return Invisibly \code{TRUE}; called for its side effect of stopping.
#'
#' @noRd
.tdt_check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                              strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x))
    stop(name, " must be a single finite number.")
  if (strict) {
    if (x <= lower || x >= upper)
      stop(name, " must be a single number in (", lower, ",", upper, ").")
  } else {
    if (x < lower || x > upper)
      stop(name, " must be a single number in [", lower, ",", upper, "].")
  }
  invisible(TRUE)
}


#' Resolve the heterozygote relative risk from a mode of inheritance
#'
#' Mirrors the mode-of-inheritance handling in
#' \code{cc_power_conditional_full()} and \code{cc_mssn_conditional_full()}:
#' multiplicative gives \eqn{R_1 = \sqrt{R_2}}, dominant gives \eqn{R_1 = R_2},
#' and recessive gives \eqn{R_1 = 1}. An explicitly supplied \code{R1}
#' overrides the mode of inheritance.
#'
#' @param R1 Numeric or \code{NULL}. Heterozygote relative risk.
#' @param R2 Numeric. Homozygote relative risk.
#' @param MOI Character. One of \code{"M"}, \code{"D"}, or \code{"Rec"}.
#'
#' @return A list with \code{R1}, \code{R2}, \code{MOI}, and \code{R1_source}.
#'
#' @noRd
.tdt_resolve_relative_risks <- function(R1, R2, MOI) {
  if (is.null(R2))
    stop("For input_mode='model_based', R2 must be supplied.")

  .tdt_check_scalar(R2, "R2", lower = 0, strict = TRUE)

  if (is.null(R1)) {
    R1 <- if (MOI == "M") sqrt(R2) else if (MOI == "D") R2 else 1
    R1_source <- paste0("derived from MOI='", MOI, "'")
  } else {
    .tdt_check_scalar(R1, "R1", lower = 0, strict = TRUE)
    R1_source <- "supplied explicitly (MOI ignored)"
  }

  list(R1 = R1, R2 = R2, MOI = MOI, R1_source = R1_source)
}


#' Expected transmission and non-transmission probabilities from a genetic model
#'
#' Implements the model-based components of Equation 1.25 of Gordon et al.
#' (2020). Penetrances follow Equations 1.6 and 1.7:
#' \deqn{f_0 = \phi_1 / (p_+^2 + 2 R_1 p_+ p_d + R_2 p_d^2),\quad
#'       f_1 = R_1 f_0,\quad f_2 = R_2 f_0.}
#' The transmission probabilities per heterozygous parent are
#' \deqn{g_T^* = p_d p_+ + \delta p_+ C / \phi_1,\quad
#'       g_{NT}^* = p_d p_+ - \delta p_d C / \phi_1,}
#' with \eqn{C = p_d f_2 + (1 - 2 p_d) f_1 - p_+ f_0} and
#' \eqn{\delta = \delta' p_d p_+}. Expected counts over \eqn{N} affected trios
#' are \eqn{ET = 2 N g_T^*} and \eqn{ENT = 2 N g_{NT}^*}.
#'
#' @param prev Numeric in (0,1). Disease prevalence \eqn{\phi_1}.
#' @param pd Numeric in (0,1). Disease-allele frequency.
#' @param R1,R2 Numeric > 0. Heterozygote and homozygote relative risks.
#' @param delta_prime Numeric in \eqn{[0,1]}. Scaled linkage disequilibrium
#'   parameter.
#'
#' @return A list with \code{gT}, \code{gNT}, penetrances, \code{C},
#'   \code{delta}, and \code{phi1_check}.
#'
#' @noRd
.tdt_model_transmission_probs <- function(prev, pd, R1, R2, delta_prime) {
  p_plus <- 1 - pd

  # Eq. 1.6 / 1.7: penetrances from prevalence and genotype relative risks.
  Z  <- p_plus^2 + 2 * R1 * pd * p_plus + R2 * pd^2
  f0 <- prev / Z
  f1 <- R1 * f0
  f2 <- R2 * f0

  if (f1 > 1 || f2 > 1)
    warning("Computed penetrances f1 or f2 exceed 1; check prev, R1, and R2.")

  # Eq. 1.25 components.
  C     <- pd * f2 + (1 - 2 * pd) * f1 - p_plus * f0
  delta <- delta_prime * pd * p_plus

  # phi1 recomputed from the penetrances; algebraically identical to prev and
  # retained as an internal consistency check.
  phi1_check <- pd^2 * f2 + 2 * pd * p_plus * f1 + p_plus^2 * f0

  gT  <- pd * p_plus + delta * p_plus * C / phi1_check
  gNT <- pd * p_plus - delta * pd     * C / phi1_check

  list(
    gT = gT,
    gNT = gNT,
    penetrances = c(f0 = f0, f1 = f1, f2 = f2),
    C = C,
    delta = delta,
    p_plus = p_plus,
    phi1_check = phi1_check
  )
}


#' Transmission-probability pipeline for the conditional TDT framework
#'
#' Produces the expected transmission and non-transmission probabilities per
#' heterozygous parent that determine the TDT non-centrality parameter. The
#' pipeline is the family-based analogue of the case-control genotype-frequency
#' pipeline: a baseline stage followed by optional heterogeneity modifiers.
#'
#' Stage 1 (implemented) builds baseline \eqn{g_T^*} and \eqn{g_{NT}^*} either
#' from a genetic model or from user-supplied \eqn{ET} and \eqn{ENT}.
#'
#' Stages 2 to 4 are reserved for locus heterogeneity, phenotype
#' misclassification, and genotype misclassification. They currently pass the
#' baseline values through unchanged and record disabled modifier slots, so that
#' the returned object layout is stable across future versions.
#'
#' @param input_mode Character. \code{"model_based"} or \code{"model_free"}.
#' @param prev,pd,R1,R2,delta_prime Model-based inputs. See
#'   \code{.tdt_model_transmission_probs()}.
#' @param ET,ENT Numeric. Model-free expected transmission and non-transmission
#'   counts.
#' @param n_trios Numeric. Number of affected trios the model-free \code{ET} and
#'   \code{ENT} correspond to.
#'
#' @return A list with baseline and final \code{gT}/\code{gNT}, a
#'   \code{model_info} block, a \code{locus_het} block, and an \code{errors}
#'   block.
#'
#' @noRd
.tdt_transmission_pipeline <- function(input_mode,
                                       prev, pd, R1, R2, MOI, delta_prime,
                                       ET, ENT, n_trios) {

  # ---- stage 1: baseline transmission probabilities ----
  if (input_mode == "model_based") {

    .tdt_check_scalar(prev, "prev", lower = 0, upper = 1, strict = TRUE)
    .tdt_check_scalar(pd, "pd", lower = 0, upper = 1, strict = TRUE)
    .tdt_check_scalar(delta_prime, "delta_prime", lower = 0, upper = 1)

    rr <- .tdt_resolve_relative_risks(R1 = R1, R2 = R2, MOI = MOI)

    mod <- .tdt_model_transmission_probs(
      prev = prev, pd = pd,
      R1 = rr$R1, R2 = rr$R2,
      delta_prime = delta_prime
    )

    gT_base  <- mod$gT
    gNT_base <- mod$gNT

    model_info <- list(
      input_mode = "model_based",
      prev = prev,
      pd = pd,
      qd = mod$p_plus,
      R1 = rr$R1,
      R2 = rr$R2,
      MOI = rr$MOI,
      R1_source = rr$R1_source,
      delta_prime = delta_prime,
      delta = mod$delta,
      C = mod$C,
      penetrances = mod$penetrances,
      phi1_check = mod$phi1_check
    )

  } else {

    if (is.null(ET) || is.null(ENT))
      stop("For input_mode='model_free', ET and ENT must both be supplied.")

    .tdt_check_scalar(ET, "ET", lower = 0)
    .tdt_check_scalar(ENT, "ENT", lower = 0)

    if (ET + ENT <= 0)
      stop("ET + ENT must be positive.")

    # ET and ENT are counts accumulated over n_trios affected trios, each
    # contributing two parental transmissions: ET = 2 * n_trios * gT.
    gT_base  <- ET  / (2 * n_trios)
    gNT_base <- ENT / (2 * n_trios)

    model_info <- list(
      input_mode = "model_free",
      ET_supplied = ET,
      ENT_supplied = ENT,
      n_trios_supplied = n_trios
    )
  }

  # ---- stage 2: locus heterogeneity (not implemented in this version) ----
  # TODO: apply the locus-heterogeneity mixture to gT/gNT, mirroring the
  # case-control modifier g_case = pi * g_case + (1 - pi) * g_ctrl. See
  # Gordon et al. (2020) Sect. 5.3.3.
  gT  <- gT_base
  gNT <- gNT_base

  locus_het_info <- list(
    enabled = FALSE,
    pi = 1,
    gT_before_locus_het = gT_base,
    gNT_before_locus_het = gNT_base,
    gT_after_locus_het = gT,
    gNT_after_locus_het = gNT
  )

  # ---- stage 3: phenotype misclassification (not implemented) ----
  # TODO: Gordon et al. (2020) Sect. 5.2.6.
  pheno_misclass_info <- list(
    enabled = FALSE,
    model = "none",
    gT_after_pheno_misclass = gT,
    gNT_after_pheno_misclass = gNT
  )

  # ---- stage 4: genotype misclassification (not implemented) ----
  # TODO: Gordon et al. (2020) Sect. 5.2.5, including type I error inflation.
  geno_misclass_info <- list(
    enabled = FALSE,
    model = "none"
  )

  list(
    gT_base = gT_base,
    gNT_base = gNT_base,
    gT = gT,
    gNT = gNT,
    model_info = model_info,
    locus_het = locus_het_info,
    errors = list(
      phenotype_misclass = pheno_misclass_info,
      genotype_misclass = geno_misclass_info
    )
  )
}


#' Validate transmission probabilities before computing the NCP
#'
#' @param gT,gNT Numeric transmission and non-transmission probabilities.
#' @param context Character. \code{"power"} or \code{"mssn"}; controls whether a
#'   null effect is a warning or an error.
#'
#' @return Invisibly \code{TRUE}.
#'
#' @noRd
.tdt_check_transmission_probs <- function(gT, gNT, context = c("power", "mssn")) {
  context <- match.arg(context)

  if (!is.finite(gT) || !is.finite(gNT))
    stop("Transmission probabilities are not finite; check inputs.")
  if (gT < 0 || gNT < 0)
    stop("Transmission probabilities cannot be negative; check pd, delta_prime, ",
         "and the relative risks.")
  if (gT + gNT <= 0)
    stop("gT + gNT must be positive; check inputs.")

  if (abs(gT - gNT) < 1e-12) {
    if (context == "mssn") {
      stop("gT equals gNT, so the non-centrality parameter is 0 and no finite ",
           "sample size attains the target power. This happens when there is ",
           "no linkage disequilibrium (delta_prime = 0) or no genotype effect ",
           "(R1 = R2 = 1).")
    } else {
      warning("gT equals gNT, so the non-centrality parameter is 0 and power ",
              "equals alpha.")
    }
  }

  invisible(TRUE)
}


#' Format helpers for clean console output
#' @noRd
.tdt_fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

#' @noRd
.tdt_fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)


#' Print the shared input block for the conditional TDT summaries
#' @noRd
.tdt_print_inputs <- function(pipe, input_mode, alpha) {
  fmt2 <- "%-32s %12s  |  %-28s %12s"
  mi <- pipe$model_info

  message(sprintf(
    fmt2,
    "Input Mode:", input_mode,
    "Significance Level (alpha):", .tdt_fmt_e(alpha, 2)
  ))

  if (input_mode == "model_based") {
    message(sprintf(
      fmt2,
      "Disease Prevalence (phi1):", .tdt_fmt_f(mi$prev, 4),
      "Risk Allele Freq (p_d):", .tdt_fmt_f(mi$pd, 4)
    ))
    message(sprintf(
      fmt2,
      "MOI:", mi$MOI,
      "Relative Risks (R1,R2):",
      paste0(.tdt_fmt_f(mi$R1, 3), ",", .tdt_fmt_f(mi$R2, 3))
    ))
    message(sprintf(
      "%-32s %12s",
      "LD parameter (delta_prime):", .tdt_fmt_f(mi$delta_prime, 4)
    ))
  } else {
    message("Model-free input: user-supplied ET and ENT")
    message(sprintf(
      fmt2,
      "Expected Transmissions (ET):", .tdt_fmt_f(mi$ET_supplied, 4),
      "Expected Non-Transm. (ENT):", .tdt_fmt_f(mi$ENT_supplied, 4)
    ))
    message(sprintf(
      "%-32s %12s",
      "Trios for supplied ET/ENT:", format(mi$n_trios_supplied)
    ))
  }

  message(sprintf("%-32s %12s", "Locus heterogeneity:", "none"))
  message(sprintf("%-32s %12s", "Phenotype misclassification:", "none"))
  message(sprintf("%-32s %12s", "Genotype misclassification:", "none"))

  invisible(NULL)
}


# ---- exported functions ------------------------------------------------------

#' Family-Based (TDT) Power for Conditional Transmission Probabilities
#'
#' Computes power for the transmission disequilibrium test at a fixed number of
#' affected trios. This is the family-based analogue of
#' \code{\link{cc_power_conditional_full}}: it accepts either a genetic model or
#' user-supplied expected transmission counts through a single generalized
#' interface, and returns a structured object with the same overall layout.
#'
#' @param n_trios Numeric \eqn{> 0}. Number of affected trios.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param input_mode Character. One of \code{"model_based"} or
#'   \code{"model_free"}. See Details.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence \eqn{\phi_1} for
#'   \code{input_mode = "model_based"}.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency for
#'   \code{input_mode = "model_based"}.
#' @param R1 Numeric \eqn{> 0} or \code{NULL}. Heterozygote relative risk. When
#'   \code{NULL}, it is derived from \code{MOI}.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk for
#'   \code{input_mode = "model_based"}.
#' @param MOI Character. Mode of inheritance used to derive \code{R1} when
#'   \code{R1} is \code{NULL}: \code{"M"} for multiplicative, \code{"D"} for
#'   dominant, or \code{"Rec"} for recessive.
#' @param delta_prime Numeric in \eqn{[0,1]}. Scaled linkage disequilibrium
#'   parameter \eqn{\delta'} between the marker and trait loci.
#' @param ET,ENT Numeric \eqn{\ge 0} for \code{input_mode = "model_free"}.
#'   Expected transmission and non-transmission counts accumulated over
#'   \code{n_trios} affected trios.
#' @param verbose Logical. If \code{TRUE}, prints a clean formatted summary.
#'
#' @details
#' The workflow is:
#' \enumerate{
#' \item Construct baseline expected transmission probabilities \eqn{g_T^*} and
#' \eqn{g_{NT}^*} per heterozygous parent.
#' \item Reserved for locus heterogeneity (not implemented in this version).
#' \item Reserved for phenotype misclassification (not implemented).
#' \item Reserved for genotype misclassification (not implemented).
#' \item Compute the non-centrality parameter and power from the resulting
#' expected counts.
#' }
#'
#' With \code{input_mode = "model_based"}, penetrances are derived from
#' \code{prev}, \code{pd}, \code{R1}, and \code{R2} using Equations 1.6 and 1.7,
#' and expected counts follow Equation 1.25:
#' \deqn{ET = 2 N [p_d p_+ + \delta p_+ C / \phi_1], \quad
#'       ENT = 2 N [p_d p_+ - \delta p_d C / \phi_1],}
#' where \eqn{C = p_d f_2 + (1 - 2 p_d) f_1 - p_+ f_0} and
#' \eqn{\delta = \delta' p_d p_+}.
#'
#' With \code{input_mode = "model_free"}, the user supplies \code{ET} and
#' \code{ENT} directly. Unlike the case-control framework, whose model-free
#' inputs are case and control genotype frequencies, the natural model-free
#' inputs here are the expected transmission counts, because they completely
#' determine the TDT non-centrality parameter. No genotype frequencies are
#' involved.
#'
#' In both modes the non-centrality parameter is
#' \deqn{\lambda = (ET - ENT)^2 / (ET + ENT),}
#' and the TDT statistic follows a central chi-square distribution with 1 degree
#' of freedom under the null and a non-central chi-square distribution with 1
#' degree of freedom under the alternative (Table 1.5).
#'
#' The \code{locus_het} and \code{errors} components of the returned object are
#' always present and currently report disabled modifiers, so that downstream
#' wrappers and plotting functions can rely on a stable layout as heterogeneity
#' models are added.
#'
#' @return An object of class \code{"tdt_power_conditional_full"}: a nested list
#' with components \code{alpha}, \code{n_trios}, \code{input_mode},
#' \code{delta_prime}, \code{locus_het}, \code{errors}, \code{model_info},
#' \code{tests}, and \code{transmissions}. \code{tests$tdt} contains the test
#' label, degrees of freedom, non-centrality parameter, internal \code{S}
#' component, and \code{power}. \code{transmissions} stores baseline and
#' observed \code{gT}/\code{gNT} probabilities and \code{ET}/\code{ENT} counts.
#'
#' @examples
#' tdt_power_conditional_full(
#'   n_trios = 10000, alpha = 0.05,
#'   input_mode = "model_based",
#'   prev = 0.05, pd = 0.25, R1 = 1, R2 = 1.1,
#'   delta_prime = 1,
#'   verbose = FALSE
#' )
#'
#' tdt_power_conditional_full(
#'   n_trios = 500, alpha = 0.05,
#'   input_mode = "model_based",
#'   prev = 0.005, pd = 0.25, R2 = 2, MOI = "D",
#'   delta_prime = 1,
#'   verbose = FALSE
#' )
#'
#' tdt_power_conditional_full(
#'   n_trios = 120, alpha = 0.05,
#'   input_mode = "model_free",
#'   ET = 140, ENT = 100,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#' Equations 1.6, 1.7, and 1.25; Table 1.5.
#'
#' Spielman, R. S., McGinnis, R. E., & Ewens, W. J. (1993). Transmission test
#' for linkage disequilibrium: the insulin gene region and insulin-dependent
#' diabetes mellitus (IDDM). \emph{American Journal of Human Genetics}, 52,
#' 506-516.
#'
#' @seealso \code{\link{tdt_mssn_conditional_full}} for the sample-size
#'   counterpart, and \code{\link{cc_power_conditional_full}} for the
#'   population-based analogue.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export
tdt_power_conditional_full <- function(
    n_trios, alpha,

    input_mode = c("model_based", "model_free"),

    # model-based inputs
    prev = NULL,
    pd   = NULL,
    R1   = NULL,
    R2   = NULL,
    MOI  = c("M", "D", "Rec"),
    delta_prime = 1,

    # model-free inputs
    ET  = NULL,
    ENT = NULL,

    verbose = TRUE
) {

  input_mode <- match.arg(input_mode)
  MOI <- match.arg(MOI)

  # ---- checks ----
  .tdt_check_scalar(n_trios, "n_trios", lower = 0, strict = TRUE)
  .tdt_check_scalar(alpha, "alpha", lower = 0, upper = 1, strict = TRUE)

  # ---- transmission pipeline ----
  pipe <- .tdt_transmission_pipeline(
    input_mode = input_mode,
    prev = prev, pd = pd, R1 = R1, R2 = R2, MOI = MOI,
    delta_prime = delta_prime,
    ET = ET, ENT = ENT, n_trios = n_trios
  )

  gT  <- pipe$gT
  gNT <- pipe$gNT

  .tdt_check_transmission_probs(gT, gNT, context = "power")

  # ---- expected counts and non-centrality parameter ----
  ET_obs  <- 2 * n_trios * gT
  ENT_obs <- 2 * n_trios * gNT

  # S is the per-trio effect size; lambda = 2 * n_trios * S. Retained for
  # validation and for the MSSN counterpart, but not printed.
  S <- (gT - gNT)^2 / (gT + gNT)
  lambda <- 2 * n_trios * S

  crit  <- qchisq(1 - alpha, df = 1)
  power <- pchisq(crit, df = 1, ncp = lambda, lower.tail = FALSE)

  # ---- output ----
  out <- list(
    alpha = alpha,
    n_trios = n_trios,
    input_mode = input_mode,
    delta_prime = if (input_mode == "model_based") delta_prime else NA_real_,
    locus_het = pipe$locus_het,
    errors = pipe$errors,
    model_info = pipe$model_info,
    tests = list(
      tdt = list(
        test = "transmission disequilibrium test",
        df = 1,
        lambda = lambda,
        S = S,
        power = power
      )
    ),
    transmissions = list(
      gT_base = pipe$gT_base,
      gNT_base = pipe$gNT_base,
      gT_obs = gT,
      gNT_obs = gNT,
      ET = ET_obs,
      ENT = ENT_obs
    )
  )

  class(out) <- "tdt_power_conditional_full"

  # ---- clean printed output ----
  if (isTRUE(verbose)) {

    message("\n--- Family-Based (TDT): Power for Fixed Number of Trios ---")
    message("Equation: 1.25  |  Test: transmission disequilibrium test (df = 1)")
    message("--------------------------------------------------------------------------")

    message(sprintf(
      "%-32s %12s  |  %-28s %12s",
      "Number of Trios (N):", format(n_trios),
      "Test Degrees of Freedom:", "1"
    ))

    .tdt_print_inputs(pipe, input_mode, alpha)

    message("--------------------------------------------------------------------------")
    message("Power")
    message(sprintf(
      "  %-16s lambda=%10.4f  |  power=%8.4f",
      "TDT:", lambda, power
    ))

    message("--------------------------------------------------------------------------")
    message("Expected transmissions and non-transmissions")
    message(sprintf("  %-28s %12.4f", "ET:", ET_obs))
    message(sprintf("  %-28s %12.4f", "ENT:", ENT_obs))
    message(sprintf("  %-28s %12.6f", "gT  (per parent):", gT))
    message(sprintf("  %-28s %12.6f", "gNT (per parent):", gNT))

    message("--------------------------------------------------------------------------")
  }

  invisible(out)
}


#' Family-Based (TDT) Minimum Sample Size for Conditional Transmission Probabilities
#'
#' Computes the minimum sample size necessary (MSSN), expressed as a number of
#' affected trios, for the transmission disequilibrium test to attain a target
#' power. This is the family-based analogue of
#' \code{\link{cc_mssn_conditional_full}} and shares its interface, workflow,
#' and returned-object layout.
#'
#' @param power Numeric in \eqn{(0,1)}. Desired target power.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param input_mode Character. One of \code{"model_based"} or
#'   \code{"model_free"}. See Details.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence \eqn{\phi_1} for
#'   \code{input_mode = "model_based"}.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency for
#'   \code{input_mode = "model_based"}.
#' @param R1 Numeric \eqn{> 0} or \code{NULL}. Heterozygote relative risk. When
#'   \code{NULL}, it is derived from \code{MOI}.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk for
#'   \code{input_mode = "model_based"}.
#' @param MOI Character. Mode of inheritance used to derive \code{R1} when
#'   \code{R1} is \code{NULL}: \code{"M"} for multiplicative, \code{"D"} for
#'   dominant, or \code{"Rec"} for recessive.
#' @param delta_prime Numeric in \eqn{[0,1]}. Scaled linkage disequilibrium
#'   parameter \eqn{\delta'} between the marker and trait loci.
#' @param ET,ENT Numeric \eqn{\ge 0} for \code{input_mode = "model_free"}.
#'   Expected transmission and non-transmission counts accumulated over
#'   \code{n_trios} affected trios.
#' @param n_trios Numeric \eqn{> 0}. Number of affected trios that the supplied
#'   \code{ET} and \code{ENT} correspond to. Required for
#'   \code{input_mode = "model_free"} and ignored otherwise. See Details.
#' @param verbose Logical. If \code{TRUE}, prints a clean formatted summary.
#'
#' @details
#' The workflow mirrors \code{\link{tdt_power_conditional_full}}: baseline
#' transmission probabilities are constructed, reserved heterogeneity stages
#' pass them through unchanged in this version, and the MSSN is computed from
#' the resulting per-trio effect size.
#'
#' Because \eqn{ET = 2 N g_T^*} and \eqn{ENT = 2 N g_{NT}^*} are both
#' proportional to the number of trios, the non-centrality parameter is linear
#' in \eqn{N}:
#' \deqn{\lambda = 2 N (g_T^* - g_{NT}^*)^2 / (g_T^* + g_{NT}^*).}
#' The MSSN therefore has a closed form and requires no simulation:
#' \deqn{N^* = \left\lceil \frac{\lambda^*}{2}
#'   \frac{g_T^* + g_{NT}^*}{(g_T^* - g_{NT}^*)^2} \right\rceil,}
#' where \eqn{\lambda^*} is the non-centrality parameter attaining the target
#' power at the given significance level with 1 degree of freedom.
#'
#' With \code{input_mode = "model_free"}, \code{ET} and \code{ENT} are absolute
#' counts that already embed a sample size, so \code{n_trios} must also be
#' supplied to recover the per-parent probabilities
#' \eqn{g_T^* = ET / (2 n_{trios})} and \eqn{g_{NT}^* = ENT / (2 n_{trios})}
#' before rescaling to \eqn{N^*}. This requirement has no case-control
#' counterpart, because model-free genotype frequencies are already normalized.
#'
#' @return An object of class \code{"tdt_mssn_conditional_full"}: a nested list
#' with components \code{alpha}, \code{target_power}, \code{input_mode},
#' \code{delta_prime}, \code{locus_het}, \code{errors}, \code{model_info},
#' \code{tests}, and \code{transmissions}. \code{tests$tdt} contains the test
#' label, degrees of freedom, target non-centrality parameter
#' \code{lambda_star}, internal \code{S} component, and \code{MSSN_trios}.
#' \code{transmissions} stores baseline and observed \code{gT}/\code{gNT}
#' probabilities and the \code{ET}/\code{ENT} counts implied at the MSSN.
#'
#' @examples
#' tdt_mssn_conditional_full(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_based",
#'   prev = 0.005, pd = 0.25, R1 = 2, R2 = 2,
#'   delta_prime = 1,
#'   verbose = FALSE
#' )
#'
#' tdt_mssn_conditional_full(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_based",
#'   prev = 0.05, pd = 0.25, R2 = 1.5, MOI = "M",
#'   delta_prime = 0.8,
#'   verbose = FALSE
#' )
#'
#' tdt_mssn_conditional_full(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   ET = 140, ENT = 100, n_trios = 120,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#' Equations 1.6, 1.7, and 1.25; Table 1.5.
#'
#' Spielman, R. S., McGinnis, R. E., & Ewens, W. J. (1993). Transmission test
#' for linkage disequilibrium: the insulin gene region and insulin-dependent
#' diabetes mellitus (IDDM). \emph{American Journal of Human Genetics}, 52,
#' 506-516.
#'
#' @seealso \code{\link{tdt_power_conditional_full}} for the power counterpart,
#'   and \code{\link{cc_mssn_conditional_full}} for the population-based
#'   analogue.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export
tdt_mssn_conditional_full <- function(
    power, alpha,

    input_mode = c("model_based", "model_free"),

    # model-based inputs
    prev = NULL,
    pd   = NULL,
    R1   = NULL,
    R2   = NULL,
    MOI  = c("M", "D", "Rec"),
    delta_prime = 1,

    # model-free inputs
    ET  = NULL,
    ENT = NULL,
    n_trios = NULL,

    verbose = TRUE
) {

  input_mode <- match.arg(input_mode)
  MOI <- match.arg(MOI)

  # ---- checks ----
  .tdt_check_scalar(power, "power", lower = 0, upper = 1, strict = TRUE)
  .tdt_check_scalar(alpha, "alpha", lower = 0, upper = 1, strict = TRUE)

  if (input_mode == "model_free") {
    if (is.null(n_trios))
      stop("For input_mode='model_free', n_trios must be supplied: it is the ",
           "number of affected trios that the given ET and ENT correspond to.")
    .tdt_check_scalar(n_trios, "n_trios", lower = 0, strict = TRUE)
  } else {
    # Any value works here; the pipeline normalizes by it only in model-free
    # mode, and MSSN output does not depend on it.
    n_trios <- 1
  }

  # ---- transmission pipeline ----
  pipe <- .tdt_transmission_pipeline(
    input_mode = input_mode,
    prev = prev, pd = pd, R1 = R1, R2 = R2, MOI = MOI,
    delta_prime = delta_prime,
    ET = ET, ENT = ENT, n_trios = n_trios
  )

  gT  <- pipe$gT
  gNT <- pipe$gNT

  .tdt_check_transmission_probs(gT, gNT, context = "mssn")

  # ---- target non-centrality parameter and MSSN ----
  lambda_star <- .tdt_chisq_ncp_target(power = power, alpha = alpha, df = 1)

  S <- (gT - gNT)^2 / (gT + gNT)

  if (!is.finite(S) || S <= 0)
    stop("TDT S <= 0; check inputs.")

  MSSN_trios <- ceiling(lambda_star / (2 * S))

  ET_at_mssn  <- 2 * MSSN_trios * gT
  ENT_at_mssn <- 2 * MSSN_trios * gNT

  # ---- output ----
  out <- list(
    alpha = alpha,
    target_power = power,
    input_mode = input_mode,
    delta_prime = if (input_mode == "model_based") delta_prime else NA_real_,
    locus_het = pipe$locus_het,
    errors = pipe$errors,
    model_info = pipe$model_info,
    tests = list(
      tdt = list(
        test = "transmission disequilibrium test",
        df = 1,
        lambda_star = lambda_star,
        S = S,
        MSSN_trios = MSSN_trios
      )
    ),
    transmissions = list(
      gT_base = pipe$gT_base,
      gNT_base = pipe$gNT_base,
      gT_obs = gT,
      gNT_obs = gNT,
      ET_at_mssn = ET_at_mssn,
      ENT_at_mssn = ENT_at_mssn
    )
  )

  class(out) <- "tdt_mssn_conditional_full"

  # ---- clean printed output ----
  if (isTRUE(verbose)) {

    message("\n--- Family-Based (TDT): Minimum Sample Size Necessary (MSSN) ---")
    message("Equation: 1.25  |  Test: transmission disequilibrium test (df = 1)")
    message("--------------------------------------------------------------------------")

    message(sprintf(
      "%-32s %12s  |  %-28s %12s",
      "Target Power:", .tdt_fmt_f(power, 3),
      "Test Degrees of Freedom:", "1"
    ))

    .tdt_print_inputs(pipe, input_mode, alpha)

    message("--------------------------------------------------------------------------")
    message("Minimum Sample Size Necessary")
    message(sprintf(
      "  %-16s MSSN_trios=%8d  |  MSSN_parents=%8d",
      "TDT:", MSSN_trios, 2 * MSSN_trios
    ))

    message("--------------------------------------------------------------------------")
    message("Expected transmissions and non-transmissions at the MSSN")
    message(sprintf("  %-28s %12.4f", "ET:", ET_at_mssn))
    message(sprintf("  %-28s %12.4f", "ENT:", ENT_at_mssn))
    message(sprintf("  %-28s %12.6f", "gT  (per parent):", gT))
    message(sprintf("  %-28s %12.6f", "gNT (per parent):", gNT))

    message("--------------------------------------------------------------------------")
  }

  invisible(out)
}
