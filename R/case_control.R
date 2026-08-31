.cc_test_components <- function(g_case, g_ctrl, k, w) {
  S_g <- sum((g_case - g_ctrl)^2 / (g_case + k * g_ctrl))
  numerator_t <- (sum(w * (g_case - g_ctrl)))^2
  denominator_t <- sum(w^2 * (g_case + k * g_ctrl)) -
    (sum(w * (g_case + k * g_ctrl)))^2 / (1 + k)

  list(
    S_g = S_g,
    numerator_t = numerator_t,
    denominator_t = denominator_t,
    S_t = numerator_t / denominator_t
  )
}

.cc_validate_test_components <- function(components) {
  if (!is.finite(components$S_g) || components$S_g <= 0)
    stop("Genotype S <= 0; check inputs.")
  if (!is.finite(components$denominator_t) || components$denominator_t <= 0)
    stop("Trend denominator <= 0; check inputs/weights.")
  if (!is.finite(components$numerator_t) || components$numerator_t <= 0)
    stop("Trend numerator is 0; implies no weighted mean difference.")
  invisible(components)
}

.cc_mssn_test_results <- function(g_case, g_ctrl, k, w,
                                  lambda_star_g, lambda_star_t,
                                  validate = TRUE) {
  components <- .cc_test_components(g_case, g_ctrl, k, w)
  null_design <- isTRUE(all(g_case == g_ctrl))
  if (isTRUE(validate) && null_design) {
    stop("No finite MSSN exists because the trend contrast is zero under this design.")
  }
  if (isTRUE(validate))
    .cc_validate_test_components(components)

  case_g <- ceiling(lambda_star_g / (k * components$S_g))
  ctrl_g <- ceiling(k * case_g)
  case_t <- ceiling(lambda_star_t / (k * components$S_t))
  ctrl_t <- ceiling(k * case_t)

  list(
    genotypes = list(
      lambda_star = lambda_star_g,
      S = components$S_g,
      MSSN_case = case_g,
      MSSN_ctrl = ctrl_g,
      MSSN_total = case_g + ctrl_g
    ),
    trend = list(
      lambda_star = lambda_star_t,
      S = components$S_t,
      numerator = components$numerator_t,
      denominator = components$denominator_t,
      MSSN_case = case_t,
      MSSN_ctrl = ctrl_t,
      MSSN_total = case_t + ctrl_t
    )
  )
}

.cc_power_test_results <- function(g_case, g_ctrl, k, w, N_case, alpha,
                                   validate = TRUE) {
  components <- .cc_test_components(g_case, g_ctrl, k, w)
  null_design <- isTRUE(all(g_case == g_ctrl))
  if (isTRUE(validate) && !null_design)
    .cc_validate_test_components(components)

  lambda_g <- k * N_case * components$S_g
  lambda_t <- k * N_case * components$S_t

  list(
    genotypes = list(
      lambda = lambda_g,
      S = components$S_g,
      power = pchisq(qchisq(1 - alpha, df = 2), df = 2, ncp = lambda_g,
                     lower.tail = FALSE)
    ),
    trend = list(
      lambda = lambda_t,
      S = components$S_t,
      numerator = components$numerator_t,
      denominator = components$denominator_t,
      power = pchisq(qchisq(1 - alpha, df = 1), df = 1, ncp = lambda_t,
                     lower.tail = FALSE)
    )
  )
}

.cc_genotype_error_active <- function(M_case, M_ctrl) {
  identity <- diag(3)
  any(M_case != identity) || any(M_ctrl != identity)
}

#' Case-Control Minimum Sample Size for Conditional Genotype Frequencies
#'
#' Computes the minimum sample size necessary (MSSN) for case-control association
#' tests using conditional genotype frequencies. The function supports
#' model-based genotype frequencies, model-free genotype frequencies, optional
#' locus heterogeneity, optional phenotype misclassification, optional genotype
#' misclassification, genotype chi-square tests, and genotype trend tests.
#'
#' @param power Numeric in \eqn{(0,1)}. Desired target power.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param input_mode Character. One of \code{"model_based"} or
#'   \code{"model_free"}. See Details.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence for
#'   \code{input_mode = "model_based"}.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency for
#'   \code{input_mode = "model_based"}.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk for
#'   \code{input_mode = "model_based"}.
#' @param MOI Character. Mode of inheritance for model-based frequencies:
#'   \code{"M"} for multiplicative, \code{"D"} for dominant, or \code{"Rec"}
#'   for recessive.
#' @param g1,g0 Numeric vectors of length 3 for \code{input_mode = "model_free"}.
#'   \code{g1} gives case genotype frequencies and \code{g0} gives control
#'   genotype frequencies, each ordered as \code{c(g0, g1, g2)} and summing to 1.
#' @param locus_het Logical. If \code{TRUE}, applies locus heterogeneity to the
#'   case genotype frequencies before phenotype and genotype misclassification.
#' @param pi Numeric in \eqn{[0,1]}. Locus-homogeneity fraction used when
#'   \code{locus_het = TRUE}. When \code{locus_het = FALSE}, \code{pi} must
#'   remain at its default value of 1.
#' @param pheno_misclass Logical. If \code{TRUE}, applies phenotype
#'   misclassification before genotype misclassification.
#' @param theta Numeric in \eqn{[0,1)}. Probability that a truly affected
#'   individual is classified as a control.
#' @param phi Numeric in \eqn{[0,1)}. Probability that a truly unaffected
#'   individual is classified as a case.
#' @param k Numeric \eqn{> 0}. Control-to-case sample size ratio
#'   \eqn{N_{ctrl} / N_{case}}.
#' @param w Numeric vector of length 3. Genotype trend-test scores. The three
#'   weights cannot all be equal.
#' @param geno_misclass Character. Genotype misclassification model:
#'   \code{"none"}, \code{"1p"}, \code{"2p"}, \code{"3p"}, or \code{"diff3p"}.
#' @param e Numeric in \eqn{[0,0.5]}. For \code{geno_misclass = "1p"}, the
#'   probability assigned to each adjacent off-diagonal genotype call; the
#'   diagonal probability is \eqn{1-2e}. Thus the \eqn{\epsilon} in textbook
#'   Eq. 2.5 equals \eqn{2e}; see Details.
#' @param e1,e2 Numeric. Two-parameter misclassification rates for
#'   \code{geno_misclass = "2p"}.
#' @param e01,e02,e03 Numeric. Non-differential three-parameter
#'   misclassification rates for \code{geno_misclass = "3p"}.
#' @param case_e01,case_e02,case_e03 Numeric. Case-specific three-parameter
#'   misclassification rates for \code{geno_misclass = "diff3p"}.
#' @param ctrl_e01,ctrl_e02,ctrl_e03 Numeric. Control-specific three-parameter
#'   misclassification rates for \code{geno_misclass = "diff3p"}.
#' @param diff_source Character. For \code{geno_misclass = "diff3p"}, one of
#'   \code{"explicit"}, \code{"case"}, or \code{"ctrl"}. See Details.
#' @param diff_multiplier Numeric \eqn{\ge 0}. Multiplier used when
#'   \code{diff_source = "case"} or \code{diff_source = "ctrl"}.
#' @param verbose Logical. If \code{TRUE}, prints a clean formatted summary.
#'
#' @details
#' The workflow is:
#' \enumerate{
#' \item Construct baseline conditional genotype frequencies for cases and
#' controls.
#' \item Optionally apply locus heterogeneity to cases as
#' \eqn{g_{case,true} = \pi g_{case,base} + (1 - \pi) g_{ctrl,base}}.
#' \item Optionally apply phenotype misclassification, where
#' \code{theta = Pr(affected -> control)} and
#' \code{phi = Pr(unaffected -> case)}.
#' \item Optionally apply genotype misclassification matrices to the resulting
#' case and control genotype frequencies.
#' \item Compute genotype chi-square and genotype trend-test MSSN values from
#' the observed genotype frequencies.
#' }
#' The defaults (\code{locus_het = FALSE}, \code{pheno_misclass = FALSE}, and
#' \code{geno_misclass = "none"}) give the ordinary no-error design. When
#' modifiers are active, they are applied sequentially in the order above and
#' therefore form one combined adjusted design.
#'
#' With \code{input_mode = "model_based"}, conditional case and control genotype
#' frequencies are derived from \code{prev}, \code{pd}, \code{R2}, and
#' \code{MOI} using Chapter 1, Section 1.4.2, Eqs. 1.6--1.7 (p. 13) of
#' Gordon, Finch, and Kim (2020). With \code{input_mode = "model_free"}, the user supplies
#' \code{g1} and \code{g0} directly; when phenotype misclassification is
#' enabled, \code{g1} is treated as the true affected genotype distribution and
#' \code{g0} is treated as the true unaffected genotype distribution.
#'
#' Phenotype misclassification requires \code{prev} in both input modes because
#' disease prevalence is used to mix the true affected and unaffected genotype
#' distributions into observed case and control genotype distributions. It is
#' applied after optional locus heterogeneity and before optional genotype
#' misclassification.
#'
#' The genotype misclassification models are:
#' \code{"none"} for identity matrices, \code{"1p"} for one symmetric error
#' rate, \code{"2p"} for adjacent homozygote/heterozygote and heterozygote
#' error rates, \code{"3p"} for a non-differential three-parameter matrix, and
#' \code{"diff3p"} for separate case and control three-parameter matrices.
#' The one-, two-, and three-parameter matrices correspond to textbook
#' Eqs. 2.5, 2.6, and 2.7 (pp. 57--58). For \code{"1p"}, package \code{e}
#' is the probability in each adjacent off-diagonal cell and the diagonal is
#' \eqn{1-2e}; textbook Eq. 2.5 uses off-diagonal \eqn{\epsilon/2}, so
#' \eqn{\epsilon=2e}.
#'
#' The genotype test uses Eq. 1.22 (p. 26). The trend test is defined by
#' Eqs. 1.20--1.21 (p. 24) and uses the NCP in Eq. 1.24 (p. 27). Locus
#' heterogeneity follows Eq. 2.16 (p. 88): \code{pi} is the homogeneous
#' fraction and \eqn{1-\pi} is the heterogeneous fraction. Its trend-test
#' construction corresponds to Eqs. 5.29a--b (pp. 287--288).
#'
#' For \code{geno_misclass = "diff3p"}, \code{diff_source = "explicit"} uses
#' the case and control error parameters exactly as supplied. With
#' \code{diff_source = "case"}, control parameters are computed by multiplying
#' the case parameters by \code{diff_multiplier}. With
#' \code{diff_source = "ctrl"}, case parameters are computed by multiplying the
#' control parameters by \code{diff_multiplier}.
#'
#' With \code{geno_misclass = "diff3p"}, returned MSSNs are nominal,
#' asymptotic values evaluated with the usual chi-square critical value.
#' Different case and control error mechanisms can distort the null
#' distribution and inflate type I error. This function does not independently
#' recalibrate that null distribution for arbitrary differential genotyping
#' error.
#'
#' Internal effect-size components \code{S} are retained in the returned object
#' for validation and debugging, but are not printed in the clean verbose
#' output.
#'
#' @return An object of class \code{"cc_mssn"}, containing:
#' \describe{
#' \item{alpha, target_power}{Requested significance level and power.}
#' \item{input_mode, k, w, locus_het}{Input mode, control-to-case ratio, trend
#' scores, and locus-heterogeneity settings.}
#' \item{errors}{Genotype-error model and matrices, plus phenotype-error
#' settings and intermediate frequencies.}
#' \item{model_info}{Model-based penetrances and risk-model inputs, or
#' model-free identifying information.}
#' \item{tests$genotypes, tests$trend}{Test label, degrees of freedom, target
#' NCP \code{lambda_star}, internal \code{S}, and case, control, and total MSSN.
#' The trend result also contains its numerator and denominator. Sample sizes
#' are numbers of individuals.}
#' \item{freqs}{Baseline, post-heterogeneity (true), post-phenotype-error, and
#' final observed case and control genotype-probability vectors.}
#' }
#'
#' @examples
#' cc_mssn(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   geno_misclass = "none",
#'   verbose = FALSE
#' )
#'
#' cc_mssn(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   verbose = FALSE
#' )
#'
#' cc_mssn(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' cc_mssn(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   prev = 0.05,
#'   pheno_misclass = TRUE, theta = 0.05, phi = 0.01,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' Gordon, D., Finch, S. J., Nothnagel, M., & Ott, J. (2002). Power and sample
#' size calculations for case-control genetic association tests when errors
#' are present: application to single nucleotide polymorphisms.
#' \emph{Human Heredity}, 54(1), 22--33. \doi{10.1159/000066696}.
#'
#' Armitage, P. (1955). Tests for linear trends in proportions and frequencies.
#' \emph{Biometrics}, 11(3), 375--386. \doi{10.2307/3001775}.
#'
#' Slager, S. L., & Schaid, D. J. (2001). Case-control studies of genetic
#' markers: power and sample size approximations for Armitage's test for trend.
#' \emph{Human Heredity}, 52(3), 149--153. \doi{10.1159/000053370}.
#'
#' Edwards, B. J., Haynes, C., Levenstien, M. A., Finch, S. J., & Gordon, D.
#' (2005). Power and sample size calculations in the presence of phenotype
#' errors for case/control genetic association studies. \emph{BMC Genetics},
#' 6, 18. \doi{10.1186/1471-2156-6-18}.
#'
#' @seealso \code{\link{cc_power}},
#' \code{\link{case_control_genotype_misclassification}},
#' \code{\link{case_control_locus_heterogeneity}}, and
#' \code{\link{case_control_phenotype_misclassification}}.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export
cc_mssn <- function(
    power, alpha,

    input_mode = c("model_based", "model_free"),

    # model-based inputs
    prev = NULL,
    pd   = NULL,
    R2   = NULL,
    MOI  = c("M", "D", "Rec"),

    # model-free inputs
    g1 = NULL,
    g0 = NULL,

    # locus heterogeneity modifier
    locus_het = FALSE,
    pi = 1,

    # phenotype misclassification modifier
    pheno_misclass = FALSE,
    theta = 0,
    phi = 0,

    k = 1,
    w = c(0, 1, 2),

    # genotype misclassification controls
    geno_misclass = c("none", "1p", "2p", "3p", "diff3p"),

    # 1p
    e  = 0,

    # 2p
    e1 = 0,
    e2 = 0,

    # 3p non-differential
    e01 = 0,
    e02 = 0,
    e03 = 0,

    # diff3p separate case and control matrices
    case_e01 = 0,
    case_e02 = 0,
    case_e03 = 0,
    ctrl_e01 = 0,
    ctrl_e02 = 0,
    ctrl_e03 = 0,

    # diff3p multiplier shortcut
    diff_source = c("explicit", "case", "ctrl"),
    diff_multiplier = 1,

    verbose = TRUE
) {

  input_mode <- match.arg(input_mode)
  MOI <- match.arg(MOI)
  geno_misclass <- match.arg(geno_misclass)
  diff_source <- match.arg(diff_source)

  # ---- local helpers ----
  chisq_ncp_target <- function(power, alpha, df) {
    crit <- qchisq(1 - alpha, df = df)
    f <- function(lambda) {
      pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - power
    }
    uniroot(f, lower = 0, upper = 1e6)$root
  }

  check_genotype_freqs <- function(g, name = "g") {
    if (!is.numeric(g) || length(g) != 3)
      stop(name, " must be a numeric vector of length 3: c(g0, g1, g2).")
    if (any(!is.finite(g)))
      stop(name, " contains non-finite values.")
    if (any(g < 0))
      stop(name, " cannot contain negative genotype frequencies.")
    if (abs(sum(g) - 1) > 1e-6)
      stop(name, " must sum to 1.")
    invisible(TRUE)
  }

  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M", "D", "Rec")) {
    MOI <- match.arg(MOI)
    p_plus <- 1 - pd
    R1 <- if (MOI == "M") sqrt(R2) else if (MOI == "D") R2 else 1

    f0 <- prev / (p_plus^2 + R1 * 2 * pd * p_plus + R2 * pd^2)
    f1 <- R1 * f0
    f2 <- R2 * f0

    g_j1 <- c(
      f0 * p_plus^2 / prev,
      f1 * 2 * pd * p_plus / prev,
      f2 * pd^2 / prev
    )

    g_j0 <- c(
      (1 - f0) * p_plus^2 / (1 - prev),
      (1 - f1) * 2 * pd * p_plus / (1 - prev),
      (1 - f2) * pd^2 / (1 - prev)
    )

    list(
      R1 = R1,
      f0 = f0,
      f1 = f1,
      f2 = f2,
      g_j1 = g_j1,
      g_j0 = g_j0
    )
  }

  cc_apply_pheno_misclass <- function(g_aff, g_unaff, prev, theta, phi) {
    check_genotype_freqs(g_aff, "g_aff")
    check_genotype_freqs(g_unaff, "g_unaff")

    case_denom <- (1 - theta) * prev + phi * (1 - prev)
    ctrl_denom <- theta * prev + (1 - phi) * (1 - prev)

    if (case_denom <= 0 || ctrl_denom <= 0)
      stop("Observed case/control denominators must be positive.")

    g_case_obs <- (g_aff * (1 - theta) * prev + g_unaff * phi * (1 - prev)) /
      case_denom

    g_ctrl_obs <- (g_aff * theta * prev + g_unaff * (1 - phi) * (1 - prev)) /
      ctrl_denom

    list(
      g_case_obs = as.numeric(g_case_obs / sum(g_case_obs)),
      g_ctrl_obs = as.numeric(g_ctrl_obs / sum(g_ctrl_obs)),
      case_denom = case_denom,
      ctrl_denom = ctrl_denom
    )
  }

  cc_misclass_matrix_1p <- function(e) {
    if (!is.numeric(e) || length(e) != 1 || e < 0 || e > 0.5)
      stop("e must be a single number in [0, 0.5].")

    M <- matrix(c(
      1 - 2*e, e,       e,
      e,       1 - 2*e, e,
      e,       e,       1 - 2*e
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10))
      stop("1p misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("1p misclassification matrix has negative entries; check e.")

    M
  }

  cc_misclass_matrix_2p <- function(e1, e2) {
    if (!is.numeric(e1) || length(e1) != 1 || e1 < 0 || e1 > 1)
      stop("e1 must be a single number in [0,1].")
    if (!is.numeric(e2) || length(e2) != 1 || e2 < 0 || e2 > 0.5)
      stop("e2 must be a single number in [0,0.5].")

    M <- matrix(c(
      1 - e1,  e1,       0,
      e2,      1 - 2*e2, e2,
      0,       e1,       1 - e1
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10))
      stop("2p misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("2p misclassification matrix has negative entries; check e1/e2.")

    M
  }

  cc_misclass_matrix_3p <- function(e01, e02, e03) {
    if (!is.numeric(e01) || length(e01) != 1 || e01 < 0 || e01 > 1)
      stop("e01 must be a single number in [0,1].")
    if (!is.numeric(e02) || length(e02) != 1 || e02 < 0 || e02 > 0.5)
      stop("e02 must be a single number in [0,0.5].")
    if (!is.numeric(e03) || length(e03) != 1 || e03 < 0 || e03 > 1)
      stop("e03 must be a single number in [0,1].")
    if (e01 + e03 > 1)
      stop("Need e01 + e03 <= 1.")

    M <- matrix(c(
      1 - (e01 + e03),  e01,            e03,
      e02,              1 - 2*e02,      e02,
      e03,              e01,            1 - (e01 + e03)
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10))
      stop("3p misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("3p misclassification matrix has negative entries; check e01/e02/e03.")

    M
  }

  cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
    if (length(g_true) != 3)
      stop("g_true must be length 3.")
    as.numeric(t(M_true_to_obs) %*% g_true)
  }

  cc_scale_3p_errors <- function(e01, e02, e03, multiplier) {
    if (!is.numeric(multiplier) || length(multiplier) != 1 || multiplier < 0)
      stop("diff_multiplier must be a single nonnegative number.")

    out <- c(
      e01 = e01 * multiplier,
      e02 = e02 * multiplier,
      e03 = e03 * multiplier
    )

    if (out["e01"] < 0 || out["e01"] > 1)
      stop("Scaled e01 is outside [0,1].")
    if (out["e02"] < 0 || out["e02"] > 0.5)
      stop("Scaled e02 is outside [0,0.5].")
    if (out["e03"] < 0 || out["e03"] > 1)
      stop("Scaled e03 is outside [0,1].")
    if ((out["e01"] + out["e03"]) > 1)
      stop("Scaled e01 + e03 > 1.")

    out
  }

  .fmt_f <- function(x, digits = 3) {
    formatC(x, format = "f", digits = digits)
  }

  .fmt_e <- function(x, digits = 2) {
    formatC(x, format = "e", digits = digits)
  }

  # ---- checks ----
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("power must be a single number in (0,1).")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number: N_ctrl / N_case.")
  if (!is.numeric(w) || length(w) != 3)
    stop("w must be a numeric vector of length 3.")
  if (length(unique(w)) == 1)
    stop("Trend weights w cannot all be equal.")
  if (!is.logical(locus_het) || length(locus_het) != 1)
    stop("locus_het must be TRUE or FALSE.")
  if (!is.numeric(pi) || length(pi) != 1 || !is.finite(pi) || pi < 0 || pi > 1)
    stop("pi must be a single number in [0,1].")
  if (!isTRUE(locus_het) && pi != 1)
    stop("pi is used only when locus_het = TRUE; set pi = 1 or enable locus heterogeneity.")
  if (!is.logical(pheno_misclass) || length(pheno_misclass) != 1)
    stop("pheno_misclass must be TRUE or FALSE.")
  if (!is.numeric(theta) || length(theta) != 1 || theta < 0 || theta >= 1)
    stop("theta must be a single number in [0,1).")
  if (!is.numeric(phi) || length(phi) != 1 || phi < 0 || phi >= 1)
    stop("phi must be a single number in [0,1).")
  if (isTRUE(pheno_misclass) && (is.null(prev) || !is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1))
    stop("When pheno_misclass=TRUE, prev must be a single number in (0,1).")
  if (!is.numeric(diff_multiplier) || length(diff_multiplier) != 1 || diff_multiplier < 0)
    stop("diff_multiplier must be a single nonnegative number.")

  # ---- determine baseline true genotype frequencies ----
  if (input_mode == "model_based") {

    if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
      stop("For input_mode='model_based', prev must be a single number in (0,1).")
    if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
      stop("For input_mode='model_based', pd must be a single number in (0,1).")
    if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
      stop("For input_mode='model_based', R2 must be a single positive number.")

    freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
    g1_base <- freqs$g_j1
    g0_base <- freqs$g_j0

    model_info <- list(
      input_mode = "model_based",
      prev = prev,
      pd = pd,
      qd = 1 - pd,
      R1 = freqs$R1,
      R2 = R2,
      MOI = MOI,
      penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2)
    )

  } else if (input_mode == "model_free") {

    if (is.null(g1) || is.null(g0))
      stop("For input_mode='model_free', g1 and g0 must both be supplied.")

    check_genotype_freqs(g1, "g1")
    check_genotype_freqs(g0, "g0")

    g1_base <- as.numeric(g1)
    g0_base <- as.numeric(g0)

    model_info <- list(
      input_mode = "model_free"
    )
  }

  # ---- apply locus heterogeneity as optional modifier ----
  if (isTRUE(locus_het)) {
    g1_true <- pi * g1_base + (1 - pi) * g0_base
    g0_true <- g0_base
  } else {
    g1_true <- g1_base
    g0_true <- g0_base
  }

  locus_het_info <- list(
    enabled = locus_het,
    pi = pi,
    g_case_before_locus_het = g1_base,
    g_ctrl_before_locus_het = g0_base,
    g_case_after_locus_het = g1_true,
    g_ctrl_after_locus_het = g0_true
  )

  # ---- apply phenotype misclassification as optional modifier ----
  if (isTRUE(pheno_misclass)) {
    pheno <- cc_apply_pheno_misclass(
      g_aff = g1_true,
      g_unaff = g0_true,
      prev = prev,
      theta = theta,
      phi = phi
    )
    g1_true <- pheno$g_case_obs
    g0_true <- pheno$g_ctrl_obs
  } else {
    pheno <- list(
      g_case_obs = g1_true,
      g_ctrl_obs = g0_true,
      case_denom = NA_real_,
      ctrl_denom = NA_real_
    )
  }

  pheno_misclass_info <- list(
    enabled = pheno_misclass,
    theta = theta,
    phi = phi,
    g_case_after_pheno_misclass = g1_true,
    g_ctrl_after_pheno_misclass = g0_true,
    case_denom = pheno$case_denom,
    ctrl_denom = pheno$ctrl_denom
  )

# ---- choose genotype misclassification model ----
  M_case <- diag(3)
  M_ctrl <- diag(3)
  misclass_info <- list(enabled = FALSE, model = "none")

  if (geno_misclass == "1p") {

    M_case <- M_ctrl <- cc_misclass_matrix_1p(e)

    misclass_info <- list(
      enabled = (e > 0),
      model = "1p_symmetric",
      e = e,
      M = M_case
    )

  } else if (geno_misclass == "2p") {

    M_case <- M_ctrl <- cc_misclass_matrix_2p(e1, e2)

    misclass_info <- list(
      enabled = (e1 > 0 || e2 > 0),
      model = "2p_hom_het",
      e1 = e1,
      e2 = e2,
      M = M_case
    )

  } else if (geno_misclass == "3p") {

    M_case <- M_ctrl <- cc_misclass_matrix_3p(e01, e02, e03)

    misclass_info <- list(
      enabled = (e01 > 0 || e02 > 0 || e03 > 0),
      model = "3p_homhet_homhom",
      e01 = e01,
      e02 = e02,
      e03 = e03,
      M = M_case
    )

  } else if (geno_misclass == "diff3p") {

    if (diff_source == "explicit") {
      case_params <- c(e01 = case_e01, e02 = case_e02, e03 = case_e03)
      ctrl_params <- c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)

    } else if (diff_source == "case") {
      case_params <- c(e01 = case_e01, e02 = case_e02, e03 = case_e03)
      ctrl_params <- cc_scale_3p_errors(case_e01, case_e02, case_e03, diff_multiplier)

    } else if (diff_source == "ctrl") {
      ctrl_params <- c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)
      case_params <- cc_scale_3p_errors(ctrl_e01, ctrl_e02, ctrl_e03, diff_multiplier)
    }

    M_case <- cc_misclass_matrix_3p(
      case_params["e01"],
      case_params["e02"],
      case_params["e03"]
    )

    M_ctrl <- cc_misclass_matrix_3p(
      ctrl_params["e01"],
      ctrl_params["e02"],
      ctrl_params["e03"]
    )

    misclass_info <- list(
      enabled = any(c(case_params, ctrl_params) > 0),
      model = "diff3p_homhet_homhom",
      diff_source = diff_source,
      diff_multiplier = diff_multiplier,
      case_params = case_params,
      ctrl_params = ctrl_params,
      M_case = M_case,
      M_ctrl = M_ctrl
    )
  }

  # ---- observed genotype frequencies after misclassification ----
  g1_obs <- cc_apply_genotype_misclass(g1_true, M_case)
  g1_obs <- g1_obs / sum(g1_obs)

  g0_obs <- cc_apply_genotype_misclass(g0_true, M_ctrl)
  g0_obs <- g0_obs / sum(g0_obs)

  # ---- target lambdas ----
  lambda_star_g <- chisq_ncp_target(power = power, alpha = alpha, df = 2)
  lambda_star_1 <- chisq_ncp_target(power = power, alpha = alpha, df = 1)

  adjusted_tests <- .cc_mssn_test_results(
    g_case = g1_obs, g_ctrl = g0_obs, k = k, w = w,
    lambda_star_g = lambda_star_g, lambda_star_t = lambda_star_1
  )
  baseline_tests <- if (isTRUE(verbose)) {
    .cc_mssn_test_results(
      g_case = g1_base, g_ctrl = g0_base, k = k, w = w,
      lambda_star_g = lambda_star_g, lambda_star_t = lambda_star_1,
      validate = FALSE
    )
  } else {
    NULL
  }

  # ---- output ----
  out <- list(
    alpha = alpha,
    target_power = power,
    input_mode = input_mode,
    k = k,
    w = w,
    locus_het = locus_het_info,
    errors = list(
      phenotype_misclass = pheno_misclass_info,
      genotype_misclass = misclass_info
    ),
    model_info = model_info,
    tests = list(
      genotypes = list(
        test = "case-control chi-square test of independence for genotypes",
        df = 2,
        lambda_star = lambda_star_g,
        S = adjusted_tests$genotypes$S,
        MSSN_case = adjusted_tests$genotypes$MSSN_case,
        MSSN_ctrl = adjusted_tests$genotypes$MSSN_ctrl,
        MSSN_total = adjusted_tests$genotypes$MSSN_total
      ),
      trend = list(
        test = "trend test for genotypes",
        df = 1,
        lambda_star = lambda_star_1,
        S = adjusted_tests$trend$S,
        numerator = adjusted_tests$trend$numerator,
        denominator = adjusted_tests$trend$denominator,
        MSSN_case = adjusted_tests$trend$MSSN_case,
        MSSN_ctrl = adjusted_tests$trend$MSSN_ctrl,
        MSSN_total = adjusted_tests$trend$MSSN_total
      )
    ),
    freqs = list(
      g_base_case = g1_base,
      g_base_ctrl = g0_base,
      g_true_case = g1_true,
      g_true_ctrl = g0_true,
      g_after_pheno_case = pheno$g_case_obs,
      g_after_pheno_ctrl = pheno$g_ctrl_obs,
      g_obs_case  = g1_obs,
      g_obs_ctrl  = g0_obs
    )
  )

  class(out) <- "cc_mssn"

  # ---- clean printed output ----
  if (isTRUE(verbose)) {
    locus_active <- isTRUE(locus_het) && pi < 1
    pheno_active <- isTRUE(pheno_misclass) && (theta > 0 || phi > 0)
    geno_active <- .cc_genotype_error_active(M_case, M_ctrl)
    any_modifier <- locus_active || pheno_active || geno_active

    format_count <- function(x) {
      if (is.finite(x)) formatC(x, format = "d", big.mark = ",") else as.character(x)
    }

    message("Case-control minimum sample size")
    message(sprintf("Target power: %.1f%%", 100 * power))
    message("")
    message("No-error design")
    message(sprintf("  Genotype test: %s cases, %s controls, %s total",
                    format_count(baseline_tests$genotypes$MSSN_case),
                    format_count(baseline_tests$genotypes$MSSN_ctrl),
                    format_count(baseline_tests$genotypes$MSSN_total)))
    message(sprintf("  Trend test: %s cases, %s controls, %s total",
                    format_count(baseline_tests$trend$MSSN_case),
                    format_count(baseline_tests$trend$MSSN_ctrl),
                    format_count(baseline_tests$trend$MSSN_total)))

    if (any_modifier) {
      message("")
      message("Adjusted design")
      message("  Active modifiers:")
      if (locus_active) {
        message(sprintf("    Locus heterogeneity: %.1f%%", 100 * (1 - pi)))
      }
      if (pheno_active) {
        message(sprintf(
          "    Phenotype misclassification: theta %.1f%%, phi %.1f%%",
          100 * theta, 100 * phi
        ))
      }
      if (geno_active) {
        geno_label <- switch(
          geno_misclass,
          `1p` = "1-parameter", `2p` = "2-parameter",
          `3p` = "3-parameter", diff3p = "differential 3-parameter"
        )
        message(paste0("    Genotype misclassification: ", geno_label))
      }
      message(sprintf("  Genotype test: %s cases, %s controls, %s total",
                      format_count(adjusted_tests$genotypes$MSSN_case),
                      format_count(adjusted_tests$genotypes$MSSN_ctrl),
                      format_count(adjusted_tests$genotypes$MSSN_total)))
      message(sprintf("  Trend test: %s cases, %s controls, %s total",
                      format_count(adjusted_tests$trend$MSSN_case),
                      format_count(adjusted_tests$trend$MSSN_ctrl),
                      format_count(adjusted_tests$trend$MSSN_total)))
    }
  }

  invisible(out)
}





#' Case-Control Power for Conditional Genotype Frequencies
#'
#' Computes power for fixed case-control sample sizes using conditional genotype
#' frequencies. The function supports model-based genotype frequencies,
#' model-free genotype frequencies, optional locus heterogeneity, optional
#' phenotype misclassification, optional genotype misclassification, genotype
#' chi-square tests and genotype trend tests.
#'
#' @param N_case Numeric \eqn{> 0}. Number of cases.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param input_mode Character. One of \code{"model_based"} or
#'   \code{"model_free"}. See Details.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence for
#'   \code{input_mode = "model_based"}.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency for
#'   \code{input_mode = "model_based"}.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk for
#'   \code{input_mode = "model_based"}.
#' @param MOI Character. Mode of inheritance for model-based frequencies:
#'   \code{"M"} for multiplicative, \code{"D"} for dominant, or \code{"Rec"}
#'   for recessive.
#' @param g1,g0 Numeric vectors of length 3 for \code{input_mode = "model_free"}.
#'   \code{g1} gives case genotype frequencies and \code{g0} gives control
#'   genotype frequencies, each ordered as \code{c(g0, g1, g2)} and summing to 1.
#' @param locus_het Logical. If \code{TRUE}, applies locus heterogeneity to the
#'   case genotype frequencies before phenotype and genotype misclassification.
#' @param pi Numeric in \eqn{[0,1]}. Locus-homogeneity fraction used when
#'   \code{locus_het = TRUE}. When \code{locus_het = FALSE}, \code{pi} must
#'   remain at its default value of 1.
#' @param pheno_misclass Logical. If \code{TRUE}, applies phenotype
#'   misclassification before genotype misclassification.
#' @param theta Numeric in \eqn{[0,1)}. Probability that a truly affected
#'   individual is classified as a control.
#' @param phi Numeric in \eqn{[0,1)}. Probability that a truly unaffected
#'   individual is classified as a case.
#' @param k Numeric \eqn{> 0}. Control-to-case sample size ratio
#'   \eqn{N_{ctrl} / N_{case}}.
#' @param w Numeric vector of length 3. Genotype trend-test scores. The three
#'   weights cannot all be equal.
#' @param geno_misclass Character. Genotype misclassification model:
#'   \code{"none"}, \code{"1p"}, \code{"2p"}, \code{"3p"}, or \code{"diff3p"}.
#' @param e Numeric in \eqn{[0,0.5]}. For \code{geno_misclass = "1p"}, the
#'   probability assigned to each adjacent off-diagonal genotype call. The
#'   textbook Eq. 2.5 parameter satisfies \eqn{\epsilon=2e}; see Details.
#' @param e1,e2 Numeric. Two-parameter misclassification rates for
#'   \code{geno_misclass = "2p"}.
#' @param e01,e02,e03 Numeric. Non-differential three-parameter
#'   misclassification rates for \code{geno_misclass = "3p"}.
#' @param case_e01,case_e02,case_e03 Numeric. Case-specific three-parameter
#'   misclassification rates for \code{geno_misclass = "diff3p"}.
#' @param ctrl_e01,ctrl_e02,ctrl_e03 Numeric. Control-specific three-parameter
#'   misclassification rates for \code{geno_misclass = "diff3p"}.
#' @param diff_source Character. For \code{geno_misclass = "diff3p"}, one of
#'   \code{"explicit"}, \code{"case"}, or \code{"ctrl"}. See Details.
#' @param diff_multiplier Numeric \eqn{\ge 0}. Multiplier used when
#'   \code{diff_source = "case"} or \code{diff_source = "ctrl"}.
#' @param verbose Logical. If \code{TRUE}, prints a clean formatted summary.
#'
#' @details
#' The workflow is:
#' \enumerate{
#' \item Construct baseline conditional genotype frequencies for cases and
#' controls.
#' \item Optionally apply locus heterogeneity to cases as
#' \eqn{g_{case,true} = \pi g_{case,base} + (1 - \pi) g_{ctrl,base}}.
#' \item Optionally apply phenotype misclassification, where
#' \code{theta = Pr(affected -> control)} and
#' \code{phi = Pr(unaffected -> case)}.
#' \item Optionally apply genotype misclassification matrices to the resulting
#' case and control genotype frequencies.
#' \item Compute genotype chi-square and genotype trend-test non-centrality
#' parameters and powers from the observed genotype frequencies.
#' }
#' The defaults (\code{locus_het = FALSE}, \code{pheno_misclass = FALSE}, and
#' \code{geno_misclass = "none"}) give the ordinary no-error design. When
#' modifiers are active, they are applied sequentially in the order above and
#' therefore form one combined adjusted design.
#'
#' With \code{input_mode = "model_based"}, conditional case and control genotype
#' frequencies are derived from \code{prev}, \code{pd}, \code{R2}, and
#' \code{MOI} using Chapter 1, Section 1.4.2, Eqs. 1.6--1.7 (p. 13) of
#' Gordon, Finch, and Kim (2020). With \code{input_mode = "model_free"}, the user supplies
#' \code{g1} and \code{g0} directly; when phenotype misclassification is
#' enabled, \code{g1} is treated as the true affected genotype distribution and
#' \code{g0} is treated as the true unaffected genotype distribution.
#'
#' Phenotype misclassification requires \code{prev} in both input modes because
#' disease prevalence is used to mix the true affected and unaffected genotype
#' distributions into observed case and control genotype distributions. It is
#' applied after optional locus heterogeneity and before optional genotype
#' misclassification.
#'
#' The genotype misclassification models are:
#' \code{"none"} for identity matrices, \code{"1p"} for one symmetric error
#' rate, \code{"2p"} for adjacent homozygote/heterozygote and heterozygote
#' error rates, \code{"3p"} for a non-differential three-parameter matrix, and
#' \code{"diff3p"} for separate case and control three-parameter matrices.
#' These matrices correspond to Eqs. 2.5--2.7 (pp. 57--58). For \code{"1p"},
#' package \code{e} is the probability in each adjacent off-diagonal cell and
#' the diagonal is \eqn{1-2e}; textbook Eq. 2.5 uses off-diagonal
#' \eqn{\epsilon/2}, so \eqn{\epsilon=2e}. The genotype and trend tests use
#' Eqs. 1.22 and 1.24 (pp. 26--27), with the trend statistic defined by
#' Eqs. 1.20--1.21 (p. 24). Locus heterogeneity follows Eq. 2.16 (p. 88), where
#' \code{pi} is the homogeneous fraction; its trend-test form is given in
#' Eqs. 5.29a--b (pp. 287--288).
#'
#' For \code{geno_misclass = "diff3p"}, \code{diff_source = "explicit"} uses
#' the case and control error parameters exactly as supplied. With
#' \code{diff_source = "case"}, control parameters are computed by multiplying
#' the case parameters by \code{diff_multiplier}. With
#' \code{diff_source = "ctrl"}, case parameters are computed by multiplying the
#' control parameters by \code{diff_multiplier}.
#'
#' With \code{geno_misclass = "diff3p"}, returned power is nominal asymptotic
#' power evaluated with the usual chi-square critical value. Different case
#' and control error mechanisms can distort the null distribution and inflate
#' type I error. This function does not independently recalibrate the null
#' distribution for arbitrary differential genotyping error.
#'
#' Internal effect-size components \code{S} are retained in the returned object
#' for validation and debugging, but are not printed in the clean verbose
#' output.
#'
#' @return An object of class \code{"cc_power"}, containing:
#' \describe{
#' \item{alpha, N_case, N_ctrl, N_total}{Significance level and numbers of
#' case, control, and total individuals.}
#' \item{input_mode, k, w, locus_het}{Input mode, control-to-case ratio, trend
#' scores, and locus-heterogeneity settings.}
#' \item{errors}{Genotype-error model and matrices, plus phenotype-error
#' settings and intermediate frequencies.}
#' \item{model_info}{Model-based penetrances and risk-model inputs, or
#' model-free identifying information.}
#' \item{tests$genotypes, tests$trend}{Test label, degrees of freedom, NCP
#' \code{lambda}, internal \code{S}, and power. The trend result also contains
#' its numerator and denominator.}
#' \item{freqs}{Baseline, post-heterogeneity (true), post-phenotype-error, and
#' final observed case and control genotype-probability vectors.}
#' }
#'
#' @examples
#' cc_power(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   geno_misclass = "none",
#'   verbose = FALSE
#' )
#'
#' cc_power(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   verbose = FALSE
#' )
#'
#' cc_power(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' cc_power(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_based",
#'   prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M",
#'   pheno_misclass = TRUE, theta = 0.05, phi = 0.01,
#'   verbose = FALSE
#' )
#'
#' cc_power(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_based",
#'   prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M",
#'   locus_het = TRUE, pi = 0.8,
#'   geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' cc_power(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   geno_misclass = "diff3p",
#'   diff_source = "explicit",
#'   case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
#'   ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002,
#'   verbose = FALSE
#' )
#'
#' cc_power(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   geno_misclass = "diff3p",
#'   diff_source = "case", diff_multiplier = 0.5,
#'   case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Kim, W. (2020).
#' \emph{Heterogeneity in Statistical Genetics: How to Assess, Address, and
#' Account for Mixtures in Association Studies}. Springer.
#' \doi{10.1007/978-3-030-61121-7}.
#'
#' Gordon, D., Finch, S. J., Nothnagel, M., & Ott, J. (2002). Power and sample
#' size calculations for case-control genetic association tests when errors
#' are present: application to single nucleotide polymorphisms.
#' \emph{Human Heredity}, 54(1), 22--33. \doi{10.1159/000066696}.
#'
#' Armitage, P. (1955). Tests for linear trends in proportions and frequencies.
#' \emph{Biometrics}, 11(3), 375--386. \doi{10.2307/3001775}.
#'
#' Slager, S. L., & Schaid, D. J. (2001). Case-control studies of genetic
#' markers: power and sample size approximations for Armitage's test for trend.
#' \emph{Human Heredity}, 52(3), 149--153. \doi{10.1159/000053370}.
#'
#' Edwards, B. J., Haynes, C., Levenstien, M. A., Finch, S. J., & Gordon, D.
#' (2005). Power and sample size calculations in the presence of phenotype
#' errors for case/control genetic association studies. \emph{BMC Genetics},
#' 6, 18. \doi{10.1186/1471-2156-6-18}.
#'
#' @seealso \code{\link{cc_mssn}},
#' \code{\link{case_control_genotype_misclassification}},
#' \code{\link{case_control_locus_heterogeneity}}, and
#' \code{\link{case_control_phenotype_misclassification}}.
#'
#' @importFrom stats pchisq qchisq
#' @export
cc_power <- function(
    N_case, alpha,

    input_mode = c("model_based", "model_free"),

    # model-based inputs
    prev = NULL,
    pd   = NULL,
    R2   = NULL,
    MOI  = c("M", "D", "Rec"),

    # model-free inputs
    g1 = NULL,
    g0 = NULL,

    # locus heterogeneity modifier
    locus_het = FALSE,
    pi = 1,

    # phenotype misclassification modifier
    pheno_misclass = FALSE,
    theta = 0,
    phi = 0,

    k = 1,
    w = c(0, 1, 2),

    # genotype misclassification controls
    geno_misclass = c("none", "1p", "2p", "3p", "diff3p"),

    # 1p
    e  = 0,

    # 2p
    e1 = 0,
    e2 = 0,

    # 3p non-differential
    e01 = 0,
    e02 = 0,
    e03 = 0,

    # diff3p separate case and control matrices
    case_e01 = 0,
    case_e02 = 0,
    case_e03 = 0,
    ctrl_e01 = 0,
    ctrl_e02 = 0,
    ctrl_e03 = 0,

    # diff3p multiplier shortcut
    diff_source = c("explicit", "case", "ctrl"),
    diff_multiplier = 1,

    verbose = TRUE
) {

  input_mode <- match.arg(input_mode)
  MOI <- match.arg(MOI)
  geno_misclass <- match.arg(geno_misclass)
  diff_source <- match.arg(diff_source)

  # ---- local helpers ----
  check_genotype_freqs <- function(g, name = "g") {
    if (!is.numeric(g) || length(g) != 3)
      stop(name, " must be a numeric vector of length 3: c(g0, g1, g2).")
    if (any(!is.finite(g)))
      stop(name, " contains non-finite values.")
    if (any(g < 0))
      stop(name, " cannot contain negative genotype frequencies.")
    if (abs(sum(g) - 1) > 1e-6)
      stop(name, " must sum to 1.")
    invisible(TRUE)
  }

  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M", "D", "Rec")) {
    MOI <- match.arg(MOI)
    p_plus <- 1 - pd
    R1 <- if (MOI == "M") sqrt(R2) else if (MOI == "D") R2 else 1

    f0 <- prev / (p_plus^2 + R1 * 2 * pd * p_plus + R2 * pd^2)
    f1 <- R1 * f0
    f2 <- R2 * f0

    g_j1 <- c(
      f0 * p_plus^2 / prev,
      f1 * 2 * pd * p_plus / prev,
      f2 * pd^2 / prev
    )

    g_j0 <- c(
      (1 - f0) * p_plus^2 / (1 - prev),
      (1 - f1) * 2 * pd * p_plus / (1 - prev),
      (1 - f2) * pd^2 / (1 - prev)
    )

    list(
      R1 = R1,
      f0 = f0,
      f1 = f1,
      f2 = f2,
      g_j1 = g_j1,
      g_j0 = g_j0
    )
  }

  cc_apply_pheno_misclass <- function(g_aff, g_unaff, prev, theta, phi) {
    check_genotype_freqs(g_aff, "g_aff")
    check_genotype_freqs(g_unaff, "g_unaff")

    case_denom <- (1 - theta) * prev + phi * (1 - prev)
    ctrl_denom <- theta * prev + (1 - phi) * (1 - prev)

    if (case_denom <= 0 || ctrl_denom <= 0)
      stop("Observed case/control denominators must be positive.")

    g_case_obs <- (g_aff * (1 - theta) * prev + g_unaff * phi * (1 - prev)) /
      case_denom

    g_ctrl_obs <- (g_aff * theta * prev + g_unaff * (1 - phi) * (1 - prev)) /
      ctrl_denom

    list(
      g_case_obs = as.numeric(g_case_obs / sum(g_case_obs)),
      g_ctrl_obs = as.numeric(g_ctrl_obs / sum(g_ctrl_obs)),
      case_denom = case_denom,
      ctrl_denom = ctrl_denom
    )
  }

  cc_misclass_matrix_1p <- function(e) {
    if (!is.numeric(e) || length(e) != 1 || e < 0 || e > 0.5)
      stop("e must be a single number in [0, 0.5].")

    M <- matrix(c(
      1 - 2*e, e,       e,
      e,       1 - 2*e, e,
      e,       e,       1 - 2*e
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10))
      stop("1p misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("1p misclassification matrix has negative entries; check e.")

    M
  }

  cc_misclass_matrix_2p <- function(e1, e2) {
    if (!is.numeric(e1) || length(e1) != 1 || e1 < 0 || e1 > 1)
      stop("e1 must be a single number in [0,1].")
    if (!is.numeric(e2) || length(e2) != 1 || e2 < 0 || e2 > 0.5)
      stop("e2 must be a single number in [0,0.5].")

    M <- matrix(c(
      1 - e1,  e1,       0,
      e2,      1 - 2*e2, e2,
      0,       e1,       1 - e1
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10))
      stop("2p misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("2p misclassification matrix has negative entries; check e1/e2.")

    M
  }

  cc_misclass_matrix_3p <- function(e01, e02, e03) {
    if (!is.numeric(e01) || length(e01) != 1 || e01 < 0 || e01 > 1)
      stop("e01 must be a single number in [0,1].")
    if (!is.numeric(e02) || length(e02) != 1 || e02 < 0 || e02 > 0.5)
      stop("e02 must be a single number in [0,0.5].")
    if (!is.numeric(e03) || length(e03) != 1 || e03 < 0 || e03 > 1)
      stop("e03 must be a single number in [0,1].")
    if (e01 + e03 > 1)
      stop("Need e01 + e03 <= 1.")

    M <- matrix(c(
      1 - (e01 + e03),  e01,            e03,
      e02,              1 - 2*e02,      e02,
      e03,              e01,            1 - (e01 + e03)
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10))
      stop("3p misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("3p misclassification matrix has negative entries; check e01/e02/e03.")

    M
  }

  cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
    if (length(g_true) != 3)
      stop("g_true must be length 3.")
    as.numeric(t(M_true_to_obs) %*% g_true)
  }

  cc_scale_3p_errors <- function(e01, e02, e03, multiplier) {
    if (!is.numeric(multiplier) || length(multiplier) != 1 || multiplier < 0)
      stop("diff_multiplier must be a single nonnegative number.")

    out <- c(
      e01 = e01 * multiplier,
      e02 = e02 * multiplier,
      e03 = e03 * multiplier
    )

    if (out["e01"] < 0 || out["e01"] > 1)
      stop("Scaled e01 is outside [0,1].")
    if (out["e02"] < 0 || out["e02"] > 0.5)
      stop("Scaled e02 is outside [0,0.5].")
    if (out["e03"] < 0 || out["e03"] > 1)
      stop("Scaled e03 is outside [0,1].")
    if ((out["e01"] + out["e03"]) > 1)
      stop("Scaled e01 + e03 > 1.")

    out
  }

  .fmt_f <- function(x, digits = 3) {
    formatC(x, format = "f", digits = digits)
  }

  .fmt_e <- function(x, digits = 2) {
    formatC(x, format = "e", digits = digits)
  }

  # ---- checks ----
  if (!is.numeric(N_case) || length(N_case) != 1 || N_case <= 0)
    stop("N_case must be a single positive number.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number: N_ctrl / N_case.")
  if (!is.numeric(w) || length(w) != 3)
    stop("w must be a numeric vector of length 3.")
  if (length(unique(w)) == 1)
    stop("Trend weights w cannot all be equal.")
  if (!is.logical(locus_het) || length(locus_het) != 1)
    stop("locus_het must be TRUE or FALSE.")
  if (!is.numeric(pi) || length(pi) != 1 || !is.finite(pi) || pi < 0 || pi > 1)
    stop("pi must be a single number in [0,1].")
  if (!isTRUE(locus_het) && pi != 1)
    stop("pi is used only when locus_het = TRUE; set pi = 1 or enable locus heterogeneity.")
  if (!is.logical(pheno_misclass) || length(pheno_misclass) != 1)
    stop("pheno_misclass must be TRUE or FALSE.")
  if (!is.numeric(theta) || length(theta) != 1 || theta < 0 || theta >= 1)
    stop("theta must be a single number in [0,1).")
  if (!is.numeric(phi) || length(phi) != 1 || phi < 0 || phi >= 1)
    stop("phi must be a single number in [0,1).")
  if (isTRUE(pheno_misclass) && (is.null(prev) || !is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1))
    stop("When pheno_misclass=TRUE, prev must be a single number in (0,1).")
  if (!is.numeric(diff_multiplier) || length(diff_multiplier) != 1 || diff_multiplier < 0)
    stop("diff_multiplier must be a single nonnegative number.")

  # ---- determine baseline true genotype frequencies ----
  if (input_mode == "model_based") {

    if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
      stop("For input_mode='model_based', prev must be a single number in (0,1).")
    if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
      stop("For input_mode='model_based', pd must be a single number in (0,1).")
    if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
      stop("For input_mode='model_based', R2 must be a single positive number.")

    freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
    g1_base <- freqs$g_j1
    g0_base <- freqs$g_j0

    model_info <- list(
      input_mode = "model_based",
      prev = prev,
      pd = pd,
      qd = 1 - pd,
      R1 = freqs$R1,
      R2 = R2,
      MOI = MOI,
      penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2)
    )

  } else if (input_mode == "model_free") {

    if (is.null(g1) || is.null(g0))
      stop("For input_mode='model_free', g1 and g0 must both be supplied.")

    check_genotype_freqs(g1, "g1")
    check_genotype_freqs(g0, "g0")

    g1_base <- as.numeric(g1)
    g0_base <- as.numeric(g0)

    model_info <- list(
      input_mode = "model_free"
    )
  }

  # ---- apply locus heterogeneity as optional modifier ----
  if (isTRUE(locus_het)) {
    g1_true <- pi * g1_base + (1 - pi) * g0_base
    g0_true <- g0_base
  } else {
    g1_true <- g1_base
    g0_true <- g0_base
  }

  locus_het_info <- list(
    enabled = locus_het,
    pi = pi,
    g_case_before_locus_het = g1_base,
    g_ctrl_before_locus_het = g0_base,
    g_case_after_locus_het = g1_true,
    g_ctrl_after_locus_het = g0_true
  )

  # ---- apply phenotype misclassification as optional modifier ----
  if (isTRUE(pheno_misclass)) {
    pheno <- cc_apply_pheno_misclass(
      g_aff = g1_true,
      g_unaff = g0_true,
      prev = prev,
      theta = theta,
      phi = phi
    )
    g1_true <- pheno$g_case_obs
    g0_true <- pheno$g_ctrl_obs
  } else {
    pheno <- list(
      g_case_obs = g1_true,
      g_ctrl_obs = g0_true,
      case_denom = NA_real_,
      ctrl_denom = NA_real_
    )
  }

  pheno_misclass_info <- list(
    enabled = pheno_misclass,
    theta = theta,
    phi = phi,
    g_case_after_pheno_misclass = g1_true,
    g_ctrl_after_pheno_misclass = g0_true,
    case_denom = pheno$case_denom,
    ctrl_denom = pheno$ctrl_denom
  )

# ---- choose genotype misclassification model ----
  M_case <- diag(3)
  M_ctrl <- diag(3)
  misclass_info <- list(enabled = FALSE, model = "none")

  if (geno_misclass == "1p") {

    M_case <- M_ctrl <- cc_misclass_matrix_1p(e)

    misclass_info <- list(
      enabled = (e > 0),
      model = "1p_symmetric",
      e = e,
      M = M_case
    )

  } else if (geno_misclass == "2p") {

    M_case <- M_ctrl <- cc_misclass_matrix_2p(e1, e2)

    misclass_info <- list(
      enabled = (e1 > 0 || e2 > 0),
      model = "2p_hom_het",
      e1 = e1,
      e2 = e2,
      M = M_case
    )

  } else if (geno_misclass == "3p") {

    M_case <- M_ctrl <- cc_misclass_matrix_3p(e01, e02, e03)

    misclass_info <- list(
      enabled = (e01 > 0 || e02 > 0 || e03 > 0),
      model = "3p_homhet_homhom",
      e01 = e01,
      e02 = e02,
      e03 = e03,
      M = M_case
    )

  } else if (geno_misclass == "diff3p") {

    if (diff_source == "explicit") {
      case_params <- c(e01 = case_e01, e02 = case_e02, e03 = case_e03)
      ctrl_params <- c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)

    } else if (diff_source == "case") {
      case_params <- c(e01 = case_e01, e02 = case_e02, e03 = case_e03)
      ctrl_params <- cc_scale_3p_errors(case_e01, case_e02, case_e03, diff_multiplier)

    } else if (diff_source == "ctrl") {
      ctrl_params <- c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)
      case_params <- cc_scale_3p_errors(ctrl_e01, ctrl_e02, ctrl_e03, diff_multiplier)
    }

    M_case <- cc_misclass_matrix_3p(
      case_params["e01"],
      case_params["e02"],
      case_params["e03"]
    )

    M_ctrl <- cc_misclass_matrix_3p(
      ctrl_params["e01"],
      ctrl_params["e02"],
      ctrl_params["e03"]
    )

    misclass_info <- list(
      enabled = any(c(case_params, ctrl_params) > 0),
      model = "diff3p_homhet_homhom",
      diff_source = diff_source,
      diff_multiplier = diff_multiplier,
      case_params = case_params,
      ctrl_params = ctrl_params,
      M_case = M_case,
      M_ctrl = M_ctrl
    )
  }

  # ---- observed genotype frequencies after misclassification ----
  g1_obs <- cc_apply_genotype_misclass(g1_true, M_case)
  g1_obs <- g1_obs / sum(g1_obs)

  g0_obs <- cc_apply_genotype_misclass(g0_true, M_ctrl)
  g0_obs <- g0_obs / sum(g0_obs)

  # ---- sample sizes ----
  N_ctrl <- k * N_case

  adjusted_tests <- .cc_power_test_results(
    g_case = g1_obs, g_ctrl = g0_obs, k = k, w = w,
    N_case = N_case, alpha = alpha
  )
  baseline_tests <- if (isTRUE(verbose)) {
    .cc_power_test_results(
      g_case = g1_base, g_ctrl = g0_base, k = k, w = w,
      N_case = N_case, alpha = alpha, validate = FALSE
    )
  } else {
    NULL
  }

  # ---- output ----
  out <- list(
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    input_mode = input_mode,
    k = k,
    w = w,
    locus_het = locus_het_info,
    errors = list(
      phenotype_misclass = pheno_misclass_info,
      genotype_misclass = misclass_info
    ),
    model_info = model_info,
    tests = list(
      genotypes = list(
        test = "case-control chi-square test of independence for genotypes",
        df = 2,
        lambda = adjusted_tests$genotypes$lambda,
        S = adjusted_tests$genotypes$S,
        power = adjusted_tests$genotypes$power
      ),
      trend = list(
        test = "trend test for genotypes",
        df = 1,
        lambda = adjusted_tests$trend$lambda,
        S = adjusted_tests$trend$S,
        numerator = adjusted_tests$trend$numerator,
        denominator = adjusted_tests$trend$denominator,
        power = adjusted_tests$trend$power
      )
    ),
    freqs = list(
      g_base_case = g1_base,
      g_base_ctrl = g0_base,
      g_true_case = g1_true,
      g_true_ctrl = g0_true,
      g_after_pheno_case = pheno$g_case_obs,
      g_after_pheno_ctrl = pheno$g_ctrl_obs,
      g_obs_case  = g1_obs,
      g_obs_ctrl  = g0_obs
    )
  )

  class(out) <- "cc_power"

  # ---- clean printed output ----
  if (isTRUE(verbose)) {
    locus_active <- isTRUE(locus_het) && pi < 1
    pheno_active <- isTRUE(pheno_misclass) && (theta > 0 || phi > 0)
    geno_active <- .cc_genotype_error_active(M_case, M_ctrl)
    any_modifier <- locus_active || pheno_active || geno_active

    message("Case-control power")
    message(sprintf("Cases: %s; controls: %s; total: %s",
                    formatC(N_case, format = "f", digits = 0, big.mark = ","),
                    formatC(N_ctrl, format = "f", digits = 0, big.mark = ","),
                    formatC(N_case + N_ctrl, format = "f", digits = 0,
                            big.mark = ",")))
    message("")
    message("No-error design")
    message(sprintf("  Genotype test power: %.1f%%",
                    100 * baseline_tests$genotypes$power))
    message(sprintf("  Trend test power: %.1f%%",
                    100 * baseline_tests$trend$power))

    if (any_modifier) {
      message("")
      message("Adjusted design")
      message("  Active modifiers:")
      if (locus_active) {
        message(sprintf("    Locus heterogeneity: %.1f%%", 100 * (1 - pi)))
      }
      if (pheno_active) {
        message(sprintf(
          "    Phenotype misclassification: theta %.1f%%, phi %.1f%%",
          100 * theta, 100 * phi
        ))
      }
      if (geno_active) {
        geno_label <- switch(
          geno_misclass,
          `1p` = "1-parameter", `2p` = "2-parameter",
          `3p` = "3-parameter", diff3p = "differential 3-parameter"
        )
        message(paste0("    Genotype misclassification: ", geno_label))
      }
      message(sprintf("  Genotype test power: %.1f%%",
                      100 * adjusted_tests$genotypes$power))
      message(sprintf("  Trend test power: %.1f%%",
                      100 * adjusted_tests$trend$power))
    }
  }

  invisible(out)
}
