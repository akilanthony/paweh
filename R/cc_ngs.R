.cc_model_genotype_frequencies <- function(pd, R2, MOI, prev) {
  MOI <- match.arg(MOI, c("M", "D", "Rec"))

  if (!is.numeric(prev) || length(prev) != 1L || !is.finite(prev) ||
      prev <= 0 || prev >= 1) {
    stop("prev must be a single finite number in (0, 1).")
  }
  if (!is.numeric(pd) || length(pd) != 1L || !is.finite(pd) ||
      pd <= 0 || pd >= 1) {
    stop("pd must be a single finite number in (0, 1).")
  }
  if (!is.numeric(R2) || length(R2) != 1L || !is.finite(R2) || R2 <= 0) {
    stop("R2 must be a single finite positive number.")
  }

  p_plus <- 1 - pd
  R1 <- if (MOI == "M") sqrt(R2) else if (MOI == "D") R2 else 1

  f0 <- prev / (p_plus^2 + R1 * 2 * pd * p_plus + R2 * pd^2)
  f1 <- R1 * f0
  f2 <- R2 * f0

  population <- c(p_plus^2, 2 * pd * p_plus, pd^2)
  case <- c(
    f0 * p_plus^2 / prev,
    f1 * 2 * pd * p_plus / prev,
    f2 * pd^2 / prev
  )
  control <- c(
    (1 - f0) * p_plus^2 / (1 - prev),
    (1 - f1) * 2 * pd * p_plus / (1 - prev),
    (1 - f2) * pd^2 / (1 - prev)
  )

  list(
    case = as.numeric(case),
    control = as.numeric(control),
    population = as.numeric(population),
    penetrances = c(f0 = f0, f1 = f1, f2 = f2),
    R1 = R1,
    R2 = R2,
    MOI = MOI,
    pd = pd,
    prev = prev
  )
}

.cc_ngs_scores_from_moi <- function(MOI) {
  MOI <- match.arg(MOI, c("M", "D", "Rec"))
  switch(
    MOI,
    M = c(0, 1, 2),
    D = c(0, 1, 1),
    Rec = c(0, 0, 1)
  )
}

.cc_ngs_apply_locus_heterogeneity <- function(
    g_case,
    g_control,
    locus_het = FALSE,
    pi = 1
) {
  .cc_ngs_validate_genotype_frequencies(g_case, "g_case")
  .cc_ngs_validate_genotype_frequencies(g_control, "g_control")
  if (!is.logical(locus_het) || length(locus_het) != 1L || is.na(locus_het)) {
    stop("locus_het must be TRUE or FALSE.")
  }
  if (!is.numeric(pi) || length(pi) != 1L || !is.finite(pi) ||
      pi < 0 || pi > 1) {
    stop("pi must be a single finite number in [0, 1].")
  }
  if (!isTRUE(locus_het) && pi != 1) {
    stop(
      "pi is used only when locus_het = TRUE; set pi = 1 or enable locus heterogeneity."
    )
  }

  effective_pi <- if (isTRUE(locus_het)) pi else 1
  adjusted <- cc_apply_locus_het(
    g_case_assoc = g_case,
    g_ctrl = g_control,
    pi = effective_pi
  )

  list(
    enabled = locus_het,
    pi = pi,
    effective_pi = effective_pi,
    g_case_before_locus_het = as.numeric(g_case),
    g_ctrl_before_locus_het = as.numeric(g_control),
    g_case_after_locus_het = as.numeric(adjusted$g_case_het),
    g_ctrl_after_locus_het = as.numeric(adjusted$g_ctrl_het)
  )
}

.cc_ngs_chisq_power <- function(lambda, alpha) {
  critical <- qchisq(1 - alpha, df = 1)
  as.numeric(pchisq(
    critical,
    df = 1,
    ncp = lambda,
    lower.tail = FALSE
  ))
}

.cc_ngs_target_ncp <- function(power, alpha) {
  if (power <= alpha) {
    return(0)
  }
  cc_chisq_ncp_target(power = power, alpha = alpha, df = 1)
}

.cc_ngs_mssn_components <- function(g_case, g_control, k, scores,
                                    lambda_target) {
  .cc_ngs_validate_genotype_frequencies(g_case, "g_case")
  .cc_ngs_validate_genotype_frequencies(g_control, "g_control")
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("k must be a single finite positive number.")
  }
  if (!is.numeric(scores) || length(scores) != 3L ||
      any(!is.finite(scores)) || length(unique(scores)) == 1L) {
    stop("scores must be a finite, nonconstant numeric vector of length 3.")
  }
  if (!is.numeric(lambda_target) || length(lambda_target) != 1L ||
      !is.finite(lambda_target) || lambda_target < 0) {
    stop("lambda_target must be a single finite nonnegative number.")
  }

  D <- sum(scores * (g_case - g_control))
  pooled <- g_case + k * g_control
  Q <- sum(scores^2 * pooled) -
    sum(scores * pooled)^2 / (1 + k)

  if (!is.finite(Q) || Q <= 0) {
    stop("Ahn trend-test Q must be finite and positive.")
  }
  if (D^2 < 1e-15) {
    if (lambda_target == 0) {
      return(list(D = D, Q = Q, N_case_continuous = 0))
    }
    stop(
      "No finite MSSN exists because the trend contrast is zero under ",
      "this design."
    )
  }

  list(
    D = D,
    Q = Q,
    N_case_continuous = lambda_target * Q / (k * D^2)
  )
}

#' Analytic Power for a Case-Control Sequencing Study
#'
#' Computes prospective asymptotic power for a model-based case-control
#' sequencing design. The calculation constructs true case and control
#' genotype probabilities, applies a fixed-depth symmetric sequencing-error
#' model with deterministic maximum-likelihood genotype calls, and evaluates
#' the Ahn/Chapman-Nam Cochran-Armitage trend-test noncentrality parameter on
#' the resulting called-genotype probabilities.
#'
#' @param N_case Numeric \eqn{> 0}. Number of cases.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk.
#' @param coverage Positive integer sequencing depth.
#' @param seq_error Symmetric per-read sequencing-error probability in
#'   \eqn{[0,0.5)}.
#' @param MOI Character mode of inheritance: \code{"M"} for multiplicative,
#'   \code{"D"} for dominant, or \code{"Rec"} for recessive.
#' @param k Numeric \eqn{> 0}. Control-to-case sample-size ratio
#'   \eqn{N_{ctrl}/N_{case}}.
#' @param verbose Logical. If \code{TRUE}, print a concise result summary.
#' @param locus_het Logical. If \code{TRUE}, apply the canonical PAWEH
#'   case-control locus-heterogeneity mixture before sequencing observation.
#' @param pi Numeric in \eqn{[0,1]}. Locus-homogeneity fraction used when
#'   \code{locus_het = TRUE}. One retains the original associated-case
#'   distribution; zero makes the case distribution equal to controls. When
#'   \code{locus_het = FALSE}, \code{pi} must remain at its default value of 1.
#'
#' @details
#' Trend scores are selected from \code{MOI}: \code{"M"} uses
#' \code{c(0,1,2)}, \code{"D"} uses \code{c(0,1,1)}, and \code{"Rec"} uses
#' \code{c(0,0,1)}. Power is the upper-tail probability beyond the central
#' one-degree-of-freedom chi-square critical value under a noncentral
#' chi-square distribution with the calculated NCP.
#'
#' Locus heterogeneity uses the same parameterization as ordinary PAWEH
#' case-control design:
#' \deqn{g_{case,H} = \pi g_{case} + (1-\pi)g_{control},}
#' with the control distribution unchanged. This biological mixture is applied
#' to true genotype probabilities before sequencing observation. Because the
#' current sequencing-error model is nondifferential, the same transition
#' matrix is applied to cases and controls, so mixing and sequencing commute by
#' matrix linearity. This identity does not extend automatically to future
#' differential case/control sequencing-error models. At \eqn{\pi=0}, the
#' case-control contrast and NCP are zero, and asymptotic power equals
#' \code{alpha}.
#'
#' This is an analytic study-design calculation applied to sequencing-derived
#' called genotypes. It is not a raw-read likelihood or EM analysis, performs
#' no downstream association testing, and uses no simulation. Even when
#' \code{seq_error = 0}, finite depth can cause call uncertainty because a true
#' heterozygote can yield reads from only one allele.
#'
#' @return Invisibly, an object of class \code{"cc_ngs_power"} containing the
#'   design inputs, trend scores, NCP and power, model information, true and
#'   called genotype frequencies, and the true-to-called transition matrix.
#'
#' @references
#' Ahn, K., Haynes, C., Kim, W., St. Fleur, R., Gordon, D., & Finch, S. J.
#' (2007). The effects of SNP genotyping errors on the power of the
#' Cochran-Armitage linear trend test for case/control association studies.
#' \emph{Annals of Human Genetics}, 71, 249--261.
#' \doi{10.1111/j.1469-1809.2006.00318.x}.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020). \emph{Heterogeneity in
#' Statistical Genetics}. Springer. \doi{10.1007/978-3-030-61121-7}.
#'
#' @examples
#' cc_ngs_power(
#'   N_case = 1000, alpha = 0.05,
#'   prev = 0.05, pd = 0.30, R2 = 1.8,
#'   coverage = 20, seq_error = 0.01,
#'   MOI = "M", verbose = FALSE
#' )
#'
#' @importFrom stats pchisq qchisq
#' @export
cc_ngs_power <- function(
    N_case, alpha,
    prev, pd, R2,
    coverage, seq_error,
    MOI = c("M", "D", "Rec"),
    k = 1,
    verbose = TRUE,
    locus_het = FALSE,
    pi = 1
) {
  MOI <- match.arg(MOI)

  if (!is.numeric(N_case) || length(N_case) != 1L || !is.finite(N_case) ||
      N_case <= 0) {
    stop("N_case must be a single finite positive number.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number in (0, 1).")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("k must be a single finite positive number: N_ctrl / N_case.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  model <- .cc_model_genotype_frequencies(
    pd = pd,
    R2 = R2,
    MOI = MOI,
    prev = prev
  )
  heterogeneity <- .cc_ngs_apply_locus_heterogeneity(
    g_case = model$case,
    g_control = model$control,
    locus_het = locus_het,
    pi = pi
  )
  scores <- .cc_ngs_scores_from_moi(MOI)
  N_ctrl <- k * N_case
  ngs <- .cc_ngs_ahn_ncp(
    g_case = heterogeneity$g_case_after_locus_het,
    g_control = heterogeneity$g_ctrl_after_locus_het,
    N_case = N_case,
    N_control = N_ctrl,
    coverage = coverage,
    seq_error = seq_error,
    scores = scores
  )

  power <- .cc_ngs_chisq_power(ngs$lambda, alpha)

  out <- list(
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    k = k,
    power = as.numeric(power),
    lambda = ngs$lambda,
    MOI = MOI,
    scores = scores,
    coverage = coverage,
    seq_error = seq_error,
    locus_het = heterogeneity,
    model_info = list(
      input_mode = "model_based",
      prev = prev,
      pd = pd,
      qd = 1 - pd,
      R1 = model$R1,
      R2 = model$R2,
      MOI = model$MOI,
      penetrances = model$penetrances
    ),
    freqs = list(
      population = model$population,
      case_true_pre_heterogeneity = model$case,
      control_true_pre_heterogeneity = model$control,
      case_true = ngs$case_true,
      control_true = ngs$control_true,
      case_called = ngs$case_called,
      control_called = ngs$control_called
    ),
    transition_matrix = ngs$E
  )
  class(out) <- "cc_ngs_power"

  if (isTRUE(verbose)) {
    print(out)
  }

  invisible(out)
}

#' @export
print.cc_ngs_power <- function(x, ...) {
  cat("Case-control sequencing trend-test power\n")
  cat(sprintf("Cases: %s; controls: %s; total: %s\n",
              formatC(x$N_case, format = "f", digits = 0, big.mark = ","),
              formatC(x$N_ctrl, format = "f", digits = 0, big.mark = ","),
              formatC(x$N_total, format = "f", digits = 0, big.mark = ",")))
  cat(sprintf("Coverage: %s; sequencing error: %.4g\n",
              formatC(x$coverage, format = "f", digits = 0), x$seq_error))
  if (isTRUE(x$locus_het$enabled) && x$locus_het$pi < 1) {
    cat(sprintf("Locus heterogeneity: %.1f%% (pi = %.4g)\n",
                100 * (1 - x$locus_het$pi), x$locus_het$pi))
  }
  cat(sprintf("MOI: %s; alpha: %.4g\n", x$MOI, x$alpha))
  cat(sprintf("NCP: %.4f; power: %.1f%%\n", x$lambda, 100 * x$power))
  invisible(x)
}

#' Analytic MSSN for a Case-Control Sequencing Study
#'
#' Computes the minimum sample size necessary (MSSN) for a model-based
#' case-control sequencing trend design. It uses the same fixed-depth,
#' symmetric sequencing-error model and deterministic maximum-likelihood
#' genotype calls as \code{\link{cc_ngs_power}}.
#'
#' @param power Numeric in \eqn{(0,1)}. Requested power.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk.
#' @param coverage Positive integer sequencing depth.
#' @param seq_error Symmetric per-read sequencing-error probability in
#'   \eqn{[0,0.5)}.
#' @param MOI Character mode of inheritance: \code{"M"} for multiplicative,
#'   \code{"D"} for dominant, or \code{"Rec"} for recessive.
#' @param k Numeric \eqn{> 0}. Planned control-to-case sample-size ratio.
#' @param verbose Logical. If \code{TRUE}, print a concise result summary.
#' @param locus_het Logical. If \code{TRUE}, apply the canonical PAWEH
#'   case-control locus-heterogeneity mixture before sequencing observation.
#' @param pi Numeric in \eqn{[0,1]}. Locus-homogeneity fraction used when
#'   \code{locus_het = TRUE}. One retains the original associated-case
#'   distribution; zero makes the case distribution equal to controls. When
#'   \code{locus_het = FALSE}, \code{pi} must remain at its default value of 1.
#'
#' @details
#' Locus heterogeneity is applied to true case genotype probabilities as
#' \eqn{g_{case,H}=\pi g_{case}+(1-\pi)g_{control}}, using the same
#' parameterization as ordinary PAWEH case-control design. Sequencing
#' observation follows this mixture. Under the current nondifferential error
#' model, applying the common transition matrix before or after forming the
#' mixture is algebraically equivalent; this need not hold for future
#' differential case/control sequencing error. When \eqn{\pi=0}, no finite
#' MSSN exists for target power greater than \code{alpha} because the trend
#' contrast is zero.
#'
#' The function numerically inverts the one-degree-of-freedom noncentral
#' chi-square distribution only to obtain the target NCP. It then solves the
#' Ahn/Chapman-Nam trend-test sample-size equation analytically. Planned cases
#' are the ceiling of the continuous requirement and planned controls are
#' \code{ceiling(k * MSSN_case)}, following the PAWEH convention. Achieved NCP
#' and power are recomputed using these actual integer sample sizes, with a
#' local boundary adjustment if required to ensure target attainment and
#' integer minimality.
#'
#' Scores are \code{c(0,1,2)} for \code{"M"}, \code{c(0,1,1)} for
#' \code{"D"}, and \code{c(0,0,1)} for \code{"Rec"}. Finite depth can cause
#' genotype-call uncertainty even when \code{seq_error = 0}, because a true
#' heterozygote can yield reads from only one allele.
#'
#' This is an analytic study-design calculation, not LTTae,NGS, a raw-read
#' latent-genotype likelihood or EM method, downstream association testing, or
#' simulation.
#'
#' @return Invisibly, an object of class \code{"cc_ngs_mssn"} containing the
#'   target and achieved power and NCP, continuous and integer sample sizes,
#'   model and sequencing inputs, true and called genotype frequencies, and
#'   the true-to-called transition matrix.
#'
#' @references
#' Ahn, K., Haynes, C., Kim, W., St. Fleur, R., Gordon, D., & Finch, S. J.
#' (2007). The effects of SNP genotyping errors on the power of the
#' Cochran-Armitage linear trend test for case/control association studies.
#' \emph{Annals of Human Genetics}, 71, 249--261.
#' \doi{10.1111/j.1469-1809.2006.00318.x}.
#'
#' Gordon, D., Finch, S. J., & Kim, W. (2020). \emph{Heterogeneity in
#' Statistical Genetics}. Springer. \doi{10.1007/978-3-030-61121-7}.
#'
#' @examples
#' cc_ngs_mssn(
#'   power = 0.80, alpha = 0.05,
#'   prev = 0.05, pd = 0.30, R2 = 1.8,
#'   coverage = 20, seq_error = 0.01,
#'   MOI = "M", verbose = FALSE
#' )
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export
cc_ngs_mssn <- function(
    power, alpha,
    prev, pd, R2,
    coverage, seq_error,
    MOI = c("M", "D", "Rec"),
    k = 1,
    verbose = TRUE,
    locus_het = FALSE,
    pi = 1
) {
  MOI <- match.arg(MOI)

  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0 || power >= 1) {
    stop("power must be a single finite number in (0, 1).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number in (0, 1).")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("k must be a single finite positive number: N_ctrl / N_case.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }

  model <- .cc_model_genotype_frequencies(
    pd = pd,
    R2 = R2,
    MOI = MOI,
    prev = prev
  )
  heterogeneity <- .cc_ngs_apply_locus_heterogeneity(
    g_case = model$case,
    g_control = model$control,
    locus_het = locus_het,
    pi = pi
  )
  scores <- .cc_ngs_scores_from_moi(MOI)
  called <- .cc_ngs_called_frequencies(
    g_case = heterogeneity$g_case_after_locus_het,
    g_control = heterogeneity$g_ctrl_after_locus_het,
    coverage = coverage,
    seq_error = seq_error
  )
  lambda_target <- .cc_ngs_target_ncp(power, alpha)
  components <- .cc_ngs_mssn_components(
    g_case = called$case_called,
    g_control = called$control_called,
    k = k,
    scores = scores,
    lambda_target = lambda_target
  )

  initial_case <- max(1, ceiling(components$N_case_continuous))
  evaluate_design <- function(N_case) {
    N_control <- ceiling(k * N_case)
    lambda <- .cc_ahn_trend_ncp(
      g_case = called$case_called,
      g_control = called$control_called,
      N_case = N_case,
      N_control = N_control,
      scores = scores
    )
    list(
      N_case = N_case,
      N_control = N_control,
      lambda = lambda,
      power = .cc_ngs_chisq_power(lambda, alpha)
    )
  }

  planned <- evaluate_design(initial_case)
  tolerance <- 1e-12
  while (planned$N_case > 1) {
    previous <- evaluate_design(planned$N_case - 1)
    if (previous$power < power - tolerance) {
      break
    }
    planned <- previous
  }
  while (planned$power < power - tolerance) {
    planned <- evaluate_design(planned$N_case + 1)
  }

  out <- list(
    power_target = power,
    alpha = alpha,
    MSSN_case = planned$N_case,
    MSSN_ctrl = planned$N_control,
    MSSN_total = planned$N_case + planned$N_control,
    N_case_continuous = components$N_case_continuous,
    achieved_power = planned$power,
    achieved_lambda = planned$lambda,
    lambda_target = lambda_target,
    initial_MSSN_case = initial_case,
    rounding_adjustment = planned$N_case - initial_case,
    k = k,
    MOI = MOI,
    scores = scores,
    coverage = coverage,
    seq_error = seq_error,
    locus_het = heterogeneity,
    model_info = list(
      input_mode = "model_based",
      prev = prev,
      pd = pd,
      qd = 1 - pd,
      R1 = model$R1,
      R2 = model$R2,
      MOI = model$MOI,
      penetrances = model$penetrances
    ),
    freqs = list(
      population = model$population,
      case_true_pre_heterogeneity = model$case,
      control_true_pre_heterogeneity = model$control,
      case_true = called$case_true,
      control_true = called$control_true,
      case_called = called$case_called,
      control_called = called$control_called
    ),
    transition_matrix = called$E
  )
  class(out) <- "cc_ngs_mssn"

  if (isTRUE(verbose)) {
    print(out)
  }

  invisible(out)
}

#' @export
print.cc_ngs_mssn <- function(x, ...) {
  cat("Case-control sequencing trend-test MSSN\n")
  cat(sprintf("Target power: %.1f%%; alpha: %.4g\n",
              100 * x$power_target, x$alpha))
  cat(sprintf("Cases: %s; controls: %s; total MSSN: %s\n",
              formatC(x$MSSN_case, format = "f", digits = 0, big.mark = ","),
              formatC(x$MSSN_ctrl, format = "f", digits = 0, big.mark = ","),
              formatC(x$MSSN_total, format = "f", digits = 0, big.mark = ",")))
  cat(sprintf("Coverage: %s; sequencing error: %.4g; MOI: %s\n",
              formatC(x$coverage, format = "f", digits = 0),
              x$seq_error, x$MOI))
  if (isTRUE(x$locus_het$enabled) && x$locus_het$pi < 1) {
    cat(sprintf("Locus heterogeneity: %.1f%% (pi = %.4g)\n",
                100 * (1 - x$locus_het$pi), x$locus_het$pi))
  }
  cat(sprintf("Achieved power: %.1f%%\n", 100 * x$achieved_power))
  invisible(x)
}
