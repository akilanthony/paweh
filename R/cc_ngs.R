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

.cc_ngs_apply_pheno_misclassification <- function(
    g_case,
    g_control,
    prev,
    pheno_misclass = FALSE,
    theta = 0,
    phi = 0
) {
  .cc_ngs_validate_genotype_frequencies(g_case, "g_case")
  .cc_ngs_validate_genotype_frequencies(g_control, "g_control")
  if (!is.logical(pheno_misclass) || length(pheno_misclass) != 1L ||
      is.na(pheno_misclass)) {
    stop("pheno_misclass must be TRUE or FALSE.")
  }
  if (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta) ||
      theta < 0 || theta >= 1) {
    stop("theta must be a single number in [0,1).")
  }
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) ||
      phi < 0 || phi >= 1) {
    stop("phi must be a single number in [0,1).")
  }
  if (!is.numeric(prev) || length(prev) != 1L || !is.finite(prev) ||
      prev <= 0 || prev >= 1) {
    stop("When pheno_misclass=TRUE, prev must be a single number in (0,1).")
  }

  transformed <- if (isTRUE(pheno_misclass)) {
    .cc_apply_pheno_misclass(
      g_aff = g_case,
      g_unaff = g_control,
      prev = prev,
      theta = theta,
      phi = phi
    )
  } else {
    list(
      g_case_obs = as.numeric(g_case),
      g_ctrl_obs = as.numeric(g_control),
      case_denom = NA_real_,
      ctrl_denom = NA_real_
    )
  }

  list(
    enabled = pheno_misclass,
    theta = theta,
    phi = phi,
    g_case_before_pheno_misclass = as.numeric(g_case),
    g_ctrl_before_pheno_misclass = as.numeric(g_control),
    g_case_after_pheno_misclass = transformed$g_case_obs,
    g_ctrl_after_pheno_misclass = transformed$g_ctrl_obs,
    case_denom = transformed$case_denom,
    ctrl_denom = transformed$ctrl_denom
  )
}

.cc_ngs_apply_genotype_misclassification <- function(
    g_case, g_control,
    geno_misclass = c("none", "1p", "2p", "3p", "diff3p"),
    e = 0, e1 = 0, e2 = 0,
    e01 = 0, e02 = 0, e03 = 0,
    case_e01 = 0, case_e02 = 0, case_e03 = 0,
    ctrl_e01 = 0, ctrl_e02 = 0, ctrl_e03 = 0,
    diff_source = c("explicit", "case", "ctrl"),
    diff_multiplier = 1
) {
  genotype <- .cc_resolve_genotype_misclassification(
    geno_misclass = geno_misclass,
    e = e, e1 = e1, e2 = e2,
    e01 = e01, e02 = e02, e03 = e03,
    case_e01 = case_e01, case_e02 = case_e02, case_e03 = case_e03,
    ctrl_e01 = ctrl_e01, ctrl_e02 = ctrl_e02, ctrl_e03 = ctrl_e03,
    diff_source = diff_source,
    diff_multiplier = diff_multiplier
  )
  if (identical(genotype$info$model, "none")) {
    info <- genotype$info
    info$M_case <- genotype$M_case
    info$M_ctrl <- genotype$M_ctrl
    return(list(
      case_final = as.numeric(g_case),
      control_final = as.numeric(g_control),
      info = info
    ))
  }
  case_final <- .cc_apply_genotype_misclass(g_case, genotype$M_case)
  case_final <- case_final / sum(case_final)
  control_final <- .cc_apply_genotype_misclass(g_control, genotype$M_ctrl)
  control_final <- control_final / sum(control_final)
  .cc_ngs_validate_genotype_frequencies(case_final, "case_final")
  .cc_ngs_validate_genotype_frequencies(control_final, "control_final")

  info <- genotype$info
  if (is.null(info$M_case)) info$M_case <- genotype$M_case
  if (is.null(info$M_ctrl)) info$M_ctrl <- genotype$M_ctrl

  list(
    case_final = case_final,
    control_final = control_final,
    info = info
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
#' genotype probabilities, applies optional locus and phenotype modifiers,
#' observes each group through its fixed-depth symmetric sequencing-error model
#' with deterministic maximum-likelihood genotype calls, applies optional
#' genotype misclassification, and evaluates the Ahn/Chapman-Nam
#' Cochran-Armitage trend-test noncentrality parameter on the final observed
#' genotype probabilities.
#'
#' @param N_case Numeric \eqn{> 0}. Number of cases.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk.
#' @param coverage Common positive integer fixed sequencing depth. May be omitted
#'   when both group-specific depths are supplied.
#' @param seq_error Symmetric per-read sequencing-error probability in
#'   \eqn{[0,0.5)}. Common shorthand; may be omitted when both group-specific
#'   errors are supplied.
#' @param case_coverage,ctrl_coverage Case/control fixed depths, respectively.
#'   Each defaults to \code{coverage}; an explicit value overrides it.
#' @param case_seq_error,ctrl_seq_error Case/control per-read error probabilities.
#'   Each defaults to \code{seq_error}; an explicit value overrides it.
#'   Effective depths and errors obey the common parameter constraints.
#' @param pheno_misclass Logical. If \code{TRUE}, apply the ordinary PAWEH
#'   case-control phenotype-misclassification model before sequencing.
#' @param theta Numeric in \eqn{[0,1)}. Probability that a truly affected
#'   individual is classified as a control.
#' @param phi Numeric in \eqn{[0,1)}. Probability that a truly unaffected
#'   individual is classified as a case.
#' @param geno_misclass Character genotype-misclassification model:
#'   \code{"none"}, \code{"1p"}, \code{"2p"}, \code{"3p"}, or
#'   \code{"diff3p"}. The selected model is applied after sequencing calls.
#' @param e Numeric in \eqn{[0,0.5]}. One-parameter symmetric error rate.
#' @param e1,e2 Numeric two-parameter homozygote-to-heterozygote and
#'   heterozygote-to-homozygote error rates.
#' @param e01,e02,e03 Numeric non-differential three-parameter error rates.
#' @param case_e01,case_e02,case_e03 Case-specific three-parameter error rates.
#' @param ctrl_e01,ctrl_e02,ctrl_e03 Control-specific three-parameter error
#'   rates.
#' @param diff_source Character source for differential three-parameter errors:
#'   \code{"explicit"}, \code{"case"}, or \code{"ctrl"}.
#' @param diff_multiplier Nonnegative multiplier used to derive the other
#'   group's rates when \code{diff_source} is \code{"case"} or \code{"ctrl"}.
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
#' with the control distribution unchanged. The full observation order is:
#' baseline genotype model, optional locus heterogeneity, optional phenotype
#' misclassification, group-specific sequencing, optional genotype
#' misclassification, then the trend-test NCP.
#' Phenotype misclassification uses the same prevalence-weighted transformation
#' as \code{\link{cc_power}}, with
#' \code{theta = Pr(affected -> control)} and
#' \code{phi = Pr(unaffected -> case)}. The resulting observed-case and
#' observed-control distributions enter sequencing. Separate row-true,
#' column-called matrices are then constructed using each observed group's
#' effective depth and error, and applied as `t(M_case) %*% g_case` and
#' `t(M_ctrl) %*% g_control`. The ordinary PAWEH 1p, 2p, 3p, or differential
#' 3p genotype-error matrices are then applied to those sequencing-called
#' frequencies with the same row-current, column-observed convention.
#' Phenotype, sequencing, and genotype errors may be active simultaneously.
#' Only equal sequencing mechanisms permit commuting
#' a shared biological mixture and sequencing. At \eqn{\pi=0}, equal mechanisms
#' with no phenotype-induced difference give zero NCP and power equal to
#' \code{alpha}.
#'
#' Differential sequencing or genotype-misclassification mechanisms can create
#' observed case/control differences even under a biological null. Results
#' remain nominal asymptotic calculations using the existing chi-square
#' critical value; the null distribution and Type I error are not separately
#' recalibrated here.
#'
#' This is an analytic study-design calculation applied to sequencing-derived
#' called genotypes. It is not a raw-read likelihood or EM analysis, performs
#' no downstream association testing, and uses no simulation. Even when
#' \code{seq_error = 0}, finite depth can cause call uncertainty because a true
#' heterozygote can yield reads from only one allele.
#'
#' @return Invisibly, an object of class \code{"cc_ngs_power"} containing the
#'   design inputs, trend scores, NCP and power, model information, staged and
#'   called genotype frequencies, phenotype- and genotype-error metadata, and
#'   the true-to-called transition matrix. The frequency list distinguishes
#'   post-sequencing calls from final post-genotype-misclassification values.
#'   The \code{sequencing} list records the four effective group parameters and
#'   \code{case_transition_matrix}/\code{ctrl_transition_matrix}. Legacy
#'   \code{coverage}/\code{seq_error} retain supplied common inputs (or NULL);
#'   \code{transition_matrix} aliases the case matrix, shared when mechanisms
#'   are equal. Use the two explicit matrices for differential designs.
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
    coverage = NULL, seq_error = NULL,
    MOI = c("M", "D", "Rec"),
    k = 1,
    verbose = TRUE,
    locus_het = FALSE,
    pi = 1,
    case_coverage = coverage,
    ctrl_coverage = coverage,
    case_seq_error = seq_error,
    ctrl_seq_error = seq_error,
    pheno_misclass = FALSE,
    theta = 0,
    phi = 0,
    geno_misclass = c("none", "1p", "2p", "3p", "diff3p"),
    e = 0,
    e1 = 0,
    e2 = 0,
    e01 = 0,
    e02 = 0,
    e03 = 0,
    case_e01 = 0,
    case_e02 = 0,
    case_e03 = 0,
    ctrl_e01 = 0,
    ctrl_e02 = 0,
    ctrl_e03 = 0,
    diff_source = c("explicit", "case", "ctrl"),
    diff_multiplier = 1
) {
  MOI <- match.arg(MOI)
  geno_misclass <- match.arg(geno_misclass)
  diff_source <- match.arg(diff_source)

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
  phenotype <- .cc_ngs_apply_pheno_misclassification(
    g_case = heterogeneity$g_case_after_locus_het,
    g_control = heterogeneity$g_ctrl_after_locus_het,
    prev = prev,
    pheno_misclass = pheno_misclass,
    theta = theta,
    phi = phi
  )
  scores <- .cc_ngs_scores_from_moi(MOI)
  N_ctrl <- k * N_case
  called <- .cc_ngs_called_frequencies(
    g_case = phenotype$g_case_after_pheno_misclass,
    g_control = phenotype$g_ctrl_after_pheno_misclass,
    coverage = coverage,
    seq_error = seq_error,
    case_coverage = case_coverage, ctrl_coverage = ctrl_coverage,
    case_seq_error = case_seq_error, ctrl_seq_error = ctrl_seq_error
  )
  genotype <- .cc_ngs_apply_genotype_misclassification(
    g_case = called$case_called,
    g_control = called$control_called,
    geno_misclass = geno_misclass,
    e = e, e1 = e1, e2 = e2,
    e01 = e01, e02 = e02, e03 = e03,
    case_e01 = case_e01, case_e02 = case_e02, case_e03 = case_e03,
    ctrl_e01 = ctrl_e01, ctrl_e02 = ctrl_e02, ctrl_e03 = ctrl_e03,
    diff_source = diff_source,
    diff_multiplier = diff_multiplier
  )
  lambda <- .cc_ahn_trend_ncp(
    g_case = genotype$case_final,
    g_control = genotype$control_final,
    N_case = N_case,
    N_control = N_ctrl,
    scores = scores
  )

  power <- .cc_ngs_chisq_power(lambda, alpha)

  out <- list(
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    k = k,
    power = as.numeric(power),
    lambda = lambda,
    MOI = MOI,
    scores = scores,
    coverage = coverage,
    seq_error = seq_error,
    locus_het = heterogeneity,
    errors = list(
      phenotype_misclass = phenotype,
      genotype_misclass = genotype$info
    ),
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
      case_post_heterogeneity = heterogeneity$g_case_after_locus_het,
      control_post_heterogeneity = heterogeneity$g_ctrl_after_locus_het,
      case_preseq = phenotype$g_case_after_pheno_misclass,
      control_preseq = phenotype$g_ctrl_after_pheno_misclass,
      case_true = called$case_true,
      control_true = called$control_true,
      case_post_sequencing = called$case_called,
      control_post_sequencing = called$control_called,
      case_final = genotype$case_final,
      control_final = genotype$control_final,
      case_called = genotype$case_final,
      control_called = genotype$control_final
    ),
    transition_matrix = called$E,
    sequencing = called$sequencing
  )
  class(out) <- "cc_ngs_power"

  if (isTRUE(verbose)) {
    .paweh_print_cc_ngs_power(out)
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
  .cc_ngs_print_sequencing(x)
  if (isTRUE(x$locus_het$enabled) && x$locus_het$pi < 1) {
    cat(sprintf("Locus heterogeneity: %.1f%% (pi = %.4g)\n",
                100 * (1 - x$locus_het$pi), x$locus_het$pi))
  }
  pheno <- x$errors$phenotype_misclass
  if (!is.null(pheno) && isTRUE(pheno$enabled)) {
    cat(sprintf("Phenotype misclassification: theta = %.4g; phi = %.4g\n",
                pheno$theta, pheno$phi))
  }
  genotype <- x$errors$genotype_misclass
  if (!is.null(genotype) && isTRUE(genotype$enabled)) {
    cat(sprintf("Genotype misclassification: %s%s\n", genotype$model,
                if (identical(genotype$model, "diff3p_homhet_homhom"))
                  " (differential)" else ""))
  }
  cat(sprintf("MOI: %s; alpha: %.4g\n", x$MOI, x$alpha))
  cat(sprintf("NCP: %.4f; power: %.1f%%\n", x$lambda, 100 * x$power))
  invisible(x)
}

#' Analytic MSSN for a Case-Control Sequencing Study
#'
#' Computes the minimum sample size necessary (MSSN) for a model-based
#' case-control sequencing trend design. It uses the same phenotype modifier,
#' group-specific fixed-depth symmetric sequencing-error models, deterministic
#' maximum-likelihood genotype calls, and genotype-misclassification stage as
#' \code{\link{cc_ngs_power}}.
#'
#' @param power Numeric in \eqn{(0,1)}. Requested power.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk.
#' @param coverage Common positive integer fixed sequencing depth. May be omitted
#'   when both group-specific depths are supplied.
#' @param seq_error Symmetric per-read sequencing-error probability in
#'   \eqn{[0,0.5)}. Common shorthand; may be omitted when both group-specific
#'   errors are supplied.
#' @param case_coverage,ctrl_coverage Case/control fixed depths, respectively.
#'   Each defaults to \code{coverage}; an explicit value overrides it.
#' @param case_seq_error,ctrl_seq_error Case/control per-read error probabilities.
#'   Each defaults to \code{seq_error}; an explicit value overrides it.
#'   Effective depths and errors obey the common parameter constraints.
#' @param pheno_misclass Logical. If \code{TRUE}, apply the ordinary PAWEH
#'   case-control phenotype-misclassification model before sequencing.
#' @param theta Numeric in \eqn{[0,1)}. Probability that a truly affected
#'   individual is classified as a control.
#' @param phi Numeric in \eqn{[0,1)}. Probability that a truly unaffected
#'   individual is classified as a case.
#' @param geno_misclass Character genotype-misclassification model:
#'   \code{"none"}, \code{"1p"}, \code{"2p"}, \code{"3p"}, or
#'   \code{"diff3p"}. The selected model is applied after sequencing calls.
#' @param e Numeric in \eqn{[0,0.5]}. One-parameter symmetric error rate.
#' @param e1,e2 Numeric two-parameter homozygote-to-heterozygote and
#'   heterozygote-to-homozygote error rates.
#' @param e01,e02,e03 Numeric non-differential three-parameter error rates.
#' @param case_e01,case_e02,case_e03 Case-specific three-parameter error rates.
#' @param ctrl_e01,ctrl_e02,ctrl_e03 Control-specific three-parameter error
#'   rates.
#' @param diff_source Character source for differential three-parameter errors:
#'   \code{"explicit"}, \code{"case"}, or \code{"ctrl"}.
#' @param diff_multiplier Nonnegative multiplier used to derive the other
#'   group's rates when \code{diff_source} is \code{"case"} or \code{"ctrl"}.
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
#' parameterization as ordinary PAWEH case-control design. The optional
#' prevalence-weighted phenotype transformation from \code{\link{cc_mssn}}
#' follows locus heterogeneity and precedes sequencing. Case sequencing settings
#' apply to individuals observed as cases after phenotype classification, and
#' control settings apply to individuals observed as controls. Genotype
#' misclassification is applied to the resulting sequencing calls. See
#' \code{\link{cc_ngs_power}} for the complete order and matrix convention.
#' When \eqn{\pi=0} and sequencing mechanisms are equal, no finite MSSN exists
#' for target power greater than \code{alpha} if the observed contrast is zero.
#' Differential mechanisms may create an observed contrast even under a
#' biological null. MSSN is a nominal asymptotic calculation using the existing
#' chi-square critical value, without separate null-distribution or Type I error
#' recalibration.
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
#'   the true-to-called transition matrix. See \code{\link{cc_ngs_power}}
#'   for sequencing metadata and legacy-field conventions.
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
    coverage = NULL, seq_error = NULL,
    MOI = c("M", "D", "Rec"),
    k = 1,
    verbose = TRUE,
    locus_het = FALSE,
    pi = 1,
    case_coverage = coverage,
    ctrl_coverage = coverage,
    case_seq_error = seq_error,
    ctrl_seq_error = seq_error,
    pheno_misclass = FALSE,
    theta = 0,
    phi = 0,
    geno_misclass = c("none", "1p", "2p", "3p", "diff3p"),
    e = 0,
    e1 = 0,
    e2 = 0,
    e01 = 0,
    e02 = 0,
    e03 = 0,
    case_e01 = 0,
    case_e02 = 0,
    case_e03 = 0,
    ctrl_e01 = 0,
    ctrl_e02 = 0,
    ctrl_e03 = 0,
    diff_source = c("explicit", "case", "ctrl"),
    diff_multiplier = 1
) {
  MOI <- match.arg(MOI)
  geno_misclass <- match.arg(geno_misclass)
  diff_source <- match.arg(diff_source)

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
  phenotype <- .cc_ngs_apply_pheno_misclassification(
    g_case = heterogeneity$g_case_after_locus_het,
    g_control = heterogeneity$g_ctrl_after_locus_het,
    prev = prev,
    pheno_misclass = pheno_misclass,
    theta = theta,
    phi = phi
  )
  scores <- .cc_ngs_scores_from_moi(MOI)
  called <- .cc_ngs_called_frequencies(
    g_case = phenotype$g_case_after_pheno_misclass,
    g_control = phenotype$g_ctrl_after_pheno_misclass,
    coverage = coverage,
    seq_error = seq_error,
    case_coverage = case_coverage, ctrl_coverage = ctrl_coverage,
    case_seq_error = case_seq_error, ctrl_seq_error = ctrl_seq_error
  )
  genotype <- .cc_ngs_apply_genotype_misclassification(
    g_case = called$case_called,
    g_control = called$control_called,
    geno_misclass = geno_misclass,
    e = e, e1 = e1, e2 = e2,
    e01 = e01, e02 = e02, e03 = e03,
    case_e01 = case_e01, case_e02 = case_e02, case_e03 = case_e03,
    ctrl_e01 = ctrl_e01, ctrl_e02 = ctrl_e02, ctrl_e03 = ctrl_e03,
    diff_source = diff_source,
    diff_multiplier = diff_multiplier
  )
  lambda_target <- .cc_ngs_target_ncp(power, alpha)
  components <- .cc_ngs_mssn_components(
    g_case = genotype$case_final,
    g_control = genotype$control_final,
    k = k,
    scores = scores,
    lambda_target = lambda_target
  )

  initial_case <- max(1, ceiling(components$N_case_continuous))
  evaluate_design <- function(N_case) {
    N_control <- ceiling(k * N_case)
    lambda <- .cc_ahn_trend_ncp(
      g_case = genotype$case_final,
      g_control = genotype$control_final,
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
    errors = list(
      phenotype_misclass = phenotype,
      genotype_misclass = genotype$info
    ),
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
      case_post_heterogeneity = heterogeneity$g_case_after_locus_het,
      control_post_heterogeneity = heterogeneity$g_ctrl_after_locus_het,
      case_preseq = phenotype$g_case_after_pheno_misclass,
      control_preseq = phenotype$g_ctrl_after_pheno_misclass,
      case_true = called$case_true,
      control_true = called$control_true,
      case_post_sequencing = called$case_called,
      control_post_sequencing = called$control_called,
      case_final = genotype$case_final,
      control_final = genotype$control_final,
      case_called = genotype$case_final,
      control_called = genotype$control_final
    ),
    transition_matrix = called$E,
    sequencing = called$sequencing
  )
  class(out) <- "cc_ngs_mssn"

  if (isTRUE(verbose)) {
    .paweh_print_cc_ngs_mssn(out)
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
  .cc_ngs_print_sequencing(x)
  cat(sprintf("MOI: %s\n", x$MOI))
  if (isTRUE(x$locus_het$enabled) && x$locus_het$pi < 1) {
    cat(sprintf("Locus heterogeneity: %.1f%% (pi = %.4g)\n",
                100 * (1 - x$locus_het$pi), x$locus_het$pi))
  }
  pheno <- x$errors$phenotype_misclass
  if (!is.null(pheno) && isTRUE(pheno$enabled)) {
    cat(sprintf("Phenotype misclassification: theta = %.4g; phi = %.4g\n",
                pheno$theta, pheno$phi))
  }
  genotype <- x$errors$genotype_misclass
  if (!is.null(genotype) && isTRUE(genotype$enabled)) {
    cat(sprintf("Genotype misclassification: %s%s\n", genotype$model,
                if (identical(genotype$model, "diff3p_homhet_homhom"))
                  " (differential)" else ""))
  }
  cat(sprintf("Achieved power: %.1f%%\n", 100 * x$achieved_power))
  invisible(x)
}

.cc_ngs_print_sequencing <- function(x) {
  z <- x$sequencing
  if (is.null(z)) {
    z <- list(case_coverage = x$coverage, ctrl_coverage = x$coverage,
              case_seq_error = x$seq_error, ctrl_seq_error = x$seq_error)
  }
  if (z$case_coverage == z$ctrl_coverage &&
      z$case_seq_error == z$ctrl_seq_error) {
    cat(sprintf("Coverage: %s; sequencing error: %.4g\n",
                z$case_coverage, z$case_seq_error))
  } else {
    cat(sprintf("Cases coverage: %s; sequencing error: %.4g\n",
                z$case_coverage, z$case_seq_error))
    cat(sprintf("Controls coverage: %s; sequencing error: %.4g\n",
                z$ctrl_coverage, z$ctrl_seq_error))
  }
}
