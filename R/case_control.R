#' Case-Control Minimum Sample Size for Conditional Genotype Frequencies
#'
#' Computes the minimum sample size necessary (MSSN) for case-control association
#' tests using conditional genotype frequencies. The function supports
#' model-based genotype frequencies, model-free genotype frequencies, optional
#' locus heterogeneity, optional genotype misclassification, genotype
#' chi-square tests, optional allelic chi-square tests, and genotype trend tests.
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
#'   case genotype frequencies before misclassification.
#' @param pi Numeric in \eqn{[0,1]}. Locus-homogeneity fraction used when
#'   \code{locus_het = TRUE}.
#' @param k Numeric \eqn{> 0}. Control-to-case sample size ratio
#'   \eqn{N_{ctrl} / N_{case}}.
#' @param w Numeric vector of length 3. Genotype trend-test scores. The three
#'   weights cannot all be equal.
#' @param include_allelic Logical. If \code{TRUE}, includes the allelic
#'   chi-square test in the returned object and verbose output.
#' @param geno_misclass Character. Genotype misclassification model:
#'   \code{"none"}, \code{"1p"}, \code{"2p"}, \code{"3p"}, or \code{"diff3p"}.
#' @param e Numeric. Symmetric one-parameter misclassification rate for
#'   \code{geno_misclass = "1p"}.
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
#' \item Optionally apply genotype misclassification matrices to the true case
#' and control genotype frequencies.
#' \item Compute genotype chi-square, optional allelic chi-square, and genotype
#' trend-test MSSN values from the observed genotype frequencies.
#' }
#'
#' With \code{input_mode = "model_based"}, conditional case and control genotype
#' frequencies are derived from \code{prev}, \code{pd}, \code{R2}, and
#' \code{MOI}. With \code{input_mode = "model_free"}, the user supplies
#' \code{g1} and \code{g0} directly.
#'
#' The genotype misclassification models are:
#' \code{"none"} for identity matrices, \code{"1p"} for one symmetric error
#' rate, \code{"2p"} for adjacent homozygote/heterozygote and heterozygote
#' error rates, \code{"3p"} for a non-differential three-parameter matrix, and
#' \code{"diff3p"} for separate case and control three-parameter matrices.
#'
#' For \code{geno_misclass = "diff3p"}, \code{diff_source = "explicit"} uses
#' the case and control error parameters exactly as supplied. With
#' \code{diff_source = "case"}, control parameters are computed by multiplying
#' the case parameters by \code{diff_multiplier}. With
#' \code{diff_source = "ctrl"}, case parameters are computed by multiplying the
#' control parameters by \code{diff_multiplier}.
#'
#' If \code{include_allelic = FALSE}, the allelic test is omitted and
#' \code{tests$alleles} is \code{NULL}. Internal effect-size components
#' \code{S} are retained in the returned object for validation and debugging,
#' but are not printed in the clean verbose output.
#'
#' @return An object of class \code{"cc_mssn_conditional_full"}: a nested list
#' with components \code{alpha}, \code{target_power}, \code{input_mode},
#' \code{k}, \code{w}, \code{include_allelic}, \code{locus_het},
#' \code{errors}, \code{model_info}, \code{tests}, and \code{freqs}.
#' \code{tests$genotypes} and \code{tests$trend} contain test labels, degrees
#' of freedom, target non-centrality parameters, internal \code{S} components,
#' and MSSN case/control/total values. \code{tests$alleles} has the same MSSN
#' fields when requested, otherwise it is \code{NULL}. \code{freqs} stores
#' baseline, true, and observed case/control genotype frequencies and observed
#' allele frequencies.
#'
#' @examples
#' cc_mssn_conditional_full(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   geno_misclass = "none",
#'   verbose = FALSE
#' )
#'
#' cc_mssn_conditional_full(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   verbose = FALSE
#' )
#'
#' cc_mssn_conditional_full(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' cc_mssn_conditional_full(
#'   power = 0.8, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   include_allelic = FALSE,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' TODO: Provide exact textbook equation numbers/pages for the case-control
#' genotype, allelic, trend, locus-heterogeneity, and genotype-misclassification
#' calculations.
#'
#' @importFrom stats pchisq qchisq uniroot
#' @export
cc_mssn_conditional_full <- function(
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
    
    k = 1,
    w = c(0, 1, 2),
    include_allelic = TRUE,
    
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
  
  cc_geno_to_allele_freqs <- function(gj) {
    check_genotype_freqs(gj, "gj")
    p <- gj[3] + 0.5 * gj[2]
    c(q = 1 - p, p = p)
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
  if (!is.logical(include_allelic) || length(include_allelic) != 1)
    stop("include_allelic must be TRUE or FALSE.")
  if (!is.logical(locus_het) || length(locus_het) != 1)
    stop("locus_het must be TRUE or FALSE.")
  if (!is.numeric(pi) || length(pi) != 1 || pi < 0 || pi > 1)
    stop("pi must be a single number in [0,1].")
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
  
  p1_obs <- cc_geno_to_allele_freqs(g1_obs)
  p0_obs <- cc_geno_to_allele_freqs(g0_obs)
  
  # ---- target lambdas ----
  lambda_star_g <- chisq_ncp_target(power = power, alpha = alpha, df = 2)
  lambda_star_1 <- chisq_ncp_target(power = power, alpha = alpha, df = 1)
  
  # ---- genotype chi-square ----
  S_g <- sum((g1_obs - g0_obs)^2 / (g1_obs + k * g0_obs))
  
  if (!is.finite(S_g) || S_g <= 0)
    stop("Genotype S <= 0; check inputs.")
  
  MSSN_case_g <- ceiling(lambda_star_g / (k * S_g))
  MSSN_ctrl_g <- ceiling(k * MSSN_case_g)
  
  # ---- allelic chi-square ----
  allelic_result <- NULL
  
  if (isTRUE(include_allelic)) {
    
    S_a <- sum((p1_obs - p0_obs)^2 / (p1_obs + k * p0_obs))
    
    if (!is.finite(S_a) || S_a <= 0)
      stop("Allele S <= 0; check inputs.")
    
    MSSN_case_a <- ceiling(lambda_star_1 / (2 * k * S_a))
    MSSN_ctrl_a <- ceiling(k * MSSN_case_a)
    
    allelic_result <- list(
      test = "case-control chi-square test of independence for alleles",
      df = 1,
      lambda_star = lambda_star_1,
      S = S_a,
      MSSN_case = MSSN_case_a,
      MSSN_ctrl = MSSN_ctrl_a,
      MSSN_total = MSSN_case_a + MSSN_ctrl_a
    )
  }
  
  # ---- trend test ----
  num_t <- (sum(w * (g1_obs - g0_obs)))^2
  
  den_t <- sum(w^2 * (g1_obs + k * g0_obs)) -
    (sum(w * (g1_obs + k * g0_obs)))^2 / (1 + k)
  
  if (!is.finite(den_t) || den_t <= 0)
    stop("Trend denominator <= 0; check inputs/weights.")
  
  if (num_t < 1e-15)
    stop("Trend numerator is approximately 0; implies no weighted mean difference.")
  
  S_t <- num_t / den_t
  
  MSSN_case_t <- ceiling(lambda_star_1 / (k * S_t))
  MSSN_ctrl_t <- ceiling(k * MSSN_case_t)
  
  # ---- output ----
  out <- list(
    alpha = alpha,
    target_power = power,
    input_mode = input_mode,
    k = k,
    w = w,
    include_allelic = include_allelic,
    locus_het = locus_het_info,
    errors = list(
      genotype_misclass = misclass_info
    ),
    model_info = model_info,
    tests = list(
      genotypes = list(
        test = "case-control chi-square test of independence for genotypes",
        df = 2,
        lambda_star = lambda_star_g,
        S = S_g,
        MSSN_case = MSSN_case_g,
        MSSN_ctrl = MSSN_ctrl_g,
        MSSN_total = MSSN_case_g + MSSN_ctrl_g
      ),
      alleles = allelic_result,
      trend = list(
        test = "trend test for genotypes",
        df = 1,
        lambda_star = lambda_star_1,
        S = S_t,
        numerator = num_t,
        denominator = den_t,
        MSSN_case = MSSN_case_t,
        MSSN_ctrl = MSSN_ctrl_t,
        MSSN_total = MSSN_case_t + MSSN_ctrl_t
      )
    ),
    freqs = list(
      g_base_case = g1_base,
      g_base_ctrl = g0_base,
      g_true_case = g1_true,
      g_true_ctrl = g0_true,
      g_obs_case  = g1_obs,
      g_obs_ctrl  = g0_obs,
      p_obs_case  = p1_obs,
      p_obs_ctrl  = p0_obs
    )
  )
  
  class(out) <- "cc_mssn_conditional_full"
  
  # ---- clean printed output ----
  if (isTRUE(verbose)) {
    
    message("\n--- Case-Control: Minimum Sample Size Necessary (MSSN) ---")
    if (isTRUE(include_allelic)) {
      message("Outputs: Genotype chi-square, allelic chi-square, and genotype trend test")
    } else {
      message("Outputs: Genotype chi-square and genotype trend test")
    }
    message("--------------------------------------------------------------------------")
    
    fmt2 <- "%-32s %12s  |  %-28s %12s"
    
    message(sprintf(
      fmt2,
      "Input Mode:", input_mode,
      "Significance Level (alpha):", .fmt_e(alpha, 2)
    ))
    
    message(sprintf(
      fmt2,
      "Target Power:", .fmt_f(power, 3),
      "Case:Control Ratio (k):", .fmt_f(k, 3)
    ))
    
    message(sprintf(
      "%-32s %12s",
      "Trend Weights (w):", paste0(w, collapse = ",")
    ))
    
    if (input_mode == "model_based") {
      message(sprintf(
        fmt2,
        "Disease Prevalence (prev):", .fmt_f(prev, 4),
        "Risk Allele Freq (p_d):", .fmt_f(pd, 4)
      ))
      message(sprintf(
        fmt2,
        "MOI:", MOI,
        "R2:", .fmt_f(R2, 4)
      ))
    } else if (input_mode == "model_free") {
      message("Model-free input: user-supplied genotype frequencies g1 and g0")
    }
    
    if (isTRUE(locus_het)) {
      message(sprintf(
        "%-32s %12s",
        "Locus heterogeneity:", paste0("enabled, pi=", .fmt_f(pi, 3))
      ))
    } else {
      message(sprintf(
        "%-32s %12s",
        "Locus heterogeneity:", "none"
      ))
    }
    
    if (geno_misclass == "none") {
      message(sprintf("%-32s %12s", "Genotype misclassification:", "none"))
    } else if (geno_misclass == "1p") {
      message(sprintf(
        "%-32s %12s",
        "Genotype misclassification:",
        paste0("1-parameter, e=", .fmt_f(e, 4))
      ))
    } else if (geno_misclass == "2p") {
      message(sprintf(
        "%-32s %12s",
        "Genotype misclassification:",
        paste0("2-parameter, e1=", .fmt_f(e1, 4), ", e2=", .fmt_f(e2, 4))
      ))
    } else if (geno_misclass == "3p") {
      message(sprintf(
        "%-32s %12s",
        "Genotype misclassification:",
        paste0(
          "3-parameter, e01=", .fmt_f(e01, 4),
          ", e02=", .fmt_f(e02, 4),
          ", e03=", .fmt_f(e03, 4)
        )
      ))
    } else if (geno_misclass == "diff3p") {
      message("Genotype misclassification:     differential 3-parameter")
      message(sprintf("  diff_source:                  %s", misclass_info$diff_source))
      message(sprintf("  diff_multiplier:              %s", .fmt_f(misclass_info$diff_multiplier, 4)))
      message(sprintf(
        "  Case parameters:              e01=%s, e02=%s, e03=%s",
        .fmt_f(misclass_info$case_params["e01"], 4),
        .fmt_f(misclass_info$case_params["e02"], 4),
        .fmt_f(misclass_info$case_params["e03"], 4)
      ))
      message(sprintf(
        "  Control parameters:           e01=%s, e02=%s, e03=%s",
        .fmt_f(misclass_info$ctrl_params["e01"], 4),
        .fmt_f(misclass_info$ctrl_params["e02"], 4),
        .fmt_f(misclass_info$ctrl_params["e03"], 4)
      ))
    }
    
    message("--------------------------------------------------------------------------")
    message("Minimum Sample Size Necessary")
    
    message(sprintf(
      "  %-16s MSSN_case=%8d  |  MSSN_ctrl=%8d  |  MSSN_total=%8d",
      "Genotypes:", MSSN_case_g, MSSN_ctrl_g, MSSN_case_g + MSSN_ctrl_g
    ))
    
    if (isTRUE(include_allelic)) {
      message(sprintf(
        "  %-16s MSSN_case=%8d  |  MSSN_ctrl=%8d  |  MSSN_total=%8d",
        "Alleles:", MSSN_case_a, MSSN_ctrl_a, MSSN_case_a + MSSN_ctrl_a
      ))
    }
    
    message(sprintf(
      "  %-16s MSSN_case=%8d  |  MSSN_ctrl=%8d  |  MSSN_total=%8d",
      "Trend:", MSSN_case_t, MSSN_ctrl_t, MSSN_case_t + MSSN_ctrl_t
    ))
    
    message("--------------------------------------------------------------------------")
    message("Observed genotype frequencies: cases vs controls")
    message(sprintf("  g0: %6.3f vs %6.3f", g1_obs[1], g0_obs[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g1_obs[2], g0_obs[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g1_obs[3], g0_obs[3]))
    
    if (isTRUE(include_allelic)) {
      message("Observed risk-allele frequencies")
      message(sprintf("  p_case: %6.3f", unname(p1_obs["p"])))
      message(sprintf("  p_ctrl: %6.3f", unname(p0_obs["p"])))
    }
    
    message("--------------------------------------------------------------------------")
  }
  
  invisible(out)
}





#' Case-Control Power for Conditional Genotype Frequencies
#'
#' Computes power for fixed case-control sample sizes using conditional genotype
#' frequencies. The function supports model-based genotype frequencies,
#' model-free genotype frequencies, optional locus heterogeneity, optional
#' genotype misclassification, genotype chi-square tests, optional allelic
#' chi-square tests, and genotype trend tests.
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
#'   case genotype frequencies before misclassification.
#' @param pi Numeric in \eqn{[0,1]}. Locus-homogeneity fraction used when
#'   \code{locus_het = TRUE}.
#' @param k Numeric \eqn{> 0}. Control-to-case sample size ratio
#'   \eqn{N_{ctrl} / N_{case}}.
#' @param w Numeric vector of length 3. Genotype trend-test scores. The three
#'   weights cannot all be equal.
#' @param include_allelic Logical. If \code{TRUE}, includes the allelic
#'   chi-square test in the returned object and verbose output.
#' @param geno_misclass Character. Genotype misclassification model:
#'   \code{"none"}, \code{"1p"}, \code{"2p"}, \code{"3p"}, or \code{"diff3p"}.
#' @param e Numeric. Symmetric one-parameter misclassification rate for
#'   \code{geno_misclass = "1p"}.
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
#' \item Optionally apply genotype misclassification matrices to the true case
#' and control genotype frequencies.
#' \item Compute genotype chi-square, optional allelic chi-square, and genotype
#' trend-test non-centrality parameters and powers from the observed genotype
#' frequencies.
#' }
#'
#' With \code{input_mode = "model_based"}, conditional case and control genotype
#' frequencies are derived from \code{prev}, \code{pd}, \code{R2}, and
#' \code{MOI}. With \code{input_mode = "model_free"}, the user supplies
#' \code{g1} and \code{g0} directly.
#'
#' The genotype misclassification models are:
#' \code{"none"} for identity matrices, \code{"1p"} for one symmetric error
#' rate, \code{"2p"} for adjacent homozygote/heterozygote and heterozygote
#' error rates, \code{"3p"} for a non-differential three-parameter matrix, and
#' \code{"diff3p"} for separate case and control three-parameter matrices.
#'
#' For \code{geno_misclass = "diff3p"}, \code{diff_source = "explicit"} uses
#' the case and control error parameters exactly as supplied. With
#' \code{diff_source = "case"}, control parameters are computed by multiplying
#' the case parameters by \code{diff_multiplier}. With
#' \code{diff_source = "ctrl"}, case parameters are computed by multiplying the
#' control parameters by \code{diff_multiplier}.
#'
#' If \code{include_allelic = FALSE}, the allelic test is omitted and
#' \code{tests$alleles} is \code{NULL}. Internal effect-size components
#' \code{S} are retained in the returned object for validation and debugging,
#' but are not printed in the clean verbose output.
#'
#' @return An object of class \code{"cc_power_conditional_full"}: a nested list
#' with components \code{alpha}, \code{N_case}, \code{N_ctrl}, \code{N_total},
#' \code{input_mode}, \code{k}, \code{w}, \code{include_allelic},
#' \code{locus_het}, \code{errors}, \code{model_info}, \code{tests}, and
#' \code{freqs}. \code{tests$genotypes} and \code{tests$trend} contain test
#' labels, degrees of freedom, non-centrality parameters, internal \code{S}
#' components, and power. \code{tests$alleles} has the same power fields when
#' requested, otherwise it is \code{NULL}. \code{freqs} stores baseline, true,
#' and observed case/control genotype frequencies and observed allele
#' frequencies.
#'
#' @examples
#' cc_power_conditional_full(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   geno_misclass = "none",
#'   verbose = FALSE
#' )
#'
#' cc_power_conditional_full(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   verbose = FALSE
#' )
#'
#' cc_power_conditional_full(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   locus_het = TRUE, pi = 0.8,
#'   geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' cc_power_conditional_full(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_based",
#'   prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M",
#'   locus_het = TRUE, pi = 0.8,
#'   geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' cc_power_conditional_full(
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
#' cc_power_conditional_full(
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
#' cc_power_conditional_full(
#'   N_case = 500, alpha = 0.05,
#'   input_mode = "model_free",
#'   g1 = c(0.25, 0.50, 0.25),
#'   g0 = c(0.36, 0.48, 0.16),
#'   include_allelic = FALSE,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' TODO: Provide exact textbook equation numbers/pages for the case-control
#' genotype, allelic, trend, locus-heterogeneity, and genotype-misclassification
#' calculations.
#'
#' @importFrom stats pchisq qchisq
#' @export
cc_power_conditional_full <- function(
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
    
    k = 1,
    w = c(0, 1, 2),
    include_allelic = TRUE,
    
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
  
  cc_geno_to_allele_freqs <- function(gj) {
    check_genotype_freqs(gj, "gj")
    p <- gj[3] + 0.5 * gj[2]
    c(q = 1 - p, p = p)
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
  if (!is.logical(include_allelic) || length(include_allelic) != 1)
    stop("include_allelic must be TRUE or FALSE.")
  if (!is.logical(locus_het) || length(locus_het) != 1)
    stop("locus_het must be TRUE or FALSE.")
  if (!is.numeric(pi) || length(pi) != 1 || pi < 0 || pi > 1)
    stop("pi must be a single number in [0,1].")
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
  
  p1_obs <- cc_geno_to_allele_freqs(g1_obs)
  p0_obs <- cc_geno_to_allele_freqs(g0_obs)
  
  # ---- sample sizes ----
  N_ctrl <- k * N_case
  
  # ---- genotype chi-square ----
  S_g <- sum((g1_obs - g0_obs)^2 / (g1_obs + k * g0_obs))
  
  if (!is.finite(S_g) || S_g <= 0)
    stop("Genotype S <= 0; check inputs.")
  
  lambda_g <- k * N_case * S_g
  crit_g <- qchisq(1 - alpha, df = 2)
  power_g <- pchisq(crit_g, df = 2, ncp = lambda_g, lower.tail = FALSE)
  
  # ---- allelic chi-square ----
  allelic_result <- NULL
  
  if (isTRUE(include_allelic)) {
    
    S_a <- sum((p1_obs - p0_obs)^2 / (p1_obs + k * p0_obs))
    
    if (!is.finite(S_a) || S_a <= 0)
      stop("Allele S <= 0; check inputs.")
    
    lambda_a <- 2 * k * N_case * S_a
    crit_1 <- qchisq(1 - alpha, df = 1)
    power_a <- pchisq(crit_1, df = 1, ncp = lambda_a, lower.tail = FALSE)
    
    allelic_result <- list(
      test = "case-control chi-square test of independence for alleles",
      df = 1,
      lambda = lambda_a,
      S = S_a,
      power = power_a
    )
  }
  
  # ---- trend test ----
  num_t <- (sum(w * (g1_obs - g0_obs)))^2
  
  den_t <- sum(w^2 * (g1_obs + k * g0_obs)) -
    (sum(w * (g1_obs + k * g0_obs)))^2 / (1 + k)
  
  if (!is.finite(den_t) || den_t <= 0)
    stop("Trend denominator <= 0; check inputs/weights.")
  
  if (num_t < 1e-15)
    stop("Trend numerator is approximately 0; implies no weighted mean difference.")
  
  S_t <- num_t / den_t
  lambda_t <- k * N_case * S_t
  
  crit_1 <- qchisq(1 - alpha, df = 1)
  power_t <- pchisq(crit_1, df = 1, ncp = lambda_t, lower.tail = FALSE)
  
  # ---- output ----
  out <- list(
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    input_mode = input_mode,
    k = k,
    w = w,
    include_allelic = include_allelic,
    locus_het = locus_het_info,
    errors = list(
      genotype_misclass = misclass_info
    ),
    model_info = model_info,
    tests = list(
      genotypes = list(
        test = "case-control chi-square test of independence for genotypes",
        df = 2,
        lambda = lambda_g,
        S = S_g,
        power = power_g
      ),
      alleles = allelic_result,
      trend = list(
        test = "trend test for genotypes",
        df = 1,
        lambda = lambda_t,
        S = S_t,
        numerator = num_t,
        denominator = den_t,
        power = power_t
      )
    ),
    freqs = list(
      g_base_case = g1_base,
      g_base_ctrl = g0_base,
      g_true_case = g1_true,
      g_true_ctrl = g0_true,
      g_obs_case  = g1_obs,
      g_obs_ctrl  = g0_obs,
      p_obs_case  = p1_obs,
      p_obs_ctrl  = p0_obs
    )
  )
  
  class(out) <- "cc_power_conditional_full"
  
  # ---- clean printed output ----
  if (isTRUE(verbose)) {
    
    message("\n--- Case-Control: Power for Fixed Sample Size ---")
    if (isTRUE(include_allelic)) {
      message("Outputs: Genotype chi-square, allelic chi-square, and genotype trend test")
    } else {
      message("Outputs: Genotype chi-square and genotype trend test")
    }
    message("--------------------------------------------------------------------------")
    
    fmt2 <- "%-32s %12s  |  %-28s %12s"
    
    message(sprintf(
      fmt2,
      "Input Mode:", input_mode,
      "Significance Level (alpha):", .fmt_e(alpha, 2)
    ))
    
    message(sprintf(
      fmt2,
      "N_case:", .fmt_f(N_case, 0),
      "N_ctrl:", .fmt_f(N_ctrl, 0)
    ))
    
    message(sprintf(
      fmt2,
      "Case:Control Ratio (k):", .fmt_f(k, 3),
      "Trend Weights (w):", paste0(w, collapse = ",")
    ))
    
    if (input_mode == "model_based") {
      message(sprintf(
        fmt2,
        "Disease Prevalence (prev):", .fmt_f(prev, 4),
        "Risk Allele Freq (p_d):", .fmt_f(pd, 4)
      ))
      message(sprintf(
        fmt2,
        "MOI:", MOI,
        "R2:", .fmt_f(R2, 4)
      ))
    } else if (input_mode == "model_free") {
      message("Model-free input: user-supplied genotype frequencies g1 and g0")
    }
    
    if (isTRUE(locus_het)) {
      message(sprintf(
        "%-32s %12s",
        "Locus heterogeneity:", paste0("enabled, pi=", .fmt_f(pi, 3))
      ))
    } else {
      message(sprintf(
        "%-32s %12s",
        "Locus heterogeneity:", "none"
      ))
    }
    
    if (geno_misclass == "none") {
      message(sprintf("%-32s %12s", "Genotype misclassification:", "none"))
    } else if (geno_misclass == "1p") {
      message(sprintf(
        "%-32s %12s",
        "Genotype misclassification:",
        paste0("1-parameter, e=", .fmt_f(e, 4))
      ))
    } else if (geno_misclass == "2p") {
      message(sprintf(
        "%-32s %12s",
        "Genotype misclassification:",
        paste0("2-parameter, e1=", .fmt_f(e1, 4), ", e2=", .fmt_f(e2, 4))
      ))
    } else if (geno_misclass == "3p") {
      message(sprintf(
        "%-32s %12s",
        "Genotype misclassification:",
        paste0(
          "3-parameter, e01=", .fmt_f(e01, 4),
          ", e02=", .fmt_f(e02, 4),
          ", e03=", .fmt_f(e03, 4)
        )
      ))
    } else if (geno_misclass == "diff3p") {
      message("Genotype misclassification:     differential 3-parameter")
      message(sprintf("  diff_source:                  %s", misclass_info$diff_source))
      message(sprintf("  diff_multiplier:              %s", .fmt_f(misclass_info$diff_multiplier, 4)))
      message(sprintf(
        "  Case parameters:              e01=%s, e02=%s, e03=%s",
        .fmt_f(misclass_info$case_params["e01"], 4),
        .fmt_f(misclass_info$case_params["e02"], 4),
        .fmt_f(misclass_info$case_params["e03"], 4)
      ))
      message(sprintf(
        "  Control parameters:           e01=%s, e02=%s, e03=%s",
        .fmt_f(misclass_info$ctrl_params["e01"], 4),
        .fmt_f(misclass_info$ctrl_params["e02"], 4),
        .fmt_f(misclass_info$ctrl_params["e03"], 4)
      ))
    }
    
    message("--------------------------------------------------------------------------")
    message("Power")
    
    message(sprintf(
      "  %-16s %12.6f",
      "Genotypes:", power_g
    ))
    
    if (isTRUE(include_allelic)) {
      message(sprintf(
        "  %-16s %12.6f",
        "Alleles:", power_a
      ))
    }
    
    message(sprintf(
      "  %-16s %12.6f",
      "Trend:", power_t
    ))
    
    message("--------------------------------------------------------------------------")
    message("Non-Centrality Parameters")
    
    message(sprintf(
      "  %-16s %12.5f  | df=%d",
      "Genotypes:", lambda_g, 2
    ))
    
    if (isTRUE(include_allelic)) {
      message(sprintf(
        "  %-16s %12.5f  | df=%d",
        "Alleles:", lambda_a, 1
      ))
    }
    
    message(sprintf(
      "  %-16s %12.5f  | df=%d",
      "Trend:", lambda_t, 1
    ))
    
    message("--------------------------------------------------------------------------")
    message("Observed genotype frequencies: cases vs controls")
    message(sprintf("  g0: %6.3f vs %6.3f", g1_obs[1], g0_obs[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g1_obs[2], g0_obs[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g1_obs[3], g0_obs[3]))
    
    if (isTRUE(include_allelic)) {
      message("Observed risk-allele frequencies")
      message(sprintf("  p_case: %6.3f", unname(p1_obs["p"])))
      message(sprintf("  p_ctrl: %6.3f", unname(p0_obs["p"])))
    }
    
    message("--------------------------------------------------------------------------")
  }
  
  invisible(out)
}
