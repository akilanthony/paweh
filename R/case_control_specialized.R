#' Specialized Case-Control Genotype Misclassification Functions
#'
#' Convenience functions for genotype-only case-control chi-square calculations
#' under conditional model-based genotype frequencies and genotype
#' misclassification. These functions are narrower than
#' \code{\link{cc_mssn}} and
#' \code{\link{cc_power}}: they compute only the genotype
#' chi-square test for one specific misclassification model.
#'
#' @param power Numeric in \eqn{(0,1)}. Desired target power for MSSN functions.
#' @param N_case Numeric \eqn{> 0}. Number of cases for power functions.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence.
#' @param pd Numeric in \eqn{(0,1)}. Disease-allele frequency.
#' @param R2 Numeric \eqn{> 0}. Homozygote relative risk.
#' @param MOI Character. Mode of inheritance: \code{"M"}, \code{"D"}, or
#'   \code{"Rec"}.
#' @param k Numeric \eqn{> 0}. Control-to-case ratio \eqn{N_{ctrl}/N_{case}}.
#' @param e Numeric. One-parameter symmetric genotype error rate.
#' @param e1,e2 Numeric. Two-parameter genotype misclassification rates.
#' @param e01,e02,e03 Numeric. Three-parameter non-differential genotype
#'   misclassification rates.
#' @param case_e01,case_e02,case_e03 Numeric. Case-specific three-parameter
#'   genotype misclassification rates for differential misclassification.
#' @param ctrl_e01,ctrl_e02,ctrl_e03 Numeric. Control-specific three-parameter
#'   genotype misclassification rates for differential misclassification.
#' @param verbose Logical. If \code{TRUE}, prints a formatted summary.
#'
#' @details
#' These are specialized convenience functions retained for users who want a
#' focused genotype-only calculation. The all-in-one case-control functions are
#' recommended when combining model-free inputs, locus heterogeneity, trend
#' tests, or multiple modifier choices in one call.
#'
#' The returned objects keep internal effect-size components such as \code{S}
#' and, for differential misclassification, per-genotype components, for
#' validation and debugging.
#'
#' @return A list with class matching the function name. MSSN functions include
#' target non-centrality parameter, case/control sample sizes, model parameters,
#' penetrances, misclassification matrices, internal \code{S} values, and true
#' and observed genotype frequencies. Power functions include case/control
#' sample sizes, non-centrality parameter, power, model parameters,
#' penetrances, misclassification matrices, internal \code{S} values, and true
#' and observed genotype frequencies.
#'
#' @examples
#' cc_chisq_mssn_genotype_misclassification_1p(
#'   power = 0.8, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
#'   MOI = "D", e = 0.02, verbose = FALSE
#' )
#'
#' cc_chisq_power_genotype_misclassification_3p(
#'   N_case = 500, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
#'   MOI = "D", e01 = 0.02, e02 = 0.01, e03 = 0.005,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' TODO: Provide exact textbook equation numbers/pages for the specialized
#' case-control genotype misclassification calculations.
#'
#' @name case_control_genotype_misclassification
#' @importFrom stats pchisq qchisq uniroot
#' @export
cc_chisq_mssn_genotype_misclassification_1p <- function(
    power, alpha, prev, pd, R2,
    MOI = c("M","D","Rec"),
    k = 1,
    e = 0,
    verbose = TRUE
) {
  MOI <- match.arg(MOI)


  chisq_ncp_target <- function(power, alpha, df) {
    crit <- qchisq(1 - alpha, df = df)
    f <- function(lambda) pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - power
    uniroot(f, lower = 0, upper = 1e6)$root
  }

  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M","D","Rec")) {
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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  cc_misclass_matrix_symmetric <- function(e) {
    if (!is.numeric(e) || length(e) != 1 || e < 0 || e > 0.5)
      stop("e must be a single number in [0, 0.5].")
    matrix(c(
      1 - 2*e, e,       e,
      e,       1 - 2*e, e,
      e,       e,       1 - 2*e
    ), nrow = 3, byrow = TRUE)
  }

  cc_apply_genotype_misclass <- function(g, M) {
    as.numeric(M %*% g)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)


  # checks
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("power must be a single number in (0,1).")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")


  # target NCP for genotype chi-square (df=2)
  lambda_star <- chisq_ncp_target(power = power, alpha = alpha, df = 2)

  # true conditional genotype freqs
  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0


  M <- cc_misclass_matrix_symmetric(e)
  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M)


  g_obs_case <- g_obs_case / sum(g_obs_case)
  g_obs_ctrl <- g_obs_ctrl / sum(g_obs_ctrl)

  # Eq 1.22 (k-form) effect-size component
  S <- sum((g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl))
  if (!is.finite(S) || S <= 0) stop("Genotype S <= 0; check inputs (or error model).")

  # MSSN
  N_case <- ceiling(lambda_star / (k * S))
  N_ctrl <- ceiling(k * N_case)

  out <- list(
    alpha = alpha,
    target_power = power,
    k = k,
    e = e,
    lambda_star = lambda_star,
    N_case = N_case,
    N_ctrl = N_ctrl,
    S = S,
    model_parameters = list(prev = prev, pd = pd, qd = 1 - pd, R1 = freqs$R1, R2 = R2, MOI = MOI),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    misclassification = list(model = "symmetric", M = M),
    freqs = list(
      g_true_case = g_true_case,
      g_true_ctrl = g_true_ctrl,
      g_obs_case  = g_obs_case,
      g_obs_ctrl  = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_mssn_genotype_misclassification_1p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): MSSN with Genotype Misclassification ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df=2")
    message("Misclassification: symmetric matrix (diag=1-2e, offdiag=e)")
    message("--------------------------------------------------------------------------")

    fmt2 <- "%-32s %12s  |  %-28s %12s"
    message(sprintf(fmt2, "Target Power:", .fmt_f(power, 3),
                    "Significance Level (alpha):", .fmt_e(alpha, 2)))

    message(sprintf(fmt2, "Disease Prevalence (prev):", .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):", .fmt_f(pd, 4)))

    message(sprintf(fmt2, "Wild-type Freq (p_plus):", .fmt_f(1 - pd, 4),
                    "MOI:", MOI))

    message(sprintf(fmt2, "Case:Control Ratio (k):", .fmt_f(k, 3),
                    "Misclass error (e):", .fmt_f(e, 4)))

    message(sprintf(fmt2, "Relative Risks (R1,R2):",
                    paste0(.fmt_f(out$model_parameters$R1, 3), ", ", .fmt_f(R2, 3)),
                    "Target NCP (lambda*):",
                    formatC(lambda_star, format = "f", digits = 5)))

    message("--------------------------------------------------------------------------")
    message("True conditional genotype frequencies (cases vs controls)")
    message(sprintf("  g0: %6.3f vs %6.3f", g_true_case[1], g_true_ctrl[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g_true_case[2], g_true_ctrl[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g_true_case[3], g_true_ctrl[3]))

    message("Observed conditional genotype frequencies (after misclassification)")
    message(sprintf("  g0~: %6.6f vs %6.6f", g_obs_case[1], g_obs_ctrl[1]))
    message(sprintf("  g1~: %6.6f vs %6.6f", g_obs_case[2], g_obs_ctrl[2]))
    message(sprintf("  g2~: %6.6f vs %6.6f", g_obs_case[3], g_obs_ctrl[3]))

    message("--------------------------------------------------------------------------")
    message("Effect-size component S and required sample sizes")
    message(sprintf("  %-16s %12.6g", "S:", S))
    message(sprintf("  %-16s N_case=%8d  |  N_ctrl=%8d", "MSSN:", N_case, N_ctrl))
    message("--------------------------------------------------------------------------")
  }

  invisible(out)
}
#' @rdname case_control_genotype_misclassification
#' @export
cc_chisq_power_genotype_misclassification_1p <- function(
    N_case, alpha, prev, pd, R2,
    MOI = c("M","D","Rec"),
    k = 1,
    e = 0,
    verbose = TRUE
) {
  MOI <- match.arg(MOI)


  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M","D","Rec")) {
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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  cc_misclass_matrix_symmetric <- function(e) {
    if (!is.numeric(e) || length(e) != 1 || e < 0 || e > 0.5)
      stop("e must be a single number in [0, 0.5].")
    matrix(c(
      1 - 2*e, e,       e,
      e,       1 - 2*e, e,
      e,       e,       1 - 2*e
    ), nrow = 3, byrow = TRUE)
  }

  cc_apply_genotype_misclass <- function(g, M) {
    as.numeric(M %*% g)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)


  if (!is.numeric(N_case) || length(N_case) != 1 || N_case <= 0)
    stop("N_case must be a single positive number.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")

  N_ctrl <- k * N_case


  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0


  M <- cc_misclass_matrix_symmetric(e)
  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M)
  g_obs_case <- g_obs_case / sum(g_obs_case)
  g_obs_ctrl <- g_obs_ctrl / sum(g_obs_ctrl)

  # Eq 1.22 (k-form) NCP
  S <- sum((g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl))
  if (!is.finite(S) || S <= 0) stop("Genotype S <= 0; check inputs (or error model).")

  lambda <- k * N_case * S

  # power (df=2)
  crit <- qchisq(1 - alpha, df = 2)
  power <- pchisq(crit, df = 2, ncp = lambda, lower.tail = FALSE)

  out <- list(
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    k = k,
    e = e,
    lambda = lambda,
    power = power,
    S = S,
    model_parameters = list(prev = prev, pd = pd, qd = 1 - pd, R1 = freqs$R1, R2 = R2, MOI = MOI),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    misclassification = list(model = "symmetric", M = M),
    freqs = list(
      g_true_case = g_true_case,
      g_true_ctrl = g_true_ctrl,
      g_obs_case  = g_obs_case,
      g_obs_ctrl  = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_power_genotype_misclassification_1p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): Power with Genotype Misclassification ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df=2")
    message("Misclassification: symmetric matrix (diag=1-2e, offdiag=e)")
    message("-------------------------------------------------------------------")

    fmt2 <- "%-32s %12s  |  %-28s %12s"
    message(sprintf(fmt2, "N_case:", formatC(N_case, format = "f", digits = 0),
                    "N_ctrl:", formatC(N_ctrl, format = "f", digits = 0)))

    message(sprintf(fmt2, "Significance Level (alpha):", .fmt_e(alpha, 2),
                    "Case:Control Ratio (k):", .fmt_f(k, 3)))

    message(sprintf(fmt2, "Disease Prevalence (prev):", .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):", .fmt_f(pd, 4)))

    message(sprintf(fmt2, "Wild-type Freq (p_plus):", .fmt_f(1 - pd, 4),
                    "MOI:", MOI))

    message(sprintf(fmt2, "Relative Risks (R1,R2):",
                    paste0(.fmt_f(out$model_parameters$R1, 3), ", ", .fmt_f(R2, 3)),
                    "Misclass error (e):", .fmt_f(e, 4)))

    message("-------------------------------------------------------------------")
    message("Observed conditional genotype frequencies (after misclassification)")
    message(sprintf("  g0~: %6.6f vs %6.6f", g_obs_case[1], g_obs_ctrl[1]))
    message(sprintf("  g1~: %6.6f vs %6.6f", g_obs_case[2], g_obs_ctrl[2]))
    message(sprintf("  g2~: %6.6f vs %6.6f", g_obs_case[3], g_obs_ctrl[3]))

    message("-------------------------------------------------------------------")
    message("Non-Centrality Parameter (lambda) and Power")
    message(sprintf("  %-16s %12.5f  | df=%d", "lambda:", lambda, 2))
    message(sprintf("  %-16s %12.6f", "power:", power))
    message("-------------------------------------------------------------------")
  }

  invisible(out)
}










# -------------------------------------------------------------------------
### two parameter genotype error model
# -------------------------------------------------------------------
#' @rdname case_control_genotype_misclassification
#' @export
cc_chisq_mssn_genotype_misclassification_2p <- function(
    power, alpha, prev, pd, R2,
    MOI = c("M","D","Rec"),
    k = 1,
    e1 = 0,
    e2 = 0,
    verbose = TRUE
) {
  MOI <- match.arg(MOI)

  # ---- local helpers ----
  chisq_ncp_target <- function(power, alpha, df) {
    crit <- qchisq(1 - alpha, df = df)
    f <- function(lambda) pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - power
    uniroot(f, lower = 0, upper = 1e6)$root
  }

  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M","D","Rec")) {
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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  # Excel orientation: rows=true, cols=observed; rows sum to 1
  cc_misclass_matrix_two_param <- function(e1, e2) {
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
      stop("Misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("Misclassification matrix has negative entries; check e1/e2.")

    M
  }

  cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
    if (length(g_true) != 3) stop("g_true must be length 3 (g0,g1,g2).")
    as.numeric(t(M_true_to_obs) %*% g_true)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)

  # ---- checks ----
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("power must be a single number in (0,1).")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")

  # ---- target NCP (df=2) ----
  lambda_star <- chisq_ncp_target(power = power, alpha = alpha, df = 2)

  # ---- true conditional genotype freqs ----
  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0

  # ---- apply 2p misclassification ----
  M <- cc_misclass_matrix_two_param(e1, e2)
  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M)

  # defensive renorm
  g_obs_case <- g_obs_case / sum(g_obs_case)
  g_obs_ctrl <- g_obs_ctrl / sum(g_obs_ctrl)

  # ---- Eq 1.22 (k-form): S ----
  S <- sum((g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl))
  if (!is.finite(S) || S <= 0) stop("Genotype S <= 0; check inputs (or error model).")

  # ---- MSSN ----
  N_case <- ceiling(lambda_star / (k * S))
  N_ctrl <- ceiling(k * N_case)

  out <- list(
    alpha = alpha,
    target_power = power,
    k = k,
    e1 = e1,
    e2 = e2,
    lambda_star = lambda_star,
    N_case = N_case,
    N_ctrl = N_ctrl,
    S = S,
    model_parameters = list(prev = prev, pd = pd, qd = 1 - pd, R1 = freqs$R1, R2 = R2, MOI = MOI),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    misclassification = list(model = "two_param_hom_het", M = M),
    freqs = list(
      g_true_case = g_true_case,
      g_true_ctrl = g_true_ctrl,
      g_obs_case  = g_obs_case,
      g_obs_ctrl  = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_mssn_genotype_misclassification_2p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): MSSN with Genotype Misclassification (2-parameter) ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df=2")
    message("Misclassification: hom<->het two-parameter (e1=hom->het, e2=het->hom)")
    message("--------------------------------------------------------------------------")

    fmt2 <- "%-32s %12s  |  %-28s %12s"
    message(sprintf(fmt2, "Target Power:", .fmt_f(power, 3),
                    "Significance Level (alpha):", .fmt_e(alpha, 2)))

    message(sprintf(fmt2, "Disease Prevalence (prev):", .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):", .fmt_f(pd, 4)))

    message(sprintf(fmt2, "Wild-type Freq (p_plus):", .fmt_f(1 - pd, 4),
                    "MOI:", MOI))

    message(sprintf(fmt2, "Case:Control Ratio (k):", .fmt_f(k, 3),
                    "Misclass (e1,e2):", paste0(.fmt_f(e1, 4), ", ", .fmt_f(e2, 4))))

    message(sprintf(fmt2, "Relative Risks (R1,R2):",
                    paste0(.fmt_f(out$model_parameters$R1, 3), ", ", .fmt_f(R2, 3)),
                    "Target NCP (lambda*):",
                    formatC(lambda_star, format = "f", digits = 5)))

    message("--------------------------------------------------------------------------")
    message("True conditional genotype frequencies (cases vs controls)")
    message(sprintf("  g0: %6.3f vs %6.3f", g_true_case[1], g_true_ctrl[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g_true_case[2], g_true_ctrl[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g_true_case[3], g_true_ctrl[3]))

    message("Observed conditional genotype frequencies (after misclassification)")
    message(sprintf("  g0~: %6.6f vs %6.6f", g_obs_case[1], g_obs_ctrl[1]))
    message(sprintf("  g1~: %6.6f vs %6.6f", g_obs_case[2], g_obs_ctrl[2]))
    message(sprintf("  g2~: %6.6f vs %6.6f", g_obs_case[3], g_obs_ctrl[3]))

    message("--------------------------------------------------------------------------")
    message("Effect-size component S and required sample sizes")
    message(sprintf("  %-16s %12.6g", "S:", S))
    message(sprintf("  %-16s N_case=%8d  |  N_ctrl=%8d", "MSSN:", N_case, N_ctrl))
    message("--------------------------------------------------------------------------")
  }

  invisible(out)
}
#' @rdname case_control_genotype_misclassification
#' @export
cc_chisq_power_genotype_misclassification_2p <- function(
    N_case, alpha, prev, pd, R2,
    MOI = c("M","D","Rec"),
    k = 1,
    e1 = 0,
    e2 = 0,
    verbose = TRUE
) {
  MOI <- match.arg(MOI)

  # ---- local helpers ----
  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M","D","Rec")) {
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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  cc_misclass_matrix_two_param <- function(e1, e2) {
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
      stop("Misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12))
      stop("Misclassification matrix has negative entries; check e1/e2.")
    M
  }

  cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
    if (length(g_true) != 3) stop("g_true must be length 3 (g0,g1,g2).")
    as.numeric(t(M_true_to_obs) %*% g_true)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)

  # ---- checks ----
  if (!is.numeric(N_case) || length(N_case) != 1 || N_case <= 0)
    stop("N_case must be a single positive number.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")

  N_ctrl <- k * N_case

  # ---- true freqs ----
  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0

  # ---- observed freqs ----
  M <- cc_misclass_matrix_two_param(e1, e2)
  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M)

  g_obs_case <- g_obs_case / sum(g_obs_case)
  g_obs_ctrl <- g_obs_ctrl / sum(g_obs_ctrl)

  # ---- Eq 1.22 (k-form): lambda = k*N_case*S ----
  S <- sum((g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl))
  if (!is.finite(S) || S <= 0) stop("Genotype S <= 0; check inputs (or error model).")

  lambda <- k * N_case * S

  crit <- qchisq(1 - alpha, df = 2)
  power <- pchisq(crit, df = 2, ncp = lambda, lower.tail = FALSE)

  out <- list(
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    k = k,
    e1 = e1,
    e2 = e2,
    lambda = lambda,
    power = power,
    S = S,
    model_parameters = list(prev = prev, pd = pd, qd = 1 - pd, R1 = freqs$R1, R2 = R2, MOI = MOI),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    misclassification = list(model = "two_param_hom_het", M = M),
    freqs = list(
      g_true_case = g_true_case,
      g_true_ctrl = g_true_ctrl,
      g_obs_case  = g_obs_case,
      g_obs_ctrl  = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_power_genotype_misclassification_2p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): Power with Genotype Misclassification (2-parameter) ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df=2")
    message("Misclassification: hom<->het two-parameter (e1=hom->het, e2=het->hom)")
    message("-------------------------------------------------------------------")

    fmt2 <- "%-32s %12s  |  %-28s %12s"
    message(sprintf(fmt2, "N_case:", formatC(N_case, format = "f", digits = 0),
                    "N_ctrl:", formatC(N_ctrl, format = "f", digits = 0)))

    message(sprintf(fmt2, "Significance Level (alpha):", .fmt_e(alpha, 2),
                    "Case:Control Ratio (k):", .fmt_f(k, 3)))

    message(sprintf(fmt2, "Disease Prevalence (prev):", .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):", .fmt_f(pd, 4)))

    message(sprintf(fmt2, "Wild-type Freq (p_plus):", .fmt_f(1 - pd, 4),
                    "MOI:", MOI))

    message(sprintf(fmt2, "Relative Risks (R1,R2):",
                    paste0(.fmt_f(out$model_parameters$R1, 3), ", ", .fmt_f(R2, 3)),
                    "Misclass (e1,e2):", paste0(.fmt_f(e1, 4), ", ", .fmt_f(e2, 4))))

    message("-------------------------------------------------------------------")
    message("Observed conditional genotype frequencies (after misclassification)")
    message(sprintf("  g0~: %6.6f vs %6.6f", g_obs_case[1], g_obs_ctrl[1]))
    message(sprintf("  g1~: %6.6f vs %6.6f", g_obs_case[2], g_obs_ctrl[2]))
    message(sprintf("  g2~: %6.6f vs %6.6f", g_obs_case[3], g_obs_ctrl[3]))

    message("-------------------------------------------------------------------")
    message("Non-Centrality Parameter (lambda) and Power")
    message(sprintf("  %-16s %12.5f  | df=%d", "lambda:", lambda, 2))
    message(sprintf("  %-16s %12.6f", "power:", power))
    message("-------------------------------------------------------------------")
  }

  invisible(out)
}



#
#
# mssn_2p <- cc_chisq_mssn_genotype_misclassification_2p(
#   power = 0.8,
#   alpha = 5e-6,
#   prev  = 0.1,
#   pd    = 0.3,
#   R2    = 2.25,
#   MOI   = "D",
#   k     = 1,
#   e1    = 0.001,
#   e2    = 0.025,
#   verbose = TRUE
# )
#
#
# mssn_2p$N_case
# mssn_2p$N_ctrl
#
#
# mssn_2p$freqs$g_obs_case
# mssn_2p$freqs$g_obs_ctrl
#
#
# pow_check <- cc_chisq_power_genotype_misclassification_2p(
#   N_case = mssn_2p$N_case,
#   alpha = 5e-6,
#   prev  = 0.1,
#   pd    = 0.3,
#   R2    = 2.25,
#   MOI   = "D",
#   k     = 1,
#   e1    = 0.001,
#   e2    = 0.025,
#   verbose = TRUE
# )
#
# pow_check$power
#



# 3 parameter model

# -------------------------------------------------------------------
# Case-Control genotype only
# Genotype misclassification: 3-parameter model (e01, e02, e03)
#' @rdname case_control_genotype_misclassification
#' @export
cc_chisq_mssn_genotype_misclassification_3p <- function(
    power, alpha, prev, pd, R2,
    MOI = c("M","D","Rec"),
    k = 1,
    e01 = 0,  # hom -> het
    e02 = 0,  # het -> hom (split equally to 0 and 2)
    e03 = 0,  # hom -> other hom
    verbose = TRUE
) {
  MOI <- match.arg(MOI)

  #helpers
  chisq_ncp_target <- function(power, alpha, df) {
    crit <- qchisq(1 - alpha, df = df)
    f <- function(lambda) pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - power
    uniroot(f, lower = 0, upper = 1e6)$root
  }

  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M","D","Rec")) {
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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  cc_misclass_matrix_3p <- function(e01, e02, e03) {
    # checks to guarantee rows sum to 1 and probabilities are valid
    if (!is.numeric(e01) || length(e01) != 1 || e01 < 0 || e01 > 1)
      stop("e01 must be a single number in [0,1].")
    if (!is.numeric(e02) || length(e02) != 1 || e02 < 0 || e02 > 0.5)
      stop("e02 must be a single number in [0,0.5] (since row1 has 1-2e02).")
    if (!is.numeric(e03) || length(e03) != 1 || e03 < 0 || e03 > 1)
      stop("e03 must be a single number in [0,1].")

    if (e01 + e03 > 1) stop("Need e01 + e03 <= 1 so hom rows remain nonnegative.")

    M <- matrix(c(
      1 - (e01 + e03),  e01,            e03,
      e02,              1 - 2*e02,      e02,
      e03,              e01,            1 - (e01 + e03)
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10)) stop("Misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12)) stop("Misclassification matrix has negative entries; check e01/e02/e03.")
    M
  }

  cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
    if (length(g_true) != 3) stop("g_true must be length 3.")
    as.numeric(t(M_true_to_obs) %*% g_true)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)

  # checks
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("power must be a single number in (0,1).")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")

  # NCP
  lambda_star <- chisq_ncp_target(power = power, alpha = alpha, df = 2)

  # ---- true freqs ----
  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0

  # ---- misclassification ----
  M <- cc_misclass_matrix_3p(e01, e02, e03)
  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M)

  g_obs_case <- g_obs_case / sum(g_obs_case)
  g_obs_ctrl <- g_obs_ctrl / sum(g_obs_ctrl)

  # ---- Eq 1.22 (k-form) effect-size component ----
  S <- sum((g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl))
  if (!is.finite(S) || S <= 0) stop("Genotype S <= 0; check inputs (or error model).")

  # ---- MSSN ----
  N_case <- ceiling(lambda_star / (k * S))
  N_ctrl <- ceiling(k * N_case)

  out <- list(
    alpha = alpha,
    target_power = power,
    k = k,
    misclassification = list(model = "3p", e01 = e01, e02 = e02, e03 = e03, M = M),
    lambda_star = lambda_star,
    MSSN_case = N_case,
    MSSN_ctrl = N_ctrl,
    S = S,
    model_parameters = list(prev = prev, pd = pd, qd = 1 - pd, R1 = freqs$R1, R2 = R2, MOI = MOI),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    freqs = list(
      g_true_case = g_true_case, g_true_ctrl = g_true_ctrl,
      g_obs_case  = g_obs_case,  g_obs_ctrl  = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_mssn_genotype_misclassification_3p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): Minimum Sample Size Necessary (MSSN) ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df=2")
    message("Misclassification: 3-parameter matrix (e01,e02,e03)")
    message("--------------------------------------------------------------------------")

    fmt2 <- "%-34s %12s  |  %-30s %12s"
    message(sprintf(fmt2, "Target Power:", .fmt_f(power, 3),
                    "Significance Level (alpha):", .fmt_e(alpha, 2)))

    message(sprintf(fmt2, "Disease Prevalence (prev):", .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):", .fmt_f(pd, 4)))

    message(sprintf(fmt2, "Wild-type Freq (p_plus):", .fmt_f(1 - pd, 4),
                    "MOI:", MOI))

    message(sprintf(fmt2, "Case:Control Ratio (k):", .fmt_f(k, 3),
                    "Target NCP (lambda*):", formatC(lambda_star, format="f", digits=5)))

    message(sprintf(fmt2, "e01 (hom->het):", .fmt_f(e01, 4),
                    "e02 (het->hom):", .fmt_f(e02, 4)))
    message(sprintf("%-34s %12s", "e03 (hom->hom):", .fmt_f(e03, 4)))

    message("--------------------------------------------------------------------------")
    message("Observed conditional genotype frequencies (after misclassification)")
    message(sprintf("  g0~: %10.6f vs %10.6f", g_obs_case[1], g_obs_ctrl[1]))
    message(sprintf("  g1~: %10.6f vs %10.6f", g_obs_case[2], g_obs_ctrl[2]))
    message(sprintf("  g2~: %10.6f vs %10.6f", g_obs_case[3], g_obs_ctrl[3]))

    message("--------------------------------------------------------------------------")
    message("Effect-size component S and MSSN")
    message(sprintf("  %-16s %12.6g", "S:", S))
    message(sprintf("  %-16s MSSN_case=%8d  |  MSSN_ctrl=%8d", "MSSN:", N_case, N_ctrl))
    message("--------------------------------------------------------------------------")
  }

  invisible(out)
}
#' @rdname case_control_genotype_misclassification
#' @export
cc_chisq_power_genotype_misclassification_3p <- function(
    N_case, alpha, prev, pd, R2,
    MOI = c("M","D","Rec"),
    k = 1,
    e01 = 0,
    e02 = 0,
    e03 = 0,
    verbose = TRUE
) {
  MOI <- match.arg(MOI)

  #helpers
  cc_conditional_geno_freqs <- function(prev, pd, R2, MOI = c("M","D","Rec")) {
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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  cc_misclass_matrix_3p <- function(e01, e02, e03) {
    if (!is.numeric(e01) || length(e01) != 1 || e01 < 0 || e01 > 1)
      stop("e01 must be a single number in [0,1].")
    if (!is.numeric(e02) || length(e02) != 1 || e02 < 0 || e02 > 0.5)
      stop("e02 must be a single number in [0,0.5].")
    if (!is.numeric(e03) || length(e03) != 1 || e03 < 0 || e03 > 1)
      stop("e03 must be a single number in [0,1].")
    if (e01 + e03 > 1) stop("Need e01 + e03 <= 1 so hom rows remain nonnegative.")

    M <- matrix(c(
      1 - (e01 + e03),  e01,            e03,
      e02,              1 - 2*e02,      e02,
      e03,              e01,            1 - (e01 + e03)
    ), nrow = 3, byrow = TRUE)

    if (any(abs(rowSums(M) - 1) > 1e-10)) stop("Misclassification matrix rows do not sum to 1.")
    if (any(M < -1e-12)) stop("Misclassification matrix has negative entries; check e01/e02/e03.")
    M
  }

  cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
    as.numeric(t(M_true_to_obs) %*% g_true)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)

  #checks
  if (!is.numeric(N_case) || length(N_case) != 1 || N_case <= 0)
    stop("N_case must be a single positive number.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")

  N_ctrl <- k * N_case

  # ---- true freqs ----
  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0

  # ---- misclassification ----
  M <- cc_misclass_matrix_3p(e01, e02, e03)
  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M); g_obs_case <- g_obs_case / sum(g_obs_case)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M); g_obs_ctrl <- g_obs_ctrl / sum(g_obs_ctrl)

  # ---- Eq 1.22 (k-form) NCP ----
  S <- sum((g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl))
  if (!is.finite(S) || S <= 0) stop("Genotype S <= 0; check inputs (or error model).")

  lambda <- k * N_case * S

  # ---- power (df=2) ----
  crit <- qchisq(1 - alpha, df = 2)
  power_out <- pchisq(crit, df = 2, ncp = lambda, lower.tail = FALSE)

  out <- list(
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    k = k,
    misclassification = list(model = "3p", e01 = e01, e02 = e02, e03 = e03, M = M),
    lambda = lambda,
    power = power_out,
    S = S,
    model_parameters = list(prev = prev, pd = pd, qd = 1 - pd, R1 = freqs$R1, R2 = R2, MOI = MOI),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    freqs = list(
      g_true_case = g_true_case, g_true_ctrl = g_true_ctrl,
      g_obs_case  = g_obs_case,  g_obs_ctrl  = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_power_genotype_misclassification_3p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): Power for Fixed Sample Size ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df=2")
    message("Misclassification: 3-parameter matrix (e01,e02,e03)")
    message("-------------------------------------------------------------------")

    fmt2 <- "%-34s %12s  |  %-30s %12s"
    message(sprintf(fmt2, "N_case:", formatC(N_case, format="f", digits=0),
                    "N_ctrl:", formatC(N_ctrl, format="f", digits=0)))

    message(sprintf(fmt2, "Significance Level (alpha):", .fmt_e(alpha, 2),
                    "Case:Control Ratio (k):", .fmt_f(k, 3)))

    message(sprintf(fmt2, "Disease Prevalence (prev):", .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):", .fmt_f(pd, 4)))

    message(sprintf(fmt2, "Wild-type Freq (p_plus):", .fmt_f(1 - pd, 4),
                    "MOI:", MOI))

    message(sprintf(fmt2, "e01 (hom->het):", .fmt_f(e01, 4),
                    "e02 (het->hom):", .fmt_f(e02, 4)))
    message(sprintf("%-34s %12s", "e03 (hom->hom):", .fmt_f(e03, 4)))

    message("-------------------------------------------------------------------")
    message("Non-Centrality Parameter (lambda) and Power")
    message(sprintf("  %-16s %12.5f  | df=%d", "lambda:", lambda, 2))
    message(sprintf("  %-16s %12.6f", "power:", power_out))
    message("-------------------------------------------------------------------")
  }

  invisible(out)
}





# Differential genotype misclass
# -------------------------------------------------------------------



cc_misclass_matrix_3p <- function(e01 = 0, e02 = 0, e03 = 0) {
  if (!is.numeric(e01) || length(e01) != 1 || e01 < 0 || e01 > 1)
    stop("e01 must be a single number in [0,1].")
  if (!is.numeric(e02) || length(e02) != 1 || e02 < 0 || e02 > 0.5)
    stop("e02 must be a single number in [0,0.5].")
  if (!is.numeric(e03) || length(e03) != 1 || e03 < 0 || e03 > 1)
    stop("e03 must be a single number in [0,1].")
  if (e01 + e03 > 1)
    stop("Need e01 + e03 <= 1 so homozygote rows remain nonnegative.")

  M <- matrix(c(
    1 - (e01 + e03),  e01,            e03,
    e02,              1 - 2 * e02,    e02,
    e03,              e01,            1 - (e01 + e03)
  ), nrow = 3, byrow = TRUE)

  if (any(abs(rowSums(M) - 1) > 1e-10))
    stop("Misclassification matrix rows do not sum to 1.")
  if (any(M < -1e-12))
    stop("Misclassification matrix has negative entries; check e01/e02/e03.")

  M
}


cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
  if (!is.numeric(g_true) || length(g_true) != 3)
    stop("g_true must be a numeric vector of length 3.")
  if (any(g_true < -1e-12))
    stop("g_true has negative entries.")
  if (abs(sum(g_true) - 1) > 1e-8)
    stop("g_true must sum to 1.")
  if (!is.matrix(M_true_to_obs) || any(dim(M_true_to_obs) != c(3, 3)))
    stop("M_true_to_obs must be a 3 x 3 matrix.")

  g_obs <- as.numeric(t(M_true_to_obs) %*% g_true)
  g_obs / sum(g_obs)
}
#' @rdname case_control_genotype_misclassification
#' @export
cc_chisq_mssn_genotype_misclassification_differential_3p <- function(
    power, alpha, prev, pd, R2,
    MOI = c("M", "D", "Rec"),
    k = 1,

    # case misclassification parameters
    case_e01 = 0, case_e02 = 0, case_e03 = 0,

    # control misclassification parameters
    ctrl_e01 = 0, ctrl_e02 = 0, ctrl_e03 = 0,

    verbose = TRUE
) {
  MOI <- match.arg(MOI)

  chisq_ncp_target <- function(power, alpha, df) {
    crit <- qchisq(1 - alpha, df = df)
    f <- function(lambda) pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - power
    uniroot(f, lower = 0, upper = 1e6)$root
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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)

  # checks
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("power must be a single number in (0,1).")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")

  lambda_star <- chisq_ncp_target(power = power, alpha = alpha, df = 2)

  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0

  M_case <- cc_misclass_matrix_3p(case_e01, case_e02, case_e03)
  M_ctrl <- cc_misclass_matrix_3p(ctrl_e01, ctrl_e02, ctrl_e03)

  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M_case)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M_ctrl)

  # Eq 1.22 (k-form)
  comp <- (g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl)
  S <- sum(comp)

  if (!is.finite(S) || S <= 0)
    stop("Genotype S <= 0; check inputs or differential misclassification settings.")

  MSSN_case <- ceiling(lambda_star / (k * S))
  MSSN_ctrl <- ceiling(k * MSSN_case)

  out <- list(
    alpha = alpha,
    target_power = power,
    df = 2,
    k = k,
    lambda_star = lambda_star,
    MSSN_case = MSSN_case,
    MSSN_ctrl = MSSN_ctrl,
    components = comp,
    S = S,
    model_parameters = list(
      prev = prev, pd = pd, qd = 1 - pd,
      R1 = freqs$R1, R2 = R2, MOI = MOI
    ),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    misclassification = list(
      model = "differential_3p",
      case_matrix = M_case,
      ctrl_matrix = M_ctrl,
      case_params = c(e01 = case_e01, e02 = case_e02, e03 = case_e03),
      ctrl_params = c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)
    ),
    freqs = list(
      g_true_case = g_true_case,
      g_true_ctrl = g_true_ctrl,
      g_obs_case = g_obs_case,
      g_obs_ctrl = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_mssn_genotype_misclassification_differential_3p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): Minimum Sample Size Necessary (MSSN) ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df = 2")
    message("Misclassification: DIFFERENTIAL 3-parameter model (separate case/control matrices)")
    message("-----------------------------------------------------------------------------------")

    fmt2 <- "%-32s %12s  |  %-28s %12s"
    message(sprintf(fmt2,
                    "Target Power:",
                    .fmt_f(power, 3),
                    "Significance Level (alpha):",
                    .fmt_e(alpha, 2)))

    message(sprintf(fmt2,
                    "Disease Prevalence (prev):",
                    .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):",
                    .fmt_f(pd, 4)))

    message(sprintf(fmt2,
                    "Wild-type Freq (p_plus):",
                    .fmt_f(1 - pd, 4),
                    "MOI:",
                    MOI))

    message(sprintf(fmt2,
                    "Case:Control Ratio (k):",
                    .fmt_f(k, 3),
                    "Target NCP (lambda*):",
                    formatC(lambda_star, format = "f", digits = 5)))

    message("-----------------------------------------------------------------------------------")
    message("Case genotype error parameters")
    message(sprintf("  e01=%0.4f   e02=%0.4f   e03=%0.4f", case_e01, case_e02, case_e03))
    message("Control genotype error parameters")
    message(sprintf("  e01=%0.4f   e02=%0.4f   e03=%0.4f", ctrl_e01, ctrl_e02, ctrl_e03))

    message("-----------------------------------------------------------------------------------")
    message("Observed genotype frequencies (cases vs controls)")
    message(sprintf("  g0~: %10.6f vs %10.6f", g_obs_case[1], g_obs_ctrl[1]))
    message(sprintf("  g1~: %10.6f vs %10.6f", g_obs_case[2], g_obs_ctrl[2]))
    message(sprintf("  g2~: %10.6f vs %10.6f", g_obs_case[3], g_obs_ctrl[3]))

    message("-----------------------------------------------------------------------------------")
    message("Eq 1.22 components")
    message(sprintf("  Component 01: %.10g", comp[1]))
    message(sprintf("  Component 02: %.10g", comp[2]))
    message(sprintf("  Component 03: %.10g", comp[3]))
    message(sprintf("  Sum(components): %.10g", S))

    message("-----------------------------------------------------------------------------------")
    message(sprintf("MSSN: N_case = %d  |  N_ctrl = %d", MSSN_case, MSSN_ctrl))
    message("-----------------------------------------------------------------------------------")
  }

  invisible(out)
}
#' @rdname case_control_genotype_misclassification
#' @export
cc_chisq_power_genotype_misclassification_differential_3p <- function(
    N_case, alpha, prev, pd, R2,
    MOI = c("M", "D", "Rec"),
    k = 1,

    case_e01 = 0, case_e02 = 0, case_e03 = 0,
    ctrl_e01 = 0, ctrl_e02 = 0, ctrl_e03 = 0,

    verbose = TRUE
) {
  MOI <- match.arg(MOI)

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

    list(R1 = R1, f0 = f0, f1 = f1, f2 = f2, g_j1 = g_j1, g_j0 = g_j0)
  }

  .fmt_f <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
  .fmt_e <- function(x, digits = 2) formatC(x, format = "e", digits = digits)

  if (!is.numeric(N_case) || length(N_case) != 1 || N_case <= 0)
    stop("N_case must be a single positive number.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")
  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1)
    stop("prev must be a single number in (0,1).")
  if (!is.numeric(pd) || length(pd) != 1 || pd <= 0 || pd >= 1)
    stop("pd must be a single number in (0,1).")
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0)
    stop("R2 must be a single positive number.")
  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number (N_ctrl / N_case).")

  N_ctrl <- k * N_case

  freqs <- cc_conditional_geno_freqs(prev = prev, pd = pd, R2 = R2, MOI = MOI)
  g_true_case <- freqs$g_j1
  g_true_ctrl <- freqs$g_j0

  M_case <- cc_misclass_matrix_3p(case_e01, case_e02, case_e03)
  M_ctrl <- cc_misclass_matrix_3p(ctrl_e01, ctrl_e02, ctrl_e03)

  g_obs_case <- cc_apply_genotype_misclass(g_true_case, M_case)
  g_obs_ctrl <- cc_apply_genotype_misclass(g_true_ctrl, M_ctrl)

  comp <- (g_obs_case - g_obs_ctrl)^2 / (g_obs_case + k * g_obs_ctrl)
  S <- sum(comp)

  if (!is.finite(S) || S <= 0)
    stop("Genotype S <= 0; check inputs or differential misclassification settings.")

  lambda <- k * N_case * S
  crit <- qchisq(1 - alpha, df = 2)
  power <- pchisq(crit, df = 2, ncp = lambda, lower.tail = FALSE)

  out <- list(
    alpha = alpha,
    df = 2,
    N_case = N_case,
    N_ctrl = N_ctrl,
    k = k,
    lambda = lambda,
    power = power,
    components = comp,
    S = S,
    model_parameters = list(
      prev = prev, pd = pd, qd = 1 - pd,
      R1 = freqs$R1, R2 = R2, MOI = MOI
    ),
    penetrances = c(f0 = freqs$f0, f1 = freqs$f1, f2 = freqs$f2),
    misclassification = list(
      model = "differential_3p",
      case_matrix = M_case,
      ctrl_matrix = M_ctrl,
      case_params = c(e01 = case_e01, e02 = case_e02, e03 = case_e03),
      ctrl_params = c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)
    ),
    freqs = list(
      g_true_case = g_true_case,
      g_true_ctrl = g_true_ctrl,
      g_obs_case = g_obs_case,
      g_obs_ctrl = g_obs_ctrl
    )
  )
  class(out) <- "cc_chisq_power_genotype_misclassification_differential_3p"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control (Conditional Model): Power with Differential Genotype Misclassification ---")
    message("Test: Genotype chi-square (Eq 1.22, k-form), df = 2")
    message("-----------------------------------------------------------------------------------")

    fmt2 <- "%-32s %12s  |  %-28s %12s"
    message(sprintf(fmt2,
                    "N_case:",
                    formatC(N_case, format = "f", digits = 0),
                    "N_ctrl:",
                    formatC(N_ctrl, format = "f", digits = 0)))

    message(sprintf(fmt2,
                    "Significance Level (alpha):",
                    .fmt_e(alpha, 2),
                    "Case:Control Ratio (k):",
                    .fmt_f(k, 3)))

    message(sprintf(fmt2,
                    "Disease Prevalence (prev):",
                    .fmt_f(prev, 4),
                    "Risk Allele Freq (p_d):",
                    .fmt_f(pd, 4)))

    message(sprintf(fmt2,
                    "Wild-type Freq (p_plus):",
                    .fmt_f(1 - pd, 4),
                    "MOI:",
                    MOI))

    message("-----------------------------------------------------------------------------------")
    message("Case genotype error parameters")
    message(sprintf("  e01=%0.4f   e02=%0.4f   e03=%0.4f", case_e01, case_e02, case_e03))
    message("Control genotype error parameters")
    message(sprintf("  e01=%0.4f   e02=%0.4f   e03=%0.4f", ctrl_e01, ctrl_e02, ctrl_e03))

    message("-----------------------------------------------------------------------------------")
    message("Observed genotype frequencies (cases vs controls)")
    message(sprintf("  g0~: %10.6f vs %10.6f", g_obs_case[1], g_obs_ctrl[1]))
    message(sprintf("  g1~: %10.6f vs %10.6f", g_obs_case[2], g_obs_ctrl[2]))
    message(sprintf("  g2~: %10.6f vs %10.6f", g_obs_case[3], g_obs_ctrl[3]))

    message("-----------------------------------------------------------------------------------")
    message(sprintf("Sum(components): %.10g", S))
    message(sprintf("lambda: %.10f", lambda))
    message(sprintf("power:  %.6f", power))
    message("-----------------------------------------------------------------------------------")
  }

  invisible(out)
}



## helper function

# locus het

cc_chisq_ncp_target <- function(power, alpha, df) {
  crit <- qchisq(1 - alpha, df = df)
  f <- function(lambda) {
    pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - power
  }
  uniroot(f, lower = 0, upper = 1e6)$root
}


cc_check_genotype_freqs <- function(g, name = "g") {
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


cc_apply_locus_het <- function(g_case_assoc, g_ctrl, pi) {

  cc_check_genotype_freqs(g_case_assoc, "g_case_assoc")
  cc_check_genotype_freqs(g_ctrl, "g_ctrl")

  if (!is.numeric(pi) || length(pi) != 1 || pi < 0 || pi > 1)
    stop("pi must be a single number in [0,1].")

  g_case_het <- pi * g_case_assoc + (1 - pi) * g_ctrl
  g_ctrl_het <- g_ctrl

  list(
    pi = pi,
    g_case_assoc = g_case_assoc,
    g_ctrl = g_ctrl,
    g_case_het = g_case_het,
    g_ctrl_het = g_ctrl_het
  )
}


cc_fmt_f <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}


cc_fmt_e <- function(x, digits = 2) {
  formatC(x, format = "e", digits = digits)
}




#' Specialized Case-Control Locus-Heterogeneity Functions
#'
#' Convenience functions for case-control association calculations under locus
#' heterogeneity. These focused functions compute genotype chi-square or
#' genotype trend calculations one test at a time. The all-in-one case-control
#' functions remain the recommended interface when multiple tests or modifiers
#' are needed in one call.
#'
#' @param power Numeric in \eqn{(0,1)}. Desired target power for MSSN functions.
#' @param N_case Numeric \eqn{> 0}. Number of cases for power functions.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param g_case_assoc Numeric vector of length 3. Case genotype frequencies
#'   under the associated-locus model, ordered as \code{c(g0, g1, g2)} and
#'   summing to 1.
#' @param g_ctrl Numeric vector of length 3. Control genotype frequencies,
#'   ordered as \code{c(g0, g1, g2)} and summing to 1.
#' @param pi Numeric in \eqn{[0,1]}. Locus-homogeneity fraction. The adjusted
#'   case genotype frequencies are
#'   \eqn{\pi g_{case,assoc} + (1 - \pi) g_{ctrl}}.
#' @param k Numeric \eqn{> 0}. Control-to-case ratio \eqn{N_{ctrl}/N_{case}}.
#' @param w Numeric vector of length 3. Trend-test genotype scores. The three
#'   weights cannot all be equal.
#' @param verbose Logical. If \code{TRUE}, prints a formatted summary.
#'
#' @return A list with class matching the function name. MSSN functions include
#' target non-centrality parameter, sample sizes, internal \code{S} values, and
#' adjusted frequencies. Power functions include sample sizes,
#' non-centrality parameter, power, internal \code{S} values, and adjusted
#' frequencies. Trend functions additionally return the trend numerator and
#' denominator.
#'
#' @examples
#' cc_chisq_mssn_locus_heterogeneity(
#'   power = 0.8, alpha = 0.05,
#'   g_case_assoc = c(0.25, 0.50, 0.25),
#'   g_ctrl = c(0.36, 0.48, 0.16),
#'   pi = 0.8,
#'   verbose = FALSE
#' )
#'
#' cc_trend_power_locus_heterogeneity(
#'   N_case = 500, alpha = 0.05,
#'   g_case_assoc = c(0.25, 0.50, 0.25),
#'   g_ctrl = c(0.36, 0.48, 0.16),
#'   pi = 0.8,
#'   verbose = FALSE
#' )
#'
#' @references
#' Gordon, D., Finch, S. J., & Nothnagel, M. (2020).
#' *Heterogeneity in Statistical Genetics*. Springer Nature.
#'
#' TODO: Provide exact textbook equation numbers/pages for the specialized
#' case-control locus-heterogeneity calculations.
#'
#' @name case_control_locus_heterogeneity
#' @importFrom stats pchisq qchisq uniroot
#' @export
cc_chisq_mssn_locus_heterogeneity <- function(
    power,
    alpha,
    g_case_assoc,
    g_ctrl,
    pi,
    k = 1,
    verbose = TRUE
) {

  # checks
  if (!is.numeric(power) || length(power) != 1 || power <= 0 || power >= 1)
    stop("power must be a single number in (0,1).")

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")

  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number: N_ctrl / N_case.")

  # apply locus heterogeneity
  het <- cc_apply_locus_het(
    g_case_assoc = g_case_assoc,
    g_ctrl = g_ctrl,
    pi = pi
  )

  g1 <- het$g_case_het
  g0 <- het$g_ctrl_het

  # target NCP, df = 2 for genotype chi-square
  lambda_star <- cc_chisq_ncp_target(
    power = power,
    alpha = alpha,
    df = 2
  )

  # genotype chi-square S component
  S_g <- sum((g1 - g0)^2 / (g1 + k * g0))

  if (!is.finite(S_g) || S_g <= 0)
    stop("Genotype S_g <= 0; check inputs or pi.")

  N_case <- ceiling(lambda_star / (k * S_g))
  N_ctrl <- ceiling(k * N_case)

  out <- list(
    test = "case-control chi-square test of independence for genotypes",
    df = 2,
    alpha = alpha,
    target_power = power,
    k = k,
    pi = pi,
    lambda_star = lambda_star,
    S = S_g,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    freqs = list(
      g_case_assoc = het$g_case_assoc,
      g_ctrl = het$g_ctrl,
      g_case_het = g1,
      g_ctrl_het = g0
    )
  )

  class(out) <- "cc_chisq_mssn_locus_heterogeneity"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control Locus Heterogeneity: MSSN for Genotype Chi-Square ---")
    message("-----------------------------------------------------------------------")
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "Target Power:", cc_fmt_f(power, 3),
                    "Significance Level (alpha):", cc_fmt_e(alpha, 2)))
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "Locus Heterogeneity (pi):", cc_fmt_f(pi, 3),
                    "Case:Control Ratio (k):", cc_fmt_f(k, 3)))
    message("-----------------------------------------------------------------------")
    message(sprintf("%-32s %12.5f", "Target NCP (lambda_star):", lambda_star))
    message("-----------------------------------------------------------------------")
    message("Required Sample Sizes")
    message(sprintf("  %-16s %8d", "N_case:", N_case))
    message(sprintf("  %-16s %8d", "N_ctrl:", N_ctrl))
    message(sprintf("  %-16s %8d", "N_total:", N_case + N_ctrl))
    message("-----------------------------------------------------------------------")
    message("Genotype frequencies after locus heterogeneity adjustment")
    message(sprintf("  g0: %6.3f vs %6.3f", g1[1], g0[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g1[2], g0[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g1[3], g0[3]))
    message("-----------------------------------------------------------------------")
  }

  invisible(out)
}

## POWER
#' @rdname case_control_locus_heterogeneity
#' @export
cc_chisq_power_locus_heterogeneity <- function(
    N_case,
    alpha,
    g_case_assoc,
    g_ctrl,
    pi,
    k = 1,
    verbose = TRUE
) {

  # checks
  if (!is.numeric(N_case) || length(N_case) != 1 || N_case <= 0)
    stop("N_case must be a single positive number.")

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
    stop("alpha must be a single number in (0,1).")

  if (!is.numeric(k) || length(k) != 1 || k <= 0)
    stop("k must be a single positive number: N_ctrl / N_case.")

  N_ctrl <- k * N_case

  # apply locus heterogeneity
  het <- cc_apply_locus_het(
    g_case_assoc = g_case_assoc,
    g_ctrl = g_ctrl,
    pi = pi
  )

  g1 <- het$g_case_het
  g0 <- het$g_ctrl_het

  # genotype chi-square NCP
  S_g <- sum((g1 - g0)^2 / (g1 + k * g0))

  if (!is.finite(S_g) || S_g <= 0)
    stop("Genotype S_g <= 0; check inputs or pi.")

  lambda_g <- k * N_case * S_g

  crit <- qchisq(1 - alpha, df = 2)
  power_g <- pchisq(crit, df = 2, ncp = lambda_g, lower.tail = FALSE)

  out <- list(
    test = "case-control chi-square test of independence for genotypes",
    df = 2,
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    k = k,
    pi = pi,
    lambda = lambda_g,
    S = S_g,
    power = power_g,
    freqs = list(
      g_case_assoc = het$g_case_assoc,
      g_ctrl = het$g_ctrl,
      g_case_het = g1,
      g_ctrl_het = g0
    )
  )

  class(out) <- "cc_chisq_power_locus_heterogeneity"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control Locus Heterogeneity: Power for Genotype Chi-Square ---")
    message("----------------------------------------------------------------------")
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "N_case:", cc_fmt_f(N_case, 0),
                    "N_ctrl:", cc_fmt_f(N_ctrl, 0)))
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "Significance Level (alpha):", cc_fmt_e(alpha, 2),
                    "Case:Control Ratio (k):", cc_fmt_f(k, 3)))
    message(sprintf("%-32s %12s", "Locus Heterogeneity (pi):", cc_fmt_f(pi, 3)))
    message("----------------------------------------------------------------------")
    message(sprintf("%-32s %12.5f", "NCP (lambda):", lambda_g))
    message(sprintf("%-32s %12.6f", "Power:", power_g))
    message("----------------------------------------------------------------------")
    message("Genotype frequencies after locus heterogeneity adjustment")
    message(sprintf("  g0: %6.3f vs %6.3f", g1[1], g0[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g1[2], g0[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g1[3], g0[3]))
    message("----------------------------------------------------------------------")
  }

  invisible(out)
}













#' @rdname case_control_locus_heterogeneity
#' @export
cc_trend_mssn_locus_heterogeneity <- function(
    power,
    alpha,
    g_case_assoc,
    g_ctrl,
    pi,
    k = 1,
    w = c(0, 1, 2),
    verbose = TRUE
) {

  # checks
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

  # apply locus heterogeneity
  het <- cc_apply_locus_het(
    g_case_assoc = g_case_assoc,
    g_ctrl = g_ctrl,
    pi = pi
  )

  g1 <- het$g_case_het
  g0 <- het$g_ctrl_het

  # target NCP, df = 1 for trend test
  lambda_star <- cc_chisq_ncp_target(
    power = power,
    alpha = alpha,
    df = 1
  )

  # trend S component
  num_t <- (sum(w * (g1 - g0)))^2

  den_t <- sum(w^2 * (g1 + k * g0)) -
    (sum(w * (g1 + k * g0)))^2 / (1 + k)

  if (!is.finite(den_t) || den_t <= 0)
    stop("Trend denominator <= 0; check inputs/weights.")

  if (num_t < 1e-15)
    stop("Trend numerator is approximately 0; implies no weighted mean difference.")

  S_t <- num_t / den_t

  N_case <- ceiling(lambda_star / (k * S_t))
  N_ctrl <- ceiling(k * N_case)

  out <- list(
    test = "trend test for genotypes in the presence of locus heterogeneity",
    df = 1,
    alpha = alpha,
    target_power = power,
    k = k,
    pi = pi,
    w = w,
    lambda_star = lambda_star,
    S = S_t,
    numerator = num_t,
    denominator = den_t,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    freqs = list(
      g_case_assoc = het$g_case_assoc,
      g_ctrl = het$g_ctrl,
      g_case_het = g1,
      g_ctrl_het = g0
    )
  )

  class(out) <- "cc_trend_mssn_locus_heterogeneity"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control Locus Heterogeneity: MSSN for Trend Test ---")
    message("-------------------------------------------------------------")
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "Target Power:", cc_fmt_f(power, 3),
                    "Significance Level (alpha):", cc_fmt_e(alpha, 2)))
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "Locus Heterogeneity (pi):", cc_fmt_f(pi, 3),
                    "Case:Control Ratio (k):", cc_fmt_f(k, 3)))
    message(sprintf("%-32s %12s", "Trend Weights (w):", paste0(w, collapse = ",")))
    message("-------------------------------------------------------------")
    message(sprintf("%-32s %12.5f", "Target NCP (lambda_star):", lambda_star))
    message("-------------------------------------------------------------")
    message("Required Sample Sizes")
    message(sprintf("  %-16s %8d", "N_case:", N_case))
    message(sprintf("  %-16s %8d", "N_ctrl:", N_ctrl))
    message(sprintf("  %-16s %8d", "N_total:", N_case + N_ctrl))
    message("-------------------------------------------------------------")
    message("Genotype frequencies after locus heterogeneity adjustment")
    message(sprintf("  g0: %6.3f vs %6.3f", g1[1], g0[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g1[2], g0[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g1[3], g0[3]))
    message("-------------------------------------------------------------")
  }

  invisible(out)
}

#' @rdname case_control_locus_heterogeneity
#' @export
cc_trend_power_locus_heterogeneity <- function(
    N_case,
    alpha,
    g_case_assoc,
    g_ctrl,
    pi,
    k = 1,
    w = c(0, 1, 2),
    verbose = TRUE
) {

  # checks
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

  N_ctrl <- k * N_case

  # apply locus heterogeneity
  het <- cc_apply_locus_het(
    g_case_assoc = g_case_assoc,
    g_ctrl = g_ctrl,
    pi = pi
  )

  g1 <- het$g_case_het
  g0 <- het$g_ctrl_het

  # trend NCP
  num_t <- (sum(w * (g1 - g0)))^2

  den_t <- sum(w^2 * (g1 + k * g0)) -
    (sum(w * (g1 + k * g0)))^2 / (1 + k)

  if (!is.finite(den_t) || den_t <= 0)
    stop("Trend denominator <= 0; check inputs/weights.")

  if (num_t < 1e-15)
    stop("Trend numerator is approximately 0; implies no weighted mean difference.")

  S_t <- num_t / den_t
  lambda_t <- k * N_case * S_t

  crit <- qchisq(1 - alpha, df = 1)
  power_t <- pchisq(crit, df = 1, ncp = lambda_t, lower.tail = FALSE)

  out <- list(
    test = "trend test for genotypes in the presence of locus heterogeneity",
    df = 1,
    alpha = alpha,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl,
    k = k,
    pi = pi,
    w = w,
    lambda = lambda_t,
    S = S_t,
    numerator = num_t,
    denominator = den_t,
    power = power_t,
    freqs = list(
      g_case_assoc = het$g_case_assoc,
      g_ctrl = het$g_ctrl,
      g_case_het = g1,
      g_ctrl_het = g0
    )
  )

  class(out) <- "cc_trend_power_locus_heterogeneity"

  if (isTRUE(verbose)) {
    message("\n--- Case-Control Locus Heterogeneity: Power for Trend Test ---")
    message("--------------------------------------------------------------")
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "N_case:", cc_fmt_f(N_case, 0),
                    "N_ctrl:", cc_fmt_f(N_ctrl, 0)))
    message(sprintf("%-32s %12s  |  %-28s %12s",
                    "Significance Level (alpha):", cc_fmt_e(alpha, 2),
                    "Case:Control Ratio (k):", cc_fmt_f(k, 3)))
    message(sprintf("%-32s %12s", "Locus Heterogeneity (pi):", cc_fmt_f(pi, 3)))
    message(sprintf("%-32s %12s", "Trend Weights (w):", paste0(w, collapse = ",")))
    message("--------------------------------------------------------------")
    message(sprintf("%-32s %12.5f", "NCP (lambda):", lambda_t))
    message(sprintf("%-32s %12.6f", "Power:", power_t))
    message("--------------------------------------------------------------")
    message("Genotype frequencies after locus heterogeneity adjustment")
    message(sprintf("  g0: %6.3f vs %6.3f", g1[1], g0[1]))
    message(sprintf("  g1: %6.3f vs %6.3f", g1[2], g0[2]))
    message(sprintf("  g2: %6.3f vs %6.3f", g1[3], g0[3]))
    message("--------------------------------------------------------------")
  }

  invisible(out)
}

#' Specialized Case-Control Phenotype Misclassification Functions
#'
#' Convenience functions for genotype chi-square case-control calculations in
#' the presence of phenotype misclassification only. These focused functions are
#' narrower than \code{\link{cc_mssn}} and
#' \code{\link{cc_power}}: they use true affected and true
#' unaffected genotype frequencies supplied directly by the user and do not
#' apply locus heterogeneity or genotype misclassification.
#'
#' @param N_case Numeric \eqn{> 0}. Number of cases for power calculations.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param g_aff Numeric vector of true affected genotype frequencies, ordered
#'   consistently across affected and unaffected groups and summing to 1.
#' @param g_unaff Numeric vector of true unaffected genotype frequencies,
#'   ordered consistently with \code{g_aff} and summing to 1.
#' @param prev Numeric in \eqn{(0,1)}. Disease prevalence.
#' @param theta Numeric in \eqn{[0,1)}. Probability that a truly affected
#'   individual is classified as a control,
#'   \code{Pr(affected -> control)}.
#' @param phi Numeric in \eqn{[0,1)}. Probability that a truly unaffected
#'   individual is classified as a case,
#'   \code{Pr(unaffected -> case)}.
#' @param k Numeric \eqn{> 0}. Control-to-case ratio \eqn{N_{ctrl}/N_{case}}.
#' @param target_power Numeric in \eqn{(0,1)}. Desired target power for MSSN
#'   calculations.
#'
#' @details
#' Phenotype misclassification is applied by mixing the true affected genotype
#' distribution \code{g_aff} and true unaffected genotype distribution
#' \code{g_unaff} using disease prevalence and the phenotype-error
#' probabilities. Here \code{theta = Pr(affected -> control)} and
#' \code{phi = Pr(unaffected -> case)}. These specialized functions apply only
#' phenotype misclassification; use the full case-control functions for
#' combined phenotype misclassification, locus heterogeneity, genotype
#' misclassification, or trend tests.
#'
#' @return
#' \code{cc_chisq_power_phenotype_misclassification()} returns a list containing sample sizes,
#' phenotype-error parameters, observed case/control genotype frequencies,
#' internal \code{S}, non-centrality parameter \code{lambda}, and power.
#' \code{cc_chisq_mssn_phenotype_misclassification()} returns a list containing target power,
#' phenotype-error parameters, observed case/control genotype frequencies,
#' internal \code{S}, target non-centrality parameter \code{lambda_star}, and
#' case/control/total MSSN.
#'
#' @examples
#' g_aff <- c((1 - 0.05)^2, 2 * 0.05 * (1 - 0.05), 0.05^2)
#' g_unaff <- c((1 - 0.15)^2, 2 * 0.15 * (1 - 0.15), 0.15^2)
#'
#' cc_chisq_power_phenotype_misclassification(
#'   N_case = 250, alpha = 0.01,
#'   g_aff = g_aff, g_unaff = g_unaff,
#'   prev = 0.05, theta = 0, phi = 0.01
#' )
#'
#' cc_chisq_mssn_phenotype_misclassification(
#'   target_power = 0.8, alpha = 0.01,
#'   g_aff = g_aff, g_unaff = g_unaff,
#'   prev = 0.05, theta = 0, phi = 0.01
#' )
#'
#' @references
#' Edwards, B. J., Haynes, C., Levenstien, M. A., Finch, S. J., & Gordon, D.
#' (2005). Power and sample size calculations in the presence of phenotype
#' errors for case/control genetic association studies. \emph{BMC Genetics},
#' 6, 18.
#'
#' @name case_control_phenotype_misclassification
#' @importFrom stats pchisq qchisq uniroot
#' @export
cc_chisq_power_phenotype_misclassification <- function(
    N_case,
    alpha = 0.05,
    g_aff,
    g_unaff,
    prev,
    theta = 0,
    phi = 0,
    k = 1
) {
  # g_aff = true affected genotype frequencies
  # g_unaff = true unaffected genotype frequencies
  # prev = disease prevalence K
  # theta = Pr(true affected classified as control)
  # phi = Pr(true unaffected classified as case)
  # k = N_ctrl / N_case

  if (length(g_aff) != length(g_unaff)) {
    stop("g_aff and g_unaff must have the same length.")
  }
  if (abs(sum(g_aff) - 1) > 1e-6) stop("g_aff must sum to 1.")
  if (abs(sum(g_unaff) - 1) > 1e-6) stop("g_unaff must sum to 1.")
  if (prev <= 0 || prev >= 1) stop("prev must be in (0,1).")
  if (theta < 0 || theta >= 1) stop("theta must be in [0,1).")
  if (phi < 0 || phi >= 1) stop("phi must be in [0,1).")
  if (N_case <= 0) stop("N_case must be positive.")
  if (k <= 0) stop("k must be positive.")

  N_ctrl <- k * N_case
  df <- length(g_aff) - 1

  case_denom <- (1 - theta) * prev + phi * (1 - prev)
  ctrl_denom <- theta * prev + (1 - phi) * (1 - prev)

  g_case_obs <- (g_aff * (1 - theta) * prev + g_unaff * phi * (1 - prev)) /
    case_denom

  g_ctrl_obs <- (g_aff * theta * prev + g_unaff * (1 - phi) * (1 - prev)) /
    ctrl_denom

  # Mitra-style noncentrality parameter for 2 x n case-control table
  S <- sum((g_case_obs - g_ctrl_obs)^2 / (g_case_obs + k * g_ctrl_obs))
  lambda <- k * N_case * S

  crit <- qchisq(1 - alpha, df = df)
  power <- pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE)

  list(
    N_case = N_case,
    N_ctrl = N_ctrl,
    alpha = alpha,
    df = df,
    prev = prev,
    theta = theta,
    phi = phi,
    g_case_obs = g_case_obs,
    g_ctrl_obs = g_ctrl_obs,
    S = S,
    lambda = lambda,
    power = power
  )
}


#' @rdname case_control_phenotype_misclassification
#' @export
cc_chisq_mssn_phenotype_misclassification <- function(
    target_power = 0.80,
    alpha = 0.05,
    g_aff,
    g_unaff,
    prev,
    theta = 0,
    phi = 0,
    k = 1
) {
  if (length(g_aff) != length(g_unaff)) {
    stop("g_aff and g_unaff must have the same length.")
  }
  if (abs(sum(g_aff) - 1) > 1e-6) stop("g_aff must sum to 1.")
  if (abs(sum(g_unaff) - 1) > 1e-6) stop("g_unaff must sum to 1.")
  if (prev <= 0 || prev >= 1) stop("prev must be in (0,1).")
  if (theta < 0 || theta >= 1) stop("theta must be in [0,1).")
  if (phi < 0 || phi >= 1) stop("phi must be in [0,1).")
  if (target_power <= 0 || target_power >= 1) {
    stop("target_power must be in (0,1).")
  }
  if (k <= 0) stop("k must be positive.")

  df <- length(g_aff) - 1

  case_denom <- (1 - theta) * prev + phi * (1 - prev)
  ctrl_denom <- theta * prev + (1 - phi) * (1 - prev)

  g_case_obs <- (g_aff * (1 - theta) * prev + g_unaff * phi * (1 - prev)) /
    case_denom

  g_ctrl_obs <- (g_aff * theta * prev + g_unaff * (1 - phi) * (1 - prev)) /
    ctrl_denom

  S <- sum((g_case_obs - g_ctrl_obs)^2 / (g_case_obs + k * g_ctrl_obs))

  if (!is.finite(S) || S <= 0) {
    stop("S is not positive; there may be no detectable genotype-frequency difference.")
  }

  crit <- qchisq(1 - alpha, df = df)

  f <- function(lambda) {
    pchisq(crit, df = df, ncp = lambda, lower.tail = FALSE) - target_power
  }

  lambda_star <- uniroot(f, lower = 0, upper = 1e6)$root

  N_case <- ceiling(lambda_star / (k * S))
  N_ctrl <- ceiling(k * N_case)

  list(
    target_power = target_power,
    alpha = alpha,
    df = df,
    prev = prev,
    theta = theta,
    phi = phi,
    g_case_obs = g_case_obs,
    g_ctrl_obs = g_ctrl_obs,
    S = S,
    lambda_star = lambda_star,
    N_case = N_case,
    N_ctrl = N_ctrl,
    N_total = N_case + N_ctrl
  )
}
