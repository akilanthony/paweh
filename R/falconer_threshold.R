#' Threshold-Selected Parameters Under the Falconer Model
#'
#' Converts the genotype-specific Falconer trait distributions into upper-tail
#' case and lower-tail control penetrances and conditional genotype frequencies.
#'
#' @inheritParams qtl_falconer_parameters
#' @param threshold_mode Either \code{"percentile"} or \code{"direct"}.
#' @param x_upper Percentage selected from the upper population tail in
#'   percentile mode.
#' @param x_lower Percentage selected from the lower population tail in
#'   percentile mode.
#' @param upper_threshold Direct upper-tail threshold in direct mode.
#' @param lower_threshold Direct lower-tail threshold in direct mode.
#' @param verbose Logical. If \code{TRUE}, prints a formatted threshold-model
#'   summary. Set to \code{FALSE} to suppress all console output.
#'
#' @details
#' Cases are selected above the upper threshold and controls below the lower
#' threshold. The middle portion is excluded, so the affected and unaffected
#' selection events are not complements. Percentile inputs refer to standard
#' normal population percentiles before conditioning on genotype.
#'
#' @return An object of class \code{"qtl_falconer_threshold_parameters"} containing
#'   thresholds, genotype-specific penetrances, selected-population
#'   prevalences, conditional genotype frequencies, and Falconer parameters.
#'
#' @examples
#' qtl_falconer_threshold_parameters(
#'   qtl_var = 0.025, tau = 0.5, pd = 0.15,
#'   x_upper = 5, x_lower = 5
#' )
#'
#' @references
#' Gordon et al. (2020), \emph{Heterogeneity in Statistical Genetics},
#' Chapter 6, Section 6.1, Equations 6.3--6.5.
#'
#' @importFrom stats pnorm qnorm
#' @export
qtl_falconer_threshold_parameters <- function(
    qtl_var,
    tau,
    pd,
    threshold_mode = c("percentile", "direct"),
    x_upper = NULL,
    x_lower = NULL,
    upper_threshold = NULL,
    lower_threshold = NULL,
    verbose = TRUE
) {
  threshold_mode <- match.arg(threshold_mode)
  .falconer_check_flag(verbose, "verbose")
  parameters <- qtl_falconer_parameters(
    qtl_var = qtl_var, tau = tau, pd = pd, verbose = FALSE
  )

  if (identical(threshold_mode, "percentile")) {
    .falconer_check_scalar(
      x_upper, "x_upper", lower = 0, upper = 100,
      lower_open = TRUE, upper_open = TRUE
    )
    .falconer_check_scalar(
      x_lower, "x_lower", lower = 0, upper = 100,
      lower_open = TRUE, upper_open = TRUE
    )
    upper_threshold <- stats::qnorm(1 - x_upper / 100)
    lower_threshold <- stats::qnorm(x_lower / 100)
  } else {
    .falconer_check_scalar(upper_threshold, "upper_threshold")
    .falconer_check_scalar(lower_threshold, "lower_threshold")
  }

  if (lower_threshold >= upper_threshold) {
    stop(
      "lower_threshold must be less than upper_threshold so the middle is excluded.",
      call. = FALSE
    )
  }

  f_affected <- 1 - stats::pnorm(
    (upper_threshold - parameters$mu) / parameters$residual_sd
  )
  f_unaffected <- stats::pnorm(
    (lower_threshold - parameters$mu) / parameters$residual_sd
  )
  names(f_affected) <- names(f_unaffected) <- names(parameters$mu)

  prev_affected <- sum(parameters$pi * f_affected)
  prev_unaffected <- sum(parameters$pi * f_unaffected)
  if (!is.finite(prev_affected) || prev_affected <= 0 ||
      !is.finite(prev_unaffected) || prev_unaffected <= 0) {
    stop("Threshold selection probabilities must be positive and finite.", call. = FALSE)
  }

  g_affected <- parameters$pi * f_affected / prev_affected
  g_unaffected <- parameters$pi * f_unaffected / prev_unaffected
  if (any(!is.finite(c(f_affected, f_unaffected, g_affected, g_unaffected)))) {
    stop("Threshold parameter calculation produced non-finite values.", call. = FALSE)
  }
  if (abs(sum(g_affected) - 1) > 1e-10 ||
      abs(sum(g_unaffected) - 1) > 1e-10) {
    stop("Conditional genotype frequencies do not sum to 1.", call. = FALSE)
  }

  out <- structure(
    list(
      threshold_mode = threshold_mode,
      x_upper = if (identical(threshold_mode, "percentile")) x_upper else NULL,
      x_lower = if (identical(threshold_mode, "percentile")) x_lower else NULL,
      upper_threshold = upper_threshold,
      lower_threshold = lower_threshold,
      thresholds = c(lower = lower_threshold, upper = upper_threshold),
      f_affected = f_affected,
      f_unaffected = f_unaffected,
      prev_affected = prev_affected,
      prev_unaffected = prev_unaffected,
      g_affected = g_affected,
      g_unaffected = g_unaffected,
      falconer = parameters
    ),
    class = "qtl_falconer_threshold_parameters"
  )

  if (isTRUE(verbose)) {
    .falconer_print_header("Falconer Threshold-Selected Trait Parameters")

    .falconer_print_section("Input Parameters")
    .falconer_print_parameter("QTL variance", out$falconer$qtl_var, 4L)
    .falconer_print_parameter(
      "Dominance/Additivity ratio (tau)", out$falconer$tau, 4L
    )
    .falconer_print_parameter(
      "Increaser allele frequency", out$falconer$pd, 4L
    )

    .falconer_print_rule()
    .falconer_print_section("Threshold Information")
    .falconer_print_parameter("Threshold mode", out$threshold_mode)
    .falconer_print_parameter(
      "Upper percentile",
      if (is.null(out$x_upper)) "direct threshold" else paste0(
        .falconer_fmt(out$x_upper, 2L), "%"
      )
    )
    .falconer_print_parameter(
      "Lower percentile",
      if (is.null(out$x_lower)) "direct threshold" else paste0(
        .falconer_fmt(out$x_lower, 2L), "%"
      )
    )
    .falconer_print_parameter("Upper threshold", out$upper_threshold, 5L)
    .falconer_print_parameter("Lower threshold", out$lower_threshold, 5L)

    .falconer_print_rule()
    .falconer_print_model(out$falconer, "Derived Falconer Parameters")

    .falconer_print_rule()
    .falconer_print_section("Genotype-Specific Quantities")
    .falconer_print_table(list(
      list(label = "Trait mean", values = out$falconer$mu, digits = 4L),
      list(label = "Affected penetrance", values = out$f_affected, digits = 4L),
      list(label = "Unaffected penetrance", values = out$f_unaffected, digits = 4L),
      list(label = "Affected genotype freq.", values = out$g_affected, digits = 4L),
      list(label = "Unaffected genotype freq.", values = out$g_unaffected, digits = 4L)
    ))

    .falconer_print_rule()
    .falconer_print_section("Selected Prevalences")
    .falconer_print_parameter(
      "Affected prevalence", out$prev_affected, 5L
    )
    .falconer_print_parameter(
      "Unaffected prevalence", out$prev_unaffected, 5L
    )
    .falconer_print_rule()
  }

  invisible(out)
}

#' Power for a Threshold-Selected Genotype Chi-Square Test
#'
#' Calculates power for the 2-df genotype chi-square test after a quantitative
#' trait is converted into upper-tail cases and lower-tail controls.
#'
#' @param N_case Positive number of selected cases.
#' @param alpha Significance level in \eqn{(0,1)}.
#' @inheritParams qtl_falconer_threshold_parameters
#' @param k Positive control-to-case ratio \eqn{N_{control}/N_{case}}.
#' @param verbose Logical. If \code{TRUE}, prints a concise summary.
#'
#' @return An object of class \code{"qtl_threshold_chisq_power"} containing
#'   selected sample sizes, power, non-centrality parameter, internal effect
#'   component \code{S}, thresholds, penetrances, prevalences, conditional
#'   genotype frequencies, and Falconer parameters.
#'
#' @examples
#' qtl_threshold_chisq_power(
#'   N_case = 126, alpha = 0.0001,
#'   qtl_var = 0.025, tau = 0.5, pd = 0.15,
#'   x_upper = 5, x_lower = 5, verbose = FALSE
#' )
#'
#' @references
#' Gordon et al. (2020), \emph{Heterogeneity in Statistical Genetics},
#' Chapter 6, Section 6.1, Equations 6.1 and 6.3--6.5.
#'
#' @export
qtl_threshold_chisq_power <- function(
    N_case,
    alpha,
    qtl_var,
    tau,
    pd,
    x_upper,
    x_lower,
    k = 1,
    verbose = TRUE
) {
  .falconer_check_scalar(N_case, "N_case", lower = 0, lower_open = TRUE)
  .falconer_check_scalar(
    alpha, "alpha", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_scalar(k, "k", lower = 0, lower_open = TRUE)
  .falconer_check_flag(verbose, "verbose")

  threshold <- qtl_falconer_threshold_parameters(
    qtl_var = qtl_var, tau = tau, pd = pd,
    x_upper = x_upper, x_lower = x_lower, verbose = FALSE
  )

  # Reuse the package's validated genotype chi-square backend with pi = 1,
  # which applies no locus-heterogeneity mixing.
  cc <- cc_chisq_power_locus_heterogeneity(
    N_case = N_case,
    alpha = alpha,
    g_case_assoc = threshold$g_affected,
    g_ctrl = threshold$g_unaffected,
    pi = 1,
    k = k,
    verbose = FALSE
  )

  out <- list(
    test = "threshold-selected case-control chi-square test of genotypes",
    df = 2,
    alpha = alpha,
    N_case = cc$N_case,
    N_control = cc$N_ctrl,
    N_total = cc$N_total,
    k = k,
    lambda = cc$lambda,
    S = cc$S,
    power = cc$power,
    thresholds = threshold$thresholds,
    penetrances = list(
      affected = threshold$f_affected,
      unaffected = threshold$f_unaffected
    ),
    prevalences = c(
      affected = threshold$prev_affected,
      unaffected = threshold$prev_unaffected
    ),
    frequencies = list(
      affected = threshold$g_affected,
      unaffected = threshold$g_unaffected
    ),
    falconer = threshold$falconer,
    threshold_parameters = threshold
  )
  class(out) <- "qtl_threshold_chisq_power"

  if (isTRUE(verbose)) {
    .falconer_print_header(
      "Falconer Threshold-Selected Trait",
      "Genotype Chi-Square Power"
    )

    .falconer_print_section("Study Design")
    .falconer_print_parameter(
      "Number of selected cases", out$N_case,
      digits = if (out$N_case == floor(out$N_case)) 0L else 2L,
      integer = out$N_case == floor(out$N_case)
    )
    .falconer_print_parameter(
      "Number of selected controls", out$N_control,
      digits = if (out$N_control == floor(out$N_control)) 0L else 2L,
      integer = out$N_control == floor(out$N_control)
    )
    .falconer_print_parameter("Control/case ratio (k)", out$k, 3L)
    .falconer_print_parameter(
      "Total selected sample", out$N_total,
      digits = if (out$N_total == floor(out$N_total)) 0L else 2L,
      integer = out$N_total == floor(out$N_total)
    )
    .falconer_print_parameter(
      "Significance level (alpha)", out$alpha, 2L, scientific = TRUE
    )

    .falconer_print_rule()
    .falconer_print_section("Threshold Selection")
    .falconer_print_parameter(
      "Upper percentile", paste0(.falconer_fmt(x_upper, 2L), "%")
    )
    .falconer_print_parameter(
      "Lower percentile", paste0(.falconer_fmt(x_lower, 2L), "%")
    )
    .falconer_print_parameter("Upper threshold", out$thresholds["upper"], 5L)
    .falconer_print_parameter("Lower threshold", out$thresholds["lower"], 5L)

    .falconer_print_rule()
    .falconer_print_model(out$falconer)

    .falconer_print_rule()
    .falconer_print_section("Genotype-Specific Quantities")
    .falconer_print_table(list(
      list(label = "Trait mean", values = out$falconer$mu, digits = 4L),
      list(label = "Affected penetrance", values = out$penetrances$affected, digits = 4L),
      list(label = "Unaffected penetrance", values = out$penetrances$unaffected, digits = 4L),
      list(label = "Case genotype frequency", values = out$frequencies$affected, digits = 4L),
      list(label = "Control genotype frequency", values = out$frequencies$unaffected, digits = 4L)
    ))

    .falconer_print_rule()
    .falconer_print_section("Selected Population")
    .falconer_print_parameter(
      "Affected prevalence", out$prevalences["affected"], 5L
    )
    .falconer_print_parameter(
      "Unaffected prevalence", out$prevalences["unaffected"], 5L
    )

    .falconer_print_rule()
    .falconer_print_section("Test Results")
    .falconer_print_parameter("Degrees of freedom", out$df, integer = TRUE)
    .falconer_print_parameter(
      "Non-centrality parameter", out$lambda, 5L
    )
    .falconer_print_parameter("Power", out$power, 6L)
    .falconer_print_rule()
  }

  invisible(out)
}

#' Minimum Sample Size for a Threshold-Selected Genotype Chi-Square Test
#'
#' Calculates the minimum selected case and control sample sizes after deriving
#' their genotype distributions from a single-trait Falconer threshold model.
#'
#' @param power Target power in \eqn{(0,1)}.
#' @inheritParams qtl_threshold_chisq_power
#'
#' @details
#' The primary MSSN is the selected case-control sample size used by the
#' association test. The separately reported screening quantities are expected
#' population counts needed to obtain those selected samples under simple
#' population sampling; they are not statistical MSSN values.
#'
#' @return An object of class \code{"qtl_threshold_chisq_mssn"} containing
#'   selected MSSN values, target NCP, internal \code{S}, thresholds,
#'   penetrances, prevalences, conditional frequencies, Falconer parameters,
#'   and separately labelled expected population screening counts.
#'
#' @examples
#' qtl_threshold_chisq_mssn(
#'   power = 0.8, alpha = 0.0001,
#'   qtl_var = 0.025, tau = 0.5, pd = 0.15,
#'   x_upper = 5, x_lower = 5, verbose = FALSE
#' )
#'
#' @export
qtl_threshold_chisq_mssn <- function(
    power,
    alpha,
    qtl_var,
    tau,
    pd,
    x_upper,
    x_lower,
    k = 1,
    verbose = TRUE
) {
  .falconer_check_scalar(
    power, "power", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_scalar(
    alpha, "alpha", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_scalar(k, "k", lower = 0, lower_open = TRUE)
  .falconer_check_flag(verbose, "verbose")

  threshold <- qtl_falconer_threshold_parameters(
    qtl_var = qtl_var, tau = tau, pd = pd,
    x_upper = x_upper, x_lower = x_lower, verbose = FALSE
  )
  cc <- cc_chisq_mssn_locus_heterogeneity(
    power = power,
    alpha = alpha,
    g_case_assoc = threshold$g_affected,
    g_ctrl = threshold$g_unaffected,
    pi = 1,
    k = k,
    verbose = FALSE
  )

  screened_cases <- cc$N_case / threshold$prev_affected
  screened_controls <- cc$N_ctrl / threshold$prev_unaffected
  out <- list(
    test = "threshold-selected case-control chi-square test of genotypes",
    df = 2,
    target_power = power,
    alpha = alpha,
    N_case = cc$N_case,
    N_control = cc$N_ctrl,
    N_total = cc$N_total,
    k = k,
    lambda_star = cc$lambda_star,
    S = cc$S,
    expected_population_screened_cases = screened_cases,
    expected_population_screened_controls = screened_controls,
    screening = list(
      expected_population_screened_cases = screened_cases,
      expected_population_screened_controls = screened_controls
    ),
    thresholds = threshold$thresholds,
    penetrances = list(
      affected = threshold$f_affected,
      unaffected = threshold$f_unaffected
    ),
    prevalences = c(
      affected = threshold$prev_affected,
      unaffected = threshold$prev_unaffected
    ),
    frequencies = list(
      affected = threshold$g_affected,
      unaffected = threshold$g_unaffected
    ),
    falconer = threshold$falconer,
    threshold_parameters = threshold
  )
  class(out) <- "qtl_threshold_chisq_mssn"

  if (isTRUE(verbose)) {
    achieved <- qtl_threshold_chisq_power(
      N_case = out$N_case,
      alpha = out$alpha,
      qtl_var = out$falconer$qtl_var,
      tau = out$falconer$tau,
      pd = out$falconer$pd,
      x_upper = x_upper,
      x_lower = x_lower,
      k = out$k,
      verbose = FALSE
    )$power

    .falconer_print_header(
      "Falconer Threshold-Selected Trait",
      "Genotype Chi-Square Minimum Sample Size"
    )

    .falconer_print_section("Target Design")
    .falconer_print_parameter("Target power", out$target_power, 3L)
    .falconer_print_parameter(
      "Significance level (alpha)", out$alpha, 2L, scientific = TRUE
    )
    .falconer_print_parameter("Control/case ratio (k)", out$k, 3L)

    .falconer_print_rule()
    .falconer_print_section("Required Selected Sample")
    .falconer_print_parameter("Cases", out$N_case, integer = TRUE)
    .falconer_print_parameter("Controls", out$N_control, integer = TRUE)
    .falconer_print_parameter(
      "Total selected sample", out$N_total, integer = TRUE
    )
    .falconer_print_parameter("Achieved power", achieved, 6L)
    .falconer_print_parameter("Required NCP", out$lambda_star, 5L)

    .falconer_print_rule()
    .falconer_print_section("Threshold Selection")
    .falconer_print_parameter(
      "Upper percentile", paste0(.falconer_fmt(x_upper, 2L), "%")
    )
    .falconer_print_parameter(
      "Lower percentile", paste0(.falconer_fmt(x_lower, 2L), "%")
    )
    .falconer_print_parameter("Upper threshold", out$thresholds["upper"], 5L)
    .falconer_print_parameter("Lower threshold", out$thresholds["lower"], 5L)

    .falconer_print_rule()
    .falconer_print_model(out$falconer)

    .falconer_print_rule()
    .falconer_print_section("Genotype-Specific Quantities")
    .falconer_print_table(list(
      list(label = "Trait mean", values = out$falconer$mu, digits = 4L),
      list(label = "Affected penetrance", values = out$penetrances$affected, digits = 4L),
      list(label = "Unaffected penetrance", values = out$penetrances$unaffected, digits = 4L),
      list(label = "Case genotype frequency", values = out$frequencies$affected, digits = 4L),
      list(label = "Control genotype frequency", values = out$frequencies$unaffected, digits = 4L)
    ))

    .falconer_print_rule()
    .falconer_print_section("Selected Population")
    .falconer_print_parameter(
      "Affected prevalence", out$prevalences["affected"], 5L
    )
    .falconer_print_parameter(
      "Unaffected prevalence", out$prevalences["unaffected"], 5L
    )

    .falconer_print_rule()
    .falconer_print_section("Expected Screening Burden")
    .falconer_print_parameter(
      "Individuals screened for selected cases",
      out$expected_population_screened_cases,
      1L
    )
    .falconer_print_parameter(
      "Individuals screened for selected controls",
      out$expected_population_screened_controls,
      1L
    )
    message("  Screening burden is separate from the statistical selected-sample MSSN.")
    .falconer_print_rule()
  }

  invisible(out)
}
