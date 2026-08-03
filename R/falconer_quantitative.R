# Internal validation helpers for the single-trait Falconer framework.

#' Format a Value for Falconer Console Output
#' @noRd
.falconer_fmt <- function(x, digits = 4L, scientific = FALSE,
                          integer = FALSE) {
  if (is.character(x)) return(x)
  if (integer) return(formatC(as.integer(x), format = "d", big.mark = ","))
  if (!scientific && is.finite(x) && abs(x) < 0.5 * 10^(-digits)) x <- 0
  formatC(
    x,
    format = if (scientific) "e" else "f",
    digits = digits
  )
}

#' Print a Falconer Console Header
#' @noRd
.falconer_print_header <- function(title, subtitle = NULL) {
  line <- paste(rep("=", 74L), collapse = "")
  message("\n", line)
  message(title)
  if (!is.null(subtitle)) message(subtitle)
  message(line)
}

#' Print a Falconer Section Heading
#' @noRd
.falconer_print_section <- function(title) {
  message("\n", title)
}

#' Print a Falconer Section Separator
#' @noRd
.falconer_print_rule <- function() {
  message(paste(rep("-", 74L), collapse = ""))
}

#' Print an Aligned Falconer Parameter
#' @noRd
.falconer_print_parameter <- function(label, value, digits = 4L,
                                       scientific = FALSE,
                                       integer = FALSE) {
  formatted <- .falconer_fmt(
    value, digits = digits, scientific = scientific, integer = integer
  )
  message(sprintf("%-44s %14s", paste0(label, ":"), formatted))
}

#' Print an Aligned Three-Genotype Falconer Vector
#' @noRd
.falconer_print_vector <- function(label, values, digits = 4L,
                                    integer = FALSE) {
  formatted <- vapply(values, function(value) {
    .falconer_fmt(value, digits = digits, integer = integer)
  }, character(1))
  message(sprintf(
    "%-28s %14s %14s %14s",
    paste0(label, ":"), formatted[1L], formatted[2L], formatted[3L]
  ))
}

#' Print a Falconer Three-Genotype Table
#' @noRd
.falconer_print_table <- function(rows) {
  message(sprintf(
    "%-28s %14s %14s %14s",
    "", "Genotype 0", "Genotype 1", "Genotype 2"
  ))
  for (row in rows) {
    .falconer_print_vector(
      row$label,
      row$values,
      digits = if (is.null(row$digits)) 4L else row$digits,
      integer = isTRUE(row$integer)
    )
  }
}

#' Print the Shared Falconer Model Section
#' @noRd
.falconer_print_model <- function(parameters, title = "Falconer Model") {
  .falconer_print_section(title)
  .falconer_print_parameter("QTL variance", parameters$qtl_var, 4L)
  .falconer_print_parameter(
    "Residual variance", parameters$residual_variance, 4L
  )
  .falconer_print_parameter("Residual SD", parameters$residual_sd, 4L)
  .falconer_print_parameter(
    "Increaser allele frequency", parameters$pd, 4L
  )
  .falconer_print_parameter(
    "Dominance/Additivity ratio (tau)", parameters$tau, 4L
  )
  .falconer_print_parameter("Additive effect (a)", parameters$a, 4L)
  .falconer_print_parameter(
    "Dominance effect (delta)", parameters$delta, 4L
  )
  .falconer_print_parameter("Centering mean (m)", parameters$m, 4L)
}

.falconer_check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                                   lower_open = FALSE, upper_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(name, " must be a single finite number.", call. = FALSE)
  }

  lower_bad <- if (lower_open) x <= lower else x < lower
  upper_bad <- if (upper_open) x >= upper else x > upper
  if (lower_bad || upper_bad) {
    interval <- paste0(
      if (lower_open) "(" else "[", lower, ", ", upper,
      if (upper_open) ")" else "]"
    )
    stop(name, " must be in ", interval, ".", call. = FALSE)
  }

  invisible(TRUE)
}

.falconer_check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

.falconer_anova_counts <- function(N, pi, count_method) {
  expected <- N * pi
  if (identical(count_method, "rounded")) round(expected) else expected
}

.falconer_anova_components <- function(N, parameters, count_method) {
  counts <- .falconer_anova_counts(N, parameters$pi, count_method)
  if (identical(count_method, "rounded") && any(counts <= 0)) {
    return(NULL)
  }

  df1 <- 2
  df2 <- N - 3
  lambda <- sum(counts * parameters$mu^2) / parameters$residual_variance
  critical <- stats::qf(1 - parameters$alpha, df1 = df1, df2 = df2)
  achieved_power <- stats::pf(
    critical,
    df1 = df1,
    df2 = df2,
    ncp = lambda,
    lower.tail = FALSE
  )

  list(
    counts = counts,
    df1 = df1,
    df2 = df2,
    lambda = lambda,
    critical = critical,
    power = achieved_power
  )
}

#' Falconer Parameters for a Single Quantitative Trait
#'
#' Calculates genotype-specific normal-mixture parameters for a standardized
#' continuous trait under the single-locus Falconer model.
#'
#' @param qtl_var Variance in the standardized trait attributable to the QTL.
#'   Must be in \eqn{(0,1)}.
#' @param tau Finite dominance-to-additivity ratio \eqn{\delta/a}.
#' @param pd Frequency of the increaser allele. Must be in \eqn{(0,1)}.
#' @param verbose Logical. If \code{TRUE}, prints a formatted parameter
#'   summary. Set to \code{FALSE} to suppress all console output.
#'
#' @details
#' Genotypes are ordered as zero, one, and two copies of the increaser allele.
#' Their Hardy-Weinberg mixing proportions are
#' \eqn{((1-p_d)^2, 2p_d(1-p_d), p_d^2)}. The population-centering constant
#' makes the genotype-weighted trait mean zero, and the residual variance is
#' \eqn{1 - V_{QTL}}.
#'
#' @return An object of class \code{"falconer_parameters"}. The list contains
#'   \code{qtl_var}, \code{tau}, \code{pd}, \code{p_plus}, additive effect
#'   \code{a}, dominance effect \code{delta}, centering constant \code{m},
#'   genotype means \code{mu}, mixing proportions \code{pi},
#'   \code{residual_variance}, \code{residual_sd}, \code{weighted_mean}, and
#'   \code{total_variance}.
#'
#' @examples
#' falconer_parameters(qtl_var = 0.025, tau = 0.5, pd = 0.15)
#'
#' @references
#' Gordon et al. (2020), \emph{Heterogeneity in Statistical Genetics},
#' Chapter 6, Section 6.1, pp. 324--325.
#'
#' @export
falconer_parameters <- function(qtl_var, tau, pd, verbose = TRUE) {
  .falconer_check_scalar(
    qtl_var, "qtl_var", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_scalar(tau, "tau")
  .falconer_check_scalar(
    pd, "pd", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_flag(verbose, "verbose")

  p_plus <- 1 - pd
  denominator <-
    2 * pd * p_plus * (1 + tau * (p_plus - pd))^2 +
    4 * (pd * p_plus * tau)^2
  a <- sqrt(qtl_var / denominator)
  delta <- tau * a
  m <- -(pd^2 * a + 2 * pd * p_plus * delta - p_plus^2 * a)
  mu <- c(m - a, m + delta, m + a)
  pi <- c(p_plus^2, 2 * pd * p_plus, pd^2)
  names(mu) <- names(pi) <- c("genotype_0", "genotype_1", "genotype_2")

  residual_variance <- 1 - qtl_var
  residual_sd <- sqrt(residual_variance)
  weighted_mean <- sum(pi * mu)
  total_variance <- residual_variance + sum(pi * (mu - weighted_mean)^2)

  computed <- c(
    p_plus = p_plus,
    denominator = denominator,
    a = a,
    delta = delta,
    m = m,
    mu,
    pi,
    residual_variance = residual_variance,
    residual_sd = residual_sd,
    weighted_mean = weighted_mean,
    total_variance = total_variance
  )
  if (any(!is.finite(computed))) {
    stop("Falconer parameter calculation produced non-finite values.", call. = FALSE)
  }
  if (abs(sum(pi) - 1) > 1e-12) {
    stop("Falconer genotype mixing proportions do not sum to 1.", call. = FALSE)
  }
  if (abs(weighted_mean) > 1e-10) {
    stop("Falconer genotype-weighted mean is not approximately zero.", call. = FALSE)
  }
  if (abs(total_variance - 1) > 1e-10) {
    stop("Falconer total trait variance is not approximately 1.", call. = FALSE)
  }

  out <- structure(
    list(
      qtl_var = qtl_var,
      tau = tau,
      pd = pd,
      p_plus = p_plus,
      a = a,
      delta = delta,
      m = m,
      mu = mu,
      pi = pi,
      residual_variance = residual_variance,
      residual_sd = residual_sd,
      weighted_mean = weighted_mean,
      total_variance = total_variance
    ),
    class = "falconer_parameters"
  )

  if (isTRUE(verbose)) {
    .falconer_print_header("Falconer Quantitative Trait Parameters")

    .falconer_print_section("Input Parameters")
    .falconer_print_parameter("QTL variance", out$qtl_var, 4L)
    .falconer_print_parameter(
      "Dominance/Additivity ratio (tau)", out$tau, 4L
    )
    .falconer_print_parameter("Increaser allele frequency", out$pd, 4L)

    .falconer_print_rule()
    .falconer_print_section("Derived Parameters")
    .falconer_print_parameter("Additive effect (a)", out$a, 4L)
    .falconer_print_parameter("Dominance effect (delta)", out$delta, 4L)
    .falconer_print_parameter("Centering mean (m)", out$m, 4L)
    .falconer_print_parameter(
      "Residual variance", out$residual_variance, 4L
    )
    .falconer_print_parameter("Residual SD", out$residual_sd, 4L)

    .falconer_print_rule()
    .falconer_print_section("Genotype-Specific Quantities")
    .falconer_print_table(list(
      list(label = "Mixing proportion", values = out$pi, digits = 4L),
      list(label = "Trait mean", values = out$mu, digits = 4L)
    ))

    .falconer_print_rule()
    .falconer_print_section("Validation")
    .falconer_print_parameter(
      "Weighted population mean", out$weighted_mean, 6L
    )
    .falconer_print_parameter("Total variance", out$total_variance, 6L)
    .falconer_print_rule()
  }

  invisible(out)
}

#' Power for One-Way ANOVA Under the Falconer Model
#'
#' Calculates power for a one-way ANOVA comparing the means of the three SNP
#' genotype groups for one continuous quantitative trait.
#'
#' @param N Integer total sample size greater than 3.
#' @param alpha Significance level in \eqn{(0,1)}.
#' @inheritParams falconer_parameters
#' @param count_method Genotype count method. \code{"rounded"} rounds expected
#'   genotype counts to the nearest integer before computing the ANOVA
#'   quantities. \code{"expected"} uses expected genotype counts directly
#'   without rounding.
#' @param verbose Logical. If \code{TRUE}, prints a concise summary.
#'
#' @return An object of class \code{"qtl_anova_power"} containing power,
#'   sample size, degrees of freedom, critical value, non-centrality parameter,
#'   genotype counts, genotype means, residual variance, and the complete
#'   Falconer parameter object.
#'
#' @examples
#' qtl_anova_power(
#'   N = 996, alpha = 0.0001, qtl_var = 0.025,
#'   tau = 0.5, pd = 0.15, verbose = FALSE
#' )
#'
#' @references
#' Gordon et al. (2020), \emph{Heterogeneity in Statistical Genetics},
#' Chapter 6, Section 6.1, Equations 6.1, 6.7, and 6.8.
#'
#' @importFrom stats pf qf
#' @export
qtl_anova_power <- function(
    N,
    alpha,
    qtl_var,
    tau,
    pd,
    count_method = c("rounded", "expected"),
    verbose = TRUE
) {
  count_method <- match.arg(count_method)
  .falconer_check_scalar(N, "N", lower = 3, lower_open = TRUE)
  if (N != floor(N)) stop("N must be an integer.", call. = FALSE)
  .falconer_check_scalar(
    alpha, "alpha", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_flag(verbose, "verbose")

  parameters <- falconer_parameters(
    qtl_var = qtl_var, tau = tau, pd = pd, verbose = FALSE
  )
  parameters$alpha <- alpha
  components <- .falconer_anova_components(N, parameters, count_method)
  if (is.null(components)) {
    stop(
      "Rounded genotype counts must all be positive for a three-group ANOVA.",
      call. = FALSE
    )
  }

  out <- list(
    test = "one-way ANOVA for a single quantitative trait by genotype",
    N = N,
    alpha = alpha,
    df1 = components$df1,
    df2 = components$df2,
    critical = components$critical,
    lambda = components$lambda,
    power = components$power,
    count_method = count_method,
    genotype_counts = components$counts,
    mu = parameters$mu,
    residual_variance = parameters$residual_variance,
    falconer = parameters
  )
  class(out) <- "qtl_anova_power"

  if (isTRUE(verbose)) {
    .falconer_print_header(
      "Falconer Quantitative Trait",
      "One-Way ANOVA Power"
    )

    .falconer_print_section("Study Design")
    .falconer_print_parameter("Total sample size", out$N, integer = TRUE)
    .falconer_print_parameter(
      "Significance level (alpha)", out$alpha, 2L, scientific = TRUE
    )
    .falconer_print_parameter(
      "Degrees of freedom (numerator, denominator)",
      sprintf("(%d, %d)", out$df1, out$df2)
    )
    .falconer_print_parameter(
      "Genotype count method",
      if (out$count_method == "rounded") "Rounded" else "Expected"
    )

    .falconer_print_rule()
    .falconer_print_model(out$falconer)

    .falconer_print_rule()
    .falconer_print_section("Genotype-Specific Quantities")
    .falconer_print_table(list(
      list(label = "Mixing proportion", values = out$falconer$pi, digits = 4L),
      list(label = "Trait mean", values = out$mu, digits = 4L),
      list(
        label = "Genotype count",
        values = out$genotype_counts,
        digits = if (out$count_method == "rounded") 0L else 2L,
        integer = out$count_method == "rounded"
      )
    ))

    .falconer_print_rule()
    .falconer_print_section("Test Results")
    .falconer_print_parameter(
      "Non-centrality parameter", out$lambda, 5L
    )
    .falconer_print_parameter("Power", out$power, 6L)
    .falconer_print_rule()
  }

  invisible(out)
}

#' Minimum Sample Size for One-Way ANOVA Under the Falconer Model
#'
#' Finds the smallest sufficient integer total sample size for the three-group
#' one-way ANOVA. Both denominator degrees of freedom and the non-centrality
#' parameter are recalculated for every candidate sample size.
#'
#' @param power Target power in \eqn{(0,1)}.
#' @inheritParams qtl_anova_power
#' @param multiple_of_three Logical. If \code{TRUE}, searches sample sizes in
#'   multiples of three. If \code{FALSE}, searches every integer sample size.
#'
#' @return An object of class \code{"qtl_anova_mssn"} containing the smallest
#'   sufficient \code{N}, achieved power, degrees of freedom, non-centrality
#'   parameter, genotype counts, and Falconer parameters.
#'
#' @examples
#' qtl_anova_mssn(
#'   power = 0.8, alpha = 0.0001, qtl_var = 0.025,
#'   tau = 0.5, pd = 0.15, verbose = FALSE
#' )
#'
#' @importFrom stats pf qf
#' @export
qtl_anova_mssn <- function(
    power,
    alpha,
    qtl_var,
    tau,
    pd,
    count_method = c("rounded", "expected"),
    multiple_of_three = TRUE,
    verbose = TRUE
) {
  count_method <- match.arg(count_method)
  .falconer_check_scalar(
    power, "power", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_scalar(
    alpha, "alpha", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  .falconer_check_flag(multiple_of_three, "multiple_of_three")
  .falconer_check_flag(verbose, "verbose")

  parameters <- falconer_parameters(
    qtl_var = qtl_var, tau = tau, pd = pd, verbose = FALSE
  )
  step <- if (multiple_of_three) 3L else 1L
  first_N <- if (multiple_of_three) 6L else 4L
  max_N <- 10000000L
  chunk_length <- 100000L
  found <- NULL

  for (chunk_start in seq.int(first_N, max_N, by = step * chunk_length)) {
    chunk_end <- min(max_N, chunk_start + step * (chunk_length - 1L))
    candidates <- seq.int(chunk_start, chunk_end, by = step)
    df2 <- candidates - 3

    if (identical(count_method, "rounded")) {
      counts <- round(outer(candidates, parameters$pi))
      valid <- rowSums(counts <= 0) == 0
      lambda <- as.vector(counts %*% (parameters$mu^2)) /
        parameters$residual_variance
    } else {
      valid <- rep(TRUE, length(candidates))
      lambda <- candidates * sum(parameters$pi * parameters$mu^2) /
        parameters$residual_variance
    }

    critical <- stats::qf(1 - alpha, df1 = 2, df2 = df2)
    candidate_power <- stats::pf(
      critical, df1 = 2, df2 = df2, ncp = lambda, lower.tail = FALSE
    )
    sufficient <- which(valid & candidate_power >= power)
    if (length(sufficient)) {
      index <- sufficient[1L]
      found <- list(
        N = candidates[index],
        power = candidate_power[index],
        lambda = lambda[index],
        critical = critical[index]
      )
      break
    }
  }

  if (is.null(found)) {
    stop(
      "No sufficient ANOVA sample size was found at or below 10,000,000.",
      call. = FALSE
    )
  }

  counts <- .falconer_anova_counts(found$N, parameters$pi, count_method)
  out <- list(
    test = "one-way ANOVA for a single quantitative trait by genotype",
    target_power = power,
    achieved_power = found$power,
    N = found$N,
    alpha = alpha,
    df1 = 2,
    df2 = found$N - 3,
    critical = found$critical,
    lambda = found$lambda,
    count_method = count_method,
    multiple_of_three = multiple_of_three,
    genotype_counts = counts,
    mu = parameters$mu,
    residual_variance = parameters$residual_variance,
    falconer = parameters
  )
  class(out) <- "qtl_anova_mssn"

  if (isTRUE(verbose)) {
    .falconer_print_header(
      "Falconer Quantitative Trait",
      "One-Way ANOVA Minimum Sample Size"
    )

    .falconer_print_section("Target Design")
    .falconer_print_parameter("Target power", out$target_power, 3L)
    .falconer_print_parameter(
      "Significance level (alpha)", out$alpha, 2L, scientific = TRUE
    )
    .falconer_print_parameter(
      "Genotype count method",
      if (out$count_method == "rounded") "Rounded" else "Expected"
    )

    .falconer_print_rule()
    .falconer_print_section("Required Sample Size")
    .falconer_print_parameter(
      "Minimum total sample size", out$N, integer = TRUE
    )
    .falconer_print_parameter("Achieved power", out$achieved_power, 6L)
    .falconer_print_parameter(
      "Degrees of freedom (numerator, denominator)",
      sprintf("(%d, %d)", out$df1, out$df2)
    )
    .falconer_print_parameter("Required NCP", out$lambda, 5L)

    .falconer_print_rule()
    .falconer_print_model(out$falconer)

    .falconer_print_rule()
    .falconer_print_section("Genotype-Specific Quantities")
    .falconer_print_table(list(
      list(
        label = "Genotype count",
        values = out$genotype_counts,
        digits = if (out$count_method == "rounded") 0L else 2L,
        integer = out$count_method == "rounded"
      ),
      list(label = "Trait mean", values = out$mu, digits = 4L),
      list(label = "Mixing proportion", values = out$falconer$pi, digits = 4L)
    ))
    .falconer_print_rule()
  }

  invisible(out)
}
