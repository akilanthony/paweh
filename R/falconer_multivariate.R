#' Internal validation for the multivariate Falconer model.
#'
#' @noRd
.falconer_mv_validate_model <- function(qtl_var, tau, pd, cor_matrix) {
  if (!is.numeric(qtl_var) || length(qtl_var) < 1L ||
      any(!is.finite(qtl_var)) || any(qtl_var <= 0 | qtl_var >= 1)) {
    stop("qtl_var must be a finite numeric vector with entries in (0, 1).", call. = FALSE)
  }
  p <- length(qtl_var)
  if (!is.numeric(tau) || length(tau) != p || any(!is.finite(tau))) {
    stop("tau must be a finite numeric vector with length(qtl_var) entries.", call. = FALSE)
  }
  .falconer_check_scalar(
    pd, "pd", lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE
  )
  if (!is.matrix(cor_matrix) || !is.numeric(cor_matrix) ||
      !identical(dim(cor_matrix), c(p, p)) || any(!is.finite(cor_matrix))) {
    stop("cor_matrix must be a finite numeric p x p matrix.", call. = FALSE)
  }
  if (!isTRUE(all.equal(cor_matrix, t(cor_matrix), tolerance = 1e-12))) {
    stop("cor_matrix must be symmetric.", call. = FALSE)
  }
  if (any(abs(diag(cor_matrix) - 1) > 1e-12)) {
    stop("cor_matrix must have a unit diagonal.", call. = FALSE)
  }
  if (any(cor_matrix < -1 - 1e-12 | cor_matrix > 1 + 1e-12)) {
    stop("All entries of cor_matrix must be valid correlations in [-1, 1].", call. = FALSE)
  }
  if (inherits(try(chol(cor_matrix), silent = TRUE), "try-error")) {
    stop("cor_matrix must be positive definite.", call. = FALSE)
  }
  invisible(p)
}

#' Construct the common multivariate Falconer model.
#'
#' The returned mean matrix has phenotypes in rows and genotypes 0, 1, and 2
#' in columns.
#'
#' @noRd
.falconer_mv_parameters <- function(qtl_var, tau, pd, cor_matrix) {
  p <- .falconer_mv_validate_model(qtl_var, tau, pd, cor_matrix)
  one_trait <- lapply(seq_len(p), function(i) {
    qtl_falconer_parameters(qtl_var[i], tau[i], pd, verbose = FALSE)
  })
  additive_variance <- vapply(one_trait, function(x) {
    2 * x$pd * x$p_plus * (x$a + x$delta * (x$p_plus - x$pd))^2
  }, numeric(1))
  dominance_variance <- vapply(one_trait, function(x) {
    (2 * x$pd * x$p_plus * x$delta)^2
  }, numeric(1))
  if (any(abs(additive_variance + dominance_variance - qtl_var) > 1e-10)) {
    stop("Additive and dominance variances do not reproduce qtl_var.", call. = FALSE)
  }

  phenotype_names <- paste0("phenotype_", seq_len(p))
  genotype_names <- paste0("genotype_", 0:2)
  mean_matrix <- do.call(rbind, lapply(one_trait, `[[`, "mu"))
  dimnames(mean_matrix) <- list(phenotype_names, genotype_names)
  residual_variance <- 1 - qtl_var
  residual_sd <- sqrt(residual_variance)
  Sigma <- diag(residual_sd, p) %*% cor_matrix %*% diag(residual_sd, p)
  dimnames(Sigma) <- dimnames(cor_matrix) <- list(phenotype_names, phenotype_names)
  if (inherits(try(chol(Sigma), silent = TRUE), "try-error")) {
    stop("The resulting residual covariance matrix must be positive definite.", call. = FALSE)
  }

  parameters <- data.frame(
    qtl_var = qtl_var,
    tau = tau,
    a = vapply(one_trait, `[[`, numeric(1), "a"),
    delta = vapply(one_trait, `[[`, numeric(1), "delta"),
    m = vapply(one_trait, `[[`, numeric(1), "m"),
    additive_variance = additive_variance,
    dominance_variance = dominance_variance,
    residual_variance = residual_variance,
    residual_sd = residual_sd,
    row.names = phenotype_names
  )
  pi <- one_trait[[1L]]$pi

  list(
    number_of_phenotypes = p,
    number_of_genotype_groups = 3L,
    pd = pd,
    p_plus = 1 - pd,
    genotype_frequencies = pi,
    parameters = parameters,
    mean_matrix = mean_matrix,
    mean_matrix_orientation = "rows = phenotypes; columns = genotypes",
    phenotype_correlation_matrix = cor_matrix,
    residual_covariance_matrix = Sigma,
    single_trait_parameters = one_trait
  )
}

#' Gordon et al. analytic Pillai quantities.
#'
#' Implements the general matrix form preceding Equation (1): Phi-star is
#' Sigma^-1 (CB)' times the inverse of C diag(1 / pi) C' times (CB), while V-star is the sum of
#' phi-star / (1 + phi-star), and lambda = N s V-star / (s - V-star).
#'
#' @noRd
.falconer_mv_pillai_components <- function(N, model) {
  p <- model$number_of_phenotypes
  g <- model$number_of_genotype_groups
  s <- min(g - 1L, p)
  C <- cbind(rep(1, g - 1L), -diag(g - 1L))
  B <- t(model$mean_matrix)
  CB <- C %*% B
  middle <- solve(C %*% diag(1 / model$genotype_frequencies) %*% t(C))
  hypothesis_per_subject <- t(CB) %*% middle %*% CB
  phi_star <- solve(model$residual_covariance_matrix, hypothesis_per_subject)
  roots <- eigen(phi_star, only.values = TRUE)$values
  if (any(abs(Im(roots)) > 1e-9)) {
    stop("Pillai characteristic roots were not real-valued.", call. = FALSE)
  }
  roots <- sort(Re(roots), decreasing = TRUE)[seq_len(s)]
  roots[roots < 0 & roots > -1e-10] <- 0
  if (any(!is.finite(roots)) || any(roots < 0)) {
    stop("Pillai characteristic roots must be finite and non-negative.", call. = FALSE)
  }
  V_star <- sum(roots / (1 + roots))
  lambda_per_subject <- s * V_star / (s - V_star)
  df1 <- (g - 1L) * p
  df2 <- s * (N - g + s - p)
  if (!is.finite(df2) || df2 <= 0) {
    stop("N is too small for positive Pillai denominator degrees of freedom.", call. = FALSE)
  }
  list(
    contrast_matrix = C,
    coefficient_matrix = B,
    contrast_means = CB,
    hypothesis_per_subject = hypothesis_per_subject,
    phi_star = phi_star,
    eigenvalues = roots,
    s = s,
    V_star = V_star,
    lambda_per_subject = lambda_per_subject,
    lambda = N * lambda_per_subject,
    numerator_df = df1,
    denominator_df = df2
  )
}

#' Compute Pillai power from a model and total sample size.
#'
#' @noRd
.falconer_mv_pillai_power <- function(N, alpha, model) {
  comp <- .falconer_mv_pillai_components(N, model)
  critical <- stats::qf(1 - alpha, comp$numerator_df, comp$denominator_df)
  power <- stats::pf(
    critical, comp$numerator_df, comp$denominator_df,
    ncp = comp$lambda, lower.tail = FALSE
  )
  c(comp, list(critical_value = critical, power = power))
}

#' Evaluate one multivariate-normal rectangle probability.
#'
#' Miwa's deterministic algorithm is used through dimension 20. Higher
#' dimensions use Genz-Bretz with a locally fixed seed and restored RNG state.
#'
#' @noRd
.falconer_mv_rectangle_probability <- function(lower, upper, mean, Sigma) {
  p <- length(mean)
  if (p == 1L) {
    probability <- stats::pnorm(upper, mean, sqrt(Sigma[1, 1])) -
      stats::pnorm(lower, mean, sqrt(Sigma[1, 1]))
    return(list(
      probability = probability, error = 0, message = "Normal completion",
      algorithm = "univariate pnorm"
    ))
  }
  algorithm_name <- if (p <= 20L) "Miwa" else "GenzBretz"
  algorithm <- if (p <= 20L) {
    mvtnorm::Miwa(steps = 128L)
  } else {
    mvtnorm::GenzBretz(maxpts = 25000L * p, abseps = 1e-8, releps = 0)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  if (p > 20L) set.seed(271828L)
  on.exit({
    if (p > 20L) {
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  value <- mvtnorm::pmvnorm(
    lower = lower, upper = upper, mean = mean, sigma = Sigma,
    algorithm = algorithm
  )
  list(
    probability = as.numeric(value),
    error = if (is.null(attr(value, "error"))) NA_real_ else attr(value, "error"),
    message = if (is.null(attr(value, "msg"))) "" else attr(value, "msg"),
    algorithm = algorithm_name
  )
}

#' Derive joint-AND threshold penetrances and conditional frequencies.
#'
#' @noRd
.falconer_mv_threshold_components <- function(model, x_upper, x_lower) {
  p <- model$number_of_phenotypes
  check_percent <- function(x, name) {
    if (!is.numeric(x) || length(x) != p || any(!is.finite(x)) ||
        any(x <= 0 | x >= 100)) {
      stop(sprintf("%s must be a finite numeric vector of length p with entries in (0, 100).", name), call. = FALSE)
    }
  }
  check_percent(x_upper, "x_upper")
  check_percent(x_lower, "x_lower")
  upper_threshold <- stats::qnorm(1 - x_upper / 100)
  lower_threshold <- stats::qnorm(x_lower / 100)
  if (any(lower_threshold >= upper_threshold)) {
    stop("Every lower threshold must be less than its corresponding upper threshold.", call. = FALSE)
  }

  affected_details <- unaffected_details <- vector("list", 3L)
  for (j in 1:3) {
    affected_details[[j]] <- .falconer_mv_rectangle_probability(
      lower = upper_threshold, upper = rep(Inf, p),
      mean = model$mean_matrix[, j], Sigma = model$residual_covariance_matrix
    )
    unaffected_details[[j]] <- .falconer_mv_rectangle_probability(
      lower = rep(-Inf, p), upper = lower_threshold,
      mean = model$mean_matrix[, j], Sigma = model$residual_covariance_matrix
    )
  }
  f_affected <- vapply(affected_details, `[[`, numeric(1), "probability")
  f_unaffected <- vapply(unaffected_details, `[[`, numeric(1), "probability")
  names(f_affected) <- names(f_unaffected) <- names(model$genotype_frequencies)
  prev_affected <- sum(model$genotype_frequencies * f_affected)
  prev_unaffected <- sum(model$genotype_frequencies * f_unaffected)
  if (any(!is.finite(c(prev_affected, prev_unaffected))) ||
      prev_affected <= 0 || prev_unaffected <= 0) {
    stop("Joint threshold selection probabilities must be positive and finite.", call. = FALSE)
  }
  g_case <- model$genotype_frequencies * f_affected / prev_affected
  g_control <- model$genotype_frequencies * f_unaffected / prev_unaffected
  if (abs(sum(g_case) - 1) > 1e-10 || abs(sum(g_control) - 1) > 1e-10) {
    stop("Conditional genotype frequencies do not sum to one.", call. = FALSE)
  }
  list(
    x_upper = x_upper, x_lower = x_lower,
    upper_threshold = upper_threshold, lower_threshold = lower_threshold,
    penetrances = list(affected = f_affected, unaffected = f_unaffected),
    prevalences = c(affected = prev_affected, unaffected = prev_unaffected),
    frequencies = list(case = g_case, control = g_control),
    integration = list(affected = affected_details, unaffected = unaffected_details),
    selection_definition = "affected: all traits >= upper thresholds; unaffected: all traits <= lower thresholds"
  )
}

#' Expected-count and sparse-cell diagnostics.
#'
#' @noRd
.falconer_mv_expected_counts <- function(N_case, N_control, frequencies) {
  counts <- rbind(
    cases = N_case * frequencies$case,
    controls = N_control * frequencies$control
  )
  sparse <- sum(counts < 1)
  list(
    expected_counts = counts,
    cells_below_one = sparse,
    percent_cells_below_one = 100 * sparse / length(counts)
  )
}

#' Print a compact labelled matrix.
#'
#' @noRd
.falconer_mv_print_matrix <- function(x, digits = 5L) {
  shown <- formatC(x, format = "f", digits = digits)
  dim(shown) <- dim(x)
  dimnames(shown) <- dimnames(x)
  message(paste(utils::capture.output(print(shown, quote = FALSE)), collapse = "\n"))
}

#' Print the common multivariate model.
#'
#' @noRd
.falconer_mv_print_model <- function(model) {
  .falconer_print_section("Falconer Model")
  .falconer_print_parameter("Increaser allele frequency", model$pd, 4L)
  .falconer_mv_print_matrix(t(model$parameters), 6L)
  .falconer_print_rule()
  .falconer_print_section("Genotype-Specific Mean Vectors")
  .falconer_mv_print_matrix(model$mean_matrix, 6L)
  .falconer_print_rule()
  .falconer_print_section("Phenotype Correlation Matrix")
  .falconer_mv_print_matrix(model$phenotype_correlation_matrix, 5L)
  .falconer_print_section("Residual Covariance Matrix")
  .falconer_mv_print_matrix(model$residual_covariance_matrix, 6L)
}

#' Multivariate Falconer Power
#'
#' Computes prospective power for a one-way MANOVA using Pillai's trace or for
#' a genotype chi-square test after joint multivariate threshold selection.
#'
#' @param N Integer total sample size for `test = "pillai"`.
#' @param N_case Selected case sample size for
#'   `test = "threshold_chisq"`.
#' @param alpha Significance level in (0, 1).
#' @param qtl_var Numeric vector of phenotype-specific QTL variances in
#'   (0, 1).
#' @param tau Numeric vector of phenotype-specific dominance/additivity ratios.
#' @param pd Shared increaser-allele frequency in (0, 1).
#' @param cor_matrix Positive-definite phenotype correlation matrix. Its order
#'   must match `qtl_var` and `tau`.
#' @param test Either `"pillai"` or `"threshold_chisq"`.
#' @param x_upper,x_lower Vectors of upper- and lower-tail percentages. A case
#'   satisfies every upper threshold (joint AND); a control satisfies every
#'   lower threshold (joint AND).
#' @param k Control/case ratio for the threshold chi-square design.
#' @param verbose Logical; whether to print a polished summary.
#'
#' @details
#' Genotype-specific means are stored in a matrix whose rows are phenotypes and
#' columns are genotypes 0, 1, and 2. All genotypes share the residual
#' covariance matrix formed by scaling `cor_matrix` by residual standard
#' deviations. Joint probabilities use `mvtnorm::pmvnorm()` with the
#' deterministic Miwa algorithm through 20 dimensions; higher dimensions use
#' reproducibly seeded Genz--Bretz integration and return its error/status.
#'
#' Pillai power follows Gordon et al.'s general matrix derivation: the
#' characteristic roots of Phi-star determine V-star and
#' lambda = N s V-star / (s - V-star). The null critical value and
#' alternative power use central and noncentral F distributions.
#'
#' @return An object of class `"qtl_multivariate_power_full"`. Both test
#'   modes retain the complete Falconer model and auditable intermediates.
#'
#' @examples
#' cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
#' qtl_multivariate_power_full(
#'   N = 4514, alpha = 5e-8, qtl_var = c(0.01, 0.005),
#'   tau = c(0, 0.5), pd = 0.25, cor_matrix = cor_matrix,
#'   test = "pillai", verbose = FALSE
#' )
#'
#' @references
#' Gordon, Finch, and Kim (2020), \emph{Heterogeneity in Statistical
#' Genetics}, Chapter 6, Sections 6.2.1--6.2.5, including Equations 6.9--6.10.
#'
#' Gordon D, Londono D, Patel P, Kim W, Finch SJ, Heiman GA (2017).
#' "An Analytic Solution to the Computation of Power and Sample Size for
#' Genetic Association Studies under a Pleiotropic Mode of Inheritance."
#' \emph{Human Heredity}, 81(4), 194--209. doi:10.1159/000457135.
#'
#' @importFrom stats pf qf qnorm pnorm
#' @importFrom mvtnorm pmvnorm Miwa GenzBretz
#' @export
qtl_multivariate_power_full <- function(
    N = NULL,
    N_case = NULL,
    alpha,
    qtl_var,
    tau,
    pd,
    cor_matrix,
    test = c("pillai", "threshold_chisq"),
    x_upper = NULL,
    x_lower = NULL,
    k = 1,
    verbose = TRUE
) {
  test <- match.arg(test)
  .falconer_check_scalar(alpha, "alpha", 0, 1, TRUE, TRUE)
  .falconer_check_flag(verbose, "verbose")
  model <- .falconer_mv_parameters(qtl_var, tau, pd, cor_matrix)

  if (identical(test, "pillai")) {
    .falconer_check_scalar(N, "N", lower = 0, lower_open = TRUE)
    if (N != floor(N)) stop("N must be an integer total sample size.", call. = FALSE)
    pillai <- .falconer_mv_pillai_power(N, alpha, model)
    out <- list(
      test = "pillai", N = N, alpha = alpha, power = pillai$power,
      number_of_phenotypes = model$number_of_phenotypes,
      number_of_genotype_groups = model$number_of_genotype_groups,
      genotype_frequencies = model$genotype_frequencies,
      expected_genotype_counts = N * model$genotype_frequencies,
      noncentrality_parameter = pillai$lambda,
      numerator_df = pillai$numerator_df,
      denominator_df = pillai$denominator_df,
      critical_value = pillai$critical_value,
      pillai = pillai,
      falconer = model
    )
  } else {
    .falconer_check_scalar(N_case, "N_case", lower = 0, lower_open = TRUE)
    .falconer_check_scalar(k, "k", lower = 0, lower_open = TRUE)
    threshold <- .falconer_mv_threshold_components(model, x_upper, x_lower)
    cc <- cc_chisq_power_locus_heterogeneity(
      N_case = N_case, alpha = alpha,
      g_case_assoc = threshold$frequencies$case,
      g_ctrl = threshold$frequencies$control,
      pi = 1, k = k, verbose = FALSE
    )
    diagnostics <- .falconer_mv_expected_counts(
      cc$N_case, cc$N_ctrl, threshold$frequencies
    )
    out <- c(list(
      test = "threshold_chisq", alpha = alpha,
      N_case = cc$N_case, N_control = cc$N_ctrl, N_total = cc$N_total,
      k = k, df = cc$df, noncentrality_parameter = cc$lambda,
      S = cc$S, power = cc$power,
      thresholds = threshold, falconer = model
    ), diagnostics)
  }
  class(out) <- "qtl_multivariate_power_full"

  if (isTRUE(verbose)) {
    subtitle <- if (test == "pillai") "Pillai MANOVA Power" else
      "Threshold-Selected Genotype Chi-Square Power"
    .falconer_print_header("Falconer Multivariate Quantitative Traits", subtitle)
    .falconer_print_section("Study Design")
    .falconer_print_parameter("Number of phenotypes", model$number_of_phenotypes, integer = TRUE)
    if (test == "pillai") .falconer_print_parameter("Total sample size", out$N, integer = TRUE)
    else {
      .falconer_print_parameter("Selected cases", out$N_case, 2L)
      .falconer_print_parameter("Selected controls", out$N_control, 2L)
      .falconer_print_parameter("Control/case ratio", out$k, 3L)
    }
    .falconer_print_parameter("Significance level (alpha)", alpha, 2L, scientific = TRUE)
    .falconer_print_rule()
    .falconer_mv_print_model(model)
    .falconer_print_rule()
    if (test == "pillai") {
      .falconer_print_section("Pillai Test Results")
      .falconer_print_parameter("Numerator df", out$numerator_df, integer = TRUE)
      .falconer_print_parameter("Denominator df", out$denominator_df, 2L)
      .falconer_print_parameter("Non-centrality parameter", out$noncentrality_parameter, 6L)
      .falconer_print_parameter("Critical value", out$critical_value, 6L)
      .falconer_print_parameter("Power", out$power, 6L)
    } else {
      .falconer_print_section("Threshold Selection")
      .falconer_mv_print_matrix(rbind(
        upper_tail_percent = out$thresholds$x_upper,
        upper_threshold = out$thresholds$upper_threshold,
        lower_tail_percent = out$thresholds$x_lower,
        lower_threshold = out$thresholds$lower_threshold
      ), 5L)
      .falconer_print_section("Joint Penetrances")
      .falconer_mv_print_matrix(rbind(
        affected = out$thresholds$penetrances$affected,
        unaffected = out$thresholds$penetrances$unaffected
      ), 6L)
      .falconer_print_parameter("Affected prevalence", out$thresholds$prevalences["affected"], 6L)
      .falconer_print_parameter("Unaffected prevalence", out$thresholds$prevalences["unaffected"], 6L)
      .falconer_print_section("Conditional Genotype Frequencies")
      .falconer_mv_print_matrix(rbind(
        cases = out$thresholds$frequencies$case,
        controls = out$thresholds$frequencies$control
      ), 6L)
      .falconer_print_section("Test Results")
      .falconer_print_parameter("Non-centrality parameter", out$noncentrality_parameter, 6L)
      .falconer_print_parameter("Power", out$power, 6L)
      .falconer_print_section("Expected Genotype Counts")
      .falconer_mv_print_matrix(out$expected_counts, 3L)
      .falconer_print_parameter("Cells with expected count < 1", out$cells_below_one, integer = TRUE)
      .falconer_print_parameter("Percent of cells with expected count < 1", out$percent_cells_below_one, 2L)
    }
    .falconer_print_rule()
  }
  invisible(out)
}

#' Multivariate Falconer Minimum Sample Size
#'
#' Finds the minimum integer total sample size for Pillai MANOVA or the minimum
#' selected case-control sample for a joint-threshold genotype chi-square test.
#'
#' @param power Target power in (0, 1).
#' @inheritParams qtl_multivariate_power_full
#'
#' @details
#' Pillai MSSN is found by an integer search because its denominator degrees of
#' freedom depend on sample size. For threshold selection, the statistical MSSN
#' is kept separate from the expected population screening burden.
#'
#' @return An object of class `"qtl_multivariate_mssn_full"` containing
#'   integer MSSN, achieved power, complete model quantities, and test-specific
#'   audit information.
#'
#' @examples
#' cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
#' qtl_multivariate_mssn_full(
#'   power = 0.95, alpha = 5e-8, qtl_var = c(0.01, 0.005),
#'   tau = c(0, 0.5), pd = 0.25, cor_matrix = cor_matrix,
#'   test = "pillai", verbose = FALSE
#' )
#'
#' @inherit qtl_multivariate_power_full references
#' @importFrom stats uniroot
#' @export
qtl_multivariate_mssn_full <- function(
    power,
    alpha,
    qtl_var,
    tau,
    pd,
    cor_matrix,
    test = c("pillai", "threshold_chisq"),
    x_upper = NULL,
    x_lower = NULL,
    k = 1,
    verbose = TRUE
) {
  test <- match.arg(test)
  .falconer_check_scalar(power, "power", 0, 1, TRUE, TRUE)
  .falconer_check_scalar(alpha, "alpha", 0, 1, TRUE, TRUE)
  .falconer_check_flag(verbose, "verbose")
  model <- .falconer_mv_parameters(qtl_var, tau, pd, cor_matrix)

  if (test == "pillai") {
    p <- model$number_of_phenotypes
    s <- min(2L, p)
    lower <- max(4L, 3L - s + p + 1L)
    power_at <- function(N) .falconer_mv_pillai_power(N, alpha, model)$power
    upper <- lower
    while (power_at(upper) < power) {
      upper <- upper * 2L
      if (upper > .Machine$integer.max / 2L) {
        stop("Unable to bracket the requested Pillai sample size.", call. = FALSE)
      }
    }
    lo <- lower
    hi <- upper
    while (lo < hi) {
      mid <- floor((lo + hi) / 2)
      if (power_at(mid) >= power) hi <- mid else lo <- mid + 1L
    }
    N <- as.integer(lo)
    achieved <- .falconer_mv_pillai_power(N, alpha, model)
    continuous_root <- stats::uniroot(
      function(n) power_at(n) - power,
      lower = max(lower, N - 1), upper = N, tol = 1e-10
    )$root
    out <- list(
      test = "pillai", target_power = power, achieved_power = achieved$power,
      alpha = alpha, N = N, historical_fractional_mssn = continuous_root,
      number_of_phenotypes = model$number_of_phenotypes,
      number_of_genotype_groups = model$number_of_genotype_groups,
      genotype_frequencies = model$genotype_frequencies,
      expected_genotype_counts = N * model$genotype_frequencies,
      noncentrality_parameter = achieved$lambda,
      numerator_df = achieved$numerator_df,
      denominator_df = achieved$denominator_df,
      critical_value = achieved$critical_value,
      pillai = achieved, falconer = model
    )
  } else {
    .falconer_check_scalar(k, "k", lower = 0, lower_open = TRUE)
    threshold <- .falconer_mv_threshold_components(model, x_upper, x_lower)
    cc <- cc_chisq_mssn_locus_heterogeneity(
      power = power, alpha = alpha,
      g_case_assoc = threshold$frequencies$case,
      g_ctrl = threshold$frequencies$control,
      pi = 1, k = k, verbose = FALSE
    )
    achieved <- cc_chisq_power_locus_heterogeneity(
      N_case = cc$N_case, alpha = alpha,
      g_case_assoc = threshold$frequencies$case,
      g_ctrl = threshold$frequencies$control,
      pi = 1, k = k, verbose = FALSE
    )
    diagnostics <- .falconer_mv_expected_counts(
      cc$N_case, cc$N_ctrl, threshold$frequencies
    )
    screened_cases <- cc$N_case / threshold$prevalences["affected"]
    screened_controls <- cc$N_ctrl / threshold$prevalences["unaffected"]
    out <- c(list(
      test = "threshold_chisq", target_power = power,
      achieved_power = achieved$power, alpha = alpha,
      N_case = cc$N_case, N_control = cc$N_ctrl, N_total = cc$N_total,
      historical_fractional_cases = cc$lambda_star / (k * cc$S),
      k = k, df = cc$df, target_noncentrality_parameter = cc$lambda_star,
      noncentrality_parameter = achieved$lambda, S = cc$S,
      expected_population_screened_cases = as.numeric(screened_cases),
      expected_population_screened_controls = as.numeric(screened_controls),
      screening = list(cases = as.numeric(screened_cases), controls = as.numeric(screened_controls)),
      thresholds = threshold, falconer = model
    ), diagnostics)
  }
  class(out) <- "qtl_multivariate_mssn_full"

  if (isTRUE(verbose)) {
    subtitle <- if (test == "pillai") "Pillai MANOVA Minimum Sample Size" else
      "Threshold-Selected Genotype Chi-Square Minimum Sample Size"
    .falconer_print_header("Falconer Multivariate Quantitative Traits", subtitle)
    .falconer_print_section("Target Design")
    .falconer_print_parameter("Target power", power, 3L)
    .falconer_print_parameter("Significance level (alpha)", alpha, 2L, scientific = TRUE)
    if (test == "threshold_chisq") .falconer_print_parameter("Control/case ratio", k, 3L)
    .falconer_print_rule()
    .falconer_mv_print_model(model)
    .falconer_print_rule()
    if (test == "pillai") {
      .falconer_print_section("Required Sample Size")
      .falconer_print_parameter("Minimum total sample size", out$N, integer = TRUE)
      .falconer_print_parameter("Achieved power", out$achieved_power, 6L)
      .falconer_print_parameter("Non-centrality parameter", out$noncentrality_parameter, 6L)
      .falconer_print_parameter("Numerator df", out$numerator_df, integer = TRUE)
      .falconer_print_parameter("Denominator df", out$denominator_df, 2L)
      .falconer_print_parameter("Critical value", out$critical_value, 6L)
    } else {
      .falconer_print_section("Threshold Selection")
      .falconer_mv_print_matrix(rbind(
        upper_tail_percent = out$thresholds$x_upper,
        upper_threshold = out$thresholds$upper_threshold,
        lower_tail_percent = out$thresholds$x_lower,
        lower_threshold = out$thresholds$lower_threshold
      ), 5L)
      .falconer_print_section("Joint Penetrances")
      .falconer_mv_print_matrix(rbind(
        affected = out$thresholds$penetrances$affected,
        unaffected = out$thresholds$penetrances$unaffected
      ), 6L)
      .falconer_print_parameter("Affected prevalence", out$thresholds$prevalences["affected"], 6L)
      .falconer_print_parameter("Unaffected prevalence", out$thresholds$prevalences["unaffected"], 6L)
      .falconer_print_section("Conditional Genotype Frequencies")
      .falconer_mv_print_matrix(rbind(
        cases = out$thresholds$frequencies$case,
        controls = out$thresholds$frequencies$control
      ), 6L)
      .falconer_print_section("Required Selected Sample")
      .falconer_print_parameter("Cases", out$N_case, integer = TRUE)
      .falconer_print_parameter("Controls", out$N_control, integer = TRUE)
      .falconer_print_parameter("Total selected sample size", out$N_total, integer = TRUE)
      .falconer_print_parameter("Non-centrality parameter", out$noncentrality_parameter, 6L)
      .falconer_print_parameter("Achieved power", out$achieved_power, 6L)
      .falconer_print_section("Expected Genotype Counts")
      .falconer_mv_print_matrix(out$expected_counts, 3L)
      .falconer_print_parameter("Cells with expected count < 1", out$cells_below_one, integer = TRUE)
      .falconer_print_parameter("Percent of cells with expected count < 1", out$percent_cells_below_one, 2L)
      .falconer_print_section("Expected Screening Burden")
      .falconer_print_parameter("Population screened for selected cases", out$expected_population_screened_cases, 1L)
      .falconer_print_parameter("Population screened for selected controls", out$expected_population_screened_controls, 1L)
    }
    .falconer_print_rule()
  }
  invisible(out)
}
