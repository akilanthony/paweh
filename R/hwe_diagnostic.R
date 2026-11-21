.hwe_validate_order <- function(labels, label) {
  if (!is.null(labels) && !identical(unname(labels), c("0", "1", "2"))) {
    stop(label, ' must be ordered as "0", "1", "2".', call. = FALSE)
  }
}

.hwe_validate_vector <- function(x, label) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != 3L ||
      any(!is.finite(x)) || any(x < 0)) {
    stop(label, " must be a vector of three finite, nonnegative numbers.",
         call. = FALSE)
  }
  .hwe_validate_order(names(x), paste(label, "names"))
}

.hwe_validate_scalar <- function(x, label, lower, upper, open = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) != 1L ||
      !is.finite(x) || x < lower || x > upper ||
      (open && (x == lower || x == upper))) {
    interval <- if (open) paste0("(", lower, ", ", upper, ")") else
      paste0("[", lower, ", ", upper, "]")
    stop(label, " must be a single finite number in ", interval, ".", call. = FALSE)
  }
}

.hwe_proportions <- function(p, q) {
  c("0" = q^2, "1" = 2 * p * q, "2" = p^2)
}

.hwe_fit <- function(g) {
  het_half <- g[[2L]] / 2
  p <- het_half + g[[3L]]
  q <- het_half + g[[1L]]
  expected <- .hwe_proportions(p, q)
  F <- NA_real_
  if (p > 0 && q > 0) {
    # Either homozygote residual equals half the heterozygote deficit.
    # Use the rarer homozygote to preserve the contrast near a boundary.
    rare_hom <- if (p <= q) 3L else 1L
    residual <- g[[rare_hom]] - expected[[rare_hom]]
    rounding <- 8 * .Machine$double.eps * max(g[[rare_hom]], expected[[rare_hom]])
    F <- if (abs(residual) <= rounding) 0 else residual / (p * q)
  }
  list(allele_frequency = p, expected_probs = expected, F = F)
}

.hwe_exact_p <- function(counts, N) {
  minor <- counts[[2L]] + 2 * min(counts[[1L]], counts[[3L]])
  n_states <- floor(minor / 2) + 1
  if (n_states > 1000000) {
    stop("Exact enumeration exceeds one million conditional tables; use method = 'chisq'.",
         call. = FALSE)
  }

  # The weight at h+2 divided by the weight at h, from the conditional
  # multinomial weights 2^h / (n0! h! n2!). Work outward from the mode in
  # log space so that even very small table probabilities can be compared.
  h_min <- minor %% 2
  h <- h_min + 2 * (seq_len(n_states - 1) - 1)
  log_ratios <- log((minor - h) * (2 * N - minor - h) / ((h + 1) * (h + 2)))
  left <- seq_along(log_ratios) <= sum(log_ratios > 0)
  log_weights <- c(rev(cumsum(-rev(log_ratios[left]))), 0, cumsum(log_ratios[!left]))
  log_observed <- log_weights[(counts[[2L]] - h_min) / 2 + 1]
  tolerance <- 1e-12 + 16 * .Machine$double.eps * abs(log_observed)
  weights <- exp(log_weights)
  sum(weights[log_weights <= log_observed + tolerance]) / sum(weights)
}

#' Test observed genotype counts for Hardy-Weinberg equilibrium
#'
#' Test one biallelic autosomal locus using counts of diploid individuals.
#' The null is Hardy-Weinberg proportions with an unknown allele frequency.
#' Pearson inference estimates that frequency from the counts. Exact inference
#' conditions on the allele total. A large p-value does not establish equilibrium.
#'
#' @param counts Numeric vector of three nonnegative whole counts, in 0, 1,
#'   and 2 copies of the same counted allele. If named, use `"0"`, `"1"`,
#'   and `"2"` in that order. The sum must be positive and no larger than
#'   `2^52`.
#' @param method One of `"exact"` (default) or `"chisq"` (uncorrected Pearson).
#'
#' @details The uncorrected Pearson statistic has an asymptotic chi-square
#'   reference distribution with one degree of freedom: three cells minus one
#'   total-count constraint and one estimated allele frequency. Expected counts
#'   below five prompt a warning for Pearson inference. The two-sided exact
#'   p-value conditions on the observed allele count and sums the probabilities
#'   of tables at least as unlikely as the observed table, including ties
#'   (Wigginton et al., 2005). Exact enumeration is limited to one million
#'   conditional tables. At a monomorphic locus there is one conditional
#'   table (exact p-value 1) and no regular Pearson test. The function warns
#'   and returns `NA` for the Pearson statistic, degrees of freedom, and `F`.
#'
#' @return A list with `N`, `counts`, `observed_probs`, `allele_frequency`,
#'   `expected_probs`, `expected_counts`, Pearson `statistic` and `df`,
#'   the selected `method` and its `p_value`, the signed heterozygote
#'   deficit `F = 1 - observed/expected`, and `sparse_expected` indicating
#'   any expected count below five. Negative `F` indicates heterozygote excess.
#'   The Pearson degrees of freedom describe `statistic`, even when the
#'   selected p-value is exact.
#'
#' @references
#' Wigginton JE, Cutler DJ, Abecasis GR (2005). A note on exact tests of
#' Hardy-Weinberg equilibrium. *American Journal of Human Genetics*
#' 76(5), 887-893. \doi{10.1086/429864}.
#'
#' @seealso [hwe_power()], [hwe_mssn()]
#' @examples
#' hwe_test(c(88, 10, 2))
#' hwe_test(c(40, 30, 30), method = "chisq")
#' @export
hwe_test <- function(counts, method = c("exact", "chisq")) {
  method <- match.arg(method)
  .hwe_validate_vector(counts, "counts")
  N <- sum(counts)
  if (any(counts != floor(counts)) || N < 1 || N > 2^52) {
    stop("counts must be whole numbers with a sum in [1, 2^52].", call. = FALSE)
  }
  names(counts) <- c("0", "1", "2")
  g <- counts / N
  fit <- .hwe_fit(g)
  expected_counts <- N * fit$expected_probs
  sparse <- any(expected_counts < 5)

  statistic <- N * fit$F^2
  df <- 1L
  if (is.na(fit$F)) {
    warning("The locus is monomorphic; the sample is uninformative for HWE.",
            call. = FALSE)
    df <- NA_integer_
  } else if (method == "chisq" && sparse) {
    warning("An expected genotype count is below 5; consider the exact HWE test.",
            call. = FALSE)
  }
  p_value <- if (method == "exact") .hwe_exact_p(counts, N) else
    stats::pchisq(statistic, df, lower.tail = FALSE)

  list(N = N, counts = counts, observed_probs = g,
       allele_frequency = fit$allele_frequency, expected_probs = fit$expected_probs,
       expected_counts = expected_counts, statistic = statistic, df = df,
       method = method, p_value = p_value, F = fit$F, sparse_expected = sparse)
}

.hwe_alternative <- function(genotype_probs, allele_frequency, misclassification_matrix) {
  true_probs <- NULL
  if (is.null(genotype_probs)) {
    if (is.null(allele_frequency) || is.null(misclassification_matrix)) {
      stop("Supply genotype_probs, or both allele_frequency and misclassification_matrix.",
           call. = FALSE)
    }
    .hwe_validate_scalar(allele_frequency, "allele_frequency", 0, 1)
    M <- misclassification_matrix
    if (!is.numeric(M) || !identical(dim(M), c(3L, 3L)) ||
        any(!is.finite(M)) || any(M < 0) || any(M > 1)) {
      stop("misclassification_matrix must be a numeric 3 by 3 matrix with finite probabilities in [0, 1].",
           call. = FALSE)
    }
    .hwe_validate_order(rownames(M), "Matrix rows (true genotypes)")
    .hwe_validate_order(colnames(M), "Matrix columns (observed genotypes)")
    totals <- rowSums(M)
    if (any(abs(totals - 1) > 1e-12)) {
      stop("misclassification_matrix rows must sum to 1.", call. = FALSE)
    }
    true_probs <- .hwe_proportions(allele_frequency, 1 - allele_frequency)
    genotype_probs <- drop(true_probs %*% (M / totals))
  } else {
    if (!is.null(allele_frequency) || !is.null(misclassification_matrix)) {
      stop("Supply genotype_probs, or both allele_frequency and misclassification_matrix.",
           call. = FALSE)
    }
    .hwe_validate_vector(genotype_probs, "genotype_probs")
    if (abs(sum(genotype_probs) - 1) > 1e-12) {
      stop("genotype_probs must sum to 1.", call. = FALSE)
    }
  }

  observed_probs <- genotype_probs / sum(genotype_probs)
  names(observed_probs) <- c("0", "1", "2")
  fit <- .hwe_fit(observed_probs)
  if (is.na(fit$F)) {
    stop("The observed alternative is monomorphic; Pearson HWE planning requires both alleles.",
         call. = FALSE)
  }
  c(list(true_genotype_probs = true_probs, observed_probs = observed_probs), fit)
}

.hwe_asymptotic_power <- function(N, F, z, alpha) {
  if (F == 0) return(alpha)
  shift <- sqrt(N) * abs(F)
  # A noncentral chi-square with one df is the square of N(shift, 1).
  # Evaluating both normal tails is stable even at large noncentrality.
  max(alpha, stats::pnorm(-z - shift) + stats::pnorm(z - shift, lower.tail = FALSE))
}

.hwe_planning_result <- function(fit, N, alpha, z, power) {
  if (is.finite(N)) {
    expected_counts <- N * fit$expected_probs
    sparse <- any(expected_counts < 5)
    if (sparse) {
      warning("A planned HWE expected count is below 5; asymptotic power can be unreliable.",
              call. = FALSE)
    }
  } else {
    expected_counts <- stats::setNames(rep(NA_real_, 3L), names(fit$expected_probs))
    sparse <- NA
  }
  lambda <- if (is.finite(N)) (sqrt(N) * fit$F)^2 else NA_real_
  c(list(N = N, power = power, alpha = alpha, df = 1L, critical_value = z^2,
         lambda = lambda, lambda_per_person = fit$F^2,
         expected_counts = expected_counts, sparse_expected = sparse), fit)
}

#' Plan power and sample size for a Pearson HWE test
#'
#' Compute prospective power or the minimum whole-number sample size for the
#' ordinary, uncorrected Pearson Hardy-Weinberg equilibrium (HWE) test.
#' This is a noncentral chi-square planning approximation for independently
#' sampled diploid individuals at one biallelic autosomal locus.
#'
#' @param N Number of analyzable individuals, a positive whole number at
#'   most `2^52`.
#' @param alpha Nominal significance level, strictly between 0 and 1.
#' @param genotype_probs Assumed observed probabilities for genotypes with 0,
#'   1, and 2 copies of the same allele. Supply a finite, nonnegative vector
#'   summing to 1, or instead supply both `allele_frequency` and
#'   `misclassification_matrix`. Names, if given, must be `"0"`, `"1"`,
#'   and `"2"` in that order.
#' @param allele_frequency Frequency of the counted allele **before** error,
#'   in `[0, 1]`. It defines the true HWE probabilities.
#' @param misclassification_matrix Numeric 3 by 3 matrix with rows for true
#'   genotype and columns for observed genotype, both in 0, 1, 2 allele-copy
#'   order. Any row or column names must be `"0"`, `"1"`, `"2"` in that order.
#'   Entry `[i + 1, j + 1]` is the probability of calling `j` copies when
#'   there are truly `i` copies. Each row sums to 1. Multiplying the true
#'   probability row vector by this matrix gives the observed probabilities.
#' @param target_power Desired rejection probability, strictly between 0 and 1.
#'
#' @details The fitted null uses the counted-allele frequency from the
#'   assumed observed probabilities `g`. Write `h` for the HWE probabilities
#'   at that frequency. The noncentrality parameter is
#'   `lambda = N * sum((g - h)^2 / h) = N * F^2`, where
#'   `F = 1 - g[2] / h[2]` uses the heterozygote probabilities. Power is the
#'   upper-tail probability under this one-degree-of-freedom
#'   noncentral chi-square at the central chi-square cutoff. No exact-test
#'   power or goodness of fit to a fixed external frequency is calculated.
#'   The fitted distribution must contain both alleles. Expected counts below
#'   five trigger a warning for finite designs. They are a diagnostic for the
#'   approximation, rather than a constraint in the sample-size search.
#'   Sample sizes above `2^52` are outside the supported whole-number range.
#'
#' @return A list with the sample size `N`, achieved `power`, `alpha`,
#'   `df`, `critical_value`, `lambda` (the noncentrality parameter),
#'   `lambda_per_person`, `expected_counts`, `sparse_expected`,
#'   `true_genotype_probs` (matrix scenario only), `observed_probs`,
#'   fitted `allele_frequency`, `expected_probs`, and `F`.
#'   `hwe_mssn()` also returns `target_power` and `previous_power` at the
#'   immediately smaller positive size when available. An equilibrium
#'   alternative has power `alpha` at every size in this approximation.
#'   For such an alternative and target above `alpha`, `N` is `Inf`
#'   and achieved power and size-dependent diagnostics are `NA`.
#'
#' @references
#' Rohlfs RV, Weir BS (2008). Distributions of Hardy-Weinberg equilibrium test
#' statistics. *Genetics* 180(3), 1609-1616.
#' \doi{10.1534/genetics.108.088005}.
#'
#' @seealso [hwe_test()]
#' @examples
#' g <- c(0.35588, 0.48424, 0.15988)
#' hwe_power(1000, genotype_probs = g)$power
#' hwe_mssn(genotype_probs = g)$N
#' @export
hwe_power <- function(N, alpha = 0.05, genotype_probs = NULL,
                      allele_frequency = NULL, misclassification_matrix = NULL) {
  .hwe_validate_scalar(N, "N", 1, 2^52)
  if (N != floor(N)) stop("N must be a whole number.", call. = FALSE)
  .hwe_validate_scalar(alpha, "alpha", 0, 1, open = TRUE)
  fit <- .hwe_alternative(genotype_probs, allele_frequency, misclassification_matrix)
  z <- stats::qnorm(log(alpha) - log(2), lower.tail = FALSE, log.p = TRUE)
  power <- .hwe_asymptotic_power(N, fit$F, z, alpha)
  .hwe_planning_result(fit, N, alpha, z, power)
}

#' @rdname hwe_power
#' @export
hwe_mssn <- function(target_power = 0.80, alpha = 0.05, genotype_probs = NULL,
                     allele_frequency = NULL, misclassification_matrix = NULL) {
  .hwe_validate_scalar(target_power, "target_power", 0, 1, open = TRUE)
  .hwe_validate_scalar(alpha, "alpha", 0, 1, open = TRUE)
  fit <- .hwe_alternative(genotype_probs, allele_frequency, misclassification_matrix)
  z <- stats::qnorm(log(alpha) - log(2), lower.tail = FALSE, log.p = TRUE)
  power_at <- function(N) .hwe_asymptotic_power(N, fit$F, z, alpha)

  N <- 1
  if (power_at(N) < target_power) {
    if (fit$F == 0) {
      res <- .hwe_planning_result(fit, Inf, alpha, z, NA_real_)
      return(c(res, list(target_power = target_power, previous_power = NA_real_)))
    }
    lo <- 1
    hi <- 2^52
    if (power_at(hi) < target_power) {
      stop("Required sample size exceeds the supported whole-number range (2^52).",
           call. = FALSE)
    }
    # The noncentral chi-square upper tail increases with N * F^2.
    while (hi - lo > 1) {
      mid <- lo + floor((hi - lo) / 2)
      if (power_at(mid) >= target_power) hi <- mid else lo <- mid
    }
    N <- hi
  }
  res <- .hwe_planning_result(fit, N, alpha, z, power_at(N))
  c(res, list(target_power = target_power,
              previous_power = if (N > 1) power_at(N - 1) else NA_real_))
}
