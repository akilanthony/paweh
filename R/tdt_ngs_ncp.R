# Private implementation of the published TDT1-NGS noncentrality parameter
# in Kim (2015), Appendix B. The calculation uses latent trio states and raw
# read-count probabilities; it does not classify genotypes and is analytic.

.tdt_ngs_validate_scalar <- function(x, name, lower, upper,
                                     lower_open = FALSE,
                                     upper_open = FALSE,
                                     integer = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (valid) {
    valid <- if (lower_open) x > lower else x >= lower
  }
  if (valid) {
    valid <- if (upper_open) x < upper else x <= upper
  }
  if (valid && integer) {
    valid <- x == floor(x)
  }
  if (!valid) {
    stop(name, " is outside its supported domain.")
  }
  invisible(TRUE)
}

.tdt_ngs_trio_states <- function() {
  score_delta <- c(
    0, -0.5, -0.5, 0.5, 0.5, 0, 0, -1, 0, 1,
    -0.5, -0.5, 0.5, 0.5, 0
  )

  data.frame(
    state = 1:15,
    father = c(2, 1, 2, 1, 2, 0, 2, 1, 1, 1, 0, 1, 0, 1, 0),
    mother = c(2, 2, 1, 2, 1, 2, 0, 1, 1, 1, 1, 0, 1, 0, 0),
    child = c(2, 1, 1, 2, 2, 1, 1, 0, 1, 2, 0, 0, 1, 1, 0),
    mating = c(
      "mu22", "mu12", "mu21", "mu12", "mu21", "mu02", "mu20",
      "mu11", "mu11", "mu11", "mu01", "mu10", "mu01", "mu10",
      "mu00"
    ),
    score_delta = score_delta,
    stringsAsFactors = FALSE
  )
}

.tdt_ngs_hwe_mating_freqs <- function(pd) {
  .tdt_ngs_validate_scalar(pd, "pd", 0, 1, TRUE, TRUE)

  genotype <- c(g0 = (1 - pd)^2, g1 = 2 * pd * (1 - pd), g2 = pd^2)
  mu <- c(
    mu22 = genotype[["g2"]] * genotype[["g2"]],
    mu12 = genotype[["g1"]] * genotype[["g2"]],
    mu21 = genotype[["g2"]] * genotype[["g1"]],
    mu02 = genotype[["g0"]] * genotype[["g2"]],
    mu20 = genotype[["g2"]] * genotype[["g0"]],
    mu11 = genotype[["g1"]] * genotype[["g1"]],
    mu01 = genotype[["g0"]] * genotype[["g1"]],
    mu10 = genotype[["g1"]] * genotype[["g0"]],
    mu00 = genotype[["g0"]] * genotype[["g0"]]
  )

  if (any(!is.finite(mu)) || any(mu <= 0) || abs(sum(mu) - 1) > 1e-12) {
    stop("HWE mating frequencies must be positive and sum to 1.")
  }
  mu
}

.tdt_ngs_validate_mating_freqs <- function(mu) {
  required <- c(
    "mu22", "mu12", "mu21", "mu02", "mu20",
    "mu11", "mu01", "mu10", "mu00"
  )
  if (!is.numeric(mu) || length(mu) != 9L || any(!is.finite(mu)) ||
      any(mu <= 0) || abs(sum(mu) - 1) > 1e-8) {
    stop("mu must contain nine positive finite mating probabilities summing to 1.")
  }
  if (is.null(names(mu)) || !identical(names(mu), required)) {
    stop("mu must use Kim's order: mu22, mu12, mu21, mu02, mu20, ",
         "mu11, mu01, mu10, mu00.")
  }
  invisible(TRUE)
}

.tdt_ngs_state_probabilities <- function(mu, t) {
  .tdt_ngs_validate_mating_freqs(mu)
  .tdt_ngs_validate_scalar(t, "t", 0, 1, TRUE, TRUE)

  transmission <- c(
    1,
    1 - t, 1 - t,
    t, t,
    1, 1,
    (1 - t)^2, 2 * t * (1 - t), t^2,
    1 - t, 1 - t,
    t, t,
    1
  )
  states <- .tdt_ngs_trio_states()
  pi <- unname(mu[states$mating]) * transmission
  names(pi) <- paste0("state_", states$state)

  if (any(!is.finite(pi)) || any(pi < 0) || abs(sum(pi) - 1) > 1e-12) {
    stop("TDT1-NGS trio-state probabilities must be finite, nonnegative, and sum to 1.")
  }
  pi
}

.tdt_ngs_read_probability <- function(genotype, epsilon0, epsilon1) {
  if (!is.numeric(genotype) || length(genotype) != 1L ||
      !is.finite(genotype) || !(genotype %in% 0:2)) {
    stop("genotype must be a single value in {0, 1, 2}.")
  }
  .tdt_ngs_validate_scalar(epsilon0, "epsilon0", 0, 1, FALSE, TRUE)
  .tdt_ngs_validate_scalar(epsilon1, "epsilon1", 0, 1, FALSE, TRUE)
  if (epsilon0 + epsilon1 >= 1) {
    stop("epsilon0 + epsilon1 must be less than 1.")
  }

  epsilon0 + ((1 - epsilon0 - epsilon1) / 2) * genotype
}

.tdt_ngs_read_likelihood <- function(x, coverage, genotypes,
                                     epsilon0, epsilon1, log = FALSE) {
  .tdt_ngs_validate_scalar(
    coverage, "coverage", 1, Inf, FALSE, FALSE, integer = TRUE
  )
  if (!is.numeric(x) || length(x) != 3L || any(!is.finite(x)) ||
      any(x < 0) || any(x > coverage) || any(x != floor(x))) {
    stop("x must contain three integer read counts between 0 and coverage.")
  }
  if (!is.numeric(genotypes) || length(genotypes) != 3L ||
      any(!is.finite(genotypes)) || any(!(genotypes %in% 0:2))) {
    stop("genotypes must contain three values in {0, 1, 2}.")
  }

  q <- vapply(
    genotypes, .tdt_ngs_read_probability, numeric(1),
    epsilon0 = epsilon0, epsilon1 = epsilon1
  )
  log_probability <- sum(stats::dbinom(x, coverage, q, log = TRUE))
  if (log) log_probability else exp(log_probability)
}

.tdt_ngs_log_sum_exp <- function(x) {
  maximum <- max(x)
  if (!is.finite(maximum)) {
    return(maximum)
  }
  maximum + log(sum(exp(x - maximum)))
}

.tdt_ngs_posterior <- function(x, coverage, pi, epsilon0, epsilon1) {
  if (!is.numeric(pi) || length(pi) != 15L || any(!is.finite(pi)) ||
      any(pi <= 0) || abs(sum(pi) - 1) > 1e-8) {
    stop("pi must contain 15 positive finite probabilities summing to 1.")
  }
  states <- .tdt_ngs_trio_states()
  log_weights <- vapply(seq_len(nrow(states)), function(i) {
    log(pi[[i]]) + .tdt_ngs_read_likelihood(
      x = x,
      coverage = coverage,
      genotypes = unlist(states[i, c("father", "mother", "child")]),
      epsilon0 = epsilon0,
      epsilon1 = epsilon1,
      log = TRUE
    )
  }, numeric(1))
  log_h <- .tdt_ngs_log_sum_exp(log_weights)
  tau <- exp(log_weights - log_h)

  if (any(!is.finite(tau)) || any(tau < 0) || abs(sum(tau) - 1) > 1e-10) {
    stop("Posterior trio-state probabilities are numerically invalid.")
  }
  names(tau) <- paste0("state_", states$state)
  tau
}

.tdt_ngs_mating_score_matrix <- function(mu) {
  .tdt_ngs_validate_mating_freqs(mu)
  states <- .tdt_ngs_trio_states()
  free_names <- names(mu)[1:8]
  score <- matrix(
    0, nrow = 15, ncol = 8,
    dimnames = list(paste0("state_", 1:15), free_names)
  )
  for (j in seq_along(free_names)) {
    score[states$mating == free_names[[j]], j] <- 1 / mu[[j]]
  }
  score[states$mating == "mu00", ] <- -1 / mu[["mu00"]]
  score
}

.tdt_ngs_read_score <- function(x, coverage, genotypes,
                                epsilon0, epsilon1) {
  likelihood_args <- list(
    x = x, coverage = coverage, genotypes = genotypes,
    epsilon0 = epsilon0, epsilon1 = epsilon1, log = TRUE
  )
  do.call(.tdt_ngs_read_likelihood, likelihood_args)

  q <- vapply(
    genotypes, .tdt_ngs_read_probability, numeric(1),
    epsilon0 = epsilon0, epsilon1 = epsilon1
  )
  if (any(q <= 0) || any(q >= 1)) {
    stop("Component read scores are defined only for interior error probabilities.")
  }
  base <- x / q - (coverage - x) / (1 - q)
  c(
    epsilon0 = sum(base * (1 - genotypes / 2)),
    epsilon1 = sum(base * (-genotypes / 2))
  )
}

.tdt_ngs_binomial_pmf_derivative <- function(x, size, probability) {
  pmf <- stats::dbinom(x, size = size, prob = probability)
  if (probability == 0) {
    derivative <- numeric(length(x))
    derivative[x == 0] <- -size
    derivative[x == 1] <- size
    return(list(pmf = pmf, derivative = derivative))
  }
  if (probability == 1) {
    derivative <- numeric(length(x))
    derivative[x == size] <- size
    derivative[x == size - 1] <- -size
    return(list(pmf = pmf, derivative = derivative))
  }
  list(
    pmf = pmf,
    derivative = pmf * (x / probability - (size - x) / (1 - probability))
  )
}

.tdt_ngs_information_matrix <- function(pd, coverage, seq_error) {
  .tdt_ngs_validate_scalar(pd, "pd", 0, 1, TRUE, TRUE)
  .tdt_ngs_validate_scalar(
    coverage, "coverage", 1, Inf, FALSE, FALSE, integer = TRUE
  )
  .tdt_ngs_validate_scalar(seq_error, "seq_error", 0, 0.5, FALSE, TRUE)

  states <- .tdt_ngs_trio_states()
  mu <- .tdt_ngs_hwe_mating_freqs(pd)
  pi <- .tdt_ngs_state_probabilities(mu, t = 0.5)
  counts <- 0:coverage
  grid <- expand.grid(
    father = counts, mother = counts, child = counts,
    KEEP.OUT.ATTRS = FALSE
  )
  member_columns <- c("father", "mother", "child")

  q <- vapply(
    0:2, .tdt_ngs_read_probability, numeric(1),
    epsilon0 = seq_error, epsilon1 = seq_error
  )
  pmf <- dq <- matrix(0, nrow = coverage + 1L, ncol = 3L)
  for (g in 0:2) {
    values <- .tdt_ngs_binomial_pmf_derivative(counts, coverage, q[[g + 1L]])
    pmf[, g + 1L] <- values$pmf
    dq[, g + 1L] <- values$derivative
  }
  dq0 <- sweep(dq, 2, 1 - (0:2) / 2, `*`)
  dq1 <- sweep(dq, 2, -(0:2) / 2, `*`)

  n_read_states <- nrow(grid)
  indices <- Map(function(column) grid[[column]] + 1L, member_columns)
  state_likelihood <- matrix(0, nrow = n_read_states, ncol = 15L)
  derivative0 <- derivative1 <- matrix(
    0, nrow = n_read_states, ncol = 15L
  )
  for (i in seq_len(nrow(states))) {
    g <- unlist(states[i, member_columns])
    p <- Map(function(index, genotype) pmf[index, genotype + 1L], indices, g)
    d0 <- Map(function(index, genotype) dq0[index, genotype + 1L], indices, g)
    d1 <- Map(function(index, genotype) dq1[index, genotype + 1L], indices, g)

    state_likelihood[, i] <- p[[1]] * p[[2]] * p[[3]]
    derivative0[, i] <-
      d0[[1]] * p[[2]] * p[[3]] +
      p[[1]] * d0[[2]] * p[[3]] +
      p[[1]] * p[[2]] * d0[[3]]
    derivative1[, i] <-
      d1[[1]] * p[[2]] * p[[3]] +
      p[[1]] * d1[[2]] * p[[3]] +
      p[[1]] * p[[2]] * d1[[3]]
  }

  weighted <- sweep(state_likelihood, 2, pi, `*`)
  h <- rowSums(weighted)
  if (any(!is.finite(h)) || any(h <= 0) || abs(sum(h) - 1) > 1e-10) {
    stop("The exact null read-count distribution is numerically invalid.")
  }
  tau <- weighted / h

  scores <- cbind(
    delta = as.numeric(tau %*% states$score_delta),
    tau %*% .tdt_ngs_mating_score_matrix(mu),
    epsilon0 = as.numeric(derivative0 %*% pi) / h,
    epsilon1 = as.numeric(derivative1 %*% pi) / h
  )
  information <- crossprod(scores * sqrt(h))
  information <- (information + t(information)) / 2

  if (any(!is.finite(information))) {
    stop("The exact TDT1-NGS information matrix contains nonfinite values.")
  }
  list(
    information_matrix = information,
    parameter_names = colnames(information),
    mu = mu,
    pi = pi,
    null_read_probability_sum = sum(h),
    score_mean = colSums(scores * h),
    coverage = coverage,
    seq_error = seq_error
  )
}

.tdt_ngs_information <- function(pd, coverage, seq_error) {
  result <- .tdt_ngs_information_matrix(pd, coverage, seq_error)
  information <- result$information_matrix
  nuisance <- information[-1, -1, drop = FALSE]
  nuisance_scale <- sqrt(diag(nuisance))
  scaled_nuisance <- nuisance / outer(nuisance_scale, nuisance_scale)
  scaled_rcond <- rcond(scaled_nuisance)

  if (any(!is.finite(nuisance_scale)) || any(nuisance_scale <= 0) ||
      !is.finite(scaled_rcond) || scaled_rcond < 1e-12) {
    stop(
      "The published TDT1-NGS nuisance information matrix is singular for ",
      "this design; efficient information cannot be computed without changing ",
      "the Appendix B parameterization."
    )
  }
  cross <- information[1, -1, drop = FALSE]
  scaled_cross <- as.numeric(cross) / nuisance_scale
  efficient <- as.numeric(
    information[1, 1] -
      scaled_cross %*% solve(scaled_nuisance, scaled_cross)
  )
  if (!is.finite(efficient) || efficient < -1e-10) {
    stop("TDT1-NGS efficient information must be finite and nonnegative.")
  }
  efficient <- max(efficient, 0)

  c(result, list(
    I_dd = unname(information[1, 1]),
    I_deta = unname(cross),
    I_etaeta = nuisance,
    nuisance_rcond = scaled_rcond,
    efficient_information = efficient
  ))
}

.tdt_ngs_ncp <- function(N, pd, R1, coverage, seq_error) {
  .tdt_ngs_validate_scalar(N, "N", 1, Inf, FALSE, FALSE, integer = TRUE)
  .tdt_ngs_validate_scalar(pd, "pd", 0, 1, TRUE, TRUE)
  .tdt_ngs_validate_scalar(R1, "R1", 0, Inf, TRUE, FALSE)
  .tdt_ngs_validate_scalar(
    coverage, "coverage", 1, Inf, FALSE, FALSE, integer = TRUE
  )
  .tdt_ngs_validate_scalar(seq_error, "seq_error", 0, 0.5, FALSE, TRUE)

  t <- R1 / (1 + R1)
  delta <- log(t / (1 - t))
  information <- .tdt_ngs_information(pd, coverage, seq_error)
  lambda <- N * delta^2 * information$efficient_information

  if (!is.finite(lambda) || lambda < 0) {
    stop("The published TDT1-NGS NCP must be finite and nonnegative.")
  }
  list(
    lambda = unname(lambda),
    N = N,
    pd = pd,
    R1 = R1,
    R2 = R1^2,
    t = t,
    delta = delta,
    coverage = coverage,
    seq_error = seq_error,
    efficient_information = information$efficient_information,
    information_matrix = information$information_matrix,
    nuisance_information = information$I_etaeta,
    nuisance_rcond = information$nuisance_rcond,
    score_mean = information$score_mean
  )
}
