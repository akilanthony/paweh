kim_states <- data.frame(
  state = 1:15,
  father = c(2, 1, 2, 1, 2, 0, 2, 1, 1, 1, 0, 1, 0, 1, 0),
  mother = c(2, 2, 1, 2, 1, 2, 0, 1, 1, 1, 1, 0, 1, 0, 0),
  child = c(2, 1, 1, 2, 2, 1, 1, 0, 1, 2, 0, 0, 1, 1, 0),
  mating = c(
    "mu22", "mu12", "mu21", "mu12", "mu21", "mu02", "mu20",
    "mu11", "mu11", "mu11", "mu01", "mu10", "mu01", "mu10",
    "mu00"
  ),
  score_delta = c(
    0, -0.5, -0.5, 0.5, 0.5, 0, 0, -1, 0, 1,
    -0.5, -0.5, 0.5, 0.5, 0
  ),
  stringsAsFactors = FALSE
)

kim_reference_mu <- function(pd) {
  g <- c((1 - pd)^2, 2 * pd * (1 - pd), pd^2)
  c(
    mu22 = g[3] * g[3], mu12 = g[2] * g[3], mu21 = g[3] * g[2],
    mu02 = g[1] * g[3], mu20 = g[3] * g[1], mu11 = g[2] * g[2],
    mu01 = g[1] * g[2], mu10 = g[2] * g[1], mu00 = g[1] * g[1]
  )
}

kim_reference_pi <- function(mu, delta) {
  t <- stats::plogis(delta)
  child_probability <- c(
    1, 1 - t, 1 - t, t, t, 1, 1,
    (1 - t)^2, 2 * t * (1 - t), t^2,
    1 - t, 1 - t, t, t, 1
  )
  unname(mu[kim_states$mating]) * child_probability
}

kim_reference_h <- function(x, coverage, delta, mu_free, epsilon0, epsilon1) {
  mu <- c(mu_free, mu00 = 1 - sum(mu_free))
  pi <- kim_reference_pi(mu, delta)
  likelihood <- vapply(seq_len(nrow(kim_states)), function(i) {
    genotypes <- unlist(kim_states[i, c("father", "mother", "child")])
    q <- epsilon0 + (1 - epsilon0 - epsilon1) * genotypes / 2
    prod(stats::dbinom(x, coverage, q))
  }, numeric(1))
  sum(pi * likelihood)
}

kim_reference_information_fd <- function(pd, coverage, seq_error,
                                         step = 1e-6) {
  mu <- kim_reference_mu(pd)
  theta <- c(delta = 0, mu[1:8], epsilon0 = seq_error,
             epsilon1 = seq_error)
  grid <- expand.grid(
    father = 0:coverage, mother = 0:coverage, child = 0:coverage,
    KEEP.OUT.ATTRS = FALSE
  )
  score <- matrix(0, nrow(grid), length(theta))
  h0 <- numeric(nrow(grid))

  evaluate <- function(x, parameter) {
    kim_reference_h(
      x, coverage, parameter[["delta"]], parameter[2:9],
      parameter[["epsilon0"]], parameter[["epsilon1"]]
    )
  }
  for (row in seq_len(nrow(grid))) {
    x <- unlist(grid[row, ])
    h0[[row]] <- evaluate(x, theta)
    for (j in seq_along(theta)) {
      upper <- lower <- theta
      upper[[j]] <- upper[[j]] + step
      lower[[j]] <- lower[[j]] - step
      score[row, j] <-
        (log(evaluate(x, upper)) - log(evaluate(x, lower))) / (2 * step)
    }
  }
  crossprod(score * sqrt(h0))
}

test_that("Kim's published 15-state ordering and TDT1 score are exact", {
  observed <- .tdt_ngs_trio_states()

  expect_identical(observed, kim_states)
  expect_identical(observed$state, 1:15)
  expect_equal(nrow(unique(observed[c("father", "mother", "child")])), 15)
  expect_true(all(unlist(observed[c("father", "mother", "child")]) %in% 0:2))

  mendelian <- function(father, mother) {
    gametes <- function(g) if (g == 0) 0 else if (g == 2) 1 else c(0, 1)
    unique(as.vector(outer(gametes(father), gametes(mother), `+`)))
  }
  for (i in seq_len(nrow(observed))) {
    expect_true(observed$child[[i]] %in%
                  mendelian(observed$father[[i]], observed$mother[[i]]))
  }
})

test_that("HWE and Kim-ordered mating frequencies are exact", {
  for (pd in c(0.15, 0.325, 0.5)) {
    observed <- .tdt_ngs_hwe_mating_freqs(pd)
    expected <- kim_reference_mu(pd)
    expect_equal(observed, expected, tolerance = 1e-15)
    expect_equal(sum(observed), 1, tolerance = 1e-15)
    expect_true(all(observed > 0))
    expect_equal(observed[["mu12"]], observed[["mu21"]])
    expect_equal(observed[["mu02"]], observed[["mu20"]])
    expect_equal(observed[["mu01"]], observed[["mu10"]])
    expect_equal(observed[["mu00"]], (1 - pd)^4)
  }
})

test_that("all 15 trio-state probabilities follow the Appendix A formulas", {
  mu <- kim_reference_mu(0.325)
  for (t in c(0.35, 0.5, 0.65)) {
    observed <- .tdt_ngs_state_probabilities(mu, t)
    expected <- kim_reference_pi(mu, stats::qlogis(t))
    expect_equal(unname(observed), expected, tolerance = 1e-15)
    expect_equal(sum(observed), 1, tolerance = 1e-15)
    expect_true(all(observed >= 0))
  }
})

test_that("null trio-state probabilities have the published transmission symmetry", {
  pi <- .tdt_ngs_state_probabilities(kim_reference_mu(0.325), 0.5)
  expect_equal(pi[[2]], pi[[4]])
  expect_equal(pi[[3]], pi[[5]])
  expect_equal(pi[[8]], pi[[10]])
  expect_equal(pi[[9]], 2 * pi[[8]])
  expect_equal(pi[[11]], pi[[13]])
  expect_equal(pi[[12]], pi[[14]])
})

test_that("the Appendix B transmission score is the delta derivative", {
  mu <- kim_reference_mu(0.325)
  step <- 1e-6
  upper <- log(kim_reference_pi(mu, step))
  lower <- log(kim_reference_pi(mu, -step))
  finite_difference <- (upper - lower) / (2 * step)

  expect_equal(finite_difference, kim_states$score_delta, tolerance = 2e-10)
})

test_that("the eight mating scores match constrained finite differences", {
  mu <- kim_reference_mu(0.325)
  score <- .tdt_ngs_mating_score_matrix(mu)
  step <- 1e-7

  for (j in 1:8) {
    upper <- lower <- mu
    upper[[j]] <- upper[[j]] + step
    upper[[9]] <- upper[[9]] - step
    lower[[j]] <- lower[[j]] - step
    lower[[9]] <- lower[[9]] + step
    finite_difference <-
      (log(kim_reference_pi(upper, 0)) -
         log(kim_reference_pi(lower, 0))) / (2 * step)
    expect_equal(finite_difference, unname(score[, j]), tolerance = 2e-6)
  }
})

test_that("multiplicative transmission gives delta equal to log R1", {
  for (R1 in c(1, 1.1, 1.2, 1.3)) {
    t <- R1 / (1 + R1)
    expect_equal(log(t / (1 - t)), log(R1), tolerance = 1e-14)
    expect_equal(t, R1^2 / (R1 + R1^2), tolerance = 1e-15)
  }
})

test_that("trio read likelihoods equal direct binomial products", {
  settings <- list(
    list(x = c(0, 2, 4), g = c(0, 1, 2), v = 4, e0 = 0.01, e1 = 0.02),
    list(x = c(3, 1, 2), g = c(1, 2, 1), v = 5, e0 = 0.03, e1 = 0.01)
  )
  for (setting in settings) {
    q <- setting$e0 + (1 - setting$e0 - setting$e1) * setting$g / 2
    expected <- prod(stats::dbinom(setting$x, setting$v, q))
    observed <- .tdt_ngs_read_likelihood(
      setting$x, setting$v, setting$g, setting$e0, setting$e1
    )
    expect_equal(observed, expected, tolerance = 1e-13)
  }
})

test_that("analytic directional read scores match finite differences", {
  x <- c(1, 3, 4)
  genotypes <- c(0, 1, 2)
  coverage <- 5
  epsilon0 <- 0.02
  epsilon1 <- 0.03
  step <- 1e-7
  analytic <- .tdt_ngs_read_score(
    x, coverage, genotypes, epsilon0, epsilon1
  )
  log_f <- function(e0, e1) .tdt_ngs_read_likelihood(
    x, coverage, genotypes, e0, e1, log = TRUE
  )
  finite_difference <- c(
    epsilon0 = (log_f(epsilon0 + step, epsilon1) -
      log_f(epsilon0 - step, epsilon1)) / (2 * step),
    epsilon1 = (log_f(epsilon0, epsilon1 + step) -
      log_f(epsilon0, epsilon1 - step)) / (2 * step)
  )
  expect_equal(analytic, finite_difference, tolerance = 2e-8)
})

test_that("posterior state probabilities satisfy Bayes' formula", {
  mu <- kim_reference_mu(0.325)
  pi <- .tdt_ngs_state_probabilities(mu, 0.5)
  x <- c(1, 2, 3)
  coverage <- 4
  epsilon0 <- 0.01
  epsilon1 <- 0.02

  likelihood <- vapply(seq_len(nrow(kim_states)), function(i) {
    g <- unlist(kim_states[i, c("father", "mother", "child")])
    q <- epsilon0 + (1 - epsilon0 - epsilon1) * g / 2
    prod(stats::dbinom(x, coverage, q))
  }, numeric(1))
  expected <- unname(pi * likelihood / sum(pi * likelihood))
  observed <- .tdt_ngs_posterior(
    x, coverage, pi, epsilon0, epsilon1
  )
  expect_equal(unname(observed), expected, tolerance = 1e-14)
  expect_equal(sum(observed), 1, tolerance = 1e-14)
  expect_true(all(observed >= 0))
})

test_that("log-sum-exp posterior evaluation remains stable at high depth", {
  pi <- .tdt_ngs_state_probabilities(kim_reference_mu(0.325), 0.5)
  for (coverage in c(44, 100)) {
    x <- c(1, floor(coverage / 2), coverage - 1)
    tau <- .tdt_ngs_posterior(x, coverage, pi, 0.005, 0.005)
    expect_true(all(is.finite(tau)))
    expect_true(all(tau >= 0))
    expect_equal(sum(tau), 1, tolerance = 1e-13)
  }
})

test_that("exact information agrees with independent finite-difference enumeration", {
  for (coverage in c(1, 2, 4)) {
    production <- .tdt_ngs_information_matrix(0.325, coverage, 0.01)
    reference <- kim_reference_information_fd(0.325, coverage, 0.01)
    dimnames(reference) <- dimnames(production$information_matrix)
    expect_equal(
      production$information_matrix, reference,
      tolerance = 2e-7
    )
    expect_equal(production$null_read_probability_sum, 1, tolerance = 1e-14)
    expect_lt(max(abs(production$score_mean)), 1e-10)
  }
})

test_that("coverage one exposes the published nuisance-rank limitation", {
  information <- .tdt_ngs_information_matrix(0.325, 1, 0.01)
  expect_equal(dim(information$information_matrix), c(11, 11))
  expect_lte(qr(information$information_matrix)$rank, 7)
  expect_error(
    .tdt_ngs_information(0.325, 1, 0.01),
    "nuisance information matrix is singular"
  )
})

test_that("efficient information is the Appendix B Schur complement", {
  observed <- .tdt_ngs_information(0.325, 4, 0.01)
  matrix <- observed$information_matrix
  nuisance <- matrix[-1, -1, drop = FALSE]
  cross <- matrix[1, -1, drop = FALSE]
  expected <- as.numeric(matrix[1, 1] - cross %*% solve(nuisance, t(cross)))

  expect_equal(observed$efficient_information, expected, tolerance = 1e-12)
  expect_equal(dim(matrix), c(11, 11))
  expect_equal(matrix, t(matrix), tolerance = 1e-14)
  expect_true(all(diag(matrix) >= 0))
})

test_that("the private NCP is exactly N delta squared efficient information", {
  result <- .tdt_ngs_ncp(5000, 0.325, 1.2, 12, 0.005)
  expect_equal(result$delta, log(1.2), tolerance = 1e-14)
  expect_equal(result$t, 1.2 / 2.2, tolerance = 1e-15)
  expect_equal(result$R2, 1.2^2, tolerance = 1e-15)
  expect_equal(
    result$lambda,
    5000 * log(1.2)^2 * result$efficient_information,
    tolerance = 1e-14
  )
})

test_that("the null NCP is zero while information remains positive", {
  result <- .tdt_ngs_ncp(5000, 0.325, 1, 12, 0.005)
  expect_identical(result$delta, 0)
  expect_identical(result$lambda, 0)
  expect_gt(result$efficient_information, 0)
})

test_that("the published NCP scales exactly linearly with trio count", {
  base <- .tdt_ngs_ncp(1000, 0.325, 1.2, 12, 0.005)$lambda
  expect_equal(
    .tdt_ngs_ncp(2000, 0.325, 1.2, 12, 0.005)$lambda,
    2 * base, tolerance = 1e-13
  )
  expect_equal(
    .tdt_ngs_ncp(5000, 0.325, 1.2, 12, 0.005)$lambda,
    5 * base, tolerance = 1e-13
  )
})

test_that("Chapter 5 parameter-domain smoke settings are finite", {
  designs <- list(
    c(pd = 0.15, R1 = 1.1, coverage = 4, seq_error = 0.005),
    c(pd = 0.325, R1 = 1.2, coverage = 12, seq_error = 0.01),
    c(pd = 0.5, R1 = 1.3, coverage = 20, seq_error = 0.005),
    c(pd = 0.325, R1 = 1.2, coverage = 44, seq_error = 0.01)
  )
  for (design in designs) {
    expect_warning(
      result <- .tdt_ngs_ncp(
        5000, design[["pd"]], design[["R1"]],
        design[["coverage"]], design[["seq_error"]]
      ),
      NA
    )
    expect_true(is.finite(result$lambda))
    expect_gte(result$lambda, 0)
    expect_true(all(is.finite(result$information_matrix)))
    expect_gte(result$efficient_information, 0)
  }
})

test_that("coverage and sequencing-error evaluations are numerically stable", {
  for (coverage in c(4, 12, 20, 44)) {
    for (seq_error in c(0, 0.005, 0.01)) {
      result <- .tdt_ngs_ncp(5000, 0.325, 1.2, coverage, seq_error)
      expect_true(is.finite(result$lambda))
      expect_true(is.finite(result$efficient_information))
      expect_gte(result$efficient_information, 0)
      expect_lt(max(abs(result$score_mean)), 1e-8)
    }
  }
})

test_that("allele-frequency reflection symmetry holds under the null information model", {
  for (pd in c(0.15, 0.325)) {
    left <- .tdt_ngs_ncp(5000, pd, 1.2, 12, 0.005)$lambda
    right <- .tdt_ngs_ncp(5000, 1 - pd, 1.2, 12, 0.005)$lambda
    expect_equal(left, right, tolerance = 1e-12)
  }
})

test_that("invalid TDT1-NGS inputs fail informatively", {
  for (pd in list(0, 1, -0.1, NA_real_, Inf, c(0.2, 0.3))) {
    expect_error(.tdt_ngs_ncp(100, pd, 1.2, 4, 0.01), "pd")
  }
  for (N in list(0, -1, 1.5, NA_real_, Inf)) {
    expect_error(.tdt_ngs_ncp(N, 0.3, 1.2, 4, 0.01), "N")
  }
  for (R1 in list(0, -1, NA_real_, Inf, c(1.1, 1.2))) {
    expect_error(.tdt_ngs_ncp(100, 0.3, R1, 4, 0.01), "R1")
  }
  for (coverage in list(0, -1, 1.5, NA_real_, Inf)) {
    expect_error(.tdt_ngs_ncp(100, 0.3, 1.2, coverage, 0.01), "coverage")
  }
  for (seq_error in list(-0.01, 0.5, NA_real_, Inf, c(0.01, 0.02))) {
    expect_error(.tdt_ngs_ncp(100, 0.3, 1.2, 4, seq_error), "seq_error")
  }
})
