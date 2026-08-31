# Primary source:
# Ahn et al. (2007), Annals of Human Genetics 71:249-261.
# DOI: 10.1111/j.1469-1809.2006.00318.x

s7_ahn_reference_ncp <- function(g_case, g_control, N_case, N_control,
                                 scores) {
  weighted_counts <- N_case * g_case + N_control * g_control
  N_case * N_control * sum(scores * (g_case - g_control))^2 /
    (
      sum(scores^2 * weighted_counts) -
        sum(scores * weighted_counts)^2 / (N_case + N_control)
    )
}

s7_integer_design_power <- function(result, N_case) {
  N_control <- ceiling(result$k * N_case)
  lambda <- .cc_ahn_trend_ncp(
    g_case = result$freqs$case_called,
    g_control = result$freqs$control_called,
    N_case = N_case,
    N_control = N_control,
    scores = result$scores
  )
  list(
    N_control = N_control,
    lambda = lambda,
    power = pchisq(
      qchisq(1 - result$alpha, df = 1),
      df = 1,
      ncp = lambda,
      lower.tail = FALSE
    )
  )
}

# Literature metadata, not a PAWEH simulation result. Ahn et al. used 10,000
# simulation replicates per setting and reported close agreement between
# asymptotic and simulated power. PAWEH does not reproduce that simulation here.
s7_ahn_validation_metadata <- list(
  simulation_replicates_per_setting = 10000L,
  reported_max_difference_percentage_points_approximately = 0.3,
  reported_regression_r_squared_exceeded = 0.999,
  paweh_runs_simulation = FALSE,
  paweh_claims_full_ahn_table_reproduction = FALSE,
  reason_full_tables_are_out_of_scope =
    "Public CC-NGS does not separately parameterize marker p2 or LD D/Dprime."
)

test_that("Ahn et al. Equation (1) is locked to a transparent reference", {
  # Ahn et al., Eq. (1)
  settings <- list(
    list(
      g_case = c(0.45, 0.40, 0.15),
      g_control = c(0.55, 0.35, 0.10),
      N_case = 1000,
      N_control = 1200,
      scores = c(0, 1, 2)
    ),
    list(
      g_case = c(0.25, 0.50, 0.25),
      g_control = c(0.36, 0.48, 0.16),
      N_case = 750,
      N_control = 500,
      scores = c(0, 1, 1)
    )
  )

  for (setting in settings) {
    expected <- do.call(s7_ahn_reference_ncp, setting)
    observed <- do.call(.cc_ahn_trend_ncp, setting)
    expect_equal(observed, expected, tolerance = 1e-14)
  }
})

test_that("Ahn row-true error algebra gives adjusted genotype probabilities", {
  # General Ahn-style orientation: rows=true genotype, columns=observed.
  E <- rbind(
    true_AA = c(observed_AA = 0.97, observed_AB = 0.02, observed_BB = 0.01),
    true_AB = c(observed_AA = 0.03, observed_AB = 0.94, observed_BB = 0.03),
    true_BB = c(observed_AA = 0.005, observed_AB = 0.015, observed_BB = 0.98)
  )
  P_A <- c(0.45, 0.40, 0.15)
  P_U <- c(0.55, 0.35, 0.10)
  P_A_star_reference <- c(
    E[1, 1] * P_A[1] + E[2, 1] * P_A[2] + E[3, 1] * P_A[3],
    E[1, 2] * P_A[1] + E[2, 2] * P_A[2] + E[3, 2] * P_A[3],
    E[1, 3] * P_A[1] + E[2, 3] * P_A[2] + E[3, 3] * P_A[3]
  )
  P_U_star_reference <- c(
    E[1, 1] * P_U[1] + E[2, 1] * P_U[2] + E[3, 1] * P_U[3],
    E[1, 2] * P_U[1] + E[2, 2] * P_U[2] + E[3, 2] * P_U[3],
    E[1, 3] * P_U[1] + E[2, 3] * P_U[2] + E[3, 3] * P_U[3]
  )
  P_A_star <- as.numeric(t(E) %*% P_A)
  P_U_star <- as.numeric(t(E) %*% P_U)

  expect_equal(unname(rowSums(E)), rep(1, 3), tolerance = 1e-15)
  expect_equal(P_A_star, P_A_star_reference, tolerance = 1e-15)
  expect_equal(P_U_star, P_U_star_reference, tolerance = 1e-15)
  expect_equal(sum(P_A_star), 1, tolerance = 1e-15)
  expect_equal(sum(P_U_star), 1, tolerance = 1e-15)

  expected_ncp <- s7_ahn_reference_ncp(
    P_A_star, P_U_star, 1000, 1200, c(0, 1, 2)
  )
  observed_ncp <- .cc_ahn_trend_ncp(
    P_A_star, P_U_star, 1000, 1200, c(0, 1, 2)
  )
  expect_equal(observed_ncp, expected_ncp, tolerance = 1e-14)
})

test_that("analytic Ahn MSSN rearrangement is locked", {
  P_A_star <- c(0.44925, 0.38725, 0.16350)
  P_U_star <- c(0.54450, 0.34150, 0.11400)
  scores <- c(0, 1, 2)
  k <- 1.4
  lambda_target <- 7.84883445809286
  D <- sum(scores * (P_A_star - P_U_star))
  Q <- sum(scores^2 * (P_A_star + k * P_U_star)) -
    sum(scores * (P_A_star + k * P_U_star))^2 / (1 + k)
  expected_cases <- lambda_target * Q / (k * D^2)
  observed <- .cc_ngs_mssn_components(
    P_A_star, P_U_star, k, scores, lambda_target
  )

  expect_equal(observed$D, D, tolerance = 1e-15)
  expect_equal(observed$Q, Q, tolerance = 1e-15)
  expect_equal(observed$N_case_continuous, expected_cases, tolerance = 1e-12)
})

test_that("professor pooled-theta form equals Ahn Equation (1)", {
  settings <- list(
    list(c(0.45, 0.40, 0.15), c(0.55, 0.35, 0.10), 1000, 1200,
         c(0, 1, 2)),
    list(c(0.25, 0.50, 0.25), c(0.36, 0.48, 0.16), 500, 1500,
         c(0, 0, 1)),
    list(c(0.70, 0.25, 0.05), c(0.60, 0.30, 0.10), 1800, 600,
         c(-1, 0, 1))
  )

  for (setting in settings) {
    g_case <- setting[[1]]
    g_control <- setting[[2]]
    N_case <- setting[[3]]
    N_control <- setting[[4]]
    scores <- setting[[5]]
    N <- N_case + N_control
    theta <- N_case / N
    pbar <- theta * g_case + (1 - theta) * g_control
    pooled <- N * theta * (1 - theta) *
      sum(scores * (g_case - g_control))^2 /
      (sum(scores^2 * pbar) - sum(scores * pbar)^2)
    ahn <- s7_ahn_reference_ncp(
      g_case, g_control, N_case, N_control, scores
    )

    expect_equal(ahn, pooled, tolerance = 1e-12)
    expect_equal(
      .cc_ahn_trend_ncp(
        g_case, g_control, N_case, N_control, scores
      ),
      pooled,
      tolerance = 1e-12
    )
  }
})

test_that("published Ahn T001 T011 and T012 score vectors are locked", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  scores <- list(
    T001 = c(0, 0, 1),
    T011 = c(0, 1, 1),
    T012 = c(0, 1, 2)
  )

  for (score in scores) {
    expected <- s7_ahn_reference_ncp(
      g_case, g_control, 1000, 1200, score
    )
    observed <- .cc_ahn_trend_ncp(
      g_case, g_control, 1000, 1200, score
    )
    expect_equal(observed, expected, tolerance = 1e-14)
  }
})

test_that("Ahn simulation conclusions are metadata, not computed PAWEH claims", {
  expect_identical(
    s7_ahn_validation_metadata$simulation_replicates_per_setting,
    10000L
  )
  expect_equal(
    s7_ahn_validation_metadata$
      reported_max_difference_percentage_points_approximately,
    0.3
  )
  expect_equal(
    s7_ahn_validation_metadata$reported_regression_r_squared_exceeded,
    0.999
  )
  expect_false(s7_ahn_validation_metadata$paweh_runs_simulation)
  expect_false(
    s7_ahn_validation_metadata$paweh_claims_full_ahn_table_reproduction
  )
  expect_match(
    s7_ahn_validation_metadata$reason_full_tables_are_out_of_scope,
    "p2.*LD"
  )
})

test_that("PAWEH zero-error sequencing transition fixtures are frozen", {
  # PAWEH sequencing bridge fixture; not a published Ahn numerical result.
  for (coverage in c(1, 2, 4, 10)) {
    tail <- 2^(-coverage)
    expected <- rbind(
      c(1, 0, 0),
      c(tail, 1 - 2 * tail, tail),
      c(0, 0, 1)
    )
    dimnames(expected) <- list(
      c("true_0", "true_1", "true_2"),
      c("called_0", "called_1", "called_2")
    )
    expect_equal(
      ngs_genotype_error_matrix(coverage, 0),
      expected,
      tolerance = 1e-15
    )
  }
})

test_that("PAWEH S4 direct and finite-depth NCP fixtures are frozen", {
  # PAWEH fixtures with complete parameter provenance; not Ahn table values.
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)

  direct <- .cc_ahn_trend_ncp(
    g_case, g_control, 1000, 1200, c(0, 1, 2)
  )
  finite <- .cc_ngs_ahn_ncp(
    g_case, g_control, 1000, 1200,
    coverage = 4, seq_error = 0, scores = c(0, 1, 2)
  )$lambda

  expect_equal(direct, 25.489186405767256, tolerance = 1e-14)
  expect_equal(finite, 23.240345866426704, tolerance = 1e-14)
})

test_that("PAWEH S5 high-depth public power fixture is frozen", {
  inputs <- list(
    N_case = 1000, alpha = 0.05,
    prev = 0.05, pd = 0.30, R2 = 1.8,
    coverage = 500, seq_error = 0.01,
    MOI = "M", k = 1.2, verbose = FALSE
  )
  ngs <- do.call(cc_ngs_power, inputs)
  canonical <- cc_power(
    N_case = inputs$N_case,
    alpha = inputs$alpha,
    input_mode = "model_based",
    prev = inputs$prev,
    pd = inputs$pd,
    R2 = inputs$R2,
    MOI = inputs$MOI,
    k = inputs$k,
    w = c(0, 1, 2),
    verbose = FALSE
  )

  expect_equal(
    canonical$tests$trend$power,
    0.99780646057077005,
    tolerance = 1e-14
  )
  expect_equal(ngs$power, 0.99780646057076983, tolerance = 1e-14)
  expect_equal(ngs$power, canonical$tests$trend$power, tolerance = 1e-12)
})

test_that("PAWEH S6 high-depth and finite-depth MSSN fixtures are frozen", {
  inputs <- list(
    power = 0.80, alpha = 0.05,
    prev = 0.05, pd = 0.30, R2 = 1.8,
    coverage = 500, seq_error = 0.01,
    MOI = "M", k = 1, verbose = FALSE
  )
  high <- do.call(cc_ngs_mssn, inputs)
  inputs$coverage <- 4
  inputs$seq_error <- 0
  finite <- do.call(cc_ngs_mssn, inputs)
  canonical <- cc_mssn(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05,
    pd = 0.30,
    R2 = 1.8,
    MOI = "M",
    k = 1,
    w = c(0, 1, 2),
    verbose = FALSE
  )

  expect_equal(canonical$tests$trend$MSSN_case, 373)
  expect_equal(high$MSSN_case, 373)
  expect_equal(finite$MSSN_case, 419)
  expect_gte(high$achieved_power + 1e-12, high$power_target)
  expect_gte(finite$achieved_power + 1e-12, finite$power_target)
  expect_lt(
    s7_integer_design_power(high, high$MSSN_case - 1)$power,
    high$power_target
  )
  expect_lt(
    s7_integer_design_power(finite, finite$MSSN_case - 1)$power,
    finite$power_target
  )
})

test_that("PAWEH analytic power-MSSN round trips hold on a modest grid", {
  # Deterministic coverage of every requested representative value without a
  # full Cartesian product. These are PAWEH regressions, not Ahn table values.
  designs <- list(
    list(MOI = "M", coverage = 2, seq_error = 0, k = 0.5,
         alpha = 0.05, power = 0.80, pd = 0.30, R2 = 1.8, prev = 0.05),
    list(MOI = "D", coverage = 4, seq_error = 0.001, k = 1,
         alpha = 5e-8, power = 0.90, pd = 0.20, R2 = 2.0, prev = 0.02),
    list(MOI = "Rec", coverage = 20, seq_error = 0.01, k = 2,
         alpha = 0.05, power = 0.90, pd = 0.40, R2 = 2.5, prev = 0.10),
    list(MOI = "M", coverage = 100, seq_error = 0.001, k = 2,
         alpha = 5e-8, power = 0.80, pd = 0.15, R2 = 2.2, prev = 0.01),
    list(MOI = "D", coverage = 20, seq_error = 0, k = 0.5,
         alpha = 0.05, power = 0.90, pd = 0.35, R2 = 1.7, prev = 0.08),
    list(MOI = "Rec", coverage = 4, seq_error = 0.01, k = 1,
         alpha = 5e-8, power = 0.80, pd = 0.45, R2 = 2.8, prev = 0.05)
  )

  for (design in designs) {
    result <- do.call(
      cc_ngs_mssn,
      c(design, list(verbose = FALSE))
    )
    achieved <- s7_integer_design_power(result, result$MSSN_case)

    expect_equal(result$MSSN_ctrl, achieved$N_control)
    expect_equal(result$achieved_lambda, achieved$lambda, tolerance = 1e-11)
    expect_equal(result$achieved_power, achieved$power, tolerance = 1e-14)
    expect_gte(achieved$power + 1e-12, design$power)
    if (result$MSSN_case > 1) {
      previous <- s7_integer_design_power(result, result$MSSN_case - 1)
      expect_lt(previous$power, design$power)
    }

    expect_equal(
      unname(rowSums(result$transition_matrix)),
      rep(1, 3),
      tolerance = 1e-14
    )
    for (g in result$freqs) {
      expect_true(all(is.finite(g)))
      expect_true(all(g >= 0 & g <= 1))
      expect_equal(sum(g), 1, tolerance = 1e-14)
    }
  }
})

test_that("PAWEH high-depth called frequencies approach true frequencies", {
  E <- ngs_genotype_error_matrix(coverage = 500, seq_error = 0.01)
  g_true <- c(0.70, 0.25, 0.05)
  g_called <- as.numeric(t(E) %*% g_true)

  expect_equal(unname(E), diag(3), tolerance = 1e-12)
  expect_equal(g_called, g_true, tolerance = 1e-12)
})

test_that("core Ahn NCP invariances remain frozen", {
  g_case <- c(0.45, 0.40, 0.15)
  g_control <- c(0.55, 0.35, 0.10)
  baseline <- .cc_ahn_trend_ncp(
    g_case, g_control, 1000, 1200, c(0, 1, 2)
  )

  expect_equal(
    .cc_ahn_trend_ncp(g_case, g_control, 2000, 2400, c(0, 1, 2)),
    2 * baseline,
    tolerance = 1e-12
  )
  expect_equal(
    .cc_ahn_trend_ncp(g_case, g_control, 1000, 1200, c(1, 2, 3)),
    baseline,
    tolerance = 1e-12
  )
  expect_equal(
    .cc_ahn_trend_ncp(g_control, g_case, 1200, 1000, c(0, 1, 2)),
    baseline,
    tolerance = 1e-14
  )
})
