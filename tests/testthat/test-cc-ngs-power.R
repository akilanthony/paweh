cc_ngs_test_args <- function(...) {
  modifyList(
    list(
      N_case = 1000,
      alpha = 0.05,
      prev = 0.05,
      pd = 0.30,
      R2 = 1.8,
      coverage = 20,
      seq_error = 0.01,
      MOI = "M",
      k = 1.2,
      verbose = FALSE
    ),
    list(...)
  )
}

test_that("S5 disease model reproduces canonical cc_power frequencies", {
  settings <- list(
    list(pd = 0.10, R2 = 1.4, MOI = "M", prev = 0.01),
    list(pd = 0.30, R2 = 1.8, MOI = "D", prev = 0.05),
    list(pd = 0.45, R2 = 2.2, MOI = "Rec", prev = 0.10),
    list(pd = 0.20, R2 = 0.8, MOI = "M", prev = 0.02)
  )

  for (setting in settings) {
    model <- do.call(.cc_model_genotype_frequencies, setting)
    canonical <- cc_power(
      N_case = 500,
      alpha = 0.05,
      input_mode = "model_based",
      prev = setting$prev,
      pd = setting$pd,
      R2 = setting$R2,
      MOI = setting$MOI,
      locus_het = FALSE,
      pheno_misclass = FALSE,
      geno_misclass = "none",
      w = .cc_ngs_scores_from_moi(setting$MOI),
      verbose = FALSE
    )

    expect_equal(model$case, canonical$freqs$g_base_case, tolerance = 1e-15)
    expect_equal(model$control, canonical$freqs$g_base_ctrl, tolerance = 1e-15)
    expect_equal(model$penetrances, canonical$model_info$penetrances,
                 tolerance = 1e-15)
    expect_equal(model$R1, canonical$model_info$R1, tolerance = 1e-15)
  }
})

test_that("population frequencies are HWE probabilities", {
  for (pd in c(0.01, 0.10, 0.30, 0.75, 0.99)) {
    model <- .cc_model_genotype_frequencies(pd, 1.8, "M", 0.05)
    expected <- c((1 - pd)^2, 2 * pd * (1 - pd), pd^2)

    expect_equal(model$population, expected, tolerance = 1e-15)
    expect_equal(sum(model$population), 1, tolerance = 1e-15)
  }
})

test_that("penetrances reproduce the supplied prevalence", {
  settings <- list(
    list(pd = 0.10, R2 = 1.4, MOI = "M", prev = 0.01),
    list(pd = 0.30, R2 = 1.8, MOI = "D", prev = 0.05),
    list(pd = 0.45, R2 = 2.2, MOI = "Rec", prev = 0.10)
  )

  for (setting in settings) {
    model <- do.call(.cc_model_genotype_frequencies, setting)
    prevalence <- sum(model$population * model$penetrances)
    expect_equal(prevalence, setting$prev, tolerance = 1e-15)
  }
})

test_that("MOI maps to the specified sequencing trend scores", {
  expect_identical(.cc_ngs_scores_from_moi("M"), c(0, 1, 2))
  expect_identical(.cc_ngs_scores_from_moi("D"), c(0, 1, 1))
  expect_identical(.cc_ngs_scores_from_moi("Rec"), c(0, 0, 1))
})

test_that("public power matches an independent sequencing composition", {
  args <- cc_ngs_test_args(
    N_case = 900,
    alpha = 0.01,
    coverage = 12,
    seq_error = 0.04,
    MOI = "D",
    k = 1.5
  )
  observed <- do.call(cc_ngs_power, args)
  model <- .cc_model_genotype_frequencies(
    args$pd, args$R2, args$MOI, args$prev
  )
  scores <- c(0, 1, 1)
  E <- ngs_genotype_error_matrix(args$coverage, args$seq_error)
  case_called <- as.numeric(t(E) %*% model$case)
  control_called <- as.numeric(t(E) %*% model$control)
  weighted_counts <- args$N_case * case_called +
    (args$k * args$N_case) * control_called
  contrast <- sum(scores * (case_called - control_called))
  denominator <- sum(scores^2 * weighted_counts) -
    sum(scores * weighted_counts)^2 /
      (args$N_case + args$k * args$N_case)
  lambda <- args$N_case * (args$k * args$N_case) *
    contrast^2 / denominator
  critical <- qchisq(1 - args$alpha, df = 1)
  power <- pchisq(critical, df = 1, ncp = lambda, lower.tail = FALSE)

  expect_equal(observed$lambda, lambda, tolerance = 1e-12)
  expect_equal(observed$power, power, tolerance = 1e-14)
  expect_equal(observed$transition_matrix, E, tolerance = 1e-15)
  expect_equal(observed$freqs$case_called, case_called, tolerance = 1e-15)
  expect_equal(observed$freqs$control_called, control_called, tolerance = 1e-15)
})

test_that("high-depth sequencing power approaches canonical trend power", {
  for (MOI in c("M", "D", "Rec")) {
    args <- cc_ngs_test_args(MOI = MOI, coverage = 500, seq_error = 0.01)
    sequencing <- do.call(cc_ngs_power, args)
    canonical <- cc_power(
      N_case = args$N_case,
      alpha = args$alpha,
      input_mode = "model_based",
      prev = args$prev,
      pd = args$pd,
      R2 = args$R2,
      MOI = MOI,
      k = args$k,
      w = .cc_ngs_scores_from_moi(MOI),
      verbose = FALSE
    )

    expect_equal(
      sequencing$power,
      canonical$tests$trend$power,
      tolerance = 1e-10
    )
    expect_equal(
      sequencing$lambda,
      canonical$tests$trend$lambda,
      tolerance = 1e-10
    )
  }
})

test_that("zero-error finite-depth power retains heterozygote uncertainty", {
  args <- cc_ngs_test_args(coverage = 4, seq_error = 0, MOI = "M")
  sequencing <- do.call(cc_ngs_power, args)
  canonical <- cc_power(
    N_case = args$N_case,
    alpha = args$alpha,
    input_mode = "model_based",
    prev = args$prev,
    pd = args$pd,
    R2 = args$R2,
    MOI = args$MOI,
    k = args$k,
    w = c(0, 1, 2),
    verbose = FALSE
  )
  tail <- 2^(-args$coverage)

  expect_equal(
    unname(sequencing$transition_matrix["true_1", ]),
    c(tail, 1 - 2 * tail, tail),
    tolerance = 1e-15
  )
  expect_lt(sequencing$power, canonical$tests$trend$power)
  expect_false(isTRUE(all.equal(
    sequencing$power,
    canonical$tests$trend$power,
    tolerance = 1e-12
  )))
})

test_that("representative coverage values approach ordinary power", {
  coverages <- c(2, 4, 20, 100)
  powers <- vapply(
    coverages,
    function(coverage) {
      do.call(
        cc_ngs_power,
        cc_ngs_test_args(coverage = coverage, seq_error = 0.01)
      )$power
    },
    numeric(1)
  )
  canonical <- cc_power(
    N_case = 1000,
    alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05,
    pd = 0.30,
    R2 = 1.8,
    MOI = "M",
    k = 1.2,
    w = c(0, 1, 2),
    verbose = FALSE
  )$tests$trend$power

  expect_true(all(is.finite(powers)))
  expect_lt(abs(powers[[4]] - canonical), abs(powers[[1]] - canonical))
  expect_equal(powers[[4]], canonical, tolerance = 1e-10)
})

test_that("lambda scales with case and control sample sizes", {
  baseline <- do.call(cc_ngs_power, cc_ngs_test_args(N_case = 600))
  doubled <- do.call(cc_ngs_power, cc_ngs_test_args(N_case = 1200))

  expect_equal(doubled$lambda, 2 * baseline$lambda, tolerance = 1e-12)
  expect_equal(doubled$N_ctrl, 2 * baseline$N_ctrl)
})

test_that("power increases with less stringent alpha", {
  powers <- vapply(
    c(0.05, 0.001, 5e-8),
    function(alpha) {
      do.call(cc_ngs_power, cc_ngs_test_args(alpha = alpha))$power
    },
    numeric(1)
  )

  expect_gte(powers[[1]], powers[[2]])
  expect_gte(powers[[2]], powers[[3]])
})

test_that("cc_ngs_power has a transparent valid return structure", {
  visible <- withVisible(do.call(cc_ngs_power, cc_ngs_test_args()))
  result <- visible$value
  expected_fields <- c(
    "alpha", "N_case", "N_ctrl", "N_total", "k", "power", "lambda",
    "MOI", "scores", "coverage", "seq_error", "model_info", "freqs",
    "transition_matrix"
  )

  expect_false(visible$visible)
  expect_identical(class(result), "cc_ngs_power")
  expect_true(all(expected_fields %in% names(result)))
  expect_true(is.finite(result$power))
  expect_true(result$power >= 0 && result$power <= 1)
  expect_true(is.finite(result$lambda))
  expect_gte(result$lambda, 0)
  expect_equal(result$N_ctrl, result$k * result$N_case)
  expect_equal(result$N_total, result$N_case + result$N_ctrl)
  expect_identical(dim(result$transition_matrix), c(3L, 3L))
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
})

test_that("verbose output is concise and optional", {
  expect_output(
    invisible(do.call(cc_ngs_power, cc_ngs_test_args(verbose = FALSE))),
    NA
  )
  expect_output(
    invisible(do.call(cc_ngs_power, cc_ngs_test_args(verbose = TRUE))),
    "Case-control sequencing trend-test power"
  )
  output <- capture.output(
    invisible(do.call(cc_ngs_power, cc_ngs_test_args(verbose = TRUE)))
  )

  expect_true(any(grepl("Coverage:", output, fixed = TRUE)))
  expect_true(any(grepl("NCP:", output, fixed = TRUE)))
  expect_false(any(grepl("transition_matrix", output, fixed = TRUE)))
  expect_lte(length(output), 5)
})

test_that("cc_ngs_power validates public inputs", {
  invalid_values <- list(
    pd = list(0, 1, -0.1, NA_real_, NaN, Inf, c(0.2, 0.3)),
    R2 = list(0, -1, NA_real_, NaN, Inf, c(1, 2)),
    prev = list(0, 1, -0.1, NA_real_, NaN, Inf, c(0.1, 0.2)),
    N_case = list(0, -1, NA_real_, NaN, Inf, c(100, 200)),
    k = list(0, -1, NA_real_, NaN, Inf, c(1, 2)),
    alpha = list(0, 1, -0.1, NA_real_, NaN, Inf, c(0.01, 0.05)),
    coverage = list(0, -1, 1.5, NA_real_, NaN, Inf, c(10, 20)),
    seq_error = list(-0.1, 0.5, 1, NA_real_, NaN, Inf, c(0.01, 0.02))
  )

  for (argument in names(invalid_values)) {
    for (value in invalid_values[[argument]]) {
      args <- cc_ngs_test_args()
      args[[argument]] <- value
      expect_error(do.call(cc_ngs_power, args), argument)
    }
  }

  expect_error(
    do.call(cc_ngs_power, cc_ngs_test_args(MOI = "additive")),
    "arg"
  )
  expect_error(
    do.call(cc_ngs_power, cc_ngs_test_args(verbose = NA)),
    "verbose"
  )
})

test_that("CC-NGS power API excludes deferred sequencing features", {
  arguments <- names(formals(cc_ngs_power))
  excluded <- c(
    "target_power", "MSSN", "pheno_misclass", "theta", "phi",
    "case_seq_error", "ctrl_seq_error", "coverage_dist", "g1", "g0",
    "raw_reads", "trios"
  )

  expect_false(any(excluded %in% arguments))
})
