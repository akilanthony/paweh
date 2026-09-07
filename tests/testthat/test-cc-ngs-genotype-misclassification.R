cc_ngs_genotype_args <- function(...) {
  modifyList(list(
    N_case = 800, alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
    coverage = 12, seq_error = 0.02, MOI = "M", k = 1.25,
    verbose = FALSE
  ), list(...))
}

test_that("CC-NGS applies each nondifferential genotype model after sequencing", {
  settings <- list(
    `1p` = list(
      params = list(e = 0.03),
      M = matrix(c(.94, .03, .03, .03, .94, .03, .03, .03, .94),
                 nrow = 3, byrow = TRUE)
    ),
    `2p` = list(
      params = list(e1 = 0.04, e2 = 0.02),
      M = matrix(c(.96, .04, 0, .02, .96, .02, 0, .04, .96),
                 nrow = 3, byrow = TRUE)
    ),
    `3p` = list(
      params = list(e01 = 0.04, e02 = 0.02, e03 = 0.01),
      M = matrix(c(.95, .04, .01, .02, .96, .02, .01, .04, .95),
                 nrow = 3, byrow = TRUE)
    )
  )

  for (model in names(settings)) {
    out <- do.call(cc_ngs_power, c(
      cc_ngs_genotype_args(geno_misclass = model), settings[[model]]$params
    ))
    M <- settings[[model]]$M
    expected_case <- as.numeric(t(M) %*% out$freqs$case_post_sequencing)
    expected_ctrl <- as.numeric(t(M) %*% out$freqs$control_post_sequencing)

    expect_equal(out$errors$genotype_misclass$M_case, M,
                 tolerance = 1e-15, ignore_attr = TRUE)
    expect_equal(out$freqs$case_final, expected_case, tolerance = 1e-15)
    expect_equal(out$freqs$control_final, expected_ctrl, tolerance = 1e-15)
    expect_identical(out$freqs$case_called, out$freqs$case_final)
    expect_identical(out$freqs$control_called, out$freqs$control_final)
    expect_identical(
      out$errors$genotype_misclass$M_case,
      out$errors$genotype_misclass$M_ctrl
    )
  }
})

test_that("differential 3p supports explicit and multiplier parameterizations", {
  explicit <- cc_ngs_power(
    N_case = 800, alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
    coverage = 10, seq_error = 0.01, verbose = FALSE,
    geno_misclass = "diff3p",
    case_e01 = 0.04, case_e02 = 0.02, case_e03 = 0.01,
    ctrl_e01 = 0.02, ctrl_e02 = 0.01, ctrl_e03 = 0.005
  )
  from_case <- cc_ngs_power(
    N_case = 800, alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
    coverage = 10, seq_error = 0.01, verbose = FALSE,
    geno_misclass = "diff3p", diff_source = "case",
    diff_multiplier = 0.5,
    case_e01 = 0.04, case_e02 = 0.02, case_e03 = 0.01
  )
  M_case <- matrix(c(.95, .04, .01, .02, .96, .02, .01, .04, .95),
                   nrow = 3, byrow = TRUE)
  M_ctrl <- matrix(c(.975, .02, .005, .01, .98, .01, .005, .02, .975),
                   nrow = 3, byrow = TRUE)

  expect_equal(explicit$errors$genotype_misclass$M_case, M_case,
               tolerance = 1e-15, ignore_attr = TRUE)
  expect_equal(explicit$errors$genotype_misclass$M_ctrl, M_ctrl,
               tolerance = 1e-15, ignore_attr = TRUE)
  expect_equal(explicit$freqs$case_final,
               as.numeric(t(M_case) %*% explicit$freqs$case_post_sequencing),
               tolerance = 1e-15)
  expect_equal(explicit$freqs$control_final,
               as.numeric(t(M_ctrl) %*%
                            explicit$freqs$control_post_sequencing),
               tolerance = 1e-15)
  expect_equal(from_case$freqs$case_final, explicit$freqs$case_final,
               tolerance = 1e-15)
  expect_equal(from_case$freqs$control_final, explicit$freqs$control_final,
               tolerance = 1e-15)
  expect_equal(from_case$lambda, explicit$lambda, tolerance = 1e-15)
  expect_identical(
    from_case$errors$genotype_misclass$model,
    "diff3p_homhet_homhom"
  )
  expect_identical(from_case$errors$genotype_misclass$diff_source, "case")
  expect_false(identical(
    from_case$errors$genotype_misclass$M_case,
    from_case$errors$genotype_misclass$M_ctrl
  ))
})

test_that("CC-NGS genotype stage follows sequencing and feeds the Ahn trend", {
  out <- cc_ngs_power(
    N_case = 700, alpha = 0.01, prev = 0.03, pd = 0.25, R2 = 2,
    case_coverage = 8, ctrl_coverage = 20,
    case_seq_error = 0.04, ctrl_seq_error = 0.01,
    pheno_misclass = TRUE, theta = 0.03, phi = 0.01,
    locus_het = TRUE, pi = 0.8,
    geno_misclass = "diff3p",
    case_e01 = 0.03, case_e02 = 0.01, case_e03 = 0.004,
    ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002,
    k = 1.3, verbose = FALSE
  )
  expected_lambda <- .cc_ahn_trend_ncp(
    out$freqs$case_final, out$freqs$control_final,
    out$N_case, out$N_ctrl, out$scores
  )

  expect_equal(out$lambda, expected_lambda, tolerance = 1e-14)
  expect_equal(
    out$freqs$case_post_sequencing,
    as.numeric(t(out$sequencing$case_transition_matrix) %*%
                 out$freqs$case_preseq),
    tolerance = 1e-15
  )
  expect_equal(
    out$freqs$case_final,
    as.numeric(t(out$errors$genotype_misclass$M_case) %*%
                 out$freqs$case_post_sequencing),
    tolerance = 1e-15
  )
})

test_that("CC-NGS MSSN uses the final genotype frequencies", {
  out <- cc_ngs_mssn(
    power = 0.8, alpha = 0.05, prev = 0.05, pd = 0.3, R2 = 1.8,
    coverage = 15, seq_error = 0.02, k = 1.2, verbose = FALSE,
    geno_misclass = "2p", e1 = 0.04, e2 = 0.02
  )
  lambda <- .cc_ahn_trend_ncp(
    out$freqs$case_final, out$freqs$control_final,
    out$MSSN_case, out$MSSN_ctrl, out$scores
  )

  expect_equal(out$achieved_lambda, lambda, tolerance = 1e-14)
  expect_identical(out$freqs$case_called, out$freqs$case_final)
  expect_true(out$errors$genotype_misclass$enabled)
})

test_that("inactive genotype defaults preserve prior CC-NGS values", {
  baseline <- do.call(cc_ngs_power, cc_ngs_genotype_args())
  explicit <- do.call(cc_ngs_power, cc_ngs_genotype_args(
    geno_misclass = "none"
  ))

  expect_identical(baseline$lambda, explicit$lambda)
  expect_identical(baseline$power, explicit$power)
  expect_identical(
    baseline$freqs$case_post_sequencing, baseline$freqs$case_called
  )
  expect_identical(
    baseline$freqs$control_post_sequencing, baseline$freqs$control_called
  )
  expect_false(baseline$errors$genotype_misclass$enabled)

  for (zero in list(
    list(geno_misclass = "1p", e = 0),
    list(geno_misclass = "2p", e1 = 0, e2 = 0),
    list(geno_misclass = "3p", e01 = 0, e02 = 0, e03 = 0)
  )) {
    observed <- do.call(cc_ngs_power, c(cc_ngs_genotype_args(), zero))
    expect_identical(observed$lambda, baseline$lambda)
    expect_identical(observed$power, baseline$power)
    expect_equal(observed$freqs$case_final, baseline$freqs$case_called,
                 tolerance = 1e-15)
    expect_equal(observed$freqs$control_final, baseline$freqs$control_called,
                 tolerance = 1e-15)
  }
})

test_that("CC-NGS and ordinary case-control share genotype model semantics", {
  ngs <- cc_ngs_power(
    N_case = 600, alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
    coverage = 12, seq_error = 0.02, k = 1.4, verbose = FALSE,
    geno_misclass = "3p", e01 = 0.04, e02 = 0.02, e03 = 0.01
  )
  ordinary <- cc_power(
    N_case = 600, alpha = 0.05, input_mode = "model_free",
    g1 = ngs$freqs$case_post_sequencing,
    g0 = ngs$freqs$control_post_sequencing,
    k = 1.4, w = ngs$scores, verbose = FALSE,
    geno_misclass = "3p", e01 = 0.04, e02 = 0.02, e03 = 0.01
  )

  expect_identical(
    ngs$errors$genotype_misclass$M,
    ordinary$errors$genotype_misclass$M
  )
  expect_equal(ngs$freqs$case_final, ordinary$freqs$g_obs_case,
               tolerance = 1e-15)
  expect_equal(ngs$freqs$control_final, ordinary$freqs$g_obs_ctrl,
               tolerance = 1e-15)
  expect_equal(ngs$lambda, ordinary$tests$trend$lambda, tolerance = 1e-12)
})

test_that("genotype model validation errors propagate through both CC-NGS APIs", {
  for (fun in list(cc_ngs_power, cc_ngs_mssn)) {
    first <- if (identical(fun, cc_ngs_power)) {
      list(N_case = 500)
    } else {
      list(power = 0.8)
    }
    common <- list(
      alpha = 0.05, prev = 0.05, pd = 0.3, R2 = 1.8,
      coverage = 10, seq_error = 0.01, verbose = FALSE
    )
    expect_error(
      do.call(fun, c(first, common, list(geno_misclass = "1p", e = 0.6))),
      "e must be a single number in [0, 0.5].", fixed = TRUE
    )
    expect_error(
      do.call(fun, c(first, common, list(
        geno_misclass = "diff3p", diff_source = "case",
        case_e01 = 0.8, case_e03 = 0.1, diff_multiplier = 2
      ))),
      "Scaled e01 is outside [0,1].", fixed = TRUE
    )
    expect_error(
      do.call(fun, c(first, common, list(geno_misclass = "2p", e1 = -0.1))),
      "e1 must be a single number in [0,1].", fixed = TRUE
    )
    expect_error(
      do.call(fun, c(first, common, list(geno_misclass = "2p", e2 = 0.6))),
      "e2 must be a single number in [0,0.5].", fixed = TRUE
    )
    expect_error(
      do.call(fun, c(first, common, list(
        geno_misclass = "3p", e01 = 0.8, e03 = 0.3
      ))),
      "Need e01 + e03 <= 1.", fixed = TRUE
    )
    expect_error(
      do.call(fun, c(first, common, list(diff_source = "invalid"))),
      "should be one of"
    )
    expect_error(
      do.call(fun, c(first, common, list(diff_multiplier = -1))),
      "diff_multiplier must be a single nonnegative number.", fixed = TRUE
    )
  }
})

test_that("appended genotype arguments preserve legacy positional calls", {
  positional <- cc_ngs_power(
    500, 0.05, 0.05, 0.3, 1.8, 10, 0.01, "M", 1.2, FALSE,
    TRUE, 0.9, 8, 12, 0.02, 0.01, TRUE, 0.02, 0.01
  )
  named <- cc_ngs_power(
    N_case = 500, alpha = 0.05, prev = 0.05, pd = 0.3, R2 = 1.8,
    coverage = 10, seq_error = 0.01, MOI = "M", k = 1.2,
    verbose = FALSE, locus_het = TRUE, pi = 0.9,
    case_coverage = 8, ctrl_coverage = 12,
    case_seq_error = 0.02, ctrl_seq_error = 0.01,
    pheno_misclass = TRUE, theta = 0.02, phi = 0.01
  )

  expect_identical(positional, named)
})

test_that("sequential and total-matrix compositions agree for both groups", {
  out <- cc_ngs_power(
    N_case = 900, alpha = 0.05, prev = 0.04, pd = 0.28, R2 = 1.9,
    case_coverage = 6, ctrl_coverage = 18,
    case_seq_error = 0.05, ctrl_seq_error = 0.005,
    geno_misclass = "diff3p",
    case_e01 = 0.04, case_e02 = 0.02, case_e03 = 0.006,
    ctrl_e01 = 0.015, ctrl_e02 = 0.008, ctrl_e03 = 0.002,
    verbose = FALSE
  )
  case_total <- out$sequencing$case_transition_matrix %*%
    out$errors$genotype_misclass$M_case
  ctrl_total <- out$sequencing$ctrl_transition_matrix %*%
    out$errors$genotype_misclass$M_ctrl

  expect_equal(out$freqs$case_final,
               as.numeric(t(case_total) %*% out$freqs$case_preseq),
               tolerance = 1e-15)
  expect_equal(out$freqs$control_final,
               as.numeric(t(ctrl_total) %*% out$freqs$control_preseq),
               tolerance = 1e-15)
  expect_false(identical(out$sequencing$case_transition_matrix,
                         out$sequencing$ctrl_transition_matrix))
  expect_false(identical(out$errors$genotype_misclass$M_case,
                         out$errors$genotype_misclass$M_ctrl))

  for (M in list(
    out$sequencing$case_transition_matrix,
    out$sequencing$ctrl_transition_matrix,
    out$errors$genotype_misclass$M_case,
    out$errors$genotype_misclass$M_ctrl,
    case_total, ctrl_total
  )) {
    expect_true(all(M >= 0 & M <= 1))
    expect_equal(unname(rowSums(M)), rep(1, 3), tolerance = 1e-14)
  }
})

test_that("effectively perfect sequencing retains the CC-NGS sequential kernel", {
  common <- list(
    prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M", k = 1.2,
    locus_het = TRUE, pi = 0.8,
    pheno_misclass = TRUE, theta = 0.02, phi = 0.01,
    geno_misclass = "3p", e01 = 0.03, e02 = 0.015, e03 = 0.004,
    verbose = FALSE
  )
  ngs <- do.call(cc_ngs_power, c(
    list(N_case = 1000, alpha = 0.05, coverage = 100, seq_error = 0),
    common
  ))
  expected_case <- .cc_apply_genotype_misclass(
    ngs$freqs$case_post_sequencing,
    ngs$errors$genotype_misclass$M_case
  )
  expected_ctrl <- .cc_apply_genotype_misclass(
    ngs$freqs$control_post_sequencing,
    ngs$errors$genotype_misclass$M_ctrl
  )
  expected_lambda <- .cc_ahn_trend_ncp(
    expected_case, expected_ctrl, 1000, 1200, c(0, 1, 2)
  )
  expected_power <- pchisq(
    qchisq(0.95, df = 1), df = 1, ncp = expected_lambda,
    lower.tail = FALSE
  )

  expect_equal(ngs$freqs$case_final, expected_case, tolerance = 1e-15)
  expect_equal(ngs$freqs$control_final, expected_ctrl, tolerance = 1e-15)
  expect_equal(ngs$lambda, expected_lambda, tolerance = 1e-12)
  expect_equal(ngs$power, expected_power, tolerance = 1e-14)
})

test_that("MSSN round trips for 1p, 3p, and differential 3p", {
  settings <- list(
    list(geno_misclass = "1p", e = 0.02),
    list(
      geno_misclass = "3p", e01 = 0.03, e02 = 0.01, e03 = 0.004,
      pheno_misclass = TRUE, theta = 0.02, phi = 0.01
    ),
    list(
      geno_misclass = "diff3p",
      case_e01 = 0.03, case_e02 = 0.01, case_e03 = 0.004,
      ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002,
      case_coverage = 8, ctrl_coverage = 20,
      case_seq_error = 0.04, ctrl_seq_error = 0.01
    )
  )

  for (setting in settings) {
    common <- c(list(
      alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
      coverage = 12, seq_error = 0.02, k = 1, verbose = FALSE
    ), setting)
    mssn <- do.call(cc_ngs_mssn, c(list(power = 0.8), common))
    achieved <- do.call(cc_ngs_power, c(list(N_case = mssn$MSSN_case), common))
    expect_gte(achieved$power + 1e-12, 0.8)
    expect_equal(achieved$lambda, mssn$achieved_lambda, tolerance = 1e-12)

    if (mssn$MSSN_case > 1) {
      previous <- do.call(cc_ngs_power, c(
        list(N_case = mssn$MSSN_case - 1), common
      ))
      expect_lt(previous$power, 0.8)
    }
  }
})

test_that("active console output identifies genotype model and differential status", {
  out <- cc_ngs_power(
    N_case = 500, alpha = 0.05, prev = 0.05, pd = 0.3, R2 = 1.8,
    coverage = 10, seq_error = 0.01,
    geno_misclass = "diff3p",
    case_e01 = 0.02, ctrl_e01 = 0.01,
    verbose = FALSE
  )
  printed <- capture.output(.paweh_print_cc_ngs_power(out), type = "message")
  expect_true(any(grepl("Differential 3-parameter", printed, fixed = TRUE)))
  expect_true(any(grepl("Differential error source", printed, fixed = TRUE)))
})
