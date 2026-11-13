cc_ngs_scenario_args <- function(...) {
  utils::modifyList(list(
    alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
    coverage = 10, seq_error = 0.01, MOI = "M", k = 1.2,
    verbose = FALSE
  ), list(...))
}

cc_ngs_all_modifiers <- function(...) {
  cc_ngs_scenario_args(
    pheno_misclass = TRUE, theta = 0.05, phi = 0.01,
    locus_het = TRUE, pi = 0.7,
    geno_misclass = "3p", e01 = 0.03, e02 = 0.01, e03 = 0.004,
    ...
  )
}

test_that("CC-NGS scenario activation is explicit and independent", {
  expected <- list(
    list(list(), "sequencing_only"),
    list(list(pheno_misclass = TRUE, theta = 0.05, phi = 0.01),
         c("sequencing_only", "phenotype_misclassification")),
    list(list(locus_het = TRUE, pi = 0.7),
         c("sequencing_only", "heterogeneity")),
    list(list(geno_misclass = "1p", e = 0.02),
         c("sequencing_only", "genotype_misclassification")),
    list(list(pheno_misclass = TRUE, theta = 0.05, phi = 0.01,
              locus_het = TRUE, pi = 0.7),
         c("sequencing_only", "phenotype_misclassification", "heterogeneity")),
    list(list(pheno_misclass = TRUE, theta = 0.05, phi = 0.01,
              geno_misclass = "1p", e = 0.02),
         c("sequencing_only", "phenotype_misclassification",
           "genotype_misclassification")),
    list(list(locus_het = TRUE, pi = 0.7,
              geno_misclass = "1p", e = 0.02),
         c("sequencing_only", "heterogeneity", "genotype_misclassification")),
    list(cc_ngs_all_modifiers(),
         c("sequencing_only", "phenotype_misclassification", "heterogeneity",
           "genotype_misclassification"))
  )
  for (entry in expected) {
    extra <- entry[[1L]]
    args <- if (length(extra) && !is.null(extra$alpha)) extra else
      utils::modifyList(cc_ngs_scenario_args(), extra)
    power <- do.call(cc_ngs_power, c(list(N_case = 800), args))
    mssn <- do.call(cc_ngs_mssn, c(list(power = 0.8), args))
    expect_identical(names(power$scenarios), entry[[2L]])
    expect_identical(names(mssn$scenarios), entry[[2L]])
  }
})

test_that("each modifier branch equals its standalone call", {
  all_args <- cc_ngs_all_modifiers(
    case_coverage = 4, ctrl_coverage = 15,
    case_seq_error = 0.04, ctrl_seq_error = 0.005
  )
  combined_power <- do.call(cc_ngs_power, c(list(N_case = 900), all_args))
  combined_mssn <- do.call(cc_ngs_mssn, c(list(power = 0.8), all_args))
  standalone <- list(
    phenotype_misclassification = cc_ngs_scenario_args(
      pheno_misclass = TRUE, theta = 0.05, phi = 0.01),
    heterogeneity = cc_ngs_scenario_args(locus_het = TRUE, pi = 0.7),
    genotype_misclassification = cc_ngs_scenario_args(
      geno_misclass = "3p", e01 = 0.03, e02 = 0.01, e03 = 0.004)
  )
  sequencing <- list(case_coverage = 4, ctrl_coverage = 15,
                     case_seq_error = 0.04, ctrl_seq_error = 0.005)
  for (name in names(standalone)) {
    args <- utils::modifyList(standalone[[name]], sequencing)
    p <- do.call(cc_ngs_power, c(list(N_case = 900), args))
    m <- do.call(cc_ngs_mssn, c(list(power = 0.8), args))
    for (field in c("freqs", "sequencing", "transition_matrix", "lambda",
                    "power", "D", "Q")) {
      expect_identical(combined_power$scenarios[[name]][[field]],
                       p$scenarios[[name]][[field]])
    }
    for (field in c("freqs", "sequencing", "transition_matrix", "D", "Q",
                    "N_case_continuous", "MSSN_case", "MSSN_ctrl",
                    "MSSN_total", "achieved_lambda", "achieved_power")) {
      expect_identical(combined_mssn$scenarios[[name]][[field]],
                       m$scenarios[[name]][[field]])
    }
  }
})

test_that("all modifier scenarios branch from common biological baseline", {
  out <- do.call(cc_ngs_power, c(list(N_case = 900), cc_ngs_all_modifiers(
    case_coverage = 4, ctrl_coverage = 15,
    case_seq_error = 0.04, ctrl_seq_error = 0.005
  )))
  base <- out$scenarios$sequencing_only
  phenotype <- out$scenarios$phenotype_misclassification
  heterogeneity <- out$scenarios$heterogeneity
  genotype <- out$scenarios$genotype_misclassification

  for (scenario in out$scenarios) {
    expect_identical(scenario$freqs$case_true_pre_heterogeneity,
                     base$freqs$case_true_pre_heterogeneity)
    expect_identical(scenario$freqs$control_true_pre_heterogeneity,
                     base$freqs$control_true_pre_heterogeneity)
    expect_identical(scenario$sequencing, base$sequencing)
  }
  expect_identical(genotype$freqs$case_preseq,
                   base$freqs$case_true_pre_heterogeneity)
  expect_identical(genotype$freqs$control_preseq,
                   base$freqs$control_true_pre_heterogeneity)
  expect_identical(genotype$freqs$case_post_sequencing,
                   base$freqs$case_post_sequencing)
  expect_identical(genotype$freqs$control_post_sequencing,
                   base$freqs$control_post_sequencing)
  expect_false(phenotype$errors$genotype_misclass$enabled)
  expect_false(phenotype$locus_het$enabled)
  expect_true(heterogeneity$locus_het$enabled)
  expect_false(heterogeneity$errors$phenotype_misclass$enabled)
  expect_false(genotype$locus_het$enabled)
  expect_false(genotype$errors$phenotype_misclass$enabled)
})

test_that("identity modifiers remain visible and equal sequencing only", {
  calls <- list(
    list(name = "phenotype_misclassification",
         args = list(pheno_misclass = TRUE, theta = 0, phi = 0)),
    list(name = "heterogeneity",
         args = list(locus_het = TRUE, pi = 1)),
    list(name = "genotype_misclassification",
         args = list(geno_misclass = "1p", e = 0))
  )
  for (call in calls) {
    args <- utils::modifyList(cc_ngs_scenario_args(), call$args)
    out <- do.call(cc_ngs_power, c(list(N_case = 700), args))
    expect_true(call$name %in% names(out$scenarios))
    expect_equal(out$scenarios[[call$name]]$lambda,
                 out$scenarios$sequencing_only$lambda, tolerance = 1e-13)
    expect_equal(out$scenarios[[call$name]]$power,
                 out$scenarios$sequencing_only$power, tolerance = 1e-14)
    expect_equal(out$scenarios[[call$name]]$freqs$case_final,
                 out$scenarios$sequencing_only$freqs$case_final,
                 tolerance = 1e-15)
  }
})

test_that("legacy fields use deterministic compatibility scenarios", {
  zero <- do.call(cc_ngs_power, c(list(N_case = 700), cc_ngs_scenario_args()))
  one <- do.call(cc_ngs_power, c(list(N_case = 700), cc_ngs_scenario_args(
    locus_het = TRUE, pi = 0.7
  )))
  many <- do.call(cc_ngs_power, c(list(N_case = 700), cc_ngs_all_modifiers()))
  expect_identical(zero$compatibility_scenario, "sequencing_only")
  expect_identical(one$compatibility_scenario, "heterogeneity")
  expect_identical(many$compatibility_scenario, "sequencing_only")
  for (out in list(zero, one, many)) {
    selected <- out$scenarios[[out$compatibility_scenario]]
    expect_identical(out$lambda, selected$lambda)
    expect_identical(out$power, selected$power)
    expect_identical(out$freqs, selected$freqs)
    expect_identical(out$sequencing, selected$sequencing)
  }
})

test_that("scenario D and Q reconstruct power and MSSN calculations", {
  power <- do.call(cc_ngs_power, c(list(N_case = 650), cc_ngs_all_modifiers()))
  mssn <- do.call(cc_ngs_mssn, c(list(power = 0.85), cc_ngs_all_modifiers()))
  for (name in names(power$scenarios)) {
    p <- power$scenarios[[name]]
    expect_equal(p$lambda, 650 * power$k * p$D^2 / p$Q,
                 tolerance = 1e-13)
    m <- mssn$scenarios[[name]]
    expect_equal(m$N_case_continuous,
                 m$lambda_target * m$Q / (mssn$k * m$D^2),
                 tolerance = 1e-12)
    achieved <- .cc_ahn_trend_ncp(
      m$freqs$case_final, m$freqs$control_final,
      m$MSSN_case, m$MSSN_ctrl, mssn$scores
    )
    expect_equal(m$achieved_lambda, achieved, tolerance = 1e-13)
    expect_gte(m$achieved_power + 1e-12, mssn$power_target)
  }
})

test_that("multi-scenario console output names every independent result", {
  out <- do.call(cc_ngs_power, c(list(N_case = 700), cc_ngs_all_modifiers()))
  text <- capture.output(.paweh_print_cc_ngs_power(out), type = "message")
  for (label in c("Sequencing only", "Phenotype misclassification",
                  "Locus heterogeneity", "Genotype misclassification")) {
    expect_true(any(grepl(label, text, fixed = TRUE)))
  }
})
