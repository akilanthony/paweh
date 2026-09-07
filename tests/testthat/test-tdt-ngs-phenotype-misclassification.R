tdt_ngs_pheno_power <- function(...) {
  args <- list(
    N = 1200, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, alpha = 0.05, verbose = FALSE,
    pheno_misclass = TRUE, prev = 0.05, pi01 = 0.08
  )
  do.call(tdt_ngs_power, utils::modifyList(args, list(...)))
}

tdt_ngs_pheno_mssn <- function(...) {
  args <- list(
    power = 0.80, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, alpha = 0.05, verbose = FALSE,
    pheno_misclass = TRUE, prev = 0.05, pi01 = 0.08
  )
  do.call(tdt_ngs_mssn, utils::modifyList(args, list(...)))
}

test_that("phenotype sensitivity creates the two exact scenarios", {
  inactive <- tdt_ngs_power(
    N = 1200, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, verbose = FALSE
  )
  active <- tdt_ngs_pheno_power()

  expect_named(inactive$scenarios, "sequencing_only")
  expect_identical(inactive$compatibility_scenario, "sequencing_only")
  expect_named(
    active$scenarios,
    c("sequencing_only", "phenotype_misclassification")
  )
  expect_identical(
    active$compatibility_scenario, "phenotype_misclassification"
  )
  expect_identical(active$lambda, active$scenarios[[2]]$lambda)
  expect_identical(active$power, active$scenarios[[2]]$power)
  expect_identical(active$scenarios[[1]]$lambda, inactive$lambda)
  expect_identical(active$information_matrix, inactive$information_matrix)
  expect_identical(
    active$scenarios[[2]]$model_info$raw_read_likelihood_modified, FALSE
  )
})

test_that("ordinary-TDT helpers exactly reconstruct the bridge factor", {
  out <- tdt_ngs_pheno_power()
  bridge <- out$scenarios$phenotype_misclassification$ordinary_tdt_bridge
  fixed <- list(
    pd = out$pd, prev = 0.05, R1 = out$R1, R2 = out$R1^2,
    delta_prime = 1, theta1 = out$pd, verbose = FALSE
  )
  t0 <- do.call(tdt_expected_transmission_probability,
                c(fixed, list(pi01 = 0)))$gT_star
  nt0 <- do.call(tdt_expected_nontransmission_probability,
                 c(fixed, list(pi01 = 0)))$gNT_star
  t1 <- do.call(tdt_expected_transmission_probability,
                c(fixed, list(pi01 = 0.08)))$gT_star
  nt1 <- do.call(tdt_expected_nontransmission_probability,
                 c(fixed, list(pi01 = 0.08)))$gNT_star
  c0 <- 2 * (t0 - nt0)^2 / (t0 + nt0)
  c1 <- 2 * (t1 - nt1)^2 / (t1 + nt1)

  expect_equal(bridge$gT_noerror, t0, tolerance = 1e-15)
  expect_equal(bridge$gNT_noerror, nt0, tolerance = 1e-15)
  expect_equal(bridge$gT_pheno, t1, tolerance = 1e-15)
  expect_equal(bridge$gNT_pheno, nt1, tolerance = 1e-15)
  expect_equal(bridge$ncp_per_trio_noerror, c0, tolerance = 1e-15)
  expect_equal(bridge$ncp_per_trio_pheno, c1, tolerance = 1e-15)
  expect_equal(bridge$attenuation_factor, c1 / c0, tolerance = 1e-15)
  expect_identical(bridge$R2, out$R1^2)
  expect_identical(bridge$delta_prime, 1)
  expect_identical(bridge$theta1, out$pd)
  expect_identical(bridge$method, "ordinary_TDT_NCP_attenuation_bridge")
  expect_identical(bridge$raw_read_likelihood_modified, FALSE)
})

test_that("phenotype NCP and power are independently reproducible", {
  out <- tdt_ngs_pheno_power()
  base <- out$scenarios$sequencing_only
  adjusted <- out$scenarios$phenotype_misclassification
  expected_lambda <- base$lambda * adjusted$attenuation_factor
  critical <- stats::qchisq(1 - out$alpha, df = 1)
  expected_power <- stats::pchisq(
    critical, df = 1, ncp = expected_lambda, lower.tail = FALSE
  )

  expect_equal(adjusted$base_lambda, base$lambda, tolerance = 1e-15)
  expect_equal(adjusted$lambda, expected_lambda, tolerance = 1e-13)
  expect_equal(adjusted$power, expected_power, tolerance = 1e-15)
  expect_lt(adjusted$lambda, base$lambda)
  expect_lt(adjusted$power, base$power)
})

test_that("MSSN applies phenotype attenuation before inversion and is minimal", {
  out <- tdt_ngs_pheno_mssn()
  base <- out$scenarios$sequencing_only
  adjusted <- out$scenarios$phenotype_misclassification
  expected_coefficient <- base$ncp_per_trio * adjusted$attenuation_factor
  expected_continuous <- adjusted$lambda_target / expected_coefficient

  expect_equal(adjusted$ncp_per_trio, expected_coefficient, tolerance = 1e-15)
  expect_equal(adjusted$N_trios_continuous, expected_continuous,
               tolerance = 1e-12)
  expect_identical(out$MSSN_trios, adjusted$MSSN_trios)
  expect_gte(adjusted$achieved_power, out$power_target)

  at_plan <- tdt_ngs_pheno_power(N = adjusted$MSSN_trios)
  expect_equal(at_plan$power, adjusted$achieved_power, tolerance = 1e-14)
  expect_gte(at_plan$power, out$power_target)
  if (adjusted$MSSN_trios > 1) {
    one_lower <- tdt_ngs_pheno_power(N = adjusted$MSSN_trios - 1)
    expect_lt(one_lower$power, out$power_target)
  }
})

test_that("bridge factor is invariant to sequencing design", {
  equal <- tdt_ngs_pheno_power()
  unequal_directional <- tdt_ngs_pheno_power(
    coverage = NULL, seq_error = NULL,
    father_coverage = 4, mother_coverage = 10, child_coverage = 18,
    epsilon0 = 0.004, epsilon1 = 0.025
  )
  e <- equal$scenarios$phenotype_misclassification
  u <- unequal_directional$scenarios$phenotype_misclassification

  expect_identical(e$attenuation_factor, u$attenuation_factor)
  expect_identical(e$ordinary_tdt_bridge, u$ordinary_tdt_bridge)
  expect_equal(u$lambda,
               u$base_lambda * u$attenuation_factor, tolerance = 1e-13)
  expect_identical(
    unequal_directional$sequencing$coverage,
    c(father = 4, mother = 10, child = 18)
  )
  expect_identical(unequal_directional$epsilon0, 0.004)
  expect_identical(unequal_directional$epsilon1, 0.025)
})

test_that("pi01 zero is an explicit identity scenario", {
  active <- tdt_ngs_pheno_power(pi01 = 0)
  base <- active$scenarios$sequencing_only
  adjusted <- active$scenarios$phenotype_misclassification
  mssn <- tdt_ngs_pheno_mssn(pi01 = 0)

  expect_identical(adjusted$attenuation_factor, 1)
  expect_identical(adjusted$lambda, base$lambda)
  expect_identical(adjusted$power, base$power)
  expect_identical(adjusted$ordinary_tdt_bridge$status, "identity_pi01_zero")
  expect_identical(
    mssn$scenarios$phenotype_misclassification$MSSN_trios,
    mssn$scenarios$sequencing_only$MSSN_trios
  )
  expect_identical(
    mssn$scenarios$phenotype_misclassification$attenuation_factor, 1
  )
})

test_that("R1 one avoids an indeterminate phenotype ratio", {
  out <- expect_no_warning(tdt_ngs_pheno_power(R1 = 1))
  adjusted <- out$scenarios$phenotype_misclassification
  expect_identical(adjusted$attenuation_factor, 1)
  expect_identical(
    adjusted$ordinary_tdt_bridge$status, "zero_baseline_effect_identity"
  )
  expect_identical(out$scenarios$sequencing_only$lambda, 0)
  expect_identical(adjusted$lambda, 0)
  numeric_bridge <- adjusted$ordinary_tdt_bridge[vapply(
    adjusted$ordinary_tdt_bridge, is.numeric, logical(1)
  )]
  expect_true(all(is.finite(unlist(numeric_bridge))))

  expect_error(
    tdt_ngs_pheno_mssn(R1 = 1),
    "No finite MSSN exists because R1 = 1"
  )
  supported <- expect_no_warning(tdt_ngs_pheno_mssn(
    R1 = 1, power = 0.01, alpha = 0.05
  ))
  expect_identical(supported$MSSN_trios, 1)
  expect_identical(
    supported$scenarios$phenotype_misclassification$attenuation_factor, 1
  )
})

test_that("phenotype bridge inputs are validated", {
  base_power <- list(
    N = 100, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, verbose = FALSE
  )
  invalid <- list(
    list(pheno_misclass = 1),
    list(pheno_misclass = NA),
    list(pheno_misclass = c(TRUE, FALSE)),
    list(pheno_misclass = TRUE),
    list(pheno_misclass = TRUE, prev = 0),
    list(pheno_misclass = TRUE, prev = 1),
    list(pheno_misclass = TRUE, prev = NA_real_),
    list(pheno_misclass = TRUE, prev = c(0.01, 0.02)),
    list(pheno_misclass = TRUE, prev = 0.05, pi01 = -0.01),
    list(pheno_misclass = TRUE, prev = 0.05, pi01 = 1),
    list(pheno_misclass = TRUE, prev = 0.05, pi01 = NA_real_),
    list(pheno_misclass = FALSE, pi01 = 0.01)
  )
  for (values in invalid) {
    expect_error(do.call(tdt_ngs_power, utils::modifyList(base_power, values)))
  }
  with_prev <- do.call(
    tdt_ngs_power, utils::modifyList(base_power, list(prev = 0.05))
  )
  without_prev <- do.call(tdt_ngs_power, base_power)
  expect_identical(with_prev$lambda, without_prev$lambda)
  expect_identical(with_prev$power, without_prev$power)
})

test_that("high-quality phenotype TDT1-NGS approaches ordinary phenotype TDT", {
  ngs <- tdt_ngs_power(
    N = 1000, pd = 0.325, R1 = 1.01,
    father_coverage = 100, mother_coverage = 100, child_coverage = 100,
    epsilon0 = 0, epsilon1 = 0,
    pheno_misclass = TRUE, prev = 0.05, pi01 = 0.05,
    alpha = 0.05, verbose = FALSE
  )
  ordinary <- tdt_power(
    N = 1000, input_mode = "model_based",
    pd = 0.325, prev = 0.05, R1 = 1.01, R2 = 1.01^2,
    alpha = 0.05, delta_prime = 1, misclass_rate = 0.05,
    verbose = FALSE
  )
  ratio <- ordinary$lambda$misclassification / ordinary$lambda$no_error

  expect_equal(
    ngs$scenarios$phenotype_misclassification$attenuation_factor,
    ratio, tolerance = 1e-12
  )
  expect_equal(ngs$scenarios$sequencing_only$lambda,
               ordinary$lambda$no_error, tolerance = 0.002)
  expect_equal(ngs$lambda, ordinary$lambda$misclassification,
               tolerance = 0.002)
})

test_that("phenotype console reports both scenarios and bridge inputs", {
  power_text <- capture.output(
    .paweh_print_tdt_ngs_power(tdt_ngs_pheno_power()), type = "message"
  )
  mssn_text <- capture.output(
    .paweh_print_tdt_ngs_mssn(tdt_ngs_pheno_mssn()), type = "message"
  )
  text <- paste(c(power_text, mssn_text), collapse = "\n")

  expect_match(text, "Sequencing only", fixed = TRUE)
  expect_match(text, "Phenotype misclassification", fixed = TRUE)
  expect_match(text, "ordinary-TDT NCP attenuation bridge", fixed = TRUE)
  expect_match(text, "Unaffected as affected (pi01)", fixed = TRUE)
  expect_match(text, "Prevalence", fixed = TRUE)
  expect_match(text, "Attenuation factor", fixed = TRUE)
  expect_match(text, "Adjusted", fixed = TRUE)
})

test_that("TDT-NGS plot consumers use the active phenotype scenario", {
  power_data <- plot_ngs_power(
    design = "tdt", coverage = c(4, 8), seq_error = 0.01,
    return_data = TRUE, N = 1200, pd = 0.25, R1 = 1.4,
    pheno_misclass = TRUE, prev = 0.05, pi01 = 0.08,
    alpha = 0.05
  )
  mssn_data <- plot_ngs_mssn(
    design = "tdt", coverage = c(4, 8), seq_error = 0.01,
    return_data = TRUE, power = 0.80, pd = 0.25, R1 = 1.4,
    pheno_misclass = TRUE, prev = 0.05, pi01 = 0.08,
    alpha = 0.05
  )
  for (i in seq_len(nrow(power_data))) {
    direct <- tdt_ngs_pheno_power(coverage = power_data$coverage[i])
    expect_equal(power_data$lambda[i], direct$lambda)
    expect_equal(power_data$power[i], direct$power)
  }
  for (i in seq_len(nrow(mssn_data))) {
    direct <- tdt_ngs_pheno_mssn(coverage = mssn_data$coverage[i])
    expect_identical(mssn_data$MSSN_trios[i], direct$MSSN_trios)
    expect_equal(mssn_data$achieved_power[i], direct$achieved_power)
  }
})
