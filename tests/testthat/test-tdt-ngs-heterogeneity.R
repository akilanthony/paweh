tdt_ngs_het_power <- function(...) {
  args <- list(
    N = 1200, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, alpha = 0.05, verbose = FALSE,
    locus_het = TRUE, prev = 0.05, heter_rate = 0.25
  )
  do.call(tdt_ngs_power, utils::modifyList(args, list(...)))
}

tdt_ngs_het_mssn <- function(...) {
  args <- list(
    power = 0.80, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, alpha = 0.05, verbose = FALSE,
    locus_het = TRUE, prev = 0.05, heter_rate = 0.25
  )
  do.call(tdt_ngs_mssn, utils::modifyList(args, list(...)))
}

test_that("heterogeneity sensitivity creates the exact scenario architecture", {
  inactive <- tdt_ngs_power(
    N = 1200, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, verbose = FALSE
  )
  active <- tdt_ngs_het_power()

  expect_named(inactive$scenarios, "sequencing_only")
  expect_identical(inactive$compatibility_scenario, "sequencing_only")
  expect_named(active$scenarios, c("sequencing_only", "heterogeneity"))
  expect_identical(active$compatibility_scenario, "heterogeneity")
  expect_identical(active$lambda, active$scenarios$heterogeneity$lambda)
  expect_identical(active$power, active$scenarios$heterogeneity$power)
  expect_identical(active$scenarios$sequencing_only$lambda, inactive$lambda)
  expect_identical(active$information_matrix, inactive$information_matrix)
  expect_identical(active$scenarios$heterogeneity$modifier, "heterogeneity")
  expect_identical(active$scenarios$heterogeneity$label,
                   "Locus heterogeneity")
})

test_that("ordinary Eq. 5.31-5.32 counts reconstruct the heterogeneity bridge", {
  out <- tdt_ngs_het_power()
  bridge <- out$scenarios$heterogeneity$ordinary_tdt_bridge
  fixed <- list(
    N_star = 1, pd = out$pd, prev = 0.05, R1 = out$R1,
    R2 = out$R1^2, delta_prime = 1, theta1 = out$pd, verbose = FALSE
  )
  homogeneous <- do.call(
    tdt_expected_transmission_counts, c(fixed, list(pi = 1))
  )
  heterogeneous <- do.call(
    tdt_expected_transmission_counts, c(fixed, list(pi = 0.75))
  )
  gT0 <- homogeneous$ET_star / 2
  gNT0 <- homogeneous$ENT_star / 2
  gT1 <- heterogeneous$ET_star / 2
  gNT1 <- heterogeneous$ENT_star / 2
  c0 <- 2 * (gT0 - gNT0)^2 / (gT0 + gNT0)
  c1 <- 2 * (gT1 - gNT1)^2 / (gT1 + gNT1)

  expect_equal(bridge$gT_noerror, gT0, tolerance = 1e-15)
  expect_equal(bridge$gNT_noerror, gNT0, tolerance = 1e-15)
  expect_equal(bridge$gT_heterogeneity, gT1, tolerance = 1e-15)
  expect_equal(bridge$gNT_heterogeneity, gNT1, tolerance = 1e-15)
  expect_equal(bridge$ncp_per_trio_noerror, c0, tolerance = 1e-15)
  expect_equal(bridge$ncp_per_trio_heterogeneity, c1, tolerance = 1e-15)
  expect_equal(bridge$attenuation_factor, c1 / c0, tolerance = 1e-15)
  expect_identical(bridge$heter_rate, 0.25)
  expect_identical(bridge$effective_pi, 0.75)
  expect_identical(bridge$R2, out$R1^2)
  expect_identical(bridge$delta_prime, 1)
  expect_identical(bridge$theta1, out$pd)
  expect_identical(bridge$method, "ordinary_TDT_NCP_attenuation_bridge")
  expect_identical(bridge$raw_read_likelihood_modified, FALSE)
})

test_that("heterogeneity NCP and power are independently reproducible", {
  out <- tdt_ngs_het_power()
  base <- out$scenarios$sequencing_only
  adjusted <- out$scenarios$heterogeneity
  expected_lambda <- base$lambda * adjusted$attenuation_factor
  critical <- stats::qchisq(1 - out$alpha, df = 1)
  expected_power <- stats::pchisq(
    critical, df = 1, ncp = expected_lambda, lower.tail = FALSE
  )

  expect_equal(adjusted$base_lambda, base$lambda, tolerance = 1e-15)
  expect_equal(adjusted$lambda, expected_lambda, tolerance = 1e-14)
  expect_equal(adjusted$power, expected_power, tolerance = 1e-15)
  expect_lt(adjusted$lambda, base$lambda)
  expect_lt(adjusted$power, base$power)
  expect_identical(
    adjusted$model_info$heterogeneity_adjustment,
    "ordinary_TDT_NCP_attenuation_bridge"
  )
  expect_identical(adjusted$model_info$raw_read_likelihood_modified, FALSE)
})

test_that("heterogeneity MSSN attenuates before inversion and is minimal", {
  out <- tdt_ngs_het_mssn()
  base <- out$scenarios$sequencing_only
  adjusted <- out$scenarios$heterogeneity
  expected_coefficient <- base$ncp_per_trio * adjusted$attenuation_factor
  expected_continuous <- adjusted$lambda_target / expected_coefficient

  expect_equal(adjusted$ncp_per_trio, expected_coefficient, tolerance = 1e-15)
  expect_equal(adjusted$N_trios_continuous, expected_continuous,
               tolerance = 1e-12)
  expect_gt(adjusted$MSSN_trios, base$MSSN_trios)
  expect_identical(out$MSSN_trios, adjusted$MSSN_trios)
  expect_gte(adjusted$achieved_power, out$power_target)

  at_plan <- tdt_ngs_het_power(N = adjusted$MSSN_trios)
  expect_equal(at_plan$power, adjusted$achieved_power, tolerance = 1e-14)
  expect_gte(at_plan$power, out$power_target)
  one_lower <- tdt_ngs_het_power(N = adjusted$MSSN_trios - 1)
  expect_lt(one_lower$power, out$power_target)
})

test_that("heterogeneity factor is invariant to sequencing design", {
  equal <- tdt_ngs_het_power()
  unequal <- tdt_ngs_het_power(
    coverage = NULL, seq_error = NULL,
    father_coverage = 4, mother_coverage = 10, child_coverage = 18,
    epsilon0 = 0.004, epsilon1 = 0.025
  )
  e <- equal$scenarios$heterogeneity
  u <- unequal$scenarios$heterogeneity

  expect_identical(e$attenuation_factor, u$attenuation_factor)
  expect_identical(e$ordinary_tdt_bridge, u$ordinary_tdt_bridge)
  expect_equal(e$lambda, e$base_lambda * e$attenuation_factor,
               tolerance = 1e-14)
  expect_equal(u$lambda, u$base_lambda * u$attenuation_factor,
               tolerance = 1e-14)
  expect_identical(u$sequencing$coverage,
                   c(father = 4, mother = 10, child = 18))
  expect_identical(u$sequencing$epsilon0, 0.004)
  expect_identical(u$sequencing$epsilon1, 0.025)
})

test_that("heter_rate zero remains a visible identity scenario", {
  power <- tdt_ngs_het_power(heter_rate = 0)
  mssn <- tdt_ngs_het_mssn(heter_rate = 0)
  pbase <- power$scenarios$sequencing_only
  padjusted <- power$scenarios$heterogeneity
  mbase <- mssn$scenarios$sequencing_only
  madjusted <- mssn$scenarios$heterogeneity

  expect_identical(padjusted$attenuation_factor, 1)
  expect_identical(padjusted$lambda, pbase$lambda)
  expect_identical(padjusted$power, pbase$power)
  expect_identical(padjusted$ordinary_tdt_bridge$effective_pi, 1)
  expect_identical(padjusted$ordinary_tdt_bridge$status,
                   "identity_heter_rate_zero")
  expect_identical(madjusted$MSSN_trios, mbase$MSSN_trios)
  expect_identical(madjusted$N_trios_continuous, mbase$N_trios_continuous)
})

test_that("phenotype and heterogeneity remain independent scenarios", {
  common_power <- list(
    N = 1200, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, alpha = 0.05, verbose = FALSE, prev = 0.05
  )
  common_mssn <- c(list(power = 0.80),
                   common_power[setdiff(names(common_power), "N")])
  phenotype_power <- do.call(
    tdt_ngs_power,
    c(common_power, list(pheno_misclass = TRUE, pi01 = 0.08))
  )
  heterogeneity_power <- do.call(
    tdt_ngs_power,
    c(common_power, list(locus_het = TRUE, heter_rate = 0.25))
  )
  both_power <- do.call(
    tdt_ngs_power,
    c(common_power, list(
      pheno_misclass = TRUE, pi01 = 0.08,
      locus_het = TRUE, heter_rate = 0.25
    ))
  )
  baseline_power <- do.call(tdt_ngs_power, common_power)

  expect_named(
    both_power$scenarios,
    c("sequencing_only", "phenotype_misclassification", "heterogeneity")
  )
  expect_identical(both_power$compatibility_scenario, "sequencing_only")
  expect_identical(both_power$scenarios$sequencing_only,
                   baseline_power$scenarios$sequencing_only)
  expect_identical(both_power$scenarios$phenotype_misclassification,
                   phenotype_power$scenarios$phenotype_misclassification)
  expect_identical(both_power$scenarios$heterogeneity,
                   heterogeneity_power$scenarios$heterogeneity)
  expect_identical(both_power$lambda,
                   both_power$scenarios$sequencing_only$lambda)
  expect_false(any(grepl("combined|plus", names(both_power$scenarios))))

  phenotype_mssn <- do.call(
    tdt_ngs_mssn,
    c(common_mssn, list(pheno_misclass = TRUE, pi01 = 0.08))
  )
  heterogeneity_mssn <- do.call(
    tdt_ngs_mssn,
    c(common_mssn, list(locus_het = TRUE, heter_rate = 0.25))
  )
  both_mssn <- do.call(
    tdt_ngs_mssn,
    c(common_mssn, list(
      pheno_misclass = TRUE, pi01 = 0.08,
      locus_het = TRUE, heter_rate = 0.25
    ))
  )
  baseline_mssn <- do.call(tdt_ngs_mssn, common_mssn)

  expect_named(
    both_mssn$scenarios,
    c("sequencing_only", "phenotype_misclassification", "heterogeneity")
  )
  expect_identical(both_mssn$compatibility_scenario, "sequencing_only")
  expect_identical(both_mssn$scenarios$sequencing_only,
                   baseline_mssn$scenarios$sequencing_only)
  expect_identical(both_mssn$scenarios$phenotype_misclassification,
                   phenotype_mssn$scenarios$phenotype_misclassification)
  expect_identical(both_mssn$scenarios$heterogeneity,
                   heterogeneity_mssn$scenarios$heterogeneity)
  expect_identical(both_mssn$MSSN_trios,
                   both_mssn$scenarios$sequencing_only$MSSN_trios)
})

test_that("both identity modifiers remain explicit with baseline compatibility", {
  power <- tdt_ngs_power(
    N = 1200, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, verbose = FALSE, prev = 0.05,
    pheno_misclass = TRUE, pi01 = 0,
    locus_het = TRUE, heter_rate = 0
  )
  mssn <- tdt_ngs_mssn(
    power = 0.80, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, verbose = FALSE, prev = 0.05,
    pheno_misclass = TRUE, pi01 = 0,
    locus_het = TRUE, heter_rate = 0
  )
  expected <- c("sequencing_only", "phenotype_misclassification",
                "heterogeneity")

  expect_named(power$scenarios, expected)
  expect_named(mssn$scenarios, expected)
  expect_identical(power$compatibility_scenario, "sequencing_only")
  expect_identical(mssn$compatibility_scenario, "sequencing_only")
  expect_true(all(vapply(power$scenarios, function(x) {
    identical(x$lambda, power$scenarios$sequencing_only$lambda)
  }, logical(1))))
  expect_true(all(vapply(mssn$scenarios, function(x) {
    identical(x$MSSN_trios, mssn$scenarios$sequencing_only$MSSN_trios)
  }, logical(1))))
})

test_that("complete heterogeneity has zero power signal and no finite MSSN", {
  out <- tdt_ngs_het_power(heter_rate = 1)
  adjusted <- out$scenarios$heterogeneity

  expect_identical(adjusted$ordinary_tdt_bridge$effective_pi, 0)
  expect_equal(adjusted$attenuation_factor, 0, tolerance = 1e-15)
  expect_equal(adjusted$ncp_per_trio, 0, tolerance = 1e-15)
  expect_equal(adjusted$lambda, 0, tolerance = 1e-15)
  expect_equal(adjusted$power, out$alpha, tolerance = 1e-15)
  expect_error(
    tdt_ngs_het_mssn(heter_rate = 1),
    "heterogeneity scenario has zero per-trio NCP",
    fixed = TRUE
  )
})

test_that("R1 one safely reports all independent null scenarios", {
  out <- expect_no_warning(tdt_ngs_power(
    N = 1200, pd = 0.25, R1 = 1, coverage = 8,
    seq_error = 0.01, alpha = 0.05, verbose = FALSE,
    pheno_misclass = TRUE, prev = 0.05, pi01 = 0.08,
    locus_het = TRUE, heter_rate = 0.25
  ))
  expect_named(out$scenarios, c(
    "sequencing_only", "phenotype_misclassification", "heterogeneity"
  ))
  expect_true(all(vapply(out$scenarios, function(x) identical(x$lambda, 0),
                         logical(1))))
  expect_equal(unname(vapply(out$scenarios, `[[`, numeric(1), "power")),
               rep(out$alpha, 3), tolerance = 1e-15)
  expect_identical(
    out$scenarios$phenotype_misclassification$ordinary_tdt_bridge$status,
    "zero_baseline_effect_identity"
  )
  expect_identical(
    out$scenarios$heterogeneity$ordinary_tdt_bridge$status,
    "zero_baseline_effect_identity"
  )
  expect_error(
    tdt_ngs_mssn(
      power = 0.80, pd = 0.25, R1 = 1, coverage = 8,
      seq_error = 0.01, alpha = 0.05, verbose = FALSE,
      pheno_misclass = TRUE, prev = 0.05, pi01 = 0.08,
      locus_het = TRUE, heter_rate = 0.25
    ),
    "R1 = 1 implies zero transmission effect",
    fixed = TRUE
  )
})

test_that("heterogeneity bridge inputs are validated", {
  base <- list(
    N = 100, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, verbose = FALSE
  )
  invalid <- list(
    list(locus_het = 1),
    list(locus_het = NA),
    list(locus_het = c(TRUE, FALSE)),
    list(locus_het = TRUE),
    list(locus_het = TRUE, prev = 0),
    list(locus_het = TRUE, prev = 1),
    list(locus_het = TRUE, prev = NA_real_),
    list(locus_het = TRUE, prev = c(0.01, 0.02)),
    list(locus_het = TRUE, prev = 0.05, heter_rate = -0.01),
    list(locus_het = TRUE, prev = 0.05, heter_rate = 1.01),
    list(locus_het = TRUE, prev = 0.05, heter_rate = NA_real_),
    list(locus_het = FALSE, heter_rate = 0.01)
  )
  for (values in invalid) {
    expect_error(do.call(tdt_ngs_power, utils::modifyList(base, values)))
  }
  with_prev <- do.call(
    tdt_ngs_power, utils::modifyList(base, list(prev = 0.05))
  )
  without_prev <- do.call(tdt_ngs_power, base)
  expect_identical(with_prev, without_prev)
})

test_that("high-quality heterogeneity TDT1-NGS approaches ordinary TDT", {
  ngs <- tdt_ngs_power(
    N = 1000, pd = 0.325, R1 = 1.01,
    father_coverage = 100, mother_coverage = 100, child_coverage = 100,
    epsilon0 = 0, epsilon1 = 0,
    locus_het = TRUE, prev = 0.05, heter_rate = 0.25,
    alpha = 0.05, verbose = FALSE
  )
  ordinary <- tdt_power(
    N = 1000, input_mode = "model_based",
    pd = 0.325, prev = 0.05, R1 = 1.01, R2 = 1.01^2,
    alpha = 0.05, delta_prime = 1, heter_rate = 0.25,
    verbose = FALSE
  )
  ratio <- ordinary$lambda$heterogeneity / ordinary$lambda$no_error

  expect_equal(ngs$scenarios$heterogeneity$attenuation_factor,
               ratio, tolerance = 1e-12)
  expect_equal(ngs$scenarios$sequencing_only$lambda,
               ordinary$lambda$no_error, tolerance = 0.002)
  expect_equal(ngs$lambda, ordinary$lambda$heterogeneity, tolerance = 0.002)
})

test_that("heterogeneity console reports independent bridge scenarios", {
  heter_power <- capture.output(
    .paweh_print_tdt_ngs_power(tdt_ngs_het_power()), type = "message"
  )
  heter_mssn <- capture.output(
    .paweh_print_tdt_ngs_mssn(tdt_ngs_het_mssn()), type = "message"
  )
  both <- tdt_ngs_power(
    N = 1200, pd = 0.25, R1 = 1.4, coverage = 8,
    seq_error = 0.01, verbose = FALSE, prev = 0.05,
    pheno_misclass = TRUE, pi01 = 0.08,
    locus_het = TRUE, heter_rate = 0.25
  )
  both_text <- capture.output(.paweh_print_tdt_ngs_power(both), type = "message")
  text <- paste(c(heter_power, heter_mssn, both_text), collapse = "\n")

  expect_match(text, "Sequencing only", fixed = TRUE)
  expect_match(text, "Locus heterogeneity", fixed = TRUE)
  expect_match(text, "Heterogeneity adjustment", fixed = TRUE)
  expect_match(text, "ordinary-TDT NCP attenuation bridge", fixed = TRUE)
  expect_match(text, "heter_rate", fixed = TRUE)
  expect_match(text, "Linked/homogeneous fraction (pi)", fixed = TRUE)
  expect_match(paste(both_text, collapse = "\n"),
               "Phenotype misclassification", fixed = TRUE)
  expect_match(paste(both_text, collapse = "\n"),
               "Modifiers are evaluated separately", fixed = TRUE)
})

test_that("NGS plot consumers honor heterogeneity compatibility selection", {
  heter_data <- plot_ngs_power(
    design = "tdt", coverage = c(4, 8), seq_error = 0.01,
    return_data = TRUE, N = 1200, pd = 0.25, R1 = 1.4,
    locus_het = TRUE, prev = 0.05, heter_rate = 0.25,
    alpha = 0.05
  )
  both_data <- plot_ngs_power(
    design = "tdt", coverage = c(4, 8), seq_error = 0.01,
    return_data = TRUE, N = 1200, pd = 0.25, R1 = 1.4,
    pheno_misclass = TRUE, pi01 = 0.08,
    locus_het = TRUE, prev = 0.05, heter_rate = 0.25,
    alpha = 0.05
  )
  heter_mssn <- plot_ngs_mssn(
    design = "tdt", coverage = c(4, 8), seq_error = 0.01,
    return_data = TRUE, power = 0.80, pd = 0.25, R1 = 1.4,
    locus_het = TRUE, prev = 0.05, heter_rate = 0.25,
    alpha = 0.05
  )
  both_mssn <- plot_ngs_mssn(
    design = "tdt", coverage = c(4, 8), seq_error = 0.01,
    return_data = TRUE, power = 0.80, pd = 0.25, R1 = 1.4,
    pheno_misclass = TRUE, pi01 = 0.08,
    locus_het = TRUE, prev = 0.05, heter_rate = 0.25,
    alpha = 0.05
  )
  for (i in seq_len(nrow(heter_data))) {
    heter <- tdt_ngs_het_power(coverage = heter_data$coverage[[i]])
    baseline <- tdt_ngs_power(
      N = 1200, pd = 0.25, R1 = 1.4,
      coverage = both_data$coverage[[i]], seq_error = 0.01,
      alpha = 0.05, verbose = FALSE
    )
    expect_equal(heter_data$lambda[[i]], heter$lambda)
    expect_equal(heter_data$power[[i]], heter$power)
    expect_equal(both_data$lambda[[i]], baseline$lambda)
    expect_equal(both_data$power[[i]], baseline$power)
    adjusted_mssn <- tdt_ngs_het_mssn(coverage = heter_mssn$coverage[[i]])
    baseline_mssn <- tdt_ngs_mssn(
      power = 0.80, pd = 0.25, R1 = 1.4,
      coverage = both_mssn$coverage[[i]], seq_error = 0.01,
      alpha = 0.05, verbose = FALSE
    )
    expect_identical(heter_mssn$MSSN_trios[[i]], adjusted_mssn$MSSN_trios)
    expect_equal(heter_mssn$achieved_power[[i]], adjusted_mssn$achieved_power)
    expect_identical(both_mssn$MSSN_trios[[i]], baseline_mssn$MSSN_trios)
    expect_equal(both_mssn$achieved_power[[i]], baseline_mssn$achieved_power)
  }
})
