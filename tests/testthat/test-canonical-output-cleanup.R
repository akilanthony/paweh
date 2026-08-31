capture_canonical_messages <- function(expr) {
  paste(capture.output(force(expr), type = "message"), collapse = "\n")
}

tdt_design_args <- list(
  input_mode = "model_based",
  pd = 0.3, prev = 0.05, R1 = 1.5, R2 = 2.25,
  alpha = 0.05, verbose = FALSE
)

cc_design_args <- list(
  alpha = 0.05, input_mode = "model_free",
  g1 = c(0.36, 0.48, 0.16),
  g0 = c(0.25, 0.50, 0.25),
  verbose = FALSE
)

test_that("TDT omitted modifiers equal explicit zero modifiers", {
  power_default <- do.call(tdt_power, c(list(N = 600), tdt_design_args))
  power_zero <- do.call(
    tdt_power,
    c(list(N = 600), tdt_design_args,
      list(misclass_rate = 0, heter_rate = 0))
  )
  expect_equal(power_default, power_zero)

  mssn_default <- do.call(
    tdt_mssn,
    c(list(target_power = 0.8), tdt_design_args)
  )
  mssn_zero <- do.call(
    tdt_mssn,
    c(list(target_power = 0.8), tdt_design_args,
      list(misclass_rate = 0, heter_rate = 0))
  )
  expect_equal(mssn_default, mssn_zero)
})

test_that("TDT verbose output shows only requested separate scenarios", {
  check_output <- function(fun, extra, present = character(), absent = character()) {
    args <- c(
      if (identical(fun, tdt_power)) list(N = 600) else list(target_power = 0.8),
      tdt_design_args[names(tdt_design_args) != "verbose"], extra,
      list(verbose = TRUE)
    )
    text <- capture_canonical_messages(do.call(fun, args))
    expect_match(text, "No-Error Design", fixed = TRUE)
    for (label in present) expect_match(text, label, fixed = TRUE)
    for (label in absent) expect_false(grepl(label, text, fixed = TRUE))
    expect_false(grepl("combined", text, ignore.case = TRUE))
  }

  for (fun in list(tdt_power, tdt_mssn)) {
    check_output(fun, list(), absent = c(
      "Phenotype Misclassification", "Locus Heterogeneity"
    ))
    check_output(fun, list(misclass_rate = 0.05),
                 present = "Phenotype Misclassification",
                 absent = "Locus Heterogeneity")
    check_output(fun, list(heter_rate = 0.1),
                 present = "Locus Heterogeneity",
                 absent = "Phenotype Misclassification")
    check_output(fun, list(misclass_rate = 0.05, heter_rate = 0.1),
                 present = c("Phenotype Misclassification", "Locus Heterogeneity"))
  }
})

test_that("case-control omitted modifiers equal explicit no-error modifiers", {
  power_default <- do.call(cc_power, c(list(N_case = 400), cc_design_args))
  power_zero <- do.call(
    cc_power,
    c(list(N_case = 400), cc_design_args,
      list(locus_het = FALSE, pi = 1, pheno_misclass = FALSE,
           theta = 0, phi = 0, geno_misclass = "none"))
  )
  expect_equal(power_default, power_zero)

  mssn_default <- do.call(cc_mssn, c(list(power = 0.8), cc_design_args))
  mssn_zero <- do.call(
    cc_mssn,
    c(list(power = 0.8), cc_design_args,
      list(locus_het = FALSE, pi = 1, pheno_misclass = FALSE,
           theta = 0, phi = 0, geno_misclass = "none"))
  )
  expect_equal(mssn_default, mssn_zero)
})

test_that("case-control verbose output reports baseline then one adjusted design", {
  base <- cc_design_args[names(cc_design_args) != "verbose"]
  for (fun in list(cc_power, cc_mssn)) {
    first <- if (identical(fun, cc_power)) list(N_case = 400) else list(power = 0.8)
    quiet_text <- capture_canonical_messages(do.call(fun, c(first, base, list(verbose = TRUE))))
    expect_match(quiet_text, "No-Error Design", fixed = TRUE)
    expect_false(grepl("Adjusted Design", quiet_text, fixed = TRUE))

    adjusted_text <- capture_canonical_messages(do.call(
      fun,
      c(first, base, list(
        prev = 0.05,
        locus_het = TRUE, pi = 0.8,
        pheno_misclass = TRUE, theta = 0.02, phi = 0.01,
        geno_misclass = "2p", e1 = 0.02, e2 = 0.01,
        verbose = TRUE
      ))
    ))
    expect_match(adjusted_text, "No-Error Design", fixed = TRUE)
    expect_match(adjusted_text, "Adjusted Design", fixed = TRUE)
    expect_match(adjusted_text, "Locus heterogeneity", fixed = TRUE)
    expect_match(adjusted_text, "Phenotype misclassification", fixed = TRUE)
    expect_match(adjusted_text, "Genotype misclassification", fixed = TRUE)
  }
})

test_that("zero-valued genotype error models remain no-error designs", {
  base <- cc_design_args[names(cc_design_args) != "verbose"]
  zero_models <- list(
    `1p` = list(geno_misclass = "1p", e = 0),
    `2p` = list(geno_misclass = "2p", e1 = 0, e2 = 0),
    `3p` = list(geno_misclass = "3p", e01 = 0, e02 = 0, e03 = 0),
    diff3p = list(
      geno_misclass = "diff3p",
      case_e01 = 0, case_e02 = 0, case_e03 = 0,
      ctrl_e01 = 0, ctrl_e02 = 0, ctrl_e03 = 0
    )
  )

  for (fun in list(cc_power, cc_mssn)) {
    first <- if (identical(fun, cc_power)) list(N_case = 400) else list(power = 0.8)
    no_error <- do.call(fun, c(first, base, list(verbose = FALSE)))

    for (model in zero_models) {
      zero_error <- do.call(fun, c(first, base, model, list(verbose = FALSE)))
      expect_equal(zero_error$tests, no_error$tests)
      expect_equal(zero_error$freqs$g_obs_case, no_error$freqs$g_obs_case)
      expect_equal(zero_error$freqs$g_obs_ctrl, no_error$freqs$g_obs_ctrl)
      expect_false(zero_error$errors$genotype_misclass$enabled)

      text <- capture_canonical_messages(
        do.call(fun, c(first, base, model, list(verbose = TRUE)))
      )
      expect_match(text, "No-Error Design", fixed = TRUE)
      expect_false(grepl("Adjusted Design", text, fixed = TRUE))
      expect_false(grepl("Genotype misclassification", text, fixed = TRUE))
    }
  }
})

test_that("nonzero effective genotype error triggers adjusted narration", {
  base <- cc_design_args[names(cc_design_args) != "verbose"]
  for (fun in list(cc_power, cc_mssn)) {
    first <- if (identical(fun, cc_power)) list(N_case = 400) else list(power = 0.8)
    text <- capture_canonical_messages(do.call(
      fun,
      c(first, base, list(geno_misclass = "1p", e = 0.01, verbose = TRUE))
    ))
    expect_match(text, "Adjusted Design", fixed = TRUE)
    expect_match(text, "Genotype misclassification", fixed = TRUE)
  }
})

test_that("TDT percent increase handles a non-finite no-error baseline", {
  out <- tdt_mssn(
    target_power = 0.8, input_mode = "model_free",
    ET = 100, ENT = 100, n_trios = 100,
    pd = 0.3, heter_rate = 0.1, verbose = FALSE
  )
  expect_true(is.infinite(out$N$no_error))
  expect_true(is.na(out$percent_increase$misclassification))
  expect_true(is.na(out$percent_increase$heterogeneity))

  text <- capture_canonical_messages(tdt_mssn(
    target_power = 0.8, input_mode = "model_free",
    ET = 100, ENT = 100, n_trios = 100,
    pd = 0.3, heter_rate = 0.1, verbose = TRUE
  ))
  expect_match(text, "Required complete trios:.*Inf")
  expect_match(text, "MSSN increase:.*Not defined")
})

test_that("canonical verbose FALSE calls emit no narrative", {
  expect_silent(do.call(tdt_power, c(list(N = 600), tdt_design_args)))
  expect_silent(do.call(tdt_mssn, c(list(target_power = 0.8), tdt_design_args)))
  expect_silent(do.call(cc_power, c(list(N_case = 400), cc_design_args)))
  expect_silent(do.call(cc_mssn, c(list(power = 0.8), cc_design_args)))
})
