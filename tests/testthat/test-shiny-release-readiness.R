test_that("Case-Control release matrix renders and reproduces canonical results", {
  designs <- list(
    power_no_modifiers = list(),
    mssn_all_modifiers = list(
      objective = "mssn", locus_het = TRUE, pi = 75,
      pheno_misclass = TRUE, theta = 2, phi = 1,
      genotype_error = TRUE, geno_misclass = "diff3p",
      case_e01 = 2, case_e02 = 1, case_e03 = 0.5,
      ctrl_e01 = 1, ctrl_e02 = 0.5, ctrl_e03 = 0.25
    ),
    direct_genotypes = list(input_mode = "model_free")
  )
  for (name in names(designs)) {
    values <- modifyList(pawh:::.pawh_cc_defaults(), designs[[name]])
    calculation <- pawh:::.pawh_cc_calculate(pawh:::.pawh_cc_snapshot(values))
    fun <- if (values$objective == "mssn") cc_mssn else cc_power
    reproduced <- do.call(fun, pawh:::.pawh_cc_repro_args(calculation))
    expect_equal(reproduced$tests, calculation$adjusted$tests, info = name)
    expect_error(pawh:::.pawh_cc_results_ui(calculation), NA, info = name)
    expect_error(pawh:::.pawh_cc_advanced_ui(calculation), NA, info = name)
    expect_error(pawh:::.pawh_cc_methods_ui(calculation), NA, info = name)
    expect_match(pawh:::.pawh_call_text(pawh:::.pawh_cc_repro_call(calculation)),
      if (values$objective == "mssn") "cc_mssn" else "cc_power", fixed = TRUE)
  }
})

test_that("TDT release matrix renders and reproduces canonical results", {
  designs <- list(
    power_no_modifiers = list(),
    mssn_both_scenarios = list(
      objective = "mssn", misclassification = TRUE, misclass_rate = 5,
      heterogeneity = TRUE, heter_rate = 20
    ),
    direct_transmissions = list(input_mode = "model_free"),
    null_effect_mssn = list(
      objective = "mssn", input_mode = "model_free", ET = 100, ENT = 100
    )
  )
  for (name in names(designs)) {
    values <- modifyList(pawh:::.pawh_tdt_defaults(), designs[[name]])
    calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
    fun <- if (values$objective == "mssn") tdt_mssn else tdt_power
    reproduced <- do.call(fun, pawh:::.pawh_tdt_repro_args(calculation))
    if (values$objective == "mssn") {
      expect_equal(reproduced$N, calculation$result$N, info = name)
    } else {
      expect_equal(reproduced$power, calculation$result$power, info = name)
    }
    expect_error(pawh:::.pawh_tdt_results_ui(calculation), NA, info = name)
    expect_error(pawh:::.pawh_tdt_advanced_ui(calculation), NA, info = name)
    expect_error(pawh:::.pawh_tdt_methods_ui(calculation), NA, info = name)
    expect_match(pawh:::.pawh_call_text(pawh:::.pawh_tdt_repro_call(calculation)),
      if (values$objective == "mssn") "tdt_mssn" else "tdt_power", fixed = TRUE)
  }
})

test_that("QTL release matrix renders and reproduces canonical results", {
  designs <- list(
    continuous_power = list(subtype = "continuous", objective = "power"),
    continuous_mssn = list(subtype = "continuous", objective = "mssn"),
    extreme_power = list(subtype = "extreme", objective = "power"),
    extreme_mssn = list(subtype = "extreme", objective = "mssn"),
    multivariate_pillai_power = list(subtype = "multivariate", objective = "power"),
    multivariate_pillai_mssn = list(subtype = "multivariate", objective = "mssn"),
    multivariate_threshold_power = list(
      subtype = "multivariate", objective = "power", mv_test = "threshold_chisq"
    ),
    multivariate_threshold_mssn = list(
      subtype = "multivariate", objective = "mssn", mv_test = "threshold_chisq"
    )
  )
  for (name in names(designs)) {
    values <- modifyList(pawh:::.pawh_qtl_defaults(), designs[[name]])
    calculation <- pawh:::.pawh_qtl_calculate(pawh:::.pawh_qtl_snapshot(values))
    fun <- get(pawh:::.pawh_qtl_function(calculation$snapshot), envir = asNamespace("pawh"))
    reproduced <- do.call(fun, pawh:::.pawh_qtl_repro_args(calculation))
    if (values$objective == "power") {
      expect_equal(reproduced$power, calculation$result$power, info = name)
    } else if (values$subtype == "continuous" ||
               values$subtype == "multivariate" && values$mv_test == "pillai") {
      expect_equal(reproduced$N, calculation$result$N, info = name)
    } else {
      expect_equal(reproduced$N_total, calculation$result$N_total, info = name)
    }
    expect_error(pawh:::.pawh_qtl_results_ui(calculation), NA, info = name)
    expect_error(pawh:::.pawh_qtl_advanced_ui(calculation), NA, info = name)
    expect_error(pawh:::.pawh_qtl_methods_ui(calculation), NA, info = name)
    expect_match(
      pawh:::.pawh_call_text(pawh:::.pawh_qtl_repro_call(calculation)),
      pawh:::.pawh_qtl_function(calculation$snapshot), fixed = TRUE
    )
  }
})
