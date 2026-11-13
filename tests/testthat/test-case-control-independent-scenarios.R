cc_scenario_g_case <- c(0.25, 0.50, 0.25)
cc_scenario_g_ctrl <- c(0.36, 0.48, 0.16)

cc_scenario_base <- function(mode = c("model_free", "model_based")) {
  mode <- match.arg(mode)
  if (mode == "model_free") {
    list(
      input_mode = mode, g1 = cc_scenario_g_case, g0 = cc_scenario_g_ctrl,
      prev = 0.05, k = 1, w = c(0, 1, 2), verbose = FALSE
    )
  } else {
    list(
      input_mode = mode, prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M",
      k = 1, w = c(0, 1, 2), verbose = FALSE
    )
  }
}

cc_scenario_modifiers <- list(
  phenotype = list(pheno_misclass = TRUE, theta = 0.05, phi = 0.01),
  heterogeneity = list(locus_het = TRUE, pi = 0.80),
  genotype = list(
    geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005
  )
)

cc_scenario_combine <- function(modifiers = names(cc_scenario_modifiers)) {
  do.call(c, unname(cc_scenario_modifiers[modifiers]))
}

cc_scenario_call <- function(fun = c("power", "mssn"), mode = "model_free",
                             modifiers = list()) {
  fun <- match.arg(fun)
  leading <- if (fun == "power") {
    list(N_case = 500, alpha = 0.05)
  } else {
    list(power = 0.80, alpha = 0.05)
  }
  do.call(
    if (fun == "power") cc_power else cc_mssn,
    c(leading, cc_scenario_base(mode), modifiers)
  )
}

test_that("ordinary CC no-error results retain the compatibility structure", {
  for (mode in c("model_based", "model_free")) {
    for (fun in c("power", "mssn")) {
      out <- cc_scenario_call(fun, mode)
      expect_named(out$scenarios, "no_error")
      expect_identical(out$compatibility_scenario, "no_error")
      expect_identical(out$tests, out$scenarios$no_error$tests)
      expect_identical(out$freqs$g_base_case,
                       out$scenarios$no_error$freqs$g_case_base)
      expect_identical(out$freqs$g_base_ctrl,
                       out$scenarios$no_error$freqs$g_ctrl_base)
    }
  }
})

test_that("each ordinary CC modifier retains its standalone result", {
  names_by_modifier <- c(
    phenotype = "phenotype_misclassification",
    heterogeneity = "heterogeneity",
    genotype = "genotype_misclassification"
  )
  for (mode in c("model_based", "model_free")) {
    for (fun in c("power", "mssn")) {
      for (modifier in names(cc_scenario_modifiers)) {
        out <- cc_scenario_call(fun, mode, cc_scenario_modifiers[[modifier]])
        scenario <- names_by_modifier[[modifier]]
        expect_named(out$scenarios, c("no_error", scenario))
        expect_identical(out$compatibility_scenario, scenario)
        expect_identical(out$tests, out$scenarios[[scenario]]$tests)
      }
    }
  }
})

test_that("all genotype-error models retain independent standalone scenarios", {
  models <- list(
    `1p` = list(geno_misclass = "1p", e = 0.02),
    `2p` = list(geno_misclass = "2p", e1 = 0.02, e2 = 0.01),
    `3p` = cc_scenario_modifiers$genotype,
    diff3p = list(
      geno_misclass = "diff3p", diff_source = "explicit",
      case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
      ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002
    )
  )
  for (fun in c("power", "mssn")) {
    for (model in models) {
      out <- cc_scenario_call(fun, modifiers = model)
      expect_named(
        out$scenarios,
        c("no_error", "genotype_misclassification")
      )
      expect_identical(
        out$tests,
        out$scenarios$genotype_misclassification$tests
      )
    }
  }
})

test_that("every pair of modifiers creates only independent scenarios", {
  pairs <- list(
    c("phenotype", "heterogeneity"),
    c("phenotype", "genotype"),
    c("heterogeneity", "genotype")
  )
  scenario_names <- c(
    phenotype = "phenotype_misclassification",
    heterogeneity = "heterogeneity",
    genotype = "genotype_misclassification"
  )
  for (fun in c("power", "mssn")) {
    for (pair in pairs) {
      modifiers <- cc_scenario_combine(pair)
      combined_call <- cc_scenario_call(fun, modifiers = modifiers)
      expect_named(
        combined_call$scenarios,
        c("no_error", unname(scenario_names[pair]))
      )
      expect_identical(combined_call$compatibility_scenario, "no_error")
      for (modifier in pair) {
        standalone <- cc_scenario_call(
          fun, modifiers = cc_scenario_modifiers[[modifier]]
        )
        scenario <- scenario_names[[modifier]]
        expect_identical(
          combined_call$scenarios[[scenario]],
          standalone$scenarios[[scenario]]
        )
      }
    }
  }
})

test_that("all three modifiers create exactly four independent scenarios", {
  modifiers <- cc_scenario_combine()
  expected <- c(
    "no_error", "phenotype_misclassification", "heterogeneity",
    "genotype_misclassification"
  )
  for (mode in c("model_based", "model_free")) {
    for (fun in c("power", "mssn")) {
      out <- cc_scenario_call(fun, mode, modifiers)
      expect_named(out$scenarios, expected)
      expect_identical(out$compatibility_scenario, "no_error")
      expect_identical(out$tests, out$scenarios$no_error$tests)
      for (modifier in names(cc_scenario_modifiers)) {
        scenario <- switch(
          modifier,
          phenotype = "phenotype_misclassification",
          heterogeneity = "heterogeneity",
          genotype = "genotype_misclassification"
        )
        standalone <- cc_scenario_call(
          fun, mode, cc_scenario_modifiers[[modifier]]
        )
        expect_identical(out$scenarios[[scenario]],
                         standalone$scenarios[[scenario]])
      }
      expect_false(any(grepl("combined|sequential", names(out$scenarios))))
    }
  }
})

test_that("all ordinary CC scenarios record the same baseline provenance", {
  out <- cc_scenario_call(
    "power", modifiers = cc_scenario_combine()
  )
  baseline <- out$scenarios$no_error$freqs
  for (scenario in out$scenarios[-1]) {
    expect_identical(scenario$freqs$g_case_input, baseline$g_case_base)
    expect_identical(scenario$freqs$g_ctrl_input, baseline$g_ctrl_base)
  }
  expect_equal(
    out$scenarios$heterogeneity$freqs$g_case,
    0.8 * baseline$g_case_base + 0.2 * baseline$g_ctrl_base,
    tolerance = 1e-15
  )
})

test_that("differential genotype error applies group matrices to baseline", {
  args <- list(
    geno_misclass = "diff3p", diff_source = "explicit",
    case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
    ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002,
    locus_het = TRUE, pi = 0.7,
    pheno_misclass = TRUE, theta = 0.04, phi = 0.01
  )
  out <- cc_scenario_call("power", modifiers = args)
  scenario <- out$scenarios$genotype_misclassification
  expect_equal(
    scenario$freqs$g_case,
    paweh:::.cc_apply_genotype_misclass(
      scenario$freqs$g_case_base, scenario$freqs$M_case
    ),
    tolerance = 1e-15
  )
  expect_equal(
    scenario$freqs$g_ctrl,
    paweh:::.cc_apply_genotype_misclass(
      scenario$freqs$g_ctrl_base, scenario$freqs$M_ctrl
    ),
    tolerance = 1e-15
  )
  expect_identical(
    scenario$modifier$model,
    "diff3p_homhet_homhom"
  )
})

test_that("identity-valued requested modifiers remain explicit scenarios", {
  identity_modifiers <- list(
    phenotype = list(pheno_misclass = TRUE, theta = 0, phi = 0),
    heterogeneity = list(locus_het = TRUE, pi = 1),
    genotype = list(geno_misclass = "1p", e = 0)
  )
  scenario_names <- c(
    phenotype = "phenotype_misclassification",
    heterogeneity = "heterogeneity",
    genotype = "genotype_misclassification"
  )
  for (modifier in names(identity_modifiers)) {
    out <- cc_scenario_call("power", modifiers = identity_modifiers[[modifier]])
    scenario <- scenario_names[[modifier]]
    expect_true(scenario %in% names(out$scenarios))
    expect_equal(
      out$scenarios[[scenario]]$tests,
      out$scenarios$no_error$tests,
      tolerance = 1e-14
    )
  }
})

test_that("inactive ordinary CC modifiers do not create scenarios", {
  out <- cc_scenario_call(
    "power",
    modifiers = list(
      pheno_misclass = FALSE, locus_het = FALSE, geno_misclass = "none"
    )
  )
  expect_named(out$scenarios, "no_error")
})

test_that("scenario-specific MSSNs attain target power", {
  modifiers <- cc_scenario_combine()
  mssn <- cc_scenario_call("mssn", modifiers = modifiers)
  for (scenario in names(mssn$scenarios)) {
    for (test in c("genotypes", "trend")) {
      n_case <- mssn$scenarios[[scenario]]$tests[[test]]$MSSN_case
      args <- c(
        list(N_case = n_case, alpha = 0.05),
        cc_scenario_base(), modifiers
      )
      power <- do.call(cc_power, args)
      expect_gte(power$scenarios[[scenario]]$tests[[test]]$power, 0.80)
      if (n_case > 1) {
        args$N_case <- n_case - 1
        lower <- do.call(cc_power, args)
        expect_lt(lower$scenarios[[scenario]]$tests[[test]]$power, 0.80)
      }
    }
  }
})

test_that("ordinary CC console output names independent scenarios", {
  messages <- capture_messages(
    cc_power(
      N_case = 500, alpha = 0.05,
      input_mode = "model_free", g1 = cc_scenario_g_case,
      g0 = cc_scenario_g_ctrl, prev = 0.05,
      locus_het = TRUE, pi = 0.8,
      pheno_misclass = TRUE, theta = 0.05, phi = 0.01,
      geno_misclass = "1p", e = 0.02
    )
  )
  messages <- paste(messages, collapse = "\n")
  expect_match(messages, "No error", fixed = TRUE)
  expect_match(messages, "Phenotype misclassification", fixed = TRUE)
  expect_match(messages, "Locus heterogeneity", fixed = TRUE)
  expect_match(messages, "Genotype misclassification", fixed = TRUE)
  expect_false(grepl("Adjusted Design|Combined", messages))
})
