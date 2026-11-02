se_case_control_fixture <- function() {
  genotypes <- rbind(
    rs1 = c(0, 0, 1, 0, 1, 1, 2, 1, 2, 2),
    rs2 = c(2, 1, 1, 1, 0, 1, 1, 0, 2, 0)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(GT = genotypes),
    colData = S4Vectors::DataFrame(status = rep(c("case", "control"), 5))
  )
}

test_that("the shared helper separates samples and orders all three frequencies", {
  se <- se_case_control_fixture()
  expected <- list(g1 = c(0.2, 0.4, 0.4), g0 = c(0.4, 0.4, 0.2))
  expect_equal(
    .cc_se_genotype_frequencies(se, "rs1", "GT", "status", "case", "control"),
    expected
  )
  expect_equal(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "control", "case"),
    list(g1 = expected$g0, g0 = expected$g1)
  )
  expect_equal(
    .cc_se_genotype_frequencies(se, "rs2", "GT", "status", "case", "control"),
    list(g1 = c(0.2, 0.4, 0.4), g0 = c(0.4, 0.6, 0))
  )

  SummarizedExperiment::colData(se)$status <- rep(c(2, 1), 5)
  expect_equal(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", 2, 1), expected
  )
  SummarizedExperiment::colData(se)$status <- factor(
    rep(c("case", "control"), 5), levels = c("control", "case")
  )
  expect_equal(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "case", "control"),
    expected
  )
})

test_that("both wrappers return the original results for the prospective design", {
  se <- se_case_control_fixture()
  for (ranged in c(FALSE, TRUE)) {
    input <- if (ranged) methods::as(se, "RangedSummarizedExperiment") else se
    if (ranged) expect_s4_class(input, "RangedSummarizedExperiment")
    common <- list(
      se = input, variant = "rs1", assay_name = "GT", group_col = "status",
      case_value = "case", control_value = "control", alpha = 0.05,
      k = 2, w = c(0, 1, 3), verbose = FALSE
    )
    direct <- list(
      input_mode = "model_free", g1 = c(0.2, 0.4, 0.4),
      g0 = c(0.4, 0.4, 0.2), alpha = 0.05, k = 2,
      w = c(0, 1, 3), verbose = FALSE
    )
    mssn <- withVisible(do.call(cc_mssn_from_se, c(common, list(power = 0.8))))
    power <- withVisible(do.call(cc_power_from_se, c(common, list(N_case = 123))))
    expect_identical(mssn$value, do.call(cc_mssn, c(direct, list(power = 0.8))))
    expect_identical(power$value, do.call(cc_power, c(direct, list(N_case = 123))))
    expect_s3_class(mssn$value, "cc_mssn")
    expect_s3_class(power$value, "cc_power")
    expect_false(mssn$visible)
    expect_false(power$visible)
    expect_equal(power$value$N_case, 123)
    expect_equal(power$value$N_ctrl, 246)
  }
})

test_that("missing genotypes are excluded separately and other groups are ignored", {
  se <- se_case_control_fixture()
  SummarizedExperiment::assay(se, "GT")[1, c(1, 2)] <- c(NA, NaN)
  expect_equal(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "case", "control"),
    list(g1 = c(0, 0.5, 0.5), g0 = c(0.25, 0.5, 0.25))
  )
  expect_identical(
    cc_power_from_se(se, 1, "GT", "status", "case", "control",
                     N_case = 100, alpha = 0.05, verbose = FALSE),
    cc_power(N_case = 100, alpha = 0.05, input_mode = "model_free",
             g1 = c(0, 0.5, 0.5), g0 = c(0.25, 0.5, 0.25), verbose = FALSE)
  )
  SummarizedExperiment::colData(se)$status[c(1, 2)] <- c(NA, "other")
  SummarizedExperiment::assay(se, "GT")[1, c(1, 2)] <- c(99, 0.3)
  expect_equal(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "case", "control"),
    list(g1 = c(0, 0.5, 0.5), g0 = c(0.25, 0.5, 0.25))
  )
  for (group in c("case", "control")) {
    empty <- se_case_control_fixture()
    selected <- SummarizedExperiment::colData(empty)$status == group
    SummarizedExperiment::assay(empty, "GT")[1, selected] <- NA_real_
    expect_error(
      .cc_se_genotype_frequencies(empty, 1, "GT", "status", "case", "control"),
      "Each group must retain at least one"
    )
  }
})

test_that("extraction rejects invalid and ambiguous selectors", {
  se <- se_case_control_fixture()
  args <- list(se, "rs1", "GT", "status", "case", "control")
  expect_error(
    .cc_se_genotype_frequencies(matrix(0, 2, 2), 1, "GT", "status", 1, 0),
    "se must be a SummarizedExperiment"
  )
  for (variant in list("absent", "", NA_character_, 0, -1, 1.5, 3, Inf,
                       NA_real_, c(1, 2), TRUE, NULL)) {
    bad <- args
    bad[2] <- list(variant)
    expect_error(do.call(.cc_se_genotype_frequencies, bad), "variant")
  }
  for (selector in list("absent", "", NA_character_, c("GT", "GT"), 1, NULL)) {
    bad <- args
    bad[3] <- list(selector)
    expect_error(do.call(.cc_se_genotype_frequencies, bad), "assay_name")
    bad <- args
    bad[4] <- list(selector)
    expect_error(do.call(.cc_se_genotype_frequencies, bad), "group_col")
  }
  rownames(se) <- c("rs1", "rs1")
  expect_error(
    .cc_se_genotype_frequencies(se, "rs1", "GT", "status", "case", "control"),
    "exactly one row name"
  )
  expect_no_error(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "case", "control")
  )
  SummarizedExperiment::assays(se) <- list(
    GT = SummarizedExperiment::assay(se), GT = SummarizedExperiment::assay(se)
  )
  expect_error(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "case", "control"),
    "exactly one assay"
  )
  se <- se_case_control_fixture()
  colnames(SummarizedExperiment::colData(se)) <- "duplicated"
  SummarizedExperiment::colData(se) <- cbind(
    SummarizedExperiment::colData(se), SummarizedExperiment::colData(se)
  )
  expect_error(
    .cc_se_genotype_frequencies(se, 1, "GT", "duplicated", "case", "control"),
    "exactly one column"
  )
})

test_that("invalid group labels and genotype calls fail clearly", {
  se <- se_case_control_fixture()
  for (label in list(NA, NA_real_, Inf, "", c("case", "control"), TRUE, NULL)) {
    expect_error(
      .cc_se_genotype_frequencies(se, 1, "GT", "status", label, "control"),
      "single non-missing numeric or character"
    )
    expect_error(
      .cc_se_genotype_frequencies(se, 1, "GT", "status", "case", label),
      "single non-missing numeric or character"
    )
  }
  expect_error(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "case", "case"),
    "must be different"
  )
  expect_error(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", "absent", "control"),
    "must identify samples"
  )
  SummarizedExperiment::colData(se)$status <- rep(TRUE, ncol(se))
  expect_error(
    .cc_se_genotype_frequencies(se, 1, "GT", "status", 1, 0), "grouping column"
  )
  for (genotype in list(-1, 3, 0.5, Inf, -Inf, "1", TRUE)) {
    bad <- se_case_control_fixture()
    SummarizedExperiment::assay(bad, "GT")[1, 1] <- genotype
    if (is.logical(genotype)) {
      SummarizedExperiment::assay(bad, "GT") <-
        SummarizedExperiment::assay(bad, "GT") > 0
    }
    expect_error(
      .cc_se_genotype_frequencies(bad, 1, "GT", "status", "case", "control"),
      "numeric hard calls"
    )
  }
})

test_that("modifiers are forwarded and conflicting or misspelled dots are rejected", {
  common <- list(
    se = se_case_control_fixture(), variant = 1, assay_name = "GT",
    group_col = "status", case_value = "case", control_value = "control",
    alpha = 0.05, k = 1, w = c(0, 1, 2), verbose = FALSE
  )
  modifiers <- list(
    locus_het = TRUE, pi = 0.8, geno_misclass = "1p", e = 0.01,
    pheno_misclass = TRUE, prev = 0.1, theta = 0.05, phi = 0.01
  )
  direct <- list(
    input_mode = "model_free", g1 = c(0.2, 0.4, 0.4),
    g0 = c(0.4, 0.4, 0.2), alpha = 0.05, verbose = FALSE
  )
  for (kind in c("mssn", "power")) {
    wrapper <- get(paste0("cc_", kind, "_from_se"))
    original <- get(paste0("cc_", kind))
    design <- if (kind == "mssn") list(power = 0.8) else list(N_case = 100)
    expect_identical(
      do.call(wrapper, c(common, design, modifiers)),
      do.call(original, c(direct, design, modifiers))
    )
    invalid_dots <- list(
      list(g1 = c(0.1, 0.3, 0.6)), list(g0 = c(0.1, 0.3, 0.6)),
      list(input_mode = "model_based"), list(inp = "model_based"),
      list(g = 1), list(pd = 0.3), list(R2 = 2), list(MOI = "D"),
      list(locus_he = TRUE), list(unknown = 1), list(0.5),
      list(pi = 0.8, pi = 0.9)
    )
    for (dots in invalid_dots) {
      expect_error(do.call(wrapper, c(common, design, dots)), "exact")
    }
  }
  expect_error(
    do.call(cc_power_from_se, c(common, list(N_case = -1))), "N_case"
  )
  expect_error(
    do.call(cc_mssn_from_se, c(common, list(power = 1))), "power"
  )
})

test_that("wrappers preserve the formatted reports and invisible returns", {
  common <- list(
    se = se_case_control_fixture(), variant = 1, assay_name = "GT",
    group_col = "status", case_value = "case", control_value = "control",
    alpha = 0.05
  )
  direct <- list(
    input_mode = "model_free", g1 = c(0.2, 0.4, 0.4),
    g0 = c(0.4, 0.4, 0.2), alpha = 0.05
  )
  for (kind in c("mssn", "power")) {
    wrapper <- get(paste0("cc_", kind, "_from_se"))
    original <- get(paste0("cc_", kind))
    design <- if (kind == "mssn") list(power = 0.8) else list(N_case = 100)
    report <- capture.output(
      result <- withVisible(do.call(wrapper, c(common, design))), type = "message"
    )
    expected <- capture.output(do.call(original, c(direct, design)), type = "message")
    expect_identical(report, expected)
    expect_equal(sum(grepl("PAWEH Case-Control Study", report, fixed = TRUE)), 1)
    expect_true(any(grepl("Genotype", report, fixed = TRUE)))
    expect_true(any(grepl("Trend", report, fixed = TRUE)))
    expect_false(result$visible)
    expect_output(suppressMessages(do.call(wrapper, c(common, design))), NA)
    expect_message(do.call(wrapper, c(common, design, list(verbose = FALSE))), NA)
    expect_output(do.call(wrapper, c(common, design, list(verbose = FALSE))), NA)
  }
})
