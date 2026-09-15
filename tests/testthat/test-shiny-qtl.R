qtl_dashboard_calculation <- function(...) {
  values <- modifyList(paweh:::.paweh_qtl_defaults(), list(...))
  snapshot <- paweh:::.paweh_qtl_snapshot(values)
  paweh:::.paweh_qtl_calculate(snapshot)
}

test_that("QTL UI is a simple input and canonical-output workspace", {
  html <- paste(as.character(paweh:::.paweh_qtl_ui("qtl")), collapse = "\n")
  expect_match(html, "Full continuous trait", fixed = TRUE)
  expect_match(html, "Extreme phenotype sampling", fixed = TRUE)
  expect_match(html, "Multiple quantitative traits", fixed = TRUE)
  expect_match(html, "Estimate power", fixed = TRUE)
  expect_match(html, "Minimum sample size", fixed = TRUE)
  expect_match(html, "Advanced assumptions", fixed = TRUE)
  expect_match(html, ">Output<", fixed = TRUE)
  expect_match(html, "qtl-verbose_output", fixed = TRUE)
  for (removed in c(
    "Sensitivity", "Visualize", "Methods", "Interpretation",
    "Advanced visualization", "Your calculated design"
  )) {
    expect_false(grepl(removed, html, fixed = TRUE), info = removed)
  }
})

test_that("QTL subtype controls retain every required input", {
  shiny::testServer(paweh:::.paweh_qtl_server, {
    session$setInputs(subtype = "continuous", objective = "power")
    session$flushReact()
    core <- paste(as.character(output$core_inputs), collapse = "\n")
    advanced <- paste(as.character(output$advanced_inputs), collapse = "\n")
    expect_match(core, "Variance explained by the QTL", fixed = TRUE)
    expect_match(core, "Dominance parameter", fixed = TRUE)
    expect_match(advanced, "Genotype-count method", fixed = TRUE)

    session$setInputs(subtype = "extreme")
    session$flushReact()
    core <- paste(as.character(output$core_inputs), collapse = "\n")
    advanced <- paste(as.character(output$advanced_inputs), collapse = "\n")
    expect_match(core, "Upper population tail selected", fixed = TRUE)
    expect_match(core, "Lower population tail selected", fixed = TRUE)
    expect_match(advanced, "Lower-tail selected per upper-tail", fixed = TRUE)

    session$setInputs(
      subtype = "multivariate", n_traits = "3",
      mv_test = "threshold_chisq"
    )
    session$flushReact()
    core <- paste(as.character(output$core_inputs), collapse = "\n")
    expect_match(core, "Joint continuous-trait test", fixed = TRUE)
    expect_match(core, "Joint extreme-selection test", fixed = TRUE)
    expect_match(core, "Trait 3", fixed = TRUE)
    expect_match(core, "Correlation: Trait 2 and Trait 3", fixed = TRUE)
    expect_match(core, "Upper population tail selected", fixed = TRUE)
  })
})

test_that("all six QTL workflows call the matching canonical function once", {
  scenarios <- list(
    list("continuous", "power", qtl_anova_power, "power", "qtl_anova_power"),
    list("continuous", "mssn", qtl_anova_mssn, "N", "qtl_anova_mssn"),
    list("extreme", "power", qtl_threshold_chisq_power, "power", "qtl_threshold_chisq_power"),
    list("extreme", "mssn", qtl_threshold_chisq_mssn, "N_total", "qtl_threshold_chisq_mssn"),
    list("multivariate", "power", qtl_multivariate_power_full, "power", "qtl_multivariate_power_full"),
    list("multivariate", "mssn", qtl_multivariate_mssn_full, "N", "qtl_multivariate_mssn_full")
  )
  for (scenario in scenarios) {
    calculation <- qtl_dashboard_calculation(
      subtype = scenario[[1]], objective = scenario[[2]]
    )
    snapshot <- calculation$snapshot
    direct <- do.call(
      scenario[[3]],
      paweh:::.paweh_qtl_call_args(snapshot, verbose = FALSE)
    )
    expect_equal(
      calculation$result[[scenario[[4]]]], direct[[scenario[[4]]]],
      info = paste(scenario[[1]], scenario[[2]])
    )
    expect_identical(
      paweh:::.paweh_qtl_function(snapshot),
      scenario[[5]]
    )
  }
})

test_that("QTL output is exactly the canonical verbose console output", {
  designs <- list(
    continuous = list(subtype = "continuous"),
    extreme = list(subtype = "extreme"),
    multivariate = list(subtype = "multivariate")
  )
  headings <- c(
    continuous = "Falconer Quantitative Trait",
    extreme = "Falconer Threshold-Selected Trait",
    multivariate = "Falconer Multivariate Quantitative Traits"
  )
  for (name in names(designs)) {
    calculation <- do.call(qtl_dashboard_calculation, designs[[name]])
    snapshot <- calculation$snapshot
    expected <- capture.output(
      invisible(do.call(
        get(
          paweh:::.paweh_qtl_function(snapshot),
          envir = asNamespace("paweh")
        ),
        paweh:::.paweh_qtl_call_args(snapshot, verbose = TRUE)
      )),
      type = "message"
    )
    expect_identical(calculation$verbose_output, paste(expected, collapse = "\n"))
    expect_match(calculation$verbose_output, headings[[name]], fixed = TRUE)
    expect_match(calculation$verbose_output, "Test Results|Required Sample Size")
  }
})

test_that("joint threshold QTL power and MSSN remain canonical", {
  for (objective in c("power", "mssn")) {
    calculation <- qtl_dashboard_calculation(
      subtype = "multivariate", objective = objective,
      mv_test = "threshold_chisq"
    )
    direct <- paweh:::.paweh_qtl_call(calculation$snapshot)
    key <- if (objective == "power") "power" else "N_total"
    expect_equal(calculation$result[[key]], direct[[key]])
    expect_identical(calculation$result$test, "threshold_chisq")
    expect_match(calculation$verbose_output, "Threshold-Selected", fixed = TRUE)
  }
})

test_that("QTL validation and canonical errors are displayed without crashing", {
  values <- paweh:::.paweh_qtl_defaults()
  values$alpha <- 1
  expect_error(paweh:::.paweh_qtl_snapshot(values), "Significance level")

  shiny::testServer(paweh:::.paweh_qtl_server, {
    session$setInputs(
      subtype = "extreme", objective = "power", N_case = 150,
      alpha = .05, pd = 30, qtl_var = 10, tau = 0,
      x_upper = 60, x_lower = 50, k = 1
    )
    session$setInputs(calculate = 1)
    expect_null(session$returned$calculation())
    expect_match(session$returned$error(), "middle is excluded", fixed = TRUE)
    expect_match(output$verbose_output, "could not be calculated", fixed = TRUE)
    expect_match(output$verbose_output, "middle is excluded", fixed = TRUE)
  })
})

test_that("QTL calculated output remains frozen until recalculation", {
  shiny::testServer(paweh:::.paweh_qtl_server, {
    session$setInputs(
      subtype = "continuous", objective = "power", N = 500,
      alpha = .05, pd = 30, qtl_var = 10, tau = 0,
      count_method = "rounded"
    )
    session$setInputs(calculate = 1)
    first <- session$returned$calculation()$result$power
    first_output <- session$returned$verbose_output()
    session$setInputs(alpha = .01)
    session$flushReact()
    expect_true(session$returned$changed())
    expect_equal(session$returned$calculation()$result$power, first)
    expect_identical(session$returned$verbose_output(), first_output)
    session$setInputs(calculate = 2)
    expect_false(session$returned$changed())
    expect_false(identical(session$returned$calculation()$result$power, first))
  })
})

test_that("multivariate correlation inputs form the canonical matrix", {
  valid <- paweh:::.paweh_qtl_snapshot(modifyList(
    paweh:::.paweh_qtl_defaults(),
    list(
      subtype = "multivariate", n_traits = 3,
      corr_1_2 = .2, corr_1_3 = .1, corr_2_3 = -.1
    )
  ))
  expect_equal(unname(diag(valid$backend_args$cor_matrix)), rep(1, 3))
  expect_true(isSymmetric(valid$backend_args$cor_matrix))

  invalid <- modifyList(paweh:::.paweh_qtl_defaults(), list(
    subtype = "multivariate", n_traits = 3,
    corr_1_2 = .9, corr_1_3 = .9, corr_2_3 = -.9
  ))
  snapshot <- paweh:::.paweh_qtl_snapshot(invalid)
  expect_error(paweh:::.paweh_qtl_calculate(snapshot), "positive definite")
})
