qtl_dashboard_calculation <- function(...) {
  values <- modifyList(paweh:::.paweh_qtl_defaults(), list(...))
  paweh:::.paweh_qtl_calculate(paweh:::.paweh_qtl_snapshot(values))
}

test_that("QTL UI exposes all production workflows", {
  html <- paste(as.character(paweh:::.paweh_qtl_ui("qtl")), collapse = "\n")
  expect_match(html, "How will the phenotype be analyzed?", fixed = TRUE)
  expect_match(html, "Full continuous trait", fixed = TRUE)
  expect_match(html, "Extreme phenotype sampling", fixed = TRUE)
  expect_match(html, "Multiple quantitative traits", fixed = TRUE)
  expect_match(html, "Estimate power", fixed = TRUE)
  expect_match(html, "Minimum sample size", fixed = TRUE)
  expect_match(html, "Advanced assumptions", fixed = TRUE)
  for (label in c("Results", "Sensitivity", "Visualize", "Methods")) {
    expect_match(html, label, fixed = TRUE)
  }
  expect_match(html, "Advanced visualization", fixed = TRUE)
  expect_false(grepl("<details[^>]* open", html))
})

test_that("QTL subtype controls use human-readable canonical semantics", {
  shiny::testServer(paweh:::.paweh_qtl_server, {
    session$setInputs(subtype = "continuous", objective = "power")
    session$flushReact()
    core <- paste(as.character(output$core_inputs), collapse = "\n")
    objective <- paste(as.character(output$objective_inputs), collapse = "\n")
    expect_match(core, "Modeled-allele frequency", fixed = TRUE)
    expect_match(core, "Variance explained by the QTL", fixed = TRUE)
    expect_match(core, "Dominance-to-additivity ratio", fixed = TRUE)
    expect_match(objective, "Total sample size", fixed = TRUE)

    session$setInputs(subtype = "extreme")
    session$flushReact()
    core <- paste(as.character(output$core_inputs), collapse = "\n")
    expect_match(core, "Upper population tail selected", fixed = TRUE)
    expect_match(core, "standardized-normal population percentiles", fixed = TRUE)

    session$setInputs(subtype = "multivariate", n_traits = "3", mv_test = "threshold_chisq")
    session$flushReact()
    core <- paste(as.character(output$core_inputs), collapse = "\n")
    expect_match(core, "Joint continuous-trait test", fixed = TRUE)
    expect_match(core, "Joint extreme-selection test", fixed = TRUE)
    expect_match(core, "Trait 3", fixed = TRUE)
    expect_match(core, "Correlation: Trait 2 and Trait 3", fixed = TRUE)
    expect_match(core, "display a selected trait pair", fixed = TRUE)
  })
})

test_that("all six QTL interfaces exactly reproduce canonical calls", {
  scenarios <- list(
    list(subtype = "continuous", objective = "power", fun = qtl_anova_power, key = "power"),
    list(subtype = "continuous", objective = "mssn", fun = qtl_anova_mssn, key = "N"),
    list(subtype = "extreme", objective = "power", fun = qtl_threshold_chisq_power, key = "power"),
    list(subtype = "extreme", objective = "mssn", fun = qtl_threshold_chisq_mssn, key = "N_total"),
    list(subtype = "multivariate", objective = "power", fun = qtl_multivariate_power_full, key = "power"),
    list(subtype = "multivariate", objective = "mssn", fun = qtl_multivariate_mssn_full, key = "N")
  )
  for (scenario in scenarios) {
    calculation <- qtl_dashboard_calculation(
      subtype = scenario$subtype, objective = scenario$objective
    )
    args <- paweh:::.paweh_qtl_repro_args(calculation)
    direct <- do.call(scenario$fun, args)
    expect_equal(calculation$result[[scenario$key]], direct[[scenario$key]])
    expect_false(any(c("subtype", "objective", "display") %in% names(args)))
    expect_match(
      paweh:::.paweh_call_text(paweh:::.paweh_qtl_repro_call(calculation)),
      paweh:::.paweh_qtl_function(calculation$snapshot), fixed = TRUE
    )
  }
})

test_that("joint threshold QTL calls also remain canonical", {
  for (objective in c("power", "mssn")) {
    calculation <- qtl_dashboard_calculation(
      subtype = "multivariate", objective = objective,
      mv_test = "threshold_chisq"
    )
    fun <- if (objective == "power") qtl_multivariate_power_full else qtl_multivariate_mssn_full
    direct <- do.call(fun, paweh:::.paweh_qtl_repro_args(calculation))
    key <- if (objective == "power") "power" else "N_total"
    expect_equal(calculation$result[[key]], direct[[key]])
    expect_identical(calculation$result$test, "threshold_chisq")
  }
})

test_that("QTL validation is friendly and correlation matrices are structured", {
  values <- paweh:::.paweh_qtl_defaults()
  values$alpha <- 1
  expect_error(paweh:::.paweh_qtl_snapshot(values), "Significance level")
  values <- paweh:::.paweh_qtl_defaults()
  values$pd <- 0
  expect_error(paweh:::.paweh_qtl_snapshot(values), "Modeled-allele frequency")
  values <- paweh:::.paweh_qtl_defaults()
  values$subtype <- "extreme"
  values$x_upper <- 60
  values$x_lower <- 50
  expect_error(paweh:::.paweh_qtl_snapshot(values), "middle is excluded")
  values <- paweh:::.paweh_qtl_defaults()
  values$subtype <- "multivariate"
  values$n_traits <- 3
  values$corr_1_2 <- values$corr_1_3 <- .9
  values$corr_2_3 <- -.9
  expect_error(paweh:::.paweh_qtl_snapshot(values), "positive definite")

  valid <- paweh:::.paweh_qtl_snapshot(modifyList(
    paweh:::.paweh_qtl_defaults(),
    list(subtype = "multivariate", n_traits = 3, corr_1_2 = .2, corr_1_3 = .1, corr_2_3 = -.1)
  ))
  expect_equal(unname(diag(valid$backend_args$cor_matrix)), rep(1, 3))
  expect_true(isSymmetric(valid$backend_args$cor_matrix))
})

test_that("QTL calculated state is frozen until recalculation", {
  shiny::testServer(paweh:::.paweh_qtl_server, {
    session$setInputs(
      subtype = "continuous", objective = "power", N = 500,
      alpha = .05, pd = 30, qtl_var = 10, tau = 0,
      count_method = "rounded"
    )
    session$setInputs(calculate = 1)
    first <- session$returned$calculation()$result$power
    session$setInputs(alpha = .01)
    session$flushReact()
    expect_true(session$returned$changed())
    expect_equal(session$returned$calculation()$result$power, first)
    expect_match(paste(as.character(output$changed_notice), collapse = "\n"),
      "Inputs have changed", fixed = TRUE)
    session$setInputs(calculate = 2)
    expect_false(session$returned$changed())
    expect_false(identical(session$returned$calculation()$result$power, first))
  })
})

test_that("QTL results and advanced details use canonical returned quantities", {
  calculations <- list(
    qtl_dashboard_calculation(subtype = "continuous"),
    qtl_dashboard_calculation(subtype = "extreme"),
    qtl_dashboard_calculation(subtype = "multivariate"),
    qtl_dashboard_calculation(subtype = "multivariate", mv_test = "threshold_chisq")
  )
  for (calculation in calculations) {
    html <- paste(as.character(paweh:::.paweh_qtl_results_ui(calculation)), collapse = "\n")
    expect_match(html, "Interpretation", fixed = TRUE)
    expect_match(html, "Advanced calculation details", fixed = TRUE)
    expect_match(html, "Model specification", fixed = TRUE)
    expect_match(html, "Genotype and model quantities", fixed = TRUE)
    expect_match(html, "Test-specific quantities", fixed = TRUE)
    expect_match(html, "Reproduce in R", fixed = TRUE)
    expect_match(html, paweh:::.paweh_qtl_function(calculation$snapshot), fixed = TRUE)
    expect_false(grepl("<details[^>]* open", html))
  }
  extreme <- paste(as.character(paweh:::.paweh_qtl_advanced_ui(calculations[[2]])), collapse = "\n")
  expect_match(extreme, "Selection quantities", fixed = TRUE)
  expect_match(extreme, "Expected upper-tail proportion", fixed = TRUE)
  multivariate <- paste(as.character(paweh:::.paweh_qtl_advanced_ui(calculations[[3]])), collapse = "\n")
  expect_match(multivariate, "Multivariate model quantities", fixed = TRUE)
  expect_match(multivariate, "Residual covariance", fixed = TRUE)
})

test_that("QTL sensitivity points are separate canonical calls", {
  calculations <- list(
    qtl_dashboard_calculation(subtype = "continuous"),
    qtl_dashboard_calculation(subtype = "extreme"),
    qtl_dashboard_calculation(subtype = "multivariate"),
    qtl_dashboard_calculation(subtype = "multivariate", mv_test = "threshold_chisq")
  )
  for (calculation in calculations) {
    key <- names(paweh:::.paweh_qtl_specs(calculation))[1]
    spec <- paweh:::.paweh_qtl_specs(calculation)[[key]]
    sensitivity <- paweh:::.paweh_qtl_sensitivity(
      calculation, key, c(spec$value * .8, spec$value * 1.2), n = 5
    )
    expect_equal(nrow(sensitivity$data), 5)
    expect_equal(sensitivity$baseline_x, spec$value)
    expect_equal(sensitivity$data$y[3], calculation$result$power)
    expect_s3_class(paweh:::.paweh_qtl_sensitivity_plot(sensitivity), "ggplot")
  }
})

test_that("QTL two-dimensional visualizations use canonical plot helpers", {
  continuous <- qtl_dashboard_calculation(subtype = "continuous")
  extreme <- qtl_dashboard_calculation(subtype = "extreme")
  multivariate <- qtl_dashboard_calculation(subtype = "multivariate")
  threshold <- qtl_dashboard_calculation(subtype = "multivariate", mv_test = "threshold_chisq")
  expect_s3_class(paweh:::.paweh_qtl_single_plot(continuous), "ggplot")
  extreme_plot <- paweh:::.paweh_qtl_single_plot(extreme)
  expect_s3_class(extreme_plot, "ggplot")
  expect_true(any(vapply(extreme_plot$layers, function(layer) inherits(layer$geom, "GeomVline"), logical(1))))
  expect_s3_class(paweh:::.paweh_qtl_multivariate_plot(multivariate, "genotype"), "ggplot")
  expect_s3_class(paweh:::.paweh_qtl_multivariate_plot(multivariate, "mixture"), "ggplot")
  expect_s3_class(paweh:::.paweh_qtl_multivariate_plot(threshold, "selection"), "ggplot")
})

test_that("three- and four-trait designs visualize a frozen selected pair", {
  calculation <- qtl_dashboard_calculation(
    subtype = "multivariate", n_traits = 4,
    corr_1_2 = .1, corr_1_3 = .05, corr_1_4 = 0,
    corr_2_3 = .1, corr_2_4 = .05, corr_3_4 = .1
  )
  args <- paweh:::.paweh_qtl_multivariate_plot_args(
    calculation, "genotype_density", grid_n = 20, trait_pair = c(2, 4)
  )
  expect_equal(args$qtl_var, calculation$snapshot$backend_args$qtl_var[c(2, 4)])
  expect_equal(args$tau, calculation$snapshot$backend_args$tau[c(2, 4)])
  expect_equal(args$cor_matrix, calculation$snapshot$backend_args$cor_matrix[c(2, 4), c(2, 4)])
  expect_s3_class(paweh:::.paweh_qtl_multivariate_plot(calculation, "genotype", c(2, 4)), "ggplot")

  shiny::testServer(paweh:::.paweh_qtl_server, {
    session$setInputs(
      subtype = "multivariate", objective = "power", mv_test = "pillai",
      n_traits = "4", N = 500, alpha = .05, pd = 30,
      mv_qtl_var_1 = 10, mv_qtl_var_2 = 5, mv_qtl_var_3 = 3, mv_qtl_var_4 = 2,
      mv_tau_1 = 0, mv_tau_2 = .5, mv_tau_3 = 0, mv_tau_4 = 0,
      corr_1_2 = .1, corr_1_3 = .05, corr_1_4 = 0,
      corr_2_3 = .1, corr_2_4 = .05, corr_3_4 = .1
    )
    session$setInputs(calculate = 1)
    controls <- paste(as.character(output$visualization_controls), collapse = "\n")
    expect_match(controls, "Horizontal trait", fixed = TRUE)
    expect_match(controls, "Vertical trait", fixed = TRUE)
    expect_null(session$returned$surface())
  })
})

test_that("QTL 3D surfaces are restrained, on-demand, and scientifically distinct", {
  skip_if_not_installed("plotly")
  calculation <- qtl_dashboard_calculation(subtype = "multivariate")
  genotype <- plotly::plotly_build(paweh:::.paweh_qtl_multivariate_surface(calculation, "genotype"))
  mixture <- plotly::plotly_build(paweh:::.paweh_qtl_multivariate_surface(calculation, "mixture"))
  genotype_names <- vapply(genotype$x$data, function(trace) trace$name %||% "", character(1))
  mixture_names <- vapply(mixture$x$data, function(trace) trace$name %||% "", character(1))
  expect_identical(genotype_names[1:3], paste("Genotype", 0:2))
  expect_identical(mixture_names[1], "Mixture density")
  expect_false(any(grepl("Genotype [0-2]", mixture_names)))
  expect_lte(as.numeric(object.size(genotype)), 2e6)

  shiny::testServer(paweh:::.paweh_qtl_server, {
    session$setInputs(
      subtype = "multivariate", objective = "power", mv_test = "pillai",
      n_traits = "2", N = 500, alpha = .05, pd = 30,
      mv_qtl_var_1 = 10, mv_qtl_var_2 = 5,
      mv_tau_1 = 0, mv_tau_2 = .5, corr_1_2 = 0
    )
    session$setInputs(calculate = 1)
    expect_null(session$returned$surface())
    session$setInputs(alpha = .01, surface_mode = "genotype", generate_surface = 1)
    expect_s3_class(session$returned$surface(), "plotly")
    expect_identical(session$returned$calculation()$snapshot$backend_args$alpha, .05)
  })
})

test_that("Gordon multivariate anchor is preserved through the dashboard", {
  calculation <- qtl_dashboard_calculation(
    subtype = "multivariate", objective = "mssn", alpha = 5e-8,
    pd = 5, mv_qtl_var_1 = 10, mv_qtl_var_2 = 5,
    mv_tau_1 = 0, mv_tau_2 = .5, corr_1_2 = 0
  )
  expect_identical(calculation$result$N, 326L)
  expect_equal(calculation$result$historical_fractional_mssn, 325.5056664947, tolerance = 1e-9)
  direct <- do.call(qtl_multivariate_mssn_full, paweh:::.paweh_qtl_repro_args(calculation))
  expect_identical(direct$N, calculation$result$N)
  expect_equal(direct$historical_fractional_mssn, calculation$result$historical_fractional_mssn)
})

test_that("QTL plots use the shared restrained genotype palette", {
  expect_identical(
    unname(paweh:::.paweh_qtl_genotype_colors()),
    c("#3F4850", "#6F879A", "#355C7D")
  )
  code <- paste(vapply(
    list(
      plot_qtl_genotype_distribution,
      plot_qtl_multivariate_contour,
      plot_qtl_multivariate_surface3d
    ), function(fun) paste(deparse(body(fun)), collapse = "\n"), character(1)
  ), collapse = "\n")
  expect_match(code, ".paweh_qtl_genotype_colors()", fixed = TRUE)
  expect_false(grepl("viridis|rainbow|#00FFFF|cyan|coral|#FF0000", code, ignore.case = TRUE))
})
