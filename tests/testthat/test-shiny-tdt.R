test_that("TDT UI exposes the complete production workflow", {
  html <- paste(as.character(pawh:::.pawh_tdt_ui("tdt")), collapse = "\n")
  expect_match(html, "Estimate power", fixed = TRUE)
  expect_match(html, "Minimum sample size", fixed = TRUE)
  expect_match(html, "Direct transmission quantities", fixed = TRUE)
  expect_match(html, "Advanced assumptions", fixed = TRUE)
  expect_match(html, "Phenotype misclassification", fixed = TRUE)
  expect_match(html, "Locus heterogeneity", fixed = TRUE)
  expect_match(html, "Sensitivity", fixed = TRUE)
  expect_match(html, "Visualize", fixed = TRUE)
  expect_match(html, "Methods", fixed = TRUE)
  expect_match(html, "<details class=\"pawh-sidebar-section\">", fixed = TRUE)
  expect_false(grepl("<details[^>]* open", html))
})

test_that("objective and input controls are mode-specific", {
  shiny::testServer(pawh:::.pawh_tdt_server, {
    session$setInputs(objective = "power", input_mode = "model_based")
    session$flushReact()
    power_ui <- paste(as.character(output$objective_inputs), collapse = "\n")
    model_ui <- paste(as.character(output$model_inputs), collapse = "\n")
    expect_match(power_ui, "Number of affected-child trios", fixed = TRUE)
    expect_false(grepl("Target power", power_ui, fixed = TRUE))
    expect_match(model_ui, "Heterozygote relative risk", fixed = TRUE)
    expect_false(grepl("Expected transmitted count", model_ui, fixed = TRUE))

    session$setInputs(objective = "mssn", input_mode = "model_free")
    session$flushReact()
    mssn_ui <- paste(as.character(output$objective_inputs), collapse = "\n")
    direct_ui <- paste(as.character(output$model_inputs), collapse = "\n")
    expect_match(mssn_ui, "Target power", fixed = TRUE)
    expect_false(grepl("Number of affected-child trios", mssn_ui, fixed = TRUE))
    expect_match(direct_ui, "Expected transmitted count", fixed = TRUE)
    expect_match(direct_ui, "represented by ET/ENT", fixed = TRUE)

    session$setInputs(misclassification = TRUE, heterogeneity = FALSE)
    session$flushReact()
    modifier_ui <- paste(as.character(output$model_inputs), collapse = "\n")
    expect_match(modifier_ui, "prevalence for misclassification", ignore.case = TRUE)
  })
})

test_that("Shiny snapshots reproduce canonical power and MSSN", {
  values <- pawh:::.pawh_tdt_defaults()
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  direct <- tdt_power(
    N = 600, input_mode = "model_based", pd = .3, prev = .05,
    R1 = 1.5, R2 = 2.25, alpha = .05, delta_prime = 1,
    misclass_rate = 0, heter_rate = 0, verbose = FALSE
  )
  expect_equal(calculation$result$power, direct$power)
  expect_equal(calculation$result$gT_star, direct$gT_star)
  expect_identical(pawh:::.pawh_tdt_scenarios(calculation), "no_error")

  values$objective <- "mssn"
  values$input_mode <- "model_free"
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  direct <- tdt_mssn(
    target_power = .8, input_mode = "model_free", ET = 140, ENT = 100,
    n_trios = 120, alpha = .05, misclass_rate = 0, heter_rate = 0,
    verbose = FALSE
  )
  expect_equal(calculation$result$N, direct$N)
  expect_equal(calculation$result$percent_increase, direct$percent_increase)
})

test_that("modifier scenarios remain separate canonical results", {
  scenarios <- list(
    misclassification = c(misclassification = TRUE, misclass_rate = 5),
    heterogeneity = c(heterogeneity = TRUE, heter_rate = 20),
    both = c(misclassification = TRUE, misclass_rate = 5, heterogeneity = TRUE, heter_rate = 20)
  )
  for (name in names(scenarios)) {
    values <- pawh:::.pawh_tdt_defaults()
    for (field in names(scenarios[[name]])) values[[field]] <- scenarios[[name]][[field]]
    calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
    direct <- tdt_power(
      N = 600, input_mode = "model_based", pd = .3, prev = .05,
      R1 = 1.5, R2 = 2.25, alpha = .05, delta_prime = 1,
      misclass_rate = if (isTRUE(values$misclassification)) .05 else 0,
      heter_rate = if (isTRUE(values$heterogeneity)) .2 else 0,
      verbose = FALSE
    )
    expect_equal(calculation$result$power, direct$power, info = name)
  }

  values <- pawh:::.pawh_tdt_defaults()
  values$misclassification <- values$heterogeneity <- TRUE
  values$misclass_rate <- 5
  values$heter_rate <- 20
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  expect_identical(
    pawh:::.pawh_tdt_scenarios(calculation),
    c("no_error", "misclassification", "heterogeneity")
  )
  html <- paste(as.character(pawh:::.pawh_tdt_results_ui(calculation)), collapse = "\n")
  expect_match(html, "No-error design", fixed = TRUE)
  expect_match(html, "Phenotype misclassification", fixed = TRUE)
  expect_match(html, "Locus heterogeneity", fixed = TRUE)
  expect_false(grepl("Combined|Adjusted combined|Misclassification \\+ heterogeneity", html))
  expect_match(html, "separate sensitivity scenarios", fixed = TRUE)
})

test_that("zero-valued modifier toggles do not create active scenarios", {
  values <- pawh:::.pawh_tdt_defaults()
  values$misclassification <- values$heterogeneity <- TRUE
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  expect_false(any(unlist(calculation$active)))
  expect_identical(pawh:::.pawh_tdt_scenarios(calculation), "no_error")
})

test_that("validation is friendly and calculated state remains frozen", {
  values <- pawh:::.pawh_tdt_defaults()
  values$alpha <- 1
  expect_error(pawh:::.pawh_tdt_snapshot(values), "Significance level")
  values <- pawh:::.pawh_tdt_defaults()
  values$input_mode <- "model_free"
  values$ET <- values$ENT <- 0
  expect_error(pawh:::.pawh_tdt_snapshot(values), "cannot both be zero")

  shiny::testServer(pawh:::.pawh_tdt_server, {
    session$setInputs(
      objective = "power", input_mode = "model_based", N = 600,
      alpha = .05, prev = 5, pd = 30, R1 = 1.5, R2 = 2.25,
      delta_prime = 1, misclassification = FALSE, heterogeneity = FALSE
    )
    session$setInputs(calculate = 1)
    first <- session$returned$calculation()$result$power$no_error
    session$setInputs(alpha = .01)
    session$flushReact()
    expect_true(session$returned$changed())
    expect_equal(session$returned$calculation()$result$power$no_error, first)
    session$setInputs(calculate = 2)
    expect_false(session$returned$changed())
    expect_false(identical(session$returned$calculation()$result$power$no_error, first))
  })
})

test_that("null-effect MSSN is displayed without invalid percentages", {
  values <- pawh:::.pawh_tdt_defaults()
  values$objective <- "mssn"
  values$input_mode <- "model_free"
  values$ET <- values$ENT <- 100
  values$misclassification <- values$heterogeneity <- TRUE
  values$misclass_rate <- 5
  values$heter_rate <- 20
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  expect_true(all(is.infinite(unlist(calculation$result$N))))
  html <- paste(as.character(pawh:::.pawh_tdt_results_ui(calculation)), collapse = "\n")
  expect_match(html, "Required affected-child trios", fixed = TRUE)
  expect_match(html, ">Inf<", fixed = TRUE)
  expect_match(html, "not defined", fixed = TRUE)
  expect_false(grepl("NaN%|Inf%", html))
})

test_that("sensitivity points are canonical calls of the frozen design", {
  values <- pawh:::.pawh_tdt_defaults()
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  sensitivity <- pawh:::.pawh_tdt_sensitivity(calculation, "pd", c(.2, .4), n = 5)
  expect_equal(nrow(sensitivity$data), 5)
  expected <- tdt_power(
    N = 600, input_mode = "model_based", pd = .3, prev = .05,
    R1 = 1.5, R2 = 2.25, alpha = .05, delta_prime = 1,
    misclass_rate = 0, heter_rate = 0, verbose = FALSE
  )
  expect_equal(sensitivity$data$y[3], expected$power$no_error)
  expect_equal(sensitivity$baseline_x, .3)

  values$objective <- "mssn"
  values$misclassification <- TRUE
  values$misclass_rate <- 5
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  expect_identical(names(pawh:::.pawh_tdt_specs(calculation))[1], "misclass_rate")
  sensitivity <- pawh:::.pawh_tdt_sensitivity(calculation, "misclass_rate", c(0, .1), n = 3)
  expected <- tdt_mssn(
    target_power = .8, input_mode = "model_based", pd = .3, prev = .05,
    R1 = 1.5, R2 = 2.25, alpha = .05, delta_prime = 1,
    misclass_rate = .05, heter_rate = 0, verbose = FALSE
  )
  middle <- sensitivity$data[sensitivity$data$x == sensitivity$data$x[3], ]
  expect_equal(middle$y, unlist(expected$N[c("no_error", "misclassification")], use.names = FALSE))
  expect_equal(sensitivity$baseline_x, .05)
})

test_that("TDT plots use shared restrained colors and redundant line styles", {
  values <- pawh:::.pawh_tdt_defaults()
  values$misclassification <- values$heterogeneity <- TRUE
  values$misclass_rate <- 5
  values$heter_rate <- 20
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  sensitivity <- pawh:::.pawh_tdt_sensitivity(calculation, "misclass_rate", c(0, .1), n = 5)
  sensitivity_plot <- pawh:::.pawh_tdt_sensitivity_plot(sensitivity)
  transmission_plot <- pawh:::.pawh_tdt_transmission_plot(calculation)
  expect_s3_class(sensitivity_plot, "ggplot")
  expect_s3_class(transmission_plot, "ggplot")
  expect_true(any(vapply(sensitivity_plot$scales$scales,
    function(scale) "colour" %in% scale$aesthetics, logical(1))))
  expect_true(any(vapply(sensitivity_plot$scales$scales,
    function(scale) "linetype" %in% scale$aesthetics, logical(1))))
  expect_true(any(vapply(transmission_plot$scales$scales,
    function(scale) "fill" %in% scale$aesthetics, logical(1))))
  expect_false(grepl("#00FFFF|cyan|coral|#FF0000", paste(deparse(sensitivity_plot), collapse = ""), ignore.case = TRUE))
})

test_that("TDT literature anchors remain unchanged", {
  buyske <- tdt_mssn(
    target_power = .80, alpha = 1e-5, pd = .50, prev = .01,
    R1 = sqrt(2.5), R2 = 2.5, delta_prime = 1,
    misclass_rate = .05, heter_rate = 0, verbose = FALSE
  )
  expect_equal(buyske$N$no_error, 545.551, tolerance = 1e-3)
  expect_equal(buyske$N$misclassification, 21400.41, tolerance = 1e-2)
  expect_equal(buyske$N$misclassification / buyske$N$no_error, 39.22715, tolerance = 1e-5)
  expect_equal(buyske$percent_increase$misclassification, 3822.715, tolerance = 1e-3)

  common <- list(
    target_power = .80, alpha = 1e-5, pd = .25, prev = .10,
    R1 = sqrt(1.5), R2 = 1.5, delta_prime = 1,
    misclass_rate = 0, verbose = FALSE
  )
  expect_equal(do.call(tdt_mssn, c(common, list(heter_rate = 0)))$N$heterogeneity, 3430.696, tolerance = 1e-3)
  expect_equal(do.call(tdt_mssn, c(common, list(heter_rate = .5)))$N$heterogeneity, 13722.79, tolerance = 1e-2)
})

test_that("TDT advanced details use separate canonical quantities", {
  values <- pawh:::.pawh_tdt_defaults()
  values$misclassification <- values$heterogeneity <- TRUE
  values$misclass_rate <- 5
  values$heter_rate <- 20
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  html <- paste(as.character(pawh:::.pawh_tdt_advanced_ui(calculation)), collapse = "\n")
  expect_match(html, "Advanced calculation details", fixed = TRUE)
  expect_false(grepl("<details[^>]* open", html))
  expect_match(html, "No-error design transmission quantities", fixed = TRUE)
  expect_match(html, "Phenotype misclassification transmission quantities", fixed = TRUE)
  expect_match(html, "Locus heterogeneity transmission quantities", fixed = TRUE)
  expect_match(html, formatC(calculation$result$gT_star$misclassification, format = "f", digits = 4), fixed = TRUE)
  expect_match(html, formatC(calculation$result$gNT_star$heterogeneity, format = "f", digits = 4), fixed = TRUE)
  expect_match(html, formatC(calculation$result$ET$no_error, format = "f", digits = 2), fixed = TRUE)
  expect_false(grepl("combined", html, ignore.case = TRUE))

  values$objective <- "mssn"
  values$input_mode <- "model_free"
  values$ET <- values$ENT <- 100
  null <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  null_html <- paste(as.character(pawh:::.pawh_tdt_advanced_ui(null)), collapse = "\n")
  expect_match(null_html, "Direct transmission quantities", fixed = TRUE)
  expect_match(null_html, ">Inf<", fixed = TRUE)
  expect_match(null_html, "not defined", fixed = TRUE)
  expect_false(grepl("NaN%|Inf%", null_html))
})

test_that("TDT reproducible calls exactly reproduce frozen calculations", {
  values <- pawh:::.pawh_tdt_defaults()
  values$misclassification <- values$heterogeneity <- TRUE
  values$misclass_rate <- 5
  values$heter_rate <- 20
  power <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  args <- pawh:::.pawh_tdt_repro_args(power)
  reproduced <- do.call(tdt_power, args)
  expect_equal(reproduced$power, power$result$power)
  expect_identical(args$misclass_rate, .05)
  expect_identical(args$heter_rate, .2)
  expect_false(any(c("objective", "display") %in% names(args)))
  expect_match(pawh:::.pawh_call_text(pawh:::.pawh_tdt_repro_call(power)), "tdt_power", fixed = TRUE)

  values$objective <- "mssn"
  mssn <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  reproduced <- do.call(tdt_mssn, pawh:::.pawh_tdt_repro_args(mssn))
  expect_equal(reproduced$N, mssn$result$N)
  expect_match(pawh:::.pawh_call_text(pawh:::.pawh_tdt_repro_call(mssn)), "tdt_mssn", fixed = TRUE)
})

test_that("shared colors and transmission labels use canonical values", {
  colors <- pawh:::.pawh_plot_colors()
  expect_identical(unname(colors[c("tdt_baseline", "tdt_misclassification", "tdt_heterogeneity")]),
    c("#3F4850", "#355C7D", "#6F879A"))
  expect_identical(unname(colors[c("transmitted", "nontransmitted")]), c("#3F4850", "#8FA1AF"))
  values <- pawh:::.pawh_tdt_defaults()
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  data <- pawh:::.pawh_tdt_transmissions(calculation)
  expect_equal(data$Probability, c(calculation$result$gT_star$no_error, calculation$result$gNT_star$no_error))
  expect_identical(data$Label, formatC(data$Probability, format = "f", digits = 3))
  plot <- pawh:::.pawh_tdt_transmission_plot(calculation)
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))))
})

test_that("TDT 3D surface is canonical, modest, and frozen", {
  skip_if_not_installed("plotly")
  values <- pawh:::.pawh_tdt_defaults()
  values$misclassification <- TRUE
  values$misclass_rate <- 5
  calculation <- pawh:::.pawh_tdt_calculate(pawh:::.pawh_tdt_snapshot(values))
  args <- pawh:::.pawh_tdt_surface_args(
    calculation, "misclassification", "pd", "misclass_rate", n = 4
  )
  expect_length(args$x_values, 4)
  expect_length(args$y_values, 4)
  expect_identical(args$alpha, .05)
  surface <- pawh:::.pawh_tdt_surface(
    calculation, "misclassification", "pd", "misclass_rate", n = 4
  )
  expect_s3_class(surface, "plotly")
  expect_equal(nrow(attr(surface, "surface_data")), 16)
  expect_identical(attr(surface, "surface_spec")$fixed_parameters$alpha, .05)

  shiny::testServer(pawh:::.pawh_tdt_server, {
    session$setInputs(
      objective = "power", input_mode = "model_based", N = 600,
      alpha = .05, prev = 5, pd = 30, R1 = 1.5, R2 = 2.25,
      delta_prime = 1, misclassification = TRUE, misclass_rate = 5,
      heterogeneity = FALSE, heter_rate = 0
    )
    session$setInputs(calculate = 1)
    expect_null(session$returned$surface())
    session$setInputs(alpha = .01)
    session$setInputs(
      surface_scenario = "misclassification", surface_x = "pd",
      surface_y = "misclass_rate", generate_surface = 1
    )
    expect_s3_class(session$returned$surface(), "plotly")
    expect_identical(
      attr(session$returned$surface(), "surface_spec")$fixed_parameters$alpha,
      .05
    )
  })
})
