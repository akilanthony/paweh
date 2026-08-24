test_that("Case-Control UI exposes the production workflow", {
  html <- paste(as.character(paweh:::.paweh_case_control_ui("cc")), collapse = "\n")
  expect_match(html, "Estimate power")
  expect_match(html, "Minimum sample size")
  expect_match(html, "Direct genotype probabilities")
  expect_match(html, "Advanced assumptions")
  expect_match(html, "Sensitivity")
  expect_match(html, "Visualize")
  expect_match(html, "Methods")
})

test_that("model and direct snapshots reproduce canonical power and MSSN", {
  v <- paweh:::.paweh_cc_defaults()
  calc <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(v))
  direct <- cc_power(N_case = 500, alpha = .05, input_mode = "model_based",
    prev = .1, pd = .3, R2 = 2, MOI = "M", k = 1, verbose = FALSE)
  expect_equal(calc$adjusted$tests, direct$tests)

  v$objective <- "mssn"; v$input_mode <- "model_free"
  calc <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(v))
  direct <- cc_mssn(power = .8, alpha = .05, input_mode = "model_free",
    g1 = c(.35,.45,.2), g0 = c(.49,.42,.09), k = 1, verbose = FALSE)
  expect_equal(calc$adjusted$tests, direct$tests)
})

test_that("advanced modifiers and zero-error narration follow canonical metadata", {
  v <- paweh:::.paweh_cc_defaults()
  v$locus_het <- TRUE; v$pi <- 75; v$pheno_misclass <- TRUE
  v$theta <- 3; v$phi <- 1; v$genotype_error <- TRUE; v$e <- 2
  calc <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(v))
  direct <- cc_power(N_case=500,alpha=.05,input_mode="model_based",prev=.1,
    pd=.3,R2=2,MOI="M",k=1,locus_het=TRUE,pi=.75,
    pheno_misclass=TRUE,theta=.03,phi=.01,geno_misclass="1p",e=.02,
    verbose=FALSE)
  expect_equal(calc$adjusted$tests, direct$tests)
  expect_true(all(unlist(calc$active)))

  for (model in c("1p","2p","3p","diff3p")) {
    z <- paweh:::.paweh_cc_defaults(); z$genotype_error <- TRUE; z$geno_misclass <- model
    zero <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(z))
    expect_false(zero$active$genotype)
    expect_equal(zero$adjusted$tests, zero$baseline$tests)
  }
})

test_that("validation is friendly and calculated state remains frozen", {
  v <- paweh:::.paweh_cc_defaults(); v$input_mode <- "model_free"; v$g1_2 <- .3
  expect_error(paweh:::.paweh_cc_snapshot(v), "sum to 1")
  shiny::testServer(paweh:::.paweh_case_control_server, {
    session$setInputs(objective="power",input_mode="model_based",N_case=500,
      alpha=.05,k=1,prev=10,pd=30,MOI="M",R2=2)
    session$setInputs(calculate=1)
    first <- session$returned$calculation()$adjusted$tests$genotypes$power
    session$setInputs(alpha=.01)
    session$flushReact()
    expect_true(session$returned$changed())
    expect_equal(session$returned$calculation()$adjusted$tests$genotypes$power, first)
  })
})

test_that("sensitivity points are direct canonical calls of frozen input", {
  v <- paweh:::.paweh_cc_defaults()
  calc <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(v))
  sens <- paweh:::.paweh_cc_sensitivity(calc, "N_case", c(400,600), n=5)
  expect_equal(nrow(sens$data), 10)
  expected <- cc_power(N_case=500,alpha=.05,input_mode="model_based",prev=.1,
    pd=.3,R2=2,MOI="M",k=1,verbose=FALSE)
  middle <- subset(sens$data, x == 500)
  expect_equal(middle$y, c(expected$tests$genotypes$power, expected$tests$trend$power))
  expect_equal(sens$baseline_x, 500)
})

test_that("Case-Control plots use semantic colors and redundant line styles", {
  v <- paweh:::.paweh_cc_defaults()
  calc <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(v))
  sens <- paweh:::.paweh_cc_sensitivity(calc, "N_case", c(400, 600), n = 5)

  sensitivity_plot <- paweh:::.paweh_cc_sensitivity_plot(sens)
  genotype_plot <- paweh:::.paweh_cc_genotype_plot(calc)

  expect_s3_class(sensitivity_plot, "ggplot")
  expect_s3_class(genotype_plot, "ggplot")
  expect_true(any(vapply(sensitivity_plot$scales$scales,
    function(scale) "colour" %in% scale$aesthetics, logical(1))))
  expect_true(any(vapply(sensitivity_plot$scales$scales,
    function(scale) "linetype" %in% scale$aesthetics, logical(1))))
  expect_true(any(vapply(genotype_plot$scales$scales,
    function(scale) "fill" %in% scale$aesthetics, logical(1))))
})

test_that("Case-Control advanced details are collapsed and canonical", {
  values <- paweh:::.paweh_cc_defaults()
  calculation <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(values))
  html <- paste(as.character(paweh:::.paweh_cc_advanced_ui(calculation)), collapse = "\n")
  expect_match(html, "Advanced calculation details", fixed = TRUE)
  expect_match(html, "Case genotype probabilities", fixed = TRUE)
  expect_match(html, formatC(calculation$adjusted$freqs$g_base_case[1], format = "f", digits = 4), fixed = TRUE)
  expect_match(html, formatC(calculation$adjusted$tests$genotypes$lambda, format = "f", digits = 4), fixed = TRUE)
  expect_false(grepl("<details[^>]* open", html))

  values$input_mode <- "model_free"
  direct <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(values))
  direct_html <- paste(as.character(paweh:::.paweh_cc_advanced_ui(direct)), collapse = "\n")
  expect_match(direct_html, "Direct genotype probabilities", fixed = TRUE)
  expect_match(direct_html, "0.3500", fixed = TRUE)
  expect_equal(direct$adjusted$freqs$g_base_case, c(.35, .45, .2))

  values <- paweh:::.paweh_cc_defaults()
  values$locus_het <- TRUE; values$pi <- 75
  values$pheno_misclass <- TRUE; values$theta <- 3; values$phi <- 1
  values$genotype_error <- TRUE; values$e <- 2
  modified <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(values))
  modifier_html <- paste(as.character(paweh:::.paweh_cc_advanced_ui(modified)), collapse = "\n")
  expect_match(modifier_html, "Disease attributable to locus", fixed = TRUE)
  expect_match(modifier_html, "Affected classified as control", fixed = TRUE)
  expect_match(modifier_html, modified$adjusted$errors$genotype_misclass$model, fixed = TRUE)
  expect_match(modifier_html, "Observed genotype probabilities after misclassification", fixed = TRUE)
})

test_that("Case-Control reproducible calls use frozen public arguments", {
  values <- paweh:::.paweh_cc_defaults()
  values$locus_het <- TRUE; values$pi <- 80
  values$pheno_misclass <- TRUE; values$theta <- 2; values$phi <- 1
  power <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(values))
  args <- paweh:::.paweh_cc_repro_args(power)
  reproduced <- do.call(cc_power, args)
  expect_equal(reproduced$tests, power$adjusted$tests)
  expect_identical(args$pi, .8)
  expect_identical(args$theta, .02)
  expect_false(any(c("target_power", "objective", "display") %in% names(args)))
  expect_match(paweh:::.paweh_call_text(paweh:::.paweh_cc_repro_call(power)), "cc_power", fixed = TRUE)

  values$objective <- "mssn"
  mssn <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(values))
  reproduced <- do.call(cc_mssn, paweh:::.paweh_cc_repro_args(mssn))
  expect_equal(reproduced$tests, mssn$adjusted$tests)
  expect_match(paweh:::.paweh_call_text(paweh:::.paweh_cc_repro_call(mssn)), "cc_mssn", fixed = TRUE)
})

test_that("power sensitivity zoom note detects only narrow ranges", {
  narrow <- list(objective = "power", data = data.frame(y = c(.991, .999)))
  broad <- list(objective = "power", data = data.frame(y = c(.1, .9)))
  mssn <- list(objective = "mssn", data = data.frame(y = c(100, 200)))
  expect_true(paweh:::.paweh_power_axis_zoomed(narrow))
  expect_false(paweh:::.paweh_power_axis_zoomed(broad))
  expect_false(paweh:::.paweh_power_axis_zoomed(mssn))
})

test_that("Case-Control Methods records the frozen analysis", {
  values <- paweh:::.paweh_cc_defaults()
  values$locus_het <- TRUE
  values$pi <- 80
  calculation <- paweh:::.paweh_cc_calculate(paweh:::.paweh_cc_snapshot(values))
  html <- paste(as.character(paweh:::.paweh_cc_methods_ui(calculation)), collapse = "\n")
  expect_match(html, "Input specification", fixed = TRUE)
  expect_match(html, "Genotype 2 x 3 chi-square", fixed = TRUE)
  expect_match(html, "Cochran-Armitage trend", fixed = TRUE)
  expect_match(html, "Locus heterogeneity", fixed = TRUE)
  expect_match(html, "cc_power()", fixed = TRUE)
})
