falconer_print_args <- list(qtl_var = 0.025, tau = 0.5, pd = 0.15)
falconer_threshold_print_args <- c(
  falconer_print_args,
  list(x_upper = 5, x_lower = 5)
)

test_that("verbose FALSE suppresses all Falconer console output", {
  expect_silent(do.call(
    qtl_falconer_parameters,
    c(falconer_print_args, list(verbose = FALSE))
  ))
  expect_silent(do.call(
    qtl_falconer_threshold_parameters,
    c(falconer_threshold_print_args, list(verbose = FALSE))
  ))
  expect_silent(qtl_anova_power(
    N = 996, alpha = 0.0001,
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    verbose = FALSE
  ))
  expect_silent(qtl_anova_mssn(
    power = 0.8, alpha = 0.0001,
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    verbose = FALSE
  ))
  expect_silent(qtl_threshold_chisq_power(
    N_case = 126, alpha = 0.0001,
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    x_upper = 5, x_lower = 5, verbose = FALSE
  ))
  expect_silent(qtl_threshold_chisq_mssn(
    power = 0.8, alpha = 0.0001,
    qtl_var = 0.025, tau = 0.5, pd = 0.15,
    x_upper = 5, x_lower = 5, verbose = FALSE
  ))
})

test_that("verbose TRUE prints expected Falconer section headings", {
  parameter_output <- paste(capture.output(
    do.call(qtl_falconer_parameters, falconer_print_args), type = "message"
  ), collapse = "\n")
  threshold_output <- paste(capture.output(
    do.call(qtl_falconer_threshold_parameters, falconer_threshold_print_args),
    type = "message"
  ), collapse = "\n")
  anova_power_output <- paste(capture.output(
    qtl_anova_power(996, 0.0001, 0.025, 0.5, 0.15), type = "message"
  ), collapse = "\n")
  anova_mssn_output <- paste(capture.output(
    qtl_anova_mssn(0.8, 0.0001, 0.025, 0.5, 0.15), type = "message"
  ), collapse = "\n")
  chisq_power_output <- paste(capture.output(
    qtl_threshold_chisq_power(126, 0.0001, 0.025, 0.5, 0.15, 5, 5),
    type = "message"
  ), collapse = "\n")
  chisq_mssn_output <- paste(capture.output(
    qtl_threshold_chisq_mssn(0.8, 0.0001, 0.025, 0.5, 0.15, 5, 5),
    type = "message"
  ), collapse = "\n")

  expect_match(
    parameter_output,
    "(?s)Input Parameters.*Derived Parameters.*Genotype-Specific Quantities.*Validation",
    perl = TRUE
  )
  expect_match(
    threshold_output,
    "(?s)Threshold Information.*Derived Falconer Parameters.*Selected Prevalences",
    perl = TRUE
  )
  expect_match(
    anova_power_output,
    "(?s)Study Design.*Falconer Model.*Test Results",
    perl = TRUE
  )
  expect_match(
    anova_mssn_output,
    "(?s)Target Design.*Required Sample Size.*Genotype-Specific Quantities",
    perl = TRUE
  )
  expect_match(anova_power_output, "Genotype count method:\\s+Rounded")
  expect_match(anova_mssn_output, "Genotype count method:\\s+Rounded")
  expect_match(
    chisq_power_output,
    "(?s)Study Design.*Threshold Selection.*Selected Population.*Test Results",
    perl = TRUE
  )
  expect_match(
    chisq_mssn_output,
    "(?s)Target Design.*Required Selected Sample.*Expected Screening Burden",
    perl = TRUE
  )
})

test_that("printing does not alter Falconer return values", {
  quiet_parameters <- do.call(
    qtl_falconer_parameters,
    c(falconer_print_args, list(verbose = FALSE))
  )
  printed_parameters <- suppressMessages(do.call(
    qtl_falconer_parameters,
    c(falconer_print_args, list(verbose = TRUE))
  ))

  quiet_threshold <- do.call(
    qtl_falconer_threshold_parameters,
    c(falconer_threshold_print_args, list(verbose = FALSE))
  )
  printed_threshold <- suppressMessages(do.call(
    qtl_falconer_threshold_parameters,
    c(falconer_threshold_print_args, list(verbose = TRUE))
  ))

  quiet_anova_power <- qtl_anova_power(
    996, 0.0001, 0.025, 0.5, 0.15, verbose = FALSE
  )
  printed_anova_power <- suppressMessages(qtl_anova_power(
    996, 0.0001, 0.025, 0.5, 0.15, verbose = TRUE
  ))

  quiet_anova_mssn <- qtl_anova_mssn(
    0.8, 0.0001, 0.025, 0.5, 0.15, verbose = FALSE
  )
  printed_anova_mssn <- suppressMessages(qtl_anova_mssn(
    0.8, 0.0001, 0.025, 0.5, 0.15, verbose = TRUE
  ))

  quiet_chisq_power <- qtl_threshold_chisq_power(
    126, 0.0001, 0.025, 0.5, 0.15, 5, 5, verbose = FALSE
  )
  printed_chisq_power <- suppressMessages(qtl_threshold_chisq_power(
    126, 0.0001, 0.025, 0.5, 0.15, 5, 5, verbose = TRUE
  ))

  quiet_chisq_mssn <- qtl_threshold_chisq_mssn(
    0.8, 0.0001, 0.025, 0.5, 0.15, 5, 5, verbose = FALSE
  )
  printed_chisq_mssn <- suppressMessages(qtl_threshold_chisq_mssn(
    0.8, 0.0001, 0.025, 0.5, 0.15, 5, 5, verbose = TRUE
  ))

  expect_equal(printed_parameters, quiet_parameters)
  expect_equal(printed_threshold, quiet_threshold)
  expect_equal(printed_anova_power, quiet_anova_power)
  expect_equal(printed_anova_mssn, quiet_anova_mssn)
  expect_equal(printed_chisq_power, quiet_chisq_power)
  expect_equal(printed_chisq_mssn, quiet_chisq_mssn)
})
