capture_paweh_report <- function(expr) {
  paste(capture.output(force(expr), type = "message"), collapse = "\n")
}

cc_console_power_args <- list(
  N_case = 500,
  alpha = 0.05,
  input_mode = "model_free",
  g1 = c(0.25, 0.50, 0.25),
  g0 = c(0.36, 0.48, 0.16),
  geno_misclass = "diff3p",
  diff_source = "case",
  diff_multiplier = 0.5,
  case_e01 = 0.02,
  case_e02 = 0.01,
  case_e03 = 0.005
)

cc_console_mssn_args <- c(
  list(power = 0.80),
  cc_console_power_args[setdiff(names(cc_console_power_args), "N_case")]
)

tdt_console_power_args <- list(
  N = 600,
  input_mode = "model_free",
  ET = 140,
  ENT = 100,
  pd = 0.30,
  prev = 0.05,
  misclass_rate = 0.01,
  heter_rate = 0.10
)

tdt_console_mssn_args <- list(
  target_power = 0.80,
  input_mode = "model_based",
  pd = 0.30,
  prev = 0.05,
  R1 = 1.5,
  R2 = 2.25
)

cc_ngs_console_power_args <- list(
  N_case = 500,
  alpha = 0.05,
  prev = 0.05,
  pd = 0.30,
  R2 = 2.25,
  coverage = 4,
  seq_error = 0.01,
  MOI = "M"
)

cc_ngs_console_mssn_args <- c(
  list(power = 0.80),
  cc_ngs_console_power_args[setdiff(names(cc_ngs_console_power_args), "N_case")]
)

tdt_ngs_console_power_args <- list(
  N = 654,
  pd = 0.35,
  R1 = 2,
  coverage = 4,
  seq_error = 0.025,
  alpha = 5e-8
)

tdt_ngs_console_mssn_args <- list(
  power = 0.90,
  pd = 0.35,
  R1 = 2,
  coverage = 4,
  seq_error = 0.025,
  alpha = 5e-8
)

test_that("ordinary CC reports are detailed and input-mode aware", {
  power_text <- capture_paweh_report(do.call(
    cc_power,
    c(cc_console_power_args, list(verbose = TRUE))
  ))
  expect_match(power_text, "PAWEH Case-Control Study", fixed = TRUE)
  expect_match(power_text, "Study Design", fixed = TRUE)
  expect_match(power_text, "No-Error Design", fixed = TRUE)
  expect_match(power_text, "Adjusted Design", fixed = TRUE)
  expect_match(power_text, "Genotype Frequencies", fixed = TRUE)
  expect_match(power_text, "Differential 3-parameter", fixed = TRUE)
  expect_match(power_text, "Model-free supplied frequencies", fixed = TRUE)
  expect_false(grepl("Disease prevalence", power_text, fixed = TRUE))

  mssn_text <- capture_paweh_report(do.call(
    cc_mssn,
    c(cc_console_mssn_args, list(verbose = TRUE))
  ))
  expect_match(mssn_text, "Minimum Sample Size Necessary", fixed = TRUE)
  expect_match(mssn_text, "Target Design", fixed = TRUE)
  expect_match(mssn_text, "MSSN inflation", fixed = TRUE)
})

test_that("ordinary TDT reports expose expected transmission information", {
  power_text <- capture_paweh_report(do.call(
    tdt_power,
    c(tdt_console_power_args, list(verbose = TRUE))
  ))
  expect_match(
    power_text, "PAWEH Transmission Disequilibrium Test", fixed = TRUE
  )
  expect_match(power_text, "Expected transmitted count (ET)", fixed = TRUE)
  expect_match(power_text, "Expected non-transmitted probability", fixed = TRUE)
  expect_match(power_text, "Phenotype Misclassification", fixed = TRUE)
  expect_match(power_text, "Locus Heterogeneity", fixed = TRUE)

  mssn_text <- capture_paweh_report(do.call(
    tdt_mssn,
    c(tdt_console_mssn_args, list(verbose = TRUE))
  ))
  expect_match(mssn_text, "Minimum Sample Size Necessary", fixed = TRUE)
  expect_match(mssn_text, "Required complete trios", fixed = TRUE)
  expect_match(mssn_text, "Expected transmitted probability", fixed = TRUE)
})

test_that("CC-NGS reports show the calling bridge without raw-read claims", {
  power_text <- capture_paweh_report(do.call(
    cc_ngs_power,
    c(cc_ngs_console_power_args, list(verbose = TRUE))
  ))
  expect_match(power_text, "PAWEH Case-Control NGS Study", fixed = TRUE)
  expect_match(power_text, "Sequencing Model", fixed = TRUE)
  expect_match(power_text, "Deterministic maximum-likelihood", fixed = TRUE)
  expect_match(power_text, "True Genotype Frequencies", fixed = TRUE)
  expect_match(power_text, "Genotype-Call Matrix", fixed = TRUE)
  expect_match(power_text, "Rows = true genotype", fixed = TRUE)
  expect_match(power_text, "Called Genotype Frequencies", fixed = TRUE)
  expect_false(grepl("Raw-read", power_text, fixed = TRUE))

  mssn_text <- capture_paweh_report(do.call(
    cc_ngs_mssn,
    c(cc_ngs_console_mssn_args, list(verbose = TRUE))
  ))
  expect_match(mssn_text, "Minimum Sample Size Necessary", fixed = TRUE)
  expect_match(mssn_text, "Continuous case requirement", fixed = TRUE)
  expect_match(mssn_text, "Achieved power", fixed = TRUE)
})

test_that("TDT1-NGS reports remain raw-read and information based", {
  power_text <- capture_paweh_report(do.call(
    tdt_ngs_power,
    c(tdt_ngs_console_power_args, list(verbose = TRUE))
  ))
  expect_match(power_text, "PAWEH TDT1-NGS", fixed = TRUE)
  expect_match(power_text, "Raw-read TDT1-NGS likelihood", fixed = TRUE)
  expect_match(power_text, "Efficient information per trio", fixed = TRUE)
  expect_false(grepl("Genotype-Call Matrix", power_text, fixed = TRUE))

  mssn_text <- capture_paweh_report(do.call(
    tdt_ngs_mssn,
    c(tdt_ngs_console_mssn_args, list(verbose = TRUE))
  ))
  expect_match(mssn_text, "Minimum Sample Size Necessary", fixed = TRUE)
  expect_match(mssn_text, "Required complete trios", fixed = TRUE)
  expect_match(mssn_text, "Total individuals", fixed = TRUE)
  expect_match(mssn_text, "Efficient information per trio", fixed = TRUE)
})

test_that("verbose FALSE suppresses all formatted reports", {
  calls <- list(
    list(cc_power, cc_console_power_args),
    list(cc_mssn, cc_console_mssn_args),
    list(tdt_power, tdt_console_power_args),
    list(tdt_mssn, tdt_console_mssn_args),
    list(cc_ngs_power, cc_ngs_console_power_args),
    list(cc_ngs_mssn, cc_ngs_console_mssn_args),
    list(tdt_ngs_power, tdt_ngs_console_power_args),
    list(tdt_ngs_mssn, tdt_ngs_console_mssn_args)
  )
  for (call in calls) {
    messages <- capture.output(
      invisible(do.call(call[[1L]], c(call[[2L]], list(verbose = FALSE)))),
      type = "message"
    )
    output <- capture.output(
      invisible(do.call(call[[1L]], c(call[[2L]], list(verbose = FALSE))))
    )
    expect_length(messages, 0L)
    expect_length(output, 0L)
  }
})

test_that("verbose presentation does not alter result objects", {
  calls <- list(
    list(cc_power, cc_console_power_args),
    list(cc_mssn, cc_console_mssn_args),
    list(tdt_power, tdt_console_power_args),
    list(tdt_mssn, tdt_console_mssn_args),
    list(cc_ngs_power, cc_ngs_console_power_args),
    list(cc_ngs_mssn, cc_ngs_console_mssn_args),
    list(tdt_ngs_power, tdt_ngs_console_power_args),
    list(tdt_ngs_mssn, tdt_ngs_console_mssn_args)
  )
  for (call in calls) {
    quiet <- do.call(call[[1L]], c(call[[2L]], list(verbose = FALSE)))
    capture.output(
      loud <- do.call(call[[1L]], c(call[[2L]], list(verbose = TRUE))),
      type = "message"
    )
    expect_identical(loud, quiet)
    expect_identical(class(loud), class(quiet))
    expect_identical(names(loud), names(quiet))
  }
})
