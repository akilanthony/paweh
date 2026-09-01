test_that("CC-NGS power plot data exactly delegate to the public API", {
  for (moi in c("M", "D", "Rec")) {
    dat <- plot_ngs_power(
      design = "cc", coverage = c(4, 2, 4), seq_error = c(0.01, 0),
      return_data = TRUE,
      N_case = 900, alpha = 0.05, prev = 0.05, pd = 0.30,
      R2 = 1.8, MOI = moi, k = 1
    )
    expect_identical(dat$coverage, c(2, 4, 2, 4))
    expect_identical(dat$seq_error, c(0, 0, 0.01, 0.01))
    direct <- mapply(function(v, e) {
      cc_ngs_power(
        N_case = 900, alpha = 0.05, prev = 0.05, pd = 0.30,
        R2 = 1.8, coverage = v, seq_error = e, MOI = moi,
        k = 1, verbose = FALSE
      )
    }, dat$coverage, dat$seq_error, SIMPLIFY = FALSE)
    expect_equal(dat$lambda, vapply(direct, `[[`, numeric(1), "lambda"))
    expect_equal(dat$power, vapply(direct, `[[`, numeric(1), "power"))
  }
})

test_that("CC-NGS plots retain heterogeneity and exact-null semantics", {
  pis <- c(1, 0.75, 0.5, 0)
  power_dat <- plot_ngs_power(
    design = "cc", coverage = c(2, 4), seq_error = 0.005,
    return_data = TRUE,
    N_case = 900, alpha = 0.05, prev = 0.05, pd = 0.30,
    R2 = 1.8, MOI = "M", k = 1, locus_het = TRUE, pi = pis
  )
  expect_equal(sort(unique(power_dat$pi), decreasing = TRUE), pis)
  expect_true(all(power_dat$lambda[power_dat$pi == 0] == 0))
  expect_equal(power_dat$power[power_dat$pi == 0], rep(0.05, 2),
               tolerance = 1e-14)

  for (i in seq_len(nrow(power_dat))) {
    direct <- cc_ngs_power(
      N_case = 900, alpha = 0.05, prev = 0.05, pd = 0.30,
      R2 = 1.8, coverage = power_dat$coverage[i],
      seq_error = power_dat$seq_error[i], MOI = "M", k = 1,
      locus_het = TRUE, pi = power_dat$pi[i], verbose = FALSE
    )
    expect_equal(power_dat$power[i], direct$power)
    expect_equal(power_dat$lambda[i], direct$lambda)
  }

  mssn_dat <- plot_ngs_mssn(
    design = "cc", coverage = c(2, 4), seq_error = 0.005,
    return_data = TRUE,
    power = 0.80, alpha = 0.05, prev = 0.05, pd = 0.30,
    R2 = 1.8, MOI = "M", k = 1, locus_het = TRUE, pi = pis
  )
  finite <- mssn_dat$pi > 0
  for (i in which(finite)) {
    direct <- cc_ngs_mssn(
      power = 0.80, alpha = 0.05, prev = 0.05, pd = 0.30,
      R2 = 1.8, coverage = mssn_dat$coverage[i],
      seq_error = mssn_dat$seq_error[i], MOI = "M", k = 1,
      locus_het = TRUE, pi = mssn_dat$pi[i], verbose = FALSE
    )
    expect_equal(mssn_dat$MSSN_case[i], direct$MSSN_case)
    expect_equal(mssn_dat$MSSN_ctrl[i], direct$MSSN_ctrl)
    expect_equal(mssn_dat$MSSN_total[i], direct$MSSN_total)
    expect_equal(mssn_dat$achieved_power[i], direct$achieved_power)
  }
  null <- mssn_dat[!finite, ]
  expect_true(all(is.na(null$MSSN_case)))
  expect_true(all(is.na(null$MSSN_ctrl)))
  expect_true(all(is.na(null$MSSN_total)))
  expect_true(all(!null$finite_mssn))
  expect_true(all(null$status == "no finite MSSN"))
})

test_that("CC-NGS strict master switch and validation errors propagate", {
  expect_error(
    plot_ngs_power(
      "cc", 4, 0.01, N_case = 900, alpha = 0.05, prev = 0.05,
      pd = 0.30, R2 = 1.8, MOI = "M", locus_het = FALSE, pi = 0.5
    ),
    "pi is used only when locus_het = TRUE; set pi = 1 or enable locus heterogeneity\\."
  )
  expect_error(
    plot_ngs_mssn(
      "cc", 4, 0.01, power = 0.8, alpha = 0.05, prev = 0.05,
      pd = 0.30, R2 = -1, MOI = "M"
    ),
    "R2 must be"
  )
})

test_that("TDT-NGS power plot data exactly delegate to the public API", {
  dat <- plot_ngs_power(
    design = "tdt", coverage = c(4, 2, 4), seq_error = c(0.005, 0.01),
    return_data = TRUE,
    N = 500, pd = 0.325, R1 = 1.2, alpha = 5e-8
  )
  expect_identical(dat$coverage, c(2, 4, 2, 4))
  for (i in seq_len(nrow(dat))) {
    direct <- tdt_ngs_power(
      N = 500, pd = 0.325, R1 = 1.2, coverage = dat$coverage[i],
      seq_error = dat$seq_error[i], alpha = 5e-8, verbose = FALSE
    )
    expect_equal(dat$lambda[i], direct$lambda)
    expect_equal(dat$power[i], direct$power)
  }
})

test_that("TDT-NGS MSSN plot data exactly delegate to the public API", {
  dat <- plot_ngs_mssn(
    design = "tdt", coverage = c(2, 4), seq_error = c(0.005, 0.01),
    return_data = TRUE,
    power = 0.80, pd = 0.325, R1 = 1.2, alpha = 5e-8
  )
  expect_true(all(dat$finite_mssn))
  expect_true(all(dat$status == "finite"))
  for (i in seq_len(nrow(dat))) {
    direct <- tdt_ngs_mssn(
      power = 0.80, pd = 0.325, R1 = 1.2,
      coverage = dat$coverage[i], seq_error = dat$seq_error[i],
      alpha = 5e-8, verbose = FALSE
    )
    expect_equal(dat$MSSN_trios[i], direct$MSSN_trios)
    expect_equal(dat$total_individuals[i], direct$total_individuals)
    expect_equal(dat$achieved_power[i], direct$achieved_power)
  }
})

test_that("sequencing plots validate grids and reject TDT heterogeneity", {
  expect_error(plot_ngs_power("tdt", 1, 0.005), "not identifiable")
  expect_error(plot_ngs_power("cc", c(2, 2.5), 0.005), "finite integers")
  expect_error(plot_ngs_power("cc", 2, c(-0.01, 0.01)), "\\[0, 0.5\\)")
  expect_error(plot_ngs_power("cc", 2, 0.01, target_power = 1),
               "target_power")
  expect_error(
    plot_ngs_power(
      "tdt", 2, 0.005, N = 500, pd = 0.325, R1 = 1.2,
      alpha = 0.05, pi = 1
    ),
    "does not support locus heterogeneity"
  )
  expect_false("pi" %in% names(formals(plot_ngs_power)))
  expect_false("locus_het" %in% names(formals(plot_ngs_mssn)))
})

test_that("sequencing ggplots expose exact data, labels, and reference line", {
  p <- plot_ngs_power(
    "cc", c(2, 4), c(0, 0.01), target_power = 0.8,
    N_case = 900, alpha = 0.05, prev = 0.05, pd = 0.30,
    R2 = 1.8, MOI = "M", k = 1
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Sequencing coverage (\u00d7)")
  expect_equal(p$labels$y, "Power")
  expect_true(all(c("coverage", "seq_error", "lambda", "power") %in%
                    names(p$data)))
  expect_gte(length(p$layers), 3L)

  m <- plot_ngs_mssn(
    "cc", c(2, 4), 0.005,
    power = 0.8, alpha = 0.05, prev = 0.05, pd = 0.30,
    R2 = 1.8, MOI = "M", locus_het = TRUE, pi = c(1, 0)
  )
  expect_s3_class(m, "ggplot")
  expect_equal(m$labels$y, "Required cases")
  expect_true(any(!m$data$finite_mssn))
  expect_true(any(is.na(m$data$MSSN_case)))
})

test_that("sequencing plotting remains deterministic and formula-free", {
  args <- list(
    design = "tdt", coverage = c(2, 4), seq_error = 0.005,
    return_data = TRUE, N = 500, pd = 0.325, R1 = 1.2, alpha = 0.05
  )
  expect_identical(do.call(plot_ngs_power, args), do.call(plot_ngs_power, args))

  plotting_functions <- list(
    plot_ngs_power, plot_ngs_mssn, .plot_ngs_validate_coverage,
    .plot_ngs_validate_seq_error, .plot_ngs_cc_settings,
    .plot_ngs_reject_tdt_heterogeneity, .plot_ngs_series, .plot_ngs_line
  )
  source <- paste(vapply(
    plotting_functions,
    function(fun) paste(deparse(body(fun)), collapse = "\n"),
    character(1)
  ), collapse = "\n")
  forbidden <- c(
    ".cc_ngs_ahn_ncp", ".tdt_ngs_ncp", "ngs_genotype_error_matrix",
    ".ngs_call_genotype_ml", "rbinom", "runif", "sample(", "set.seed"
  )
  expect_false(any(vapply(forbidden, grepl, logical(1), x = source,
                          fixed = TRUE)))
})
