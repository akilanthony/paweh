test_that("Gordon 2002 DSB genotype-error MSSNs are reproduced", {
  # Gordon et al. Human Heredity 2002;54(1):22-33.
  # doi:10.1159/000066696, Figure 1b model-free DSB example.
  g_case <- c(0.36, 0.48, 0.16)
  g_control <- c(0.25, 0.50, 0.25)

  no_error <- cc_mssn(
    power = 0.95, alpha = 0.05,
    input_mode = "model_free", g1 = g_case, g0 = g_control, k = 1,
    geno_misclass = "2p", e1 = 0, e2 = 0, verbose = FALSE
  )
  dsb_10 <- cc_mssn(
    power = 0.95, alpha = 0.05,
    input_mode = "model_free", g1 = g_case, g0 = g_control, k = 1,
    geno_misclass = "2p", e1 = 0.10, e2 = 0.05, verbose = FALSE
  )

  expect_equal(no_error$tests$genotypes$MSSN_case, 387)
  expect_equal(dsb_10$tests$genotypes$MSSN_case, 477)
})

test_that("Buyske 2009 TDT phenotype-error inflation is reproduced", {
  # Buyske et al. Human Heredity 2009;67(4):287-292.
  # doi:10.1159/000194981, approximately 39-fold required-sample ratio.
  out <- tdt_mssn(
    target_power = 0.80, alpha = 1e-5,
    pd = 0.50, prev = 0.01, R1 = sqrt(2.5), R2 = 2.5,
    delta_prime = 1, misclass_rate = 0.05, heter_rate = 0,
    verbose = FALSE
  )

  expect_equal(out$N$no_error, 545.551, tolerance = 1e-3)
  expect_equal(out$N$misclassification, 21400.41, tolerance = 1e-2)
  expect_equal(
    out$N$misclassification / out$N$no_error,
    39.22715,
    tolerance = 1e-5
  )
})

test_that("Chen 2009 TDT locus-heterogeneity MSSNs are reproduced", {
  # Chen et al. Stat Appl Genet Mol Biol 2009;8:Article 44.
  # doi:10.2202/1544-6115.1501.
  common <- list(
    target_power = 0.80, alpha = 1e-5,
    pd = 0.25, prev = 0.10, R1 = sqrt(1.5), R2 = 1.5,
    delta_prime = 1, misclass_rate = 0, verbose = FALSE
  )
  no_heterogeneity <- do.call(tdt_mssn, c(common, list(heter_rate = 0)))
  half_heterogeneous <- do.call(tdt_mssn, c(common, list(heter_rate = 0.50)))

  expect_equal(
    no_heterogeneity$N$heterogeneity,
    3430.696,
    tolerance = 1e-3
  )
  expect_equal(
    half_heterogeneous$N$heterogeneity,
    13722.79,
    tolerance = 1e-2
  )
})

test_that("Gordon 2017 multivariate Pillai MSSN is reproduced", {
  # Gordon et al. Human Heredity 2017;81(4):194-209.
  # doi:10.1159/000457135.
  out <- qtl_multivariate_mssn_full(
    power = 0.80, alpha = 5e-8,
    qtl_var = c(0.10, 0.05), tau = c(0, 0.50), pd = 0.05,
    cor_matrix = diag(2), test = "pillai", verbose = FALSE
  )

  expect_identical(out$N, 326L)
  expect_equal(out$historical_fractional_mssn, 325.5057, tolerance = 1e-3)
})
