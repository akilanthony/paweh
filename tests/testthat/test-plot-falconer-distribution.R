test_that("single-trait density plot uses the validated Falconer backend", {
  dat <- plot_qtl_genotype_distribution(
    qtl_var = 0.5, tau = 0, pd = 0.25,
    type = "density", scale = "frequency", return_data = TRUE
  )
  model <- attr(dat, "falconer_model")

  expect_s3_class(dat, "data.frame")
  expect_equal(model$a, 1.154700538, tolerance = 1e-9)
  expect_equal(model$delta, 0)
  expect_equal(model$m, 0.577350269, tolerance = 1e-9)
  expect_equal(unname(model$mu), c(-0.577350269, 0.577350269, 1.732050808), tolerance = 1e-9)
  expect_equal(model$residual_variance, 0.5)
  expect_equal(model$residual_sd, sqrt(0.5))
  expect_equal(unname(model$pi), c(0.5625, 0.375, 0.0625))
  expect_equal(unique(dat$genotype_frequency), c(0.5625, 0.375, 0.0625))
  expect_equal(dat$value, dat$density * dat$genotype_frequency)
})

test_that("small-QTL reference values produce nearly coincident means", {
  dat <- plot_qtl_genotype_distribution(
    qtl_var = 0.0005, tau = 0, pd = 0.25, return_data = TRUE
  )
  model <- attr(dat, "falconer_model")

  expect_equal(model$a, 0.0365148371670111, tolerance = 1e-12)
  expect_equal(model$m, 0.0182574185835055, tolerance = 1e-12)
  expect_equal(
    unname(model$mu),
    c(-0.0182574185835055, 0.0182574185835055, 0.0547722557505166),
    tolerance = 1e-12
  )
  expect_equal(model$residual_variance, 0.9995)
  expect_equal(model$residual_sd, 0.999749968742185, tolerance = 1e-12)
  expect_lt(diff(range(model$mu)) / model$residual_sd, 0.08)
})

test_that("single-trait plots return customizable ggplot objects", {
  density_plot <- plot_qtl_genotype_distribution(0.5, 0, 0.25)
  histogram_plot <- plot_qtl_genotype_distribution(
    0.5, 0, 0.25, type = "histogram", scale = "frequency",
    n = 500, seed = 42, show_means = FALSE
  )

  expect_s3_class(density_plot, "ggplot")
  expect_s3_class(histogram_plot, "ggplot")
  expect_s3_class(density_plot + ggplot2::labs(subtitle = "custom"), "ggplot")
  expect_false(any(vapply(
    histogram_plot$layers, function(x) inherits(x$geom, "GeomVline"), logical(1)
  )))
})

test_that("density mode is analytic and does not alter RNG state", {
  set.seed(918)
  before <- .Random.seed
  plot_qtl_genotype_distribution(0.5, 0, 0.25, type = "density")
  expect_identical(.Random.seed, before)
})

test_that("histogram simulation is reproducible and locally seeded", {
  set.seed(812)
  before <- .Random.seed
  one <- plot_qtl_genotype_distribution(
    0.5, 0, 0.25, type = "histogram", n = 1000,
    seed = 99, return_data = TRUE
  )
  expect_identical(.Random.seed, before)
  two <- plot_qtl_genotype_distribution(
    0.5, 0, 0.25, type = "histogram", n = 1000,
    seed = 99, return_data = TRUE
  )
  expect_equal(one, two)

  large <- plot_qtl_genotype_distribution(
    0.5, 0, 0.25, type = "histogram", n = 20000,
    seed = 101, return_data = TRUE
  )
  observed <- prop.table(table(large$genotype))
  expect_equal(as.numeric(observed), c(0.5625, 0.375, 0.0625), tolerance = 0.01)
})

test_that("single-trait visualization inputs are validated", {
  expect_error(plot_qtl_genotype_distribution(0, 0, 0.25), "qtl_var")
  expect_error(plot_qtl_genotype_distribution(0.5, Inf, 0.25), "tau")
  expect_error(plot_qtl_genotype_distribution(0.5, 0, 1), "pd")
  expect_error(plot_qtl_genotype_distribution(0.5, 0, 0.25, type = "bars"), "arg")
  expect_error(plot_qtl_genotype_distribution(0.5, 0, 0.25, scale = "counts"), "arg")
  expect_error(
    plot_qtl_genotype_distribution(0.5, 0, 0.25, type = "histogram", n = 1.5),
    "positive integer"
  )
  expect_error(
    plot_qtl_genotype_distribution(0.5, 0, 0.25, type = "histogram", seed = -1),
    "seed"
  )
})
