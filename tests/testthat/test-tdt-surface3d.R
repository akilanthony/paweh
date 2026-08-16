test_that("generalized TDT surfaces reproduce the four prior conventions", {
  skip_if_not_installed("plotly")

  common <- list(
    x = "pd", x_values = c(0.2, 0.3),
    prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 1
  )

  power_misclassification <- do.call(
    plot_tdt_surface3d,
    c(common, list(
      metric = "power", scenario = "misclassification",
      y = "misclass_rate", y_values = c(0, 0.02),
      N = 300, heter_rate = 0
    ))
  )
  power_heterogeneity <- do.call(
    plot_tdt_surface3d,
    c(common, list(
      metric = "power", scenario = "heterogeneity",
      y = "heter_rate", y_values = c(0, 0.1),
      N = 300, misclass_rate = 0
    ))
  )
  mssn_misclassification <- do.call(
    plot_tdt_surface3d,
    c(common, list(
      metric = "mssn", scenario = "misclassification",
      y = "misclass_rate", y_values = c(0, 0.02),
      target_power = 0.8, heter_rate = 0
    ))
  )
  mssn_heterogeneity <- do.call(
    plot_tdt_surface3d,
    c(common, list(
      metric = "mssn", scenario = "heterogeneity",
      y = "heter_rate", y_values = c(0, 0.1),
      target_power = 0.8, misclass_rate = 0
    ))
  )

  expect_equal(
    attr(power_misclassification, "surface_matrix"),
    matrix(c(
      0.839945430684816, 0.569641846610049,
      0.911508069379948, 0.662179072750112
    ), nrow = 2),
    tolerance = 1e-12
  )
  expect_equal(
    attr(power_heterogeneity, "surface_matrix"),
    matrix(c(
      0.839945430684816, 0.757667235233636,
      0.911508069379948, 0.845921067355492
    ), nrow = 2),
    tolerance = 1e-12
  )
  expect_equal(
    attr(mssn_misclassification, "surface_matrix"),
    matrix(c(270, 517, 215, 417), nrow = 2)
  )
  expect_equal(
    attr(mssn_heterogeneity, "surface_matrix"),
    matrix(c(270, 334, 215, 266), nrow = 2)
  )
})


test_that("surface scalar extraction maps canonical result components", {
  result <- list(
    power = list(no_error = 0.91, misclassification = 0.73),
    N = list(no_error = 212.4, heterogeneity = 348.8)
  )

  expect_identical(
    genmixr:::.tdt_surface_extract_value(result, "power", "no_error"),
    0.91
  )
  expect_identical(
    genmixr:::.tdt_surface_extract_value(result, "power", "misclassification"),
    0.73
  )
  expect_identical(
    genmixr:::.tdt_surface_extract_value(result, "mssn", "no_error"),
    212.4
  )
  expect_identical(
    genmixr:::.tdt_surface_extract_value(result, "mssn", "heterogeneity"),
    348.8
  )
})


test_that("surface supports generalized axes and preserves Cartesian order", {
  skip_if_not_installed("plotly")

  plot <- plot_tdt_surface3d(
    metric = "power", scenario = "no_error",
    x = "prev", y = "R1",
    x_values = c(0.02, 0.05, 0.08),
    y_values = c(1.2, 1.6),
    N = 350, pd = 0.25, R2 = 2.1
  )
  data <- attr(plot, "surface_data")

  expect_s3_class(plot, "plotly")
  expect_s3_class(plot, "htmlwidget")
  expect_named(
    data,
    c("prev", "R1", "metric_value", "raw_metric_value", "metric", "scenario")
  )
  expect_equal(data$prev, rep(c(0.02, 0.05, 0.08), 2))
  expect_equal(data$R1, rep(c(1.2, 1.6), each = 3))
  expect_identical(dim(attr(plot, "surface_matrix")), c(2L, 3L))
  expect_equal(data$metric_value, data$raw_metric_value)
})


test_that("surface forwards fixed and swept arguments to canonical backends", {
  skip_if_not_installed("plotly")

  plot <- plot_tdt_surface3d(
    metric = "power", scenario = "misclassification",
    x = "alpha", y = "misclass_rate",
    x_values = c(0.01, 0.04), y_values = c(0, 0.03),
    N = 425, pd = 0.27, prev = 0.04, R1 = 1.4, R2 = 2.05,
    delta_prime = 0.8, heter_rate = 0.15
  )
  data <- attr(plot, "surface_data")
  direct <- tdt_power(
    N = 425, input_mode = "model_based",
    pd = 0.27, prev = 0.04, R1 = 1.4, R2 = 2.05,
    alpha = data$alpha[[4]], delta_prime = 0.8,
    misclass_rate = data$misclass_rate[[4]], heter_rate = 0.15,
    verbose = FALSE
  )$power$misclassification

  expect_equal(data$raw_metric_value[[4]], direct)
})


test_that("new model-parameter axis combinations produce finite grids", {
  skip_if_not_installed("plotly")

  designs <- list(
    list(
      metric = "power", scenario = "no_error",
      x = "R1", y = "R2",
      x_values = c(1.2, 1.5), y_values = c(1.8, 2.2)
    ),
    list(
      metric = "mssn", scenario = "no_error",
      x = "R1", y = "R2",
      x_values = c(1.2, 1.5), y_values = c(1.8, 2.2)
    ),
    list(
      metric = "power", scenario = "no_error",
      x = "pd", y = "delta_prime",
      x_values = c(0.2, 0.35), y_values = c(0.5, 0.9)
    ),
    list(
      metric = "mssn", scenario = "no_error",
      x = "pd", y = "delta_prime",
      x_values = c(0.2, 0.35), y_values = c(0.5, 0.9)
    )
  )

  for (design in designs) {
    plot <- do.call(plot_tdt_surface3d, design)
    values <- attr(plot, "surface_data")$metric_value
    expect_length(values, 4)
    expect_true(all(is.finite(values)))
    expect_gt(length(unique(values)), 1)
    expect_identical(dim(attr(plot, "surface_matrix")), c(2L, 2L))
  }
})


test_that("modifier surfaces have expected local monotonic direction", {
  skip_if_not_installed("plotly")

  modifier_surface <- function(metric, scenario, y, y_values) {
    plot_tdt_surface3d(
      metric = metric, scenario = scenario,
      x = "pd", y = y,
      x_values = c(0.25, 0.35), y_values = y_values,
      N = 400, target_power = 0.8, ceiling_N = FALSE
    )
  }

  power_misclassification <- modifier_surface(
    "power", "misclassification", "misclass_rate", c(0, 0.03, 0.06)
  )
  mssn_misclassification <- modifier_surface(
    "mssn", "misclassification", "misclass_rate", c(0, 0.03, 0.06)
  )
  power_heterogeneity <- modifier_surface(
    "power", "heterogeneity", "heter_rate", c(0, 0.1, 0.2)
  )
  mssn_heterogeneity <- modifier_surface(
    "mssn", "heterogeneity", "heter_rate", c(0, 0.1, 0.2)
  )

  expect_true(all(apply(
    attr(power_misclassification, "surface_matrix"), 2, function(z) diff(z) <= 0
  )))
  expect_true(all(apply(
    attr(mssn_misclassification, "surface_matrix"), 2, function(z) diff(z) >= 0
  )))
  expect_true(all(apply(
    attr(power_heterogeneity, "surface_matrix"), 2, function(z) diff(z) <= 0
  )))
  expect_true(all(apply(
    attr(mssn_heterogeneity, "surface_matrix"), 2, function(z) diff(z) >= 0
  )))
})


test_that("MSSN ceiling control retains raw canonical results", {
  skip_if_not_installed("plotly")

  raw_plot <- plot_tdt_surface3d(
    metric = "mssn", scenario = "heterogeneity",
    x = "delta_prime", y = "heter_rate",
    x_values = c(0.7, 0.9), y_values = c(0.05, 0.2),
    target_power = 0.85, ceiling_N = FALSE
  )
  ceiling_plot <- plot_tdt_surface3d(
    metric = "mssn", scenario = "heterogeneity",
    x = "delta_prime", y = "heter_rate",
    x_values = c(0.7, 0.9), y_values = c(0.05, 0.2),
    target_power = 0.85, ceiling_N = TRUE
  )
  raw_data <- attr(raw_plot, "surface_data")
  ceiling_data <- attr(ceiling_plot, "surface_data")

  expect_equal(raw_data$metric_value, raw_data$raw_metric_value)
  expect_equal(ceiling_data$raw_metric_value, raw_data$raw_metric_value)
  expect_equal(ceiling_data$metric_value, ceiling(raw_data$raw_metric_value))
})


test_that("parameter-specific default grids are used", {
  skip_if_not_installed("plotly")

  plot <- plot_tdt_surface3d(
    metric = "power", scenario = "misclassification",
    x = "pd", y = "misclass_rate"
  )
  spec <- attr(plot, "surface_spec")

  expect_equal(spec$x_values, seq(0.10, 0.50, length.out = 20))
  expect_equal(spec$y_values, seq(0, 0.20, length.out = 20))
  expect_identical(dim(attr(plot, "surface_matrix")), c(20L, 20L))
})


test_that("surface metadata and hover labels identify the selected design", {
  skip_if_not_installed("plotly")

  plot <- plot_tdt_surface3d(
    metric = "mssn", scenario = "heterogeneity",
    x = "prev", y = "heter_rate",
    x_values = c(0.03, 0.05), y_values = c(0, 0.2)
  )
  spec <- attr(plot, "surface_spec")
  plot_attrs <- plot$x$attrs[[1L]]
  layout <- plot$x$layoutAttrs[[1L]]$scene

  expect_identical(spec$metric, "mssn")
  expect_identical(spec$scenario, "heterogeneity")
  expect_match(plot_attrs$hovertemplate, "Disease prevalence", fixed = TRUE)
  expect_match(plot_attrs$hovertemplate, "MSSN (affected trios)", fixed = TRUE)
  expect_identical(layout$xaxis$title, "Disease prevalence")
  expect_identical(layout$yaxis$title, "Locus heterogeneity rate (1 - pi)")
  expect_identical(layout$zaxis$title, "MSSN (affected trios)")
})


test_that("surface validation rejects ambiguous or inactive designs", {
  skip_if_not_installed("plotly")

  expect_error(plot_tdt_surface3d(metric = "risk"), "arg")
  expect_error(plot_tdt_surface3d(scenario = "combined"), "arg")
  expect_error(
    plot_tdt_surface3d(x = "pd", y = "pd", x_values = 1:2, y_values = 1:2),
    "distinct"
  )
  expect_error(
    plot_tdt_surface3d(x = "N", y = "pd"),
    "x must be one of"
  )
  expect_error(
    plot_tdt_surface3d(
      scenario = "misclassification", x = "pd", y = "heter_rate"
    ),
    "inactive"
  )
  expect_error(
    plot_tdt_surface3d(
      scenario = "no_error", x = "pd", y = "misclass_rate"
    ),
    "inactive"
  )
  expect_error(
    plot_tdt_surface3d(x_values = 0.2, y_values = c(0, 0.1)),
    "x_values must contain at least two"
  )
  expect_error(
    plot_tdt_surface3d(x_values = c(0.2, NA), y_values = c(0, 0.1)),
    "finite numeric"
  )
  expect_error(
    plot_tdt_surface3d(x_values = c(0.2, 0.2), y_values = c(0, 0.1)),
    "must be distinct"
  )
  expect_error(
    plot_tdt_surface3d(x_values = c(0, 0.2), y_values = c(0, 0.1)),
    "outside its valid range"
  )
  expect_error(
    plot_tdt_surface3d(
      x_values = c(0.2, 0.3), y_values = c(0, 0.1), prev = 1
    ),
    "prev has an invalid fixed value"
  )
  expect_error(
    plot_tdt_surface3d(
      metric = "mssn", x_values = c(0.2, 0.3),
      y_values = c(0, 0.1), target_power = 1
    ),
    "target_power has an invalid fixed value"
  )
  expect_error(
    plot_tdt_surface3d(
      x_values = c(0.2, 0.3), y_values = c(0, 0.1), ceiling_N = NA
    ),
    "ceiling_N"
  )
})


test_that("backend errors identify the failing grid point", {
  skip_if_not_installed("plotly")

  expect_error(
    plot_tdt_surface3d(
      metric = "mssn", scenario = "no_error",
      x = "pd", y = "prev",
      x_values = c(0.2, 0.3), y_values = c(0.8, 0.9),
      R1 = 10, R2 = 20
    ),
    "TDT surface evaluation failed at pd = 0.2, prev = 0.8"
  )
})


test_that("surface reports a clear optional dependency error", {
  testthat::local_mocked_bindings(
    .tdt_surface_plotly_available = function() FALSE,
    .package = "genmixr"
  )
  expect_error(
    plot_tdt_surface3d(),
    "Package 'plotly' is required for this plot.",
    fixed = TRUE
  )
})
