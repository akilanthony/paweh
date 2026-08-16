mv_plot3d_args <- list(
  qtl_var = c(0.95, 0.92),
  tau = c(0, 0.5),
  pd = 0.5,
  cor_matrix = matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE),
  x_upper = c(10, 10),
  x_lower = c(15, 15)
)

test_that("3D density surface returns a Plotly widget with validated data", {
  skip_if_not_installed("plotly")
  p <- do.call(plot_qtl_multivariate_surface3d, c(
    mv_plot3d_args, list(surface = "density", grid_n = 20)
  ))
  dat <- attr(p, "plot_data")
  model <- attr(p, "falconer_model")
  direct <- do.call(pawh:::.falconer_mv_parameters, mv_plot3d_args[1:4])

  expect_s3_class(p, "plotly")
  expect_s3_class(p, "htmlwidget")
  expect_equal(model$mean_matrix, direct$mean_matrix)
  expect_equal(model$residual_covariance_matrix, direct$residual_covariance_matrix)
  expect_equal(model$genotype_frequencies, direct$genotype_frequencies)
  expect_equal(nrow(dat), 20^2)
  expect_true(all(is.finite(dat$value)))
  expect_true(all(dat$value >= 0))
  expect_equal(dat$z_value, dat$value)

  built <- plotly::plotly_build(p)
  expect_true(any(vapply(built$x$data, `[[`, character(1), "type") == "surface"))
})

test_that("3D CDF surface is bounded and monotone", {
  skip_if_not_installed("plotly")
  grid_n <- 22L
  p <- do.call(plot_qtl_multivariate_surface3d, c(
    mv_plot3d_args, list(surface = "cdf", grid_n = grid_n)
  ))
  dat <- attr(p, "plot_data")
  z <- matrix(dat$value, nrow = grid_n, ncol = grid_n)

  expect_true(all(is.finite(dat$value)))
  expect_true(all(dat$value >= 0 & dat$value <= 1))
  expect_true(all(apply(z, 2, function(column) diff(column) >= -1e-10)))
  expect_true(all(apply(z, 1, function(row) diff(row) >= -1e-10)))
})

test_that("3D and 2D views use identical surface grids and thresholds", {
  skip_if_not_installed("plotly")
  p3d <- do.call(plot_qtl_multivariate_surface3d, c(
    mv_plot3d_args, list(surface = "density", grid_n = 18)
  ))
  d2 <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot3d_args,
    list(surface = "density", grid_n = 18, return_data = TRUE)
  ))
  threshold <- attr(p3d, "thresholds")

  expect_equal(attr(p3d, "plot_data")$value, d2$value)
  expect_equal(threshold, attr(d2, "thresholds"))
  expect_equal(threshold$upper_threshold, rep(stats::qnorm(0.9), 2))
  expect_equal(threshold$lower_threshold, rep(stats::qnorm(0.15), 2))
  expect_named(attr(p3d, "threshold_overlays"), c("Affected", "Unaffected"))

  affected <- attr(p3d, "threshold_overlays")$Affected
  unaffected <- attr(p3d, "threshold_overlays")$Unaffected
  expect_true(all(affected$x >= threshold$upper_threshold[1]))
  expect_true(all(affected$y >= threshold$upper_threshold[2]))
  expect_true(all(unaffected$x <= threshold$lower_threshold[1]))
  expect_true(all(unaffected$y <= threshold$lower_threshold[2]))
  expect_equal(affected$z, rep(0, 5))
  expect_equal(unaffected$z, rep(0, 5))
})

test_that("3D mean markers and threshold overlays can be disabled", {
  skip_if_not_installed("plotly")
  full <- do.call(plot_qtl_multivariate_surface3d, c(
    mv_plot3d_args, list(surface = "cdf", grid_n = 18)
  ))
  plain <- do.call(plot_qtl_multivariate_surface3d, c(
    mv_plot3d_args,
    list(surface = "cdf", grid_n = 18, show_means = FALSE,
         show_thresholds = FALSE, show_labels = FALSE)
  ))

  expect_equal(nrow(attr(full, "mean_markers")), 3)
  expect_true(all(c("Genotype 0", "Genotype 1", "Genotype 2") %in%
                    attr(full, "mean_markers")$genotype))
  expect_null(attr(plain, "mean_markers"))
  expect_null(attr(plain, "threshold_overlays"))
  expect_s3_class(plain, "plotly")
})

test_that("normalized z scaling preserves raw values", {
  skip_if_not_installed("plotly")
  p <- do.call(plot_qtl_multivariate_surface3d, c(
    mv_plot3d_args,
    list(surface = "density", grid_n = 18, z_scale = "normalized")
  ))
  dat <- attr(p, "plot_data")

  expect_equal(range(dat$z_value), c(0, 1), tolerance = 1e-12)
  expect_true(all(dat$value >= 0))
  expect_identical(attr(p, "z_scale"), "normalized")
  expect_true(all(attr(p, "mean_markers")$z_value >= 0 &
                    attr(p, "mean_markers")$z_value <= 1))
})

test_that("3D surface validates its inputs", {
  skip_if_not_installed("plotly")
  expect_error(
    plot_qtl_multivariate_surface3d(
      qtl_var = 0.5, tau = 0, pd = 0.5, cor_matrix = matrix(1)
    ),
    "exactly two"
  )
  expect_error(
    plot_qtl_multivariate_surface3d(
      qtl_var = c(0.5, 0.4), tau = c(0, 0), pd = 0.5,
      cor_matrix = diag(3)
    ),
    "2 x 2"
  )
  expect_error(
    plot_qtl_multivariate_surface3d(
      qtl_var = c(0.5, 0.4), tau = c(0, 0), pd = 0.5,
      cor_matrix = matrix(c(1, 1.2, 1.2, 1), 2)
    ),
    "correlations"
  )
  expect_error(
    do.call(plot_qtl_multivariate_surface3d, c(
      mv_plot3d_args, list(surface = "height")
    )),
    "arg"
  )
  expect_error(
    do.call(plot_qtl_multivariate_surface3d, c(
      mv_plot3d_args, list(z_scale = "log")
    )),
    "arg"
  )
  expect_error(
    do.call(plot_qtl_multivariate_surface3d, c(
      mv_plot3d_args[1:4], list(x_upper = c(10, 10))
    )),
    "supplied together"
  )
})

test_that("3D plotting leaves validated multivariate results unchanged", {
  skip_if_not_installed("plotly")
  backend_args <- list(
    power = 0.95, alpha = 5e-8,
    qtl_var = c(0.01, 0.005), tau = c(0, 0.5), pd = 0.25,
    cor_matrix = matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE),
    test = "threshold_chisq", x_upper = c(10, 10), x_lower = c(15, 15),
    k = 1, verbose = FALSE
  )
  before <- do.call(qtl_multivariate_mssn_full, backend_args)
  invisible(do.call(plot_qtl_multivariate_surface3d, c(
    mv_plot3d_args, list(surface = "cdf", grid_n = 16)
  )))
  after <- do.call(qtl_multivariate_mssn_full, backend_args)

  expect_equal(unclass(after), unclass(before))
  expect_equal(after$historical_fractional_cases, 448.804484, tolerance = 1e-6)
  expect_identical(after$N_case, 449)
  expect_identical(after$N_control, 449)
  expect_identical(after$N_total, 898)
})
