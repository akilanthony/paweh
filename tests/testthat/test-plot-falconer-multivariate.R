mv_plot_args <- list(
  qtl_var = c(0.95, 0.92),
  tau = c(0, 0.5),
  pd = 0.5,
  cor_matrix = matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE),
  x_upper = c(10, 10),
  x_lower = c(15, 15)
)

test_that("multivariate density plot reuses validated model quantities", {
  dat <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args, list(surface = "density", grid_n = 20, return_data = TRUE)
  ))
  model <- attr(dat, "falconer_model")
  direct <- do.call(pawh:::.falconer_mv_parameters, mv_plot_args[1:4])

  expect_s3_class(dat, "data.frame")
  expect_equal(model$mean_matrix, direct$mean_matrix)
  expect_equal(model$residual_covariance_matrix, direct$residual_covariance_matrix)
  expect_equal(model$genotype_frequencies, direct$genotype_frequencies)
  expect_equal(sum(model$genotype_frequencies), 1)
  expect_true(all(is.finite(dat$value)))
  expect_true(all(dat$value >= 0))
  expect_equal(nrow(dat), 20^2)
  expect_equal(
    dat$value,
    as.numeric(attr(dat, "component_values") %*% model$genotype_frequencies)
  )
})

test_that("genotype-density data are unweighted conditional distributions", {
  grid_n <- 61L
  dat <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args,
    list(surface = "genotype_density", grid_n = grid_n, return_data = TRUE)
  ))
  model <- attr(dat, "falconer_model")

  expect_equal(nrow(dat), 3 * grid_n^2)
  expect_identical(levels(dat$genotype), paste("Genotype", 0:2))
  expect_equal(as.numeric(sort(table(dat$genotype))), rep(grid_n^2, 3))
  expect_equal(dat$value, dat$conditional_density)
  expect_true(all(is.finite(dat$conditional_density)))
  expect_true(all(dat$conditional_density >= 0))

  for (j in 1:3) {
    component <- dat[dat$genotype == paste("Genotype", j - 1L), ]
    direct <- mvtnorm::dmvnorm(
      cbind(component$phenotype_1, component$phenotype_2),
      mean = model$mean_matrix[, j],
      sigma = model$residual_covariance_matrix
    )
    expect_equal(component$conditional_density, direct)

    dx <- diff(sort(unique(component$phenotype_1)))[1]
    dy <- diff(sort(unique(component$phenotype_2)))[1]
    expect_equal(sum(component$conditional_density) * dx * dy, 1,
                 tolerance = 0.006)
  }
})

test_that("conditional densities locate genotype means with correct axes", {
  dat <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args[1:4],
    list(surface = "genotype_density", grid_n = 75, return_data = TRUE)
  ))
  model <- attr(dat, "falconer_model")

  for (j in 1:3) {
    component <- dat[dat$genotype == paste("Genotype", j - 1L), ]
    nearest <- which.min(
      (component$phenotype_1 - model$mean_matrix[1, j])^2 +
        (component$phenotype_2 - model$mean_matrix[2, j])^2
    )
    direct <- mvtnorm::dmvnorm(
      c(component$phenotype_1[nearest], component$phenotype_2[nearest]),
      mean = model$mean_matrix[, j],
      sigma = model$residual_covariance_matrix
    )
    expect_equal(component$conditional_density[nearest], direct)
    expect_gte(component$conditional_density[nearest],
               0.97 * max(component$conditional_density))
  }
})

test_that("rare and overlapping genotypes remain separate conditional curves", {
  rare_args <- list(
    qtl_var = c(0.10, 0.05), tau = c(0, 0.5), pd = 0.01,
    cor_matrix = matrix(c(1, 0.4, 0.4, 1), 2, byrow = TRUE)
  )
  dat <- do.call(plot_qtl_multivariate_contour, c(
    rare_args,
    list(surface = "genotype_density", grid_n = 40, return_data = TRUE)
  ))
  model <- attr(dat, "falconer_model")

  expect_equal(unname(model$genotype_frequencies), c(0.9801, 0.0198, 0.0001))
  expect_identical(levels(dat$genotype), paste("Genotype", 0:2))
  expect_equal(as.numeric(table(dat$genotype)), rep(40^2, 3))
  expect_true(all(is.finite(dat$conditional_density)))
  expect_true(all(dat$conditional_density >= 0))

  peak_at_mean <- vapply(1:3, function(j) {
    mvtnorm::dmvnorm(
      model$mean_matrix[, j], mean = model$mean_matrix[, j],
      sigma = model$residual_covariance_matrix
    )
  }, numeric(1))
  expect_equal(peak_at_mean, rep(peak_at_mean[1], 3))
})

test_that("near-fixed alleles and coincident means keep three finite components", {
  dat <- plot_qtl_multivariate_contour(
    qtl_var = c(0.10, 0.08), tau = c(1, -1), pd = 0.99,
    cor_matrix = matrix(c(1, -0.35, -0.35, 1), 2, byrow = TRUE),
    surface = "genotype_density", grid_n = 30, return_data = TRUE
  )
  model <- attr(dat, "falconer_model")

  expect_equal(as.numeric(table(dat$genotype)), rep(30^2, 3))
  expect_true(all(is.finite(dat$conditional_density)))
  expect_true(all(dat$conditional_density >= 0))
  expect_silent(chol(model$residual_covariance_matrix))
  expect_equal(model$mean_matrix[1, 2], model$mean_matrix[1, 3])
  expect_equal(model$mean_matrix[2, 1], model$mean_matrix[2, 2])
})

test_that("genotype-density mode retains joint threshold classifications", {
  dat <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args,
    list(surface = "genotype_density", grid_n = 25, return_data = TRUE)
  ))
  threshold <- attr(dat, "thresholds")

  expect_true(all(dat$affected ==
    (dat$phenotype_1 >= threshold$upper_threshold[1] &
       dat$phenotype_2 >= threshold$upper_threshold[2])))
  expect_true(all(dat$unaffected ==
    (dat$phenotype_1 <= threshold$lower_threshold[1] &
       dat$phenotype_2 <= threshold$lower_threshold[2])))
  expect_false(any(dat$affected & dat$unaffected))
})

test_that("multivariate CDF is bounded and agrees with direct probabilities", {
  dat <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args, list(surface = "cdf", grid_n = 18, return_data = TRUE)
  ))
  model <- attr(dat, "falconer_model")

  expect_true(all(is.finite(dat$value)))
  expect_true(all(dat$value >= 0 & dat$value <= 1))
  i <- which.min(dat$phenotype_1^2 + dat$phenotype_2^2)
  direct <- sum(vapply(1:3, function(j) {
    model$genotype_frequencies[j] * as.numeric(mvtnorm::pmvnorm(
      lower = c(-Inf, -Inf),
      upper = c(dat$phenotype_1[i], dat$phenotype_2[i]),
      mean = model$mean_matrix[, j],
      sigma = model$residual_covariance_matrix,
      algorithm = mvtnorm::Miwa(steps = 512L)
    ))
  }, numeric(1)))
  expect_equal(dat$value[i], direct, tolerance = 2e-5)
})

test_that("Chapter 6.2 thresholds and joint-AND regions are explicit", {
  dat <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args, list(surface = "density", grid_n = 25, return_data = TRUE)
  ))
  threshold <- attr(dat, "thresholds")

  expect_equal(threshold$upper_threshold, rep(stats::qnorm(0.9), 2))
  expect_equal(threshold$lower_threshold, rep(stats::qnorm(0.15), 2))
  expect_true(all(
    dat$phenotype_1[dat$affected] >= threshold$upper_threshold[1] &
      dat$phenotype_2[dat$affected] >= threshold$upper_threshold[2]
  ))
  expect_true(all(
    dat$phenotype_1[dat$unaffected] <= threshold$lower_threshold[1] &
      dat$phenotype_2[dat$unaffected] <= threshold$lower_threshold[2]
  ))
  expect_false(any(dat$affected & dat$unaffected))
  expect_true(any(dat$selection == "Not selected"))
})

test_that("multivariate contour modes return customizable plots", {
  density_plot <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args, list(surface = "density", grid_n = 20)
  ))
  cdf_plot <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args, list(surface = "cdf", grid_n = 20)
  ))
  plain_plot <- do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args[1:4],
    list(surface = "density", grid_n = 20, show_thresholds = FALSE,
         show_means = FALSE, show_labels = FALSE)
  ))

  expect_s3_class(density_plot, "ggplot")
  expect_s3_class(cdf_plot, "ggplot")
  expect_s3_class(cdf_plot + ggplot2::theme_minimal(), "ggplot")
  expect_null(attr(plain_plot, "thresholds"))
  expect_false(any(vapply(
    plain_plot$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)
  )))
})

test_that("multivariate contour inputs are validated", {
  expect_error(
    plot_qtl_multivariate_contour(
      qtl_var = 0.5, tau = 0, pd = 0.5, cor_matrix = matrix(1)
    ),
    "exactly two"
  )
  expect_error(
    plot_qtl_multivariate_contour(
      qtl_var = c(0.5, 0.4), tau = c(0, 0), pd = 0.5,
      cor_matrix = diag(3)
    ),
    "2 x 2"
  )
  expect_error(
    plot_qtl_multivariate_contour(
      qtl_var = c(0.5, 0.4), tau = c(0, 0), pd = 0.5,
      cor_matrix = matrix(c(1, 1.2, 1.2, 1), 2)
    ),
    "correlations"
  )
  expect_error(
    do.call(plot_qtl_multivariate_contour, c(
      mv_plot_args[1:4], list(x_upper = c(10, 10))
    )),
    "supplied together"
  )
  expect_error(
    do.call(plot_qtl_multivariate_contour, c(mv_plot_args, list(grid_n = 5))),
    "at least 10"
  )
})

test_that("plotting leaves validated Chapter 6.2 results unchanged", {
  backend_args <- list(
    power = 0.95, alpha = 5e-8,
    qtl_var = c(0.01, 0.005), tau = c(0, 0.5), pd = 0.25,
    cor_matrix = matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE),
    test = "threshold_chisq", x_upper = c(10, 10), x_lower = c(15, 15),
    k = 1, verbose = FALSE
  )
  before <- do.call(qtl_multivariate_mssn_full, backend_args)
  invisible(do.call(plot_qtl_multivariate_contour, c(
    mv_plot_args, list(surface = "cdf", grid_n = 15)
  )))
  after <- do.call(qtl_multivariate_mssn_full, backend_args)

  expect_equal(unclass(after), unclass(before))
  expect_equal(after$historical_fractional_cases, 448.804484, tolerance = 1e-6)
  expect_identical(after$N_case, 449)
  expect_identical(after$N_control, 449)
  expect_identical(after$N_total, 898)
})
