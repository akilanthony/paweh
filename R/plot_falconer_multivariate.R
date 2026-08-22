# Bivariate Falconer mixture and genotype-conditional visualization.

#' Evaluate one bivariate-normal density over plotting points.
#'
#' @noRd
.falconer_mv_plot_density <- function(points, mean, Sigma) {
  mvtnorm::dmvnorm(
    cbind(points$phenotype_1, points$phenotype_2),
    mean = mean, sigma = Sigma
  )
}

#' Evaluate one bivariate-normal CDF over plotting points.
#'
#' This uses mvtnorm's vectorized Genz transformation. In two dimensions the
#' integration is one-dimensional; shared deterministic midpoint weights give
#' a smooth, reproducible plotting surface without calling pmvnorm separately
#' for every grid cell.
#'
#' @noRd
.falconer_mv_plot_cdf <- function(points, mean, Sigma, M = 1024L) {
  upper <- rbind(points$phenotype_1, points$phenotype_2)
  lower <- matrix(-Inf, nrow = 2L, ncol = nrow(points))
  L <- t(chol(Sigma))
  packed_chol <- mvtnorm::ltMatrices(
    L[lower.tri(L, diag = TRUE)], diag = TRUE, byrow = FALSE
  )
  midpoint_weights <- matrix((seq_len(M) - 0.5) / M, nrow = 1L)
  log_probability <- mvtnorm::lpmvnorm(
    lower = lower, upper = upper, mean = mean, chol = packed_chol,
    logLik = FALSE, M = M, w = midpoint_weights
  )
  pmin(1, pmax(0, exp(log_probability)))
}

#' Plot Two-Phenotype Falconer Density or CDF Contours
#'
#' Visualizes three genotype-specific bivariate normal distributions for exactly
#' two quantitative phenotypes. It can show the marginal mixture density, the
#' mixture lower-tail CDF, or the three genotype-conditional density contour
#' families, together with the joint threshold regions used by the Chapter 6.2
#' selected-sampling design.
#'
#' @param qtl_var Numeric vector of two phenotype-specific QTL variances.
#' @param tau Numeric vector of two phenotype-specific dominance/additivity
#'   ratios.
#' @param pd Shared increaser-allele frequency in `(0, 1)`.
#' @param cor_matrix A positive-definite `2 x 2` phenotype correlation matrix.
#' @param x_upper,x_lower Optional length-two vectors of upper- and lower-tail
#'   selection percentages. Supply both or neither.
#' @param surface One of `"density"`, `"genotype_density"`, or `"cdf"`.
#' @param show_thresholds Logical; display joint threshold regions when
#'   thresholds are supplied.
#' @param show_means Logical; mark the three genotype-specific mean vectors.
#' @param show_labels Logical; label the affected and unaffected joint regions.
#' @param grid_n Integer number of grid points along each phenotype axis.
#' @param title Optional plot title.
#' @param return_data Logical; return plotting data instead of the ggplot. The
#'   validated model and threshold details are retained as attributes. Mixture
#'   density and CDF modes retain their existing one-row-per-grid-point form.
#'   Genotype-density mode returns long-form data with one row per grid point
#'   and genotype.
#'
#' @details
#' In density mode, each grid value is the marginal, genotype-weighted mixture
#' `sum(pi[j] * f[j](y1, y2))`, where each `f[j]` is a bivariate-normal
#' density. In CDF mode, it is the lower-tail mixture probability
#' `sum(pi[j] * P(Y1 <= y1, Y2 <= y2 | G = j))`. Density and CDF surfaces
#' therefore represent different mathematical quantities.
#'
#' Genotype-density mode displays each conditional density `f[j](y1, y2)`
#' separately and does not multiply by genotype frequency. Each conditional
#' density integrates to one, including for rare genotypes. Three genotype
#' components do not imply that their weighted mixture has three distinct
#' modes.
#'
#' When thresholds are supplied, affected subjects occupy only the joint
#' upper-right region `Y1 >= TU1 AND Y2 >= TU2`. Unaffected subjects occupy
#' only the joint lower-left region `Y1 <= TL1 AND Y2 <= TL2`. Subjects in all
#' other regions are not selected. The statistical backend can support more
#' traits where documented; this conventional contour visualization is
#' intentionally restricted to two.
#'
#' @return A \code{ggplot} object, or plotting data when `return_data = TRUE`.
#'   Genotype-density data contain `phenotype_1`, `phenotype_2`, `genotype`,
#'   `conditional_density`, and `value` (equal to `conditional_density`).
#'   Returned data retain the validated model and threshold details as
#'   attributes.
#'
#' @examples
#' cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
#' plot_qtl_multivariate_contour(
#'   qtl_var = c(0.95, 0.92), tau = c(0, 0.5), pd = 0.5,
#'   cor_matrix = cor_matrix,
#'   x_upper = c(10, 10), x_lower = c(15, 15),
#'   surface = "cdf", grid_n = 30
#' )
#'
#' @seealso [qtl_multivariate_power_full()],
#' [qtl_multivariate_mssn_full()], and
#' [plot_qtl_multivariate_surface3d()].
#'
#' @export
plot_qtl_multivariate_contour <- function(
    qtl_var,
    tau,
    pd,
    cor_matrix,
    x_upper = NULL,
    x_lower = NULL,
    surface = c("density", "genotype_density", "cdf"),
    show_thresholds = TRUE,
    show_means = TRUE,
    show_labels = TRUE,
    grid_n = 150L,
    title = NULL,
    return_data = FALSE
) {
  surface <- match.arg(surface)
  .falconer_check_flag(show_thresholds, "show_thresholds")
  .falconer_check_flag(show_means, "show_means")
  .falconer_check_flag(show_labels, "show_labels")
  .falconer_check_flag(return_data, "return_data")
  if (length(qtl_var) != 2L || length(tau) != 2L) {
    stop("plot_qtl_multivariate_contour() requires exactly two phenotypes.", call. = FALSE)
  }
  if (!is.matrix(cor_matrix) || !identical(dim(cor_matrix), c(2L, 2L))) {
    stop("cor_matrix must be a 2 x 2 matrix for this visualization.", call. = FALSE)
  }
  if (!is.numeric(grid_n) || length(grid_n) != 1L || !is.finite(grid_n) ||
      grid_n != floor(grid_n) || grid_n < 10L) {
    stop("grid_n must be an integer of at least 10.", call. = FALSE)
  }
  if (xor(is.null(x_upper), is.null(x_lower))) {
    stop("x_upper and x_lower must be supplied together.", call. = FALSE)
  }

  model <- .falconer_mv_parameters(qtl_var, tau, pd, cor_matrix)
  threshold <- if (is.null(x_upper)) NULL else {
    .falconer_mv_threshold_components(model, x_upper, x_lower)
  }
  axis_sd <- sqrt(diag(model$residual_covariance_matrix))
  axis_min <- apply(model$mean_matrix, 1L, min) - 4 * axis_sd
  axis_max <- apply(model$mean_matrix, 1L, max) + 4 * axis_sd
  if (!is.null(threshold)) {
    axis_min <- pmin(axis_min, threshold$lower_threshold - 0.5 * axis_sd)
    axis_max <- pmax(axis_max, threshold$upper_threshold + 0.5 * axis_sd)
  }
  grid <- expand.grid(
    phenotype_1 = seq(axis_min[1L], axis_max[1L], length.out = as.integer(grid_n)),
    phenotype_2 = seq(axis_min[2L], axis_max[2L], length.out = as.integer(grid_n))
  )

  component <- vapply(1:3, function(j) {
    if (surface %in% c("density", "genotype_density")) {
      .falconer_mv_plot_density(
        grid, model$mean_matrix[, j], model$residual_covariance_matrix
      )
    } else {
      .falconer_mv_plot_cdf(
        grid, model$mean_matrix[, j], model$residual_covariance_matrix
      )
    }
  }, numeric(nrow(grid)))
  colnames(component) <- paste0("genotype_", 0:2)
  grid$value <- as.numeric(component %*% model$genotype_frequencies)
  grid$surface <- surface
  if (!is.null(threshold)) {
    grid$affected <- grid$phenotype_1 >= threshold$upper_threshold[1L] &
      grid$phenotype_2 >= threshold$upper_threshold[2L]
    grid$unaffected <- grid$phenotype_1 <= threshold$lower_threshold[1L] &
      grid$phenotype_2 <= threshold$lower_threshold[2L]
    grid$selection <- ifelse(
      grid$affected, "Affected",
      ifelse(grid$unaffected, "Unaffected", "Not selected")
    )
  }
  attr(grid, "falconer_model") <- model
  attr(grid, "thresholds") <- threshold
  attr(grid, "component_values") <- component

  plot_grid <- grid
  if (identical(surface, "genotype_density")) {
    genotype_labels <- paste("Genotype", 0:2)
    plot_grid <- do.call(rbind, lapply(seq_along(genotype_labels), function(j) {
      out <- grid
      out$genotype <- factor(genotype_labels[j], levels = genotype_labels)
      out$conditional_density <- component[, j]
      out$value <- out$conditional_density
      out
    }))
    rownames(plot_grid) <- NULL
    attr(plot_grid, "falconer_model") <- model
    attr(plot_grid, "thresholds") <- threshold
    attr(plot_grid, "component_values") <- component
  }
  if (isTRUE(return_data)) return(plot_grid)

  if (identical(surface, "genotype_density")) {
    palette <- c("#0072B2", "#D55E00", "#009E73")
    p <- ggplot2::ggplot(
      plot_grid,
      ggplot2::aes(
        x = .data$phenotype_1, y = .data$phenotype_2,
        z = .data$conditional_density,
        colour = .data$genotype, linetype = .data$genotype
      )
    ) +
      ggplot2::geom_contour(bins = 9L, linewidth = 0.72) +
      ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
      ggplot2::scale_linetype_manual(
        values = c("solid", "dashed", "dotdash"), drop = FALSE
      ) +
      ggplot2::labs(
        colour = "Genotype", linetype = "Genotype",
        title = if (is.null(title)) {
          "Bivariate Falconer genotype-conditional densities"
        } else title
      )
  } else {
    p <- ggplot2::ggplot(
      plot_grid,
      ggplot2::aes(x = .data$phenotype_1, y = .data$phenotype_2, z = .data$value)
    ) +
      ggplot2::geom_contour_filled(bins = 14L) +
      ggplot2::scale_fill_viridis_d(option = "C", direction = -1) +
      ggplot2::labs(
        fill = if (identical(surface, "density")) "Mixture density" else "Mixture CDF",
        title = if (is.null(title)) {
          paste("Bivariate Falconer mixture", surface)
        } else title
      )
  }
  p <- p +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(x = "Phenotype 1 value", y = "Phenotype 2 value") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  if (isTRUE(show_thresholds) && !is.null(threshold)) {
    tu <- threshold$upper_threshold
    tl <- threshold$lower_threshold
    regions <- data.frame(
      region = c("Affected", "Unaffected"),
      xmin = c(tu[1L], axis_min[1L]), xmax = c(axis_max[1L], tl[1L]),
      ymin = c(tu[2L], axis_min[2L]), ymax = c(axis_max[2L], tl[2L])
    )
    p <- p +
      ggplot2::geom_rect(
        data = regions,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ymin, ymax = .data$ymax),
        inherit.aes = FALSE, fill = c("#D55E00", "#0072B2"),
        alpha = 0.10, colour = NA
      ) +
      ggplot2::annotate(
        "segment", x = tu[1L], xend = axis_max[1L],
        y = tu[2L], yend = tu[2L], linewidth = 0.65
      ) +
      ggplot2::annotate(
        "segment", x = tu[1L], xend = tu[1L],
        y = tu[2L], yend = axis_max[2L], linewidth = 0.65
      ) +
      ggplot2::annotate(
        "segment", x = axis_min[1L], xend = tl[1L],
        y = tl[2L], yend = tl[2L], linewidth = 0.65
      ) +
      ggplot2::annotate(
        "segment", x = tl[1L], xend = tl[1L],
        y = axis_min[2L], yend = tl[2L], linewidth = 0.65
      )
    attr(p, "selection_regions") <- regions
    if (isTRUE(show_labels)) {
      labels <- data.frame(
        label = c("Affected", "Unaffected"),
        x = c(tu[1L] + 0.62 * (axis_max[1L] - tu[1L]),
              axis_min[1L] + 0.38 * (tl[1L] - axis_min[1L])),
        y = c(tu[2L] + 0.82 * (axis_max[2L] - tu[2L]),
              axis_min[2L] + 0.18 * (tl[2L] - axis_min[2L]))
      )
      p <- p + ggplot2::geom_label(
        data = labels,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        inherit.aes = FALSE, size = 3.5, linewidth = 0.2,
        fill = "white", alpha = 0.88
      )
    }
  }

  if (isTRUE(show_means)) {
    mean_data <- data.frame(
      phenotype_1 = model$mean_matrix[1L, ],
      phenotype_2 = model$mean_matrix[2L, ],
      genotype = factor(paste("Genotype", 0:2), levels = paste("Genotype", 0:2))
    )
    if (identical(surface, "genotype_density")) {
      p <- p +
        ggplot2::geom_point(
          data = mean_data,
          ggplot2::aes(
            x = .data$phenotype_1, y = .data$phenotype_2,
            shape = .data$genotype, colour = .data$genotype
          ),
          inherit.aes = FALSE, size = 3, stroke = 1, fill = "white"
        ) +
        ggplot2::scale_shape_manual(values = c(21, 22, 24), drop = FALSE) +
        ggplot2::labs(shape = "Genotype")
    } else {
      p <- p +
        ggplot2::geom_point(
          data = mean_data,
          ggplot2::aes(x = .data$phenotype_1, y = .data$phenotype_2,
                       shape = .data$genotype),
          inherit.aes = FALSE, size = 2.7, stroke = 0.8,
          colour = "black", fill = "white"
        ) +
        ggplot2::scale_shape_manual(values = c(21, 22, 24), drop = FALSE) +
        ggplot2::labs(shape = "Genotype mean")
    }
  }
  attr(p, "falconer_model") <- model
  attr(p, "thresholds") <- threshold
  attr(p, "plot_data") <- plot_grid
  p
}
