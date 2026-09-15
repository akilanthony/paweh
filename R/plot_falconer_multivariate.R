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

#' Validate one optional two-dimensional plotting range.
#'
#' @noRd
.falconer_mv_plot_limit <- function(x, name) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) || x[1L] >= x[2L]) {
    stop(name, " must be a finite numeric vector of length 2 with lower < upper.",
         call. = FALSE)
  }
  unname(x)
}

#' Construct analytic probability-content ellipses for bivariate components.
#'
#' @noRd
.falconer_mv_probability_ellipses <- function(model, component_probs,
                                               points_per_ellipse = 181L) {
  decomposition <- eigen(model$residual_covariance_matrix, symmetric = TRUE)
  covariance_root <- decomposition$vectors %*%
    diag(sqrt(decomposition$values), nrow = 2L)
  theta <- seq(0, 2 * pi, length.out = points_per_ellipse)
  unit_circle <- rbind(cos(theta), sin(theta))
  genotype_labels <- paste("Genotype", 0:2)
  band_lower <- c(0, component_probs[-length(component_probs)])
  band_labels <- sprintf(
    "%g-%g%%", 100 * band_lower, 100 * component_probs
  )

  out <- do.call(rbind, lapply(seq_along(genotype_labels), function(j) {
    do.call(rbind, lapply(seq_along(component_probs), function(k) {
      probability <- component_probs[k]
      coordinates <- model$mean_matrix[, j] +
        sqrt(stats::qchisq(probability, df = 2L)) *
        (covariance_root %*% unit_circle)
      data.frame(
        phenotype_1 = coordinates[1L, ],
        phenotype_2 = coordinates[2L, ],
        genotype = genotype_labels[j],
        component_probability = probability,
        band_lower_probability = band_lower[k],
        component_band = band_labels[k],
        ellipse_group = paste(j, format(probability, digits = 15), sep = ":")
      )
    }))
  }))
  rownames(out) <- NULL
  out$genotype <- factor(out$genotype, levels = genotype_labels)
  out$component_band <- factor(out$component_band, levels = band_labels)
  out
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
#' @param xlim,ylim Optional finite increasing length-two vectors controlling
#'   the displayed phenotype ranges. They affect only the plotting grid and
#'   visible threshold geometry, not model parameters, thresholds, mixture
#'   weights, or other statistical calculations.
#' @param component_probs Numeric enclosed-probability levels for the analytic
#'   bivariate-normal probability regions drawn when
#'   `surface = "genotype_density"`.
#'   Values must be unique and strictly between zero and one. The default draws
#'   nested 0--50%, 50--80%, and 80--95% conditional probability bands.
#'
#' @details
#' In density mode, each grid value is the marginal, genotype-weighted mixture
#' `sum(pi[j] * f[j](y1, y2))`, where each `f[j]` is a bivariate-normal
#' density. In CDF mode, it is the lower-tail mixture probability
#' `sum(pi[j] * P(Y1 <= y1, Y2 <= y2 | G = j))`. Density and CDF surfaces
#' therefore represent different mathematical quantities.
#'
#' Genotype-density mode displays each conditional distribution
#' `f[j](y1, y2)` separately and does not multiply by genotype frequency. Each
#' conditional distribution integrates to one, including for rare genotypes.
#' Its default display uses filled nested 0--50%, 50--80%, and 80--95%
#' probability regions, with the center visually strongest and the outer band
#' lightest. Region boundaries are analytic probability-content ellipses
#' satisfying the corresponding bivariate-normal Mahalanobis-distance equation;
#' they are not arbitrary raw-density contour breaks. Three genotype components
#' do not imply that their weighted mixture has three distinct modes.
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
    return_data = FALSE,
    xlim = NULL,
    ylim = NULL,
    component_probs = c(0.50, 0.80, 0.95)
) {
  surface <- match.arg(surface)
  .falconer_check_flag(show_thresholds, "show_thresholds")
  .falconer_check_flag(show_means, "show_means")
  .falconer_check_flag(show_labels, "show_labels")
  .falconer_check_flag(return_data, "return_data")
  xlim <- .falconer_mv_plot_limit(xlim, "xlim")
  ylim <- .falconer_mv_plot_limit(ylim, "ylim")
  if (!is.numeric(component_probs) || length(component_probs) < 1L ||
      any(!is.finite(component_probs)) ||
      any(component_probs <= 0 | component_probs >= 1) ||
      anyDuplicated(component_probs)) {
    stop("component_probs must contain unique finite numeric values strictly between 0 and 1.",
         call. = FALSE)
  }
  component_probs <- sort(unname(component_probs))
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
  if (!is.null(xlim)) {
    axis_min[1L] <- xlim[1L]
    axis_max[1L] <- xlim[2L]
  }
  if (!is.null(ylim)) {
    axis_min[2L] <- ylim[1L]
    axis_max[2L] <- ylim[2L]
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
    ellipse_data <- .falconer_mv_probability_ellipses(model, component_probs)
    attr(plot_grid, "ellipse_data") <- ellipse_data
    attr(plot_grid, "component_probs") <- component_probs
  }
  if (isTRUE(return_data)) return(plot_grid)

  if (identical(surface, "genotype_density")) {
    palette <- unname(.paweh_qtl_genotype_colors())
    p <- ggplot2::ggplot(
      ellipse_data,
      ggplot2::aes(
        x = .data$phenotype_1, y = .data$phenotype_2,
        group = .data$ellipse_group,
        colour = .data$genotype, linetype = .data$genotype
      )
    )
    fill_alpha <- seq(0.46, 0.18, length.out = length(component_probs))
    for (k in rev(seq_along(component_probs))) {
      region_data <- ellipse_data[
        ellipse_data$component_probability == component_probs[k],
        , drop = FALSE
      ]
      p <- p + ggplot2::geom_polygon(
        data = region_data,
        ggplot2::aes(fill = .data$genotype),
        alpha = fill_alpha[k], colour = NA
      )
    }
    p <- p +
      ggplot2::geom_path(linewidth = 0.52, alpha = 0.82) +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
      ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
      ggplot2::scale_linetype_manual(
        values = c("solid", "dashed", "dotdash"), drop = FALSE
      ) +
      ggplot2::labs(
        fill = "Genotype", colour = "Genotype", linetype = "Genotype",
        title = if (is.null(title)) {
          "Bivariate Falconer genotype-conditional distributions"
        } else title
      )
  } else {
    p <- ggplot2::ggplot(
      plot_grid,
      ggplot2::aes(x = .data$phenotype_1, y = .data$phenotype_2, z = .data$value)
    ) +
      ggplot2::geom_contour_filled(
        ggplot2::aes(fill = ggplot2::after_stat(.data$level_mid)), bins = 14L
      ) +
      ggplot2::scale_fill_gradientn(
        colours = c("#E8ECEF", "#8FA1AF", "#355C7D", "#3F4850"),
        guide = ggplot2::guide_colorbar()
      ) +
      ggplot2::labs(
        fill = if (identical(surface, "density")) "Mixture density" else "Mixture CDF",
        title = if (is.null(title)) {
          paste("Bivariate Falconer mixture", surface)
        } else title
      )
  }
  p <- p +
    ggplot2::coord_equal(
      xlim = c(axis_min[1L], axis_max[1L]),
      ylim = c(axis_min[2L], axis_max[2L]),
      expand = FALSE
    ) +
    ggplot2::labs(x = "Phenotype 1 value", y = "Phenotype 2 value") +
    .paweh_plot_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = if (identical(surface, "genotype_density")) "top" else "right"
    )

  if (isTRUE(show_thresholds) && !is.null(threshold)) {
    tu <- threshold$upper_threshold
    tl <- threshold$lower_threshold
    regions <- data.frame(
      region = c("Affected", "Unaffected"),
      xmin = c(max(tu[1L], axis_min[1L]), axis_min[1L]),
      xmax = c(axis_max[1L], min(tl[1L], axis_max[1L])),
      ymin = c(max(tu[2L], axis_min[2L]), axis_min[2L]),
      ymax = c(axis_max[2L], min(tl[2L], axis_max[2L]))
    )
    regions <- regions[
      regions$xmin < regions$xmax & regions$ymin < regions$ymax,
      , drop = FALSE
    ]
    region_colours <- c(Affected = "#A8844F", Unaffected = "#355C7D")
    if (nrow(regions) > 0L) {
      p <- p +
        ggplot2::geom_rect(
          data = regions,
          ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                       ymin = .data$ymin, ymax = .data$ymax),
          inherit.aes = FALSE, fill = unname(region_colours[regions$region]),
          alpha = 0.10, colour = NA
        )
    }
    segments <- list()
    if (tu[2L] >= axis_min[2L] && tu[2L] <= axis_max[2L] &&
        max(tu[1L], axis_min[1L]) < axis_max[1L]) {
      segments[[length(segments) + 1L]] <- data.frame(
        x = max(tu[1L], axis_min[1L]), xend = axis_max[1L],
        y = tu[2L], yend = tu[2L]
      )
    }
    if (tu[1L] >= axis_min[1L] && tu[1L] <= axis_max[1L] &&
        max(tu[2L], axis_min[2L]) < axis_max[2L]) {
      segments[[length(segments) + 1L]] <- data.frame(
        x = tu[1L], xend = tu[1L],
        y = max(tu[2L], axis_min[2L]), yend = axis_max[2L]
      )
    }
    if (tl[2L] >= axis_min[2L] && tl[2L] <= axis_max[2L] &&
        axis_min[1L] < min(tl[1L], axis_max[1L])) {
      segments[[length(segments) + 1L]] <- data.frame(
        x = axis_min[1L], xend = min(tl[1L], axis_max[1L]),
        y = tl[2L], yend = tl[2L]
      )
    }
    if (tl[1L] >= axis_min[1L] && tl[1L] <= axis_max[1L] &&
        axis_min[2L] < min(tl[2L], axis_max[2L])) {
      segments[[length(segments) + 1L]] <- data.frame(
        x = tl[1L], xend = tl[1L],
        y = axis_min[2L], yend = min(tl[2L], axis_max[2L])
      )
    }
    threshold_segments <- if (length(segments) > 0L) {
      do.call(rbind, segments)
    } else {
      data.frame(x = numeric(), xend = numeric(), y = numeric(), yend = numeric())
    }
    if (nrow(threshold_segments) > 0L) {
      p <- p + ggplot2::geom_segment(
        data = threshold_segments,
        ggplot2::aes(
          x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend
        ),
        inherit.aes = FALSE, linewidth = 0.65
      )
    }
    attr(p, "selection_regions") <- regions
    attr(p, "threshold_segments") <- threshold_segments
    if (isTRUE(show_labels) && nrow(regions) > 0L) {
      labels <- data.frame(
        label = regions$region,
        x = (regions$xmin + regions$xmax) / 2,
        y = (regions$ymin + regions$ymax) / 2
      )
      p <- p + ggplot2::geom_label(
        data = labels,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        inherit.aes = FALSE, size = 3.5, linewidth = 0.2,
        fill = "white", alpha = 0.88
      )
      attr(p, "selection_labels") <- labels
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
  if (identical(surface, "genotype_density")) {
    attr(p, "ellipse_data") <- ellipse_data
    attr(p, "component_probs") <- component_probs
  }
  p
}
