# Interactive 3D visualization for bivariate Falconer distributions.

#' Plot an Interactive 3D Two-Phenotype Falconer Surface
#'
#' Creates an interactive Plotly visualization of three genotype-specific
#' bivariate-normal distributions. It can display either one population-mixture
#' surface or three separate genotype-conditional surfaces. The plot displays
#' exactly two quantitative phenotypes and complements
#' [plot_qtl_multivariate_contour()], which provides a static 2D view.
#'
#' @param qtl_var Numeric vector of two phenotype-specific QTL variances.
#' @param tau Numeric vector of two phenotype-specific dominance/additivity
#'   ratios.
#' @param pd Shared increaser-allele frequency in `(0, 1)`.
#' @param cor_matrix A positive-definite `2 x 2` phenotype correlation matrix.
#' @param x_upper,x_lower Optional length-two vectors of upper- and lower-tail
#'   selection percentages. Supply both or neither.
#' @param surface One of `"density"`, `"genotype_density"`, or `"cdf"`.
#' @param show_means Logical; show the three genotype-specific mean vectors at
#'   their corresponding surface heights.
#' @param show_thresholds Logical; draw the two joint threshold rectangles on
#'   the `z = 0` base plane when thresholds are supplied.
#' @param show_labels Logical; label the affected and unaffected rectangles.
#' @param grid_n Integer number of grid points along each phenotype axis.
#' @param z_scale Either `"raw"` for the calculated surface values or
#'   `"normalized"` to linearly rescale the displayed heights to `[0, 1]`.
#' @param title Optional plot title.
#'
#' @details
#' Density mode plots the marginal genotype-weighted mixture
#' `sum(pi[j] * f[j](y1, y2))`, where each `f[j]` is a bivariate-normal
#' density. CDF mode plots the lower-tail mixture
#' `sum(pi[j] * P(Y1 <= y1, Y2 <= y2 | G = j))`. Thus, the density surface
#' can have genotype-related peaks, but a three-component mixture need not have
#' three distinct modes. The CDF surface is cumulative and monotone.
#'
#' Genotype-density mode adds three Plotly surface traces, one for each
#' conditional density `f[j](y1, y2)`, without genotype-frequency weighting or
#' summation. This is the multivariate analogue of
#' [plot_qtl_genotype_distribution()]. Normalizing the z axis changes only its
#' displayed scale; raw calculated values remain available in the plot's
#' long-form `plot_data` attribute.
#'
#' Affected subjects are defined only by the joint upper-right region
#' `Y1 >= TU1 AND Y2 >= TU2`. Unaffected subjects are defined only by the
#' joint lower-left region `Y1 <= TL1 AND Y2 <= TL2`. Other phenotype
#' combinations are not selected. The statistical backend can handle more
#' traits where documented, but a surface over phenotype space is restricted
#' here to two.
#'
#' This function requires the optional package \pkg{plotly}.
#'
#' @return An object inheriting from \code{"plotly"} and \code{"htmlwidget"}.
#'   Validated model quantities, raw and displayed
#'   surface data, thresholds, and mean-marker data are attached as attributes.
#'
#' @seealso [plot_qtl_multivariate_contour()],
#' [qtl_multivariate_power_full()], and [qtl_multivariate_mssn_full()].
#'
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
#'   p <- plot_qtl_multivariate_surface3d(
#'     qtl_var = c(0.1, 0.08), tau = c(0, 0.5), pd = 0.3,
#'     cor_matrix = cor_matrix,
#'     x_upper = c(10, 10), x_lower = c(15, 15),
#'     surface = "density", grid_n = 20
#'   )
#' }
#'
#' @export
plot_qtl_multivariate_surface3d <- function(
    qtl_var,
    tau,
    pd,
    cor_matrix,
    x_upper = NULL,
    x_lower = NULL,
    surface = c("density", "genotype_density", "cdf"),
    show_means = TRUE,
    show_thresholds = TRUE,
    show_labels = TRUE,
    grid_n = 80L,
    z_scale = c("raw", "normalized"),
    title = NULL
) {
  surface <- match.arg(surface)
  z_scale <- match.arg(z_scale)
  .falconer_check_flag(show_means, "show_means")
  .falconer_check_flag(show_thresholds, "show_thresholds")
  .falconer_check_flag(show_labels, "show_labels")
  if (length(qtl_var) != 2L || length(tau) != 2L) {
    stop("plot_qtl_multivariate_surface3d() requires exactly two phenotypes.",
         call. = FALSE)
  }
  if (!is.matrix(cor_matrix) || !identical(dim(cor_matrix), c(2L, 2L))) {
    stop("cor_matrix must be a 2 x 2 matrix for this visualization.", call. = FALSE)
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for this plot.", call. = FALSE)
  }

  grid <- plot_qtl_multivariate_contour(
    qtl_var = qtl_var,
    tau = tau,
    pd = pd,
    cor_matrix = cor_matrix,
    x_upper = x_upper,
    x_lower = x_lower,
    surface = surface,
    show_thresholds = show_thresholds,
    show_means = FALSE,
    show_labels = FALSE,
    grid_n = grid_n,
    return_data = TRUE
  )
  model <- attr(grid, "falconer_model")
  threshold <- attr(grid, "thresholds")
  x_values <- sort(unique(grid$phenotype_1))
  y_values <- sort(unique(grid$phenotype_2))

  raw_range <- range(grid$value)
  raw_span <- diff(raw_range)
  scale_z <- function(z) {
    if (identical(z_scale, "raw")) return(z)
    if (raw_span <= .Machine$double.eps) return(rep(0, length(z)))
    pmin(1, pmax(0, (z - raw_range[1L]) / raw_span))
  }
  grid$z_value <- scale_z(grid$value)

  z_title <- if (identical(z_scale, "normalized")) {
    "Normalized surface"
  } else if (identical(surface, "genotype_density")) {
    "Conditional density"
  } else if (identical(surface, "density")) {
    "Mixture density"
  } else {
    "Mixture CDF"
  }

  p <- plotly::plot_ly()
  if (identical(surface, "genotype_density")) {
    palette <- unname(.pawh_qtl_genotype_colors())
    genotype_levels <- levels(grid$genotype)
    for (j in seq_along(genotype_levels)) {
      component_grid <- grid[grid$genotype == genotype_levels[j], , drop = FALSE]
      z_matrix <- matrix(
        component_grid$z_value,
        nrow = length(x_values), ncol = length(y_values)
      )
      hover_template <- paste(
        genotype_levels[j],
        "Phenotype 1: %{x:.4f}",
        "Phenotype 2: %{y:.4f}",
        paste0(z_title, ": %{z:.6f}"),
        "<extra></extra>", sep = "<br>"
      )
      p <- plotly::add_surface(
        p,
        x = x_values, y = y_values, z = t(z_matrix),
        colorscale = list(c(0, palette[j]), c(1, palette[j])),
        opacity = 0.62,
        hovertemplate = hover_template,
        contours = list(z = list(show = TRUE, usecolormap = TRUE,
                                 project = list(z = TRUE))),
        name = genotype_levels[j],
        showscale = FALSE,
        showlegend = TRUE
      )
    }
  } else {
    z_matrix <- matrix(
      grid$z_value, nrow = length(x_values), ncol = length(y_values)
    )
    hover_template <- paste(
      "Phenotype 1: %{x:.4f}",
      "Phenotype 2: %{y:.4f}",
      if (identical(z_scale, "raw")) {
        paste0(if (identical(surface, "density")) "Mixture density" else "Mixture CDF",
               ": %{z:.6f}")
      } else "Normalized surface: %{z:.6f}",
      "<extra></extra>", sep = "<br>"
    )
    p <- plotly::add_surface(
      p,
      x = x_values,
      y = y_values,
      z = t(z_matrix),
      colorscale = list(c(0, "#E8ECEF"), c(0.5, "#8FA1AF"), c(1, "#3F4850")),
      reversescale = FALSE,
      opacity = 0.90,
      hovertemplate = hover_template,
      colorbar = list(title = list(text = z_title)),
      contours = list(z = list(show = TRUE, usecolormap = TRUE,
                               project = list(z = TRUE))),
      name = z_title,
      showscale = TRUE
    )
  }

  threshold_overlays <- NULL
  if (isTRUE(show_thresholds) && !is.null(threshold)) {
    tu <- threshold$upper_threshold
    tl <- threshold$lower_threshold
    axis_min <- c(min(x_values), min(y_values))
    axis_max <- c(max(x_values), max(y_values))
    threshold_overlays <- list(
      Affected = data.frame(
        x = c(tu[1L], axis_max[1L], axis_max[1L], tu[1L], tu[1L]),
        y = c(tu[2L], tu[2L], axis_max[2L], axis_max[2L], tu[2L]),
        z = 0
      ),
      Unaffected = data.frame(
        x = c(axis_min[1L], tl[1L], tl[1L], axis_min[1L], axis_min[1L]),
        y = c(axis_min[2L], axis_min[2L], tl[2L], tl[2L], axis_min[2L]),
        z = 0
      )
    )
    overlay_colours <- c(Affected = "#A8844F", Unaffected = "#355C7D")
    for (region in names(threshold_overlays)) {
      rectangle <- threshold_overlays[[region]]
      p <- plotly::add_trace(
        p,
        x = rectangle$x, y = rectangle$y, z = rectangle$z,
        type = "scatter3d", mode = "lines",
        line = list(color = unname(overlay_colours[region]), width = 7),
        name = paste(region, "joint region"),
        hovertemplate = paste0(region, " joint-AND boundary<extra></extra>")
      )
    }
    if (isTRUE(show_labels)) {
      p <- plotly::add_trace(
        p,
        x = c(mean(c(tu[1L], axis_max[1L])),
              mean(c(axis_min[1L], tl[1L]))),
        y = c(mean(c(tu[2L], axis_max[2L])),
              mean(c(axis_min[2L], tl[2L]))),
        z = c(0, 0),
        text = c("Affected", "Unaffected"),
        type = "scatter3d", mode = "text",
        textfont = list(color = unname(overlay_colours), size = 12),
        hoverinfo = "skip", showlegend = FALSE,
        name = "Threshold labels"
      )
    }
  }

  mean_markers <- NULL
  if (isTRUE(show_means)) {
    mean_markers <- data.frame(
      phenotype_1 = model$mean_matrix[1L, ],
      phenotype_2 = model$mean_matrix[2L, ],
      genotype = paste("Genotype", 0:2)
    )
    mean_components <- vapply(1:3, function(j) {
      if (surface %in% c("density", "genotype_density")) {
        .falconer_mv_plot_density(
          mean_markers, model$mean_matrix[, j],
          model$residual_covariance_matrix
        )
      } else {
        .falconer_mv_plot_cdf(
          mean_markers, model$mean_matrix[, j],
          model$residual_covariance_matrix
        )
      }
    }, numeric(3L))
    mean_markers$value <- if (identical(surface, "genotype_density")) {
      diag(mean_components)
    } else {
      as.numeric(mean_components %*% model$genotype_frequencies)
    }
    mean_markers$z_value <- scale_z(mean_markers$value)
    mean_markers$hover_text <- sprintf(
      "%s<br>Phenotype 1: %.4f<br>Phenotype 2: %.4f<br>%s: %.6f",
      mean_markers$genotype,
      mean_markers$phenotype_1,
      mean_markers$phenotype_2,
      z_title,
      mean_markers$z_value
    )
    p <- plotly::add_trace(
      p,
      x = mean_markers$phenotype_1,
      y = mean_markers$phenotype_2,
      z = mean_markers$z_value,
      text = mean_markers$hover_text,
      type = "scatter3d", mode = "markers",
      marker = list(
        size = 5,
        color = c("#0072B2", "#D55E00", "#009E73"),
        symbol = c("circle", "square", "diamond"),
        line = list(color = "white", width = 1.5)
      ),
      hoverinfo = "text",
      name = "Genotype means"
    )
  }

  p <- plotly::layout(
    p,
    title = list(text = if (is.null(title)) {
      if (identical(surface, "genotype_density")) {
        "Bivariate Falconer genotype-conditional density surfaces"
      } else {
        paste("Bivariate Falconer mixture", surface, "surface")
      }
    } else title),
    scene = list(
      xaxis = list(title = "Phenotype 1 value", range = range(x_values)),
      yaxis = list(title = "Phenotype 2 value", range = range(y_values)),
      zaxis = list(title = z_title, rangemode = "tozero"),
      aspectmode = "cube",
      camera = list(eye = list(x = 1.45, y = -1.45, z = 1.05))
    ),
    legend = list(orientation = "h", x = 0, y = -0.08),
    margin = list(l = 10, r = 10, b = 45, t = 55)
  )

  attr(p, "falconer_model") <- model
  attr(p, "thresholds") <- threshold
  attr(p, "plot_data") <- grid
  attr(p, "mean_markers") <- mean_markers
  attr(p, "threshold_overlays") <- threshold_overlays
  attr(p, "surface") <- surface
  attr(p, "z_scale") <- z_scale
  p
}
