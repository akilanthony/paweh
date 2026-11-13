# Single-trait Falconer genotype-distribution visualization.

#' Evaluate an expression with an optional local random seed.
#'
#' @noRd
.falconer_plot_with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("seed must be NULL or a non-negative integer.", call. = FALSE)
  }
  withr::with_seed(as.integer(seed), force(code))
}

#' Plot Falconer Genotype-Specific Quantitative-Trait Distributions
#'
#' Displays the three normal quantitative-trait distributions implied by a
#' validated single-trait Falconer model. Large QTL variance separates the
#' genotype distributions relative to their residual variation; small QTL
#' variance makes them overlap heavily.
#'
#' @param qtl_var QTL variance in `(0, 1)`.
#' @param tau Dominance-to-additivity ratio.
#' @param pd Increaser-allele frequency in `(0, 1)`.
#' @param type One of `"density"` for theoretical normal curves,
#'   `"histogram"` for a simulated population, or `"binned"` for a
#'   deterministic histogram-style display calculated from Normal CDF
#'   differences.
#' @param scale Either `"density"` or `"frequency"`. For theoretical curves,
#'   frequency scaling weights each density by its Hardy--Weinberg genotype
#'   frequency. For histograms, it selects normalized densities or counts.
#' @param n Total population size simulated in histogram mode.
#' @param seed Optional non-negative integer used locally for reproducible
#'   histogram simulation. The caller's random-number state is restored.
#' @param show_means Logical; mark the three theoretical genotype means.
#' @param verbose Logical; passed to the validated Falconer parameter backend.
#' @param title Optional plot title.
#' @param return_data Logical; return the plotting data frame instead of the
#'   ggplot object. The Falconer model is retained in its `falconer_model`
#'   attribute.
#' @param bins Positive integer number of equal-width bins used when
#'   `type = "binned"`. This does not affect the existing density or simulated
#'   histogram modes.
#'
#' @details
#' Conditional on genotype `j`, the plotted trait follows a normal
#' distribution with the corresponding Falconer mean and common residual
#' variance `1 - qtl_var`. Genotypes 0, 1, and 2 have Hardy--Weinberg
#' frequencies `(1 - pd)^2`, `2 * pd * (1 - pd)`, and `pd^2`. Genotypes may
#' also be interpreted as `bb`, `Bb`, and `BB`, respectively.
#'
#' Density and binned modes are analytic and perform no simulation. In binned
#' mode, each conditional bin probability is the difference between Normal CDF
#' values at the upper and lower bin boundaries. Histogram mode first samples
#' genotypes from their Hardy--Weinberg frequencies and then samples trait
#' values from the corresponding conditional normal distributions.
#' The x-axis is the quantitative-trait value. In analytic density mode, the
#' y-axis is genotype-specific density when `scale = "density"` and
#' population-weighted density when `scale = "frequency"`. Binned mode uses
#' conditional bin probability divided by bin width, with an additional
#' genotype-frequency multiplier under `scale = "frequency"`. In simulated
#' histogram mode, the y-axis is normalized density or bin count for those
#' respective scales.
#'
#' @return A \code{ggplot} object, or a plotting data frame when
#'   `return_data = TRUE`.
#'
#' @examples
#' # A large QTL effect visibly separates the genotype distributions.
#' plot_qtl_genotype_distribution(
#'   qtl_var = 0.5, tau = 0, pd = 0.25, scale = "frequency"
#' )
#'
#' # A very small QTL effect produces almost complete overlap.
#' plot_qtl_genotype_distribution(
#'   qtl_var = 0.0005, tau = 0, pd = 0.25, scale = "frequency"
#' )
#'
#' @seealso [qtl_falconer_parameters()], [plot_qtl_anova_power()], and
#' [plot_qtl_threshold_chisq_power()].
#'
#' @export
plot_qtl_genotype_distribution <- function(
    qtl_var,
    tau,
    pd,
    type = c("density", "histogram", "binned"),
    scale = c("density", "frequency"),
    n = 3000,
    seed = NULL,
    show_means = TRUE,
    verbose = FALSE,
    title = NULL,
    return_data = FALSE,
    bins = 30L
) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  .falconer_check_flag(show_means, "show_means")
  .falconer_check_flag(verbose, "verbose")
  .falconer_check_flag(return_data, "return_data")
  model <- qtl_falconer_parameters(qtl_var, tau, pd, verbose = verbose)

  genotype_labels <- paste("Genotype", 0:2)
  means <- unname(model$mu)
  weights <- unname(model$pi)
  residual_sd <- model$residual_sd
  mean_data <- data.frame(
    genotype = factor(genotype_labels, levels = genotype_labels),
    mean = means
  )

  if (identical(type, "density")) {
    x <- seq(min(means) - 4 * residual_sd,
             max(means) + 4 * residual_sd, length.out = 501L)
    dat <- do.call(rbind, lapply(seq_along(genotype_labels), function(j) {
      density <- stats::dnorm(x, mean = means[j], sd = residual_sd)
      data.frame(
        trait_value = x,
        value = if (identical(scale, "frequency")) weights[j] * density else density,
        density = density,
        genotype_frequency = weights[j],
        theoretical_mean = means[j],
        genotype = genotype_labels[j]
      )
    }))
    dat$genotype <- factor(dat$genotype, levels = genotype_labels)
  } else if (identical(type, "histogram")) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n != floor(n) || n < 1) {
      stop("n must be a positive integer in histogram mode.", call. = FALSE)
    }
    simulated <- .falconer_plot_with_seed(seed, {
      genotype_index <- sample.int(3L, size = as.integer(n), replace = TRUE,
                                   prob = weights)
      data.frame(
        trait_value = stats::rnorm(
          n, mean = means[genotype_index], sd = residual_sd
        ),
        genotype = factor(genotype_labels[genotype_index], levels = genotype_labels),
        theoretical_mean = means[genotype_index],
        genotype_frequency = weights[genotype_index]
      )
    })
    dat <- simulated
  } else {
    if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
        bins != floor(bins) || bins < 1L) {
      stop("bins must be a positive integer in binned mode.", call. = FALSE)
    }
    support <- c(
      min(means) - 4 * residual_sd,
      max(means) + 4 * residual_sd
    )
    breaks <- seq(support[1L], support[2L], length.out = as.integer(bins) + 1L)
    bin_width <- diff(breaks)
    dat <- do.call(rbind, lapply(seq_along(genotype_labels), function(j) {
      bin_probability <- stats::pnorm(
        breaks[-1L], mean = means[j], sd = residual_sd
      ) - stats::pnorm(
        breaks[-length(breaks)], mean = means[j], sd = residual_sd
      )
      weighted_bin_probability <- weights[j] * bin_probability
      data.frame(
        trait_value = (breaks[-1L] + breaks[-length(breaks)]) / 2,
        bin_lower = breaks[-length(breaks)],
        bin_upper = breaks[-1L],
        bin_width = bin_width,
        bin_probability = bin_probability,
        weighted_bin_probability = weighted_bin_probability,
        value = if (identical(scale, "frequency")) {
          weighted_bin_probability / bin_width
        } else {
          bin_probability / bin_width
        },
        genotype_frequency = weights[j],
        theoretical_mean = means[j],
        genotype = genotype_labels[j]
      )
    }))
    dat$genotype <- factor(dat$genotype, levels = genotype_labels)
  }

  attr(dat, "falconer_model") <- model
  attr(dat, "plot_type") <- type
  attr(dat, "plot_scale") <- scale
  if (isTRUE(return_data)) return(dat)

  palette <- unname(.paweh_qtl_genotype_colors())
  if (identical(type, "density")) {
    y_label <- if (identical(scale, "density")) {
      "Genotype-specific density"
    } else "Population-weighted density"
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = .data$trait_value, y = .data$value,
        colour = .data$genotype, linetype = .data$genotype
      )
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
      ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dotdash"), drop = FALSE) +
      ggplot2::labs(colour = "Genotype", linetype = "Genotype")
  } else if (identical(type, "histogram")) {
    y_label <- if (identical(scale, "density")) "Density" else "Frequency"
    y_mapping <- if (identical(scale, "density")) {
      ggplot2::aes(y = ggplot2::after_stat(.data$density))
    } else ggplot2::aes(y = ggplot2::after_stat(.data$count))
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = .data$trait_value, fill = .data$genotype)
    ) +
      ggplot2::geom_histogram(
        mapping = y_mapping, bins = 45L, position = "identity",
        alpha = 0.48, colour = "white", linewidth = 0.15
      ) +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
      ggplot2::labs(fill = "Genotype")
  } else {
    y_label <- if (identical(scale, "density")) {
      "Conditional histogram density"
    } else {
      "Population-weighted histogram density"
    }
    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = .data$trait_value, y = .data$value,
        width = .data$bin_width, fill = .data$genotype
      )
    ) +
      ggplot2::geom_col(
        position = "identity", alpha = 0.48,
        colour = "white", linewidth = 0.15
      ) +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
      ggplot2::labs(fill = "Genotype")
  }

  if (isTRUE(show_means)) {
    p <- p + ggplot2::geom_vline(
      data = mean_data,
      ggplot2::aes(xintercept = .data$mean, linetype = .data$genotype),
      colour = "grey20", linewidth = 0.55, show.legend = FALSE
    )
  }
  p <- p +
    ggplot2::labs(
      x = "Quantitative-trait value", y = y_label,
      title = if (is.null(title)) "Falconer genotype-specific trait distributions" else title
    ) +
    .paweh_plot_theme() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "top")
  attr(p, "falconer_model") <- model
  attr(p, "plot_data") <- dat
  p
}
