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
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
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
#' @param type Either `"density"` for theoretical normal curves or
#'   `"histogram"` for a simulated population.
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
#'
#' @details
#' Conditional on genotype `j`, the plotted trait follows a normal
#' distribution with the corresponding Falconer mean and common residual
#' variance `1 - qtl_var`. Genotypes 0, 1, and 2 have Hardy--Weinberg
#' frequencies `(1 - pd)^2`, `2 * pd * (1 - pd)`, and `pd^2`. Genotypes may
#' also be interpreted as `bb`, `Bb`, and `BB`, respectively.
#'
#' Density mode is analytic and performs no simulation. Histogram mode first
#' samples genotypes from their Hardy--Weinberg frequencies and then samples
#' trait values from the corresponding conditional normal distributions.
#'
#' @return A ggplot object, or a plotting data frame when
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
#' @export
plot_qtl_genotype_distribution <- function(
    qtl_var,
    tau,
    pd,
    type = c("density", "histogram"),
    scale = c("density", "frequency"),
    n = 3000,
    seed = NULL,
    show_means = TRUE,
    verbose = FALSE,
    title = NULL,
    return_data = FALSE
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
  } else {
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
  }

  attr(dat, "falconer_model") <- model
  attr(dat, "plot_type") <- type
  attr(dat, "plot_scale") <- scale
  if (isTRUE(return_data)) return(dat)

  palette <- c("#0072B2", "#D55E00", "#009E73")
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
  } else {
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
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
  attr(p, "falconer_model") <- model
  attr(p, "plot_data") <- dat
  p
}
