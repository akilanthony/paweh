# Internal visual system for PAWEH plots.

.paweh_color_values <- function() {
  c(
    charcoal = "#3F4850",
    navy = "#355C7D",
    slate_blue = "#6F879A",
    muted_teal = "#4F7C78",
    muted_burgundy = "#7A4F5B",
    muted_olive = "#6D7456",
    light_slate = "#8FA1AF",
    reference = "#7A848C",
    grid = "#E7EBEE",
    axis_text = "#5F6B76"
  )
}

.paweh_fill_values <- function() {
  unname(.paweh_color_values()[c(
    "navy", "muted_teal", "muted_burgundy",
    "slate_blue", "muted_olive", "charcoal"
  )])
}

.paweh_linetype_values <- function(n) {
  rep(c("solid", "longdash", "dotdash", "dotted", "twodash", "dashed"),
      length.out = n)
}

.paweh_plot_colors <- function() {
  palette <- .paweh_color_values()
  c(
    cases = palette[["charcoal"]], controls = palette[["slate_blue"]],
    genotype = palette[["charcoal"]], trend = palette[["navy"]],
    tdt_baseline = palette[["charcoal"]],
    tdt_misclassification = palette[["navy"]],
    tdt_heterogeneity = palette[["slate_blue"]],
    transmitted = palette[["charcoal"]],
    nontransmitted = palette[["light_slate"]],
    baseline = "#C7CDD2", adjusted = palette[["charcoal"]],
    reference = palette[["reference"]]
  )
}

.paweh_qtl_genotype_colors <- function() {
  palette <- .paweh_color_values()
  c(
    `Genotype 0` = palette[["charcoal"]],
    `Genotype 1` = palette[["slate_blue"]],
    `Genotype 2` = palette[["navy"]]
  )
}

.paweh_plot_theme <- function(base_size = 11) {
  palette <- .paweh_color_values()
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = base_size + 3, face = "bold", color = palette[["charcoal"]],
        lineheight = 1.05, margin = ggplot2::margin(b = 7)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.25, color = palette[["axis_text"]],
        lineheight = 1.1, margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(
        size = base_size, color = palette[["charcoal"]]
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 1, color = palette[["axis_text"]]
      ),
      legend.title = ggplot2::element_text(
        size = base_size - 0.25, face = "bold", color = palette[["charcoal"]]
      ),
      legend.text = ggplot2::element_text(
        size = base_size - 1, color = palette[["axis_text"]]
      ),
      panel.grid.major = ggplot2::element_line(
        color = palette[["grid"]], linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(10, 12, 10, 10)
    )
}
