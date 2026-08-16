# Plot wrappers for the single-trait Falconer quantitative framework.

.falconer_plot_labels <- c(
  N = "Total sample size",
  N_case = "Number of selected cases",
  alpha = "Significance level",
  power = "Target power",
  qtl_var = "QTL variance",
  tau = "Dominance-to-additivity ratio",
  pd = "Increaser allele frequency",
  x_upper = "Upper-tail selection percentage",
  x_lower = "Lower-tail selection percentage",
  k = "Control-to-case ratio"
)

.falconer_plot_label <- function(name, override = NULL) {
  if (!is.null(override)) return(override)
  value <- unname(.falconer_plot_labels[name])
  if (is.na(value)) name else value
}

.falconer_plot_validate <- function(x_var, x_values, allowed, context,
                                     return_data) {
  .plot_check_x_values(x_values)
  .plot_check_x_var(x_var, allowed, context)
  .falconer_check_flag(return_data, "return_data")
  invisible(TRUE)
}

.falconer_plot_sweep <- function(fun, x_var, x_values, extractor, args0) {
  args0$verbose <- FALSE
  vapply(x_values, function(value) {
    args <- args0
    args[[x_var]] <- value
    out <- do.call(fun, args)
    extractor(out)
  }, numeric(1))
}

.falconer_plot_result <- function(x_var, x_values, y, y_name, return_data,
                                   x_label, y_label, title) {
  dat <- data.frame(x = x_values, y = y)
  if (isTRUE(return_data)) {
    names(dat) <- c(x_var, y_name)
    return(dat)
  }
  .plot_make_line(
    dat,
    x_label = x_label,
    y_label = y_label,
    title = title
  )
}

#' Plot One-Way ANOVA Power Under the Falconer Model
#'
#' Sweeps one model or design parameter and repeatedly calls
#' \code{\link{qtl_anova_power}()} for a single continuous quantitative trait.
#'
#' @param x_var Parameter to vary. One of \code{"N"}, \code{"alpha"},
#'   \code{"qtl_var"}, \code{"tau"}, or \code{"pd"}.
#' @param x_values Numeric vector of at least two finite x-axis values.
#' @param title Optional plot-title override.
#' @param x_label Optional x-axis-label override.
#' @param y_label Optional y-axis-label override.
#' @param return_data Logical. If \code{TRUE}, returns plotting data rather
#'   than a ggplot object.
#' @param ... Fixed arguments passed to \code{qtl_anova_power()}.
#'
#' @details The x-axis is \code{x_values} for the selected \code{x_var}; the
#' y-axis is one-way ANOVA power from \code{qtl_anova_power()}. Every argument
#' in \code{...} remains fixed while the selected parameter is swept.
#'
#' @return A \code{ggplot} object, or a data frame with the swept parameter and power
#'   if \code{return_data = TRUE}.
#'
#' @examples
#' plot_qtl_anova_power(
#'   x_var = "N", x_values = c(600, 800, 1000),
#'   alpha = 0.0001, qtl_var = 0.025, tau = 0.5, pd = 0.15
#' )
#'
#' @seealso \code{\link{qtl_anova_power}},
#' \code{\link{plot_qtl_anova_mssn}}.
#'
#' @export
plot_qtl_anova_power <- function(
    x_var,
    x_values,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    return_data = FALSE,
    ...
) {
  allowed <- c("N", "alpha", "qtl_var", "tau", "pd")
  .falconer_plot_validate(
    x_var, x_values, allowed, "plot_qtl_anova_power()", return_data
  )
  y <- .falconer_plot_sweep(
    qtl_anova_power, x_var, x_values,
    function(out) out$power,
    list(...)
  )
  .falconer_plot_result(
    x_var, x_values, y, "power", return_data,
    x_label = .falconer_plot_label(x_var, x_label),
    y_label = if (is.null(y_label)) "Power" else y_label,
    title = if (is.null(title)) {
      paste("QTL ANOVA power vs", .falconer_plot_label(x_var))
    } else title
  )
}

#' Plot One-Way ANOVA Minimum Sample Size Under the Falconer Model
#'
#' Sweeps one model or design parameter and repeatedly calls
#' \code{\link{qtl_anova_mssn}()}.
#'
#' @param x_var Parameter to vary. One of \code{"power"}, \code{"alpha"},
#'   \code{"qtl_var"}, \code{"tau"}, or \code{"pd"}.
#' @inheritParams plot_qtl_anova_power
#' @param ... Fixed arguments passed to \code{qtl_anova_mssn()}.
#'
#' @details The x-axis is \code{x_values} for the selected \code{x_var}; the
#' y-axis is the minimum total sample size returned by \code{qtl_anova_mssn()}.
#' Target power and other arguments in \code{...} remain fixed unless selected
#' as \code{x_var}.
#'
#' @return A \code{ggplot} object, or a data frame with the swept parameter and
#'   \code{required_N} if \code{return_data = TRUE}.
#'
#' @examples
#' plot_qtl_anova_mssn(
#'   x_var = "qtl_var", x_values = c(0.015, 0.025, 0.05),
#'   power = 0.8, alpha = 0.0001, tau = 0.5, pd = 0.15
#' )
#'
#' @seealso \code{\link{qtl_anova_mssn}},
#' \code{\link{plot_qtl_anova_power}}.
#'
#' @export
plot_qtl_anova_mssn <- function(
    x_var,
    x_values,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    return_data = FALSE,
    ...
) {
  allowed <- c("power", "alpha", "qtl_var", "tau", "pd")
  .falconer_plot_validate(
    x_var, x_values, allowed, "plot_qtl_anova_mssn()", return_data
  )
  y <- .falconer_plot_sweep(
    qtl_anova_mssn, x_var, x_values,
    function(out) out$N,
    list(...)
  )
  .falconer_plot_result(
    x_var, x_values, y, "required_N", return_data,
    x_label = .falconer_plot_label(x_var, x_label),
    y_label = if (is.null(y_label)) "Required total sample size" else y_label,
    title = if (is.null(title)) {
      paste("QTL ANOVA sample size vs", .falconer_plot_label(x_var))
    } else title
  )
}

#' Plot Threshold-Selected Genotype Chi-Square Power
#'
#' Sweeps one model, selection, or design parameter and repeatedly calls
#' \code{\link{qtl_threshold_chisq_power}()}.
#'
#' @param x_var Parameter to vary. One of \code{"N_case"}, \code{"alpha"},
#'   \code{"qtl_var"}, \code{"tau"}, \code{"pd"}, \code{"x_upper"},
#'   \code{"x_lower"}, or \code{"k"}.
#' @inheritParams plot_qtl_anova_power
#' @param ... Fixed arguments passed to \code{qtl_threshold_chisq_power()}.
#'
#' @details The x-axis is \code{x_values} for the selected model, selection,
#' or design parameter. The y-axis is power for the two-df genotype chi-square
#' test after upper-tail case and lower-tail control selection. Arguments in
#' \code{...} remain fixed while \code{x_var} is swept.
#'
#' @return A \code{ggplot} object, or a data frame with the swept parameter and power
#'   if \code{return_data = TRUE}.
#'
#' @examples
#' plot_qtl_threshold_chisq_power(
#'   x_var = "N_case", x_values = c(75, 100, 125, 150),
#'   alpha = 0.0001, qtl_var = 0.025, tau = 0.5, pd = 0.15,
#'   x_upper = 5, x_lower = 5, k = 1
#' )
#'
#' @seealso \code{\link{qtl_threshold_chisq_power}},
#' \code{\link{plot_qtl_threshold_chisq_mssn}}.
#'
#' @export
plot_qtl_threshold_chisq_power <- function(
    x_var,
    x_values,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    return_data = FALSE,
    ...
) {
  allowed <- c(
    "N_case", "alpha", "qtl_var", "tau", "pd", "x_upper", "x_lower", "k"
  )
  .falconer_plot_validate(
    x_var, x_values, allowed,
    "plot_qtl_threshold_chisq_power()", return_data
  )
  y <- .falconer_plot_sweep(
    qtl_threshold_chisq_power, x_var, x_values,
    function(out) out$power,
    list(...)
  )
  .falconer_plot_result(
    x_var, x_values, y, "power", return_data,
    x_label = .falconer_plot_label(x_var, x_label),
    y_label = if (is.null(y_label)) "Power" else y_label,
    title = if (is.null(title)) {
      paste(
        "Threshold-selected genotype chi-square power vs",
        .falconer_plot_label(x_var)
      )
    } else title
  )
}

#' Plot Threshold-Selected Genotype Chi-Square Minimum Sample Size
#'
#' Sweeps one model, selection, or design parameter and repeatedly calls
#' \code{\link{qtl_threshold_chisq_mssn}()}.
#'
#' @param x_var Parameter to vary. One of \code{"power"}, \code{"alpha"},
#'   \code{"qtl_var"}, \code{"tau"}, \code{"pd"}, \code{"x_upper"},
#'   \code{"x_lower"}, or \code{"k"}.
#' @param sample_size Selected-sample result to plot: \code{"total"},
#'   \code{"case"}, or \code{"control"}. Population screening counts are
#'   deliberately excluded because they are not statistical MSSN values.
#' @inheritParams plot_qtl_anova_power
#' @param ... Fixed arguments passed to \code{qtl_threshold_chisq_mssn()}.
#'
#' @details The x-axis is \code{x_values} for the selected model, selection,
#' or design parameter. The y-axis is the selected case, control, or total MSSN
#' requested by \code{sample_size}; it is not a source-population screening
#' count. Other arguments remain fixed while \code{x_var} is swept.
#'
#' @return A \code{ggplot} object, or a data frame with the swept parameter and the
#'   selected MSSN result if \code{return_data = TRUE}.
#'
#' @examples
#' plot_qtl_threshold_chisq_mssn(
#'   x_var = "qtl_var", x_values = c(0.015, 0.025, 0.05),
#'   power = 0.8, alpha = 0.0001, tau = 0.5, pd = 0.15,
#'   x_upper = 5, x_lower = 5, k = 1
#' )
#'
#' @seealso \code{\link{qtl_threshold_chisq_mssn}},
#' \code{\link{plot_qtl_threshold_chisq_power}}.
#'
#' @export
plot_qtl_threshold_chisq_mssn <- function(
    x_var,
    x_values,
    sample_size = c("total", "case", "control"),
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    return_data = FALSE,
    ...
) {
  sample_size <- match.arg(sample_size)
  allowed <- c("power", "alpha", "qtl_var", "tau", "pd", "x_upper", "x_lower", "k")
  .falconer_plot_validate(
    x_var, x_values, allowed,
    "plot_qtl_threshold_chisq_mssn()", return_data
  )
  extractor <- switch(
    sample_size,
    total = function(out) out$N_total,
    case = function(out) out$N_case,
    control = function(out) out$N_control
  )
  y <- .falconer_plot_sweep(
    qtl_threshold_chisq_mssn, x_var, x_values,
    extractor,
    list(...)
  )
  size_label <- switch(
    sample_size,
    total = "Required total selected sample size",
    case = "Required selected cases",
    control = "Required selected controls"
  )
  .falconer_plot_result(
    x_var, x_values, y, paste0("MSSN_", sample_size), return_data,
    x_label = .falconer_plot_label(x_var, x_label),
    y_label = if (is.null(y_label)) size_label else y_label,
    title = if (is.null(title)) {
      paste(
        "Threshold-selected genotype chi-square sample size vs",
        .falconer_plot_label(x_var)
      )
    } else title
  )
}
