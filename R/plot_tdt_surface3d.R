#' Plot a Generalized 3D TDT Sensitivity Surface
#'
#' Builds an interactive sensitivity surface for model-based TDT power or
#' minimum sample size. Any two supported model parameters can be swept while
#' the remaining parameters are held fixed. Every grid point is evaluated by
#' the canonical [tdt_power()] or [tdt_mssn()] backend.
#'
#' @param metric Character. Surface height: `"power"` (default) or `"mssn"`.
#' @param scenario Character. One of `"misclassification"` (default),
#'   `"heterogeneity"`, or `"no_error"`.
#' @param x,y Character scalars naming two distinct swept parameters. Supported
#'   axes are `"pd"`, `"prev"`, `"R1"`, `"R2"`, `"alpha"`,
#'   `"delta_prime"`, `"misclass_rate"`, and `"heter_rate"`.
#' @param x_values,y_values Numeric grid values for `x` and `y`. Each must
#'   contain at least two distinct, finite values. If `NULL`, documented
#'   parameter-specific defaults are used; see Details.
#' @param N Numeric scalar greater than zero. Number of affected trios, used
#'   when `metric = "power"`.
#' @param target_power Numeric scalar in `(0, 1)`, used when
#'   `metric = "mssn"`.
#' @param pd Numeric in `(0, 1)`. Fixed risk-allele frequency when `pd` is not
#'   an axis.
#' @param prev Numeric in `(0, 1)`. Fixed disease prevalence when `prev` is not
#'   an axis.
#' @param R1,R2 Positive numeric scalars. Fixed genotype relative risks when
#'   they are not axes.
#' @param alpha Numeric in `(0, 1)`. Fixed significance level when `alpha` is
#'   not an axis.
#' @param delta_prime Numeric in `[0, 1]`. Fixed positive linkage-disequilibrium
#'   scale when `delta_prime` is not an axis.
#' @param misclass_rate Numeric in `[0, 1)`. Fixed phenotype
#'   misclassification rate when it is not an axis.
#' @param heter_rate Numeric in `[0, 1)`. Fixed locus-heterogeneity rate when
#'   it is not an axis.
#' @param ceiling_N Logical scalar. For MSSN surfaces, whether to plot the
#'   integer ceiling of the canonical sample-size result (default `TRUE`). The
#'   unrounded result remains available in the attached surface data.
#' @param title Optional character scalar. If `NULL`, an informative title is
#'   generated from `metric` and `scenario`.
#'
#' @details
#' The supported axes and their default grids are:
#' `pd = seq(0.10, 0.50, length.out = 20)`,
#' `prev = seq(0.01, 0.10, length.out = 20)`,
#' `R1 = seq(1.1, 2, length.out = 20)`,
#' `R2 = seq(1.2, 3, length.out = 20)`,
#' `alpha = seq(0.01, 0.10, length.out = 20)`,
#' `delta_prime = seq(0.25, 1, length.out = 20)`,
#' `misclass_rate = seq(0, 0.20, length.out = 20)`, and
#' `heter_rate = seq(0, 0.50, length.out = 20)`.
#' The package's implemented LD convention is
#' `D = delta_prime * pd * (1 - pd)`, where `delta_prime` is the proportion of
#' maximum positive disequilibrium: `0` represents linkage equilibrium and `1`
#' represents maximum positive LD under the model assumptions. Negative LD
#' would require a different, allele-frequency-dependent normalization and is
#' not represented by this parameterization.
#'
#' The surface is deliberately model-based. `N` and `target_power` define the
#' design calculation and therefore remain fixed rather than becoming axes.
#' The model-free inputs `ET`, `ENT`, and `n_trios`, as well as nonnumeric
#' controls such as `verbose`, are outside this visualization interface.
#' MSSN is measured in affected-child trios. `heter_rate` is the heterogeneous
#' fraction, `1 - pi`.
#'
#' The canonical TDT backends report separate scenarios, not a combined-error
#' scenario. Consequently, `misclass_rate` is an active axis only for
#' `scenario = "misclassification"`, and `heter_rate` is active only for
#' `scenario = "heterogeneity"`. An inactive modifier axis is rejected rather
#' than silently producing a flat, misleading surface.
#'
#' This function requires the optional `plotly` package. The plotting layer
#' introduces no independent statistical approximation: it only constructs the
#' grid, delegates calculations to [tdt_power()] or [tdt_mssn()], extracts the
#' selected scenario, and renders the result.
#'
#' @return An object inheriting from `"plotly"` and `"htmlwidget"`. The
#'   Cartesian grid is attached as `attr(x, "surface_data")`, with dynamically
#'   named axis columns, `metric_value` (the plotted height), and
#'   `raw_metric_value`. The plotted z matrix is attached as
#'   `attr(x, "surface_matrix")`, and calculation settings as
#'   `attr(x, "surface_spec")`.
#'
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   power_surface <- plot_tdt_surface3d(
#'     metric = "power", scenario = "misclassification",
#'     x = "pd", y = "misclass_rate",
#'     x_values = c(0.2, 0.3, 0.4),
#'     y_values = c(0, 0.03, 0.06),
#'     N = 600
#'   )
#'
#'   mssn_surface <- plot_tdt_surface3d(
#'     metric = "mssn", scenario = "heterogeneity",
#'     x = "prev", y = "heter_rate",
#'     x_values = c(0.02, 0.05, 0.08),
#'     y_values = c(0, 0.1, 0.2),
#'     target_power = 0.8
#'   )
#' }
#'
#' @seealso [tdt_power()], [tdt_mssn()], [plot_tdt_power()],
#'   [plot_tdt_mssn()]
#' @export
plot_tdt_surface3d <- function(
    metric = c("power", "mssn"),
    scenario = c("misclassification", "heterogeneity", "no_error"),
    x = "pd",
    y = "misclass_rate",
    x_values = NULL,
    y_values = NULL,
    N = 600,
    target_power = 0.80,
    pd = 0.30,
    prev = 0.05,
    R1 = 1.5,
    R2 = 2.25,
    alpha = 0.05,
    delta_prime = 1,
    misclass_rate = 0.01,
    heter_rate = 0.10,
    ceiling_N = TRUE,
    title = NULL
) {
  if (!.tdt_surface_plotly_available()) {
    stop("Package 'plotly' is required for this plot.", call. = FALSE)
  }

  metric <- match.arg(metric)
  scenario <- match.arg(scenario)
  axis_specs <- .tdt_surface_axis_specs()
  axis_names <- names(axis_specs)

  .tdt_surface_validate_axis_name(x, "x", axis_names)
  .tdt_surface_validate_axis_name(y, "y", axis_names)
  if (identical(x, y)) {
    stop("x and y must name two distinct parameters.", call. = FALSE)
  }

  active_modifier <- switch(
    scenario,
    misclassification = "misclass_rate",
    heterogeneity = "heter_rate",
    no_error = character(0)
  )
  modifier_axes <- intersect(c(x, y), c("misclass_rate", "heter_rate"))
  inactive_axes <- setdiff(modifier_axes, active_modifier)
  if (length(inactive_axes) > 0L) {
    stop(
      sprintf(
        "Axis '%s' is inactive for scenario = '%s'.",
        inactive_axes[[1L]], scenario
      ),
      call. = FALSE
    )
  }

  if (is.null(x_values)) x_values <- axis_specs[[x]]$default
  if (is.null(y_values)) y_values <- axis_specs[[y]]$default
  .tdt_surface_validate_grid(x_values, "x_values", x, axis_specs[[x]]$valid)
  .tdt_surface_validate_grid(y_values, "y_values", y, axis_specs[[y]]$valid)

  fixed <- list(
    pd = pd, prev = prev, R1 = R1, R2 = R2, alpha = alpha,
    delta_prime = delta_prime, misclass_rate = misclass_rate,
    heter_rate = heter_rate
  )
  for (parameter in setdiff(axis_names, c(x, y))) {
    .tdt_surface_validate_scalar(
      fixed[[parameter]], parameter, axis_specs[[parameter]]$valid
    )
  }
  if (metric == "power") {
    .tdt_surface_validate_scalar(N, "N", function(value) value > 0)
  } else {
    .tdt_surface_validate_scalar(
      target_power, "target_power", function(value) value > 0 && value < 1
    )
  }
  if (!is.logical(ceiling_N) || length(ceiling_N) != 1L || is.na(ceiling_N)) {
    stop("ceiling_N must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))) {
    stop("title must be NULL or a non-missing character scalar.", call. = FALSE)
  }

  surface_data <- expand.grid(
    x_value = x_values,
    y_value = y_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  names(surface_data)[1:2] <- c(x, y)

  raw_values <- numeric(nrow(surface_data))
  backend_warnings <- character(0)
  for (i in seq_len(nrow(surface_data))) {
    point <- fixed
    point[[x]] <- surface_data[[x]][[i]]
    point[[y]] <- surface_data[[y]][[i]]

    raw_values[[i]] <- tryCatch(
      withCallingHandlers(
        .tdt_surface_backend_value(
          metric = metric,
          scenario = scenario,
          point = point,
          N = N,
          target_power = target_power
        ),
        warning = function(warning) {
          backend_warnings <<- c(backend_warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) {
        stop(
          sprintf(
            "TDT surface evaluation failed at %s = %s, %s = %s: %s",
            x, format(surface_data[[x]][[i]]),
            y, format(surface_data[[y]][[i]]),
            conditionMessage(error)
          ),
          call. = FALSE
        )
      }
    )
  }
  if (length(backend_warnings) > 0L) {
    warning(
      paste(unique(backend_warnings), collapse = " "),
      " (Repeated backend warnings were consolidated.)",
      call. = FALSE
    )
  }

  plotted_values <- raw_values
  if (metric == "mssn" && isTRUE(ceiling_N)) {
    plotted_values <- ceiling(plotted_values)
  }
  surface_data$metric_value <- plotted_values
  surface_data$raw_metric_value <- raw_values
  surface_data$metric <- metric
  surface_data$scenario <- scenario
  z <- matrix(
    plotted_values,
    nrow = length(y_values),
    ncol = length(x_values),
    byrow = TRUE
  )

  x_label <- axis_specs[[x]]$label
  y_label <- axis_specs[[y]]$label
  z_label <- if (metric == "power") "Power" else "MSSN (affected trios)"
  scenario_label <- gsub("_", " ", scenario, fixed = TRUE)
  if (is.null(title)) {
    title <- sprintf("TDT %s surface: %s", toupper(metric), scenario_label)
  }
  hover <- paste0(
    x_label, ": %{x:.4g}<br>",
    y_label, ": %{y:.4g}<br>",
    z_label, ": %{z:.4g}<br>",
    "Scenario: ", scenario_label, "<extra></extra>"
  )

  plot <- plotly::plot_ly(
    x = x_values,
    y = y_values,
    z = z,
    type = "surface",
    hovertemplate = hover
  )
  plot <- plotly::layout(
    plot,
    title = list(text = title),
    scene = list(
      xaxis = list(title = x_label),
      yaxis = list(title = y_label),
      zaxis = list(title = z_label)
    )
  )

  attr(plot, "surface_data") <- surface_data
  attr(plot, "surface_matrix") <- z
  attr(plot, "surface_spec") <- list(
    metric = metric,
    scenario = scenario,
    x = x,
    y = y,
    x_values = x_values,
    y_values = y_values,
    fixed_parameters = fixed[setdiff(names(fixed), c(x, y))],
    N = N,
    target_power = target_power,
    ceiling_N = ceiling_N
  )
  plot
}

.tdt_surface_plotly_available <- function() {
  requireNamespace("plotly", quietly = TRUE)
}

.tdt_surface_axis_specs <- function() {
  list(
    pd = list(
      label = "Risk-allele frequency (p_d)",
      default = seq(0.10, 0.50, length.out = 20),
      valid = function(value) value > 0 & value < 1
    ),
    prev = list(
      label = "Disease prevalence",
      default = seq(0.01, 0.10, length.out = 20),
      valid = function(value) value > 0 & value < 1
    ),
    R1 = list(
      label = "Heterozygote relative risk (R1)",
      default = seq(1.1, 2, length.out = 20),
      valid = function(value) value > 0
    ),
    R2 = list(
      label = "Homozygote relative risk (R2)",
      default = seq(1.2, 3, length.out = 20),
      valid = function(value) value > 0
    ),
    alpha = list(
      label = "Significance level (alpha)",
      default = seq(0.01, 0.10, length.out = 20),
      valid = function(value) value > 0 & value < 1
    ),
    delta_prime = list(
      label = "LD scale (D')",
      default = seq(0.25, 1, length.out = 20),
      valid = function(value) value >= 0 & value <= 1
    ),
    misclass_rate = list(
      label = "Phenotype misclassification rate (pi01)",
      default = seq(0, 0.20, length.out = 20),
      valid = function(value) value >= 0 & value < 1
    ),
    heter_rate = list(
      label = "Locus heterogeneity rate (1 - pi)",
      default = seq(0, 0.50, length.out = 20),
      valid = function(value) value >= 0 & value < 1
    )
  )
}

.tdt_surface_validate_axis_name <- function(value, argument, supported) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% supported) {
    stop(
      sprintf(
        "%s must be one of: %s.",
        argument, paste(supported, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

.tdt_surface_validate_grid <- function(values, argument, parameter, valid) {
  if (!is.numeric(values) || length(values) < 2L || any(!is.finite(values))) {
    stop(
      sprintf(
        "%s must contain at least two finite numeric values.",
        argument
      ),
      call. = FALSE
    )
  }
  if (anyDuplicated(values)) {
    stop(sprintf("Grid values for '%s' must be distinct.", parameter),
         call. = FALSE)
  }
  if (!all(valid(values))) {
    stop(sprintf("Grid values for '%s' are outside its valid range.", parameter),
         call. = FALSE)
  }
}

.tdt_surface_validate_scalar <- function(value, parameter, valid) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      !isTRUE(valid(value))) {
    stop(sprintf("%s has an invalid fixed value.", parameter), call. = FALSE)
  }
}

.tdt_surface_backend_value <- function(
    metric, scenario, point, N, target_power
) {
  common <- c(
    point,
    list(input_mode = "model_based", verbose = FALSE)
  )
  if (metric == "power") {
    result <- do.call(tdt_power, c(list(N = N), common))
  } else {
    result <- do.call(tdt_mssn, c(list(target_power = target_power), common))
  }
  value <- .tdt_surface_extract_value(result, metric, scenario)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("The canonical backend returned a non-finite metric value.",
         call. = FALSE)
  }
  value
}

.tdt_surface_extract_value <- function(result, metric, scenario) {
  component <- if (identical(metric, "power")) "power" else "N"
  result[[component]][[scenario]]
}
