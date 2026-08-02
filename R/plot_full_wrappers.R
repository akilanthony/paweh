# -------------------------------------------------------------------
# Plot wrappers for full case-control and TDT functions
# Version 4: clear labels and polished ggplot styling
# -------------------------------------------------------------------
# These functions call the existing full backend functions repeatedly while
# sweeping one x-axis parameter at a time.
#
# Required backend functions:
#   cc_power_conditional_full()
#   cc_mssn_conditional_full()
#   tdt_power_full()
#   tdt_required_trios_full()
#
# Requires ggplot2 for plotting.
# -------------------------------------------------------------------

.plot_require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }
}

.plot_drop_helper_args <- function(args) {
  helper_names <- c(
    "theta_base", "phi_base",
    "e_base", "e1_base", "e2_base", "e01_base", "e02_base", "e03_base",
    "case_e01_base", "case_e02_base", "case_e03_base",
    "ctrl_e01_base", "ctrl_e02_base", "ctrl_e03_base",
    "misclass_rate_base"
  )
  args[setdiff(names(args), helper_names)]
}

.plot_get_arg <- function(args, name, default = NULL) {
  if (!is.null(args[[name]])) args[[name]] else default
}

.plot_check_x_values <- function(x_values) {
  if (!is.numeric(x_values) || length(x_values) < 2 || any(!is.finite(x_values))) {
    stop("x_values must be a numeric vector of length at least 2 with finite values.")
  }
  invisible(TRUE)
}

.plot_check_x_var <- function(x_var, allowed, context) {
  if (!is.character(x_var) || length(x_var) != 1 || !nzchar(x_var)) {
    stop("x_var must be a single non-empty character string.")
  }
  if (!x_var %in% allowed) {
    stop(
      "Unsupported x_var for ", context, ": '", x_var, "'. Supported values are: ",
      paste(allowed, collapse = ", "),
      "."
    )
  }
  invisible(TRUE)
}

.plot_allowed_cc_power_x <- c(
  "N_case", "N_ctrl", "N_total", "alpha", "prev", "pd", "R2", "k", "pi",
  "locus_het_rate", "theta", "phi", "pheno_error_multiplier", "e", "e1",
  "e2", "e01", "e02", "e03", "case_e01", "case_e02", "case_e03",
  "ctrl_e01", "ctrl_e02", "ctrl_e03", "geno_error_multiplier",
  "diff_multiplier"
)

.plot_allowed_cc_mssn_x <- setdiff(
  c("power", .plot_allowed_cc_power_x),
  c("N_case", "N_ctrl", "N_total")
)

.plot_label <- function(name) {
  labels <- c(
    # sample size / study design
    N = "Number of affected trios",
    N_case = "Number of cases",
    N_ctrl = "Number of controls",
    N_total = "Total sample size",
    k = "Control-to-case ratio",
    power = "Target power",
    target_power = "Target power",
    alpha = "Significance level",

    # disease model
    prev = "Disease prevalence",
    pd = "Disease allele frequency",
    R1 = "Heterozygote relative risk",
    R2 = "Homozygote relative risk",
    MOI = "Mode of inheritance",
    delta_prime = "LD scale parameter (D')",

    # locus heterogeneity
    pi = "Locus homogeneity fraction (pi)",
    locus_het_rate = "Locus heterogeneity rate",
    heter_rate = "Heterogeneity rate",

    # phenotype misclassification
    theta = "Phenotype misclassification: affected to control",
    phi = "Phenotype misclassification: unaffected to case",
    misclass_rate = "Phenotype misclassification rate",
    pheno_error_multiplier = "Phenotype error multiplier",

    # genotype misclassification
    e = "Genotype error rate",
    e1 = "Genotype error rate e1",
    e2 = "Genotype error rate e2",
    e01 = "Genotype error e01",
    e02 = "Genotype error e02",
    e03 = "Genotype error e03",
    case_e01 = "Case genotype error e01",
    case_e02 = "Case genotype error e02",
    case_e03 = "Case genotype error e03",
    ctrl_e01 = "Control genotype error e01",
    ctrl_e02 = "Control genotype error e02",
    ctrl_e03 = "Control genotype error e03",
    geno_error_multiplier = "Genotype error multiplier",
    diff_multiplier = "Differential genotype error multiplier"
  )

  out <- unname(labels[name])
  if (is.na(out)) name else out
}

.plot_test_label <- function(test) {
  labels <- c(
    genotypes = "Genotype chi-square",
    trend = "Trend test"
  )
  out <- unname(labels[test])
  if (is.na(out)) test else out
}

.plot_scenario_label <- function(scenario) {
  labels <- c(
    no_error = "No error",
    misclassification = "Misclassification",
    heterogeneity = "Heterogeneity"
  )
  out <- unname(labels[scenario])
  if (is.na(out)) scenario else out
}

.plot_sample_size_label <- function(sample_size) {
  labels <- c(
    case = "Required cases",
    control = "Required controls",
    total = "Required total sample size"
  )
  out <- unname(labels[sample_size])
  if (is.na(out)) sample_size else out
}

.plot_axis_label <- function(name, override = NULL) {
  if (!is.null(override)) override else .plot_label(name)
}

.plot_title <- function(prefix, x_var, suffix = NULL, override = NULL) {
  if (!is.null(override)) return(override)
  title <- paste(prefix, "vs", .plot_label(x_var))
  if (!is.null(suffix) && nzchar(suffix)) title <- paste0(title, " (", suffix, ")")
  title
}

.plot_pretty_data <- function(dat) {
  if ("group" %in% names(dat)) {
    dat$group <- vapply(dat$group, .plot_test_label, character(1))
  }
  if ("scenario" %in% names(dat)) {
    dat$scenario <- vapply(dat$scenario, .plot_scenario_label, character(1))
  }
  dat
}

.plot_genmixr_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(size = base_size, margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}


.plot_apply_cc_x <- function(args, x_var, x) {
  # Special user-friendly aliases and multipliers first.
  if (x_var == "locus_het_rate") {
    args$locus_het <- TRUE
    args$pi <- 1 - x
    return(args)
  }

  if (x_var == "pheno_error_multiplier") {
    theta_base <- .plot_get_arg(args, "theta_base", .plot_get_arg(args, "theta", 0))
    phi_base   <- .plot_get_arg(args, "phi_base",   .plot_get_arg(args, "phi",   0))
    args$pheno_misclass <- TRUE
    args$theta <- theta_base * x
    args$phi   <- phi_base * x
    return(args)
  }

  if (x_var == "geno_error_multiplier") {
    geno_misclass <- .plot_get_arg(args, "geno_misclass", "none")

    if (geno_misclass == "1p") {
      e_base <- .plot_get_arg(args, "e_base", .plot_get_arg(args, "e", 0))
      args$e <- e_base * x

    } else if (geno_misclass == "2p") {
      e1_base <- .plot_get_arg(args, "e1_base", .plot_get_arg(args, "e1", 0))
      e2_base <- .plot_get_arg(args, "e2_base", .plot_get_arg(args, "e2", 0))
      args$e1 <- e1_base * x
      args$e2 <- e2_base * x

    } else if (geno_misclass == "3p") {
      e01_base <- .plot_get_arg(args, "e01_base", .plot_get_arg(args, "e01", 0))
      e02_base <- .plot_get_arg(args, "e02_base", .plot_get_arg(args, "e02", 0))
      e03_base <- .plot_get_arg(args, "e03_base", .plot_get_arg(args, "e03", 0))
      args$e01 <- e01_base * x
      args$e02 <- e02_base * x
      args$e03 <- e03_base * x

    } else if (geno_misclass == "diff3p") {
      args$case_e01 <- .plot_get_arg(args, "case_e01_base", .plot_get_arg(args, "case_e01", 0)) * x
      args$case_e02 <- .plot_get_arg(args, "case_e02_base", .plot_get_arg(args, "case_e02", 0)) * x
      args$case_e03 <- .plot_get_arg(args, "case_e03_base", .plot_get_arg(args, "case_e03", 0)) * x
      args$ctrl_e01 <- .plot_get_arg(args, "ctrl_e01_base", .plot_get_arg(args, "ctrl_e01", 0)) * x
      args$ctrl_e02 <- .plot_get_arg(args, "ctrl_e02_base", .plot_get_arg(args, "ctrl_e02", 0)) * x
      args$ctrl_e03 <- .plot_get_arg(args, "ctrl_e03_base", .plot_get_arg(args, "ctrl_e03", 0)) * x

    } else {
      stop("x_var='geno_error_multiplier' requires geno_misclass to be one of '1p', '2p', '3p', or 'diff3p'.")
    }

    return(args)
  }

  if (x_var == "N_ctrl") {
    if (is.null(args$N_case)) stop("For x_var='N_ctrl', supply a fixed N_case.")
    args$k <- x / args$N_case
    return(args)
  }

  if (x_var == "N_total") {
    k <- .plot_get_arg(args, "k", 1)
    args$N_case <- x / (1 + k)
    return(args)
  }

  # Direct sweep.
  args[[x_var]] <- x
  args
}

.plot_cc_extract_power <- function(out, test) {
  if (test == "genotypes") return(out$tests$genotypes$power)
  if (test == "trend") return(out$tests$trend$power)
  stop("Unknown test: ", test)
}

.plot_cc_extract_mssn <- function(out, test, sample_size) {
  result <- switch(
    test,
    genotypes = out$tests$genotypes,
    trend = out$tests$trend,
    stop("Unknown test: ", test)
  )
  if (is.null(result)) stop("Requested test output is NULL.")

  switch(
    sample_size,
    case = result$MSSN_case,
    control = result$MSSN_ctrl,
    total = result$MSSN_total,
    stop("Unknown sample_size: ", sample_size)
  )
}

.plot_make_line <- function(dat, x_label, y_label, title, subtitle = NULL) {
  .plot_require_ggplot2()
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line(linewidth = 1.05, lineend = "round") +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    .plot_genmixr_theme()
}

.plot_make_line_grouped <- function(dat, x_label, y_label, title, group_label, subtitle = NULL) {
  .plot_require_ggplot2()
  dat <- .plot_pretty_data(dat)
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$group, group = .data$group)) +
    ggplot2::geom_line(linewidth = 1.05, lineend = "round") +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      color = group_label
    ) +
    .plot_genmixr_theme()
}

#' Plot Case-Control Power from the Full Case-Control Function
#'
#' Sweeps one x-axis parameter and repeatedly calls
#' \code{\link{cc_power_conditional_full}()} to plot case-control power.
#' The wrapper supports both \code{input_mode = "model_based"} and
#' \code{input_mode = "model_free"}.
#'
#' @param x_var Character. Parameter to vary on the x-axis.
#' @param x_values Numeric vector of x-axis values.
#' @param test One of \code{"genotypes"} or \code{"trend"}.
#' @param input_mode One of \code{"model_based"} or \code{"model_free"}.
#' @param compare_tests Logical. If TRUE, plot genotype and trend tests
#'   together and ignore \code{test}.
#' @param title Optional character title override.
#' @param x_label Optional x-axis label override.
#' @param y_label Optional y-axis label override.
#' @param return_data Logical. If TRUE, return the data frame instead of a ggplot.
#' @param ... Arguments passed to \code{cc_power_conditional_full()}.
#'
#' @details
#' Supported \code{x_var} values include fixed sample-size variables
#' \code{"N_case"}, \code{"N_ctrl"}, and \code{"N_total"}; model and design
#' variables \code{"alpha"}, \code{"prev"}, \code{"pd"}, \code{"R2"}, and
#' \code{"k"}; locus-heterogeneity variables \code{"pi"} and
#' \code{"locus_het_rate"}; phenotype-misclassification variables
#' \code{"theta"}, \code{"phi"}, and \code{"pheno_error_multiplier"}; and
#' genotype-misclassification variables \code{"e"}, \code{"e1"}, \code{"e2"},
#' \code{"e01"}, \code{"e02"}, \code{"e03"}, case/control differential
#' three-parameter error rates, \code{"geno_error_multiplier"}, and
#' \code{"diff_multiplier"}.
#'
#' When \code{x_var = "pheno_error_multiplier"}, baseline phenotype-error
#' values are read from \code{theta_base} and \code{phi_base}, falling back to
#' \code{theta} and \code{phi} if the baseline arguments are not supplied. When
#' \code{x_var = "geno_error_multiplier"}, baseline genotype-error values are
#' read from \code{e_base}; \code{e1_base} and \code{e2_base}; \code{e01_base},
#' \code{e02_base}, and \code{e03_base}; or case/control baseline parameters
#' for \code{geno_misclass = "diff3p"}.
#'
#' @return A ggplot object, or a data frame if \code{return_data = TRUE}.
#'
#' @examples
#' \dontrun{
#' g_aff <- c((1 - 0.05)^2, 2 * 0.05 * (1 - 0.05), 0.05^2)
#' g_unaff <- c((1 - 0.15)^2, 2 * 0.15 * (1 - 0.15), 0.15^2)
#'
#' cc_plot_power(
#'   x_var = "phi",
#'   x_values = seq(0, 0.10, by = 0.01),
#'   test = "genotypes",
#'   input_mode = "model_free",
#'   N_case = 250,
#'   alpha = 0.01,
#'   g1 = g_aff,
#'   g0 = g_unaff,
#'   prev = 0.05,
#'   pheno_misclass = TRUE,
#'   theta = 0,
#'   k = 1
#' )
#'
#' cc_plot_power(
#'   x_var = "phi",
#'   x_values = seq(0, 0.10, by = 0.01),
#'   input_mode = "model_free",
#'   compare_tests = TRUE,
#'   N_case = 250,
#'   alpha = 0.01,
#'   g1 = g_aff,
#'   g0 = g_unaff,
#'   prev = 0.05,
#'   pheno_misclass = TRUE,
#'   theta = 0,
#'   k = 1
#' )
#' }
#'
#' @export
cc_plot_power <- function(
    x_var,
    x_values,
    test = c("genotypes", "trend"),
    input_mode = c("model_based", "model_free"),
    compare_tests = FALSE,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    return_data = FALSE,
    ...
) {
  input_mode <- match.arg(input_mode)
  .plot_check_x_values(x_values)
  .plot_check_x_var(x_var, .plot_allowed_cc_power_x, "cc_plot_power()")

  if (!is.logical(compare_tests) || length(compare_tests) != 1) {
    stop("compare_tests must be TRUE or FALSE.")
  }

  tests <- if (isTRUE(compare_tests)) {
    c("genotypes", "trend")
  } else {
    match.arg(test)
  }

  args0 <- list(...)
  args0$input_mode <- input_mode
  args0$verbose <- FALSE

  dat_list <- lapply(tests, function(test_i) {
    y <- vapply(x_values, function(x) {
      args <- .plot_apply_cc_x(args0, x_var, x)
      args <- .plot_drop_helper_args(args)
      out <- do.call(cc_power_conditional_full, args)
      .plot_cc_extract_power(out, test_i)
    }, numeric(1))

    data.frame(
      x = x_values,
      group = test_i,
      y = y,
      stringsAsFactors = FALSE
    )
  })

  plot_dat <- do.call(rbind, dat_list)

  if (isTRUE(return_data)) {
    out_dat <- .plot_pretty_data(plot_dat)
    names(out_dat) <- c(x_var, "test", "power")
    return(out_dat)
  }

  if (isTRUE(compare_tests)) {
    return(.plot_make_line_grouped(
      plot_dat,
      x_label = .plot_axis_label(x_var, x_label),
      y_label = if (!is.null(y_label)) y_label else "Power",
      title = .plot_title("Case-control power", x_var, override = title),
      group_label = "Case-control test"
    ))
  }

  .plot_make_line(
    data.frame(x = plot_dat$x, y = plot_dat$y),
    x_label = .plot_axis_label(x_var, x_label),
    y_label = if (!is.null(y_label)) y_label else "Power",
    title = .plot_title("Case-control power", x_var, suffix = .plot_test_label(tests), override = title)
  )
}

#' Plot Case-Control Minimum Sample Size from the Full Case-Control Function
#'
#' Sweeps one x-axis parameter and repeatedly calls
#' \code{\link{cc_mssn_conditional_full}()} to plot the minimum sample size
#' necessary (MSSN) for case-control tests.
#'
#' @param x_var Character. Parameter to vary on the x-axis.
#' @param x_values Numeric vector of x-axis values.
#' @param test One of \code{"genotypes"} or \code{"trend"}.
#' @param input_mode One of \code{"model_based"} or \code{"model_free"}.
#' @param sample_size One of \code{"total"}, \code{"case"}, or \code{"control"}.
#' @param compare_tests Logical. If TRUE, plot genotype and trend tests
#'   together and ignore \code{test}.
#' @param title Optional character title override.
#' @param x_label Optional x-axis label override.
#' @param y_label Optional y-axis label override.
#' @param return_data Logical. If TRUE, return the data frame instead of a ggplot.
#' @param ... Arguments passed to \code{cc_mssn_conditional_full()}.
#'
#' @details
#' The supported \code{test} values are \code{"genotypes"} and \code{"trend"}.
#' \code{sample_size} selects whether to plot required cases, controls, or
#' total sample size. \code{compare_tests = TRUE} plots both tests together.
#'
#' Supported \code{x_var} values are the same heterogeneity and
#' misclassification variables documented for \code{\link{cc_plot_power}()},
#' plus \code{"power"} for target power. Fixed sample-size variables such as
#' \code{"N_case"} are not valid because sample size is the MSSN output.
#'
#' Multiplier behavior matches \code{\link{cc_plot_power}()}:
#' \code{pheno_error_multiplier} multiplies baseline \code{theta_base} and
#' \code{phi_base}; \code{geno_error_multiplier} multiplies the corresponding
#' baseline genotype-error parameters for the selected genotype
#' misclassification model.
#'
#' @return A ggplot object, or a data frame if \code{return_data = TRUE}.
#'
#' @examples
#' \dontrun{
#' cc_plot_mssn(
#'   x_var = "geno_error_multiplier",
#'   x_values = seq(0, 3, by = 0.5),
#'   test = "trend",
#'   input_mode = "model_based",
#'   power = 0.80,
#'   alpha = 0.05,
#'   prev = 0.10,
#'   pd = 0.25,
#'   R2 = 2,
#'   MOI = "M",
#'   geno_misclass = "3p",
#'   e01_base = 0.02,
#'   e02_base = 0.01,
#'   e03_base = 0.005,
#'   k = 1
#' )
#' }
#'
#' @export
cc_plot_mssn <- function(
    x_var,
    x_values,
    test = c("genotypes", "trend"),
    input_mode = c("model_based", "model_free"),
    sample_size = c("total", "case", "control"),
    compare_tests = FALSE,
    title = NULL,
    x_label = NULL,
    y_label = NULL,
    return_data = FALSE,
    ...
) {
  input_mode <- match.arg(input_mode)
  sample_size <- match.arg(sample_size)
  .plot_check_x_values(x_values)
  .plot_check_x_var(x_var, .plot_allowed_cc_mssn_x, "cc_plot_mssn()")

  if (!is.logical(compare_tests) || length(compare_tests) != 1) {
    stop("compare_tests must be TRUE or FALSE.")
  }

  if (x_var %in% c("N_case", "N_ctrl", "N_total")) {
    stop("Sample-size x variables are not valid for MSSN plots because MSSN is the output. Use x_var='power' to vary target power.")
  }

  tests <- if (isTRUE(compare_tests)) {
    c("genotypes", "trend")
  } else {
    match.arg(test)
  }

  args0 <- list(...)
  args0$input_mode <- input_mode
  args0$verbose <- FALSE

  dat_list <- lapply(tests, function(test_i) {
    y <- vapply(x_values, function(x) {
      args <- .plot_apply_cc_x(args0, x_var, x)
      args <- .plot_drop_helper_args(args)
      out <- do.call(cc_mssn_conditional_full, args)
      .plot_cc_extract_mssn(out, test_i, sample_size)
    }, numeric(1))

    data.frame(
      x = x_values,
      group = test_i,
      y = y,
      stringsAsFactors = FALSE
    )
  })

  plot_dat <- do.call(rbind, dat_list)

  if (isTRUE(return_data)) {
    out_dat <- .plot_pretty_data(plot_dat)
    names(out_dat) <- c(x_var, "test", paste0("MSSN_", sample_size))
    return(out_dat)
  }

  if (isTRUE(compare_tests)) {
    return(.plot_make_line_grouped(
      plot_dat,
      x_label = .plot_axis_label(x_var, x_label),
      y_label = if (!is.null(y_label)) y_label else .plot_sample_size_label(sample_size),
      title = .plot_title("Case-control sample size", x_var, override = title),
      group_label = "Case-control test"
    ))
  }

  .plot_make_line(
    data.frame(x = plot_dat$x, y = plot_dat$y),
    x_label = .plot_axis_label(x_var, x_label),
    y_label = if (!is.null(y_label)) y_label else .plot_sample_size_label(sample_size),
    title = .plot_title("Case-control sample size", x_var, suffix = .plot_test_label(tests), override = title)
  )
}
