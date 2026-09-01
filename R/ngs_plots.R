# Plot wrappers for the validated public sequencing design APIs. Statistical
# calculations remain in cc_ngs_power()/cc_ngs_mssn() and
# tdt_ngs_power()/tdt_ngs_mssn().

.plot_ngs_validate_coverage <- function(coverage, design) {
  if (!is.numeric(coverage) || length(coverage) < 1L ||
      any(!is.finite(coverage)) || any(coverage != floor(coverage))) {
    stop("coverage must be a non-empty vector of finite integers.")
  }
  minimum <- if (design == "tdt") 2 else 1
  if (any(coverage < minimum)) {
    if (design == "tdt" && any(coverage == 1)) {
      stop(
        "coverage = 1 is unsupported: TDT1-NGS efficient information is not ",
        "identifiable under the implemented nuisance model."
      )
    }
    stop("coverage must contain integers greater than or equal to ", minimum, ".")
  }
  sort(unique(as.numeric(coverage)))
}

.plot_ngs_validate_seq_error <- function(seq_error) {
  if (!is.numeric(seq_error) || length(seq_error) < 1L ||
      any(!is.finite(seq_error)) || any(seq_error < 0) ||
      any(seq_error >= 0.5)) {
    stop("seq_error must be a non-empty finite numeric vector in [0, 0.5).")
  }
  sort(unique(as.numeric(seq_error)))
}

.plot_ngs_validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

.plot_ngs_cc_settings <- function(args) {
  locus_het <- if (is.null(args$locus_het)) FALSE else args$locus_het
  pi <- if (is.null(args$pi)) 1 else args$pi
  args$locus_het <- NULL
  args$pi <- NULL

  .plot_ngs_validate_flag(locus_het, "locus_het")
  if (!is.numeric(pi) || length(pi) < 1L || any(!is.finite(pi)) ||
      any(pi < 0) || any(pi > 1)) {
    stop("pi must be a non-empty finite numeric vector in [0, 1].")
  }
  pi <- sort(unique(as.numeric(pi)), decreasing = TRUE)
  if (!isTRUE(locus_het) && any(pi != 1)) {
    stop(
      "pi is used only when locus_het = TRUE; set pi = 1 or enable locus heterogeneity."
    )
  }

  list(args = args, locus_het = locus_het, pi = pi)
}

.plot_ngs_reject_tdt_heterogeneity <- function(args) {
  unsupported <- intersect(c("locus_het", "pi"), names(args))
  if (length(unsupported) > 0L) {
    stop(
      "TDT1-NGS plotting does not support locus heterogeneity; remove ",
      paste(unsupported, collapse = " and "), "."
    )
  }
  invisible(TRUE)
}

.plot_ngs_series <- function(dat, design) {
  if (design == "cc") {
    dat$series <- paste0(
      "Error ", format(dat$seq_error, trim = TRUE, scientific = FALSE),
      "; pi ", format(dat$pi, trim = TRUE, scientific = FALSE)
    )
  } else {
    dat$series <- paste0(
      "Error ", format(dat$seq_error, trim = TRUE, scientific = FALSE)
    )
  }
  dat$series <- factor(dat$series, levels = unique(dat$series))
  dat
}

.plot_ngs_line <- function(dat, y, y_label, title, target_power = NULL) {
  .plot_require_ggplot2()
  plot_dat <- dat
  plot_dat$y_value <- plot_dat[[y]]
  multiple <- length(unique(plot_dat$series)) > 1L

  if (multiple) {
    p <- ggplot2::ggplot(
      plot_dat,
      ggplot2::aes(
        x = .data$coverage,
        y = .data$y_value,
        color = .data$series,
        linetype = .data$series,
        group = .data$series
      )
    )
  } else {
    p <- ggplot2::ggplot(
      plot_dat,
      ggplot2::aes(x = .data$coverage, y = .data$y_value)
    )
  }

  p <- p +
    ggplot2::geom_line(linewidth = 1.05, lineend = "round", na.rm = TRUE) +
    ggplot2::geom_point(size = 2.2, alpha = 0.9, na.rm = TRUE)
  if (!is.null(target_power)) {
    p <- p + ggplot2::geom_hline(
      yintercept = target_power,
      linetype = "dashed",
      color = "grey40"
    )
  }

  subtitle <- if ("finite_mssn" %in% names(plot_dat) &&
                  any(!plot_dat$finite_mssn)) {
    "No-finite-MSSN designs are retained in plot data with MSSN = NA."
  } else {
    NULL
  }

  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Sequencing coverage (\u00d7)",
      y = y_label,
      color = "Design",
      linetype = "Design"
    ) +
    .plot_paweh_theme()
}

#' Plot Power Sensitivity for Sequencing Designs
#'
#' Plots prospective analytic power against equal fixed sequencing coverage
#' for either a case-control sequencing trend design or a TDT1-NGS design.
#' Exact values are obtained by repeated calls to the corresponding validated
#' public design function.
#'
#' @param design Either \code{"cc"} for \code{\link{cc_ngs_power}()} or
#'   \code{"tdt"} for \code{\link{tdt_ngs_power}()}.
#' @param coverage Non-empty numeric vector of integer fixed sequencing depths.
#'   TDT1-NGS requires coverage of at least 2.
#' @param seq_error Non-empty numeric vector of symmetric per-read sequencing
#'   error probabilities in \eqn{[0,0.5)}. Multiple values produce separate
#'   curves.
#' @param target_power Optional finite scalar in \eqn{(0,1)}. If supplied, a
#'   dashed horizontal reference line is added without changing calculations.
#' @param return_data Logical. If \code{TRUE}, return the exact plotting data
#'   rather than a ggplot object.
#' @param ... Fixed arguments passed to \code{cc_ngs_power()} or
#'   \code{tdt_ngs_power()}. For case-control designs,
#'   \code{locus_het = TRUE} permits a vector of \code{pi} values. TDT1-NGS
#'   does not accept \code{locus_het} or \code{pi}.
#'
#' @details
#' Case-control output retains coverage, sequencing error, locus-homogeneity
#' fraction, sample sizes, NCP, and power. At \code{pi = 0}, the exact null
#' point is retained with NCP zero and power equal to alpha. TDT1-NGS uses raw
#' read-count probabilities and does not use the case-control hard-call model.
#'
#' Coverage is fixed and equal for the relevant study members. These plots do
#' not model variable/BGE coverage distributions, cost optimization, or
#' sample-specific depth. They introduce no simulation and delegate all
#' statistical calculations to public PAWEH sequencing design APIs. Locus
#' heterogeneity is currently available only for case-control designs.
#'
#' @return A ggplot object with exact results in \code{plot$data}, or a data
#'   frame when \code{return_data = TRUE}.
#'
#' @examples
#' plot_ngs_power(
#'   design = "cc", coverage = c(4, 12, 20), seq_error = c(0, 0.01),
#'   N_case = 1000, alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
#'   MOI = "M", verbose = FALSE
#' )
#'
#' @seealso \code{\link{plot_ngs_mssn}}, \code{\link{cc_ngs_power}},
#'   \code{\link{tdt_ngs_power}}
#' @export
plot_ngs_power <- function(
    design = c("cc", "tdt"),
    coverage,
    seq_error,
    target_power = NULL,
    return_data = FALSE,
    ...
) {
  design <- match.arg(design)
  coverage <- .plot_ngs_validate_coverage(coverage, design)
  seq_error <- .plot_ngs_validate_seq_error(seq_error)
  .plot_ngs_validate_flag(return_data, "return_data")
  if (!is.null(target_power) &&
      (!is.numeric(target_power) || length(target_power) != 1L ||
       !is.finite(target_power) || target_power <= 0 || target_power >= 1)) {
    stop("target_power must be NULL or a single finite number in (0, 1).")
  }

  args <- list(...)
  args$verbose <- NULL
  if (design == "cc") {
    settings <- .plot_ngs_cc_settings(args)
    args <- settings$args
    grid <- expand.grid(
      coverage = coverage,
      seq_error = seq_error,
      pi = settings$pi,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    rows <- lapply(seq_len(nrow(grid)), function(i) {
      call_args <- c(args, list(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        locus_het = settings$locus_het,
        pi = grid$pi[i],
        verbose = FALSE
      ))
      out <- do.call(cc_ngs_power, call_args)
      data.frame(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        pi = grid$pi[i],
        N_case = out$N_case,
        N_ctrl = out$N_ctrl,
        lambda = out$lambda,
        power = out$power
      )
    })
  } else {
    .plot_ngs_reject_tdt_heterogeneity(args)
    grid <- expand.grid(
      coverage = coverage,
      seq_error = seq_error,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    rows <- lapply(seq_len(nrow(grid)), function(i) {
      call_args <- c(args, list(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        verbose = FALSE
      ))
      out <- do.call(tdt_ngs_power, call_args)
      data.frame(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        N = out$N,
        lambda = out$lambda,
        power = out$power
      )
    })
  }

  dat <- do.call(rbind, rows)
  dat <- .plot_ngs_series(dat, design)
  if (isTRUE(return_data)) {
    dat$series <- NULL
    return(dat)
  }

  .plot_ngs_line(
    dat = dat,
    y = "power",
    y_label = "Power",
    title = if (design == "cc") {
      "Case-control sequencing power vs coverage"
    } else {
      "TDT1-NGS power vs coverage"
    },
    target_power = target_power
  )
}

#' Plot MSSN Sensitivity for Sequencing Designs
#'
#' Plots analytic minimum sample size necessary (MSSN) against equal fixed
#' sequencing coverage for either a case-control sequencing trend design or a
#' TDT1-NGS design. Calculations delegate to the corresponding public MSSN API.
#'
#' @param design Either \code{"cc"} for \code{\link{cc_ngs_mssn}()} or
#'   \code{"tdt"} for \code{\link{tdt_ngs_mssn}()}.
#' @param coverage Non-empty numeric vector of integer fixed sequencing depths.
#'   TDT1-NGS requires coverage of at least 2.
#' @param seq_error Non-empty numeric vector of symmetric per-read sequencing
#'   error probabilities in \eqn{[0,0.5)}. Multiple values produce separate
#'   curves.
#' @param return_data Logical. If \code{TRUE}, return the exact plotting data
#'   rather than a ggplot object.
#' @param ... Fixed arguments passed to \code{cc_ngs_mssn()} or
#'   \code{tdt_ngs_mssn()}. For case-control designs,
#'   \code{locus_het = TRUE} permits a vector of \code{pi} values. TDT1-NGS
#'   does not accept
#'   heterogeneity arguments.
#'
#' @details
#' The primary case-control y-axis is required cases; returned plot data also
#' retain required controls, total MSSN, and achieved power. The primary
#' TDT1-NGS y-axis is required complete trios; returned data also retain total
#' individuals and achieved power.
#'
#' At the case-control exact-null boundary \code{pi = 0}, no finite MSSN exists
#' when requested power exceeds alpha. The row is retained with MSSN and
#' achieved power set to \code{NA}, \code{finite_mssn = FALSE}, and status
#' \code{"no finite MSSN"}. No artificial finite value is plotted. Only this
#' expected scientific condition is converted; invalid inputs and unexpected
#' calculation errors propagate unchanged.
#'
#' Coverage is fixed and equal for the relevant study members. These plots do
#' not model variable/BGE coverage distributions, cost optimization, or
#' sample-specific depth. They use no simulation. Locus heterogeneity is
#' currently available only for case-control designs.
#'
#' @return A ggplot object with exact results in \code{plot$data}, or a data
#'   frame when \code{return_data = TRUE}.
#'
#' @examples
#' plot_ngs_mssn(
#'   design = "tdt", coverage = c(4, 12, 20), seq_error = 0.005,
#'   power = 0.80, pd = 0.325, R1 = 1.2, alpha = 5e-8,
#'   verbose = FALSE
#' )
#'
#' @seealso \code{\link{plot_ngs_power}}, \code{\link{cc_ngs_mssn}},
#'   \code{\link{tdt_ngs_mssn}}
#' @export
plot_ngs_mssn <- function(
    design = c("cc", "tdt"),
    coverage,
    seq_error,
    return_data = FALSE,
    ...
) {
  design <- match.arg(design)
  coverage <- .plot_ngs_validate_coverage(coverage, design)
  seq_error <- .plot_ngs_validate_seq_error(seq_error)
  .plot_ngs_validate_flag(return_data, "return_data")

  args <- list(...)
  args$verbose <- NULL
  if (design == "cc") {
    settings <- .plot_ngs_cc_settings(args)
    args <- settings$args
    grid <- expand.grid(
      coverage = coverage,
      seq_error = seq_error,
      pi = settings$pi,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    rows <- lapply(seq_len(nrow(grid)), function(i) {
      call_args <- c(args, list(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        locus_het = settings$locus_het,
        pi = grid$pi[i],
        verbose = FALSE
      ))
      safe <- .plot_safe_mssn_call(cc_ngs_mssn, call_args)
      out <- safe$result
      data.frame(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        pi = grid$pi[i],
        MSSN_case = if (safe$finite_mssn) out$MSSN_case else NA_real_,
        MSSN_ctrl = if (safe$finite_mssn) out$MSSN_ctrl else NA_real_,
        MSSN_total = if (safe$finite_mssn) out$MSSN_total else NA_real_,
        achieved_power = if (safe$finite_mssn) out$achieved_power else NA_real_,
        finite_mssn = safe$finite_mssn,
        status = safe$status,
        stringsAsFactors = FALSE
      )
    })
  } else {
    .plot_ngs_reject_tdt_heterogeneity(args)
    grid <- expand.grid(
      coverage = coverage,
      seq_error = seq_error,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    rows <- lapply(seq_len(nrow(grid)), function(i) {
      call_args <- c(args, list(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        verbose = FALSE
      ))
      safe <- .plot_safe_mssn_call(tdt_ngs_mssn, call_args)
      out <- safe$result
      data.frame(
        coverage = grid$coverage[i],
        seq_error = grid$seq_error[i],
        MSSN_trios = if (safe$finite_mssn) out$MSSN_trios else NA_real_,
        total_individuals = if (safe$finite_mssn) out$total_individuals else NA_real_,
        achieved_power = if (safe$finite_mssn) out$achieved_power else NA_real_,
        finite_mssn = safe$finite_mssn,
        status = safe$status,
        stringsAsFactors = FALSE
      )
    })
  }

  dat <- do.call(rbind, rows)
  dat <- .plot_ngs_series(dat, design)
  if (isTRUE(return_data)) {
    dat$series <- NULL
    return(dat)
  }

  .plot_ngs_line(
    dat = dat,
    y = if (design == "cc") "MSSN_case" else "MSSN_trios",
    y_label = if (design == "cc") "Required cases" else "Required trios",
    title = if (design == "cc") {
      "Case-control sequencing MSSN vs coverage"
    } else {
      "TDT1-NGS MSSN vs coverage"
    }
  )
}
