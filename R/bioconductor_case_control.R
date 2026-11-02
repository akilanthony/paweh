.cc_se_genotype_frequencies <- function(se, variant, assay_name, group_col,
                                        case_value, control_value) {
  if (!methods::is(se, "SummarizedExperiment"))
    stop("se must be a SummarizedExperiment or RangedSummarizedExperiment object.")

  if (!is.character(assay_name) || length(assay_name) != 1L ||
      is.na(assay_name) || !nzchar(assay_name))
    stop("assay_name must be a single non-empty assay name.")
  assay_names <- SummarizedExperiment::assayNames(se)
  if (sum(assay_names == assay_name, na.rm = TRUE) != 1L)
    stop("assay_name must identify exactly one assay in se.")

  if (is.character(variant) && length(variant) == 1L &&
      !is.na(variant) && nzchar(variant)) {
    variant_index <- which(rownames(se) == variant)
    if (length(variant_index) != 1L)
      stop("variant must match exactly one row name in se.")
  } else if (is.numeric(variant) && length(variant) == 1L &&
             is.finite(variant) && variant == floor(variant) &&
             variant >= 1 && variant <= nrow(se)) {
    variant_index <- variant
  } else {
    stop("variant must be a unique row name or a valid positive integer row index.")
  }

  sample_data <- SummarizedExperiment::colData(se)
  if (!is.character(group_col) || length(group_col) != 1L ||
      is.na(group_col) || !nzchar(group_col))
    stop("group_col must be a single non-empty column name.")
  if (sum(colnames(sample_data) == group_col, na.rm = TRUE) != 1L)
    stop("group_col must identify exactly one column in colData(se).")
  groups <- sample_data[[group_col]]
  if ((!is.numeric(groups) && !is.character(groups) && !is.factor(groups)) ||
      !is.null(dim(groups)))
    stop("The grouping column must be a numeric, character, or factor vector.")

  for (label in list(case_value, control_value)) {
    if ((!is.numeric(label) && !is.character(label)) || length(label) != 1L ||
        is.na(label) || (is.numeric(label) && !is.finite(label)) ||
        (is.character(label) && !nzchar(label)))
      stop("case_value and control_value must be single non-missing numeric or character labels.")
  }
  # Compare explicit labels, including factor levels, without using factor codes.
  groups <- as.character(groups)
  case_value <- as.character(case_value)
  control_value <- as.character(control_value)
  if (case_value == control_value)
    stop("case_value and control_value must be different.")
  case_samples <- which(!is.na(groups) & groups == case_value)
  control_samples <- which(!is.na(groups) & groups == control_value)
  if (!length(case_samples) || !length(control_samples))
    stop("Both case_value and control_value must identify samples in colData(se).")

  # Subset before realizing the row so matrix-like assays stay inexpensive.
  selected_samples <- c(case_samples, control_samples)
  genotype_row <- SummarizedExperiment::assay(se, assay_name)[
    variant_index, selected_samples, drop = FALSE
  ]
  genotypes <- as.vector(as.matrix(genotype_row))
  if (!is.numeric(genotypes) ||
      any(!is.na(genotypes) & !(genotypes %in% 0:2)))
    stop("Selected genotypes must be numeric hard calls coded 0, 1, or 2, with NA for missing calls.")

  case_genotypes <- genotypes[seq_along(case_samples)]
  control_genotypes <- genotypes[length(case_samples) + seq_along(control_samples)]
  case_genotypes <- case_genotypes[!is.na(case_genotypes)]
  control_genotypes <- control_genotypes[!is.na(control_genotypes)]
  if (!length(case_genotypes) || !length(control_genotypes))
    stop("Each group must retain at least one non-missing genotype call.")

  list(
    g1 = tabulate(case_genotypes + 1L, nbins = 3L) / length(case_genotypes),
    g0 = tabulate(control_genotypes + 1L, nbins = 3L) / length(control_genotypes)
  )
}

#' Case-Control Minimum Sample Size from SummarizedExperiment Genotypes
#'
#' Estimates case and control genotype frequencies for one variant in a
#' SummarizedExperiment object, then calls \code{cc_mssn()} in
#' \code{input_mode = "model_free"}. This convenience wrapper uses the existing
#' genotype chi-square and genotype trend-test methods without changing them.
#'
#' @param se A \code{SummarizedExperiment} or
#'   \code{RangedSummarizedExperiment} object with variants in rows and samples
#'   in columns.
#' @param variant A single unique row name or positive integer row index.
#'   Duplicate row names can only be selected by row index.
#' @param assay_name A single assay name identifying the genotype hard calls.
#'   Must be supplied explicitly.
#' @param group_col A single column name in \code{colData(se)} containing
#'   numeric, character, or factor case/control status.
#' @param case_value,control_value Single numeric or character labels identifying
#'   cases and controls. Both must be supplied, must differ, and must occur in
#'   the grouping column. Factor levels are compared by label, not integer code;
#'   labels and group values are compared as character strings.
#' @param power Numeric in \eqn{(0,1)}. Desired target power.
#' @param alpha Numeric in \eqn{(0,1)}. Significance level.
#' @param k Numeric \eqn{> 0}. Planned control-to-case sample size ratio.
#'   Defaults to 1; it is independent of the source dataset's group sizes.
#' @param w Numeric vector of length 3. Genotype trend-test scores in genotype
#'   order 0, 1, 2. Defaults to \code{c(0, 1, 2)}.
#' @param verbose Logical. If \code{TRUE} (default), prints the original paweh
#'   formatted summary once. If \code{FALSE}, suppresses the summary.
#' @param ... Optional study-design modifiers accepted by \code{cc_mssn()}:
#'   locus heterogeneity, phenotype misclassification (including \code{prev}),
#'   and genotype misclassification settings. Use exact, unique argument names.
#'   Arguments already supplied by this wrapper, including \code{g1},
#'   \code{g0}, and \code{input_mode}, cannot be overridden. Model-based risk
#'   inputs \code{pd}, \code{R2}, and \code{MOI} are not accepted.
#'
#' @details
#' The selected assay must contain numeric diploid biallelic hard-call
#' genotypes: 0, 1, or 2 copies of the same specified allele for every sample.
#' Fractional dosages, character genotype calls, and expression counts are not
#' supported. Users must select a genotype assay; numeric values alone cannot
#' establish the biological meaning of an assay.
#'
#' Samples with missing group status or labels other than \code{case_value}
#' and \code{control_value} are excluded. Missing genotype calls (\code{NA}
#' or \code{NaN}) are excluded separately within each selected group, and each
#' group must retain at least one usable call. Only the selected variant and
#' groups are checked for invalid genotype values. Frequencies are ordered as
#' \code{c(P(G=0), P(G=1), P(G=2))}, including zero for an absent genotype.
#'
#' The source dataset only estimates genotype frequencies. The study-design
#' arguments specify the prospective design. Numerical and model-specific
#' validation is performed by the corresponding original paweh function.
#'
#' @return The unchanged \code{cc_mssn} object returned invisibly by
#'   \code{\link{cc_mssn}}, including genotype and trend-test results,
#'   frequencies, and any requested sensitivity scenarios.
#' @seealso \code{\link{cc_mssn}}, \code{\link{cc_power_from_se}}
#' @examples
#' genotypes <- matrix(c(0, 1, 1, 2, 0, 0, 1, 2), nrow = 1,
#'                     dimnames = list("rs_example", NULL))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(GT = genotypes),
#'   colData = S4Vectors::DataFrame(
#'     status = rep(c("case", "control"), each = 4)
#'   )
#' )
#' cc_mssn_from_se(
#'   se, variant = "rs_example", assay_name = "GT", group_col = "status",
#'   case_value = "case", control_value = "control",
#'   power = 0.8, alpha = 0.05
#' )
#' @export
cc_mssn_from_se <- function(se, variant, assay_name, group_col,
                            case_value, control_value, power, alpha,
                            k = 1, w = c(0, 1, 2), verbose = TRUE, ...) {
  design_args <- list(...)
  allowed_args <- setdiff(names(formals(cc_mssn)), c(
    "power", "alpha", "input_mode", "g1", "g0", "k", "w", "verbose",
    "pd", "R2", "MOI"
  ))
  if (length(design_args) &&
      (is.null(names(design_args)) || any(!names(design_args) %in% allowed_args) ||
       anyDuplicated(names(design_args))))
    stop("... must contain uniquely named study-design modifiers using exact cc_mssn argument names; g1, g0, and input_mode cannot be overridden.")

  frequencies <- .cc_se_genotype_frequencies(
    se, variant, assay_name, group_col, case_value, control_value
  )
  invisible(do.call(cc_mssn, c(
    list(power = power, alpha = alpha, input_mode = "model_free",
         k = k, w = w, verbose = verbose),
    frequencies, design_args
  )))
}

#' Case-Control Power from SummarizedExperiment Genotypes
#'
#' Estimates case and control genotype frequencies for one variant in a
#' SummarizedExperiment object, then calls \code{cc_power()} in
#' \code{input_mode = "model_free"}. This convenience wrapper uses the existing
#' genotype chi-square and genotype trend-test methods without changing them.
#'
#' @inheritParams cc_mssn_from_se
#' @param N_case Numeric \eqn{> 0}. Planned number of cases, independent of the
#'   source dataset's group sizes. Planned controls are determined by \code{k}.
#' @param ... Optional study-design modifiers accepted by \code{cc_power()}:
#'   locus heterogeneity, phenotype misclassification (including \code{prev}),
#'   and genotype misclassification settings. Use exact, unique argument names.
#'   Arguments already supplied by this wrapper, including \code{g1},
#'   \code{g0}, and \code{input_mode}, cannot be overridden. Model-based risk
#'   inputs \code{pd}, \code{R2}, and \code{MOI} are not accepted.
#' @inherit cc_mssn_from_se details
#' @return The unchanged \code{cc_power} object returned invisibly by
#'   \code{\link{cc_power}}, including genotype and trend-test results,
#'   frequencies, and any requested sensitivity scenarios.
#' @seealso \code{\link{cc_power}}, \code{\link{cc_mssn_from_se}}
#' @examples
#' genotypes <- matrix(c(0, 1, 1, 2, 0, 0, 1, 2), nrow = 1,
#'                     dimnames = list("rs_example", NULL))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(GT = genotypes),
#'   colData = S4Vectors::DataFrame(
#'     status = rep(c(1, 0), each = 4)
#'   )
#' )
#' cc_power_from_se(
#'   se, variant = 1, assay_name = "GT", group_col = "status",
#'   case_value = 1, control_value = 0,
#'   N_case = 500, alpha = 0.05, k = 2
#' )
#' @export
cc_power_from_se <- function(se, variant, assay_name, group_col,
                             case_value, control_value, N_case, alpha,
                             k = 1, w = c(0, 1, 2), verbose = TRUE, ...) {
  design_args <- list(...)
  allowed_args <- setdiff(names(formals(cc_power)), c(
    "N_case", "alpha", "input_mode", "g1", "g0", "k", "w", "verbose",
    "pd", "R2", "MOI"
  ))
  if (length(design_args) &&
      (is.null(names(design_args)) || any(!names(design_args) %in% allowed_args) ||
       anyDuplicated(names(design_args))))
    stop("... must contain uniquely named study-design modifiers using exact cc_power argument names; g1, g0, and input_mode cannot be overridden.")

  frequencies <- .cc_se_genotype_frequencies(
    se, variant, assay_name, group_col, case_value, control_value
  )
  invisible(do.call(cc_power, c(
    list(N_case = N_case, alpha = alpha, input_mode = "model_free",
         k = k, w = w, verbose = verbose),
    frequencies, design_args
  )))
}
