#' Validate a True-to-Observed Genotype Matrix
#'
#' @param M Numeric three-by-three matrix with true genotypes in rows and
#'   observed genotypes in columns, both ordered 0, 1, 2.
#' @param tolerance Tolerance for row sums, not for detecting zero entries.
#' @details An exactly all-zero row is replaced by its corresponding identity
#'   row, assuming correct classification for that genotype. One warning is
#'   issued per repaired matrix. Small nonzero entries are never treated as
#'   zero, and other malformed rows are rejected without normalization.
#' @return The validated matrix, retaining dimensions and dimnames.
#' @keywords internal
.validate_genotype_misclassification_matrix <- function(M, tolerance = 1e-10) {
  if (!is.matrix(M) || !is.numeric(M) || !identical(dim(M), c(3L, 3L))) {
    stop("Genotype misclassification matrix must be a numeric 3 x 3 matrix.")
  }
  if (any(!is.finite(M)) || any(M < 0) || any(M > 1)) {
    stop("Genotype misclassification probabilities must be finite and in [0, 1].")
  }
  zero_rows <- which(rowSums(M == 0) == ncol(M))
  if (length(zero_rows)) {
    M[zero_rows, ] <- 0
    M[cbind(zero_rows, zero_rows)] <- 1
  }
  if (any(abs(rowSums(M) - 1) > tolerance)) {
    stop("Genotype misclassification matrix rows must sum to 1.")
  }
  if (length(zero_rows)) {
    warning("All-zero genotype misclassification rows detected; replaced with corresponding identity rows.",
            call. = FALSE)
  }
  M
}

.cc_misclass_matrix_1p <- function(e) {
  if (!is.numeric(e) || length(e) != 1 || e < 0 || e > 0.5)
    stop("e must be a single number in [0, 0.5].")

  M <- matrix(c(
    1 - 2 * e, e, e,
    e, 1 - 2 * e, e,
    e, e, 1 - 2 * e
  ), nrow = 3, byrow = TRUE)

  .validate_genotype_misclassification_matrix(M)
}

.cc_misclass_matrix_2p <- function(e1, e2) {
  if (!is.numeric(e1) || length(e1) != 1 || e1 < 0 || e1 > 1)
    stop("e1 must be a single number in [0,1].")
  if (!is.numeric(e2) || length(e2) != 1 || e2 < 0 || e2 > 0.5)
    stop("e2 must be a single number in [0,0.5].")

  M <- matrix(c(
    1 - e1, e1, 0,
    e2, 1 - 2 * e2, e2,
    0, e1, 1 - e1
  ), nrow = 3, byrow = TRUE)

  .validate_genotype_misclassification_matrix(M)
}

.cc_misclass_matrix_3p <- function(e01, e02, e03) {
  if (!is.numeric(e01) || length(e01) != 1 || e01 < 0 || e01 > 1)
    stop("e01 must be a single number in [0,1].")
  if (!is.numeric(e02) || length(e02) != 1 || e02 < 0 || e02 > 0.5)
    stop("e02 must be a single number in [0,0.5].")
  if (!is.numeric(e03) || length(e03) != 1 || e03 < 0 || e03 > 1)
    stop("e03 must be a single number in [0,1].")
  if (e01 + e03 > 1)
    stop("Need e01 + e03 <= 1.")

  M <- matrix(c(
    1 - (e01 + e03), e01, e03,
    e02, 1 - 2 * e02, e02,
    e03, e01, 1 - (e01 + e03)
  ), nrow = 3, byrow = TRUE)

  .validate_genotype_misclassification_matrix(M)
}

.cc_apply_genotype_misclass <- function(g_true, M_true_to_obs) {
  if (length(g_true) != 3)
    stop("g_true must be length 3.")
  as.numeric(
    t(.validate_genotype_misclassification_matrix(M_true_to_obs)) %*% g_true
  )
}

.cc_scale_3p_errors <- function(e01, e02, e03, multiplier) {
  if (!is.numeric(multiplier) || length(multiplier) != 1 || multiplier < 0)
    stop("diff_multiplier must be a single nonnegative number.")

  out <- c(
    e01 = e01 * multiplier,
    e02 = e02 * multiplier,
    e03 = e03 * multiplier
  )

  if (out["e01"] < 0 || out["e01"] > 1)
    stop("Scaled e01 is outside [0,1].")
  if (out["e02"] < 0 || out["e02"] > 0.5)
    stop("Scaled e02 is outside [0,0.5].")
  if (out["e03"] < 0 || out["e03"] > 1)
    stop("Scaled e03 is outside [0,1].")
  if ((out["e01"] + out["e03"]) > 1)
    stop("Scaled e01 + e03 > 1.")

  out
}

.cc_resolve_genotype_misclassification <- function(
    geno_misclass = c("none", "1p", "2p", "3p", "diff3p"),
    e = 0, e1 = 0, e2 = 0,
    e01 = 0, e02 = 0, e03 = 0,
    case_e01 = 0, case_e02 = 0, case_e03 = 0,
    ctrl_e01 = 0, ctrl_e02 = 0, ctrl_e03 = 0,
    diff_source = c("explicit", "case", "ctrl"),
    diff_multiplier = 1
) {
  geno_misclass <- match.arg(geno_misclass)
  diff_source <- match.arg(diff_source)
  if (!is.numeric(diff_multiplier) || length(diff_multiplier) != 1 ||
      diff_multiplier < 0) {
    stop("diff_multiplier must be a single nonnegative number.")
  }

  M_case <- diag(3)
  M_ctrl <- diag(3)
  info <- list(enabled = FALSE, model = "none")

  if (geno_misclass == "1p") {
    M_case <- M_ctrl <- .cc_misclass_matrix_1p(e)
    info <- list(enabled = e > 0, model = "1p_symmetric", e = e, M = M_case)
  } else if (geno_misclass == "2p") {
    M_case <- M_ctrl <- .cc_misclass_matrix_2p(e1, e2)
    info <- list(
      enabled = e1 > 0 || e2 > 0, model = "2p_hom_het",
      e1 = e1, e2 = e2, M = M_case
    )
  } else if (geno_misclass == "3p") {
    M_case <- M_ctrl <- .cc_misclass_matrix_3p(e01, e02, e03)
    info <- list(
      enabled = e01 > 0 || e02 > 0 || e03 > 0,
      model = "3p_homhet_homhom", e01 = e01, e02 = e02, e03 = e03,
      M = M_case
    )
  } else if (geno_misclass == "diff3p") {
    if (diff_source == "explicit") {
      case_params <- c(e01 = case_e01, e02 = case_e02, e03 = case_e03)
      ctrl_params <- c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)
    } else if (diff_source == "case") {
      case_params <- c(e01 = case_e01, e02 = case_e02, e03 = case_e03)
      ctrl_params <- .cc_scale_3p_errors(
        case_e01, case_e02, case_e03, diff_multiplier
      )
    } else {
      ctrl_params <- c(e01 = ctrl_e01, e02 = ctrl_e02, e03 = ctrl_e03)
      case_params <- .cc_scale_3p_errors(
        ctrl_e01, ctrl_e02, ctrl_e03, diff_multiplier
      )
    }
    M_case <- .cc_misclass_matrix_3p(
      case_params["e01"], case_params["e02"], case_params["e03"]
    )
    M_ctrl <- .cc_misclass_matrix_3p(
      ctrl_params["e01"], ctrl_params["e02"], ctrl_params["e03"]
    )
    info <- list(
      enabled = any(c(case_params, ctrl_params) > 0),
      model = "diff3p_homhet_homhom",
      diff_source = diff_source,
      diff_multiplier = diff_multiplier,
      case_params = case_params,
      ctrl_params = ctrl_params,
      M_case = M_case,
      M_ctrl = M_ctrl
    )
  }

  list(M_case = M_case, M_ctrl = M_ctrl, info = info)
}
