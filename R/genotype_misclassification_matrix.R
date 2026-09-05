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
