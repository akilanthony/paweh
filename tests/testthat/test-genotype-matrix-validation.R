test_that("exact-zero genotype rows become identity rows", {
  valid <- matrix(c(.98, .02, 0, .005, .99, .005, 0, .02, .98),
                  3, byrow = TRUE)
  dimnames(valid) <- list(paste0("true_", 0:2), paste0("observed_", 0:2))
  for (rows in list(3L, 1L, 2L, c(1L, 3L), 1:3)) {
    input <- expected <- valid
    input[rows, ] <- 0
    expected[rows, ] <- diag(3)[rows, , drop = FALSE]
    expect_warning(actual <- .validate_genotype_misclassification_matrix(input),
                   "All-zero genotype misclassification rows")
    expect_identical(actual, expected)
  }
  expect_silent(expect_identical(.validate_genotype_misclassification_matrix(diag(3)), diag(3)))
  expect_silent(expect_identical(.validate_genotype_misclassification_matrix(valid), valid))
  tiny <- diag(3)
  tiny[1, ] <- c(1 - 1e-14, 1e-14, 0)
  expect_identical(.validate_genotype_misclassification_matrix(tiny), tiny)
})

test_that("nonzero malformed genotype rows are never repaired", {
  for (row in list(c(.2, .3, .1), rep(1e-15, 3))) {
    bad <- diag(3)
    bad[1, ] <- row
    expect_error(.validate_genotype_misclassification_matrix(bad), "rows must sum to 1")
  }
  for (value in c(-.1, -1e-15, NA_real_, NaN, Inf, -Inf, 1.1)) {
    bad <- diag(3)
    bad[1, 2] <- value
    expect_error(.validate_genotype_misclassification_matrix(bad), "finite and in")
  }
  for (bad in list(matrix(0, 2, 3), rep(0, 9), matrix("0", 3, 3))) {
    expect_error(.validate_genotype_misclassification_matrix(bad), "numeric 3 x 3")
  }
})

test_that("CC NGS power accepts zero rows at the matrix boundary", {
  E <- matrix(c(.98, .02, 0, .005, .99, .005, 0, 0, 0), 3, byrow = TRUE)
  testthat::local_mocked_bindings(ngs_genotype_error_matrix = function(...) E)
  calculate <- function() cc_ngs_power(
    N_case = 900, alpha = .05, prev = .05, pd = .3, R2 = 1.8,
    coverage = 10, seq_error = .01, verbose = FALSE
  )
  expect_warning(repaired <- calculate(), "All-zero genotype")
  E[3, 3] <- 1
  expect_silent(explicit <- calculate())
  expect_equal(repaired, explicit)
  expect_true(is.finite(repaired$power))
  g <- c(.6, .3, .1)
  expect_equal(cc_apply_genotype_misclass(g, E), as.numeric(t(E) %*% g))
})
