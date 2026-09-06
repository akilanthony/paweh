ngs_diff_args <- function(...) {
  c(list(alpha = 0.05, prev = 0.05, pd = 0.3, R2 = 1.8,
         coverage = 10, seq_error = 0.01, verbose = FALSE), list(...))
}

test_that("equal explicit mechanisms preserve all legacy fields exactly", {
  for (MOI in c("M", "D", "Rec")) {
    for (pi in c(1, 0.6)) {
      for (fun in list(cc_ngs_power, cc_ngs_mssn)) {
        size <- if (identical(fun, cc_ngs_power)) list(N_case = 1000) else
          list(power = 0.8)
        args <- c(size, ngs_diff_args(MOI = MOI, locus_het = TRUE, pi = pi))
        common <- do.call(fun, args)
        explicit <- do.call(fun, c(args, list(
          case_coverage = 10, ctrl_coverage = 10,
          case_seq_error = 0.01, ctrl_seq_error = 0.01
        )))
        expect_identical(explicit, common)
      }
    }
  }
})

test_that("differential designs match independent matrix compositions", {
  settings <- list(c(8, 8, 0.01, 0.08), c(4, 12, 0.01, 0.01),
                   c(4, 12, 0.08, 0.01), c(1, 4, 0, 0.499))
  for (v in settings) {
    for (MOI in c("M", "D", "Rec")) {
      for (pi in c(1, 0.6)) {
        args <- ngs_diff_args(MOI = MOI, k = 1.3, locus_het = TRUE, pi = pi,
                             case_coverage = v[1], ctrl_coverage = v[2],
                             case_seq_error = v[3], ctrl_seq_error = v[4])
        p <- do.call(cc_ngs_power, c(list(N_case = 700), args))
        m <- do.call(cc_ngs_mssn, c(list(power = 0.8), args))
        model <- .cc_model_genotype_frequencies(0.3, 1.8, MOI, 0.05)
        E_case <- ngs_genotype_error_matrix(v[1], v[3])
        E_ctrl <- ngs_genotype_error_matrix(v[2], v[4])
        gc <- as.numeric((pi * model$case + (1 - pi) * model$control) %*% E_case)
        gu <- as.numeric(model$control %*% E_ctrl)
        for (out in list(p, m)) {
          expect_equal(out$freqs$case_called, gc, tolerance = 1e-14)
          expect_equal(out$freqs$control_called, gu, tolerance = 1e-14)
          expect_identical(out$sequencing$case_transition_matrix, E_case)
          expect_identical(out$sequencing$ctrl_transition_matrix, E_ctrl)
          for (E in list(out$sequencing$case_transition_matrix,
                         out$sequencing$ctrl_transition_matrix)) {
            expect_true(all(is.finite(E) & E >= 0 & E <= 1))
            expect_equal(unname(rowSums(E)), rep(1, 3), tolerance = 1e-12)
          }
        }
        expected <- .cc_ahn_trend_ncp(gc, gu, 700, 910, p$scores)
        expect_equal(p$lambda, expected, tolerance = 1e-12)
        expect_equal(p$power, pchisq(qchisq(.95, 1), 1, expected,
                                    lower.tail = FALSE), tolerance = 1e-12)
        D <- sum(m$scores * (gc - gu))
        pooled <- gc + 1.3 * gu
        Q <- sum(m$scores^2 * pooled) - sum(m$scores * pooled)^2 / 2.3
        expect_equal(m$N_case_continuous, m$lambda_target * Q / (1.3 * D^2),
                     tolerance = 1e-12)
        design_power <- function(n) {
          lambda <- .cc_ahn_trend_ncp(gc, gu, n, ceiling(1.3 * n), m$scores)
          pchisq(qchisq(.95, 1), 1, lambda, lower.tail = FALSE)
        }
        expect_equal(m$achieved_power, design_power(m$MSSN_case))
        expect_gte(m$achieved_power, .8 - 1e-12)
        if (m$MSSN_case > 1) expect_lt(design_power(m$MSSN_case - 1), .8)
      }
    }
  }
})

test_that("each override falls back independently and common inputs may be omitted", {
  values <- list(case_coverage = 4, ctrl_coverage = 12,
                 case_seq_error = .08, ctrl_seq_error = .02)
  for (fun in list(cc_ngs_power, cc_ngs_mssn)) {
    size <- if (identical(fun, cc_ngs_power)) list(N_case = 1000) else list(power = .8)
    for (name in names(values)) {
      args <- c(size, ngs_diff_args(), values[name])
      out <- do.call(fun, args)
      expected <- list(case_coverage = 10, ctrl_coverage = 10,
                       case_seq_error = .01, ctrl_seq_error = .01)
      expected[name] <- values[name]
      expect_identical(out$sequencing[names(expected)], expected)
    }
    args <- c(size, ngs_diff_args(), values)
    out <- do.call(fun, args)
    args$coverage <- NULL
    args$seq_error <- NULL
    only <- do.call(fun, args)
    expect_identical(only$sequencing, out$sequencing)
    expect_identical(only$freqs, out$freqs)
    args$coverage <- "unused"
    args$seq_error <- "unused"
    expect_identical(do.call(fun, args)$sequencing, out$sequencing)
  }
})

test_that("invalid effective group parameters reuse sequencing validation", {
  for (fun in list(cc_ngs_power, cc_ngs_mssn)) {
    size <- if (identical(fun, cc_ngs_power)) list(N_case = 1000) else list(power = .8)
    for (name in c("case_coverage", "ctrl_coverage", "case_seq_error", "ctrl_seq_error")) {
      bad <- if (grepl("coverage", name)) list(0, -1, 1.5) else list(-.1, .5, 1)
      bad <- c(bad, list(NULL, NA_real_, NaN, Inf, "1", c(1, 2)))
      for (value in bad) {
        args <- c(size, ngs_diff_args(), setNames(list(value), name))
        group <- if (startsWith(name, "case")) "case" else "control"
        parameter <- if (grepl("coverage", name)) "coverage" else "seq_error"
        expect_error(do.call(fun, args), paste0(group, " sequencing: .*", parameter))
      }
    }
  }
})

test_that("differential mechanisms can create a nominal contrast at the biological null", {
  args <- ngs_diff_args(MOI = "D", locus_het = TRUE, pi = 0,
                       case_coverage = 2, ctrl_coverage = 20)
  p <- do.call(cc_ngs_power, c(list(N_case = 1000), args))
  expect_identical(p$freqs$case_true, p$freqs$control_true)
  expect_gt(p$lambda, 0)
  expect_gt(p$power, p$alpha)
  m <- do.call(cc_ngs_mssn, c(list(power = .8), args))
  expect_true(is.finite(m$MSSN_case))
  for (out in list(p, m)) {
    expect_output(print(out), "Cases coverage: 2")
    expect_output(print(out), "Controls coverage: 20")
  }
  args$verbose <- TRUE
  output <- capture.output(do.call(cc_ngs_power, c(list(N_case = 1000), args)),
                           type = "message")
  expect_true(any(grepl("Control coverage", output)))
  expect_true(any(grepl("Control Sequencing-Derived", output)))
})

test_that("original positional argument order is preserved", {
  p <- cc_ngs_power(1000, .05, .05, .3, 1.8, 10, .01, "D", 1.2, FALSE, TRUE, .6)
  m <- cc_ngs_mssn(.8, .05, .05, .3, 1.8, 10, .01, "D", 1.2, FALSE, TRUE, .6)
  args <- ngs_diff_args(MOI = "D", k = 1.2, locus_het = TRUE, pi = .6)
  expect_identical(p, do.call(cc_ngs_power, c(list(N_case = 1000), args)))
  expect_identical(m, do.call(cc_ngs_mssn, c(list(power = .8), args)))
})

test_that("pre-extension common-mechanism numerical fixtures remain unchanged", {
  # Captured from the repository implementation before this extension.
  old <- dget(test_path("fixtures", "cc-ngs-common-baseline.R"))
  for (i in seq_len(nrow(old))) {
    z <- old[i, ]
    args <- list(alpha = .05, prev = .05, pd = .3, R2 = 1.8, k = 1.2,
                 verbose = FALSE, locus_het = TRUE, MOI = z$MOI,
                 coverage = z$coverage, seq_error = z$seq_error, pi = z$pi)
    p <- do.call(cc_ngs_power, c(list(N_case = 1000), args))
    m <- do.call(cc_ngs_mssn, c(list(power = .8), args))
    for (name in c("lambda", "power")) {
      expect_equal(p[[name]], z[[name]], tolerance = 1e-14)
    }
    for (name in c("N_case_continuous", "achieved_power")) {
      expect_equal(m[[name]], z[[name]], tolerance = 1e-14)
    }
    expect_identical(m$MSSN_case, z$MSSN_case)
    expect_identical(m$MSSN_ctrl, z$MSSN_ctrl)
  }
})
