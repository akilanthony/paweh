cc_ngs_pheno_args <- function(...) {
  modifyList(
    list(
      alpha = 0.05, prev = 0.05, pd = 0.30, R2 = 1.8,
      coverage = 10, seq_error = 0.01, MOI = "M", k = 1,
      verbose = FALSE
    ),
    list(...)
  )
}

cc_ngs_pheno_ordinary <- function(args, N_case = 1000) {
  cc_power(
    N_case = N_case,
    alpha = args$alpha,
    input_mode = "model_based",
    prev = args$prev,
    pd = args$pd,
    R2 = args$R2,
    MOI = args$MOI,
    locus_het = if (is.null(args$locus_het)) FALSE else args$locus_het,
    pi = if (is.null(args$pi)) 1 else args$pi,
    pheno_misclass = args$pheno_misclass,
    theta = args$theta,
    phi = args$phi,
    k = args$k,
    w = .cc_ngs_scores_from_moi(args$MOI),
    geno_misclass = "none",
    verbose = FALSE
  )
}

test_that("omitted and inactive phenotype settings preserve CC-NGS exactly", {
  settings <- list(
    list(),
    list(locus_het = TRUE, pi = 0.6),
    list(case_coverage = 4, ctrl_coverage = 12),
    list(case_seq_error = 0.04, ctrl_seq_error = 0.005),
    list(locus_het = TRUE, pi = 0.6, case_coverage = 4,
         ctrl_coverage = 12, case_seq_error = 0.04,
         ctrl_seq_error = 0.005)
  )
  for (MOI in c("M", "D", "Rec")) {
    for (setting in settings) {
      args <- do.call(cc_ngs_pheno_args, c(list(MOI = MOI), setting))
      for (fun in list(cc_ngs_power, cc_ngs_mssn)) {
        size <- if (identical(fun, cc_ngs_power)) {
          list(N_case = 1000)
        } else {
          list(power = 0.8)
        }
        omitted <- do.call(fun, c(size, args))
        inactive <- do.call(
          fun,
          c(size, args, list(pheno_misclass = FALSE, theta = 0, phi = 0))
        )
        expect_identical(inactive, omitted)
      }
    }
  }
})

test_that("active zero phenotype error is a numerical identity", {
  args <- cc_ngs_pheno_args(
    locus_het = TRUE, pi = 0.6,
    case_coverage = 4, ctrl_coverage = 12,
    case_seq_error = 0.04, ctrl_seq_error = 0.005
  )
  for (fun in list(cc_ngs_power, cc_ngs_mssn)) {
    size <- if (identical(fun, cc_ngs_power)) {
      list(N_case = 1000)
    } else {
      list(power = 0.8)
    }
    none <- do.call(fun, c(size, args))
    zero <- do.call(
      fun,
      c(size, args, list(pheno_misclass = TRUE, theta = 0, phi = 0))
    )
    expect_identical(zero$freqs, none$freqs)
    expect_identical(zero$sequencing, none$sequencing)
    if (identical(fun, cc_ngs_power)) {
      expect_identical(zero$lambda, none$lambda)
      expect_identical(zero$power, none$power)
    } else {
      expect_identical(zero$N_case_continuous, none$N_case_continuous)
      expect_identical(zero$MSSN_case, none$MSSN_case)
      expect_identical(zero$MSSN_ctrl, none$MSSN_ctrl)
      expect_identical(zero$achieved_power, none$achieved_power)
    }
  }
})

test_that("CC-NGS pre-sequencing frequencies equal ordinary case-control", {
  for (MOI in c("M", "D", "Rec")) {
    for (setting in list(
      list(locus_het = FALSE, pi = 1),
      list(locus_het = TRUE, pi = 0.65)
    )) {
      args <- do.call(
        cc_ngs_pheno_args,
        c(list(MOI = MOI, pheno_misclass = TRUE,
               theta = 0.08, phi = 0.015), setting)
      )
      ngs <- do.call(cc_ngs_power, c(list(N_case = 1000), args))
      ordinary <- cc_ngs_pheno_ordinary(args)

      expect_identical(ngs$freqs$case_preseq,
                       ordinary$freqs$g_true_case)
      expect_identical(ngs$freqs$control_preseq,
                       ordinary$freqs$g_true_ctrl)
      expect_identical(
        ngs$errors$phenotype_misclass$g_case_after_pheno_misclass,
        ordinary$errors$phenotype_misclass$g_case_after_pheno_misclass
      )
      expect_identical(
        ngs$errors$phenotype_misclass$g_ctrl_after_pheno_misclass,
        ordinary$errors$phenotype_misclass$g_ctrl_after_pheno_misclass
      )
      expect_identical(
        ngs$errors$phenotype_misclass$case_denom,
        ordinary$errors$phenotype_misclass$case_denom
      )
      expect_identical(
        ngs$errors$phenotype_misclass$ctrl_denom,
        ordinary$errors$phenotype_misclass$ctrl_denom
      )
    }
  }
})

test_that("phenotype and sequencing composition matches manual matrices", {
  settings <- list(
    equal = list(),
    error = list(case_seq_error = 0.04, ctrl_seq_error = 0.005),
    coverage = list(case_coverage = 4, ctrl_coverage = 12),
    both = list(case_coverage = 4, ctrl_coverage = 12,
                case_seq_error = 0.04, ctrl_seq_error = 0.005),
    all_modifiers = list(
      locus_het = TRUE, pi = 0.6,
      case_coverage = 4, ctrl_coverage = 12,
      case_seq_error = 0.04, ctrl_seq_error = 0.005
    )
  )
  for (setting in settings) {
    args <- do.call(
      cc_ngs_pheno_args,
      c(list(MOI = "D", k = 1.3, pheno_misclass = TRUE,
             theta = 0.08, phi = 0.015), setting)
    )
    out <- do.call(cc_ngs_power, c(list(N_case = 700), args))
    E_case <- ngs_genotype_error_matrix(
      out$sequencing$case_coverage, out$sequencing$case_seq_error
    )
    E_ctrl <- ngs_genotype_error_matrix(
      out$sequencing$ctrl_coverage, out$sequencing$ctrl_seq_error
    )
    case_called <- as.numeric(t(E_case) %*% out$freqs$case_preseq)
    ctrl_called <- as.numeric(t(E_ctrl) %*% out$freqs$control_preseq)
    lambda <- .cc_ahn_trend_ncp(
      case_called, ctrl_called, 700, 910, out$scores
    )

    expect_equal(out$freqs$case_called, case_called, tolerance = 1e-15)
    expect_equal(out$freqs$control_called, ctrl_called, tolerance = 1e-15)
    expect_equal(out$lambda, lambda, tolerance = 1e-12)
    expect_equal(
      out$power,
      pchisq(qchisq(0.95, 1), 1, lambda, lower.tail = FALSE),
      tolerance = 1e-14
    )
  }
})

test_that("heterogeneity precedes phenotype mixing and sequencing", {
  args <- cc_ngs_pheno_args(
    MOI = "Rec", locus_het = TRUE, pi = 0.55,
    pheno_misclass = TRUE, theta = 0.12, phi = 0.02,
    case_coverage = 3, ctrl_coverage = 15,
    case_seq_error = 0.06, ctrl_seq_error = 0.002
  )
  out <- do.call(cc_ngs_power, c(list(N_case = 800), args))
  model <- .cc_model_genotype_frequencies(
    args$pd, args$R2, args$MOI, args$prev
  )
  case_het <- args$pi * model$case + (1 - args$pi) * model$control
  pheno <- .cc_apply_pheno_misclass(
    case_het, model$control, args$prev, args$theta, args$phi
  )
  E_case <- ngs_genotype_error_matrix(
    args$case_coverage, args$case_seq_error
  )
  E_ctrl <- ngs_genotype_error_matrix(
    args$ctrl_coverage, args$ctrl_seq_error
  )

  expect_equal(out$freqs$case_post_heterogeneity, case_het,
               tolerance = 1e-15)
  expect_identical(out$freqs$control_post_heterogeneity, model$control)
  expect_identical(out$freqs$case_preseq, pheno$g_case_obs)
  expect_identical(out$freqs$control_preseq, pheno$g_ctrl_obs)
  expect_equal(out$freqs$case_called,
               as.numeric(t(E_case) %*% pheno$g_case_obs),
               tolerance = 1e-15)
  expect_equal(out$freqs$control_called,
               as.numeric(t(E_ctrl) %*% pheno$g_ctrl_obs),
               tolerance = 1e-15)
})

test_that("phenotype-misclassified MSSNs attain target and are minimal", {
  settings <- list(
    equal = list(),
    differential = list(
      locus_het = TRUE, pi = 0.7,
      case_coverage = 4, ctrl_coverage = 15,
      case_seq_error = 0.04, ctrl_seq_error = 0.002
    )
  )
  for (setting in settings) {
    args <- do.call(
      cc_ngs_pheno_args,
      c(list(pheno_misclass = TRUE, theta = 0.05, phi = 0.01), setting)
    )
    mssn <- do.call(cc_ngs_mssn, c(list(power = 0.8), args))
    achieved <- do.call(
      cc_ngs_power,
      c(list(N_case = mssn$MSSN_case), args)
    )
    expect_equal(mssn$achieved_power, achieved$power, tolerance = 1e-14)
    expect_gte(achieved$power + 1e-12, 0.8)
    if (mssn$MSSN_case > 1) {
      previous <- do.call(
        cc_ngs_power,
        c(list(N_case = mssn$MSSN_case - 1), args)
      )
      expect_lt(previous$power, 0.8)
    }
  }
})

test_that("CC-NGS phenotype inputs follow ordinary case-control validation", {
  invalid_flags <- list(NA, 0, 1, "TRUE", c(TRUE, FALSE))
  invalid_rates <- list(-0.1, 1, NA_real_, NaN, Inf, "0.1", c(0.1, 0.2))
  for (fun in list(cc_ngs_power, cc_ngs_mssn)) {
    size <- if (identical(fun, cc_ngs_power)) {
      list(N_case = 1000)
    } else {
      list(power = 0.8)
    }
    for (value in invalid_flags) {
      expect_error(
        do.call(fun, c(size, cc_ngs_pheno_args(pheno_misclass = value))),
        "pheno_misclass must be TRUE or FALSE"
      )
    }
    for (name in c("theta", "phi")) {
      for (value in invalid_rates) {
        args <- cc_ngs_pheno_args(pheno_misclass = TRUE)
        args[[name]] <- value
        expect_error(
          do.call(fun, c(size, args)),
          paste0(name, " must be a single number in \\[0,1\\)")
        )
      }
    }
  }
})

test_that("high-depth phenotype CC-NGS agrees with ordinary trend design", {
  args <- cc_ngs_pheno_args(
    coverage = 100, seq_error = 0,
    locus_het = TRUE, pi = 0.7,
    pheno_misclass = TRUE, theta = 0.05, phi = 0.01
  )
  ngs <- do.call(cc_ngs_power, c(list(N_case = 1000), args))
  ordinary <- cc_ngs_pheno_ordinary(args)
  expect_equal(ngs$freqs$case_preseq, ordinary$freqs$g_true_case,
               tolerance = 1e-15)
  expect_equal(ngs$freqs$control_preseq, ordinary$freqs$g_true_ctrl,
               tolerance = 1e-15)
  expect_equal(ngs$lambda, ordinary$tests$trend$lambda, tolerance = 1e-12)
  expect_equal(ngs$power, ordinary$tests$trend$power, tolerance = 1e-12)
})

test_that("phenotype metadata and reporting are explicit only when active", {
  inactive <- do.call(
    cc_ngs_power,
    c(list(N_case = 1000), cc_ngs_pheno_args())
  )
  active <- do.call(
    cc_ngs_power,
    c(list(N_case = 1000), cc_ngs_pheno_args(
      pheno_misclass = TRUE, theta = 0.05, phi = 0.01
    ))
  )
  expect_false(inactive$errors$phenotype_misclass$enabled)
  expect_true(active$errors$phenotype_misclass$enabled)
  expect_identical(active$errors$phenotype_misclass$theta, 0.05)
  expect_identical(active$errors$phenotype_misclass$phi, 0.01)
  inactive_output <- capture.output(print(inactive))
  expect_false(any(grepl("Phenotype misclassification", inactive_output,
                         fixed = TRUE)))
  expect_output(print(active),
                "Phenotype misclassification: theta = 0.05; phi = 0.01")

  verbose <- capture.output(do.call(
    cc_ngs_power,
    c(list(N_case = 1000), cc_ngs_pheno_args(
      pheno_misclass = TRUE, theta = 0.05, phi = 0.01, verbose = TRUE
    ))
  ), type = "message")
  expect_true(any(grepl("Phenotype misclassification", verbose, fixed = TRUE)))
  expect_true(any(grepl("theta", verbose, fixed = TRUE)))
  expect_true(any(grepl("phi", verbose, fixed = TRUE)))
})

test_that("phenotype arguments append after the differential sequencing API", {
  expected_tail <- c(
    "case_coverage", "ctrl_coverage", "case_seq_error", "ctrl_seq_error",
    "pheno_misclass", "theta", "phi"
  )
  expect_identical(tail(names(formals(cc_ngs_power)), 7), expected_tail)
  expect_identical(tail(names(formals(cc_ngs_mssn)), 7), expected_tail)

  old_power <- cc_ngs_power(
    1000, 0.05, 0.05, 0.3, 1.8, 10, 0.01,
    "D", 1.2, FALSE, TRUE, 0.6, 4, 12, 0.04, 0.005
  )
  named_power <- do.call(cc_ngs_power, c(
    list(N_case = 1000),
    cc_ngs_pheno_args(
      MOI = "D", k = 1.2, locus_het = TRUE, pi = 0.6,
      case_coverage = 4, ctrl_coverage = 12,
      case_seq_error = 0.04, ctrl_seq_error = 0.005
    )
  ))
  expect_identical(old_power, named_power)

  old_mssn <- cc_ngs_mssn(
    0.8, 0.05, 0.05, 0.3, 1.8, 10, 0.01,
    "D", 1.2, FALSE, TRUE, 0.6, 4, 12, 0.04, 0.005
  )
  named_mssn <- do.call(cc_ngs_mssn, c(
    list(power = 0.8),
    cc_ngs_pheno_args(
      MOI = "D", k = 1.2, locus_het = TRUE, pi = 0.6,
      case_coverage = 4, ctrl_coverage = 12,
      case_seq_error = 0.04, ctrl_seq_error = 0.005
    )
  ))
  expect_identical(old_mssn, named_mssn)
})
