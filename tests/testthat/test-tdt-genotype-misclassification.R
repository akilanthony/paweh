tdt_genotype_anchor <- list(
  pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
  delta_prime = 1, alpha = 0.05
)

test_that("Chapter 5 TDT genotype-error matrix is exact and validated", {
  identity <- .tdt_genotype_misclassification_matrix(0)
  expect_identical(unname(identity), diag(3))

  for (e in c(0, 0.0025, 0.025, 0.05, 0.5)) {
    M <- .tdt_genotype_misclassification_matrix(e)
    expect_equal(unname(rowSums(M)), rep(1, 3), tolerance = 1e-15)
    expect_true(all(M >= 0))
    expect_identical(unname(M[1, 3]), 0)
    expect_identical(unname(M[3, 1]), 0)
  }

  expect_error(.tdt_genotype_misclassification_matrix(-0.001), "0, 0.5")
  expect_error(.tdt_genotype_misclassification_matrix(0.501), "0, 0.5")
  expect_error(.tdt_genotype_misclassification_matrix(NA_real_), "0, 0.5")
})

test_that("15 by 27 trio transition kernel is stochastic", {
  for (e in c(0, 0.025, 0.5)) {
    kernel <- .tdt_genotype_transition_kernel(e)
    expect_identical(dim(kernel), c(15L, 27L))
    expect_true(all(is.finite(kernel)))
    expect_true(all(kernel >= 0))
    expect_equal(unname(rowSums(kernel)), rep(1, 15), tolerance = 1e-14)
  }
})

test_that("Mendelian consistency explicitly classifies all 27 trio states", {
  states <- .tdt_all_trio_states()
  observed <- .tdt_mendelian_consistent(
    states$father, states$mother, states$child
  )
  expected_labels <- c(
    "000", "010", "011", "021", "100", "101", "110", "111", "112",
    "121", "122", "201", "211", "212", "222"
  )

  expect_length(observed, 27)
  expect_equal(sum(observed), 15)
  expect_setequal(states$label[observed], expected_labels)
})

test_that("observed trio probabilities are finite, nonnegative, and normalized", {
  true_probability <- .tdt_genotype_true_trio_probabilities(
    pd = 0.3, prev = 0.05, R1 = 1.5, R2 = 2.25,
    delta_prime = 0.8
  )
  out <- .tdt_genotype_error_expected_counts(
    true_probability, e = 0.025, N = 500
  )

  expect_length(out$observed_probability, 27)
  expect_true(all(is.finite(out$observed_probability)))
  expect_true(all(out$observed_probability >= 0))
  expect_equal(sum(out$observed_probability), 1, tolerance = 1e-14)
  expect_equal(
    out$retained_probability + out$mendelian_inconsistent_probability,
    1, tolerance = 1e-14
  )
})

test_that("deterministic H1 e > 0 calculation matches an independent enumeration", {
  pd <- 0.30
  prev <- 0.05
  R1 <- 1.5
  R2 <- 2.25
  delta_prime <- 0.8
  e <- 0.025
  N <- 500
  alpha <- 0.05

  # Independent test-side enumeration of ordered two-locus parental
  # haplotypes, transmissions, genotype error, and retained observed trios.
  qd <- 1 - pd
  D <- delta_prime * pd * qd
  hap <- data.frame(
    disease = c(1L, 1L, 0L, 0L),
    marker = c(1L, 0L, 1L, 0L),
    probability = c(pd^2 + D, pd * qd - D, qd * pd - D, qd^2 + D)
  )
  Z <- qd^2 + 2 * pd * qd * R1 + pd^2 * R2
  penetrance <- c(prev / Z, R1 * prev / Z, R2 * prev / Z)
  true_reference <- numeric(15)
  names(true_reference) <- c(
    "222", "121", "211", "122", "212", "021", "201", "110",
    "111", "112", "010", "100", "011", "101", "000"
  )

  for (fh1 in 1:4) for (fh2 in 1:4) {
    for (mh1 in 1:4) for (mh2 in 1:4) {
      parental <- prod(hap$probability[c(fh1, fh2, mh1, mh2)])
      for (ft in c(fh1, fh2)) for (mt in c(mh1, mh2)) {
        label <- paste0(
          hap$marker[[fh1]] + hap$marker[[fh2]],
          hap$marker[[mh1]] + hap$marker[[mh2]],
          hap$marker[[ft]] + hap$marker[[mt]]
        )
        disease_genotype <- hap$disease[[ft]] + hap$disease[[mt]]
        true_reference[[label]] <- true_reference[[label]] +
          parental * 0.25 * penetrance[[disease_genotype + 1L]]
      }
    }
  }
  true_reference <- true_reference / sum(true_reference)

  M <- matrix(c(
    1 - e, e, 0,
    e, 1 - 2 * e, e,
    0, e, 1 - e
  ), nrow = 3, byrow = TRUE)
  observed_states <- expand.grid(
    father = 0:2, mother = 0:2, child = 0:2,
    KEEP.OUT.ATTRS = FALSE
  )
  observed_reference <- apply(observed_states, 1, function(observed) {
    sum(vapply(names(true_reference), function(label) {
      true <- as.integer(strsplit(label, "", fixed = TRUE)[[1]])
      true_reference[[label]] *
        M[true[[1]] + 1L, observed[[1]] + 1L] *
        M[true[[2]] + 1L, observed[[2]] + 1L] *
        M[true[[3]] + 1L, observed[[3]] + 1L]
    }, numeric(1)))
  })
  observed_labels <- apply(observed_states, 1, paste0, collapse = "")
  names(observed_reference) <- observed_labels

  test_mc <- function(father, mother, child) {
    fg <- if (father == 0) 0 else if (father == 2) 1 else 0:1
    mg <- if (mother == 0) 0 else if (mother == 2) 1 else 0:1
    child %in% as.vector(outer(fg, mg, `+`))
  }
  retained <- mapply(
    test_mc,
    observed_states$father,
    observed_states$mother,
    observed_states$child
  )
  retained_probability <- sum(observed_reference[retained])
  p_star <- observed_reference[retained] / retained_probability
  ET_reference <- N * sum(
    p_star[c("011", "101", "111", "112", "112", "122", "212")]
  )
  ENT_reference <- N * sum(
    p_star[c("010", "100", "111", "110", "110", "121", "211")]
  )
  lambda_reference <-
    (ET_reference - ENT_reference)^2 / (ET_reference + ENT_reference)
  power_reference <- stats::pchisq(
    stats::qchisq(1 - alpha, df = 1),
    df = 1, ncp = lambda_reference, lower.tail = FALSE
  )

  internal <- .tdt_genotype_error_model_counts(
    pd, prev, R1, R2, delta_prime, e, N
  )
  public <- tdt_power(
    N = N, pd = pd, prev = prev, R1 = R1, R2 = R2,
    delta_prime = delta_prime, alpha = alpha,
    effect = "genotype_misclassification",
    genotype_misclassification_rate = e,
    verbose = FALSE
  )

  expect_equal(
    unname(internal$true_probability), unname(true_reference),
    tolerance = 1e-14
  )
  expect_equal(
    unname(internal$observed_probability), unname(observed_reference),
    tolerance = 1e-14
  )
  expect_equal(
    internal$retained_probability, retained_probability, tolerance = 1e-14
  )
  expect_equal(
    public$ET$genotype_misclassification, ET_reference, tolerance = 1e-12
  )
  expect_equal(
    public$ENT$genotype_misclassification, ENT_reference, tolerance = 1e-12
  )
  expect_equal(
    public$lambda$genotype_misclassification,
    lambda_reference,
    tolerance = 1e-12
  )
  expect_equal(
    public$power$genotype_misclassification,
    power_reference,
    tolerance = 1e-12
  )
})

test_that("zero genotype error reproduces ordinary TDT counts, power, and MSSN", {
  power_args <- c(
    list(N = 500, verbose = FALSE),
    tdt_genotype_anchor
  )
  ordinary_power <- do.call(tdt_power, power_args)
  genotype_power <- do.call(
    tdt_power,
    c(power_args, list(
      effect = "genotype_misclassification",
      genotype_misclassification_rate = 0
    ))
  )

  expect_equal(
    genotype_power$ET$genotype_misclassification,
    ordinary_power$ET$no_error,
    tolerance = 1e-12
  )
  expect_equal(
    genotype_power$ENT$genotype_misclassification,
    ordinary_power$ENT$no_error,
    tolerance = 1e-12
  )
  expect_equal(
    genotype_power$power$genotype_misclassification,
    ordinary_power$power$no_error,
    tolerance = 1e-12
  )
  expect_equal(genotype_power$genotype_misclassification$retained_trio_probability, 1)
  expect_equal(
    genotype_power$genotype_misclassification$mendelian_inconsistent_probability,
    0, tolerance = 1e-14
  )

  mssn_args <- c(
    list(target_power = 0.8, verbose = FALSE),
    tdt_genotype_anchor[names(tdt_genotype_anchor) != "alpha"],
    list(alpha = tdt_genotype_anchor$alpha)
  )
  ordinary_mssn <- do.call(tdt_mssn, mssn_args)
  genotype_mssn <- do.call(
    tdt_mssn,
    c(mssn_args, list(
      effect = "genotype_misclassification",
      genotype_misclassification_rate = 0
    ))
  )
  expect_equal(
    genotype_mssn$N$genotype_misclassification,
    ordinary_mssn$N$no_error,
    tolerance = 1e-10
  )
  expect_identical(
    ceiling(genotype_mssn$N$genotype_misclassification),
    ceiling(ordinary_mssn$N$no_error)
  )
})

test_that("zero-error null diagnostic equals nominal alpha", {
  out <- .tdt_genotype_error_null_diagnostic(
    pd = 0.25, e = 0, N = 5000, alpha = 5e-8
  )
  expect_equal(out$ET0, out$ENT0, tolerance = 1e-12)
  expect_equal(out$lambda_null, 0, tolerance = 1e-14)
  expect_lt(abs(out$actual_alpha - 5e-8), 1e-14)
  expect_lt(abs(out$alpha_inflation_ratio - 1), 1e-8)
})

test_that("all nine Gordon, Finch, and Kim Table 5.12 values are reproduced", {
  # The source's 40.4039 inflation ratio implies log10 = 1.6064; its printed
  # 1.6061 value is internally inconsistent and is not encoded here.
  target <- data.frame(
    pd = rep(c(0.45, 0.25, 0.05), each = 3),
    e = rep(c(0.0025, 0.025, 0.05), 3),
    alpha_star = c(
      5.15e-8, 2.62e-7, 2.02e-6,
      1.01e-7, 0.0010, 0.1797,
      1.33e-6, 0.9074, 1.0000
    ),
    alpha_tolerance = c(
      0.01e-8, 0.01e-7, 0.01e-6,
      0.01e-7, 0.00005, 0.00005,
      0.01e-6, 0.00005, 0.00005
    ),
    log10_ratio = c(
      0.0129, 0.7189, 1.6064,
      0.3046, 4.2804, 6.5555,
      1.4234, 7.2588, 7.3010
    ),
    ratio = c(
      1.0301, 5.2348, 40.4039,
      2.0166, 19073, NA,
      26.5098, NA, 2e7
    )
  )

  for (i in seq_len(nrow(target))) {
    out <- .tdt_genotype_error_null_diagnostic(
      pd = target$pd[[i]], e = target$e[[i]],
      N = 5000, alpha = 5e-8
    )
    info <- paste("pd =", target$pd[[i]], "e =", target$e[[i]])
    expect_lt(
      abs(out$actual_alpha - target$alpha_star[[i]]),
      target$alpha_tolerance[[i]]
    )
    expect_equal(
      round(out$log10_alpha_inflation_ratio, 4),
      target$log10_ratio[[i]], info = info
    )
    if (!is.na(target$ratio[[i]])) {
      expect_lt(
        abs(out$alpha_inflation_ratio - target$ratio[[i]]),
        max(5e-5, 5e-5 * target$ratio[[i]])
      )
    }
  }
  expect_gt(
    .tdt_genotype_error_null_diagnostic(0.25, 0.05)$alpha_inflation_ratio,
    3.5e6
  )
  expect_gt(
    .tdt_genotype_error_null_diagnostic(0.05, 0.025)$alpha_inflation_ratio,
    18e6
  )
})

test_that("genotype-error MSSN inverts the canonical power calculation", {
  args <- list(
    pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
    alpha = 0.05, delta_prime = 0.8,
    effect = "genotype_misclassification",
    genotype_misclassification_rate = 0.025,
    verbose = FALSE
  )
  mssn <- do.call(tdt_mssn, c(list(target_power = 0.8), args))
  required <- ceiling(mssn$N$genotype_misclassification)
  at_required <- do.call(tdt_power, c(list(N = required), args))
  below_required <- do.call(tdt_power, c(list(N = required - 1), args))

  expect_gte(at_required$power$genotype_misclassification, 0.8)
  expect_lt(below_required$power$genotype_misclassification, 0.8)
})

test_that("effect selector rejects unsupported combinations and model-free use", {
  base <- list(
    N = 500, pd = 0.3, prev = 0.05, R1 = 1.5, R2 = 2.25,
    verbose = FALSE
  )
  expect_error(
    do.call(tdt_power, c(base, list(
      effect = "genotype_misclassification",
      genotype_misclassification_rate = 0.01,
      misclass_rate = 0.01
    ))),
    "cannot be combined"
  )
  expect_error(
    do.call(tdt_power, c(base, list(
      effect = "genotype_misclassification",
      genotype_misclassification_rate = 0.01,
      heter_rate = 0.01
    ))),
    "cannot be combined"
  )
  expect_error(
    do.call(tdt_power, c(base, list(
      genotype_misclassification_rate = 0.01
    ))),
    "Set effect"
  )
  expect_error(
    tdt_power(
      N = 500, input_mode = "model_free", ET = 140, ENT = 100,
      effect = "genotype_misclassification",
      genotype_misclassification_rate = 0.01,
      verbose = FALSE
    ),
    "requires input_mode = 'model_based'"
  )
})
