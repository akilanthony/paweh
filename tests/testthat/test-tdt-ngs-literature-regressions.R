# S11 literature-validation freeze for the TDT1-NGS implementation.
#
# Provenance tiers used below:
# 1. Kim (2015), Appendix B: lambda = N * delta^2 * I_eff, with I_eff
#    obtained from the nuisance-adjusted observed read-count information.
# 2. Gordon, Finch, and Kim (2020), Chapter 5: the 108-setting parameter
#    domain reconstructed here. The accessible source material did not expose
#    an exact numerical power table, so no computed value below is described
#    as a published table value.
# 3. Fixed numeric values are PAWEH/published-kernel regression fixtures. They
#    are equation-derived outputs of the frozen S8 kernel, not quotations of
#    textbook results.

s11_chapter5_grid <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    design <- expand.grid(
      pd = c(0.15, 0.325, 0.5),
      R1 = c(1.1, 1.2, 1.3),
      seq_error = c(0.005, 0.01),
      coverage = c(4L, 12L, 20L, 28L, 36L, 44L),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    design <- design[order(
      design$pd, design$R1, design$seq_error, design$coverage
    ), , drop = FALSE]
    rownames(design) <- NULL

    started <- proc.time()[["elapsed"]]
    fits <- lapply(seq_len(nrow(design)), function(i) {
      tdt_ngs_power(
        N = 5000,
        pd = design$pd[[i]],
        R1 = design$R1[[i]],
        coverage = design$coverage[[i]],
        seq_error = design$seq_error[[i]],
        alpha = 5e-8,
        verbose = FALSE
      )
    })
    elapsed <- proc.time()[["elapsed"]] - started

    cache <<- transform(
      design,
      R2 = R1^2,
      N = 5000,
      alpha = 5e-8,
      t = vapply(fits, `[[`, numeric(1), "t"),
      delta = vapply(fits, `[[`, numeric(1), "delta"),
      efficient_information = vapply(
        fits, `[[`, numeric(1), "efficient_information"
      ),
      lambda = vapply(fits, `[[`, numeric(1), "lambda"),
      power = vapply(fits, `[[`, numeric(1), "power")
    )
    attr(cache, "elapsed") <<- elapsed
    cache
  }
})

s11_grid_groups <- function(grid, fields) {
  split(grid, interaction(grid[fields], drop = TRUE, lex.order = TRUE))
}

test_that("the complete Chapter 5 parameter-domain reconstruction is valid", {
  expect_no_warning(grid <- s11_chapter5_grid())

  expect_equal(nrow(grid), 108L)
  expect_equal(
    names(grid),
    c(
      "pd", "R1", "seq_error", "coverage", "R2", "N", "alpha",
      "t", "delta", "efficient_information", "lambda", "power"
    )
  )
  expect_identical(unique(grid$N), 5000)
  expect_identical(unique(grid$alpha), 5e-8)
  expect_equal(sort(unique(grid$pd)), c(0.15, 0.325, 0.5))
  expect_equal(sort(unique(grid$R1)), c(1.1, 1.2, 1.3))
  expect_equal(sort(unique(grid$seq_error)), c(0.005, 0.01))
  expect_equal(sort(unique(grid$coverage)), c(4L, 12L, 20L, 28L, 36L, 44L))

  expect_true(all(is.finite(grid$efficient_information)))
  expect_true(all(grid$efficient_information > 0))
  expect_true(all(is.finite(grid$lambda)))
  expect_true(all(grid$lambda >= 0))
  expect_true(all(is.finite(grid$power)))
  expect_true(all(grid$power >= 0 & grid$power <= 1))
  expect_equal(grid$R2, grid$R1^2, tolerance = 1e-15)
  expect_equal(grid$t, grid$R1 / (1 + grid$R1), tolerance = 1e-15)
  expect_equal(grid$delta, log(grid$R1), tolerance = 2e-15)
  expect_true(is.finite(attr(grid, "elapsed")))
  expect_gt(attr(grid, "elapsed"), 0)
})

test_that("all 108 Chapter 5 settings are deterministic", {
  grid <- s11_chapter5_grid()

  repeated <- lapply(seq_len(nrow(grid)), function(i) {
    tdt_ngs_power(
      N = grid$N[[i]],
      pd = grid$pd[[i]],
      R1 = grid$R1[[i]],
      coverage = grid$coverage[[i]],
      seq_error = grid$seq_error[[i]],
      alpha = grid$alpha[[i]],
      verbose = FALSE
    )
  })
  expect_equal(
    vapply(repeated, `[[`, numeric(1), "efficient_information"),
    grid$efficient_information,
    tolerance = 1e-14
  )
  expect_equal(
    vapply(repeated, `[[`, numeric(1), "lambda"),
    grid$lambda,
    tolerance = 1e-13
  )
  expect_equal(
    vapply(repeated, `[[`, numeric(1), "power"),
    grid$power,
    tolerance = 1e-14
  )
})

test_that("the complete sorted grid has compact deterministic fingerprints", {
  grid <- s11_chapter5_grid()
  observed <- vapply(
    grid[c("efficient_information", "lambda", "power")],
    function(x) {
      c(
        min = min(x), max = max(x), median = stats::median(x),
        sum = sum(x), weighted = sum(seq_along(x) * x)
      )
    },
    numeric(5)
  )
  expected <- cbind(
    efficient_information = c(
      min = 0.087668776226323672,
      max = 0.24999989503950212,
      median = 0.21762937616783018,
      sum = 20.466289529119898,
      weighted = 1270.1049143208272
    ),
    lambda = c(
      min = 3.9819291306025248,
      max = 86.043722962176460,
      median = 33.911634502516470,
      sum = 3791.7276431670875,
      weighted = 259766.11701321282
    ),
    power = c(
      min = 0.00027429727333599030,
      max = 0.99993452502669244,
      median = 0.64286808003387641,
      sum = 53.621273548249874,
      weighted = 3678.0342417914781
    )
  )
  expect_equal(observed, expected, tolerance = 2e-10)
})

test_that("established S8 coverage fixtures remain locked", {
  grid <- s11_chapter5_grid()
  fixture <- subset(
    grid,
    pd == 0.325 & R1 == 1.2 & seq_error == 0.005 &
      coverage %in% c(4L, 12L, 20L, 44L)
  )
  fixture <- fixture[order(fixture$coverage), ]

  expect_equal(
    fixture$efficient_information,
    c(
      0.166296146079271,
      0.216362516007223,
      0.219166852960543,
      0.219374904964103
    ),
    tolerance = 2e-13
  )
  expect_equal(
    fixture$lambda,
    c(
      27.6393757408911,
      35.9606943225104,
      36.4267912500960,
      36.4613706894614
    ),
    tolerance = 2e-12
  )
})

test_that("established S8 sequencing-error fixtures remain locked", {
  expected <- data.frame(
    seq_error = c(0, 0.005, 0.01),
    efficient_information = c(
      0.219190883720878,
      0.216362516007223,
      0.214412722729257
    ),
    lambda = c(
      36.4307853006492,
      35.9606943225104,
      35.6366274677014
    )
  )

  for (i in seq_len(nrow(expected))) {
    fit <- tdt_ngs_power(
      N = 5000, pd = 0.325, R1 = 1.2,
      coverage = 12, seq_error = expected$seq_error[[i]],
      alpha = 5e-8, verbose = FALSE
    )
    expect_equal(
      fit$efficient_information,
      expected$efficient_information[[i]],
      tolerance = 2e-13
    )
    expect_equal(fit$lambda, expected$lambda[[i]], tolerance = 2e-12)
  }
})

test_that("selected powers independently match the one-df chi-square tail", {
  grid <- s11_chapter5_grid()
  selected <- grid[c(1L, 18L, 37L, 54L, 73L, 90L, 108L), ]
  critical <- stats::qchisq(1 - 5e-8, df = 1)
  expected <- stats::pchisq(
    critical,
    df = 1,
    ncp = selected$lambda,
    lower.tail = FALSE
  )
  expect_equal(selected$power, expected, tolerance = 1e-14)
})

test_that("coverage response is increasing and plateaus on this exact grid", {
  grid <- s11_chapter5_grid()
  groups <- s11_grid_groups(grid, c("pd", "R1", "seq_error"))
  differences <- lapply(groups, function(x) {
    x <- x[order(x$coverage), ]
    diff(x$power)
  })

  # This is a regression property of the reconstructed 108-setting domain,
  # not a claim of universal monotonicity for every possible TDT1-NGS design.
  expect_true(all(vapply(
    differences, function(x) all(x >= -1e-14), logical(1)
  )))
  expect_true(all(vapply(groups, function(x) {
    x$power[x$coverage == 44] >= x$power[x$coverage == 4]
  }, logical(1))))
  expect_true(all(vapply(groups, function(x) {
    low_gain <- x$power[x$coverage == 12] - x$power[x$coverage == 4]
    high_gain <- x$power[x$coverage == 44] - x$power[x$coverage == 36]
    low_gain > high_gain
  }, logical(1))))
})

test_that("relative-risk response is ordered on the Chapter 5 grid", {
  grid <- s11_chapter5_grid()
  groups <- s11_grid_groups(grid, c("pd", "coverage", "seq_error"))
  expect_true(all(vapply(groups, function(x) {
    x <- x[order(x$R1), ]
    all(diff(x$power) >= -1e-14)
  }, logical(1))))
})

test_that("the paired sequencing-error effect is quantified, not overstated", {
  grid <- s11_chapter5_grid()
  groups <- s11_grid_groups(grid, c("pd", "R1", "coverage"))
  differences <- vapply(groups, function(x) {
    x$power[x$seq_error == 0.005] - x$power[x$seq_error == 0.01]
  }, numeric(1))

  expect_true(all(differences >= -1e-14))
  expect_equal(min(differences), 1.9766465131354494e-09, tolerance = 2e-13)
  expect_equal(max(differences), 0.07120978181762605, tolerance = 2e-12)
  expect_equal(stats::median(differences), 2.8938827627882297e-05,
               tolerance = 2e-12)
  expect_equal(mean(differences), 0.0045439451293689632,
               tolerance = 2e-12)
})

test_that("disease-allele effects and reflection symmetry remain explicit", {
  grid <- s11_chapter5_grid()
  means <- tapply(grid$power, grid$pd, mean)
  expect_true(means[["0.15"]] < means[["0.325"]])
  expect_true(means[["0.325"]] < means[["0.5"]])

  for (pd in c(0.15, 0.325)) {
    left <- tdt_ngs_power(
      N = 5000, pd = pd, R1 = 1.2,
      coverage = 28, seq_error = 0.005,
      alpha = 5e-8, verbose = FALSE
    )
    right <- tdt_ngs_power(
      N = 5000, pd = 1 - pd, R1 = 1.2,
      coverage = 28, seq_error = 0.005,
      alpha = 5e-8, verbose = FALSE
    )
    expect_equal(left$lambda, right$lambda, tolerance = 2e-12)
    expect_equal(left$power, right$power, tolerance = 2e-13)
  }
})

test_that("representative models converge to stable high-depth limits", {
  designs <- list(
    list(pd = 0.15, R1 = 1.1, seq_error = 0.005),
    list(pd = 0.325, R1 = 1.2, seq_error = 0.005),
    list(pd = 0.5, R1 = 1.3, seq_error = 0.01)
  )

  for (design in designs) {
    fits <- lapply(c(20L, 28L, 36L, 44L, 100L), function(coverage) {
      do.call(
        tdt_ngs_power,
        c(
          list(N = 5000, coverage = coverage, alpha = 5e-8,
               verbose = FALSE),
          design
        )
      )
    })
    information <- vapply(
      fits, `[[`, numeric(1), "efficient_information"
    )
    lambda <- vapply(fits, `[[`, numeric(1), "lambda")
    power <- vapply(fits, `[[`, numeric(1), "power")

    expect_true(all(diff(information) >= -1e-13))
    expect_true(all(diff(lambda) >= -1e-12))
    expect_true(all(diff(power) >= -1e-13))
    expect_lt(abs(information[[5]] - information[[4]]),
              abs(information[[2]] - information[[1]]))
    expect_lt(abs(lambda[[5]] - lambda[[4]]),
              abs(lambda[[2]] - lambda[[1]]))
    expect_lt(abs(power[[5]] - power[[4]]),
              abs(power[[2]] - power[[1]]))
  }
})

test_that("selected Chapter 5 MSSNs are analytic, attainable, and minimal", {
  designs <- list(
    list(power = 0.80, pd = 0.15, R1 = 1.1,
         coverage = 4L, seq_error = 0.005),
    list(power = 0.90, pd = 0.325, R1 = 1.2,
         coverage = 12L, seq_error = 0.01),
    list(power = 0.80, pd = 0.5, R1 = 1.3,
         coverage = 20L, seq_error = 0.005),
    list(power = 0.90, pd = 0.325, R1 = 1.2,
         coverage = 44L, seq_error = 0.005)
  )

  for (design in designs) {
    mssn <- do.call(
      tdt_ngs_mssn,
      c(design, list(alpha = 5e-8, verbose = FALSE))
    )
    expected_continuous <- mssn$lambda_target /
      (log(design$R1)^2 * mssn$efficient_information)
    achieved <- do.call(
      tdt_ngs_power,
      c(
        list(N = mssn$MSSN_trios, alpha = 5e-8, verbose = FALSE),
        design[c("pd", "R1", "coverage", "seq_error")]
      )
    )

    expect_equal(mssn$N_trios_continuous, expected_continuous,
                 tolerance = 2e-10)
    expect_equal(mssn$MSSN_trios, ceiling(expected_continuous))
    expect_identical(mssn$rounding_adjustment, 0)
    expect_equal(mssn$achieved_power, achieved$power, tolerance = 1e-14)
    expect_gte(achieved$power, design$power)

    if (mssn$MSSN_trios > 1) {
      previous <- do.call(
        tdt_ngs_power,
        c(
          list(N = mssn$MSSN_trios - 1L, alpha = 5e-8,
               verbose = FALSE),
          design[c("pd", "R1", "coverage", "seq_error")]
        )
      )
      expect_lt(previous$power, design$power)
    }
  }
})

test_that("the literature freeze introduces no alternate statistical engine", {
  power_source <- paste(deparse(body(tdt_ngs_power)), collapse = "\n")
  mssn_source <- paste(deparse(body(tdt_ngs_mssn)), collapse = "\n")

  expect_match(power_source, ".tdt_ngs_ncp(", fixed = TRUE)
  expect_match(mssn_source, ".tdt_ngs_information(", fixed = TRUE)
  expect_match(mssn_source, ".tdt_ngs_target_ncp(", fixed = TRUE)
  expect_false(grepl("ngs_genotype_error_matrix", power_source, fixed = TRUE))
  expect_false(grepl("ngs_genotype_error_matrix", mssn_source, fixed = TRUE))
})
