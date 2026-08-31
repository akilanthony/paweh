# Internal console-presentation helpers shared by the canonical study-design
# interfaces. These functions format quantities already computed by the public
# functions; they do not perform statistical calculations.

.paweh_console_width <- 74L

.paweh_console_fmt <- function(x, digits = 4L, scientific = FALSE,
                               integer = FALSE) {
  if (is.character(x)) {
    return(x)
  }
  if (length(x) != 1L) {
    return(paste(x, collapse = ", "))
  }
  if (is.na(x)) {
    return("NA")
  }
  if (!is.finite(x)) {
    return(as.character(x))
  }
  if (isTRUE(integer)) {
    return(formatC(x, format = "f", digits = 0L, big.mark = ","))
  }
  if (!scientific && abs(x) < 0.5 * 10^(-digits)) {
    x <- 0
  }
  formatC(
    x,
    format = if (scientific) "e" else "f",
    digits = digits
  )
}

.paweh_console_header <- function(title, subtitle) {
  line <- paste(rep("=", .paweh_console_width), collapse = "")
  message("\n", line)
  message(title)
  message(subtitle)
  message(line)
}

.paweh_console_section <- function(title) {
  message("\n", title)
}

.paweh_console_rule <- function() {
  message(paste(rep("-", .paweh_console_width), collapse = ""))
}

.paweh_console_parameter <- function(label, value, digits = 4L,
                                     scientific = FALSE,
                                     integer = FALSE) {
  shown <- .paweh_console_fmt(
    value,
    digits = digits,
    scientific = scientific,
    integer = integer
  )
  message(sprintf("%-50s %18s", paste0(label, ":"), shown))
}

.paweh_console_vector <- function(x, digits = 4L) {
  paste(vapply(
    x,
    .paweh_console_fmt,
    character(1),
    digits = digits
  ), collapse = ", ")
}

.paweh_console_genotype_table <- function(rows, digits = 5L) {
  message(sprintf(
    "%-27s %14s %14s %14s",
    "", "Genotype 0", "Genotype 1", "Genotype 2"
  ))
  for (row in rows) {
    shown <- vapply(
      row$values,
      .paweh_console_fmt,
      character(1),
      digits = if (is.null(row$digits)) digits else row$digits
    )
    message(sprintf(
      "%-27s %14s %14s %14s",
      paste0(row$label, ":"), shown[1L], shown[2L], shown[3L]
    ))
  }
}

.paweh_console_transition_matrix <- function(E, digits = 6L) {
  message("Rows = true genotype; columns = called genotype")
  message(sprintf(
    "%-20s %16s %16s %16s",
    "", "Called 0", "Called 1", "Called 2"
  ))
  for (i in seq_len(3L)) {
    shown <- vapply(
      E[i, ],
      .paweh_console_fmt,
      character(1),
      digits = digits
    )
    message(sprintf(
      "%-20s %16s %16s %16s",
      paste0("True genotype ", i - 1L, ":"),
      shown[1L], shown[2L], shown[3L]
    ))
  }
}

.paweh_input_mode_label <- function(input_mode) {
  if (identical(input_mode, "model_based")) {
    "Model-based genetic parameters"
  } else {
    "Model-free supplied frequencies/counts"
  }
}

.paweh_moi_label <- function(MOI) {
  switch(
    MOI,
    M = "Multiplicative",
    D = "Dominant",
    Rec = "Recessive",
    MOI
  )
}

.paweh_cc_modifier_flags <- function(x) {
  list(
    locus = isTRUE(x$locus_het$enabled) && x$locus_het$pi < 1,
    phenotype = isTRUE(x$errors$phenotype_misclass$enabled) &&
      (x$errors$phenotype_misclass$theta > 0 ||
         x$errors$phenotype_misclass$phi > 0),
    genotype = isTRUE(x$errors$genotype_misclass$enabled)
  )
}

.paweh_print_cc_model <- function(x) {
  .paweh_console_section("Genetic Model")
  .paweh_console_parameter(
    "Input mode", .paweh_input_mode_label(x$input_mode)
  )
  if (identical(x$input_mode, "model_based")) {
    .paweh_console_parameter("Disease prevalence", x$model_info$prev, 4L)
    .paweh_console_parameter(
      "Disease/risk allele frequency", x$model_info$pd, 4L
    )
    .paweh_console_parameter(
      "Mode of inheritance", .paweh_moi_label(x$model_info$MOI)
    )
    .paweh_console_parameter("Relative risk R1", x$model_info$R1, 4L)
    .paweh_console_parameter("Relative risk R2", x$model_info$R2, 4L)
    .paweh_console_parameter(
      "Penetrances (f0, f1, f2)",
      .paweh_console_vector(x$model_info$penetrances, 5L)
    )
  } else {
    .paweh_console_parameter(
      "Frequency source", "User-supplied case/control genotype frequencies"
    )
  }
  .paweh_console_parameter(
    "Trend weights",
    paste(.paweh_console_fmt(x$w, digits = 4L), collapse = ", ")
  )
}

.paweh_print_cc_modifiers <- function(x, flags) {
  if (!any(unlist(flags))) {
    return(invisible(NULL))
  }
  .paweh_console_section("Active Modifiers")
  if (flags$locus) {
    .paweh_console_parameter("Locus heterogeneity", "Enabled")
    .paweh_console_parameter("Locus-homogeneity fraction (pi)", x$locus_het$pi, 4L)
  }
  if (flags$phenotype) {
    pheno <- x$errors$phenotype_misclass
    .paweh_console_parameter("Phenotype misclassification", "Enabled")
    .paweh_console_parameter("Affected-to-control rate (theta)", pheno$theta, 4L)
    .paweh_console_parameter("Unaffected-to-case rate (phi)", pheno$phi, 4L)
  }
  if (flags$genotype) {
    error <- x$errors$genotype_misclass
    model <- switch(
      error$model,
      `1p_symmetric` = "1-parameter",
      `2p_hom_het` = "2-parameter",
      `3p_homhet_homhom` = "3-parameter",
      `diff3p_homhet_homhom` = "Differential 3-parameter",
      error$model
    )
    .paweh_console_parameter("Genotype misclassification", model)
    if (identical(error$model, "1p_symmetric")) {
      .paweh_console_parameter("Error parameter e", error$e, 4L)
    } else if (identical(error$model, "2p_hom_het")) {
      .paweh_console_parameter(
        "Error parameters (e1, e2)",
        .paweh_console_vector(c(error$e1, error$e2), 4L)
      )
    } else if (identical(error$model, "3p_homhet_homhom")) {
      .paweh_console_parameter(
        "Error parameters (e01, e02, e03)",
        .paweh_console_vector(c(error$e01, error$e02, error$e03), 4L)
      )
    } else if (identical(error$model, "diff3p_homhet_homhom")) {
      .paweh_console_parameter("Differential error source", error$diff_source)
      .paweh_console_parameter("Differential multiplier", error$diff_multiplier, 4L)
      .paweh_console_parameter(
        "Case errors (e01, e02, e03)",
        .paweh_console_vector(error$case_params, 4L)
      )
      .paweh_console_parameter(
        "Control errors (e01, e02, e03)",
        .paweh_console_vector(error$ctrl_params, 4L)
      )
    }
  }
  invisible(NULL)
}

.paweh_print_cc_frequencies <- function(x, flags) {
  .paweh_console_section("Genotype Frequencies")
  rows <- list(
    list(label = "No-error cases", values = x$freqs$g_base_case),
    list(label = "No-error controls", values = x$freqs$g_base_ctrl)
  )
  if (any(unlist(flags))) {
    rows <- c(rows, list(
      list(label = "Adjusted cases", values = x$freqs$g_obs_case),
      list(label = "Adjusted controls", values = x$freqs$g_obs_ctrl)
    ))
  }
  .paweh_console_genotype_table(rows)
  .paweh_console_parameter(
    "No-error risk-allele frequency (cases)",
    x$freqs$g_base_case[2L] / 2 + x$freqs$g_base_case[3L], 5L
  )
  .paweh_console_parameter(
    "No-error risk-allele frequency (controls)",
    x$freqs$g_base_ctrl[2L] / 2 + x$freqs$g_base_ctrl[3L], 5L
  )
  if (any(unlist(flags))) {
    .paweh_console_parameter(
      "Adjusted risk-allele frequency (cases)",
      x$freqs$g_obs_case[2L] / 2 + x$freqs$g_obs_case[3L], 5L
    )
    .paweh_console_parameter(
      "Adjusted risk-allele frequency (controls)",
      x$freqs$g_obs_ctrl[2L] / 2 + x$freqs$g_obs_ctrl[3L], 5L
    )
  }
}

.paweh_print_cc_power <- function(x, baseline_tests) {
  flags <- .paweh_cc_modifier_flags(x)
  adjusted <- any(unlist(flags))
  .paweh_console_header("PAWEH Case-Control Study", "Power Analysis")
  .paweh_console_section("Study Design")
  .paweh_console_parameter("Cases", x$N_case, integer = TRUE)
  .paweh_console_parameter("Controls", x$N_ctrl, integer = TRUE)
  .paweh_console_parameter("Total sample size", x$N_total, integer = TRUE)
  .paweh_console_parameter("Control-to-case ratio", x$k, 3L)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_print_cc_model(x)
  .paweh_print_cc_modifiers(x, flags)
  .paweh_console_rule()
  .paweh_console_section("No-Error Design")
  .paweh_console_parameter("Genotype-test NCP", baseline_tests$genotypes$lambda, 6L)
  .paweh_console_parameter("Genotype-test power", baseline_tests$genotypes$power, 6L)
  .paweh_console_parameter("Trend-test NCP", baseline_tests$trend$lambda, 6L)
  .paweh_console_parameter("Trend-test power", baseline_tests$trend$power, 6L)
  if (adjusted) {
    .paweh_console_rule()
    .paweh_console_section("Adjusted Design")
    .paweh_console_parameter("Genotype-test NCP", x$tests$genotypes$lambda, 6L)
    .paweh_console_parameter("Genotype-test power", x$tests$genotypes$power, 6L)
    .paweh_console_parameter(
      "Genotype-test absolute power loss",
      baseline_tests$genotypes$power - x$tests$genotypes$power, 6L
    )
    .paweh_console_parameter("Trend-test NCP", x$tests$trend$lambda, 6L)
    .paweh_console_parameter("Trend-test power", x$tests$trend$power, 6L)
    .paweh_console_parameter(
      "Trend-test absolute power loss",
      baseline_tests$trend$power - x$tests$trend$power, 6L
    )
  }
  .paweh_console_rule()
  .paweh_print_cc_frequencies(x, flags)
  .paweh_console_rule()
}

.paweh_console_inflation <- function(adjusted, baseline) {
  if (!is.finite(adjusted) || !is.finite(baseline) || baseline <= 0) {
    return("Not defined")
  }
  paste0(.paweh_console_fmt(100 * (adjusted / baseline - 1), 2L), "%")
}

.paweh_console_percent_or_undefined <- function(x) {
  if (length(x) != 1L || !is.finite(x)) {
    return("Not defined")
  }
  paste0(.paweh_console_fmt(x, 2L), "%")
}

.paweh_print_cc_mssn_test <- function(label, test) {
  .paweh_console_parameter(paste(label, "target NCP"), test$lambda_star, 6L)
  .paweh_console_parameter(paste(label, "required cases"), test$MSSN_case,
                           integer = TRUE)
  .paweh_console_parameter(paste(label, "required controls"), test$MSSN_ctrl,
                           integer = TRUE)
  .paweh_console_parameter(paste(label, "total MSSN"), test$MSSN_total,
                           integer = TRUE)
}

.paweh_print_cc_mssn <- function(x, baseline_tests) {
  flags <- .paweh_cc_modifier_flags(x)
  adjusted <- any(unlist(flags))
  .paweh_console_header(
    "PAWEH Case-Control Study", "Minimum Sample Size Necessary"
  )
  .paweh_console_section("Target Design")
  .paweh_console_parameter("Target power", x$target_power, 3L)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_console_parameter("Control-to-case ratio", x$k, 3L)
  .paweh_print_cc_model(x)
  .paweh_print_cc_modifiers(x, flags)
  .paweh_console_rule()
  .paweh_console_section("No-Error Design")
  .paweh_print_cc_mssn_test("Genotype test", baseline_tests$genotypes)
  .paweh_print_cc_mssn_test("Trend test", baseline_tests$trend)
  if (adjusted) {
    .paweh_console_rule()
    .paweh_console_section("Adjusted Design")
    .paweh_print_cc_mssn_test("Genotype test", x$tests$genotypes)
    .paweh_console_parameter(
      "Genotype-test MSSN inflation",
      .paweh_console_inflation(
        x$tests$genotypes$MSSN_total,
        baseline_tests$genotypes$MSSN_total
      )
    )
    .paweh_print_cc_mssn_test("Trend test", x$tests$trend)
    .paweh_console_parameter(
      "Trend-test MSSN inflation",
      .paweh_console_inflation(
        x$tests$trend$MSSN_total,
        baseline_tests$trend$MSSN_total
      )
    )
  }
  .paweh_console_rule()
  .paweh_print_cc_frequencies(x, flags)
  .paweh_console_rule()
}

.paweh_print_tdt_model <- function(x) {
  p <- x$model_parameters
  .paweh_console_section("Genetic Model")
  .paweh_console_parameter(
    "Input mode", .paweh_input_mode_label(x$input_mode)
  )
  if (identical(x$input_mode, "model_based")) {
    .paweh_console_parameter("Disease/risk allele frequency", p$pd, 4L)
    .paweh_console_parameter("Disease prevalence", p$prev, 4L)
    .paweh_console_parameter("Relative risk R1", p$R1, 4L)
    .paweh_console_parameter("Relative risk R2", p$R2, 4L)
  } else {
    .paweh_console_parameter(
      "Transmission source", "User-supplied expected ET and ENT"
    )
    if (!is.null(p$pd) && is.finite(p$pd)) {
      .paweh_console_parameter("Modeled-allele frequency", p$pd, 4L)
    }
    if (!is.null(p$prev) && is.finite(p$prev)) {
      .paweh_console_parameter("Disease prevalence for modifier", p$prev, 4L)
    }
  }
  .paweh_console_parameter("LD scale (delta_prime)", p$delta_prime, 4L)
}

.paweh_print_tdt_scenario_power <- function(x, scenario, title,
                                            show_loss = FALSE) {
  .paweh_console_section(title)
  .paweh_console_parameter("NCP", x$lambda[[scenario]], 6L)
  .paweh_console_parameter("Power", x$power[[scenario]], 6L)
  if (show_loss) {
    .paweh_console_parameter(
      "Absolute power loss", x$power_loss[[scenario]], 6L
    )
  }
  .paweh_console_parameter("Expected transmitted count (ET)", x$ET[[scenario]], 3L)
  .paweh_console_parameter(
    "Expected non-transmitted count (ENT)", x$ENT[[scenario]], 3L
  )
  .paweh_console_parameter("Expected transmitted probability (gT*)",
                           x$gT_star[[scenario]], 6L)
  .paweh_console_parameter("Expected non-transmitted probability (gNT*)",
                           x$gNT_star[[scenario]], 6L)
}

.paweh_print_tdt_power <- function(x) {
  p <- x$model_parameters
  .paweh_console_header(
    "PAWEH Transmission Disequilibrium Test", "Power Analysis"
  )
  .paweh_console_section("Study Design")
  .paweh_console_parameter("Affected-child trios", x$N, integer = TRUE)
  .paweh_console_parameter("Total individuals", 3 * x$N, integer = TRUE)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_print_tdt_model(x)
  .paweh_console_rule()
  .paweh_print_tdt_scenario_power(x, "no_error", "No-Error Design")
  if (p$misclass_rate > 0) {
    .paweh_console_rule()
    .paweh_console_section("Phenotype Misclassification")
    .paweh_console_parameter("Misclassification rate", p$misclass_rate, 4L)
    .paweh_console_parameter("NCP", x$lambda$misclassification, 6L)
    .paweh_console_parameter("Power", x$power$misclassification, 6L)
    .paweh_console_parameter("Absolute power loss",
                             x$power_loss$misclassification, 6L)
    .paweh_console_parameter("Expected transmitted count (ET)",
                             x$ET$misclassification, 3L)
    .paweh_console_parameter("Expected non-transmitted count (ENT)",
                             x$ENT$misclassification, 3L)
    .paweh_console_parameter("Expected transmitted probability (gT*)",
                             x$gT_star$misclassification, 6L)
    .paweh_console_parameter("Expected non-transmitted probability (gNT*)",
                             x$gNT_star$misclassification, 6L)
  }
  if (p$heter_rate > 0) {
    .paweh_console_rule()
    .paweh_console_section("Locus Heterogeneity")
    .paweh_console_parameter("Heterogeneity rate", p$heter_rate, 4L)
    .paweh_console_parameter("NCP", x$lambda$heterogeneity, 6L)
    .paweh_console_parameter("Power", x$power$heterogeneity, 6L)
    .paweh_console_parameter("Absolute power loss",
                             x$power_loss$heterogeneity, 6L)
    .paweh_console_parameter("Expected transmitted count (ET)",
                             x$ET$heterogeneity, 3L)
    .paweh_console_parameter("Expected non-transmitted count (ENT)",
                             x$ENT$heterogeneity, 3L)
    .paweh_console_parameter("Expected transmitted probability (gT*)",
                             x$gT_star$heterogeneity, 6L)
    .paweh_console_parameter("Expected non-transmitted probability (gNT*)",
                             x$gNT_star$heterogeneity, 6L)
  }
  .paweh_console_rule()
}

.paweh_print_tdt_mssn_scenario <- function(x, scenario, title,
                                           show_increase = FALSE) {
  .paweh_console_section(title)
  .paweh_console_parameter("Continuous trio requirement", x$N[[scenario]], 3L)
  .paweh_console_parameter("Required complete trios", ceiling(x$N[[scenario]]),
                           integer = TRUE)
  .paweh_console_parameter("Total individuals", 3 * ceiling(x$N[[scenario]]),
                           integer = TRUE)
  if (show_increase) {
    .paweh_console_parameter(
      "MSSN increase",
      .paweh_console_percent_or_undefined(x$percent_increase[[scenario]])
    )
    .paweh_console_parameter(
      "Power at no-error trio requirement",
      x$power_at_N_no_error[[scenario]], 6L
    )
    .paweh_console_parameter(
      "Absolute power loss at no-error requirement",
      x$power_loss_at_N_no_error[[scenario]], 6L
    )
  }
  .paweh_console_parameter("Expected transmitted probability (gT*)",
                           x$gT_star[[scenario]], 6L)
  .paweh_console_parameter("Expected non-transmitted probability (gNT*)",
                           x$gNT_star[[scenario]], 6L)
}

.paweh_print_tdt_mssn <- function(x) {
  p <- x$model_parameters
  .paweh_console_header(
    "PAWEH Transmission Disequilibrium Test",
    "Minimum Sample Size Necessary"
  )
  .paweh_console_section("Target Design")
  .paweh_console_parameter("Target power", x$target_power, 3L)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_console_parameter("Target NCP", x$lambda_star, 6L)
  .paweh_print_tdt_model(x)
  .paweh_console_rule()
  .paweh_print_tdt_mssn_scenario(x, "no_error", "No-Error Design")
  if (p$misclass_rate > 0) {
    .paweh_console_rule()
    .paweh_console_section("Phenotype Misclassification")
    .paweh_console_parameter("Misclassification rate", p$misclass_rate, 4L)
    .paweh_console_parameter("Continuous trio requirement",
                             x$N$misclassification, 3L)
    .paweh_console_parameter("Required complete trios",
                             ceiling(x$N$misclassification), integer = TRUE)
    .paweh_console_parameter("Total individuals",
                             3 * ceiling(x$N$misclassification), integer = TRUE)
    .paweh_console_parameter(
      "MSSN increase",
      .paweh_console_percent_or_undefined(x$percent_increase$misclassification)
    )
    .paweh_console_parameter("Power at no-error trio requirement",
                             x$power_at_N_no_error$misclassification, 6L)
    .paweh_console_parameter("Absolute power loss",
                             x$power_loss_at_N_no_error$misclassification, 6L)
    .paweh_console_parameter("Expected transmitted probability (gT*)",
                             x$gT_star$misclassification, 6L)
    .paweh_console_parameter("Expected non-transmitted probability (gNT*)",
                             x$gNT_star$misclassification, 6L)
  }
  if (p$heter_rate > 0) {
    .paweh_console_rule()
    .paweh_console_section("Locus Heterogeneity")
    .paweh_console_parameter("Heterogeneity rate", p$heter_rate, 4L)
    .paweh_console_parameter("Continuous trio requirement",
                             x$N$heterogeneity, 3L)
    .paweh_console_parameter("Required complete trios",
                             ceiling(x$N$heterogeneity), integer = TRUE)
    .paweh_console_parameter("Total individuals",
                             3 * ceiling(x$N$heterogeneity), integer = TRUE)
    .paweh_console_parameter(
      "MSSN increase",
      .paweh_console_percent_or_undefined(x$percent_increase$heterogeneity)
    )
    .paweh_console_parameter("Power at no-error trio requirement",
                             x$power_at_N_no_error$heterogeneity, 6L)
    .paweh_console_parameter("Absolute power loss",
                             x$power_loss_at_N_no_error$heterogeneity, 6L)
    .paweh_console_parameter("Expected transmitted probability (gT*)",
                             x$gT_star$heterogeneity, 6L)
    .paweh_console_parameter("Expected non-transmitted probability (gNT*)",
                             x$gNT_star$heterogeneity, 6L)
  }
  .paweh_console_rule()
}

.paweh_print_cc_ngs_model <- function(x) {
  .paweh_console_section("Genetic Model")
  .paweh_console_parameter("Disease prevalence", x$model_info$prev, 4L)
  .paweh_console_parameter("Disease/risk allele frequency", x$model_info$pd, 4L)
  .paweh_console_parameter("Mode of inheritance", .paweh_moi_label(x$MOI))
  .paweh_console_parameter("Relative risk R1", x$model_info$R1, 4L)
  .paweh_console_parameter("Relative risk R2", x$model_info$R2, 4L)
  .paweh_console_parameter(
    "Trend scores", .paweh_console_vector(x$scores, 2L)
  )
  if (isTRUE(x$locus_het$enabled)) {
    .paweh_console_parameter("Locus heterogeneity", "Enabled")
    .paweh_console_parameter("Locus-homogeneity fraction (pi)",
                             x$locus_het$pi, 4L)
  }
}

.paweh_print_cc_ngs_sequencing <- function(x) {
  .paweh_console_section("Sequencing Model")
  .paweh_console_parameter("Fixed coverage", x$coverage, integer = TRUE)
  .paweh_console_parameter("Per-read sequencing error", x$seq_error, 5L)
  .paweh_console_parameter(
    "Genotype calling", "Deterministic maximum-likelihood calling"
  )
  .paweh_console_parameter(
    "Observation model", "Sequencing-derived called genotypes"
  )
}

.paweh_print_cc_ngs_frequencies <- function(x) {
  .paweh_console_section("True Genotype Frequencies")
  .paweh_console_genotype_table(list(
    list(label = "Cases", values = x$freqs$case_true),
    list(label = "Controls", values = x$freqs$control_true)
  ))
  .paweh_console_rule()
  .paweh_console_section("Sequencing-Derived Genotype-Call Matrix")
  .paweh_console_transition_matrix(x$transition_matrix)
  .paweh_console_rule()
  .paweh_console_section("Called Genotype Frequencies")
  .paweh_console_genotype_table(list(
    list(label = "Cases", values = x$freqs$case_called),
    list(label = "Controls", values = x$freqs$control_called)
  ))
}

.paweh_print_cc_ngs_power <- function(x) {
  .paweh_console_header("PAWEH Case-Control NGS Study", "Power Analysis")
  .paweh_console_section("Study Design")
  .paweh_console_parameter("Cases", x$N_case, integer = TRUE)
  .paweh_console_parameter("Controls", x$N_ctrl, integer = TRUE)
  .paweh_console_parameter("Total sample size", x$N_total, integer = TRUE)
  .paweh_console_parameter("Control-to-case ratio", x$k, 3L)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_print_cc_ngs_model(x)
  .paweh_print_cc_ngs_sequencing(x)
  .paweh_console_rule()
  .paweh_print_cc_ngs_frequencies(x)
  .paweh_console_rule()
  .paweh_console_section("Power Result")
  .paweh_console_parameter("Noncentrality parameter (NCP)", x$lambda, 6L)
  .paweh_console_parameter("Power", x$power, 6L)
  .paweh_console_rule()
}

.paweh_print_cc_ngs_mssn <- function(x) {
  .paweh_console_header(
    "PAWEH Case-Control NGS Study", "Minimum Sample Size Necessary"
  )
  .paweh_console_section("Target Design")
  .paweh_console_parameter("Target power", x$power_target, 3L)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_console_parameter("Control-to-case ratio", x$k, 3L)
  .paweh_print_cc_ngs_model(x)
  .paweh_print_cc_ngs_sequencing(x)
  .paweh_console_rule()
  .paweh_print_cc_ngs_frequencies(x)
  .paweh_console_rule()
  .paweh_console_section("Required Sample Size")
  .paweh_console_parameter("Continuous case requirement",
                           x$N_case_continuous, 3L)
  .paweh_console_parameter("Required cases", x$MSSN_case, integer = TRUE)
  .paweh_console_parameter("Required controls", x$MSSN_ctrl, integer = TRUE)
  .paweh_console_parameter("Total MSSN", x$MSSN_total, integer = TRUE)
  .paweh_console_parameter("Target NCP", x$lambda_target, 6L)
  .paweh_console_parameter("Achieved NCP", x$achieved_lambda, 6L)
  .paweh_console_parameter("Achieved power", x$achieved_power, 6L)
  .paweh_console_rule()
}

.paweh_print_tdt_ngs_model <- function(x) {
  .paweh_console_section("Genetic Model")
  .paweh_console_parameter("Disease/risk allele frequency", x$pd, 4L)
  .paweh_console_parameter("Relative risk R1", x$R1, 4L)
  .paweh_console_parameter("Relative risk R2", x$R2, 4L)
  .paweh_console_parameter("Inheritance model", "Multiplicative")
  .paweh_console_parameter("Transmission parameter (t)", x$t, 6L)
  .paweh_console_parameter("Transmission log-odds (delta)", x$delta, 6L)
}

.paweh_print_tdt_ngs_sequencing <- function(x) {
  .paweh_console_section("Sequencing Model")
  .paweh_console_parameter("Equal fixed coverage", x$coverage, integer = TRUE)
  .paweh_console_parameter("Per-read sequencing error", x$seq_error, 5L)
  .paweh_console_parameter("Analysis", "Raw-read TDT1-NGS likelihood")
}

.paweh_print_tdt_ngs_power <- function(x) {
  .paweh_console_header("PAWEH TDT1-NGS", "Power Analysis")
  .paweh_console_section("Study Design")
  .paweh_console_parameter("Complete affected-child trios", x$N, integer = TRUE)
  .paweh_console_parameter("Total individuals", 3 * x$N, integer = TRUE)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_print_tdt_ngs_model(x)
  .paweh_print_tdt_ngs_sequencing(x)
  .paweh_console_rule()
  .paweh_console_section("Information and Power")
  .paweh_console_parameter("Efficient information per trio",
                           x$efficient_information, 8L)
  .paweh_console_parameter("Noncentrality parameter (NCP)", x$lambda, 6L)
  .paweh_console_parameter("Power", x$power, 6L)
  .paweh_console_rule()
}

.paweh_print_tdt_ngs_mssn <- function(x) {
  .paweh_console_header("PAWEH TDT1-NGS", "Minimum Sample Size Necessary")
  .paweh_console_section("Target Design")
  .paweh_console_parameter("Target power", x$power_target, 3L)
  .paweh_console_parameter("Significance level (alpha)", x$alpha, 2L, TRUE)
  .paweh_console_rule()
  .paweh_console_section("Required Sample Size")
  .paweh_console_parameter("Continuous trio requirement",
                           x$N_trios_continuous, 3L)
  .paweh_console_parameter("Required complete trios", x$MSSN_trios,
                           integer = TRUE)
  .paweh_console_parameter("Total individuals", x$total_individuals,
                           integer = TRUE)
  .paweh_console_parameter("Target NCP", x$lambda_target, 6L)
  .paweh_console_parameter("Achieved NCP", x$achieved_lambda, 6L)
  .paweh_console_parameter("Achieved power", x$achieved_power, 6L)
  .paweh_print_tdt_ngs_model(x)
  .paweh_print_tdt_ngs_sequencing(x)
  .paweh_console_rule()
  .paweh_console_section("Information")
  .paweh_console_parameter("Efficient information per trio",
                           x$efficient_information, 8L)
  .paweh_console_parameter("NCP per trio", x$ncp_per_trio, 8L)
  .paweh_console_rule()
}
