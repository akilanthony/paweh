# Quantitative-trait dashboard; statistical work delegates to canonical QTL functions.

.paweh_qtl_defaults <- function() {
  list(
    subtype = "continuous", objective = "power", N = 500, N_case = 150,
    target_power = 80, alpha = 0.05, pd = 30, qtl_var = 10, tau = 0,
    count_method = "rounded", multiple_of_three = TRUE,
    x_upper = 10, x_lower = 10, k = 1,
    n_traits = 2, mv_test = "pillai",
    mv_qtl_var_1 = 10, mv_qtl_var_2 = 5, mv_qtl_var_3 = 3, mv_qtl_var_4 = 2,
    mv_tau_1 = 0, mv_tau_2 = 0.5, mv_tau_3 = 0, mv_tau_4 = 0,
    mv_x_upper_1 = 10, mv_x_upper_2 = 10, mv_x_upper_3 = 10, mv_x_upper_4 = 10,
    mv_x_lower_1 = 10, mv_x_lower_2 = 10, mv_x_lower_3 = 10, mv_x_lower_4 = 10,
    corr_1_2 = 0, corr_1_3 = 0, corr_1_4 = 0,
    corr_2_3 = 0, corr_2_4 = 0, corr_3_4 = 0
  )
}

.paweh_qtl_values <- function(input) {
  shiny::reactiveValuesToList(input)
  defaults <- .paweh_qtl_defaults()
  values <- lapply(names(defaults), function(name) {
    if (is.null(input[[name]])) defaults[[name]] else input[[name]]
  })
  stats::setNames(values, names(defaults))
}

.paweh_qtl_num <- function(x, label, lower, upper, open = FALSE) {
  invalid <- !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    x < lower || x > upper || (open && (x == lower || x == upper))
  if (invalid) stop(label, " is outside its allowed range.", call. = FALSE)
}

.paweh_qtl_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    class = "paweh-workspace paweh-qtl-workspace",
    .paweh_page_heading(
      "Quantitative Trait study design",
      "Design continuous, extreme-phenotype, or joint multiple-trait studies."
    ),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = "Design setup", open = "desktop", width = 280, resizable = FALSE,
        shiny::radioButtons(ns("subtype"), "How will the phenotype be analyzed?", c(
          `Full continuous trait` = "continuous",
          `Extreme phenotype sampling` = "extreme",
          `Multiple quantitative traits` = "multivariate"
        )),
        shiny::uiOutput(ns("subtype_help")),
        shiny::radioButtons(ns("objective"), "Objective", c(
          `Estimate power` = "power", `Minimum sample size` = "mssn"
        )),
        shiny::uiOutput(ns("objective_inputs")),
        shiny::uiOutput(ns("core_inputs")),
        shiny::numericInput(ns("alpha"), "Significance level (alpha)", 0.05, 1e-8, 0.999999),
        shiny::tags$details(
          class = "paweh-sidebar-section",
          shiny::tags$summary(shiny::strong("Advanced assumptions")),
          shiny::uiOutput(ns("advanced_inputs"))
        ),
        shiny::actionButton(ns("calculate"), "Calculate study design", class = "btn-primary paweh-calculate"),
        shiny::uiOutput(ns("changed_notice")),
        bslib::card(
          class = "paweh-design-summary", bslib::card_header("Your calculated design"),
          bslib::card_body(shiny::uiOutput(ns("design_summary")))
        )
      ),
      bslib::navset_card_tab(
        id = ns("section"),
        bslib::nav_panel("Results", shiny::uiOutput(ns("results"))),
        bslib::nav_panel(
          "Sensitivity", shiny::uiOutput(ns("sensitivity_controls")),
          shiny::uiOutput(ns("sensitivity_message")),
          shiny::uiOutput(ns("sensitivity_plot_container"))
        ),
        bslib::nav_panel(
          "Visualize", shiny::uiOutput(ns("visualize_intro")),
          shiny::uiOutput(ns("visualization_controls")),
          shiny::uiOutput(ns("visualization_plot_container")),
          shiny::tags$details(
            class = "paweh-advanced-visualization",
            shiny::tags$summary("Advanced visualization"),
            shiny::div(
              class = "paweh-advanced-body",
              shiny::uiOutput(ns("surface_controls")),
              shiny::uiOutput(ns("surface_message")),
              shiny::uiOutput(ns("surface_container"))
            )
          )
        ),
        bslib::nav_panel("Methods", shiny::uiOutput(ns("methods")))
      )
    )
  )
}

.paweh_qtl_correlation_matrix <- function(values, p) {
  matrix_value <- diag(p)
  if (p > 1L) for (i in seq_len(p - 1L)) for (j in seq.int(i + 1L, p)) {
    value <- values[[paste0("corr_", i, "_", j)]]
    .paweh_qtl_num(value, paste("Correlation for traits", i, "and", j), -1, 1)
    matrix_value[i, j] <- matrix_value[j, i] <- value
  }
  dimnames(matrix_value) <- list(paste("Trait", seq_len(p)), paste("Trait", seq_len(p)))
  matrix_value
}

.paweh_qtl_snapshot <- function(values) {
  if (!values$subtype %in% c("continuous", "extreme", "multivariate")) {
    stop("Choose a valid quantitative-trait workflow.", call. = FALSE)
  }
  if (!values$objective %in% c("power", "mssn")) stop("Choose a valid objective.", call. = FALSE)
  .paweh_qtl_num(values$alpha, "Significance level", 0, 1, TRUE)
  .paweh_qtl_num(values$pd, "Modeled-allele frequency", 0, 100, TRUE)
  if (values$objective == "mssn") {
    .paweh_qtl_num(values$target_power, "Target power", 0, 100, TRUE)
    objective_value <- values$target_power / 100
  }

  if (values$subtype == "continuous") {
    .paweh_qtl_num(values$qtl_var, "QTL variance explained", 0, 100, TRUE)
    .paweh_qtl_num(values$tau, "Dominance parameter", -Inf, Inf)
    if (values$objective == "power") {
      .paweh_qtl_num(values$N, "Sample size", 3, Inf, TRUE)
      if (values$N != floor(values$N)) stop("Sample size must be an integer.", call. = FALSE)
      objective_value <- values$N
    }
    args <- list(
      alpha = values$alpha, qtl_var = values$qtl_var / 100,
      tau = values$tau, pd = values$pd / 100,
      count_method = values$count_method
    )
    if (values$objective == "mssn") args$multiple_of_three <- isTRUE(values$multiple_of_three)
  } else if (values$subtype == "extreme") {
    .paweh_qtl_num(values$qtl_var, "QTL variance explained", 0, 100, TRUE)
    .paweh_qtl_num(values$tau, "Dominance parameter", -Inf, Inf)
    .paweh_qtl_num(values$x_upper, "Upper population tail", 0, 100, TRUE)
    .paweh_qtl_num(values$x_lower, "Lower population tail", 0, 100, TRUE)
    .paweh_qtl_num(values$k, "Lower-to-upper selected ratio", 0, Inf, TRUE)
    qtl_falconer_threshold_parameters(
      qtl_var = values$qtl_var / 100, tau = values$tau,
      pd = values$pd / 100, x_upper = values$x_upper,
      x_lower = values$x_lower, verbose = FALSE
    )
    if (values$objective == "power") {
      .paweh_qtl_num(values$N_case, "Upper-tail selected sample", 0, Inf, TRUE)
      objective_value <- values$N_case
    }
    args <- list(
      alpha = values$alpha, qtl_var = values$qtl_var / 100,
      tau = values$tau, pd = values$pd / 100,
      x_upper = values$x_upper, x_lower = values$x_lower, k = values$k
    )
  } else {
    p <- as.integer(values$n_traits)
    if (!p %in% 2:4) stop("Number of traits must be between 2 and 4.", call. = FALSE)
    qtl_var <- tau <- x_upper <- x_lower <- numeric(p)
    for (i in seq_len(p)) {
      qtl_var[i] <- values[[paste0("mv_qtl_var_", i)]] / 100
      tau[i] <- values[[paste0("mv_tau_", i)]]
      x_upper[i] <- values[[paste0("mv_x_upper_", i)]]
      x_lower[i] <- values[[paste0("mv_x_lower_", i)]]
      .paweh_qtl_num(qtl_var[i], paste("Trait", i, "QTL variance"), 0, 1, TRUE)
      .paweh_qtl_num(tau[i], paste("Trait", i, "dominance parameter"), -Inf, Inf)
    }
    cor_matrix <- .paweh_qtl_correlation_matrix(values, p)
    .falconer_mv_validate_model(qtl_var, tau, values$pd / 100, cor_matrix)
    if (!values$mv_test %in% c("pillai", "threshold_chisq")) stop("Choose a valid joint test.", call. = FALSE)
    if (values$mv_test == "threshold_chisq") {
      for (i in seq_len(p)) {
        .paweh_qtl_num(x_upper[i], paste("Trait", i, "upper tail"), 0, 100, TRUE)
        .paweh_qtl_num(x_lower[i], paste("Trait", i, "lower tail"), 0, 100, TRUE)
      }
      .paweh_qtl_num(values$k, "Lower-to-upper selected ratio", 0, Inf, TRUE)
    }
    if (values$objective == "power") {
      objective_name <- if (values$mv_test == "pillai") "Total sample size" else "Upper-tail selected sample"
      objective_raw <- if (values$mv_test == "pillai") values$N else values$N_case
      .paweh_qtl_num(objective_raw, objective_name, 0, Inf, TRUE)
      if (values$mv_test == "pillai" && objective_raw != floor(objective_raw)) {
        stop("Total sample size must be an integer.", call. = FALSE)
      }
      objective_value <- objective_raw
    }
    args <- list(
      alpha = values$alpha, qtl_var = qtl_var, tau = tau,
      pd = values$pd / 100, cor_matrix = cor_matrix, test = values$mv_test
    )
    if (values$mv_test == "threshold_chisq") {
      args <- c(args, list(x_upper = x_upper, x_lower = x_lower, k = values$k))
    }
  }
  list(
    subtype = values$subtype, objective = values$objective,
    objective_value = objective_value, backend_args = args, display = values
  )
}

.paweh_qtl_function <- function(snapshot) {
  switch(snapshot$subtype,
    continuous = if (snapshot$objective == "power") "qtl_anova_power" else "qtl_anova_mssn",
    extreme = if (snapshot$objective == "power") "qtl_threshold_chisq_power" else "qtl_threshold_chisq_mssn",
    multivariate = if (snapshot$objective == "power") "qtl_multivariate_power_full" else "qtl_multivariate_mssn_full"
  )
}

.paweh_qtl_call_args <- function(snapshot, args = snapshot$backend_args) {
  if (snapshot$objective == "power") {
    if (snapshot$subtype == "continuous" ||
        snapshot$subtype == "multivariate" && args$test == "pillai") args$N <- snapshot$objective_value else args$N_case <- snapshot$objective_value
  } else {
    args$power <- snapshot$objective_value
  }
  args$verbose <- FALSE
  args
}

.paweh_qtl_call <- function(snapshot, args = snapshot$backend_args) {
  do.call(get(.paweh_qtl_function(snapshot), mode = "function"), .paweh_qtl_call_args(snapshot, args))
}

.paweh_qtl_calculate <- function(snapshot) {
  list(snapshot = snapshot, result = .paweh_qtl_call(snapshot))
}

.paweh_qtl_sig <- function(values) serialize(values, NULL)
.paweh_qtl_pct <- function(x, digits = 1) .paweh_format_percent(x, digits)
.paweh_qtl_count <- function(x) .paweh_format_count(x, round_up = TRUE)
.paweh_qtl_number <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
.paweh_qtl_subtype_labels <- c(
  continuous = "Full continuous trait", extreme = "Extreme phenotype sampling",
  multivariate = "Multiple quantitative traits"
)

.paweh_qtl_model_summary <- function(calculation) {
  s <- calculation$snapshot
  v <- s$display
  if (s$subtype == "continuous") paste0(
    "Full continuous trait | allele frequency ", v$pd,
    "% | variance explained ", v$qtl_var, "% | dominance ", v$tau
  ) else if (s$subtype == "extreme") paste0(
    "Extreme phenotype sampling | allele frequency ", v$pd,
    "% | upper/lower tails ", v$x_upper, "%/", v$x_lower, "%"
  ) else paste0(
    "Multiple quantitative traits | ", v$n_traits, " traits | allele frequency ",
    v$pd, "% | ", if (v$mv_test == "pillai") "joint continuous-trait test" else "joint extreme-selection test"
  )
}

.paweh_qtl_result_rows <- function(calculation) {
  s <- calculation$snapshot
  r <- calculation$result
  if (s$objective == "power") {
    rows <- list(.paweh_summary_row("Expected power", .paweh_qtl_pct(r$power, 1)))
    if (s$subtype == "continuous" || s$subtype == "multivariate" && r$test == "pillai") {
      rows <- c(rows, list(.paweh_summary_row("Sample size", .paweh_qtl_count(r$N))))
    } else rows <- c(rows, list(
      .paweh_summary_row("Upper-tail selected", .paweh_qtl_count(r$N_case)),
      .paweh_summary_row("Lower-tail selected", .paweh_qtl_count(r$N_control)),
      .paweh_summary_row("Total selected sample", .paweh_qtl_count(r$N_total))
    ))
  } else {
    if (s$subtype == "continuous" || s$subtype == "multivariate" && r$test == "pillai") {
      rows <- list(.paweh_summary_row("Required sample size", .paweh_qtl_count(r$N)))
    } else rows <- list(
      .paweh_summary_row("Required upper-tail selected", .paweh_qtl_count(r$N_case)),
      .paweh_summary_row("Required lower-tail selected", .paweh_qtl_count(r$N_control)),
      .paweh_summary_row("Required total selected", .paweh_qtl_count(r$N_total))
    )
    rows <- c(rows, list(
      .paweh_summary_row("Target power", .paweh_qtl_pct(r$target_power, 1)),
      if (!is.null(r$achieved_power)) .paweh_summary_row("Achieved power", .paweh_qtl_pct(r$achieved_power, 1))
    ))
  }
  rows
}

.paweh_qtl_test_label <- function(calculation) {
  s <- calculation$snapshot
  if (s$subtype == "continuous") "Quantitative-trait ANOVA" else if (s$subtype == "extreme") {
    "Threshold-selected genotype chi-square"
  } else if (calculation$result$test == "pillai") "Pillai's trace MANOVA" else {
    "Joint threshold-selected genotype chi-square"
  }
}

.paweh_qtl_interpretation <- function(calculation) {
  s <- calculation$snapshot
  r <- calculation$result
  design <- unname(.paweh_qtl_subtype_labels[s$subtype])
  if (s$objective == "power") {
    n_text <- if (s$subtype == "continuous" || s$subtype == "multivariate" && r$test == "pillai") {
      paste(.paweh_qtl_count(r$N), "individuals")
    } else paste(.paweh_qtl_count(r$N_total), "selected individuals")
    paste0("Under the specified ", tolower(design), " assumptions, a sample of ",
      n_text, " has an estimated power of ", .paweh_qtl_pct(r$power, 1),
      " using ", .paweh_qtl_test_label(calculation), ".")
  } else {
    n_text <- if (s$subtype == "continuous" || s$subtype == "multivariate" && r$test == "pillai") {
      paste(.paweh_qtl_count(r$N), "individuals")
    } else paste(.paweh_qtl_count(r$N_total), "selected individuals")
    paste0("Under the specified ", tolower(design), " assumptions, approximately ",
      n_text, " are required to achieve ", .paweh_qtl_pct(r$target_power, 0),
      " power using ", .paweh_qtl_test_label(calculation), ".")
  }
}

.paweh_qtl_repro_args <- function(calculation) {
  .paweh_qtl_call_args(calculation$snapshot)
}
.paweh_qtl_repro_call <- function(calculation) {
  .paweh_repro_call(.paweh_qtl_function(calculation$snapshot), .paweh_qtl_repro_args(calculation))
}

.paweh_qtl_named_rows <- function(values, prefix, digits = 4) {
  lapply(seq_along(values), function(i) {
    .paweh_summary_row(paste(prefix, i), .paweh_qtl_number(values[[i]], digits))
  })
}
.paweh_qtl_genotype_rows <- function(values, label, digits = 4) {
  lapply(seq_along(values), function(i) {
    .paweh_summary_row(paste(label, i - 1L), .paweh_qtl_number(values[[i]], digits))
  })
}
.paweh_qtl_matrix_rows <- function(matrix_value, label, digits = 4) {
  rows <- list()
  for (i in seq_len(nrow(matrix_value))) for (j in seq_len(ncol(matrix_value))) {
    rows[[length(rows) + 1L]] <- .paweh_summary_row(
      paste(label, i, j, sep = " "), .paweh_qtl_number(matrix_value[i, j], digits)
    )
  }
  rows
}

.paweh_qtl_advanced_ui <- function(calculation) {
  s <- calculation$snapshot
  r <- calculation$result
  falconer <- r$falconer
  if (s$subtype != "multivariate") {
    model_rows <- list(
      .paweh_summary_row("Modeled-allele frequency", .paweh_qtl_number(falconer$pd)),
      .paweh_summary_row("QTL variance", .paweh_qtl_number(falconer$qtl_var)),
      .paweh_summary_row("Dominance parameter", .paweh_qtl_number(falconer$tau)),
      .paweh_summary_row("Additive effect", .paweh_qtl_number(falconer$a)),
      .paweh_summary_row("Dominance effect", .paweh_qtl_number(falconer$delta)),
      .paweh_summary_row("Residual variance", .paweh_qtl_number(falconer$residual_variance))
    )
    genotype_rows <- c(
      .paweh_qtl_genotype_rows(falconer$pi, "Genotype frequency"),
      .paweh_qtl_genotype_rows(falconer$mu, "Genotype mean")
    )
  } else {
    model_rows <- c(list(
      .paweh_summary_row("Number of traits", r$falconer$number_of_phenotypes),
      .paweh_summary_row("Modeled-allele frequency", .paweh_qtl_number(r$falconer$pd))
    ),
      .paweh_qtl_named_rows(r$falconer$parameters$qtl_var, "Trait QTL variance"),
      .paweh_qtl_named_rows(r$falconer$parameters$tau, "Trait dominance")
    )
    genotype_rows <- c(
      .paweh_qtl_genotype_rows(r$falconer$genotype_frequencies, "Genotype frequency"),
      .paweh_qtl_matrix_rows(r$falconer$mean_matrix, "Mean: trait/genotype")
    )
  }

  if (s$subtype == "continuous") {
    test_rows <- list(
      .paweh_summary_row("Numerator df", r$df1), .paweh_summary_row("Denominator df", r$df2),
      .paweh_summary_row("Non-centrality parameter", .paweh_qtl_number(r$lambda)),
      .paweh_summary_row("Genotype count method", r$count_method)
    )
    selection <- NULL
  } else if (s$subtype == "extreme") {
    test_rows <- list(
      .paweh_summary_row("Degrees of freedom", r$df),
      .paweh_summary_row("Non-centrality parameter", .paweh_qtl_number(if (s$objective == "power") r$lambda else r$lambda_star)),
      .paweh_summary_row("Internal effect component S", .paweh_qtl_number(r$S))
    )
    selection <- .paweh_detail_section("Selection quantities", c(list(
      .paweh_summary_row("Lower standardized threshold", .paweh_qtl_number(r$thresholds[["lower"]])),
      .paweh_summary_row("Upper standardized threshold", .paweh_qtl_number(r$thresholds[["upper"]])),
      .paweh_summary_row("Expected lower-tail proportion", .paweh_qtl_pct(r$prevalences[["unaffected"]], 2)),
      .paweh_summary_row("Expected upper-tail proportion", .paweh_qtl_pct(r$prevalences[["affected"]], 2))
    ),
      .paweh_qtl_genotype_rows(r$frequencies$affected, "Upper-tail genotype frequency"),
      .paweh_qtl_genotype_rows(r$frequencies$unaffected, "Lower-tail genotype frequency"),
      if (s$objective == "mssn") list(
        .paweh_summary_row("Expected population screened for upper tail", .paweh_qtl_count(r$expected_population_screened_cases)),
        .paweh_summary_row("Expected population screened for lower tail", .paweh_qtl_count(r$expected_population_screened_controls))
      )
    ))
  } else {
    if (r$test == "pillai") test_rows <- list(
      .paweh_summary_row("Numerator df", r$numerator_df),
      .paweh_summary_row("Denominator df", .paweh_qtl_number(r$denominator_df, 2)),
      .paweh_summary_row("Non-centrality parameter", .paweh_qtl_number(r$noncentrality_parameter)),
      .paweh_summary_row("Pillai V-star", .paweh_qtl_number(r$pillai$V_star)),
      if (s$objective == "mssn") .paweh_summary_row("Historical fractional MSSN", .paweh_qtl_number(r$historical_fractional_mssn, 6))
    ) else test_rows <- list(
      .paweh_summary_row("Degrees of freedom", r$df),
      .paweh_summary_row("Non-centrality parameter", .paweh_qtl_number(r$noncentrality_parameter)),
      .paweh_summary_row("Cells with expected count below 1", r$cells_below_one),
      .paweh_summary_row("Integration method", paste(unique(vapply(
        c(r$thresholds$integration$affected, r$thresholds$integration$unaffected),
        `[[`, "", "algorithm"
      )), collapse = ", "))
    )
    selection <- .paweh_detail_section("Multivariate model quantities", c(
      .paweh_qtl_matrix_rows(r$falconer$phenotype_correlation_matrix, "Correlation"),
      .paweh_qtl_matrix_rows(r$falconer$residual_covariance_matrix, "Residual covariance"),
      if (r$test == "threshold_chisq") list(
        .paweh_summary_row("Selection definition", r$thresholds$selection_definition),
        .paweh_summary_row("Expected upper-tail proportion", .paweh_qtl_pct(r$thresholds$prevalences[["affected"]], 2)),
        .paweh_summary_row("Expected lower-tail proportion", .paweh_qtl_pct(r$thresholds$prevalences[["unaffected"]], 2))
      )
    ))
  }
  .paweh_advanced_details_ui(
    .paweh_detail_section("Model specification", model_rows),
    .paweh_detail_section("Genotype and model quantities", genotype_rows),
    .paweh_detail_section("Test-specific quantities", test_rows),
    selection,
    .paweh_reproduce_ui(.paweh_qtl_repro_call(calculation))
  )
}

.paweh_qtl_results_ui <- function(calculation) {
  shiny::tagList(
    shiny::div(class = "paweh-model-specification", .paweh_qtl_model_summary(calculation)),
    bslib::card(
      class = "paweh-result-card",
      bslib::card_header(.paweh_qtl_test_label(calculation)),
      bslib::card_body(shiny::div(class = "paweh-summary-grid", .paweh_qtl_result_rows(calculation)))
    ),
    shiny::div(
      class = "paweh-interpretation", shiny::h4("Interpretation"),
      shiny::p(.paweh_qtl_interpretation(calculation))
    ),
    .paweh_qtl_advanced_ui(calculation)
  )
}

.paweh_qtl_specs <- function(calculation) {
  s <- calculation$snapshot
  a <- s$backend_args
  specs <- list()
  if (s$subtype == "multivariate") {
    for (i in seq_along(a$qtl_var)) specs[[paste0("qtl_var_", i)]] <- list(
      label = paste("Trait", i, "variance explained"), min = 0.005, max = 0.5,
      value = a$qtl_var[i], type = "vector", key = "qtl_var", index = i
    )
    for (i in seq_along(a$tau)) specs[[paste0("tau_", i)]] <- list(
      label = paste("Trait", i, "dominance parameter"), min = -1, max = 1,
      value = a$tau[i], type = "vector", key = "tau", index = i
    )
    if (length(a$qtl_var) == 2L) specs$correlation <- list(
      label = "Phenotype correlation", min = -0.9, max = 0.9,
      value = a$cor_matrix[1, 2], type = "correlation"
    )
  } else {
    specs$qtl_var <- list(label = "Variance explained by the QTL", min = 0.005, max = 0.5, value = a$qtl_var, key = "qtl_var")
    specs$tau <- list(label = "Dominance parameter", min = -1, max = 1, value = a$tau, key = "tau")
  }
  specs$pd <- list(label = "Modeled-allele frequency", min = 0.01, max = 0.99, value = a$pd, key = "pd")
  specs$alpha <- list(label = "Significance level", min = 1e-5, max = 0.2, value = a$alpha, key = "alpha")
  if (s$subtype == "extreme") {
    specs$x_upper <- list(
      label = "Upper population-tail percentage", min = 0.5,
      max = min(49, 99 - a$x_lower), value = a$x_upper, key = "x_upper"
    )
    specs$x_lower <- list(
      label = "Lower population-tail percentage", min = 0.5,
      max = min(49, 99 - a$x_upper), value = a$x_lower, key = "x_lower"
    )
  }
  if (s$subtype == "multivariate" && a$test == "threshold_chisq") {
    for (i in seq_along(a$x_upper)) specs[[paste0("x_upper_", i)]] <- list(
      label = paste("Trait", i, "upper-tail percentage"), min = 0.5,
      max = min(49, 99 - a$x_lower[i]), value = a$x_upper[i],
      type = "vector", key = "x_upper", index = i
    )
  }
  objective_key <- if (s$objective == "power") {
    if (s$subtype == "continuous" || s$subtype == "multivariate" && a$test == "pillai") "N" else "N_case"
  } else "power"
  specs[[objective_key]] <- list(
    label = if (s$objective == "mssn") "Target power" else if (objective_key == "N") "Sample size" else "Upper-tail selected sample",
    min = if (s$objective == "mssn") 0.5 else max(4, s$objective_value / 4),
    max = if (s$objective == "mssn") 0.99 else s$objective_value * 2,
    value = s$objective_value, objective = TRUE
  )
  specs
}

.paweh_qtl_sensitivity <- function(calculation, parameter, range, n = 40L) {
  spec <- .paweh_qtl_specs(calculation)[[parameter]]
  if (is.null(spec)) stop("Unavailable sensitivity parameter.", call. = FALSE)
  grid <- seq(range[1], range[2], length.out = n)
  data <- do.call(rbind, lapply(grid, function(value) {
    snapshot <- calculation$snapshot
    args <- snapshot$backend_args
    if (isTRUE(spec$objective)) snapshot$objective_value <- value else if (identical(spec$type, "vector")) {
      args[[spec$key]][spec$index] <- value
    } else if (identical(spec$type, "correlation")) {
      args$cor_matrix[1, 2] <- args$cor_matrix[2, 1] <- value
    } else args[[spec$key]] <- value
    result <- tryCatch(.paweh_qtl_call(snapshot, args), error = function(error) NULL)
    y <- if (is.null(result)) NA_real_ else if (snapshot$objective == "power") result$power else if (
      snapshot$subtype == "continuous" || snapshot$subtype == "multivariate" && result$test == "pillai"
    ) result$N else result$N_total
    data.frame(x = value, y = y)
  }))
  plot_data <- data
  plot_data$y[!is.finite(plot_data$y)] <- NA_real_
  list(
    data = data, plot_data = plot_data, label = spec$label,
    baseline_x = spec$value, objective = calculation$snapshot$objective,
    has_nonfinite = any(!is.finite(data$y))
  )
}

.paweh_qtl_sensitivity_plot <- function(sensitivity) {
  colors <- .paweh_plot_colors()
  ggplot2::ggplot(sensitivity$plot_data, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_line(color = unname(colors["trend"]), linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_vline(
      xintercept = sensitivity$baseline_x, color = unname(colors["reference"]),
      linetype = "dotted", linewidth = 0.6
    ) +
    ggplot2::labs(
      x = sensitivity$label,
      y = if (sensitivity$objective == "power") "Power" else "Required sample size"
    ) + .paweh_plot_theme()
}

.paweh_qtl_single_plot <- function(calculation) {
  a <- calculation$snapshot$backend_args
  plot <- plot_qtl_genotype_distribution(
    qtl_var = a$qtl_var, tau = a$tau, pd = a$pd,
    type = "density", scale = "frequency", show_means = TRUE,
    verbose = FALSE
  )
  if (calculation$snapshot$subtype == "extreme") {
    thresholds <- calculation$result$thresholds
    plot <- plot +
      ggplot2::annotate(
        "rect", xmin = -Inf, xmax = thresholds[["lower"]], ymin = -Inf, ymax = Inf,
        fill = "#355C7D", alpha = 0.08
      ) +
      ggplot2::annotate(
        "rect", xmin = thresholds[["upper"]], xmax = Inf, ymin = -Inf, ymax = Inf,
        fill = "#A8844F", alpha = 0.08
      ) +
      ggplot2::geom_vline(
        xintercept = unname(thresholds), color = "#7A848C", linetype = "dotted"
      ) +
      ggplot2::labs(subtitle = "Lower selected tail | unselected middle | upper selected tail")
  }
  plot
}

.paweh_qtl_trait_pair <- function(calculation, x = NULL, y = NULL) {
  p <- length(calculation$snapshot$backend_args$qtl_var)
  if (p < 2L) stop("Multivariate visualization requires at least two traits.", call. = FALSE)
  pair <- suppressWarnings(as.integer(c(x, y)))
  if (length(pair) != 2L || anyNA(pair) || any(!pair %in% seq_len(p)) || pair[1] == pair[2]) {
    pair <- 1:2
  }
  pair
}

.paweh_qtl_multivariate_plot_args <- function(calculation, surface, grid_n = 70L, trait_pair = NULL) {
  a <- calculation$snapshot$backend_args
  if (is.null(trait_pair)) trait_pair <- .paweh_qtl_trait_pair(calculation)
  trait_pair <- .paweh_qtl_trait_pair(calculation, trait_pair[1], trait_pair[2])
  args <- list(
    qtl_var = a$qtl_var[trait_pair], tau = a$tau[trait_pair], pd = a$pd,
    cor_matrix = a$cor_matrix[trait_pair, trait_pair, drop = FALSE],
    surface = surface, grid_n = grid_n
  )
  if (a$test == "threshold_chisq" || surface == "cdf") {
    if (!is.null(a$x_upper)) args <- c(args, list(
      x_upper = a$x_upper[trait_pair], x_lower = a$x_lower[trait_pair]
    ))
  }
  attr(args, "trait_pair") <- trait_pair
  args
}

.paweh_qtl_multivariate_plot <- function(calculation, mode, trait_pair = NULL) {
  surface <- switch(mode,
    genotype = "genotype_density", mixture = "density", selection = "cdf"
  )
  args <- .paweh_qtl_multivariate_plot_args(calculation, surface, trait_pair = trait_pair)
  pair <- attr(args, "trait_pair")
  attr(args, "trait_pair") <- NULL
  do.call(plot_qtl_multivariate_contour, args) + ggplot2::labs(
    x = paste("Trait", pair[1], "value"), y = paste("Trait", pair[2], "value")
  )
}

.paweh_qtl_multivariate_surface <- function(calculation, mode, grid_n = 30L, trait_pair = NULL) {
  surface <- if (identical(mode, "genotype")) "genotype_density" else "density"
  args <- .paweh_qtl_multivariate_plot_args(calculation, surface, grid_n, trait_pair)
  pair <- attr(args, "trait_pair")
  attr(args, "trait_pair") <- NULL
  args$z_scale <- "raw"
  plotly::layout(
    do.call(plot_qtl_multivariate_surface3d, args),
    scene = list(
      xaxis = list(title = paste("Trait", pair[1], "value")),
      yaxis = list(title = paste("Trait", pair[2], "value"))
    )
  )
}

.paweh_qtl_summary_ui <- function(calculation) {
  s <- calculation$snapshot
  v <- s$display
  rows <- list(
    .paweh_summary_row("Workflow", unname(.paweh_qtl_subtype_labels[s$subtype])),
    .paweh_summary_row("Objective", if (s$objective == "power") "Estimate power" else paste("Minimum sample size at", .paweh_qtl_pct(s$objective_value))),
    .paweh_summary_row("Alpha", format(v$alpha)),
    .paweh_summary_row("Allele frequency", paste0(v$pd, "%"))
  )
  if (s$subtype != "multivariate") rows <- c(rows, list(
    .paweh_summary_row("Variance explained", paste0(v$qtl_var, "%")),
    .paweh_summary_row("Dominance", format(v$tau))
  )) else rows <- c(rows, list(
    .paweh_summary_row("Traits", v$n_traits),
    .paweh_summary_row("Joint test", if (v$mv_test == "pillai") "Continuous-trait MANOVA" else "Extreme-selection chi-square")
  ))
  shiny::div(class = "paweh-summary-grid", rows)
}

.paweh_qtl_methods_ui <- function(calculation = NULL) {
  if (is.null(calculation)) return(.paweh_empty_ui("Methods"))
  s <- calculation$snapshot
  input_specification <- if (s$subtype == "continuous") {
    "Continuous phenotype"
  } else if (s$subtype == "extreme") {
    "Upper and lower population tails"
  } else if (s$backend_args$test == "pillai") {
    "Correlated continuous phenotypes"
  } else {
    "Joint upper and lower population tails"
  }
  shiny::tagList(
    shiny::h3("Methods"),
    shiny::div(
      class = "paweh-summary-grid",
      .paweh_summary_row("Study design", unname(.paweh_qtl_subtype_labels[s$subtype])),
      .paweh_summary_row("Objective", if (s$objective == "power") "Power" else "Minimum sample size"),
      .paweh_summary_row("Input specification", input_specification),
      .paweh_summary_row("Statistical method", .paweh_qtl_test_label(calculation)),
      if (s$subtype == "multivariate") .paweh_summary_row("Traits", length(s$backend_args$qtl_var)),
      .paweh_summary_row("Canonical function", paste0(.paweh_qtl_function(s), "()"))
    ),
    shiny::p(
      "The canonical Falconer model supplies genotype-specific phenotype quantities. Sensitivity points are separate canonical calls using the frozen design."
    ),
    shiny::a(
      href = "https://akilanthony.github.io/pawh/articles/paweh-04-quantitative-trait-study-design.html",
      target = "_blank", rel = "noopener", "Read the Quantitative Trait vignette"
    )
  )
}

.paweh_qtl_trait_inputs <- function(ns, p, threshold = FALSE) {
  controls <- list()
  defaults <- .paweh_qtl_defaults()
  for (i in seq_len(p)) controls <- c(controls, list(
    shiny::h6(paste("Trait", i)),
    shiny::numericInput(
      ns(paste0("mv_qtl_var_", i)), "Variance explained (%)",
      defaults[[paste0("mv_qtl_var_", i)]], 0.01, 99.99, 0.1
    ),
    shiny::numericInput(
      ns(paste0("mv_tau_", i)), "Dominance parameter",
      defaults[[paste0("mv_tau_", i)]], NA, NA, 0.1
    ),
    if (threshold) shiny::numericInput(
      ns(paste0("mv_x_upper_", i)), "Upper population tail selected (%)",
      defaults[[paste0("mv_x_upper_", i)]], 0.01, 99.99, 0.5
    ),
    if (threshold) shiny::numericInput(
      ns(paste0("mv_x_lower_", i)), "Lower population tail selected (%)",
      defaults[[paste0("mv_x_lower_", i)]], 0.01, 99.99, 0.5
    )
  ))
  shiny::tagList(controls)
}

.paweh_qtl_correlation_inputs <- function(ns, p) {
  controls <- list()
  defaults <- .paweh_qtl_defaults()
  for (i in seq_len(p - 1L)) for (j in seq.int(i + 1L, p)) {
    id <- paste0("corr_", i, "_", j)
    controls[[length(controls) + 1L]] <- shiny::numericInput(
      ns(id), paste("Correlation: Trait", i, "and Trait", j), defaults[[id]], -0.99, 0.99, 0.05
    )
  }
  shiny::tagList(controls)
}

.paweh_qtl_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    state <- shiny::reactiveValues(
      calculation = NULL, error = NULL, signature = NULL,
      sensitivity = NULL, surface = NULL, surface_error = NULL
    )
    values <- shiny::reactive(.paweh_qtl_values(input))
    changed <- shiny::reactive(
      !is.null(state$signature) && !identical(state$signature, .paweh_qtl_sig(values()))
    )

    output$subtype_help <- shiny::renderUI({
      text <- switch(if (is.null(input$subtype)) "continuous" else input$subtype,
        continuous = "Analyze the measured quantitative phenotype directly.",
        extreme = "Select individuals from the upper and lower phenotype tails.",
        multivariate = "Analyze correlated quantitative phenotypes jointly."
      )
      shiny::p(class = "text-muted", text)
    })
    output$objective_inputs <- shiny::renderUI({
      objective <- if (is.null(input$objective)) "power" else input$objective
      subtype <- if (is.null(input$subtype)) "continuous" else input$subtype
      test <- if (is.null(input$mv_test)) "pillai" else input$mv_test
      if (objective == "mssn") {
        shiny::sliderInput(ns("target_power"), "Target power (%)", 50, 99, 80)
      } else if (subtype == "continuous" || subtype == "multivariate" && test == "pillai") {
        shiny::numericInput(ns("N"), "Total sample size", 500, 4, NA, 10)
      } else {
        shiny::numericInput(ns("N_case"), "Upper-tail selected sample", 150, 1, NA, 5)
      }
    })
    output$core_inputs <- shiny::renderUI({
      subtype <- if (is.null(input$subtype)) "continuous" else input$subtype
      common <- shiny::tagList(
        shiny::numericInput(ns("pd"), "Modeled-allele frequency (%)", 30, 0.01, 99.99, 0.1),
        shiny::tags$small(class = "text-muted", "Population frequency of the allele defining the three genotype groups.")
      )
      if (subtype == "continuous") shiny::tagList(
        common,
        shiny::numericInput(ns("qtl_var"), "Variance explained by the QTL (%)", 10, 0.01, 99.99, 0.1),
        shiny::tags$small(class = "text-muted", "Proportion of total trait variance attributable to the modeled QTL."),
        shiny::numericInput(ns("tau"), "Dominance parameter", 0, NA, NA, 0.1),
        shiny::tags$small(class = "text-muted", "Dominance-to-additivity ratio in the canonical Falconer model.")
      ) else if (subtype == "extreme") shiny::tagList(
        common,
        shiny::numericInput(ns("qtl_var"), "Variance explained by the QTL (%)", 10, 0.01, 99.99, 0.1),
        shiny::numericInput(ns("tau"), "Dominance parameter", 0, NA, NA, 0.1),
        shiny::numericInput(ns("x_upper"), "Upper population tail selected (%)", 10, 0.01, 99.99, 0.5),
        shiny::numericInput(ns("x_lower"), "Lower population tail selected (%)", 10, 0.01, 99.99, 0.5),
        shiny::tags$small(class = "text-muted", "Tail percentages define standardized-normal population percentiles; the middle is excluded.")
      ) else {
        p <- if (is.null(input$n_traits)) 2L else as.integer(input$n_traits)
        test <- if (is.null(input$mv_test)) "pillai" else input$mv_test
        shiny::tagList(
          common,
          shiny::selectInput(ns("n_traits"), "Number of quantitative traits", stats::setNames(2:4, 2:4)),
          shiny::radioButtons(ns("mv_test"), "Joint analysis", c(
            `Joint continuous-trait test (Pillai MANOVA)` = "pillai",
            `Joint extreme-selection test` = "threshold_chisq"
          )),
          .paweh_qtl_trait_inputs(ns, p, threshold = test == "threshold_chisq"),
          shiny::h6("Phenotype correlations"),
          .paweh_qtl_correlation_inputs(ns, p),
          shiny::tags$small(class = "text-muted", "Pairwise phenotype correlations; the assembled matrix must be positive definite."),
          if (p > 2L) shiny::tags$small(
            class = "text-muted",
            "Statistical calculations support up to four traits; contour and 3D views display a selected trait pair."
          )
        )
      }
    })
    output$advanced_inputs <- shiny::renderUI({
      subtype <- if (is.null(input$subtype)) "continuous" else input$subtype
      if (subtype == "continuous") shiny::tagList(
        shiny::selectInput(ns("count_method"), "Genotype-count method", c(Rounded = "rounded", Expected = "expected")),
        if (identical(input$objective, "mssn")) shiny::checkboxInput(ns("multiple_of_three"), "Restrict sample size to multiples of three", TRUE)
      ) else shiny::numericInput(ns("k"), "Lower-tail selected per upper-tail selected", 1, 0.01, NA, 0.1)
    })

    shiny::observeEvent(input$calculate, {
      state$error <- NULL; state$sensitivity <- NULL
      state$surface <- NULL; state$surface_error <- NULL
      submitted <- values()
      tryCatch({
        state$calculation <- .paweh_qtl_calculate(.paweh_qtl_snapshot(submitted))
        state$signature <- .paweh_qtl_sig(submitted)
      }, error = function(error) {
        state$error <- paste("This design could not be calculated.", conditionMessage(error))
      })
    }, ignoreInit = TRUE)

    output$changed_notice <- shiny::renderUI(if (changed()) shiny::div(
      class = "paweh-changed-notice", role = "status",
      "Inputs have changed. Recalculate to update results."
    ))
    output$design_summary <- shiny::renderUI(if (!is.null(state$error)) {
      shiny::div(class = "paweh-error", role = "alert", state$error)
    } else if (is.null(state$calculation)) {
      shiny::p(class = "text-muted", "Choose assumptions and select Calculate.")
    } else .paweh_qtl_summary_ui(state$calculation))
    output$results <- shiny::renderUI(if (!is.null(state$error)) {
      shiny::div(class = "paweh-error", role = "alert", state$error)
    } else if (is.null(state$calculation)) {
      .paweh_empty_ui("Results")
    } else .paweh_qtl_results_ui(state$calculation))

    output$sensitivity_controls <- shiny::renderUI({
      if (is.null(state$calculation)) return(.paweh_empty_ui("Sensitivity"))
      specs <- .paweh_qtl_specs(state$calculation)
      choices <- stats::setNames(names(specs), vapply(specs, `[[`, "", "label"))
      spec <- specs[[1L]]
      shiny::div(
        class = "paweh-sensitivity-controls",
        shiny::selectInput(ns("sensitivity_parameter"), "Parameter", choices),
        shiny::sliderInput(ns("sensitivity_range"), "Range", spec$min, spec$max,
          pmax(spec$min, pmin(spec$max, c(spec$value * 0.75, spec$value * 1.25)))),
        shiny::actionButton(ns("run_sensitivity"), "Run sensitivity analysis")
      )
    })
    shiny::observeEvent(input$sensitivity_parameter, {
      shiny::req(state$calculation)
      spec <- .paweh_qtl_specs(state$calculation)[[input$sensitivity_parameter]]
      shiny::req(spec)
      value <- pmax(spec$min, pmin(spec$max, c(spec$value * 0.75, spec$value * 1.25)))
      shiny::updateSliderInput(session, "sensitivity_range", min = spec$min, max = spec$max,
        value = value, step = (spec$max - spec$min) / 100)
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$run_sensitivity, {
      shiny::req(state$calculation, input$sensitivity_parameter, input$sensitivity_range)
      state$sensitivity <- .paweh_qtl_sensitivity(
        state$calculation, input$sensitivity_parameter, input$sensitivity_range
      )
    }, ignoreInit = TRUE)
    output$sensitivity_message <- shiny::renderUI(if (!is.null(state$calculation) && is.null(state$sensitivity)) {
      shiny::p(class = "text-muted", "Each point is a canonical calculation of the frozen design.")
    } else if (!is.null(state$sensitivity) && state$sensitivity$has_nonfinite) {
      shiny::div(class = "paweh-caution", "Some explored designs are invalid or non-finite and are omitted from the line geometry.")
    } else if (.paweh_power_axis_zoomed(state$sensitivity)) {
      shiny::p(class = "paweh-zoom-note", "Y-axis is zoomed to show variation in power.")
    })
    output$sensitivity_plot_container <- shiny::renderUI(if (!is.null(state$sensitivity)) {
      shiny::plotOutput(ns("sensitivity_plot"), height = "430px")
    })
    output$sensitivity_plot <- shiny::renderPlot({
      shiny::req(state$sensitivity); .paweh_qtl_sensitivity_plot(state$sensitivity)
    })

    output$visualize_intro <- shiny::renderUI(if (is.null(state$calculation)) {
      .paweh_empty_ui("Visualize")
    } else if (state$calculation$snapshot$subtype == "continuous") shiny::div(
      class = "paweh-visual-intro", shiny::h3("Genotype distributions"),
      shiny::p("Population-weighted phenotype distributions returned by the canonical Falconer plotting pathway.")
    ) else if (state$calculation$snapshot$subtype == "extreme") shiny::div(
      class = "paweh-visual-intro", shiny::h3("Threshold selection"),
      shiny::p("Canonical genotype distributions with the returned lower and upper selection thresholds.")
    ) else shiny::div(
      class = "paweh-visual-intro", shiny::h3("Multiple-trait distributions"),
      shiny::p("Statistical calculations support up to four traits. Two-dimensional and 3D visualizations display a selected pair from the frozen design.")
    ))
    output$visualization_controls <- shiny::renderUI({
      if (is.null(state$calculation) || state$calculation$snapshot$subtype != "multivariate") return(NULL)
      p <- length(state$calculation$snapshot$backend_args$qtl_var)
      trait_choices <- stats::setNames(as.character(seq_len(p)), paste("Trait", seq_len(p)))
      x <- if (!is.null(input$visual_trait_x) && input$visual_trait_x %in% unname(trait_choices)) input$visual_trait_x else "1"
      y_choices <- trait_choices[unname(trait_choices) != x]
      y <- if (!is.null(input$visual_trait_y) && input$visual_trait_y %in% unname(y_choices)) input$visual_trait_y else unname(y_choices)[1]
      choices <- c(`Genotype distributions` = "genotype", `Overall population` = "mixture")
      if (state$calculation$snapshot$backend_args$test == "threshold_chisq") choices <- c(choices, `Selection regions` = "selection")
      shiny::tagList(
        if (p > 2L) shiny::div(
          class = "paweh-sensitivity-controls",
          shiny::selectInput(ns("visual_trait_x"), "Horizontal trait", trait_choices, selected = x),
          shiny::selectInput(ns("visual_trait_y"), "Vertical trait", y_choices, selected = y)
        ),
        shiny::radioButtons(ns("visual_mode"), "Visualization", choices, inline = TRUE)
      )
    })
    output$visualization_plot_container <- shiny::renderUI(if (!is.null(state$calculation)) {
      shiny::plotOutput(ns("visualization_plot"), height = "460px")
    })
    output$visualization_plot <- shiny::renderPlot({
      shiny::req(state$calculation)
      if (state$calculation$snapshot$subtype == "multivariate") {
        .paweh_qtl_multivariate_plot(
          state$calculation,
          if (is.null(input$visual_mode)) "genotype" else input$visual_mode,
          .paweh_qtl_trait_pair(state$calculation, input$visual_trait_x, input$visual_trait_y)
        )
      } else .paweh_qtl_single_plot(state$calculation)
    })

    output$surface_controls <- shiny::renderUI({
      if (is.null(state$calculation)) return(shiny::p(class = "text-muted", "Calculate a design first."))
      if (state$calculation$snapshot$subtype != "multivariate") {
        return(shiny::p(class = "text-muted", "On-demand 3D visualization is available for multivariate designs."))
      }
      if (!requireNamespace("plotly", quietly = TRUE)) return(shiny::p(class = "text-muted", "Install the suggested plotly package to generate this optional visualization."))
      shiny::tagList(
        shiny::radioButtons(ns("surface_mode"), "3D view", c(
          `Genotype distributions` = "genotype", `Overall population` = "mixture"
        ), inline = TRUE),
        shiny::p(class = "text-muted", if (identical(input$surface_mode, "mixture")) {
          "The overall population surface is the genotype-frequency-weighted mixture. Three genotype groups do not necessarily produce three distinct mixture peaks."
        } else "Each surface represents the phenotype distribution conditional on genotype."),
        shiny::actionButton(ns("generate_surface"), "Generate 3D visualization")
      )
    })
    shiny::observeEvent(input$generate_surface, {
      shiny::req(state$calculation)
      state$surface_error <- NULL
      tryCatch({
        state$surface <- shiny::withProgress(
          message = "Generating QTL visualization", value = 0.5,
          expr = .paweh_qtl_multivariate_surface(
            state$calculation,
            if (is.null(input$surface_mode)) "genotype" else input$surface_mode,
            trait_pair = .paweh_qtl_trait_pair(state$calculation, input$visual_trait_x, input$visual_trait_y)
          )
        )
      }, error = function(error) state$surface_error <- paste("The 3D visualization could not be generated.", conditionMessage(error)))
    }, ignoreInit = TRUE)
    shiny::observeEvent(list(input$visual_trait_x, input$visual_trait_y), {
      state$surface <- NULL
      state$surface_error <- NULL
    }, ignoreInit = TRUE)
    output$surface_message <- shiny::renderUI(if (!is.null(state$surface_error)) {
      shiny::div(class = "paweh-error", role = "alert", state$surface_error)
    } else if (!is.null(state$calculation) && is.null(state$surface)) {
      shiny::p(class = "text-muted", "The Plotly visualization is generated only on request from the frozen design.")
    })
    output$surface_container <- shiny::renderUI(if (!is.null(state$surface)) {
      plotly::plotlyOutput(ns("surface_plot"), height = "500px")
    })
    if (requireNamespace("plotly", quietly = TRUE)) output$surface_plot <- plotly::renderPlotly({
      shiny::req(state$surface); state$surface
    })
    output$methods <- shiny::renderUI(.paweh_qtl_methods_ui(state$calculation))

    list(
      calculation = shiny::reactive(state$calculation), changed = changed,
      error = shiny::reactive(state$error), sensitivity = shiny::reactive(state$sensitivity),
      surface = shiny::reactive(state$surface), surface_error = shiny::reactive(state$surface_error)
    )
  })
}
