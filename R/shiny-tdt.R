# TDT / Family dashboard; all statistical work delegates to tdt_power()/tdt_mssn().

.pawh_tdt_defaults <- function() {
  list(
    objective = "power", input_mode = "model_based", N = 600,
    target_power = 80, alpha = 0.05, prev = 5, pd = 30,
    R1 = 1.5, R2 = 2.25, delta_prime = 1,
    ET = 140, ENT = 100, n_trios = 120,
    misclassification = FALSE, misclass_rate = 0,
    heterogeneity = FALSE, heter_rate = 0
  )
}

.pawh_tdt_values <- function(input) {
  shiny::reactiveValuesToList(input)
  defaults <- .pawh_tdt_defaults()
  values <- lapply(names(defaults), function(name) {
    if (is.null(input[[name]])) defaults[[name]] else input[[name]]
  })
  stats::setNames(values, names(defaults))
}

.pawh_tdt_num <- function(x, label, lower, upper, open = FALSE) {
  invalid <- !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    x < lower || x > upper || (open && (x == lower || x == upper))
  if (invalid) stop(label, " is outside its allowed range.", call. = FALSE)
}

.pawh_tdt_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    class = "pawh-workspace",
    .pawh_page_heading(
      "TDT / Family study design",
      "Estimate power or required affected-child trios for a transmission disequilibrium study."
    ),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = "Design setup", open = "desktop", width = 280,
        resizable = FALSE,
        shiny::radioButtons(ns("objective"), "Objective", c(
          `Estimate power` = "power", `Minimum sample size` = "mssn"
        )),
        shiny::uiOutput(ns("objective_inputs")),
        shiny::radioButtons(ns("input_mode"), "Input specification", c(
          `Genetic model` = "model_based",
          `Direct transmission quantities` = "model_free"
        )),
        shiny::uiOutput(ns("model_inputs")),
        shiny::numericInput(
          ns("alpha"), "Significance level (alpha)", 0.05,
          min = 1e-06, max = 0.999999
        ),
        shiny::tags$details(
          class = "pawh-sidebar-section",
          shiny::tags$summary(shiny::strong("Advanced assumptions")),
          shiny::checkboxInput(ns("misclassification"), "Phenotype misclassification"),
          shiny::conditionalPanel(
            sprintf("input['%s']", ns("misclassification")),
            shiny::numericInput(
              ns("misclass_rate"),
              "Control-to-affected misclassification rate (%)",
              0, min = 0, max = 99.9, step = 0.1
            ),
            shiny::tags$small(
              class = "text-muted",
              "The canonical model's pi01 rate: controls incorrectly classified as affected."
            )
          ),
          shiny::checkboxInput(ns("heterogeneity"), "Locus heterogeneity"),
          shiny::conditionalPanel(
            sprintf("input['%s']", ns("heterogeneity")),
            shiny::numericInput(
              ns("heter_rate"), "Locus heterogeneity rate (%)",
              0, min = 0, max = 99.9, step = 0.1
            ),
            shiny::tags$small(
              class = "text-muted",
              "Proportion of affected trios whose disease is not attributable to the modeled locus."
            )
          )
        ),
        shiny::actionButton(
          ns("calculate"), "Calculate study design",
          class = "btn-primary pawh-calculate"
        ),
        shiny::uiOutput(ns("changed_notice")),
        bslib::card(
          class = "pawh-design-summary",
          bslib::card_header("Your calculated design"),
          bslib::card_body(shiny::uiOutput(ns("design_summary")))
        )
      ),
      bslib::navset_card_tab(
        id = ns("section"),
        bslib::nav_panel("Results", shiny::uiOutput(ns("results"))),
        bslib::nav_panel(
          "Sensitivity", shiny::uiOutput(ns("sensitivity_controls")),
          shiny::uiOutput(ns("sensitivity_message")),
          shiny::plotOutput(ns("sensitivity_plot"), height = "430px")
        ),
        bslib::nav_panel(
          "Visualize", shiny::uiOutput(ns("visualize_intro")),
          shiny::plotOutput(ns("transmission_plot"), height = "430px")
        ),
        bslib::nav_panel("Methods", shiny::uiOutput(ns("methods")))
      )
    )
  )
}

.pawh_tdt_snapshot <- function(values) {
  if (!values$objective %in% c("power", "mssn")) stop("Choose a valid objective.", call. = FALSE)
  if (!values$input_mode %in% c("model_based", "model_free")) stop("Choose a valid input specification.", call. = FALSE)
  .pawh_tdt_num(values$alpha, "Significance level", 0, 1, TRUE)
  if (values$objective == "power") {
    .pawh_tdt_num(values$N, "Number of affected-child trios", 0, Inf, TRUE)
    objective_value <- values$N
  } else {
    .pawh_tdt_num(values$target_power, "Target power", 0, 100, TRUE)
    objective_value <- values$target_power / 100
  }

  args <- list(alpha = values$alpha, input_mode = values$input_mode)
  if (values$input_mode == "model_based") {
    .pawh_tdt_num(values$prev, "Disease prevalence", 0, 100, TRUE)
    .pawh_tdt_num(values$pd, "Modeled-allele frequency", 0, 100, TRUE)
    .pawh_tdt_num(values$R1, "Heterozygote relative risk", 0, Inf, TRUE)
    .pawh_tdt_num(values$R2, "Homozygote relative risk", 0, Inf, TRUE)
    .pawh_tdt_num(values$delta_prime, "D-prime", 0, 1)
    args <- c(args, list(
      prev = values$prev / 100, pd = values$pd / 100,
      R1 = values$R1, R2 = values$R2,
      delta_prime = values$delta_prime
    ))
  } else {
    .pawh_tdt_num(values$ET, "Expected transmitted count", 0, Inf)
    .pawh_tdt_num(values$ENT, "Expected non-transmitted count", 0, Inf)
    if (values$ET + values$ENT <= 0) stop("Expected transmitted and non-transmitted counts cannot both be zero.", call. = FALSE)
    args <- c(args, list(ET = values$ET, ENT = values$ENT))
    if (values$objective == "mssn") {
      .pawh_tdt_num(values$n_trios, "Represented affected-child trios", 0, Inf, TRUE)
      args$n_trios <- values$n_trios
    }
    if (isTRUE(values$misclassification) || isTRUE(values$heterogeneity)) {
      .pawh_tdt_num(values$pd, "Modeled-allele frequency", 0, 100, TRUE)
      args$pd <- values$pd / 100
    }
    if (isTRUE(values$misclassification)) {
      .pawh_tdt_num(values$prev, "Disease prevalence", 0, 100, TRUE)
      args$prev <- values$prev / 100
    }
  }

  .pawh_tdt_num(values$misclass_rate, "Phenotype misclassification rate", 0, 100)
  .pawh_tdt_num(values$heter_rate, "Locus heterogeneity rate", 0, 100)
  misclass_rate <- if (isTRUE(values$misclassification)) values$misclass_rate / 100 else 0
  heter_rate <- if (isTRUE(values$heterogeneity)) values$heter_rate / 100 else 0
  if (misclass_rate >= 1 || heter_rate >= 1) stop("Misclassification and heterogeneity rates must be below 100%.", call. = FALSE)
  args <- c(args, list(misclass_rate = misclass_rate, heter_rate = heter_rate))

  list(
    objective = values$objective, input_mode = values$input_mode,
    objective_value = objective_value, backend_args = args, display = values
  )
}

.pawh_tdt_call <- function(snapshot, args = snapshot$backend_args) {
  args$verbose <- FALSE
  if (snapshot$objective == "power") {
    args$N <- snapshot$objective_value
    do.call(tdt_power, args)
  } else {
    args$target_power <- snapshot$objective_value
    do.call(tdt_mssn, args)
  }
}

.pawh_tdt_calculate <- function(snapshot) {
  result <- .pawh_tdt_call(snapshot)
  list(
    snapshot = snapshot, result = result,
    active = list(
      misclassification = snapshot$backend_args$misclass_rate > 0,
      heterogeneity = snapshot$backend_args$heter_rate > 0
    )
  )
}

.pawh_tdt_sig <- function(values) serialize(values, NULL)
.pawh_tdt_pct <- function(x, digits = 1) {
  if (!is.finite(x)) return("not defined")
  paste0(formatC(100 * x, format = "f", digits = digits), "%")
}
.pawh_tdt_rate <- function(x, digits = 1) paste0(formatC(100 * x, format = "f", digits = digits), "%")
.pawh_tdt_count <- function(x, plan = FALSE) {
  if (is.infinite(x)) return("Inf")
  if (!is.finite(x)) return("not defined")
  value <- if (plan) ceiling(x) else x
  formatC(value, format = if (plan) "d" else "fg", digits = if (plan) 0 else 7, big.mark = ",")
}
.pawh_tdt_scenario_labels <- c(
  no_error = "No-error design",
  misclassification = "Phenotype misclassification",
  heterogeneity = "Locus heterogeneity"
)
.pawh_tdt_scenarios <- function(calculation) {
  c(
    "no_error",
    if (calculation$active$misclassification) "misclassification",
    if (calculation$active$heterogeneity) "heterogeneity"
  )
}

.pawh_tdt_result_card <- function(calculation, scenario) {
  snapshot <- calculation$snapshot
  result <- calculation$result
  rate_row <- if (scenario == "misclassification") {
    .pawh_summary_row("Misclassification rate", .pawh_tdt_rate(snapshot$backend_args$misclass_rate))
  } else if (scenario == "heterogeneity") {
    .pawh_summary_row("Heterogeneity rate", .pawh_tdt_rate(snapshot$backend_args$heter_rate))
  }
  if (snapshot$objective == "power") {
    body <- shiny::tagList(
      rate_row,
      .pawh_summary_row("Expected power", .pawh_tdt_pct(result$power[[scenario]], 1)),
      .pawh_summary_row("Affected-child trios", .pawh_tdt_count(result$N, plan = TRUE)),
      if (scenario != "no_error") .pawh_summary_row(
        "Absolute power loss",
        paste0(formatC(100 * result$power_loss[[scenario]], format = "f", digits = 1), " percentage points")
      )
    )
  } else {
    baseline <- result$N$no_error
    percent <- if (scenario == "no_error") NULL else result$percent_increase[[scenario]]
    body <- shiny::tagList(
      rate_row,
      .pawh_summary_row("Required affected-child trios", .pawh_tdt_count(result$N[[scenario]], plan = TRUE)),
      .pawh_summary_row("Target power", .pawh_tdt_pct(result$target_power, 1)),
      if (scenario != "no_error") .pawh_summary_row(
        "Increase", if (is.finite(result$N[[scenario]] - baseline)) {
          paste0("+", .pawh_tdt_count(ceiling(result$N[[scenario]]) - ceiling(baseline), plan = TRUE), " trios")
        } else "not defined"
      ),
      if (scenario != "no_error") .pawh_summary_row(
        "Percent increase",
        if (is.finite(percent)) paste0(formatC(percent, format = "f", digits = 1), "%") else "not defined"
      )
    )
  }
  bslib::card(
    class = "pawh-result-card",
    bslib::card_header(unname(.pawh_tdt_scenario_labels[scenario])),
    bslib::card_body(shiny::div(class = "pawh-summary-grid", body))
  )
}

.pawh_tdt_interpretation <- function(calculation) {
  snapshot <- calculation$snapshot
  result <- calculation$result
  scenarios <- .pawh_tdt_scenarios(calculation)
  if (snapshot$objective == "power") {
    text <- paste0(
      "Under the specified assumptions, a study with ", .pawh_tdt_count(result$N, TRUE),
      " affected-child trios has an estimated no-error power of ",
      .pawh_tdt_pct(result$power$no_error, 1), "."
    )
    additions <- vapply(setdiff(scenarios, "no_error"), function(scenario) {
      paste0(" The ", tolower(.pawh_tdt_scenario_labels[[scenario]]),
        " scenario has estimated power of ", .pawh_tdt_pct(result$power[[scenario]], 1), ".")
    }, "")
  } else {
    text <- paste0(
      "Under the no-error model, approximately ", .pawh_tdt_count(result$N$no_error, TRUE),
      " affected-child trios are required to achieve ", .pawh_tdt_pct(result$target_power, 0), " power."
    )
    additions <- vapply(setdiff(scenarios, "no_error"), function(scenario) {
      rate <- snapshot$backend_args[[if (scenario == "misclassification") "misclass_rate" else "heter_rate"]]
      paste0(" At a ", .pawh_tdt_rate(rate, 1), " ",
        if (scenario == "misclassification") "phenotype misclassification" else "locus heterogeneity",
        " rate, the required sample size is ", .pawh_tdt_count(result$N[[scenario]], TRUE), " trios.")
    }, "")
  }
  paste0(text, paste0(additions, collapse = ""))
}

.pawh_tdt_results_ui <- function(calculation) {
  scenarios <- .pawh_tdt_scenarios(calculation)
  shiny::tagList(
    bslib::layout_column_wrap(
      width = "300px",
      lapply(scenarios, function(scenario) .pawh_tdt_result_card(calculation, scenario))
    ),
    if (all(unlist(calculation$active))) shiny::div(
      class = "pawh-caution",
      "Phenotype misclassification and locus heterogeneity are evaluated as separate sensitivity scenarios in the current TDT model."
    ),
    shiny::div(
      class = "pawh-interpretation", shiny::h4("Interpretation"),
      shiny::p(.pawh_tdt_interpretation(calculation))
    )
  )
}

.pawh_tdt_specs <- function(calculation) {
  snapshot <- calculation$snapshot
  values <- snapshot$display
  specs <- list()
  if (calculation$active$misclassification) specs$misclass_rate <- list(
    label = "Phenotype misclassification rate", min = 0,
    max = min(0.5, max(0.1, snapshot$backend_args$misclass_rate * 2)),
    value = snapshot$backend_args$misclass_rate, key = "misclass_rate"
  )
  if (calculation$active$heterogeneity) specs$heter_rate <- list(
    label = "Locus heterogeneity rate", min = 0,
    max = min(0.9, max(0.25, snapshot$backend_args$heter_rate * 2)),
    value = snapshot$backend_args$heter_rate, key = "heter_rate"
  )
  if (snapshot$input_mode == "model_based") {
    specs$pd <- list(label = "Modeled-allele frequency", min = 0.01, max = 0.99, value = values$pd / 100, key = "pd")
    specs$prev <- list(label = "Disease prevalence", min = 0.001, max = 0.5, value = values$prev / 100, key = "prev")
    specs$R1 <- list(label = "Heterozygote relative risk", min = 0.1, max = max(4, values$R1 * 2), value = values$R1, key = "R1")
    specs$R2 <- list(label = "Homozygote relative risk", min = 0.1, max = max(5, values$R2 * 2), value = values$R2, key = "R2")
    specs$delta_prime <- list(label = "D-prime", min = 0, max = 1, value = values$delta_prime, key = "delta_prime")
  } else {
    specs$ET <- list(label = "Expected transmitted count", min = 0, max = max(10, values$ET * 2), value = values$ET, key = "ET")
    specs$ENT <- list(label = "Expected non-transmitted count", min = 0, max = max(10, values$ENT * 2), value = values$ENT, key = "ENT")
    if (snapshot$objective == "mssn") specs$n_trios <- list(label = "Represented affected-child trios", min = 1, max = max(10, values$n_trios * 2), value = values$n_trios, key = "n_trios")
  }
  specs$alpha <- list(label = "Significance level", min = 1e-04, max = 0.2, value = values$alpha, key = "alpha")
  if (snapshot$objective == "power") {
    specs$N <- list(label = "Number of affected-child trios", min = max(10, snapshot$objective_value / 4), max = snapshot$objective_value * 2, value = snapshot$objective_value, objective = TRUE)
  } else {
    specs$target_power <- list(label = "Target power", min = 0.5, max = 0.99, value = snapshot$objective_value, objective = TRUE)
  }
  specs
}

.pawh_tdt_sensitivity <- function(calculation, parameter, range, n = 40L) {
  spec <- .pawh_tdt_specs(calculation)[[parameter]]
  if (is.null(spec)) stop("Unavailable sensitivity parameter.", call. = FALSE)
  scenarios <- .pawh_tdt_scenarios(calculation)
  grid <- seq(range[1], range[2], length.out = n)
  data <- do.call(rbind, lapply(grid, function(value) {
    snapshot <- calculation$snapshot
    args <- snapshot$backend_args
    if (isTRUE(spec$objective)) snapshot$objective_value <- value else args[[spec$key]] <- value
    result <- tryCatch(.pawh_tdt_call(snapshot, args), error = function(e) NULL)
    y <- if (is.null(result)) rep(NA_real_, length(scenarios)) else if (snapshot$objective == "power") {
      unlist(result$power[scenarios], use.names = FALSE)
    } else unlist(result$N[scenarios], use.names = FALSE)
    data.frame(
      x = value, Scenario = unname(.pawh_tdt_scenario_labels[scenarios]),
      scenario = scenarios, y = y
    )
  }))
  plot_data <- data
  plot_data$y <- ifelse(is.finite(plot_data$y), plot_data$y, NA_real_)
  list(
    data = data, plot_data = plot_data,
    label = spec$label, baseline_x = spec$value,
    objective = calculation$snapshot$objective,
    has_nonfinite = any(!is.finite(data$y))
  )
}

.pawh_tdt_sensitivity_plot <- function(sensitivity) {
  colors <- .pawh_plot_colors()
  values <- c(
    `No-error design` = unname(colors["genotype"]),
    `Phenotype misclassification` = unname(colors["trend"]),
    `Locus heterogeneity` = unname(colors["controls"])
  )
  lines <- c(
    `No-error design` = "solid", `Phenotype misclassification` = "dashed",
    `Locus heterogeneity` = "dotdash"
  )
  ggplot2::ggplot(sensitivity$plot_data, ggplot2::aes(
    .data$x, .data$y, color = .data$Scenario, linetype = .data$Scenario
  )) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_vline(
      xintercept = sensitivity$baseline_x, color = unname(colors["reference"]),
      linetype = "dotted", linewidth = 0.6
    ) +
    ggplot2::scale_color_manual(values = values) +
    ggplot2::scale_linetype_manual(values = lines) +
    ggplot2::labs(
      x = sensitivity$label,
      y = if (sensitivity$objective == "power") "Power" else "Required affected-child trios",
      color = NULL, linetype = NULL
    ) +
    .pawh_plot_theme() + ggplot2::theme(legend.position = "top")
}

.pawh_tdt_transmissions <- function(calculation) {
  scenarios <- .pawh_tdt_scenarios(calculation)
  result <- calculation$result
  do.call(rbind, lapply(scenarios, function(scenario) data.frame(
    Scenario = unname(.pawh_tdt_scenario_labels[scenario]),
    Quantity = c("Transmitted", "Non-transmitted"),
    Probability = c(result$gT_star[[scenario]], result$gNT_star[[scenario]])
  )))
}

.pawh_tdt_transmission_plot <- function(calculation) {
  colors <- .pawh_plot_colors()
  ggplot2::ggplot(.pawh_tdt_transmissions(calculation), ggplot2::aes(
    .data$Quantity, .data$Probability, fill = .data$Quantity
  )) +
    ggplot2::geom_col(width = 0.68, alpha = 0.92) +
    ggplot2::facet_wrap(~Scenario) +
    ggplot2::scale_fill_manual(values = c(
      Transmitted = unname(colors["genotype"]),
      `Non-transmitted` = unname(colors["controls"])
    )) +
    ggplot2::labs(x = NULL, y = "Expected transmission probability", fill = NULL) +
    .pawh_plot_theme() +
    ggplot2::theme(
      legend.position = "top", panel.grid.major.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 10, face = "plain", color = "#3F4850")
    )
}

.pawh_tdt_summary_ui <- function(calculation) {
  snapshot <- calculation$snapshot
  values <- snapshot$display
  modifiers <- c(
    if (calculation$active$misclassification) paste0("Misclassification ", .pawh_tdt_rate(snapshot$backend_args$misclass_rate)),
    if (calculation$active$heterogeneity) paste0("Heterogeneity ", .pawh_tdt_rate(snapshot$backend_args$heter_rate))
  )
  common <- list(
    .pawh_summary_row("Study type", "Affected-child trio / TDT"),
    .pawh_summary_row("Objective", if (snapshot$objective == "power") {
      paste("Estimate power with", .pawh_tdt_count(snapshot$objective_value, TRUE), "trios")
    } else paste("Minimum sample size at", .pawh_tdt_pct(snapshot$objective_value))),
    .pawh_summary_row("Inputs", if (snapshot$input_mode == "model_based") "Genetic model" else "Direct transmission quantities"),
    .pawh_summary_row("Alpha", format(values$alpha))
  )
  specific <- if (snapshot$input_mode == "model_based") list(
    .pawh_summary_row("Disease prevalence", paste0(values$prev, "%")),
    .pawh_summary_row("Allele frequency", paste0(values$pd, "%")),
    .pawh_summary_row("Relative risks", paste0("R1 ", values$R1, "; R2 ", values$R2)),
    .pawh_summary_row("D-prime", format(values$delta_prime))
  ) else list(
    .pawh_summary_row("Transmission counts", paste0("ET ", values$ET, "; ENT ", values$ENT)),
    if (snapshot$objective == "mssn") .pawh_summary_row("Represented trios", .pawh_tdt_count(values$n_trios, TRUE))
  )
  shiny::div(
    class = "pawh-summary-grid", common, specific,
    .pawh_summary_row("Modifiers", if (length(modifiers)) paste(modifiers, collapse = "; ") else "None")
  )
}

.pawh_tdt_methods_ui <- function(calculation = NULL) {
  if (is.null(calculation)) return(.pawh_placeholder_ui("Methods", "Calculate a design to record its analysis specification."))
  scenarios <- unname(.pawh_tdt_scenario_labels[.pawh_tdt_scenarios(calculation)])
  shiny::tagList(
    shiny::h3("Methods"),
    shiny::div(
      class = "pawh-summary-grid",
      .pawh_summary_row("Study design", "Affected-child trio / TDT"),
      .pawh_summary_row("Objective", if (calculation$snapshot$objective == "power") "Power" else "Minimum sample size"),
      .pawh_summary_row("Input specification", if (calculation$snapshot$input_mode == "model_based") "Genetic model" else "Direct transmission quantities"),
      .pawh_summary_row("Statistical method", "Transmission Disequilibrium Test"),
      .pawh_summary_row("Scenarios", paste(scenarios, collapse = "; "))
    ),
    if (all(unlist(calculation$active))) shiny::p(
      class = "pawh-caution",
      "The current implementation evaluates phenotype misclassification and locus heterogeneity separately rather than as a joint combined model."
    ),
    shiny::p(
      "Calculations use ", shiny::code(if (calculation$snapshot$objective == "power") "tdt_power()" else "tdt_mssn()"),
      ". Sensitivity points are separate canonical calls using the frozen design."
    ),
    shiny::a(
      href = "https://akilanthony.github.io/pawh/articles/pawh-03-tdt-study-design.html",
      target = "_blank", rel = "noopener", "Read the TDT study-design vignette"
    )
  )
}

.pawh_tdt_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    state <- shiny::reactiveValues(calculation = NULL, error = NULL, signature = NULL, sensitivity = NULL)
    values <- shiny::reactive(.pawh_tdt_values(input))
    changed <- shiny::reactive(!is.null(state$signature) && !identical(state$signature, .pawh_tdt_sig(values())))

    output$objective_inputs <- shiny::renderUI(if (is.null(input$objective) || input$objective == "power") {
      shiny::numericInput(ns("N"), "Number of affected-child trios", 600, min = 1, step = 10)
    } else {
      shiny::sliderInput(ns("target_power"), "Target power (%)", min = 50, max = 99, value = 80)
    })
    output$model_inputs <- shiny::renderUI({
      direct <- identical(input$input_mode, "model_free")
      if (!direct) shiny::tagList(
        shiny::numericInput(ns("prev"), "Disease prevalence (%)", 5, 0.01, 99.99, 0.1),
        shiny::tags$small(class = "text-muted", "Population prevalence used by the canonical penetrance model."),
        shiny::numericInput(ns("pd"), "Modeled-allele frequency (%)", 30, 0.01, 99.99, 0.1),
        shiny::numericInput(ns("R1"), "Heterozygote relative risk (R1)", 1.5, 0.01, NA, 0.05),
        shiny::numericInput(ns("R2"), "Homozygote relative risk (R2)", 2.25, 0.01, NA, 0.05),
        shiny::tags$small(class = "text-muted", "R1 and R2 are entered explicitly; the canonical TDT API does not define inheritance shortcuts."),
        shiny::numericInput(ns("delta_prime"), "Linkage disequilibrium (D-prime)", 1, 0, 1, 0.05)
      ) else shiny::tagList(
        shiny::numericInput(ns("ET"), "Expected transmitted count (ET)", 140, 0, NA, 1),
        shiny::numericInput(ns("ENT"), "Expected non-transmitted count (ENT)", 100, 0, NA, 1),
        if (identical(input$objective, "mssn")) shiny::numericInput(ns("n_trios"), "Affected-child trios represented by ET/ENT", 120, 1, NA, 1),
        shiny::tags$small(class = "text-muted", "Enter the canonical expected transmission and non-transmission counts, not genotype probabilities."),
        if (isTRUE(input$misclassification) || isTRUE(input$heterogeneity)) shiny::numericInput(ns("pd"), "Modeled-allele frequency for modifier scenario (%)", 30, 0.01, 99.99, 0.1),
        if (isTRUE(input$misclassification)) shiny::numericInput(ns("prev"), "Disease prevalence for misclassification scenario (%)", 5, 0.01, 99.99, 0.1)
      )
    })

    shiny::observeEvent(input$calculate, {
      state$error <- NULL
      state$sensitivity <- NULL
      submitted <- values()
      tryCatch({
        state$calculation <- .pawh_tdt_calculate(.pawh_tdt_snapshot(submitted))
        state$signature <- .pawh_tdt_sig(submitted)
      }, error = function(error) {
        state$error <- paste("This design could not be calculated.", conditionMessage(error))
      })
    }, ignoreInit = TRUE)

    output$changed_notice <- shiny::renderUI(if (changed()) shiny::div(
      class = "pawh-changed-notice", role = "status",
      "Inputs have changed. Recalculate to update results."
    ))
    output$design_summary <- shiny::renderUI(if (!is.null(state$error)) {
      shiny::div(class = "pawh-error", role = "alert", state$error)
    } else if (is.null(state$calculation)) {
      shiny::p(class = "text-muted", "Choose assumptions and select Calculate.")
    } else .pawh_tdt_summary_ui(state$calculation))
    output$results <- shiny::renderUI(if (!is.null(state$error)) {
      shiny::div(class = "pawh-error", role = "alert", state$error)
    } else if (is.null(state$calculation)) {
      .pawh_placeholder_ui("Results", "Set up a design and select Calculate.")
    } else .pawh_tdt_results_ui(state$calculation))
    output$sensitivity_controls <- shiny::renderUI({
      if (is.null(state$calculation)) return(.pawh_placeholder_ui("Sensitivity", "Calculate a design first."))
      specs <- .pawh_tdt_specs(state$calculation)
      choices <- stats::setNames(names(specs), vapply(specs, `[[`, "", "label"))
      spec <- specs[[1L]]
      shiny::div(
        class = "pawh-sensitivity-controls",
        shiny::selectInput(ns("sensitivity_parameter"), "Parameter", choices),
        shiny::sliderInput(ns("sensitivity_range"), "Range", spec$min, spec$max,
          pmax(spec$min, pmin(spec$max, c(spec$value * 0.75, spec$value * 1.25)))),
        shiny::actionButton(ns("run_sensitivity"), "Run sensitivity analysis")
      )
    })
    shiny::observeEvent(input$sensitivity_parameter, {
      shiny::req(state$calculation)
      spec <- .pawh_tdt_specs(state$calculation)[[input$sensitivity_parameter]]
      shiny::req(spec)
      value <- pmax(spec$min, pmin(spec$max, c(spec$value * 0.75, spec$value * 1.25)))
      shiny::updateSliderInput(session, "sensitivity_range", min = spec$min, max = spec$max,
        value = value, step = (spec$max - spec$min) / 100)
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$run_sensitivity, {
      shiny::req(state$calculation, input$sensitivity_parameter, input$sensitivity_range)
      state$sensitivity <- .pawh_tdt_sensitivity(state$calculation, input$sensitivity_parameter, input$sensitivity_range)
    }, ignoreInit = TRUE)
    output$sensitivity_message <- shiny::renderUI(if (!is.null(state$calculation) && is.null(state$sensitivity)) {
      shiny::p(class = "text-muted", "Each point is a canonical calculation of the frozen design.")
    } else if (!is.null(state$sensitivity) && state$sensitivity$has_nonfinite) {
      shiny::div(class = "pawh-caution", "Some explored designs have infinite or undefined required sample size and are omitted from the line geometry.")
    })
    output$sensitivity_plot <- shiny::renderPlot({
      shiny::req(state$sensitivity)
      .pawh_tdt_sensitivity_plot(state$sensitivity)
    })
    output$visualize_intro <- shiny::renderUI(if (is.null(state$calculation)) {
      .pawh_placeholder_ui("Visualize", "Calculate a design first.")
    } else shiny::div(
      class = "pawh-visual-intro", shiny::h3("Transmission imbalance"),
      shiny::p("Expected transmission and non-transmission probabilities returned by the canonical calculation.")
    ))
    output$transmission_plot <- shiny::renderPlot({
      shiny::req(state$calculation)
      .pawh_tdt_transmission_plot(state$calculation)
    })
    output$methods <- shiny::renderUI(.pawh_tdt_methods_ui(state$calculation))

    list(
      calculation = shiny::reactive(state$calculation), changed = changed,
      error = shiny::reactive(state$error), sensitivity = shiny::reactive(state$sensitivity)
    )
  })
}
