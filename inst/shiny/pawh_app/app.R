if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The 'shiny' package is required to run this app.", call. = FALSE)
}

if (!requireNamespace("pawh", quietly = TRUE)) {
  stop("The 'pawh' package must be installed to run this app.", call. = FALSE)
}

app_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
app_dir <- if (!is.null(app_file)) {
  normalizePath(dirname(app_file), mustWork = FALSE)
} else {
  getwd()
}
pkg_root <- normalizePath(file.path(app_dir, "../../.."), mustWork = FALSE)
local_pkg_available <- file.exists(file.path(pkg_root, "DESCRIPTION")) &&
  file.exists(file.path(pkg_root, "R", "case_control.R"))

required_exports <- c(
  "cc_power",
  "cc_mssn",
  "plot_cc_power",
  "plot_cc_mssn",
  "tdt_power",
  "tdt_mssn",
  "plot_tdt_power",
  "plot_tdt_mssn"
)

missing_exports <- setdiff(required_exports, getNamespaceExports("pawh"))
if (local_pkg_available &&
    length(missing_exports) > 0 &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_root, quiet = TRUE)
  missing_exports <- setdiff(required_exports, getNamespaceExports("pawh"))
}

if (length(missing_exports) > 0) {
  stop(
    paste(
      "The loaded pawh package is missing required Shiny app exports:",
      paste(missing_exports, collapse = ", "),
      "Install this branch or run devtools::load_all() before launching the app."
    ),
    call. = FALSE
  )
}

library(shiny)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

test_label_map <- c(
  genotypes = "Genotype chi-square",
  trend = "Trend test"
)

scenario_label_map <- c(
  no_error = "No error",
  misclassification = "Misclassification",
  heterogeneity = "Heterogeneity",
  auto = "Auto"
)

choice_card <- function(id, title, subtitle = NULL) {
  actionButton(
    inputId = id,
    label = tagList(
      div(class = "choice-title", title),
      if (!is.null(subtitle)) div(class = "choice-subtitle", subtitle)
    ),
    class = "choice-card"
  )
}

card_grid <- function(...) {
  div(class = "choice-grid", ...)
}

section_card <- function(title, ..., class = NULL) {
  div(
    class = paste(c("panel-card", class), collapse = " "),
    h3(title),
    ...
  )
}

nav_controls <- function() {
  div(
    class = "nav-controls",
    actionButton("back_btn", "Back", class = "secondary-btn"),
    actionButton("reset_btn", "Start over", class = "secondary-btn")
  )
}

advanced_panel <- function(title, ...) {
  tags$details(
    class = "advanced-panel",
    tags$summary(HTML(paste0("&#9881; ", title))),
    div(class = "advanced-body", ...)
  )
}

fmt_num <- function(x, digits = 4) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return(NA_character_)
  }
  if (!is.finite(x)) {
    return(as.character(x))
  }
  if (abs(x) >= 1000) {
    formatC(x, format = "f", digits = 0, big.mark = ",")
  } else {
    formatC(x, format = "f", digits = digits)
  }
}

nonempty_tests <- function(x) {
  if (is.null(x) || length(x) == 0) {
    stop("Select at least one case-control test.")
  }
  x
}

make_range <- function(min_value, max_value, step_size) {
  if (!is.numeric(min_value) || !is.numeric(max_value) || !is.numeric(step_size) ||
      any(!is.finite(c(min_value, max_value, step_size)))) {
    stop("Plot range values must be finite numbers.")
  }
  if (step_size <= 0) {
    stop("Step size must be positive.")
  }
  if (max_value <= min_value) {
    stop("Maximum must be greater than minimum.")
  }
  values <- seq(min_value, max_value, by = step_size)
  if (length(values) < 2) {
    stop("The plot range must produce at least two x-axis values.")
  }
  values
}

safe_call <- function(expr) {
  tryCatch(
    list(result = force(expr), error = NULL),
    error = function(e) list(result = NULL, error = conditionMessage(e))
  )
}

render_error_or_table <- function(payload, table) {
  if (!is.null(payload$error)) {
    return(data.frame(Message = payload$error))
  }
  table
}

tdt_power_table <- function(out) {
  data.frame(
    Scenario = c("No error", "Misclassification", "Heterogeneity"),
    Power = vapply(out$power, fmt_num, character(1), digits = 5),
    Lambda = vapply(out$lambda, fmt_num, character(1), digits = 4),
    `Power loss` = c(
      NA_character_,
      fmt_num(out$power_loss$misclassification, 5),
      fmt_num(out$power_loss$heterogeneity, 5)
    ),
    check.names = FALSE
  )
}

tdt_mssn_table <- function(out) {
  data.frame(
    Scenario = c("No error", "Misclassification", "Heterogeneity"),
    `Required trios` = vapply(out$N, fmt_num, character(1), digits = 1),
    `Percent increase` = c(
      NA_character_,
      fmt_num(out$percent_increase$misclassification, 2),
      fmt_num(out$percent_increase$heterogeneity, 2)
    ),
    `Power at no-error N` = vapply(out$power_at_N_no_error, fmt_num, character(1), digits = 5),
    check.names = FALSE
  )
}

cc_power_table <- function(out, selected_tests) {
  rows <- lapply(selected_tests, function(test) {
    result <- out$tests[[test]]
    if (is.null(result)) {
      return(NULL)
    }
    data.frame(
      Test = unname(test_label_map[test]),
      Power = fmt_num(result$power, 5),
      Lambda = fmt_num(result$lambda, 4),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

cc_mssn_table <- function(out, selected_tests) {
  rows <- lapply(selected_tests, function(test) {
    result <- out$tests[[test]]
    if (is.null(result)) {
      return(NULL)
    }
    data.frame(
      Test = unname(test_label_map[test]),
      MSSN_case = result$MSSN_case,
      MSSN_ctrl = result$MSSN_ctrl,
      MSSN_total = result$MSSN_total,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

mode_specific_cc_ui <- function(input_mode) {
  if (identical(input_mode, "model_based")) {
    tagList(
      sliderInput("cc_prev", "Disease prevalence (prev)", min = 0.001, max = 0.50, value = 0.05, step = 0.001),
      sliderInput("cc_pd", "Risk-allele frequency (pd)", min = 0.01, max = 0.99, value = 0.30, step = 0.01),
      numericInput("cc_R2", "Homozygote relative risk (R2)", 1.80, min = 0.0001, step = 0.1),
      selectInput("cc_MOI", "Mode of inheritance (MOI)", choices = c("M", "D", "Rec"))
    )
  } else {
    tagList(
      div(class = "input-note", "If phenotype misclassification is enabled, g1 is true affected and g0 is true unaffected."),
      fluidRow(
        column(
          6,
          h4("g1 frequencies"),
          numericInput("cc_g1_0", "g1: genotype 0", 0.25, min = 0, max = 1, step = 0.01),
          numericInput("cc_g1_1", "g1: genotype 1", 0.50, min = 0, max = 1, step = 0.01),
          numericInput("cc_g1_2", "g1: genotype 2", 0.25, min = 0, max = 1, step = 0.01)
        ),
        column(
          6,
          h4("g0 frequencies"),
          numericInput("cc_g0_0", "g0: genotype 0", 0.36, min = 0, max = 1, step = 0.01),
          numericInput("cc_g0_1", "g0: genotype 1", 0.48, min = 0, max = 1, step = 0.01),
          numericInput("cc_g0_2", "g0: genotype 2", 0.16, min = 0, max = 1, step = 0.01)
        )
      ),
      sliderInput(
        "cc_prev_model_free",
        "Disease prevalence for phenotype misclassification",
        min = 0.001,
        max = 0.50,
        value = 0.05,
        step = 0.001
      )
    )
  }
}

cc_advanced_ui <- function() {
  advanced_panel(
    "Advanced heterogeneity and error settings",
    checkboxInput("cc_locus_het", "Enable locus heterogeneity", FALSE),
    sliderInput("cc_locus_het_rate", "Locus heterogeneity rate (1 - pi)", min = 0, max = 1, value = 0.20, step = 0.01),
    checkboxInput("cc_pheno_misclass", "Enable phenotype misclassification", FALSE),
    sliderInput("cc_theta", "theta: affected classified as control", min = 0, max = 0.50, value = 0, step = 0.01),
    sliderInput("cc_phi", "phi: unaffected classified as case", min = 0, max = 0.50, value = 0, step = 0.01),
    selectInput("cc_geno_misclass", "Genotype misclassification model", c("none", "1p", "2p", "3p", "diff3p")),
    uiOutput("cc_geno_error_ui")
  )
}

tdt_common_inputs <- function(include_N = FALSE, include_target_power = FALSE) {
  tagList(
    if (include_N) numericInput("tdt_N", "Number of affected trios (N)", 600, min = 1, step = 25),
    if (include_target_power) sliderInput("tdt_target_power", "Target power", min = 0.50, max = 0.99, value = 0.80, step = 0.01),
    sliderInput("tdt_pd", "Risk-allele frequency (pd)", min = 0.01, max = 0.99, value = 0.30, step = 0.01),
    sliderInput("tdt_prev", "Disease prevalence (prev)", min = 0.001, max = 0.50, value = 0.05, step = 0.001),
    numericInput("tdt_R1", "Heterozygote relative risk (R1)", 1.50, min = 0.0001, step = 0.1),
    numericInput("tdt_R2", "Homozygote relative risk (R2)", 2.25, min = 0.0001, step = 0.1),
    selectInput("tdt_alpha", "Significance level (alpha)", choices = c(0.10, 0.05, 0.01, 0.001), selected = 0.05),
    sliderInput("tdt_delta_prime", "LD scale parameter (delta_prime)", min = 0, max = 1, value = 1, step = 0.01),
    sliderInput("tdt_misclass_rate", "Phenotype misclassification rate", min = 0, max = 0.20, value = 0.01, step = 0.005),
    sliderInput("tdt_heter_rate", "Heterogeneity rate (1 - pi)", min = 0, max = 0.80, value = 0.10, step = 0.01)
  )
}

plot_settings_ui <- function(prefix) {
  advanced_panel(
    "Advanced plot settings",
    textInput(paste0(prefix, "_title"), "Custom title", ""),
    textInput(paste0(prefix, "_x_label"), "Custom x-axis label", ""),
    textInput(paste0(prefix, "_y_label"), "Custom y-axis label", "")
  )
}

alpha_gate_ui <- function() {
  div(
    class = "alpha-gate",
    div(
      class = "alpha-box",
      textInput("alpha_key", NULL, placeholder = "Enter alpha key"),
      actionButton("alpha_key_submit", "Enter", class = "primary-btn alpha-submit"),
      div(class = "alpha-error", textOutput("alpha_error", inline = TRUE))
    )
  )
}

landing_ui <- function() {
  tagList(
    div(class = "hero-card",
        h2("Choose a study design"),
        p("Start with the study family. The app will reveal only the settings needed for the selected workflow.")),
    card_grid(
      choice_card("choose_tdt", "TDT", "Transmission Disequilibrium Test workflows"),
      choice_card("choose_cc", "Case-Control", "Case-control power, sample size, and plots")
    )
  )
}

tdt_choice_ui <- function() {
  tagList(
    nav_controls(),
    div(class = "hero-card",
        h2("TDT workflow"),
        p("Choose whether to compute power, required trios, or make an exploratory plot.")),
    card_grid(
      choice_card("tdt_power_card", "Power", "Compute TDT power for a fixed number of affected trios"),
      choice_card("tdt_mssn_card", "Sample Size", "Compute required affected trios for a target power"),
      choice_card("tdt_plots_card", "Plots", "Sweep a parameter and plot the full TDT backend")
    )
  )
}

cc_setup_ui <- function(input_mode) {
  tagList(
    nav_controls(),
    div(class = "hero-card",
        h2("Case-Control setup"),
        p("Select the input type and one or more tests before choosing a workflow.")),
    h3("Input type"),
    card_grid(
      choice_card("cc_mode_model_based", "Model-based", "Use prevalence, allele frequency, relative risk, and inheritance mode"),
      choice_card("cc_mode_model_free", "Model-free", "Supply genotype frequencies directly")
    ),
    div(
      class = "setup-panel",
      strong("Current input type: "),
      span(if (is.null(input_mode)) "Not selected" else if (input_mode == "model_based") "Model-based" else "Model-free")
    ),
    checkboxGroupInput(
      "cc_tests",
      "Tests to include",
      choices = c(
        "Genotype chi-square" = "genotypes",
        "Trend test" = "trend"
      ),
      selected = c("genotypes", "trend")
    ),
    actionButton("cc_setup_continue", "Continue", class = "primary-btn")
  )
}

cc_choice_ui <- function(input_mode, selected_tests) {
  tagList(
    nav_controls(),
    div(class = "hero-card",
        h2("Case-Control workflow"),
        p(paste(
          "Input:",
          if (input_mode == "model_based") "model-based" else "model-free",
          "| Tests:",
          paste(unname(test_label_map[selected_tests]), collapse = ", ")
        ))),
    card_grid(
      choice_card("cc_power_card", "Power", "Compute power for fixed case-control sample sizes"),
      choice_card("cc_mssn_card", "Sample Size", "Compute minimum sample size for a target power"),
      choice_card("cc_plots_card", "Plots", "Sweep a parameter using the full case-control plotting wrappers")
    )
  )
}

tdt_power_ui <- function() {
  tagList(
    nav_controls(),
    div(class = "workflow-grid",
        section_card(
          "TDT power inputs",
          tdt_common_inputs(include_N = TRUE),
          selectInput(
            "tdt_power_plot_scenario",
            "Plot scenario",
            choices = c("No error" = "no_error", "Misclassification" = "misclassification", "Heterogeneity" = "heterogeneity"),
            selected = "no_error"
          ),
          actionButton("tdt_power_analyze", "Analyze", class = "primary-btn")
        ),
        section_card(
          "Results",
          tableOutput("tdt_power_table"),
          plotOutput("tdt_power_plot", height = "360px")
        )
    )
  )
}

tdt_mssn_ui <- function() {
  tagList(
    nav_controls(),
    div(class = "workflow-grid",
        section_card(
          "TDT sample-size inputs",
          tdt_common_inputs(include_target_power = TRUE),
          selectInput(
            "tdt_mssn_plot_scenario",
            "Plot scenario",
            choices = c("Misclassification" = "misclassification", "Heterogeneity" = "heterogeneity"),
            selected = "heterogeneity"
          ),
          actionButton("tdt_mssn_analyze", "Analyze", class = "primary-btn")
        ),
        section_card(
          "Results",
          tableOutput("tdt_mssn_table"),
          plotOutput("tdt_mssn_plot", height = "360px")
        )
    )
  )
}

tdt_plots_ui <- function() {
  tagList(
    nav_controls(),
    div(class = "workflow-grid",
        section_card(
          "TDT plot controls",
          selectInput("tdt_plot_type", "Plot type", c("Power plot" = "power", "Sample-size plot" = "mssn")),
          uiOutput("tdt_plot_xvar_ui"),
          fluidRow(
            column(4, numericInput("tdt_plot_xmin", "Minimum", 0, step = 0.01)),
            column(4, numericInput("tdt_plot_xmax", "Maximum", 0.50, step = 0.01)),
            column(4, numericInput("tdt_plot_xstep", "Step size", 0.05, min = 0.0001, step = 0.01))
          ),
          selectInput(
            "tdt_plot_scenario",
            "Scenario",
            choices = c("Auto" = "auto", "No error" = "no_error", "Misclassification" = "misclassification", "Heterogeneity" = "heterogeneity"),
            selected = "auto"
          ),
          conditionalPanel(
            "input.tdt_plot_type == 'power'",
            numericInput("tdt_plot_N", "Number of affected trios (N)", 600, min = 1, step = 25)
          ),
          conditionalPanel(
            "input.tdt_plot_type == 'mssn'",
            sliderInput("tdt_plot_target_power", "Target power", min = 0.50, max = 0.99, value = 0.80, step = 0.01)
          ),
          sliderInput("tdt_plot_pd", "Risk-allele frequency (pd)", min = 0.01, max = 0.99, value = 0.30, step = 0.01),
          sliderInput("tdt_plot_prev", "Disease prevalence (prev)", min = 0.001, max = 0.50, value = 0.05, step = 0.001),
          numericInput("tdt_plot_R1", "Heterozygote relative risk (R1)", 1.50, min = 0.0001, step = 0.1),
          numericInput("tdt_plot_R2", "Homozygote relative risk (R2)", 2.25, min = 0.0001, step = 0.1),
          selectInput("tdt_plot_alpha", "Significance level (alpha)", choices = c(0.10, 0.05, 0.01, 0.001), selected = 0.05),
          sliderInput("tdt_plot_delta_prime", "LD scale parameter (delta_prime)", min = 0, max = 1, value = 1, step = 0.01),
          sliderInput("tdt_plot_misclass_rate", "Phenotype misclassification rate", min = 0, max = 0.20, value = 0.01, step = 0.005),
          sliderInput("tdt_plot_heter_rate", "Heterogeneity rate (1 - pi)", min = 0, max = 0.80, value = 0.10, step = 0.01),
          plot_settings_ui("tdt_plot"),
          actionButton("tdt_plot_generate", "Generate Plot", class = "primary-btn")
        ),
        section_card("Plot", plotOutput("tdt_free_plot", height = "460px"), tableOutput("tdt_plot_error"))
    )
  )
}

cc_power_ui <- function(input_mode, selected_tests) {
  tagList(
    nav_controls(),
    div(class = "workflow-grid",
        section_card(
          "Case-control power inputs",
          numericInput("cc_N_case", "Number of cases (N_case)", 250, min = 1, step = 25),
          selectInput("cc_alpha", "Significance level (alpha)", choices = c(0.10, 0.05, 0.01, 0.001), selected = 0.05),
          numericInput("cc_k", "Control-to-case ratio (k)", 1, min = 0.001, step = 0.1),
          mode_specific_cc_ui(input_mode),
          fluidRow(
            column(4, numericInput("cc_w0", "Trend weight w0", 0, step = 1)),
            column(4, numericInput("cc_w1", "Trend weight w1", 1, step = 1)),
            column(4, numericInput("cc_w2", "Trend weight w2", 2, step = 1))
          ),
          cc_advanced_ui(),
          div(class = "input-note", paste("Selected tests:", paste(unname(test_label_map[selected_tests]), collapse = ", "))),
          actionButton("cc_power_analyze", "Analyze", class = "primary-btn")
        ),
        section_card(
          "Results",
          tableOutput("cc_power_table"),
          plotOutput("cc_power_plot", height = "360px")
        )
    )
  )
}

cc_mssn_ui <- function(input_mode, selected_tests) {
  tagList(
    nav_controls(),
    div(class = "workflow-grid",
        section_card(
          "Case-control sample-size inputs",
          sliderInput("cc_target_power", "Target power", min = 0.50, max = 0.99, value = 0.80, step = 0.01),
          selectInput("cc_alpha", "Significance level (alpha)", choices = c(0.10, 0.05, 0.01, 0.001), selected = 0.05),
          numericInput("cc_k", "Control-to-case ratio (k)", 1, min = 0.001, step = 0.1),
          mode_specific_cc_ui(input_mode),
          fluidRow(
            column(4, numericInput("cc_w0", "Trend weight w0", 0, step = 1)),
            column(4, numericInput("cc_w1", "Trend weight w1", 1, step = 1)),
            column(4, numericInput("cc_w2", "Trend weight w2", 2, step = 1))
          ),
          cc_advanced_ui(),
          div(class = "input-note", paste("Selected tests:", paste(unname(test_label_map[selected_tests]), collapse = ", "))),
          actionButton("cc_mssn_analyze", "Analyze", class = "primary-btn")
        ),
        section_card(
          "Results",
          tableOutput("cc_mssn_table"),
          plotOutput("cc_mssn_plot", height = "360px")
        )
    )
  )
}

cc_plots_ui <- function(input_mode, selected_tests) {
  tagList(
    nav_controls(),
    div(class = "workflow-grid",
        section_card(
          "Case-control plot controls",
          selectInput("cc_plot_type", "Plot type", c("Power plot" = "power", "Sample-size plot" = "mssn")),
          uiOutput("cc_plot_xvar_ui"),
          fluidRow(
            column(4, numericInput("cc_plot_xmin", "Minimum", 0, step = 0.01)),
            column(4, numericInput("cc_plot_xmax", "Maximum", 0.10, step = 0.01)),
            column(4, numericInput("cc_plot_xstep", "Step size", 0.01, min = 0.0001, step = 0.01))
          ),
          conditionalPanel(
            "input.cc_plot_type == 'power'",
            numericInput("cc_plot_N_case", "Number of cases (N_case)", 250, min = 1, step = 25)
          ),
          conditionalPanel(
            "input.cc_plot_type == 'mssn'",
            sliderInput("cc_plot_target_power", "Target power", min = 0.50, max = 0.99, value = 0.80, step = 0.01)
          ),
          selectInput("cc_plot_alpha", "Significance level (alpha)", choices = c(0.10, 0.05, 0.01, 0.001), selected = 0.05),
          numericInput("cc_plot_k", "Control-to-case ratio (k)", 1, min = 0.001, step = 0.1),
          mode_specific_cc_ui(input_mode),
          fluidRow(
            column(4, numericInput("cc_w0", "Trend weight w0", 0, step = 1)),
            column(4, numericInput("cc_w1", "Trend weight w1", 1, step = 1)),
            column(4, numericInput("cc_w2", "Trend weight w2", 2, step = 1))
          ),
          cc_advanced_ui(),
          plot_settings_ui("cc_plot"),
          div(class = "input-note", paste("Selected tests:", paste(unname(test_label_map[selected_tests]), collapse = ", "))),
          actionButton("cc_plot_generate", "Generate Plot", class = "primary-btn")
        ),
        section_card("Plot", plotOutput("cc_free_plot", height = "460px"), tableOutput("cc_plot_error"))
    )
  )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background: #f5f7fb;
        color: #172b4d;
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
      }
      .app-shell {
        max-width: 1240px;
        margin: 0 auto;
        padding: 28px 20px 56px;
      }
      .app-header {
        display: flex;
        justify-content: space-between;
        gap: 18px;
        align-items: flex-start;
        margin-bottom: 20px;
      }
      .app-header h1 {
        margin: 0 0 8px;
        font-weight: 750;
        letter-spacing: 0;
      }
      .alpha-gate {
        min-height: 82vh;
        display: flex;
        align-items: center;
        justify-content: center;
      }
      .alpha-box {
        width: min(340px, 92vw);
      }
      .alpha-box .form-group {
        margin-bottom: 10px;
      }
      .alpha-box input {
        height: 46px;
        border-radius: 10px;
        border: 1px solid #c8d2df;
        font-size: 16px;
      }
      .alpha-submit {
        width: 100%;
        margin-top: 0 !important;
      }
      .alpha-error {
        min-height: 24px;
        margin-top: 12px;
        color: #b42318;
        font-weight: 700;
        text-align: center;
      }
      .breadcrumb-line {
        color: #52627a;
        font-size: 15px;
      }
      .hero-card, .panel-card, .setup-panel {
        background: #ffffff;
        border: 1px solid #d9e2ec;
        border-radius: 16px;
        box-shadow: 0 8px 24px rgba(23, 43, 77, 0.06);
      }
      .hero-card {
        padding: 26px 30px;
        margin-bottom: 22px;
      }
      .hero-card h2, .panel-card h3 {
        margin-top: 0;
      }
      .choice-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
        gap: 18px;
      }
      .choice-card {
        width: 100%;
        min-height: 150px;
        border: 1px solid #d9e2ec !important;
        border-radius: 16px !important;
        padding: 28px !important;
        text-align: center;
        cursor: pointer;
        background: #ffffff !important;
        color: #172b4d !important;
        transition: transform 0.18s ease, box-shadow 0.18s ease, border-color 0.18s ease;
        white-space: normal;
      }
      .choice-card:hover {
        transform: scale(1.03);
        box-shadow: 0 10px 24px rgba(0,0,0,0.12);
        border-color: #2c7be5 !important;
      }
      .choice-card.selected {
        border-color: #2c7be5 !important;
        box-shadow: 0 10px 24px rgba(44,123,229,0.18);
      }
      .choice-title {
        font-size: 26px;
        font-weight: 750;
        margin-bottom: 10px;
      }
      .choice-subtitle {
        font-size: 15px;
        color: #52627a;
        line-height: 1.35;
      }
      .workflow-grid {
        display: grid;
        grid-template-columns: minmax(320px, 430px) minmax(360px, 1fr);
        gap: 20px;
        align-items: start;
      }
      .panel-card {
        padding: 22px;
      }
      .setup-panel {
        padding: 16px 18px;
        margin: 16px 0;
      }
      .nav-controls {
        display: flex;
        gap: 10px;
        margin-bottom: 16px;
      }
      .primary-btn, .secondary-btn {
        border-radius: 10px !important;
        padding: 10px 16px !important;
      }
      .primary-btn {
        background: #2c7be5 !important;
        border-color: #2c7be5 !important;
        color: #ffffff !important;
        font-weight: 650;
        margin-top: 12px;
      }
      .secondary-btn {
        background: #ffffff !important;
        border: 1px solid #c8d2df !important;
        color: #29405f !important;
      }
      .advanced-panel {
        margin: 18px 0;
        padding: 12px 14px;
        border: 1px solid #d9e2ec;
        border-radius: 12px;
        background: #fbfcfe;
      }
      .advanced-panel summary {
        cursor: pointer;
        font-weight: 700;
      }
      .advanced-body {
        padding-top: 12px;
      }
      .input-note {
        font-size: 14px;
        color: #52627a;
        margin: 10px 0 14px;
      }
      table {
        background: #ffffff;
      }
      .shiny-output-error {
        color: #b42318;
        white-space: normal;
      }
      @media (max-width: 900px) {
        .app-header, .workflow-grid {
          display: block;
        }
        .panel-card {
          margin-bottom: 18px;
        }
      }
    "))
  ),
  div(
    class = "app-shell",
    uiOutput("header_ui"),
    uiOutput("main_ui")
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    authorized = FALSE,
    alpha_error = "",
    design = NULL,
    analysis = NULL,
    cc_input_mode = NULL,
    cc_setup_done = FALSE,
    cc_tests = c("genotypes", "trend")
  )

  reset_state <- function() {
    rv$design <- NULL
    rv$analysis <- NULL
    rv$cc_input_mode <- NULL
    rv$cc_setup_done <- FALSE
    rv$cc_tests <- c("genotypes", "trend")
  }

  observeEvent(input$alpha_key_submit, {
    if (identical(input$alpha_key, "alpha_pawh")) {
      rv$authorized <- TRUE
      rv$alpha_error <- ""
      updateTextInput(session, "alpha_key", value = "")
    } else {
      rv$alpha_error <- "incorrect"
    }
  })

  observeEvent(input$alpha_key, {
    if (nzchar(rv$alpha_error)) {
      rv$alpha_error <- ""
    }
  }, ignoreInit = TRUE)

  observeEvent(input$choose_tdt, {
    rv$design <- "tdt"
    rv$analysis <- NULL
  })

  observeEvent(input$choose_cc, {
    rv$design <- "cc"
    rv$analysis <- NULL
    rv$cc_setup_done <- FALSE
  })

  observeEvent(input$tdt_power_card, rv$analysis <- "power")
  observeEvent(input$tdt_mssn_card, rv$analysis <- "mssn")
  observeEvent(input$tdt_plots_card, rv$analysis <- "plots")
  observeEvent(input$cc_power_card, rv$analysis <- "power")
  observeEvent(input$cc_mssn_card, rv$analysis <- "mssn")
  observeEvent(input$cc_plots_card, rv$analysis <- "plots")

  observeEvent(input$cc_mode_model_based, rv$cc_input_mode <- "model_based")
  observeEvent(input$cc_mode_model_free, rv$cc_input_mode <- "model_free")

  observeEvent(input$cc_setup_continue, {
    if (is.null(rv$cc_input_mode)) {
      showNotification("Choose model-based or model-free input first.", type = "error")
      return()
    }
    tests <- input$cc_tests
    if (is.null(tests) || length(tests) == 0) {
      showNotification("Select at least one case-control test.", type = "error")
      return()
    }
    rv$cc_tests <- tests
    rv$cc_setup_done <- TRUE
  })

  observeEvent(input$back_btn, {
    if (!is.null(rv$analysis)) {
      rv$analysis <- NULL
    } else if (identical(rv$design, "cc") && isTRUE(rv$cc_setup_done)) {
      rv$cc_setup_done <- FALSE
    } else {
      rv$design <- NULL
      rv$cc_setup_done <- FALSE
    }
  })

  observeEvent(input$reset_btn, reset_state())

  output$alpha_error <- renderText(rv$alpha_error)

  output$header_ui <- renderUI({
    if (!isTRUE(rv$authorized)) {
      return(NULL)
    }
    div(
      class = "app-header",
      div(
        h1("Genetic Association Study Power Calculator"),
        div(class = "breadcrumb-line", textOutput("breadcrumb", inline = TRUE))
      )
    )
  })

  output$breadcrumb <- renderText({
    parts <- c("Study Design")
    if (!is.null(rv$design)) {
      parts <- c(parts, if (rv$design == "tdt") "TDT" else "Case-Control")
    }
    if (identical(rv$design, "cc") && !is.null(rv$cc_input_mode)) {
      parts <- c(parts, if (rv$cc_input_mode == "model_based") "Model-Based" else "Model-Free")
    }
    if (!is.null(rv$analysis)) {
      parts <- c(parts, switch(rv$analysis, power = "Power", mssn = "Sample Size", plots = "Plots"))
    }
    paste(parts, collapse = " > ")
  })

  output$main_ui <- renderUI({
    if (!isTRUE(rv$authorized)) {
      return(alpha_gate_ui())
    }
    if (is.null(rv$design)) {
      return(landing_ui())
    }
    if (identical(rv$design, "tdt")) {
      if (is.null(rv$analysis)) {
        return(tdt_choice_ui())
      }
      return(switch(
        rv$analysis,
        power = tdt_power_ui(),
        mssn = tdt_mssn_ui(),
        plots = tdt_plots_ui()
      ))
    }
    if (!isTRUE(rv$cc_setup_done)) {
      return(cc_setup_ui(rv$cc_input_mode))
    }
    if (is.null(rv$analysis)) {
      return(cc_choice_ui(rv$cc_input_mode, rv$cc_tests))
    }
    switch(
      rv$analysis,
      power = cc_power_ui(rv$cc_input_mode, rv$cc_tests),
      mssn = cc_mssn_ui(rv$cc_input_mode, rv$cc_tests),
      plots = cc_plots_ui(rv$cc_input_mode, rv$cc_tests)
    )
  })

  output$cc_geno_error_ui <- renderUI({
    model <- input$cc_geno_misclass
    if (is.null(model) || model == "none") {
      return(div(class = "input-note", "No genotype error parameters needed."))
    }
    if (model == "1p") {
      return(sliderInput("cc_e", "e", min = 0, max = 0.50, value = 0.02, step = 0.005))
    }
    if (model == "2p") {
      return(tagList(
        sliderInput("cc_e1", "e1", min = 0, max = 1, value = 0.01, step = 0.005),
        sliderInput("cc_e2", "e2", min = 0, max = 0.50, value = 0.02, step = 0.005)
      ))
    }
    if (model == "3p") {
      return(tagList(
        sliderInput("cc_e01", "e01", min = 0, max = 1, value = 0.02, step = 0.005),
        sliderInput("cc_e02", "e02", min = 0, max = 0.50, value = 0.01, step = 0.005),
        sliderInput("cc_e03", "e03", min = 0, max = 1, value = 0.005, step = 0.005)
      ))
    }
    tagList(
      div(
        class = "input-note",
        strong("Important: differential genotype error is not null-calibrated."),
        "Reported power is nominal asymptotic power and reported MSSN is nominal MSSN, using the usual chi-square critical value. Different case/control error mechanisms can distort the null distribution and inflate type I error; pawh does not independently recalibrate that distribution."
      ),
      selectInput("cc_diff_source", "diff_source", c("explicit", "case", "ctrl")),
      numericInput("cc_diff_multiplier", "diff_multiplier", 1, min = 0, step = 0.05),
      fluidRow(
        column(6, h4("Case errors"),
               numericInput("cc_case_e01", "case_e01", 0.02, min = 0, max = 1, step = 0.005),
               numericInput("cc_case_e02", "case_e02", 0.01, min = 0, max = 0.5, step = 0.005),
               numericInput("cc_case_e03", "case_e03", 0.005, min = 0, max = 1, step = 0.005)),
        column(6, h4("Control errors"),
               numericInput("cc_ctrl_e01", "ctrl_e01", 0.01, min = 0, max = 1, step = 0.005),
               numericInput("cc_ctrl_e02", "ctrl_e02", 0.005, min = 0, max = 0.5, step = 0.005),
               numericInput("cc_ctrl_e03", "ctrl_e03", 0.002, min = 0, max = 1, step = 0.005))
      )
    )
  })

  output$tdt_plot_xvar_ui <- renderUI({
    choices <- if (identical(input$tdt_plot_type, "mssn")) {
      c("target_power", "alpha", "pd", "prev", "R1", "R2", "delta_prime", "misclass_rate", "heter_rate", "locus_het_rate")
    } else {
      c("N", "alpha", "pd", "prev", "R1", "R2", "delta_prime", "misclass_rate", "heter_rate", "locus_het_rate")
    }
    selectInput("tdt_plot_x_var", "x_var", choices = choices, selected = if ("heter_rate" %in% choices) "heter_rate" else choices[[1]])
  })

  output$cc_plot_xvar_ui <- renderUI({
    choices <- if (identical(input$cc_plot_type, "mssn")) {
      c("power", "alpha", "prev", "pd", "R2", "k", "pi", "locus_het_rate",
        "theta", "phi", "pheno_error_multiplier", "e", "e1", "e2", "e01",
        "e02", "e03", "geno_error_multiplier", "diff_multiplier")
    } else {
      c("N_case", "N_ctrl", "N_total", "alpha", "prev", "pd", "R2", "k", "pi",
        "locus_het_rate", "theta", "phi", "pheno_error_multiplier", "e", "e1",
        "e2", "e01", "e02", "e03", "geno_error_multiplier", "diff_multiplier")
    }
    selectInput("cc_plot_x_var", "x_var", choices = choices, selected = if (identical(input$cc_plot_type, "mssn")) "power" else "N_case")
  })

  tdt_args <- reactive({
    list(
      pd = input$tdt_pd,
      prev = input$tdt_prev,
      R1 = input$tdt_R1,
      R2 = input$tdt_R2,
      alpha = as.numeric(input$tdt_alpha),
      delta_prime = input$tdt_delta_prime,
      misclass_rate = input$tdt_misclass_rate,
      heter_rate = input$tdt_heter_rate
    )
  })

  tdt_power_result <- eventReactive(input$tdt_power_analyze, {
    args <- c(list(N = input$tdt_N), tdt_args(), list(verbose = FALSE))
    safe_call(do.call(pawh::tdt_power, args))
  }, ignoreInit = TRUE)

  output$tdt_power_table <- renderTable({
    payload <- tdt_power_result()
    validate(need(!is.null(payload), "Click Analyze to run TDT power."))
    render_error_or_table(payload, tdt_power_table(payload$result))
  }, striped = TRUE, bordered = TRUE)

  output$tdt_power_plot <- renderPlot({
    payload <- tdt_power_result()
    validate(
      need(!is.null(payload), ""),
      need(is.null(payload$error), payload$error)
    )
    n <- input$tdt_N
    x_values <- unique(round(seq(max(1, 0.25 * n), 2 * n, length.out = 24)))
    print(pawh::plot_tdt_power(
      x_var = "N",
      x_values = x_values,
      scenario = input$tdt_power_plot_scenario,
      N = n,
      pd = input$tdt_pd,
      prev = input$tdt_prev,
      R1 = input$tdt_R1,
      R2 = input$tdt_R2,
      alpha = as.numeric(input$tdt_alpha),
      delta_prime = input$tdt_delta_prime,
      misclass_rate = input$tdt_misclass_rate,
      heter_rate = input$tdt_heter_rate
    ))
  })

  tdt_mssn_result <- eventReactive(input$tdt_mssn_analyze, {
    args <- c(list(target_power = input$tdt_target_power), tdt_args(), list(verbose = FALSE))
    safe_call(do.call(pawh::tdt_mssn, args))
  }, ignoreInit = TRUE)

  output$tdt_mssn_table <- renderTable({
    payload <- tdt_mssn_result()
    validate(need(!is.null(payload), "Click Analyze to run TDT sample size."))
    render_error_or_table(payload, tdt_mssn_table(payload$result))
  }, striped = TRUE, bordered = TRUE)

  output$tdt_mssn_plot <- renderPlot({
    payload <- tdt_mssn_result()
    validate(
      need(!is.null(payload), ""),
      need(is.null(payload$error), payload$error)
    )
    scenario <- input$tdt_mssn_plot_scenario
    x_var <- if (scenario == "misclassification") "misclass_rate" else "heter_rate"
    print(pawh::plot_tdt_mssn(
      x_var = x_var,
      x_values = seq(0, if (x_var == "misclass_rate") 0.20 else 0.80, length.out = 24),
      scenario = scenario,
      target_power = input$tdt_target_power,
      pd = input$tdt_pd,
      prev = input$tdt_prev,
      R1 = input$tdt_R1,
      R2 = input$tdt_R2,
      alpha = as.numeric(input$tdt_alpha),
      delta_prime = input$tdt_delta_prime,
      misclass_rate = input$tdt_misclass_rate,
      heter_rate = input$tdt_heter_rate
    ))
  })

  tdt_plot_result <- eventReactive(input$tdt_plot_generate, {
    safe_call({
      x_values <- make_range(input$tdt_plot_xmin, input$tdt_plot_xmax, input$tdt_plot_xstep)
      common <- list(
        x_var = input$tdt_plot_x_var,
        x_values = x_values,
        scenario = input$tdt_plot_scenario,
        pd = input$tdt_plot_pd,
        prev = input$tdt_plot_prev,
        R1 = input$tdt_plot_R1,
        R2 = input$tdt_plot_R2,
        alpha = as.numeric(input$tdt_plot_alpha),
        delta_prime = input$tdt_plot_delta_prime,
        misclass_rate = input$tdt_plot_misclass_rate,
        heter_rate = input$tdt_plot_heter_rate,
        title = if (nzchar(input$tdt_plot_title)) input$tdt_plot_title else NULL,
        x_label = if (nzchar(input$tdt_plot_x_label)) input$tdt_plot_x_label else NULL,
        y_label = if (nzchar(input$tdt_plot_y_label)) input$tdt_plot_y_label else NULL
      )
      if (identical(input$tdt_plot_type, "power")) {
        do.call(pawh::plot_tdt_power, c(common, list(N = input$tdt_plot_N)))
      } else {
        do.call(pawh::plot_tdt_mssn, c(common, list(target_power = input$tdt_plot_target_power)))
      }
    })
  }, ignoreInit = TRUE)

  output$tdt_free_plot <- renderPlot({
    payload <- tdt_plot_result()
    validate(
      need(!is.null(payload), "Click Generate Plot."),
      need(is.null(payload$error), payload$error)
    )
    print(payload$result)
  })

  output$tdt_plot_error <- renderTable({
    payload <- tdt_plot_result()
    if (is.null(payload) || is.null(payload$error)) {
      return(NULL)
    }
    data.frame(Message = payload$error)
  }, bordered = TRUE)

  validate_cc_freqs <- function(g1, g0) {
    if (any(!is.finite(c(g1, g0)))) {
      stop("Genotype frequencies must be finite.")
    }
    if (any(c(g1, g0) < 0)) {
      stop("Genotype frequencies cannot be negative.")
    }
    if (abs(sum(g1) - 1) > 1e-6) {
      stop("g1 genotype frequencies must sum to 1.")
    }
    if (abs(sum(g0) - 1) > 1e-6) {
      stop("g0 genotype frequencies must sum to 1.")
    }
  }

  cc_args <- reactive({
    selected_tests <- nonempty_tests(rv$cc_tests)
    args <- list(
      alpha = as.numeric(input$cc_alpha %||% input$cc_plot_alpha),
      input_mode = rv$cc_input_mode,
      k = input$cc_k %||% input$cc_plot_k,
      w = c(input$cc_w0, input$cc_w1, input$cc_w2),
      locus_het = isTRUE(input$cc_locus_het),
      pi = 1 - (input$cc_locus_het_rate %||% 0),
      pheno_misclass = isTRUE(input$cc_pheno_misclass),
      theta = input$cc_theta %||% 0,
      phi = input$cc_phi %||% 0,
      geno_misclass = input$cc_geno_misclass %||% "none",
      e = input$cc_e %||% 0,
      e1 = input$cc_e1 %||% 0,
      e2 = input$cc_e2 %||% 0,
      e01 = input$cc_e01 %||% 0,
      e02 = input$cc_e02 %||% 0,
      e03 = input$cc_e03 %||% 0,
      case_e01 = input$cc_case_e01 %||% 0,
      case_e02 = input$cc_case_e02 %||% 0,
      case_e03 = input$cc_case_e03 %||% 0,
      ctrl_e01 = input$cc_ctrl_e01 %||% 0,
      ctrl_e02 = input$cc_ctrl_e02 %||% 0,
      ctrl_e03 = input$cc_ctrl_e03 %||% 0,
      diff_source = input$cc_diff_source %||% "explicit",
      diff_multiplier = input$cc_diff_multiplier %||% 1
    )
    if (identical(rv$cc_input_mode, "model_based")) {
      args$prev <- input$cc_prev
      args$pd <- input$cc_pd
      args$R2 <- input$cc_R2
      args$MOI <- input$cc_MOI
    } else {
      g1 <- c(input$cc_g1_0, input$cc_g1_1, input$cc_g1_2)
      g0 <- c(input$cc_g0_0, input$cc_g0_1, input$cc_g0_2)
      validate_cc_freqs(g1, g0)
      args$g1 <- g1
      args$g0 <- g0
      args$prev <- input$cc_prev_model_free
    }
    args
  })

  cc_power_result <- eventReactive(input$cc_power_analyze, {
    safe_call({
      args <- c(cc_args(), list(N_case = input$cc_N_case, verbose = FALSE))
      do.call(pawh::cc_power, args)
    })
  }, ignoreInit = TRUE)

  output$cc_power_table <- renderTable({
    payload <- cc_power_result()
    validate(need(!is.null(payload), "Click Analyze to run case-control power."))
    render_error_or_table(payload, cc_power_table(payload$result, rv$cc_tests))
  }, striped = TRUE, bordered = TRUE)

  output$cc_power_plot <- renderPlot({
    payload <- cc_power_result()
    validate(
      need(!is.null(payload), ""),
      need(is.null(payload$error), payload$error)
    )
    n <- input$cc_N_case
    args <- cc_args()
    x_values <- unique(round(seq(max(1, 0.25 * n), 2 * n, length.out = 24)))
    plot_args <- c(
      list(
        x_var = "N_case",
        x_values = x_values,
        compare_tests = length(rv$cc_tests) > 1,
        test = rv$cc_tests[[1]],
        N_case = n
      ),
      args
    )
    print(do.call(pawh::plot_cc_power, plot_args))
  })

  cc_mssn_result <- eventReactive(input$cc_mssn_analyze, {
    safe_call({
      args <- c(cc_args(), list(power = input$cc_target_power, verbose = FALSE))
      do.call(pawh::cc_mssn, args)
    })
  }, ignoreInit = TRUE)

  output$cc_mssn_table <- renderTable({
    payload <- cc_mssn_result()
    validate(need(!is.null(payload), "Click Analyze to run case-control sample size."))
    render_error_or_table(payload, cc_mssn_table(payload$result, rv$cc_tests))
  }, striped = TRUE, bordered = TRUE)

  output$cc_mssn_plot <- renderPlot({
    payload <- cc_mssn_result()
    validate(
      need(!is.null(payload), ""),
      need(is.null(payload$error), payload$error)
    )
    args <- cc_args()
    plot_args <- c(
      list(
        x_var = "power",
        x_values = seq(0.50, 0.95, length.out = 24),
        compare_tests = length(rv$cc_tests) > 1,
        test = rv$cc_tests[[1]],
        power = input$cc_target_power
      ),
      args
    )
    print(do.call(pawh::plot_cc_mssn, plot_args))
  })

  cc_plot_result <- eventReactive(input$cc_plot_generate, {
    safe_call({
      selected_tests <- nonempty_tests(rv$cc_tests)
      x_values <- make_range(input$cc_plot_xmin, input$cc_plot_xmax, input$cc_plot_xstep)
      args <- cc_args()
      args$alpha <- as.numeric(input$cc_plot_alpha)
      args$k <- input$cc_plot_k
      common <- c(
        list(
          x_var = input$cc_plot_x_var,
          x_values = x_values,
          compare_tests = length(selected_tests) > 1,
          test = selected_tests[[1]],
          title = if (nzchar(input$cc_plot_title)) input$cc_plot_title else NULL,
          x_label = if (nzchar(input$cc_plot_x_label)) input$cc_plot_x_label else NULL,
          y_label = if (nzchar(input$cc_plot_y_label)) input$cc_plot_y_label else NULL
        ),
        args
      )
      if (identical(input$cc_plot_type, "power")) {
        do.call(pawh::plot_cc_power, c(common, list(N_case = input$cc_plot_N_case)))
      } else {
        do.call(pawh::plot_cc_mssn, c(common, list(power = input$cc_plot_target_power)))
      }
    })
  }, ignoreInit = TRUE)

  output$cc_free_plot <- renderPlot({
    payload <- cc_plot_result()
    validate(
      need(!is.null(payload), "Click Generate Plot."),
      need(is.null(payload$error), payload$error)
    )
    print(payload$result)
  })

  output$cc_plot_error <- renderTable({
    payload <- cc_plot_result()
    if (is.null(payload) || is.null(payload$error)) {
      return(NULL)
    }
    data.frame(Message = payload$error)
  }, bordered = TRUE)
}

shinyApp(ui, server)
