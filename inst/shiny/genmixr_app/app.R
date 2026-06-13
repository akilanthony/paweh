if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The 'shiny' package is required to run this app.", call. = FALSE)
}
if (!requireNamespace("bslib", quietly = TRUE)) {
  stop("The 'bslib' package is required to run this app.", call. = FALSE)
}
if (!requireNamespace("DT", quietly = TRUE)) {
  stop("The 'DT' package is required to run this app.", call. = FALSE)
}
if (!requireNamespace("genmixr", quietly = TRUE)) {
  stop("The 'genmixr' package must be installed to run this app.", call. = FALSE)
}

library(shiny)

fmt_num <- function(x, digits = 5) {
  ifelse(is.na(x), NA, signif(as.numeric(x), digits))
}

capture_function_output <- function(fun, args) {
  printed <- capture.output(
    result <- do.call(fun, args),
    type = "message"
  )
  list(result = result, printed = printed)
}

case_control_summary <- function(out) {
  is_power <- inherits(out, "cc_power_conditional_full")
  rows <- list()

  add_power_row <- function(label, test) {
    if (is.null(test)) return(NULL)
    data.frame(
      Test = label,
      Power = fmt_num(test$power),
      MSSN_case = NA_real_,
      MSSN_ctrl = NA_real_,
      MSSN_total = NA_real_,
      N_case = fmt_num(out$N_case),
      N_ctrl = fmt_num(out$N_ctrl),
      N_total = fmt_num(out$N_total),
      Lambda = fmt_num(test$lambda),
      stringsAsFactors = FALSE
    )
  }

  add_mssn_row <- function(label, test) {
    if (is.null(test)) return(NULL)
    data.frame(
      Test = label,
      Power = NA_real_,
      MSSN_case = fmt_num(test$MSSN_case),
      MSSN_ctrl = fmt_num(test$MSSN_ctrl),
      MSSN_total = fmt_num(test$MSSN_total),
      N_case = NA_real_,
      N_ctrl = NA_real_,
      N_total = NA_real_,
      Lambda = fmt_num(test$lambda_star),
      stringsAsFactors = FALSE
    )
  }

  if (is_power) {
    rows <- list(
      add_power_row("Genotypes", out$tests$genotypes),
      add_power_row("Alleles", out$tests$alleles),
      add_power_row("Trend", out$tests$trend)
    )
  } else {
    rows <- list(
      add_mssn_row("Genotypes", out$tests$genotypes),
      add_mssn_row("Alleles", out$tests$alleles),
      add_mssn_row("Trend", out$tests$trend)
    )
  }

  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

case_control_genotype_freqs <- function(out) {
  data.frame(
    Genotype = c("g0", "g1", "g2"),
    Case = fmt_num(out$freqs$g_obs_case),
    Control = fmt_num(out$freqs$g_obs_ctrl),
    stringsAsFactors = FALSE
  )
}

case_control_allele_freqs <- function(out) {
  if (is.null(out$freqs$p_obs_case) || is.null(out$freqs$p_obs_ctrl)) {
    return(data.frame(
      Allele = character(),
      Case = numeric(),
      Control = numeric()
    ))
  }

  data.frame(
    Allele = c("q", "p"),
    Case = fmt_num(out$freqs$p_obs_case[c("q", "p")]),
    Control = fmt_num(out$freqs$p_obs_ctrl[c("q", "p")]),
    stringsAsFactors = FALSE
  )
}

tdt_summary <- function(out, analysis_type) {
  values <- unlist(out, use.names = TRUE)
  data.frame(
    Metric = names(values),
    Value = fmt_num(values),
    stringsAsFactors = FALSE
  )
}

tdt_plot_args <- function(input) {
  list(
    pd = input$tdt_plot_pd,
    prev = input$tdt_plot_prev,
    R1 = input$tdt_plot_R1,
    R2 = input$tdt_plot_R2,
    alpha = input$tdt_plot_alpha,
    delta_prime = input$tdt_plot_delta_prime,
    N = input$tdt_plot_N,
    misclass_seq = seq(
      input$tdt_plot_misclass_range[1],
      input$tdt_plot_misclass_range[2],
      length.out = input$tdt_plot_points
    ),
    heter_fixed = input$tdt_plot_heter_fixed,
    title = input$tdt_plot_title
  )
}

tdt_plot_call_text <- function(args) {
  paste0(
    "genmixr::tdt_plot_power_misclassification(\n",
    "  N = ", args$N, ",\n",
    "  pd = ", args$pd, ",\n",
    "  prev = ", args$prev, ",\n",
    "  R1 = ", args$R1, ",\n",
    "  R2 = ", args$R2, ",\n",
    "  alpha = ", args$alpha, ",\n",
    "  delta_prime = ", args$delta_prime, ",\n",
    "  misclass_seq = seq(", min(args$misclass_seq), ", ", max(args$misclass_seq),
    ", length.out = ", length(args$misclass_seq), "),\n",
    "  heter_fixed = ", args$heter_fixed, ",\n",
    "  title = ", deparse(args$title), "\n",
    ")"
  )
}

call_tdt_power_misclassification_plot <- function(args) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  genmixr::tdt_plot_power_misclassification(
    pd = args$pd,
    prev = args$prev,
    R1 = args$R1,
    R2 = args$R2,
    alpha = args$alpha,
    delta_prime = args$delta_prime,
    N = args$N,
    misclass_seq = args$misclass_seq,
    heter_fixed = args$heter_fixed,
    title = args$title
  )
}

ui <- bslib::page_navbar(
  title = "genmixr draft app",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  bslib::nav_panel(
    "Help",
    bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::card(
        bslib::card_header("Draft Shiny interface"),
        p("This is a draft Shiny interface for genetic association power and sample size functions in genmixr."),
        p("The app is intended for local exploration and teaching. It prioritizes runnable examples and simple output over final visual design."),
        tags$ul(
          tags$li("Use the Case-control tab for unified case-control power and MSSN calculations."),
          tags$li("Use the TDT tab for a small first-draft interface to selected TDT functions."),
          tags$li("Inputs are validated with friendly error messages where possible.")
        )
      ),
      bslib::card(
        bslib::card_header("About equations"),
        p("TODO: Finalize textbook equation and page references after verification."),
        p("No equations are implemented in the app itself; the app calls exported genmixr package functions.")
      )
    )
  ),

  bslib::nav_panel(
    "Case-control",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 360,
        selectInput("cc_analysis", "Analysis type", c("Power", "MSSN")),
        selectInput("cc_input_mode", "input_mode", c("model_free", "model_based")),

        conditionalPanel(
          "input.cc_input_mode == 'model_free'",
          tags$strong("Model-free genotype frequencies"),
          numericInput("cc_g1_0", "Case g0", 0.25, min = 0, max = 1, step = 0.01),
          numericInput("cc_g1_1", "Case g1", 0.50, min = 0, max = 1, step = 0.01),
          numericInput("cc_g1_2", "Case g2", 0.25, min = 0, max = 1, step = 0.01),
          numericInput("cc_g0_0", "Control g0", 0.36, min = 0, max = 1, step = 0.01),
          numericInput("cc_g0_1", "Control g1", 0.48, min = 0, max = 1, step = 0.01),
          numericInput("cc_g0_2", "Control g2", 0.16, min = 0, max = 1, step = 0.01)
        ),

        conditionalPanel(
          "input.cc_input_mode == 'model_based'",
          tags$strong("Model-based inputs"),
          numericInput("cc_prev", "prev", 0.05, min = 0.0001, max = 0.999, step = 0.01),
          numericInput("cc_pd", "pd", 0.30, min = 0.0001, max = 0.999, step = 0.01),
          numericInput("cc_R2", "R2", 1.80, min = 0.0001, step = 0.1),
          selectInput("cc_MOI", "MOI", c("M", "D", "Rec"))
        ),

        conditionalPanel(
          "input.cc_analysis == 'Power'",
          numericInput("cc_N_case", "N_case", 200, min = 1, step = 10)
        ),
        conditionalPanel(
          "input.cc_analysis == 'MSSN'",
          numericInput("cc_target_power", "Target power", 0.80, min = 0.001, max = 0.999, step = 0.01)
        ),
        numericInput("cc_alpha", "alpha", 0.05, min = 1e-10, max = 0.999, step = 0.01),
        numericInput("cc_k", "k", 1, min = 0.001, step = 0.1),

        tags$strong("Trend weights"),
        numericInput("cc_w0", "w0", 0, step = 1),
        numericInput("cc_w1", "w1", 1, step = 1),
        numericInput("cc_w2", "w2", 2, step = 1),
        checkboxInput("cc_include_allelic", "include_allelic", TRUE),

        tags$strong("Locus heterogeneity"),
        checkboxInput("cc_locus_het", "locus_het", FALSE),
        numericInput("cc_pi", "pi", 0.80, min = 0, max = 1, step = 0.05),

        tags$strong("Genotype misclassification"),
        selectInput("cc_geno_misclass", "geno_misclass", c("none", "1p", "2p", "3p", "diff3p")),
        numericInput("cc_e", "e", 0.02, min = 0, max = 0.5, step = 0.005),
        numericInput("cc_e1", "e1", 0.01, min = 0, max = 1, step = 0.005),
        numericInput("cc_e2", "e2", 0.02, min = 0, max = 0.5, step = 0.005),
        numericInput("cc_e01", "e01", 0.02, min = 0, max = 1, step = 0.005),
        numericInput("cc_e02", "e02", 0.01, min = 0, max = 0.5, step = 0.005),
        numericInput("cc_e03", "e03", 0.005, min = 0, max = 1, step = 0.005),

        conditionalPanel(
          "input.cc_geno_misclass == 'diff3p'",
          tags$strong("Differential 3p controls"),
          selectInput("cc_diff_source", "diff_source", c("explicit", "case", "ctrl")),
          numericInput("cc_diff_multiplier", "diff_multiplier", 0.50, min = 0, step = 0.05),
          numericInput("cc_case_e01", "case_e01", 0.02, min = 0, max = 1, step = 0.005),
          numericInput("cc_case_e02", "case_e02", 0.01, min = 0, max = 0.5, step = 0.005),
          numericInput("cc_case_e03", "case_e03", 0.005, min = 0, max = 1, step = 0.005),
          numericInput("cc_ctrl_e01", "ctrl_e01", 0.01, min = 0, max = 1, step = 0.005),
          numericInput("cc_ctrl_e02", "ctrl_e02", 0.005, min = 0, max = 0.5, step = 0.005),
          numericInput("cc_ctrl_e03", "ctrl_e03", 0.002, min = 0, max = 1, step = 0.005)
        ),
        actionButton("run_cc", "Run case-control")
      ),

      h3("Case-control results"),
      verbatimTextOutput("cc_error"),
      h4("Clean function output"),
      verbatimTextOutput("cc_printed"),
      h4("Key results"),
      DT::DTOutput("cc_summary"),
      h4("Observed genotype frequencies"),
      DT::DTOutput("cc_genotypes"),
      h4("Risk allele frequencies"),
      DT::DTOutput("cc_alleles")
    )
  ),

  bslib::nav_panel(
    "TDT",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 330,
        selectInput("tdt_analysis", "TDT analysis", c("Power from ET/ENT", "Required trios")),
        conditionalPanel(
          "input.tdt_analysis == 'Power from ET/ENT'",
          numericInput("tdt_ET", "Expected transmissions (ET)", 140, min = 0.001, step = 1),
          numericInput("tdt_ENT", "Expected non-transmissions (ENT)", 100, min = 0.001, step = 1)
        ),
        conditionalPanel(
          "input.tdt_analysis == 'Required trios'",
          numericInput("tdt_target_power", "Target power", 0.80, min = 0.001, max = 0.999, step = 0.01),
          numericInput("tdt_pd", "pd", 0.25, min = 0.0001, max = 0.999, step = 0.01),
          numericInput("tdt_prev", "prev", 0.005, min = 0.0001, max = 0.999, step = 0.001),
          numericInput("tdt_R1", "R1", 2, min = 0.0001, step = 0.1),
          numericInput("tdt_R2", "R2", 2, min = 0.0001, step = 0.1),
          numericInput("tdt_delta_prime", "delta_prime", 1, step = 0.1),
          numericInput("tdt_pi", "pi", 1, min = 0, max = 1, step = 0.05)
        ),
        numericInput("tdt_alpha", "alpha", 0.05, min = 1e-10, max = 0.999, step = 0.01),
        actionButton("run_tdt", "Run TDT")
      ),

      h3("TDT results"),
      verbatimTextOutput("tdt_error"),
      h4("Clean function output"),
      verbatimTextOutput("tdt_printed"),
      h4("Key results"),
      DT::DTOutput("tdt_summary")
    )
  ),

  bslib::nav_panel(
    "TDT Plots",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 330,
        numericInput("tdt_plot_N", "Number of affected trios (N)", 600, min = 1, step = 50),
        numericInput("tdt_plot_pd", "Disease allele frequency (pd)", 0.30, min = 0.0001, max = 0.999, step = 0.01),
        numericInput("tdt_plot_prev", "Disease prevalence (prev)", 0.05, min = 0.0001, max = 0.999, step = 0.01),
        numericInput("tdt_plot_R1", "Relative risk R1", 1.50, min = 0.0001, step = 0.1),
        numericInput("tdt_plot_R2", "Relative risk R2", 2.25, min = 0.0001, step = 0.1),
        numericInput("tdt_plot_alpha", "alpha", 0.05, min = 1e-10, max = 0.999, step = 0.01),
        numericInput("tdt_plot_delta_prime", "delta_prime", 1, step = 0.1),
        sliderInput("tdt_plot_misclass_range", "Misclassification range (pi01)", min = 0, max = 0.50, value = c(0, 0.15), step = 0.01),
        numericInput("tdt_plot_points", "Grid points", 16, min = 2, max = 101, step = 1),
        numericInput("tdt_plot_heter_fixed", "Fixed heterogeneity rate", 0, min = 0, max = 0.999, step = 0.01),
        textInput("tdt_plot_title", "Plot title", "TDT power vs misclassification (pi01)"),
        downloadButton("download_tdt_plot", "Download plot")
      ),

      h3("TDT power vs phenotype misclassification"),
      p("This draft plot calls genmixr::tdt_plot_power_misclassification() and varies the phenotype misclassification rate pi01 while holding heterogeneity fixed."),
      p("TODO: Finalize textbook equation/page references after verification."),
      verbatimTextOutput("tdt_plot_error"),
      plotOutput("tdt_power_misclass_plot", height = "480px"),
      h4("Show function call"),
      verbatimTextOutput("tdt_plot_call")
    )
  )
)

server <- function(input, output, session) {
  cc_run <- eventReactive(input$run_cc, {
    args <- list(
      alpha = input$cc_alpha,
      input_mode = input$cc_input_mode,
      locus_het = input$cc_locus_het,
      pi = input$cc_pi,
      k = input$cc_k,
      w = c(input$cc_w0, input$cc_w1, input$cc_w2),
      include_allelic = input$cc_include_allelic,
      geno_misclass = input$cc_geno_misclass,
      e = input$cc_e,
      e1 = input$cc_e1,
      e2 = input$cc_e2,
      e01 = input$cc_e01,
      e02 = input$cc_e02,
      e03 = input$cc_e03,
      case_e01 = input$cc_case_e01,
      case_e02 = input$cc_case_e02,
      case_e03 = input$cc_case_e03,
      ctrl_e01 = input$cc_ctrl_e01,
      ctrl_e02 = input$cc_ctrl_e02,
      ctrl_e03 = input$cc_ctrl_e03,
      diff_source = input$cc_diff_source,
      diff_multiplier = input$cc_diff_multiplier,
      verbose = TRUE
    )

    if (identical(input$cc_input_mode, "model_free")) {
      args$g1 <- c(input$cc_g1_0, input$cc_g1_1, input$cc_g1_2)
      args$g0 <- c(input$cc_g0_0, input$cc_g0_1, input$cc_g0_2)
    } else {
      args$prev <- input$cc_prev
      args$pd <- input$cc_pd
      args$R2 <- input$cc_R2
      args$MOI <- input$cc_MOI
    }

    if (identical(input$cc_analysis, "Power")) {
      args$N_case <- input$cc_N_case
      fun <- genmixr::cc_power_conditional_full
    } else {
      args$power <- input$cc_target_power
      fun <- genmixr::cc_mssn_conditional_full
    }

    tryCatch(
      capture_function_output(fun, args),
      error = function(e) list(error = conditionMessage(e))
    )
  }, ignoreInit = FALSE)

  output$cc_error <- renderText({
    payload <- cc_run()
    if (!is.null(payload$error)) payload$error else ""
  })

  output$cc_printed <- renderPrint({
    payload <- cc_run()
    validate(need(is.null(payload$error), payload$error))
    cat(paste(payload$printed, collapse = "\n"))
  })

  output$cc_summary <- DT::renderDT({
    payload <- cc_run()
    validate(need(is.null(payload$error), payload$error))
    DT::datatable(case_control_summary(payload$result), rownames = FALSE, options = list(dom = "t"))
  })

  output$cc_genotypes <- DT::renderDT({
    payload <- cc_run()
    validate(need(is.null(payload$error), payload$error))
    DT::datatable(case_control_genotype_freqs(payload$result), rownames = FALSE, options = list(dom = "t"))
  })

  output$cc_alleles <- DT::renderDT({
    payload <- cc_run()
    validate(need(is.null(payload$error), payload$error))
    tab <- case_control_allele_freqs(payload$result)
    validate(need(nrow(tab) > 0, "Allelic output is unavailable when include_allelic = FALSE."))
    DT::datatable(tab, rownames = FALSE, options = list(dom = "t"))
  })

  tdt_run <- eventReactive(input$run_tdt, {
    if (identical(input$tdt_analysis, "Power from ET/ENT")) {
      fun <- genmixr::tdt_power_from_ET_ENT
      args <- list(ET = input$tdt_ET, ENT = input$tdt_ENT, alpha = input$tdt_alpha)
    } else {
      fun <- genmixr::tdt_required_trios
      args <- list(
        power = input$tdt_target_power,
        alpha = input$tdt_alpha,
        df = 1,
        pd = input$tdt_pd,
        prev = input$tdt_prev,
        R1 = input$tdt_R1,
        R2 = input$tdt_R2,
        delta_prime = input$tdt_delta_prime,
        pi = input$tdt_pi
      )
    }

    tryCatch(
      capture_function_output(fun, args),
      error = function(e) list(error = conditionMessage(e))
    )
  }, ignoreInit = FALSE)

  output$tdt_error <- renderText({
    payload <- tdt_run()
    if (!is.null(payload$error)) payload$error else ""
  })

  output$tdt_printed <- renderPrint({
    payload <- tdt_run()
    validate(need(is.null(payload$error), payload$error))
    cat(paste(payload$printed, collapse = "\n"))
  })

  output$tdt_summary <- DT::renderDT({
    payload <- tdt_run()
    validate(need(is.null(payload$error), payload$error))
    DT::datatable(tdt_summary(payload$result, input$tdt_analysis), rownames = FALSE, options = list(dom = "t"))
  })

  tdt_plot_payload <- reactive({
    args <- tdt_plot_args(input)
    tryCatch(
      {
        plot_obj <- suppressMessages(suppressWarnings(
          call_tdt_power_misclassification_plot(args)
        ))
        list(args = args, plot = plot_obj)
      },
      error = function(e) list(args = args, error = conditionMessage(e))
    )
  })

  output$tdt_plot_error <- renderText({
    payload <- tdt_plot_payload()
    if (!is.null(payload$error)) payload$error else ""
  })

  output$tdt_power_misclass_plot <- renderPlot({
    payload <- tdt_plot_payload()
    validate(need(is.null(payload$error), payload$error))
    validate(need(inherits(payload$plot, "ggplot"), "The plotting function did not return a ggplot object."))
    print(payload$plot)
  })

  output$tdt_plot_call <- renderText({
    tdt_plot_call_text(tdt_plot_args(input))
  })

  output$download_tdt_plot <- downloadHandler(
    filename = function() {
      "tdt-power-misclassification.png"
    },
    content = function(file) {
      payload <- tdt_plot_payload()
      validate(need(is.null(payload$error), payload$error))
      ggplot2::ggsave(file, plot = payload$plot, width = 7, height = 5, dpi = 150)
    }
  )
}

shinyApp(ui, server)
