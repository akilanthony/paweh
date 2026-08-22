# Case-control dashboard; all statistical work delegates to cc_power()/cc_mssn().
.pawh_cc_defaults <- function() list(objective = "power", input_mode = "model_based", N_case = 500, target_power = 80, alpha = .05, k = 1, prev = 10, pd = 30, MOI = "M", R2 = 2, g1_0 = .35, g1_1 = .45, g1_2 = .2, g0_0 = .49, g0_1 = .42, g0_2 = .09, locus_het = FALSE, pi = 100, pheno_misclass = FALSE, theta = 0, phi = 0, genotype_error = FALSE, geno_misclass = "1p", e = 0, e1 = 0, e2 = 0, e01 = 0, e02 = 0, e03 = 0, case_e01 = 0, case_e02 = 0, case_e03 = 0, ctrl_e01 = 0, ctrl_e02 = 0, ctrl_e03 = 0)
.pawh_cc_values <- function(input) {
  shiny::reactiveValuesToList(input)
  d <- .pawh_cc_defaults()
  x <- lapply(names(d), function(n) if (is.null(input[[n]])) d[[n]] else input[[n]])
  names(x) <- names(d)
  x
}
.pawh_cc_num <- function(x, label, lo, hi, open = FALSE) if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < lo || x > hi || (open && (x == lo || x == hi))) stop(label, " is outside its allowed range.", call. = FALSE)
.pawh_cc_direct_ui <- function(ns) shiny::tagList(shiny::h5("Cases"), shiny::numericInput(ns("g1_0"), "0 modeled alleles", .35, 0, 1, .01), shiny::numericInput(ns("g1_1"), "1 modeled allele", .45, 0, 1, .01), shiny::numericInput(ns("g1_2"), "2 modeled alleles", .2, 0, 1, .01), shiny::h5("Controls"), shiny::numericInput(ns("g0_0"), "0 modeled alleles", .49, 0, 1, .01), shiny::numericInput(ns("g0_1"), "1 modeled allele", .42, 0, 1, .01), shiny::numericInput(ns("g0_2"), "2 modeled alleles", .09, 0, 1, .01), shiny::tags$small(class = "text-muted", "Each group must sum to 1."))
.pawh_cc_error_ui <- function(ns, m) {
  p <- function(id, label, max = 100) shiny::numericInput(ns(id), label, 0, 0, max, .1)
  switch(m,
    `1p` = p("e", "Symmetric error (%)", 50),
    `2p` = shiny::tagList(p("e1", "Homozygote to heterozygote (%)"), p("e2", "Heterozygote to either homozygote (%)", 50)),
    `3p` = shiny::tagList(p("e01", "Homozygote to heterozygote (%)"), p("e02", "Heterozygote to either homozygote (%)", 50), p("e03", "Opposite homozygote (%)")),
    diff3p = shiny::tagList(shiny::h6("Cases"), p("case_e01", "Homozygote to heterozygote (%)"), p("case_e02", "Heterozygote to either homozygote (%)", 50), p("case_e03", "Opposite homozygote (%)"), shiny::h6("Controls"), p("ctrl_e01", "Homozygote to heterozygote (%)"), p("ctrl_e02", "Heterozygote to either homozygote (%)", 50), p("ctrl_e03", "Opposite homozygote (%)"), shiny::div(class = "pawh-caution", "Differential error can affect type I error; results are nominal."))
  )
}
.pawh_case_control_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(.pawh_page_heading("Case-Control study design", "Estimate power or minimum sample size from a genetic model or direct genotype probabilities."), bslib::layout_sidebar(sidebar = bslib::sidebar(title = "Design setup", open = "desktop", shiny::radioButtons(ns("objective"), "Objective", c("Estimate power" = "power", "Minimum sample size" = "mssn")), shiny::uiOutput(ns("objective_inputs")), shiny::radioButtons(ns("input_mode"), "Genotype inputs", c("Genetic model" = "model_based", "Direct genotype probabilities" = "model_free")), shiny::uiOutput(ns("genotype_inputs")), shiny::numericInput(ns("alpha"), "Significance level (alpha)", .05, .000001, .999999), shiny::numericInput(ns("k"), "Controls per case", 1, .01, NA, .1), shiny::tags$details(class = "pawh-sidebar-section", shiny::tags$summary(shiny::strong("Advanced assumptions")), shiny::checkboxInput(ns("locus_het"), "Locus heterogeneity"), shiny::conditionalPanel(sprintf("input['%s']", ns("locus_het")), shiny::sliderInput(ns("pi"), "Disease attributable to modeled locus (%)", 0, 100, 100)), shiny::checkboxInput(ns("pheno_misclass"), "Phenotype misclassification"), shiny::conditionalPanel(sprintf("input['%s']", ns("pheno_misclass")), shiny::uiOutput(ns("phenotype_inputs"))), shiny::checkboxInput(ns("genotype_error"), "Genotype misclassification"), shiny::conditionalPanel(sprintf("input['%s']", ns("genotype_error")), shiny::selectInput(ns("geno_misclass"), "Error model", c("One-parameter" = "1p", "Two-parameter" = "2p", "Three-parameter" = "3p", "Differential case/control" = "diff3p")), shiny::uiOutput(ns("error_inputs")))), shiny::actionButton(ns("calculate"), "Calculate", class = "btn-primary w-100"), shiny::uiOutput(ns("changed_notice")), bslib::card(class = "pawh-design-summary", bslib::card_header("Your calculated design"), bslib::card_body(shiny::uiOutput(ns("design_summary"))))), bslib::navset_card_tab(id = ns("section"), bslib::nav_panel("Results", shiny::uiOutput(ns("results"))), bslib::nav_panel("Sensitivity", shiny::uiOutput(ns("sensitivity_controls")), shiny::uiOutput(ns("sensitivity_message")), shiny::plotOutput(ns("sensitivity_plot"), height = "430px")), bslib::nav_panel("Visualize", shiny::uiOutput(ns("visualize_intro")), shiny::plotOutput(ns("genotype_plot"), height = "430px")), bslib::nav_panel("Methods", shiny::uiOutput(ns("methods"))))))
}
.pawh_cc_snapshot <- function(v) {
  .pawh_cc_num(v$alpha, "Significance level", 0, 1, TRUE)
  .pawh_cc_num(v$k, "Controls per case", 0, Inf, TRUE)
  if (v$objective == "power") .pawh_cc_num(v$N_case, "Cases", 0, Inf, TRUE) else .pawh_cc_num(v$target_power, "Target power", 0, 100, TRUE)
  a <- list(alpha = v$alpha, input_mode = v$input_mode, k = v$k, w = c(0, 1, 2))
  if (v$input_mode == "model_based") {
    .pawh_cc_num(v$prev, "Prevalence", 0, 100, TRUE)
    .pawh_cc_num(v$pd, "Allele frequency", 0, 100, TRUE)
    .pawh_cc_num(v$R2, "Relative risk", 0, Inf, TRUE)
    a <- c(a, list(prev = v$prev / 100, pd = v$pd / 100, R2 = v$R2, MOI = v$MOI))
  } else if (v$input_mode == "model_free") {
    g1 <- unlist(v[c("g1_0", "g1_1", "g1_2")])
    g0 <- unlist(v[c("g0_0", "g0_1", "g0_2")])
    if (any(c(g1, g0) < 0) || any(c(g1, g0) > 1) || abs(sum(g1) - 1) > 1e-8 || abs(sum(g0) - 1) > 1e-8) stop("Case and control genotype probabilities must each sum to 1.", call. = FALSE)
    a <- c(a, list(g1 = g1, g0 = g0))
  } else {
    stop("Choose a valid genotype input method.", call. = FALSE)
  }
  pct <- function(n, max = 100) {
    .pawh_cc_num(v[[n]], n, 0, max)
    v[[n]] / 100
  }
  if (v$pheno_misclass && v$input_mode == "model_free") {
    .pawh_cc_num(v$prev, "Prevalence", 0, 100, TRUE)
    a$prev <- v$prev / 100
  }
  m <- if (v$genotype_error) v$geno_misclass else "none"
  ea <- as.list(stats::setNames(rep(0, 12), c("e", "e1", "e2", "e01", "e02", "e03", "case_e01", "case_e02", "case_e03", "ctrl_e01", "ctrl_e02", "ctrl_e03")))
  if (m == "1p") ea$e <- pct("e", 50)
  if (m == "2p") {
    ea$e1 <- pct("e1")
    ea$e2 <- pct("e2", 50)
  }
  if (m == "3p") {
    ea$e01 <- pct("e01")
    ea$e02 <- pct("e02", 50)
    ea$e03 <- pct("e03")
  }
  if (m == "diff3p") for (n in names(ea)[7:12]) ea[[n]] <- pct(n, if (grepl("e02$", n)) 50 else 100)
  if (m == "3p" && ea$e01 + ea$e03 > 1 || m == "diff3p" && (ea$case_e01 + ea$case_e03 > 1 || ea$ctrl_e01 + ea$ctrl_e03 > 1)) stop("Homozygote errors together cannot exceed 100%.", call. = FALSE)
  a <- c(a, list(locus_het = v$locus_het, pi = if (v$locus_het) pct("pi") else 1, pheno_misclass = v$pheno_misclass, theta = if (v$pheno_misclass) pct("theta", 99.99) else 0, phi = if (v$pheno_misclass) pct("phi", 99.99) else 0, geno_misclass = m), ea, list(diff_source = "explicit", diff_multiplier = 1))
  list(objective = v$objective, input_mode = v$input_mode, objective_value = if (v$objective == "power") v$N_case else v$target_power / 100, backend_args = a, display = v)
}
.pawh_cc_call <- function(s, a = s$backend_args) {
  a$verbose <- FALSE
  if (s$objective == "power") {
    a$N_case <- s$objective_value
    do.call(cc_power, a)
  } else {
    a$power <- s$objective_value
    do.call(cc_mssn, a)
  }
}
.pawh_cc_calculate <- function(s) {
  x <- .pawh_cc_call(s)
  a <- s$backend_args
  a$locus_het <- FALSE
  a$pi <- 1
  a$pheno_misclass <- FALSE
  a$theta <- a$phi <- 0
  a$geno_misclass <- "none"
  b <- .pawh_cc_call(s, a)
  list(snapshot = s, adjusted = x, baseline = b, active = list(locus = isTRUE(x$locus_het$enabled) && x$locus_het$pi < 1, phenotype = isTRUE(x$errors$phenotype_misclass$enabled) && (x$errors$phenotype_misclass$theta > 0 || x$errors$phenotype_misclass$phi > 0), genotype = isTRUE(x$errors$genotype_misclass$enabled)))
}
.pawh_cc_sig <- function(v) serialize(v, NULL)
.pawh_cc_pct <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
.pawh_cc_count <- function(x) formatC(x, format = "d", big.mark = ",")
.pawh_cc_result_data <- function(r, o) if (o == "power") data.frame(Test = c("Genotype chi-square", "Trend"), Power = c(r$tests$genotypes$power, r$tests$trend$power)) else data.frame(Test = c("Genotype chi-square", "Trend"), Cases = c(r$tests$genotypes$MSSN_case, r$tests$trend$MSSN_case), Controls = c(r$tests$genotypes$MSSN_ctrl, r$tests$trend$MSSN_ctrl), Total = c(r$tests$genotypes$MSSN_total, r$tests$trend$MSSN_total))
.pawh_cc_results_ui <- function(c) {
  o <- c$snapshot$objective
  card <- function(r, t) {
    z <- .pawh_cc_result_data(r, o)
    rows <- lapply(seq_len(nrow(z)), function(i) shiny::tags$tr(lapply(seq_along(z), function(j) shiny::tags$td(if (j == 1) z[i, j] else if (o == "power") .pawh_cc_pct(z[i, j], 2) else .pawh_cc_count(z[i, j])))))
    bslib::card(bslib::card_header(t), bslib::card_body(shiny::tags$table(class = "table table-striped", shiny::tags$thead(shiny::tags$tr(lapply(names(z), shiny::tags$th))), shiny::tags$tbody(rows))))
  }
  on <- any(unlist(c$active))
  z <- .pawh_cc_result_data(c$adjusted, o)
  txt <- if (o == "power") paste("Power is", .pawh_cc_pct(z$Power[1]), "for the genotype test and", .pawh_cc_pct(z$Power[2]), "for the trend test.") else paste("The tests require", paste(.pawh_cc_count(z$Total), collapse = " and "), "total participants; plan for the prespecified analysis.")
  shiny::tagList(if (on) bslib::layout_column_wrap(width = "360px", card(c$baseline, "No-error design"), card(c$adjusted, "Adjusted design")) else card(c$adjusted, "No-error design"), shiny::div(class = "pawh-interpretation", shiny::h4("Interpretation"), shiny::p(txt)))
}
.pawh_cc_specs <- function(c) {
  s <- c$snapshot
  v <- s$display
  x <- list(alpha = list(label = "Significance level", min = .0001, max = .2, value = v$alpha, key = "alpha"), k = list(label = "Controls per case", min = .25, max = max(4, v$k * 2), value = v$k, key = "k"))
  if (s$objective == "power") x$N_case <- list(label = "Number of cases", min = max(10, s$objective_value / 4), max = s$objective_value * 2, value = s$objective_value, objective = TRUE) else x$power <- list(label = "Target power", min = .5, max = .99, value = s$objective_value, objective = TRUE)
  if (s$input_mode == "model_based") {
    x$prev <- list(label = "Disease prevalence", min = .01, max = .99, value = v$prev / 100, key = "prev")
    x$pd <- list(label = "Allele frequency", min = .01, max = .99, value = v$pd / 100, key = "pd")
    x$R2 <- list(label = "Homozygote relative risk", min = .2, max = max(4, v$R2 * 2), value = v$R2, key = "R2")
  }
  if (c$active$locus) x$pi <- list(label = "Attributable proportion", min = .05, max = 1, value = s$backend_args$pi, key = "pi")
  if (c$active$phenotype) {
    x$theta <- list(label = "Affected-to-control error", min = 0, max = .25, value = s$backend_args$theta, key = "theta")
    x$phi <- list(label = "Unaffected-to-case error", min = 0, max = .25, value = s$backend_args$phi, key = "phi")
  }
  if (c$active$genotype && s$backend_args$geno_misclass == "1p") x$e <- list(label = "Symmetric genotype error", min = 0, max = .25, value = s$backend_args$e, key = "e")
  x
}
.pawh_cc_sensitivity <- function(c, p, range, n = 40L) {
  sp <- .pawh_cc_specs(c)[[p]]
  if (is.null(sp)) stop("Unavailable parameter.")
  d <- do.call(rbind, lapply(seq(range[1], range[2], length.out = n), function(x) {
    s <- c$snapshot
    a <- s$backend_args
    if (isTRUE(sp$objective)) s$objective_value <- x else a[[sp$key]] <- x
    r <- tryCatch(.pawh_cc_call(s, a), error = function(e) NULL)
    y <- if (is.null(r)) c(NA, NA) else if (s$objective == "power") c(r$tests$genotypes$power, r$tests$trend$power) else c(r$tests$genotypes$MSSN_total, r$tests$trend$MSSN_total)
    data.frame(x = x, Test = c("Genotype chi-square", "Trend"), y = y)
  }))
  list(data = d, label = sp$label, baseline_x = sp$value, objective = c$snapshot$objective)
}
.pawh_cc_sensitivity_plot <- function(x) {
  ggplot2::ggplot(x$data, ggplot2::aes(.data$x, .data$y, color = .data$Test, linetype = .data$Test)) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = x$baseline_x, linetype = 3) +
    ggplot2::labs(x = x$label, y = if (x$objective == "power") "Power" else "Minimum total sample size", color = "Test", linetype = "Test") +
    ggplot2::theme_minimal(base_size = 13)
}
.pawh_cc_freqs <- function(c) {
  on <- any(unlist(c$active))
  rs <- if (on) list(c$baseline, c$adjusted) else list(c$adjusted)
  sc <- if (on) c("No-error design", "Adjusted design") else "No-error design"
  do.call(rbind, lapply(seq_along(rs), function(i) rbind(data.frame(Scenario = sc[i], Group = "Cases", Genotype = factor(0:2), Probability = rs[[i]]$freqs$g_obs_case), data.frame(Scenario = sc[i], Group = "Controls", Genotype = factor(0:2), Probability = rs[[i]]$freqs$g_obs_ctrl))))
}
.pawh_cc_genotype_plot <- function(c) {
  ggplot2::ggplot(.pawh_cc_freqs(c), ggplot2::aes(.data$Genotype, .data$Probability, fill = .data$Group)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~Scenario) +
    ggplot2::labs(x = "Number of modeled alleles", y = "Genotype probability", fill = "Group") +
    ggplot2::theme_minimal(base_size = 13)
}
.pawh_case_control_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    st <- shiny::reactiveValues(calculation = NULL, error = NULL, signature = NULL, sensitivity = NULL)
    values <- shiny::reactive(.pawh_cc_values(input))
    changed <- shiny::reactive(!is.null(st$signature) && !identical(st$signature, .pawh_cc_sig(values())))
    output$objective_inputs <- shiny::renderUI(if (is.null(input$objective) || input$objective == "power") shiny::numericInput(ns("N_case"), "Number of cases", 500, 1, NA, 10) else shiny::sliderInput(ns("target_power"), "Target power (%)", 50, 99, 80))
    output$genotype_inputs <- shiny::renderUI(if (is.null(input$input_mode) || input$input_mode == "model_based") shiny::tagList(shiny::numericInput(ns("prev"), "Disease prevalence (%)", 10, .01, 99.99, .1), shiny::numericInput(ns("pd"), "Modeled-allele frequency (%)", 30, .01, 99.99, .1), shiny::selectInput(ns("MOI"), "Inheritance model", c("Multiplicative" = "M", "Dominant" = "D", "Recessive" = "Rec")), shiny::numericInput(ns("R2"), "Homozygote relative risk", 2, .01, NA, .1)) else .pawh_cc_direct_ui(ns))
    output$phenotype_inputs <- shiny::renderUI(shiny::tagList(if (identical(input$input_mode, "model_free")) shiny::numericInput(ns("prev"), "Disease prevalence (%)", 10, .01, 99.99, .1), shiny::numericInput(ns("theta"), "Affected classified as control (%)", 0, 0, 99.9, .1), shiny::numericInput(ns("phi"), "Unaffected classified as case (%)", 0, 0, 99.9, .1)))
    output$error_inputs <- shiny::renderUI(.pawh_cc_error_ui(ns, if (is.null(input$geno_misclass)) "1p" else input$geno_misclass))
    shiny::observeEvent(input$calculate,
      {
        st$error <- NULL
        st$sensitivity <- NULL
        v <- values()
        tryCatch(
          {
            st$calculation <- .pawh_cc_calculate(.pawh_cc_snapshot(v))
            st$signature <- .pawh_cc_sig(v)
          },
          error = function(e) st$error <- paste("This design could not be calculated.", conditionMessage(e))
        )
      },
      ignoreInit = TRUE
    )
    output$changed_notice <- shiny::renderUI(if (changed()) shiny::div(class = "pawh-changed-notice", role = "status", "Inputs changed. Results still show the last calculated design; select Calculate to update."))
    output$design_summary <- shiny::renderUI(if (!is.null(st$error)) {
      shiny::div(class = "pawh-error", role = "alert", st$error)
    } else if (is.null(st$calculation)) {
      shiny::p(class = "text-muted", "Choose assumptions and select Calculate.")
    } else {
      c <- st$calculation
      s <- c$snapshot
      m <- names(which(unlist(c$active)))
      shiny::tags$ul(shiny::tags$li(if (s$objective == "power") paste(.pawh_cc_count(s$objective_value), "cases; estimate power") else paste("Target power", .pawh_cc_pct(s$objective_value))), shiny::tags$li(if (s$input_mode == "model_based") "Genetic model inputs" else "Direct genotype probabilities"), shiny::tags$li(if (length(m)) paste("Active modifiers:", paste(m, collapse = ", ")) else "No active modifiers"))
    })
    output$results <- shiny::renderUI(if (!is.null(st$error)) shiny::div(class = "pawh-error", role = "alert", st$error) else if (is.null(st$calculation)) .pawh_placeholder_ui("Results", "Set up a design and select Calculate.") else .pawh_cc_results_ui(st$calculation))
    output$sensitivity_controls <- shiny::renderUI({
      if (is.null(st$calculation)) {
        return(.pawh_placeholder_ui("Sensitivity", "Calculate a design first."))
      }
      sp <- .pawh_cc_specs(st$calculation)
      ch <- stats::setNames(names(sp), vapply(sp, `[[`, "", "label"))
      z <- sp[[1]]
      shiny::tagList(shiny::selectInput(ns("sensitivity_parameter"), "Parameter", ch), shiny::sliderInput(ns("sensitivity_range"), "Range", z$min, z$max, pmax(z$min, pmin(z$max, c(z$value * .75, z$value * 1.25)))), shiny::actionButton(ns("run_sensitivity"), "Run sensitivity analysis"))
    })
    shiny::observeEvent(input$sensitivity_parameter, {
      shiny::req(st$calculation)
      z <- .pawh_cc_specs(st$calculation)[[input$sensitivity_parameter]]
      shiny::req(z)
      value <- pmax(z$min, pmin(z$max, c(z$value * .75, z$value * 1.25)))
      shiny::updateSliderInput(session, "sensitivity_range",
        min = z$min, max = z$max, value = value,
        step = (z$max - z$min) / 100
      )
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$run_sensitivity,
      {
        shiny::req(st$calculation, input$sensitivity_parameter, input$sensitivity_range)
        st$sensitivity <- .pawh_cc_sensitivity(st$calculation, input$sensitivity_parameter, input$sensitivity_range)
      },
      ignoreInit = TRUE
    )
    output$sensitivity_message <- shiny::renderUI(if (!is.null(st$calculation) && is.null(st$sensitivity)) shiny::p(class = "text-muted", "Each point is a canonical calculation of the frozen design."))
    output$sensitivity_plot <- shiny::renderPlot({
      shiny::req(st$sensitivity)
      .pawh_cc_sensitivity_plot(st$sensitivity)
    })
    output$visualize_intro <- shiny::renderUI(if (is.null(st$calculation)) .pawh_placeholder_ui("Visualize", "Calculate a design first.") else shiny::tagList(shiny::h3("Genotype distributions"), shiny::p("Probabilities returned by the canonical calculation.")))
    output$genotype_plot <- shiny::renderPlot({
      shiny::req(st$calculation)
      .pawh_cc_genotype_plot(st$calculation)
    })
    output$methods <- shiny::renderUI(shiny::tagList(shiny::h3("Methods"), shiny::p("The dashboard calls cc_power() or cc_mssn(). Baselines and sensitivity points are separate canonical calls."), shiny::a(href = "https://akilanthony.github.io/pawh/articles/pawh-02-case-control-study-design.html", target = "_blank", rel = "noopener", "Read the Case-Control vignette")))
    list(calculation = shiny::reactive(st$calculation), changed = changed, error = shiny::reactive(st$error), sensitivity = shiny::reactive(st$sensitivity))
  })
}
