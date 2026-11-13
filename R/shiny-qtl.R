# Quantitative-trait dashboard; a thin wrapper around the canonical QTL functions.

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
      "Enter a continuous, extreme-phenotype, or multiple-trait design and inspect the canonical verbose output."
    ),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        title = "Design setup", open = "desktop", width = 280,
        resizable = FALSE,
        shiny::radioButtons(
          ns("subtype"), "How will the phenotype be analyzed?", c(
            `Full continuous trait` = "continuous",
            `Extreme phenotype sampling` = "extreme",
            `Multiple quantitative traits` = "multivariate"
          )
        ),
        shiny::uiOutput(ns("subtype_help")),
        shiny::radioButtons(ns("objective"), "Objective", c(
          `Estimate power` = "power", `Minimum sample size` = "mssn"
        )),
        shiny::uiOutput(ns("objective_inputs")),
        shiny::uiOutput(ns("core_inputs")),
        shiny::numericInput(
          ns("alpha"), "Significance level (alpha)", 0.05,
          min = 1e-8, max = 0.999999
        ),
        shiny::tags$details(
          class = "paweh-sidebar-section",
          shiny::tags$summary(shiny::strong("Advanced assumptions")),
          shiny::uiOutput(ns("advanced_inputs"))
        ),
        shiny::actionButton(
          ns("calculate"), "Calculate study design",
          class = "btn-primary paweh-calculate"
        ),
        shiny::uiOutput(ns("changed_notice"))
      ),
      bslib::card(
        bslib::card_header(shiny::h3(class = "h5 mb-0", "Output")),
        bslib::card_body(shiny::verbatimTextOutput(ns("verbose_output")))
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
  labels <- paste("Trait", seq_len(p))
  dimnames(matrix_value) <- list(labels, labels)
  matrix_value
}

.paweh_qtl_snapshot <- function(values) {
  if (!values$subtype %in% c("continuous", "extreme", "multivariate")) {
    stop("Choose a valid quantitative-trait workflow.", call. = FALSE)
  }
  if (!values$objective %in% c("power", "mssn")) {
    stop("Choose a valid objective.", call. = FALSE)
  }
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
      if (values$N != floor(values$N)) {
        stop("Sample size must be an integer.", call. = FALSE)
      }
      objective_value <- values$N
    }
    args <- list(
      alpha = values$alpha, qtl_var = values$qtl_var / 100,
      tau = values$tau, pd = values$pd / 100,
      count_method = values$count_method
    )
    if (values$objective == "mssn") {
      args$multiple_of_three <- isTRUE(values$multiple_of_three)
    }
  } else if (values$subtype == "extreme") {
    .paweh_qtl_num(values$qtl_var, "QTL variance explained", 0, 100, TRUE)
    .paweh_qtl_num(values$tau, "Dominance parameter", -Inf, Inf)
    .paweh_qtl_num(values$x_upper, "Upper population tail", 0, 100, TRUE)
    .paweh_qtl_num(values$x_lower, "Lower population tail", 0, 100, TRUE)
    .paweh_qtl_num(values$k, "Lower-to-upper selected ratio", 0, Inf, TRUE)
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
    p <- suppressWarnings(as.integer(values$n_traits))
    if (length(p) != 1L || is.na(p) || !p %in% 2:4) {
      stop("Number of traits must be between 2 and 4.", call. = FALSE)
    }
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
    if (!values$mv_test %in% c("pillai", "threshold_chisq")) {
      stop("Choose a valid joint test.", call. = FALSE)
    }
    if (values$mv_test == "threshold_chisq") {
      for (i in seq_len(p)) {
        .paweh_qtl_num(x_upper[i], paste("Trait", i, "upper tail"), 0, 100, TRUE)
        .paweh_qtl_num(x_lower[i], paste("Trait", i, "lower tail"), 0, 100, TRUE)
      }
      .paweh_qtl_num(values$k, "Lower-to-upper selected ratio", 0, Inf, TRUE)
    }
    if (values$objective == "power") {
      objective_name <- if (values$mv_test == "pillai") {
        "Total sample size"
      } else {
        "Upper-tail selected sample"
      }
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
      args <- c(args, list(
        x_upper = x_upper, x_lower = x_lower, k = values$k
      ))
    }
  }
  list(
    subtype = values$subtype, objective = values$objective,
    objective_value = objective_value, backend_args = args, display = values
  )
}

.paweh_qtl_function <- function(snapshot) {
  switch(snapshot$subtype,
    continuous = if (snapshot$objective == "power") {
      "qtl_anova_power"
    } else {
      "qtl_anova_mssn"
    },
    extreme = if (snapshot$objective == "power") {
      "qtl_threshold_chisq_power"
    } else {
      "qtl_threshold_chisq_mssn"
    },
    multivariate = if (snapshot$objective == "power") {
      "qtl_multivariate_power_full"
    } else {
      "qtl_multivariate_mssn_full"
    }
  )
}

.paweh_qtl_call_args <- function(
    snapshot, args = snapshot$backend_args, verbose = FALSE
) {
  if (snapshot$objective == "power") {
    if (snapshot$subtype == "continuous" ||
        snapshot$subtype == "multivariate" && args$test == "pillai") {
      args$N <- snapshot$objective_value
    } else {
      args$N_case <- snapshot$objective_value
    }
  } else {
    args$power <- snapshot$objective_value
  }
  args$verbose <- verbose
  args
}

.paweh_qtl_call <- function(
    snapshot, args = snapshot$backend_args, verbose = FALSE
) {
  fun <- get(.paweh_qtl_function(snapshot), mode = "function")
  do.call(fun, .paweh_qtl_call_args(snapshot, args, verbose = verbose))
}

.paweh_qtl_calculate <- function(snapshot) {
  result <- NULL
  warnings <- character()
  captured <- utils::capture.output({
    invisible(result <- withCallingHandlers(
      .paweh_qtl_call(snapshot, verbose = TRUE),
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ))
  }, type = "message")
  if (!length(captured)) captured <- utils::capture.output(print(result))
  if (length(warnings)) {
    captured <- c(captured, "", "Warnings:", paste0("- ", unique(warnings)))
  }
  list(
    snapshot = snapshot,
    result = result,
    verbose_output = paste(captured, collapse = "\n")
  )
}

.paweh_qtl_sig <- function(values) serialize(values, NULL)

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
      ns(id), paste("Correlation: Trait", i, "and Trait", j),
      defaults[[id]], -0.99, 0.99, 0.05
    )
  }
  shiny::tagList(controls)
}

.paweh_qtl_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    state <- shiny::reactiveValues(
      calculation = NULL, error = NULL, signature = NULL
    )
    values <- shiny::reactive(.paweh_qtl_values(input))
    changed <- shiny::reactive(
      !is.null(state$signature) &&
        !identical(state$signature, .paweh_qtl_sig(values()))
    )

    output$subtype_help <- shiny::renderUI({
      text <- switch(
        if (is.null(input$subtype)) "continuous" else input$subtype,
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
        shiny::sliderInput(
          ns("target_power"), "Target power (%)", 50, 99, 80
        )
      } else if (subtype == "continuous" ||
                 subtype == "multivariate" && test == "pillai") {
        shiny::numericInput(ns("N"), "Total sample size", 500, 4, NA, 10)
      } else {
        shiny::numericInput(
          ns("N_case"), "Upper-tail selected sample", 150, 1, NA, 5
        )
      }
    })
    output$core_inputs <- shiny::renderUI({
      subtype <- if (is.null(input$subtype)) "continuous" else input$subtype
      common <- shiny::tagList(
        shiny::numericInput(
          ns("pd"), "Modeled-allele frequency (%)", 30, 0.01, 99.99, 0.1
        ),
        shiny::tags$small(
          class = "text-muted",
          "Population frequency of the allele defining the three genotype groups."
        )
      )
      if (subtype == "continuous") {
        shiny::tagList(
          common,
          shiny::numericInput(
            ns("qtl_var"), "Variance explained by the QTL (%)",
            10, 0.01, 99.99, 0.1
          ),
          shiny::tags$small(
            class = "text-muted",
            "Proportion of total trait variance attributable to the modeled QTL."
          ),
          shiny::numericInput(
            ns("tau"), "Dominance parameter", 0, NA, NA, 0.1
          ),
          shiny::tags$small(
            class = "text-muted",
            "Dominance-to-additivity ratio in the canonical Falconer model."
          )
        )
      } else if (subtype == "extreme") {
        shiny::tagList(
          common,
          shiny::numericInput(
            ns("qtl_var"), "Variance explained by the QTL (%)",
            10, 0.01, 99.99, 0.1
          ),
          shiny::numericInput(ns("tau"), "Dominance parameter", 0, NA, NA, 0.1),
          shiny::numericInput(
            ns("x_upper"), "Upper population tail selected (%)",
            10, 0.01, 99.99, 0.5
          ),
          shiny::numericInput(
            ns("x_lower"), "Lower population tail selected (%)",
            10, 0.01, 99.99, 0.5
          ),
          shiny::tags$small(
            class = "text-muted",
            "Tail percentages define standardized-normal population percentiles; the middle is excluded."
          )
        )
      } else {
        p <- if (is.null(input$n_traits)) 2L else as.integer(input$n_traits)
        test <- if (is.null(input$mv_test)) "pillai" else input$mv_test
        shiny::tagList(
          common,
          shiny::selectInput(
            ns("n_traits"), "Number of quantitative traits",
            stats::setNames(2:4, 2:4)
          ),
          shiny::radioButtons(ns("mv_test"), "Joint analysis", c(
            `Joint continuous-trait test (Pillai MANOVA)` = "pillai",
            `Joint extreme-selection test` = "threshold_chisq"
          )),
          .paweh_qtl_trait_inputs(ns, p, threshold = test == "threshold_chisq"),
          shiny::h6("Phenotype correlations"),
          .paweh_qtl_correlation_inputs(ns, p),
          shiny::tags$small(
            class = "text-muted",
            "Pairwise phenotype correlations; the assembled matrix must be positive definite."
          )
        )
      }
    })
    output$advanced_inputs <- shiny::renderUI({
      subtype <- if (is.null(input$subtype)) "continuous" else input$subtype
      if (subtype == "continuous") {
        shiny::tagList(
          shiny::selectInput(
            ns("count_method"), "Genotype-count method",
            c(Rounded = "rounded", Expected = "expected")
          ),
          if (identical(input$objective, "mssn")) shiny::checkboxInput(
            ns("multiple_of_three"),
            "Restrict sample size to multiples of three", TRUE
          )
        )
      } else if (subtype == "extreme" ||
                 identical(input$mv_test, "threshold_chisq")) {
        shiny::numericInput(
          ns("k"), "Lower-tail selected per upper-tail selected",
          1, 0.01, NA, 0.1
        )
      }
    })

    shiny::observeEvent(input$calculate, {
      state$error <- NULL
      submitted <- values()
      tryCatch({
        state$calculation <- .paweh_qtl_calculate(
          .paweh_qtl_snapshot(submitted)
        )
        state$signature <- .paweh_qtl_sig(submitted)
      }, error = function(error) {
        state$calculation <- NULL
        state$error <- paste(
          "This design could not be calculated.", conditionMessage(error)
        )
      })
    }, ignoreInit = TRUE)

    output$changed_notice <- shiny::renderUI({
      if (changed()) shiny::div(
        class = "paweh-changed-notice", role = "status",
        "Inputs have changed. Recalculate to update output."
      )
    })
    output$verbose_output <- shiny::renderText({
      if (!is.null(state$error)) return(state$error)
      if (is.null(state$calculation)) {
        return("Choose a workflow and select Calculate study design.")
      }
      state$calculation$verbose_output
    })

    list(
      calculation = shiny::reactive(state$calculation),
      changed = changed,
      error = shiny::reactive(state$error),
      verbose_output = shiny::reactive({
        if (is.null(state$calculation)) NULL else state$calculation$verbose_output
      })
    )
  })
}
