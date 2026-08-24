.paweh_test_ui_html <- function() {
  old_options <- options(sass.cache = FALSE)
  on.exit(options(old_options), add = TRUE)
  paste(as.character(paweh:::.paweh_app_ui()), collapse = "\n")
}

test_that("paweh_app constructs an app without launching it", {
  app <- paweh_app()

  expect_s3_class(app, "shiny.appobj")
  old_app_name <- paste0("paw", "h_app")
  expect_true("paweh_app" %in% getNamespaceExports("paweh"))
  expect_false(old_app_name %in% getNamespaceExports("paweh"))
  expect_false(exists(old_app_name, envir = asNamespace("paweh"), inherits = FALSE))
  expect_false(grepl(
    "runApp",
    paste(deparse(body(paweh_app)), collapse = "\n"),
    fixed = TRUE
  ))
})

test_that("dashboard shell exposes study and workspace navigation", {
  html <- .paweh_test_ui_html()

  expect_match(html, ">paweh<", fixed = TRUE)
  expect_false(grepl(paste0(">paw", "h<"), html, fixed = TRUE))
  expect_match(html, "Power and Sample Size for Genetic Association Studies", fixed = TRUE)
  expect_match(html, "Case-Control", fixed = TRUE)
  expect_match(html, "TDT / Family", fixed = TRUE)
  expect_match(html, "Quantitative Trait", fixed = TRUE)

  expect_match(html, "Full continuous trait", fixed = TRUE)
  expect_match(html, "Extreme phenotype sampling", fixed = TRUE)
  expect_match(html, "Multiple quantitative traits", fixed = TRUE)

  for (label in c("Results", "Sensitivity", "Visualize", "Methods")) {
    expect_match(html, label, fixed = TRUE)
  }
})

test_that("dashboard uses compact scientific presentation helpers", {
  colors <- paweh:::.paweh_plot_colors()
  expect_identical(unname(colors[c("cases", "controls")]), c("#3F4850", "#6F879A"))
  expect_identical(unname(colors[c("genotype", "trend")]), c("#3F4850", "#355C7D"))
  expect_identical(unname(colors[c("baseline", "adjusted", "reference")]),
    c("#C7CDD2", "#3F4850", "#7A848C"))
  expect_s3_class(paweh:::.paweh_plot_theme(), "theme")

  css <- paste(as.character(paweh:::.paweh_dashboard_theme()), collapse = "\n")
  expect_match(css, "#F7F8FA", fixed = TRUE)
  expect_match(css, ".paweh-study-card { min-height: 0", fixed = TRUE)
  expect_match(css, ":focus-visible", fixed = TRUE)
  expect_false(grepl("\\.shiny-plot-output\\s*:empty", css))
  expect_match(css, "grid-template-columns: minmax(0, 1fr)", fixed = TRUE)
  expect_false(grepl("#00FFFF|cyan|coral|#FF0000", css, ignore.case = TRUE))
})

test_that("pre-calculation states use consistent concise guidance", {
  expect_match(as.character(paweh:::.paweh_empty_ui("Results")), "Calculate study design", fixed = TRUE)
  expect_match(as.character(paweh:::.paweh_empty_ui("Sensitivity")), "before exploring sensitivity", fixed = TRUE)
  expect_match(as.character(paweh:::.paweh_empty_ui("Visualize")), "study-specific visualizations", fixed = TRUE)
  expect_match(as.character(paweh:::.paweh_empty_ui("Methods")), "analysis specification", fixed = TRUE)
})

test_that("plot outputs are created only after their required state exists", {
  cc_server <- getFromNamespace(".paweh_case_control_server", "paweh")
  tdt_server <- getFromNamespace(".paweh_tdt_server", "paweh")
  qtl_server <- getFromNamespace(".paweh_qtl_server", "paweh")
  shiny::testServer(cc_server, {
    expect_null(output$sensitivity_plot_container)
    expect_null(output$genotype_plot_container)
    expect_match(paste(as.character(output$sensitivity_controls), collapse = "\n"), "before exploring sensitivity", fixed = TRUE)
    expect_match(paste(as.character(output$visualize_intro), collapse = "\n"), "study-specific visualizations", fixed = TRUE)
  })
  shiny::testServer(tdt_server, {
    expect_null(output$sensitivity_plot_container)
    expect_null(output$transmission_plot_container)
    expect_match(paste(as.character(output$sensitivity_controls), collapse = "\n"), "before exploring sensitivity", fixed = TRUE)
    expect_match(paste(as.character(output$visualize_intro), collapse = "\n"), "study-specific visualizations", fixed = TRUE)
  })
  shiny::testServer(qtl_server, {
    expect_null(output$sensitivity_plot_container)
    expect_null(output$visualization_plot_container)
    expect_match(paste(as.character(output$sensitivity_controls), collapse = "\n"), "before exploring sensitivity", fixed = TRUE)
    expect_match(paste(as.character(output$visualize_intro), collapse = "\n"), "study-specific visualizations", fixed = TRUE)
  })
})

test_that("shared display formatting handles non-finite values clearly", {
  expect_identical(paweh:::.paweh_format_percent(c(NA, Inf, .8)), c("not defined", "not defined", "80.0%"))
  expect_identical(paweh:::.paweh_format_count(c(NA, Inf, -Inf, 12)), c("not defined", "Inf", "-Inf", "12"))
})

test_that("home study cards request the correct navigation targets", {
  selected <- character()
  shiny::testServer(
    paweh:::.paweh_home_server,
    args = list(navigate = function(value) selected <<- c(selected, value)),
    {
      session$setInputs(case_control = 1)
      session$flushReact()
      session$setInputs(tdt = 1)
      session$flushReact()
      session$setInputs(qtl = 1)
      session$flushReact()
    }
  )

  expect_identical(selected, c("case_control", "tdt", "qtl"))
})

test_that("all study workspaces are implemented", {
  html <- .paweh_test_ui_html()

  expect_match(html, "Case-Control study design", fixed = TRUE)
  expect_match(html, "TDT / Family study design", fixed = TRUE)
  expect_match(html, "Quantitative Trait study design", fixed = TRUE)
  expect_match(html, "Your calculated design", fixed = TRUE)
  expect_match(html, "Advanced assumptions", fixed = TRUE)
  expect_false(grepl("No calculations are performed by this skeleton", html, fixed = TRUE))
  expect_false(grepl("forthcoming", html, ignore.case = TRUE))
})

test_that("obsolete app entry point and statistical logic are absent", {
  package_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  expect_false(file.exists(file.path(package_root, "inst", "shiny", "paweh_app", "app.R")))
  expect_false(file.exists(file.path(package_root, "inst", "shiny", "paweh_app", "README.md")))

  shiny_files <- list.files(
    file.path(package_root, "R"), pattern = "^shiny-.*\\.R$", full.names = TRUE
  )
  code_lines <- unlist(lapply(shiny_files, readLines, warn = FALSE))
  code <- paste(code_lines[!grepl("^\\s*#", code_lines)], collapse = "\n")
  expect_false(grepl("shiny::runApp", code, fixed = TRUE))
  expect_false(grepl("pchisq\\s*\\(|qchisq\\s*\\(|dmvnorm\\s*\\(", code))
})
