.pawh_test_ui_html <- function() {
  old_options <- options(sass.cache = FALSE)
  on.exit(options(old_options), add = TRUE)
  paste(as.character(pawh:::.pawh_app_ui()), collapse = "\n")
}

test_that("pawh_app constructs an app without launching it", {
  app <- pawh_app()

  expect_s3_class(app, "shiny.appobj")
  expect_false(grepl(
    "runApp",
    paste(deparse(body(pawh_app)), collapse = "\n"),
    fixed = TRUE
  ))
})

test_that("dashboard shell exposes study and workspace navigation", {
  html <- .pawh_test_ui_html()

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

test_that("home study cards request the correct navigation targets", {
  selected <- character()
  shiny::testServer(
    pawh:::.pawh_home_server,
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

test_that("unimplemented study shells remain transparent placeholders", {
  html <- .pawh_test_ui_html()

  expect_match(html, "No calculations are performed by this skeleton", fixed = TRUE)
  expect_match(html, "Case-Control study design", fixed = TRUE)
  expect_match(html, "Your calculated design", fixed = TRUE)
  expect_match(html, "Advanced assumptions", fixed = TRUE)
  expect_match(html, "disabled", fixed = TRUE)
})

test_that("obsolete app entry point and statistical logic are absent", {
  package_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  expect_false(file.exists(file.path(package_root, "inst", "shiny", "pawh_app", "app.R")))
  expect_false(file.exists(file.path(package_root, "inst", "shiny", "pawh_app", "README.md")))

  shiny_files <- list.files(
    file.path(package_root, "R"), pattern = "^shiny-.*\\.R$", full.names = TRUE
  )
  code_lines <- unlist(lapply(shiny_files, readLines, warn = FALSE))
  code <- paste(code_lines[!grepl("^\\s*#", code_lines)], collapse = "\n")
  expect_false(grepl("shiny::runApp", code, fixed = TRUE))
  expect_false(grepl("pchisq\\s*\\(|qchisq\\s*\\(|dmvnorm\\s*\\(", code))
})
