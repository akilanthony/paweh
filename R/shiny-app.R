#' Construct the paweh Shiny Dashboard
#'
#' Creates the modular `bslib` dashboard for `paweh` study design. The
#' Case-Control, TDT / Family, and Quantitative Trait workspaces provide
#' canonical power, minimum-sample-size, sensitivity, and study-specific visual
#' results. The concise results are complemented by collapsed calculation
#' details and reproducible canonical R calls. Dashboard calculations delegate
#' to the same canonical package functions as the programmatic R interface.
#'
#' `paweh_app()` returns the application object without launching it. This makes
#' the app safe to construct in package code, tests, and deployment tooling.
#' Launch it explicitly with [shiny::runApp()].
#'
#' @return A Shiny application object inheriting from `shiny.appobj`.
#'
#' @examples
#' app <- paweh_app()
#' inherits(app, "shiny.appobj")
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
#'
#' @export
paweh_app <- function() {
  shiny::shinyApp(
    ui = .paweh_app_ui(),
    server = function(input, output, session) {
      .paweh_home_server("home", function(selected) {
        bslib::nav_select("main_nav", selected = selected, session = session)
      })
      .paweh_case_control_server("case_control")
      .paweh_tdt_server("tdt")
      .paweh_qtl_server("qtl")
    }
  )
}
