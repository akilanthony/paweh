#' Construct the pawh Shiny Dashboard
#'
#' Creates the modular `bslib` dashboard for `pawh` study design. The
#' Case-Control, TDT / Family, and Quantitative Trait workspaces provide
#' canonical power, minimum-sample-size, sensitivity, and study-specific visual
#' results. The concise results are complemented by collapsed calculation
#' details and reproducible canonical R calls.
#'
#' `pawh_app()` returns the application object without launching it. This makes
#' the app safe to construct in package code, tests, and deployment tooling.
#' Launch it explicitly with [shiny::runApp()].
#'
#' @return A Shiny application object inheriting from `shiny.appobj`.
#'
#' @examples
#' app <- pawh_app()
#' inherits(app, "shiny.appobj")
#' if (interactive()) {
#'   shiny::runApp(app)
#' }
#'
#' @export
pawh_app <- function() {
  shiny::shinyApp(
    ui = .pawh_app_ui(),
    server = function(input, output, session) {
      .pawh_home_server("home", function(selected) {
        bslib::nav_select("main_nav", selected = selected, session = session)
      })
      .pawh_case_control_server("case_control")
      .pawh_tdt_server("tdt")
      .pawh_qtl_server("qtl")
    }
  )
}
