#' Construct the pawh Shiny Dashboard
#'
#' Creates the modular `bslib` dashboard for `pawh` study design. The current
#' application provides the landing page, study navigation, and transparent
#' placeholders for future Case-Control, TDT / Family, and Quantitative Trait
#' workflows. It does not yet perform dashboard calculations.
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
    }
  )
}
